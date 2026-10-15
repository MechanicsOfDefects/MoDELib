/* This file is part of MoDELib, distributed under GNU GPL v2. */
#include <SurfaceMeshActor.h>
#include <GlidePlaneFactory.h>
#include <QGridLayout>
#include <QSignalBlocker>
#include <sstream>
#include <exception>
#include <cstddef>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkProperty.h>
#include <vtkLightCollection.h>

#include <DDqtControls.h>

namespace model
{
    SurfaceMeshActor::SurfaceMeshActor(vtkGenericOpenGLRenderWindow *window, vtkRenderer *rendererIn,
                             const DefectiveCrystal<3> &crystal)
        : renderWindow(window), renderer(rendererIn), defectiveCrystal(crystal),
          showSurface(new QCheckBox("Show deformed surface", this)), meshSizeEdit(new QLineEdit(this)),
          amplification(new QDoubleSpinBox(this)), showEdges(new QCheckBox("Show triangle edges", this)),
          status(new QLabel(this)), clipSurface(new QCheckBox("Clip surface", this)),
          clipCenterEdit(new QLineEdit(this)), clipNormalEdit(new QLineEdit("1 0 0", this)),
          clipStatus(new QLabel(this)), clipPositionSlider(new QSlider(Qt::Horizontal, this)),
          clipOffsetLabel(new QLabel("Offset: 0 b", this)),
          excludeFaces(new controls::GroupBox("Exclude Faces", this)), polyData(vtkSmartPointer<vtkPolyData>::New()),
          clipPlane(vtkSmartPointer<vtkPlane>::New()), clipper(vtkSmartPointer<vtkClipPolyData>::New()),
          shadingMesh(vtkSmartPointer<vtkCleanPolyData>::New()),
          normals(vtkSmartPointer<vtkPolyDataNormals>::New()), mapper(vtkSmartPointer<vtkPolyDataMapper>::New()),
          actor(vtkSmartPointer<vtkActor>::New()), reliefLight(vtkSmartPointer<vtkLight>::New())
    {
        meshSizeEdit->setObjectName("surfaceMeshArea");
        amplification->setObjectName("surfaceDisplacementAmplification");
        showSurface->setObjectName("showSurfaceMesh");
        const double area =
            std::max(1.0, (crystal.ddBase.mesh.xMax() - crystal.ddBase.mesh.xMin()).squaredNorm() / 10000.0);
        meshSizeEdit->setText(QString::number(area, 'g', 8));
        meshSizeEdit->setToolTip(
            "Maximum reference triangle area in b², as in the Fields tab. Smaller values resolve finer steps.");
        amplification->setRange(0.0, 1.e6);
        amplification->setDecimals(3);
        amplification->setValue(1.0);
        amplification->setToolTip("Displayed position = reference position + amplification × displacement. Zero shows "
                                  "the reference surface.");
        auto *layout = new QGridLayout(this);
        auto *surfaceBox = controls::box(
            "Surface", this, "Show the displaced exterior surface, including sharp slip-step walls.", showSurface);
        auto *surfaceLayout = controls::grid(surfaceBox);
        layout->addWidget(surfaceBox, 0, 0);
        controls::field(surfaceLayout, 0, "Mesh area (b²)", meshSizeEdit,
                        "Maximum reference triangle area in b squared. Smaller values give finer detail and more "
                        "computation. Press Enter to rebuild.");
        controls::field(surfaceLayout, 1, "Amplification", amplification,
                        "Displayed position equals reference position plus this factor times displacement. Zero shows "
                        "the reference shape; one is physical scale.");
        surfaceLayout->addWidget(showEdges, 2, 0, 1, 2);
        showEdges->setToolTip("Draw triangle edges over the shaded surface to inspect the mesh resolution.");
        auto *physicsBox =
            controls::box("Physics", this, "Select which simulated microstructures contribute to displacement.");
        surfaceLayout->addWidget(physicsBox, 3, 0, 1, 2);
        int row = 0;
        for(const auto &microstructure : defectiveCrystal.microstructures())
        {
            auto *check = new QCheckBox(QString::fromStdString(microstructure->tag), this);
            check->setChecked(true);
            microstructuresCheck.push_back(check);
            controls::grid(physicsBox)->addWidget(check, row++, 0, 1, 2);
            check->setToolTip("Include this physics contribution in the surface displacement. Changes recompute the "
                              "visible surface.");
            connect(check, &QCheckBox::toggled, this,
                    [this](bool)
                    {
                        dirty = true;
                        validDisplacement = false;
                        if(frameAvailable && showSurface->isChecked())
                        {
                            compute();
                        }
                        render();
                    });
        }
        const Eigen::Vector3d center = 0.5 * (crystal.ddBase.mesh.xMin() + crystal.ddBase.mesh.xMax());
        clipCenterEdit->setText(
            QString("%1 %2 %3").arg(center.x(), 0, 'g', 16).arg(center.y(), 0, 'g', 16).arg(center.z(), 0, 'g', 16));
        clipCenterEdit->setObjectName("surfaceClipCenter");
        clipNormalEdit->setObjectName("surfaceClipNormal");
        clipSurface->setObjectName("clipSurfaceMesh");
        clipCenterEdit->setToolTip("Global x y z in b. Press Enter to apply.");
        clipNormalEdit->setToolTip("Global nx ny nz. Keeps the side toward the normal; reverse it to keep the other "
                                   "side. Press Enter to apply.");
        auto *clipBox = controls::box(
            "Clipping", this, "Trim the displaced surface with a plane. The side toward the normal is retained.",
            clipSurface);
        surfaceLayout->addWidget(clipBox, 5, 0, 1, 2);
        auto *clipLayout = controls::grid(clipBox);
        controls::field(clipLayout, 0, "Center (b)", clipCenterEdit,
                        "Global x y z coordinates of the reference clip-plane center, in b. Press Enter to apply and "
                        "reset the offset.");
        controls::field(clipLayout, 1, "Normal", clipNormalEdit,
                        "Global nx ny nz components of a nonzero normal. It is normalized automatically; reverse it to "
                        "keep the other half-space. Press Enter to apply.");
        row = 2;
        appliedClipCenter = center;
        clipTravel = (crystal.ddBase.mesh.xMax() - crystal.ddBase.mesh.xMin()).norm();
        clipPositionSlider->setObjectName("surfaceClipPosition");
        clipPositionSlider->setRange(-1000, 1000);
        clipPositionSlider->setValue(0);
        clipPositionSlider->setPageStep(100);
        clipPositionSlider->setEnabled(false);
        clipPositionSlider->setToolTip("Move along the applied normal. Zero returns to the entered center; range is "
                                       "one box diagonal in each direction.");
        clipLayout->addWidget(clipOffsetLabel, row++, 0, 1, 2);
        clipLayout->addWidget(clipPositionSlider, row++, 0, 1, 2);
        connect(clipPositionSlider, &QSlider::valueChanged, this, &SurfaceMeshActor::moveClipPlane);
        clipStatus->setWordWrap(true);
        clipLayout->addWidget(clipStatus, row++, 0, 1, 2);
        clipPlane->SetOrigin(center.data());
        clipPlane->SetNormal(1.0, 0.0, 0.0);
        clipper->SetInputData(polyData);
        clipper->SetClipFunction(clipPlane);
        clipper->SetValue(0.0);
        clipper->InsideOutOff();
        connect(clipSurface, &QCheckBox::toggled, this, &SurfaceMeshActor::updateClipping);
        connect(clipCenterEdit, &QLineEdit::returnPressed, this, &SurfaceMeshActor::updateClipping);
        connect(clipNormalEdit, &QLineEdit::returnPressed, this, &SurfaceMeshActor::updateClipping);
        excludeFaces->setObjectName("surfaceExcludeFaces");
        excludeFaces->setCheckable(true);
        excludeFaces->setChecked(false);
        excludeFaces->setToolTip("Omit selected external mesh faces before triangulation. IDs match Mesh > Face IDs. Uncheck this group to restore all faces while retaining the selection.");
        auto *exclusionLayout = new QGridLayout(excludeFaces);
        for(const auto &region : crystal.ddBase.mesh.regions())
        {
            for(const auto &face : region.second->faces())
            {
                const size_t id = face.second->sID;
                if(face.second->isExternal() && !faceExclusionChecks.count(id))
                {
                    auto *check = new QCheckBox(QString("Face %1").arg(id), excludeFaces);
                    check->setObjectName(QString("surfaceExcludeFace_%1").arg(id));
                    check->setToolTip(QString("Exclude external face %1 and all of its boundary triangles and slip-step walls. This is the ID displayed in the Mesh tab.").arg(id));
                    faceExclusionChecks.emplace(id, check);
                }
            }
        }
        int faceIndex = 0;
        for(const auto &entry : faceExclusionChecks)
        {
            exclusionLayout->addWidget(entry.second, faceIndex / 2, faceIndex % 2);
            ++faceIndex;
            connect(entry.second, &QCheckBox::toggled, this, &SurfaceMeshActor::rebuild);
        }
        surfaceLayout->addWidget(excludeFaces, 6, 0, 1, 2);
        connect(excludeFaces, &QGroupBox::toggled, this, &SurfaceMeshActor::rebuild);
        status->setWordWrap(true);
        layout->addWidget(status, 1, 0);
        layout->setRowStretch(2, 1);
        shadingMesh->SetInputData(polyData);
        // Slip-trace cutting duplicates polygon vertices. Weld only exactly
        // coincident display points so smooth regions share point normals;
        // never use a distance tolerance that could close a small slip step.
        shadingMesh->SetTolerance(0.0);
        shadingMesh->PointMergingOn();
        shadingMesh->ConvertPolysToLinesOff();
        shadingMesh->ConvertLinesToPointsOff();
        shadingMesh->ConvertStripsToPolysOff();
        normals->SetInputConnection(shadingMesh->GetOutputPort());
        normals->ComputePointNormalsOn();
        normals->SplittingOn();
        normals->SetFeatureAngle(30.0);
        normals->ComputeCellNormalsOn();
        normals->ConsistencyOn();
        normals->AutoOrientNormalsOff();
        mapper->SetInputConnection(normals->GetOutputPort());
        mapper->ScalarVisibilityOff();
        actor->SetMapper(mapper);
        actor->GetProperty()->SetColor(0.65, 0.65, 0.65);
        actor->GetProperty()->SetAmbient(0.25);
        actor->GetProperty()->SetDiffuse(0.75);
        actor->GetProperty()->SetSpecular(0.0);
        auto *shadingBox = controls::box("Shading", this, "Adjust the surface material lighting without recomputing displacement.");
        surfaceLayout->addWidget(shadingBox, 4, 0, 1, 2);
        auto *ambient = new QDoubleSpinBox(shadingBox);
        auto *diffuse = new QDoubleSpinBox(shadingBox);
        auto *specular = new QDoubleSpinBox(shadingBox);
        ambient->setObjectName("surfaceAmbient");
        diffuse->setObjectName("surfaceDiffuse");
        specular->setObjectName("surfaceSpecular");
        for(auto *coefficient : {ambient, diffuse, specular})
        {
            coefficient->setRange(0.0, 1.0);
            coefficient->setDecimals(2);
            coefficient->setSingleStep(0.05);
        }
        ambient->setValue(actor->GetProperty()->GetAmbient());
        diffuse->setValue(actor->GetProperty()->GetDiffuse());
        specular->setValue(actor->GetProperty()->GetSpecular());
        controls::field(controls::grid(shadingBox), 0, "Ambient", ambient,
                        "Uniform illumination from 0 to 1. Increase to soften dark regions and reduce orientation contrast. Changes apply immediately.");
        controls::field(controls::grid(shadingBox), 1, "Diffuse", diffuse,
                        "Directional matte lighting from 0 to 1. Higher values emphasize surface slope relative to the light. Changes apply immediately.");
        controls::field(controls::grid(shadingBox), 2, "Specular", specular,
                        "Reflective highlight strength from 0 to 1. Zero gives a matte surface; higher values emphasize highlights. Changes apply immediately.");
        const auto updateShading = [this, ambient, diffuse, specular]()
        {
            actor->GetProperty()->SetAmbient(ambient->value());
            actor->GetProperty()->SetDiffuse(diffuse->value());
            actor->GetProperty()->SetSpecular(specular->value());
            render();
        };
        for(auto *coefficient : {ambient, diffuse, specular})
        {
            connect(coefficient, &QDoubleSpinBox::valueChanged, this, updateShading);
        }
        actor->GetProperty()->SetInterpolationToPhong();
        actor->GetProperty()->SetEdgeColor(0.15, 0.15, 0.15);
        actor->SetVisibility(false);
        renderer->AddActor(actor);
        // An oblique camera-relative light reveals relief as the view rotates.
        // Add it only while visible so hidden surfaces do not suppress VTK's
        // automatic default light or change the lighting of the other actors.
        reliefLight->SetLightTypeToCameraLight();
        reliefLight->SetPosition(-1.0, 1.0, 0.8);
        reliefLight->SetFocalPoint(0.0, 0.0, 0.0);
        reliefLight->PositionalOff();
        reliefLight->SetIntensity(0.8);
        status->setText("Enable the surface after loading a frame.");
        connect(meshSizeEdit, &QLineEdit::returnPressed, this, &SurfaceMeshActor::rebuild);
        connect(showSurface, &QCheckBox::toggled, this,
                [this](bool enabled)
                {
                    if(enabled && frameAvailable)
                    {
                        compute();
                    }
                    actor->SetVisibility(enabled && validDisplacement);
                    render();
                });
        connect(showEdges, &QCheckBox::toggled, this,
                [this](bool enabled)
                {
                    actor->GetProperty()->SetEdgeVisibility(enabled);
                    render();
                });
        connect(amplification, &QDoubleSpinBox::valueChanged, this,
                [this](double)
                {
                    if(validDisplacement)
                    {
                        try
                        {
                            applyDisplacement();
                        }
                        catch(const std::exception &error)
                        {
                            reportError(error);
                        }
                    }
                    render();
                });
    }

    SurfaceMeshActor::~SurfaceMeshActor()
    {
        renderer->RemoveLight(reliefLight);
        renderer->RemoveActor(actor);
    }

    void SurfaceMeshActor::updateLighting()
    {
        const bool attached = renderer->GetLights()->IsItemPresent(reliefLight) != 0;
        if(actor->GetVisibility() && !attached)
        {
            renderer->AddLight(reliefLight);
        }
        else if(!actor->GetVisibility() && attached)
        {
            renderer->RemoveLight(reliefLight);
        }
    }

    void SurfaceMeshActor::render()
    {
        updateLighting();
        if(renderWindow->GetReadyForRendering())
        {
            renderWindow->Render();
        }
    }

    void SurfaceMeshActor::moveClipPlane()
    {
        const double offset = clipTravel * clipPositionSlider->value() / 1000.0;
        const Eigen::Vector3d center = appliedClipCenter + offset * appliedClipNormal;
        clipPlane->SetOrigin(center.data());
        clipOffsetLabel->setText(QString("Offset: %1 b").arg(offset, 0, 'g', 6));
        if(validDisplacement)
        {
            normals->Update();
            renderer->ResetCameraClippingRange();
        }
        render();
    }

    void SurfaceMeshActor::updateClipping()
    {
        try
        {
            if(clipSurface->isChecked())
            {
                const auto readVector = [](const QLineEdit *edit) -> Eigen::Vector3d
                {
                    std::istringstream stream(edit->text().toStdString());
                    Eigen::Vector3d value;
                    std::string extra;
                    if(!(stream >> value.x() >> value.y() >> value.z()) || (stream >> extra) || !value.allFinite())
                    {
                        throw std::invalid_argument("Enter three finite numbers for each clip-plane vector.");
                    }
                    return value;
                };
                const Eigen::Vector3d center = readVector(clipCenterEdit);
                Eigen::Vector3d normal = readVector(clipNormalEdit);
                const double scale = normal.cwiseAbs().maxCoeff();
                if(scale == 0.0)
                {
                    throw std::invalid_argument("Clip-plane normal must be nonzero.");
                }
                normal /= scale;
                normal.normalize();
                if(center != appliedClipCenter || normal != appliedClipNormal)
                {
                    const QSignalBlocker blocker(clipPositionSlider);
                    clipPositionSlider->setValue(0);
                }
                appliedClipCenter = center;
                appliedClipNormal = normal;
                clipPlane->SetNormal(normal.data());
                shadingMesh->SetInputConnection(clipper->GetOutputPort());
            }
            else
            {
                shadingMesh->SetInputData(polyData);
            }
            clippingEnabled = clipSurface->isChecked();
            clipPositionSlider->setEnabled(clippingEnabled);
            clipStatus->clear();
            moveClipPlane();
        }
        catch(const std::exception &error)
        {
            const QSignalBlocker blocker(clipSurface);
            clipSurface->setChecked(clippingEnabled);
            clipStatus->setText(QString::fromStdString(error.what()) + " Previous clipping settings retained.");
        }
    }

    void SurfaceMeshActor::reportError(const std::exception &error)
    {
        validDisplacement = false;
        dirty = true;
        actor->SetVisibility(false);
        updateLighting();
        status->setText(QString::fromStdString(error.what()));
    }

    std::set<size_t> SurfaceMeshActor::excludedFaceIDs() const
    {
        std::set<size_t> result;
        if(excludeFaces->isChecked())
        {
            for(const auto &entry : faceExclusionChecks)
            {
                if(entry.second->isChecked())
                {
                    result.insert(entry.first);
                }
            }
        }
        return result;
    }

    void SurfaceMeshActor::rebuild()
    {
        try
        {
            bool ok = false;
            const double area = meshSizeEdit->text().toDouble(&ok);
            if(!ok)
            {
                throw std::invalid_argument("Enter a positive maximum triangle area in b²");
            }
            auto replacement = SurfaceMeshGeometry::build(defectiveCrystal.ddBase.mesh, area, excludedFaceIDs());
            referenceGeometry = std::move(replacement);
            geometry = referenceGeometry;
            displacement.clear();
            validDisplacement = false;
            dirty = true;
            if(showSurface->isChecked() && frameAvailable)
            {
                compute();
            }
            else
            {
                status->setText(QString("%1 reference nodes; %2 triangles. Enable surface to compute displacement.")
                                    .arg(geometry.referencePoints.size())
                                    .arg(geometry.triangles.size()));
            }
        }
        catch(const std::exception &error)
        {
            reportError(error);
        }
        render();
    }

    void SurfaceMeshActor::compute()
    {
        try
        {
            if(!dirty)
            {
                return;
            }
            if(referenceGeometry.referencePoints.empty())
            {
                bool ok = false;
                const double area = meshSizeEdit->text().toDouble(&ok);
                if(!ok)
                {
                    throw std::invalid_argument("Enter a positive maximum triangle area in b²");
                }
                referenceGeometry = SurfaceMeshGeometry::build(defectiveCrystal.ddBase.mesh, area, excludedFaceIDs());
            }
            std::vector<SurfaceMeshGeometry::Trace> traces;
            for(const auto &entry : defectiveCrystal.ddBase.glidePlaneFactory.glidePlanes())
            {
                const auto plane = entry.second.lock();
                if(!plane)
                {
                    continue;
                }
                SurfaceMeshGeometry::Trace trace{plane->P, plane->unitNormal, {}};
                for(const auto &segment : plane->meshIntersections)
                {
                    for(const auto *face : segment->faces)
                    {
                        if(face->isExternal())
                        {
                            trace.faces.insert(face->simplices().begin(), face->simplices().end());
                        }
                    }
                }
                if(!trace.faces.empty())
                {
                    traces.push_back(std::move(trace));
                }
            }
            geometry = referenceGeometry.split(traces);
            // Snapshot UI selections and FEM parents before starting worker threads.
            // The current configuration stays fixed until this synchronous call returns.
            std::vector<const MicrostructureBase<3> *> selected;
            for(size_t m = 0; m < microstructuresCheck.size(); ++m)
            {
                if(microstructuresCheck[m]->isChecked())
                {
                    selected.push_back(defectiveCrystal.microstructures()[m].get());
                }
            }
            std::vector<const MicrostructureBase<3>::ElementType *> elements(geometry.referencePoints.size(), nullptr);
            if(defectiveCrystal.ddBase.fe)
            {
                for(size_t k = 0; k < elements.size(); ++k)
                {
                    elements[k] = &defectiveCrystal.ddBase.fe->elements().at(geometry.parents[k]->xID);
                }
            }
            std::vector<Eigen::Vector3d> next(geometry.referencePoints.size());
            std::exception_ptr failure;
            const auto count = static_cast<std::ptrdiff_t>(next.size());
            // Each worker owns one result slot. Qt, VTK, geometry construction,
            // and rendering remain on the calling GUI thread.
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(count >= 256 && !selected.empty())
#endif
            for(std::ptrdiff_t k = 0; k < count; ++k)
            {
                try
                {
                    const auto selectedDisplacement = [&](const Eigen::Vector3d &point) -> Eigen::Vector3d
                    {
                        Eigen::Vector3d value = Eigen::Vector3d::Zero();
                        for(const auto *microstructure : selected)
                        {
                            value += microstructure->displacement(point, nullptr, elements[k], geometry.parents[k]);
                        }
                        return value;
                    };
                    Eigen::Vector3d u = selectedDisplacement(geometry.samplePoints[k]);
                    const Eigen::Vector3d offset = geometry.samplePoints[k] - geometry.referencePoints[k];
                    if(offset.squaredNorm() > 0.0)
                    {
                        // Both one-sided samples belong to the same worker.
                        u = 2.0 * u - selectedDisplacement(geometry.referencePoints[k] + 2.0 * offset);
                    }
                    if(!u.allFinite())
                    {
                        throw std::runtime_error("Nonfinite displacement on surface; adjust mesh size");
                    }
                    next[k] = u;
                }
                catch(...)
                {
                    // Exceptions cannot escape an OpenMP worker. Preserve the
                    // first failure and report it after all workers have joined.
#ifdef _OPENMP
#pragma omp critical(surfaceMeshDisplacementFailure)
#endif
                    {
                        if(!failure)
                        {
                            failure = std::current_exception();
                        }
                    }
                }
            }
            if(failure)
            {
                std::rethrow_exception(failure);
            }
            displacement = std::move(next);
            applyDisplacement();
            dirty = false;
            validDisplacement = true;
            actor->SetVisibility(showSurface->isChecked());
            updateLighting();
            status->setText(QString("%1 nodes; %2 surface triangles; %3 wall triangles.")
                                .arg(geometry.referencePoints.size())
                                .arg(geometry.triangles.size())
                                .arg(geometry.walls.size()));
        }
        catch(const std::exception &error)
        {
            reportError(error);
        }
    }

    void SurfaceMeshActor::applyDisplacement()
    {
        const auto moved = geometry.displaced(displacement, amplification->value());
        auto points = vtkSmartPointer<vtkPoints>::New();
        points->SetDataTypeToDouble();
        for(const auto &point : moved)
        {
            points->InsertNextPoint(point.data());
        }
        auto cells = vtkSmartPointer<vtkCellArray>::New();
        for(const auto &triangle : geometry.triangles)
        {
            const vtkIdType ids[3] = {static_cast<vtkIdType>(triangle[0]), static_cast<vtkIdType>(triangle[1]),
                                      static_cast<vtkIdType>(triangle[2])};
            cells->InsertNextCell(3, ids);
        }
        for(const auto &triangle : geometry.walls)
        {
            const auto &a = moved[triangle[0]];
            const auto &b = moved[triangle[1]];
            const auto &c = moved[triangle[2]];
            if((b - a).cross(c - a).squaredNorm() > 0.0)
            {
                const vtkIdType ids[3] = {static_cast<vtkIdType>(triangle[0]), static_cast<vtkIdType>(triangle[1]),
                                          static_cast<vtkIdType>(triangle[2])};
                cells->InsertNextCell(3, ids);
            }
        }
        polyData->SetPoints(points);
        polyData->SetPolys(cells);
        polyData->Modified();
        normals->Update();
        renderer->ResetCameraClippingRange();
    }

    void SurfaceMeshActor::updateConfiguration()
    {
        frameAvailable = true;
        dirty = true;
        validDisplacement = false;
        if(showSurface->isChecked())
        {
            compute();
        }
    }
} // namespace model
