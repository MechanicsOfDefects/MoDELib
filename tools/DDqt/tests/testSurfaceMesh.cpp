#include <DislocationDynamics/DislocationDynamicsBase.h>
#include <SurfaceMeshActor.h>
#include <FrameCameraRotation.h>
#include <vtkLightCollection.h>
#include <DislocationDynamics/DislocationLoopPatches.h>
#include <QApplication>
#include <vtkCellData.h>
#include <vtkProperty.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <filesystem>
#include <fstream>
#include <chrono>
#include <iostream>
#include <limits>
#ifdef _OPENMP
#include <omp.h>
#endif

using V = Eigen::Vector3d;
using M = Eigen::Matrix3d;
void check(bool condition, const char *message)
{
    if(!condition)
    {
        throw std::runtime_error(message);
    }
}

int main(int argc, char **argv)
{
    QApplication app(argc, argv);
    try
    {
        {
            auto camera = vtkSmartPointer<vtkCamera>::New();
            const V boxCenter(1, 2, 3);
            camera->SetFocalPoint(101, 202, 303);
            camera->SetPosition(111, 202, 303);
            camera->SetViewUp(0, 0, 1);
            model::FrameCameraRotation rotation;
            rotation.apply(camera, boxCenter, V(0, 0, 2), 0.5, 180);
            check(V(camera->GetPosition()).isApprox(V(1, 12, 3)), "Frame rotation ignored runID or axis normalization");
            rotation.apply(camera, boxCenter, V::UnitZ(), 0.5, 360);
            check(V(camera->GetPosition()).isApprox(V(-9, 2, 3)), "Sparse frame rotation is incorrect");
            rotation.apply(camera, boxCenter, V::UnitZ(), 0.5, 180);
            rotation.apply(camera, boxCenter, V::UnitZ(), 0.5, 180);
            check(V(camera->GetPosition()).isApprox(V(1, 12, 3)), "Backward or repeated frame accumulated rotation");
            rotation.apply(camera, boxCenter, V::UnitZ(), 0.0, 180);
            check(V(camera->GetPosition()).isApprox(V(11, 2, 3)), "Disabling rotation did not restore the view");
            rotation.apply(camera, boxCenter, V::UnitY(), -0.5, 180);
            check(V(camera->GetPosition()).isApprox(V(1, 2, 13)), "Axis change or negative angle is incorrect");
            check(V(camera->GetViewUp()).isApprox(V(-1, 0, 0)), "Camera up vector did not rotate");
            camera->Dolly(2.0);
            rotation.apply(camera, boxCenter, V::UnitY(), -0.5, 360);
            check(std::abs(camera->GetDistance() - 5.0) < 1.e-10, "Frame rotation lost interactive zoom");
            check(V(camera->GetFocalPoint()).isApprox(boxCenter), "Rotation did not target the box center");
            const V offset = V(camera->GetPosition()) - boxCenter;
            const V pan(7, 8, 9);
            const V pannedPosition = boxCenter + offset + pan;
            const V pannedFocal = boxCenter + pan;
            camera->SetPosition(pannedPosition.data());
            camera->SetFocalPoint(pannedFocal.data());
            rotation.apply(camera, boxCenter, V::UnitY(), -0.5, 360);
            check(V(camera->GetFocalPoint()).isApprox(boxCenter) &&
                  V(camera->GetPosition()).isApprox(boxCenter + offset),
                  "Manual panning moved the automatic rotation center");
            rotation.apply(camera, boxCenter, V::UnitY(), -0.5, 0);
            check(V(camera->GetFocalPoint()).isApprox(boxCenter), "Frame zero lost the box center");
        }
        {
            auto camera = vtkSmartPointer<vtkCamera>::New();
            camera->SetPosition(10,0,0);
            camera->SetFocalPoint(0,0,0);
            camera->SetViewUp(0,0,1);
            camera->ParallelProjectionOn();
            model::FrameCameraRotation motion;
            motion.apply(camera,V::Zero(),V::UnitZ(),0,10,1);
            check(V(camera->GetPosition()).norm()<1.e-10, "Camera did not reach center");
            check(camera->GetDistance()>9.99 && !camera->GetParallelProjection(), "Fly-through collapsed view or retained parallel projection");
            motion.apply(camera,V::Zero(),V::UnitZ(),0,15,1);
            check(V(camera->GetPosition()).isApprox(V(-5,0,0)), "Camera did not pass through center");
            motion.apply(camera,V::Zero(),V::UnitZ(),18,5,1);
            motion.apply(camera,V::Zero(),V::UnitZ(),18,5,1);
            check(V(camera->GetPosition()).isApprox(V(0,5,0)), "Combined rotation/travel or repeated frame incorrect");
            motion.apply(camera,V::Zero(),V::UnitZ(),0,15,0);
            check(V(camera->GetPosition()).isApprox(V(10,0,0)), "Disabling camera motion did not restore position");
        }
        check(argc == 2 || argc == 3, "Expected repository root and optional --without-dd");
        const bool withDD=argc==2;
        const std::filesystem::path root(argv[1]);
        const auto dir =
            std::filesystem::temp_directory_path() /
            ("model-surface-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        for(const auto *name : {"inputFiles", "config", "F"})
        {
            std::filesystem::create_directories(dir / name);
        }
        const auto input = dir / "inputFiles";
        std::filesystem::copy_file(root / "Library/Materials/W.txt", input / "W.txt");
        std::filesystem::copy_file(root / "Library/DislocationDynamics/DD.txt", input / "DD.txt");
        std::filesystem::copy_file(root / "Library/Meshes/unitCube24.msh", input / "mesh.msh");
        {
            std::ofstream f(input / "DefectiveCrystal.txt");
            f << "physics=" << (withDD ? "ElasticDeformation DislocationDynamics" : "ElasticDeformation") << ";\nuseFEM=0;\nNsteps=1;\nmaxResolveSteps=0;\ndtMax=1;\n"
              << "outputFrequency=1;\noutputBinary=0;\nstartAtTimeStep=0;\n";
            std::ofstream p(input / "polycrystal.txt");
            p << "materialFile=W.txt;\nmeshFile=mesh.msh;\nabsoluteTemperature=300;\n"
              << "C2G1=1 0 0 0 1 0 0 0 1;\nF=100 0 0 0 100 0 0 0 100;\nX0=0 0 0;\nperiodicFaceIDs=;\n";
            std::ofstream e(input / "ElasticDeformation.txt");
            e << "useElasticDeformationFEM=0;\ninertiaReliefPenaltyFactor=1e-3;\n"
              << "ExternalStress0=0 0 .01 0 0 0;\nExternalStressRate=0 0 .01 0 0 0;\n"
              << "ExternalStrain0=0 0 0 0 0 0;\nExternalStrainRate=0 0 0 0 0 0;\n"
              << "stiffnessRatio=0 0 0 0 0 0;\n";
        }
        {
            model::DefectiveCrystalBase<3> base(dir.string());
            const auto coarse = model::SurfaceMeshGeometry::build(base.mesh, 2500);
            const auto fine = model::SurfaceMeshGeometry::build(base.mesh, 625);
            std::set<size_t> allExternalIDs;
            std::set<const model::Simplex<3, 2> *> firstFaceSimplices;
            size_t firstFaceID = 0;
            for(const auto &region : base.mesh.regions())
            {
                for(const auto &face : region.second->faces())
                {
                    if(face.second->isExternal())
                    {
                        if(allExternalIDs.empty())
                        {
                            firstFaceID = face.second->sID;
                            firstFaceSimplices = face.second->simplices();
                        }
                        allExternalIDs.insert(face.second->sID);
                    }
                }
            }
            auto selectedIDs=allExternalIDs;
            selectedIDs.erase(firstFaceID);
            const auto excluded = model::SurfaceMeshGeometry::build(base.mesh, 2500, selectedIDs);
            check(excluded.triangles.size() == coarse.triangles.size() - firstFaceSimplices.size(),
                  "Exclusion did not remove the Mesh-tab face's actual simplices");
            for(const auto *face : excluded.boundaryFaces)
            {
                check(!firstFaceSimplices.count(face), "Excluded simplex remains in the surface");
            }
            check(model::SurfaceMeshGeometry::build(base.mesh, 2500, {}).triangles.empty(),
                  "Excluding every face did not produce an empty surface");
            check(coarse.triangles.size() == 24 && fine.triangles.size() == 96, "Incorrect boundary refinement");
            // Nonaligned, aligned, and intersecting traces must produce true seams.
            for(const double offset : {0.0, 7.0})
            {
                model::SurfaceMeshGeometry::Trace xTrace{V(offset, 0, 0), V::UnitX(), {}};
                model::SurfaceMeshGeometry::Trace yTrace{V(0, -9, 0), V::UnitY(), {}};
                for(const auto *face : fine.boundaryFaces)
                {
                    xTrace.faces.insert(face);
                    yTrace.faces.insert(face);
                }
                const auto cut = fine.split({xTrace, yTrace, xTrace});
                check(!cut.walls.empty(), "Missing slip walls");
                // Exercise the actual solid-angle branch jump, including its
                // relative on-plane tolerance, rather than only a synthetic step.
                model::DislocationLoopPatches<3> patch(nullptr);
                patch.globalPatches()[nullptr] = {V(offset, -100, -100), V(offset, 100, -100), V(offset, 100, 100),
                                                  V(offset, -100, 100)};
                bool checkedJump = false;
                for(const auto &wall : cut.walls)
                {
                    const size_t a = wall[1], b = wall[2];
                    if((cut.referencePoints[a] - cut.referencePoints[b]).norm() < 1.e-8 &&
                       std::abs(cut.referencePoints[a].x() - offset) < 1.e-8 &&
                       (cut.samplePoints[a].x() - offset) * (cut.samplePoints[b].x() - offset) < 0.0)
                    {
                        const double jump =
                            std::abs(patch.solidAngle(cut.samplePoints[a]) - patch.solidAngle(cut.samplePoints[b])) /
                            (4 * std::acos(-1.0));
                        check(jump > 0.9 && jump < 1.1, "One-sided samples lost the solid-angle displacement jump");
                        checkedJump = true;
                    }
                }
                check(checkedJump, "No paired trace samples tested");
                double cutArea = 0;
                for(const auto &t : cut.triangles)
                {
                    const V a = cut.referencePoints[t[0]], b = cut.referencePoints[t[1]], c = cut.referencePoints[t[2]];
                    const V cross = (b - a).cross(c - a);
                    cutArea += cross.norm() / 2;
                    check(cross.dot((a + b + c) / 3) > 0, "Cut reversed or degenerated surface triangle");
                    for(const auto &trace : {xTrace, yTrace})
                    {
                        double lo = 0, hi = 0;
                        for(const auto id : t)
                        {
                            const double d = trace.normal.dot(cut.referencePoints[id] - trace.point);
                            lo = std::min(lo, d);
                            hi = std::max(hi, d);
                        }
                        check(lo > -1.e-8 || hi < 1.e-8, "Triangle bridges a slip trace");
                    }
                }
                check(std::abs(cutArea - 60000) < 1.e-7, "Trace cuts changed exterior area");
                std::vector<V> jump;
                for(const auto &sample : cut.samplePoints)
                {
                    jump.push_back(V(0, 0, (sample.x() > offset ? 1.0 : 0.0) + (sample.y() > -9 ? 2.0 : 0.0)));
                }
                const auto moved = cut.displaced(jump, 10);
                double wallArea = 0;
                double topWallArea = 0;
                for(const auto &t : cut.walls)
                {
                    const V a = cut.referencePoints[t[0]], b = cut.referencePoints[t[1]], c = cut.referencePoints[t[2]];
                    check((b - a).cross(c - a).norm() < 1.e-8, "Wall has finite reference width");
                    const double area = (moved[t[1]] - moved[t[0]]).cross(moved[t[2]] - moved[t[0]]).norm() / 2;
                    wallArea += area;
                    if(std::abs(a.z() - 50) < 1.e-8 && std::abs(b.z() - 50) < 1.e-8 && std::abs(c.z() - 50) < 1.e-8)
                    {
                        topWallArea += area;
                    }
                }
                check(wallArea > 0, "Step walls collapsed after displacement");
                check(std::abs(topWallArea - 3000) < 1.e-6, "Missing or duplicate wall segments at trace crossings");
            }
            std::map<std::pair<size_t, size_t>, size_t> edges;
            double area = 0;
            for(const auto &t : fine.triangles)
            {
                const V a = fine.referencePoints[t[0]], b = fine.referencePoints[t[1]], c = fine.referencePoints[t[2]];
                const V cross = (b - a).cross(c - a);
                check(cross.dot((a + b + c) / 3) > 0, "Inward surface triangle");
                check(cross.norm() / 2 <= 625 * (1 + 1.e-12), "Triangle area exceeds requested limit");
                area += cross.norm() / 2;
                for(int k = 0; k < 3; ++k)
                {
                    const auto endpoints = std::minmax(t[k], t[(k + 1) % 3]);
                    ++edges[{endpoints.first, endpoints.second}];
                }
            }
            check(std::abs(area - 60000) < 1.e-8, "Incorrect exterior area or internal faces included");
            for(const auto &edge : edges)
            {
                check(edge.second == 2, "Crack or duplicate triangles in reference surface");
            }
            for(const auto &p : fine.referencePoints)
            {
                check(std::abs(p.cwiseAbs().maxCoeff() - 50) < 1.e-10, "Surface point off mesh boundary");
            }
            std::vector<V> step;
            for(const auto &p : fine.referencePoints)
            {
                step.push_back(V(0, 0, p.x() > 0 ? 1 : 0));
            }
            const auto stepped = fine.displaced(step, 10);
            const auto reference = fine.displaced(step, 0);
            for(size_t k = 0; k < step.size(); ++k)
            {
                check((stepped[k] - fine.referencePoints[k]).isApprox(10 * step[k]), "Step amplitude lost");
                check(reference[k].isApprox(fine.referencePoints[k]), "Zero amplification did not restore surface");
            }
            for(double badArea : {0.0, -1.0, 1.e-20})
            {
                bool rejected = false;
                try
                {
                    model::SurfaceMeshGeometry::build(base.mesh, badArea);
                }
                catch(const std::invalid_argument &)
                {
                    rejected = true;
                }
                check(rejected, "Invalid or excessive refinement accepted");
            }
            model::DefectiveCrystal<3> crystal(base);
            model::DefectiveCrystalConfigIO<3> config((dir / "config").string());
            crystal.initializeConfiguration(config);
            auto window = vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
            window->SetReadyForRendering(false);
            auto renderer = vtkSmartPointer<vtkRenderer>::New();
            window->AddRenderer(renderer);
            const int originalLights = renderer->GetLights()->GetNumberOfItems();
            {
                model::SurfaceMeshActor surface(window, renderer, crystal);
                auto *ambient = surface.findChild<QDoubleSpinBox*>("surfaceAmbient");
                auto *diffuse = surface.findChild<QDoubleSpinBox*>("surfaceDiffuse");
                auto *specular = surface.findChild<QDoubleSpinBox*>("surfaceSpecular");
                check(ambient && diffuse && specular, "Missing shading controls");
                check(ambient->value() == 0.25 && diffuse->value() == 0.75 && specular->value() == 0.0,
                      "Incorrect shading defaults");
                const auto meshTime = surface.polyData->GetMTime();
                ambient->setValue(0.4);
                diffuse->setValue(0.5);
                specular->setValue(0.2);
                check(surface.actor->GetProperty()->GetAmbient() == 0.4 &&
                      surface.actor->GetProperty()->GetDiffuse() == 0.5 &&
                      surface.actor->GetProperty()->GetSpecular() == 0.2,
                      "Shading controls did not update the material immediately");
                check(surface.polyData->GetMTime() == meshTime, "Shading change rebuilt the surface");
                surface.meshSizeEdit->setText("625");
                surface.externalFaces->setChecked(true);
                surface.updateConfiguration();
                check(surface.polyData->GetNumberOfPoints() == 0, "Hidden surface computed eagerly");
                // Cut polygons have separate vertex IDs, even away from a step.
                // Rendering must smooth coincident neighbors without welding a tiny gap.
                auto sample = vtkSmartPointer<vtkPolyData>::New();
                auto samplePoints = vtkSmartPointer<vtkPoints>::New();
                samplePoints->SetDataTypeToDouble();
                for(const V& point : {V(0,0,0), V(1,0,0), V(0,1,0),
                                      V(1,0,0), V(0,0,0), V(0,-1,0.02),
                                      V(0,0,1.e-9), V(1,0,1.e-9), V(0,1,1.e-9)})
                {
                    samplePoints->InsertNextPoint(point.data());
                }
                auto sampleCells = vtkSmartPointer<vtkCellArray>::New();
                for(vtkIdType first = 0; first < 9; first += 3)
                {
                    const vtkIdType ids[3] = {first, first + 1, first + 2};
                    sampleCells->InsertNextCell(3, ids);
                }
                sample->SetPoints(samplePoints);
                sample->SetPolys(sampleCells);
                surface.shadingMesh->SetInputData(sample);
                surface.normals->Update();
                check(surface.shadingMesh->GetOutput()->GetNumberOfPoints() == 7,
                      "Shading failed to merge duplicates or closed a small step");
                auto* shaded = surface.normals->GetOutput();
                bool smooth = false;
                for(vtkIdType k = 0; k < shaded->GetNumberOfPoints(); ++k)
                {
                    double xyz[3], normal[3];
                    shaded->GetPoint(k, xyz);
                    if(V(xyz).norm() == 0.0)
                    {
                        shaded->GetPointData()->GetNormals()->GetTuple(k, normal);
                        smooth = normal[1] > 0.001 && normal[1] < 0.019;
                    }
                }
                check(smooth, "Shading retained a normal jump across nearly coplanar triangles");
                surface.shadingMesh->SetInputData(surface.polyData);
                surface.showSurface->setChecked(true);
                check(surface.actor->GetVisibility() && surface.polyData->GetNumberOfCells() == 96,
                      "Surface tab did not build after enable");
                check(surface.normals->GetOutput()->GetCellData()->GetNormals() != nullptr, "Missing updated normals");
                check(renderer->GetLights()->IsItemPresent(surface.reliefLight), "Missing surface relief light");
                for(size_t k = 0; k < fine.referencePoints.size(); ++k)
                {
                    double xyz[3];
                    surface.polyData->GetPoint(k, xyz);
                    const V expected =
                        fine.referencePoints[k] +
                        surface.amplification->value() *
                            crystal.displacement(fine.referencePoints[k], nullptr, nullptr, fine.parents[k]);
                    check(V(xyz[0], xyz[1], xyz[2]).isApprox(expected, 1.e-12),
                          "Surface did not use displacement field");
                }
                check(surface.microstructuresCheck.size() == crystal.microstructures().size(),
                      "Missing physics selectors");
                auto *physics = surface.microstructuresCheck.front();
                check(physics->isChecked() && physics->text().toStdString() == crystal.microstructures().front()->tag,
                      "Physics selector does not match Fields defaults");
                physics->setChecked(false);
                for(size_t k = 0; k < fine.referencePoints.size(); ++k)
                {
                    double xyz[3];
                    surface.polyData->GetPoint(k, xyz);
                    check(V(xyz[0], xyz[1], xyz[2]).isApprox(fine.referencePoints[k]),
                          "Unselected physics still contributes displacement");
                }
                surface.updateConfiguration();
                check(!physics->isChecked(), "Frame update lost physics selection");
                surface.showSurface->setChecked(false);
                physics->setChecked(true);
                surface.showSurface->setChecked(true);
                for(size_t k = 0; k < fine.referencePoints.size(); ++k)
                {
                    double xyz[3];
                    surface.polyData->GetPoint(k, xyz);
                    const V expected =
                        fine.referencePoints[k] +
                        surface.amplification->value() *
                            crystal.displacement(fine.referencePoints[k], nullptr, nullptr, fine.parents[k]);
                    check(V(xyz[0], xyz[1], xyz[2]).isApprox(expected, 1.e-12),
                          "Hidden physics selection did not refresh on enable");
                }
                double first[3], second[3];
                surface.polyData->GetPoint(0, first);
                base.simulationParameters.totalTime = 1.0;
                surface.updateConfiguration();
                surface.polyData->GetPoint(0, second);
                check((V(first[0], first[1], first[2]) - V(second[0], second[1], second[2])).norm() > 0,
                      "Frame update did not recompute surface displacement");
                surface.updateConfiguration();
                surface.polyData->GetPoint(0, first);
                check(V(first[0], first[1], first[2]).isApprox(V(second[0], second[1], second[2])),
                      "Repeated frame accumulated deformation");
                surface.amplification->setValue(0);
                surface.polyData->GetPoint(0, first);
                check(V(first[0], first[1], first[2]).isApprox(fine.referencePoints[0]),
                      "Amplification control failed");
                check(!surface.clipSurface->isChecked(), "Clipping should be disabled initially");
                surface.clipSurface->setChecked(true);
                const auto checkClipped = [&](const V &center, const V &normal)
                {
                    surface.normals->Update();
                    auto *output = surface.normals->GetOutput();
                    check(output->GetNumberOfPoints() > 0, "Clipping unexpectedly removed the entire surface");
                    for(vtkIdType k = 0; k < output->GetNumberOfPoints(); ++k)
                    {
                        double xyz[3];
                        output->GetPoint(k, xyz);
                        check(normal.dot(V(xyz[0], xyz[1], xyz[2]) - center) >= -1.e-6,
                              "Rendered surface contains points on discarded side");
                    }
                };
                checkClipped(V::Zero(), V::UnitX());
                surface.clipCenterEdit->setText("10 0 0");
                surface.clipNormalEdit->setText("-2 0 0");
                surface.clipNormalEdit->returnPressed();
                checkClipped(V(10, 0, 0), -V::UnitX());
                surface.updateConfiguration();
                checkClipped(V(10, 0, 0), -V::UnitX());
                check(surface.clipPositionSlider->isEnabled(), "Clip slider not enabled with clipping");
                surface.clipPositionSlider->setValue(100);
                const double offset = (base.mesh.xMax() - base.mesh.xMin()).norm() * 0.1;
                const V shifted = V(10, 0, 0) - offset * V::UnitX();
                check(V(surface.clipPlane->GetOrigin()).isApprox(shifted), "Slider moved off the clip normal");
                check(surface.clipCenterEdit->text() == "10 0 0", "Slider overwrote reference center");
                checkClipped(shifted, -V::UnitX());
                surface.updateConfiguration();
                check(V(surface.clipPlane->GetOrigin()).isApprox(shifted), "Frame update lost slider offset");
                surface.clipSurface->setChecked(false);
                check(!surface.clipPositionSlider->isEnabled(), "Disabled clip slider remains active");
                surface.clipSurface->setChecked(true);
                check(V(surface.clipPlane->GetOrigin()).isApprox(shifted), "Clip toggle lost slider offset");
                surface.clipNormalEdit->setText("1 1 0");
                surface.clipNormalEdit->returnPressed();
                check(surface.clipPositionSlider->value() == 0, "New plane did not reset slider offset");
                surface.clipPositionSlider->setValue(-100);
                const V diagonalNormal = V(1, 1, 0).normalized();
                check(V(surface.clipPlane->GetOrigin()).isApprox(V(10, 0, 0) - offset * diagonalNormal),
                      "Slider failed for oblique normal");
                surface.clipNormalEdit->setText("-2 0 0");
                surface.clipNormalEdit->returnPressed();
                surface.clipNormalEdit->setText("0 0 0");
                surface.updateClipping();
                check(!surface.clipStatus->text().isEmpty() && surface.clipSurface->isChecked(),
                      "Invalid clipping normal was not rejected");
                checkClipped(V(10, 0, 0), -V::UnitX());
                surface.clipSurface->setChecked(false);
                check(surface.normals->GetOutput()->GetNumberOfCells() == 96,
                      "Disabling clipping did not restore complete surface");

                surface.showSurface->setChecked(false);
                surface.updateConfiguration();
                check(!surface.actor->GetVisibility(), "Hidden surface remains visible");
                check(!renderer->GetLights()->IsItemPresent(surface.reliefLight),
                      "Hidden surface changed scene lighting");
                surface.showSurface->setChecked(true);
                check(surface.actor->GetVisibility(), "Surface failed to refresh after re-enable");
                surface.meshSizeEdit->setText("invalid");
                surface.rebuild();
                check(!surface.actor->GetVisibility(), "Invalid mesh input left stale surface visible");
                surface.meshSizeEdit->setText("2500");
                surface.rebuild();
                check(surface.actor->GetVisibility() && surface.polyData->GetNumberOfCells() == 24,
                      "Surface did not recover after correcting mesh input");
            }
            if(withDD)
            {
                auto plane = base.dislocationDynamics().glidePlaneFactory.getFromKey(
                    model::GlidePlaneKey<3>(V::Zero(), base.poly.grain(1)->reciprocalLatticeDirection(V(1, 1, 0))));
                model::SurfaceMeshActor surface(window, renderer, crystal);
                surface.meshSizeEdit->setText("625");
                surface.externalFaces->setChecked(true);
                surface.updateConfiguration();
                surface.showSurface->setChecked(true);
                check(surface.actor->GetVisibility(), "Factory trace extraction failed");
                check(surface.polyData->GetNumberOfPoints() > static_cast<vtkIdType>(fine.referencePoints.size()),
                      "Live glide plane did not split surface");
                plane.reset();
                surface.updateConfiguration();
                check(surface.polyData->GetNumberOfPoints() == static_cast<vtkIdType>(fine.referencePoints.size()),
                      "Expired glide plane cuts survived frame update");
            }
            if(withDD)
            {
#ifdef _OPENMP
                const int savedThreads = omp_get_max_threads();
                const int savedDynamic = omp_get_dynamic();
                omp_set_dynamic(0);
                omp_set_num_threads(1);
#endif
                // Enough samples to exercise the parallel branch, including
                // one-sided samples where a live glide plane cuts the surface.
                auto plane = base.dislocationDynamics().glidePlaneFactory.getFromKey(
                    model::GlidePlaneKey<3>(V::Zero(), base.poly.grain(1)->reciprocalLatticeDirection(V(1, 1, 0))));
                model::SurfaceMeshActor surface(window, renderer, crystal);
                surface.meshSizeEdit->setText("39.0625");
                surface.externalFaces->setChecked(true);
                surface.updateConfiguration();
                surface.showSurface->setChecked(true);
                check(surface.actor->GetVisibility() && surface.polyData->GetNumberOfPoints() >= 256,
                      "Parallel sampling fixture failed");
                std::vector<V> serial;
                for(vtkIdType k = 0; k < surface.polyData->GetNumberOfPoints(); ++k)
                {
                    double xyz[3];
                    surface.polyData->GetPoint(k, xyz);
                    serial.emplace_back(xyz[0], xyz[1], xyz[2]);
                }
#ifdef _OPENMP
                omp_set_num_threads(4);
#endif
                surface.updateConfiguration();
                check(surface.actor->GetVisibility() &&
                          surface.polyData->GetNumberOfPoints() == static_cast<vtkIdType>(serial.size()),
                      "Parallel update failed");
                for(size_t k = 0; k < serial.size(); ++k)
                {
                    double xyz[3];
                    surface.polyData->GetPoint(k, xyz);
                    check(V(xyz[0], xyz[1], xyz[2]).isApprox(serial[k], 1.e-12),
                          "Serial and parallel displacement differ");
                }
                const double savedTime = base.simulationParameters.totalTime;
                base.simulationParameters.totalTime = std::numeric_limits<double>::quiet_NaN();
                surface.updateConfiguration();
                check(!surface.actor->GetVisibility() && surface.status->text().contains("Nonfinite displacement"),
                      "Worker failure was not reported safely");
                base.simulationParameters.totalTime = savedTime;
                surface.updateConfiguration();
                check(surface.actor->GetVisibility(), "Surface did not recover after worker failure");
                surface.meshSizeEdit->setText("2500");
                surface.rebuild();
                auto *faceCheck = surface.findChild<QCheckBox*>(QString("surfaceExternalFace_%1").arg(firstFaceID));
                check(faceCheck != nullptr, "Missing external face selection");
                check(faceCheck->isChecked() && surface.externalFaces->isChecked() && !surface.internalFaces->isChecked(),
                      "Incorrect default face selections");
                faceCheck->setChecked(false);
                check(surface.polyData->GetNumberOfCells() == static_cast<vtkIdType>(excluded.triangles.size()),
                      "Face checkbox did not rebuild the surface");
                surface.updateConfiguration();
                check(surface.polyData->GetNumberOfCells() == static_cast<vtkIdType>(excluded.triangles.size()),
                      "Frame update lost face exclusions");
                for(const auto id : allExternalIDs)
                {
                    surface.findChild<QCheckBox*>(QString("surfaceExternalFace_%1").arg(id))->setChecked(false);
                }
                check(surface.polyData->GetNumberOfCells() == 0, "Selecting all faces left stale geometry");
                surface.externalFaces->setChecked(false);
                surface.externalFaces->setChecked(true);
                check(surface.polyData->GetNumberOfCells()==0,"External group lost the empty selection");
                faceCheck->setChecked(true);
                surface.externalFaces->setChecked(false);
                check(surface.polyData->GetNumberOfCells()==0,"Disabled external group still renders faces");
                surface.externalFaces->setChecked(true);
                check(faceCheck->isChecked() && surface.polyData->GetNumberOfCells()==static_cast<vtkIdType>(firstFaceSimplices.size()),
                      "External group did not restore the individual selection");
#ifdef _OPENMP
                omp_set_num_threads(savedThreads);
                omp_set_dynamic(savedDynamic);
#endif
            }
            check(renderer->GetActors()->GetNumberOfItems() == 0, "Surface actor leaked after destruction");
            check(renderer->GetLights()->GetNumberOfItems() == originalLights,
                  "Surface relief light leaked after destruction");
        }
        std::filesystem::remove_all(dir);
        std::cout
            << "PASS: surface mesh topology, refinement, displacement steps, frame updates, controls and lifecycle\n";
        return 0;
    }
    catch(const std::exception &error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
