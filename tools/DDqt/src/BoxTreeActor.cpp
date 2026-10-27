#include <BoxTreeActor.h>
#include <GlidePlane.h>
#include <QCheckBox>
#include <QComboBox>
#include <DDqtControls.h>
#include <QSpinBox>
#include <QSlider>
#include <QDoubleSpinBox>
#include <QSignalBlocker>
#include <vtkActor.h>
#include <vtkCellArray.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <algorithm>
#include <map>

namespace model
{
    BoxTreeActor::BoxTreeActor(vtkGenericOpenGLRenderWindow* window, vtkRenderer* ren,
                              std::shared_ptr<BoxTree<3>> boxes, PlaneProvider provider)
        : planeProvider(std::move(provider)),
          cutPolydata(vtkSmartPointer<vtkPolyData>::New()),
          cutMapper(vtkSmartPointer<vtkPolyDataMapper>::New()), cutActor(vtkSmartPointer<vtkActor>::New()),
          renderWindow(window), renderer(ren), tree(std::move(boxes)),
          polydata(vtkSmartPointer<vtkPolyData>::New()),
          mapper(vtkSmartPointer<vtkPolyDataMapper>::New()), actor(vtkSmartPointer<vtkActor>::New())
    {
        auto* layout = new QVBoxLayout(this);
        group = controls::box("Show boxes", this,
            "Display the shared box hierarchy used for loop-FMM geometry. No field kernels are evaluated.");
        group->setObjectName("boxTreeGroup");
        group->setCheckable(true);
        group->setChecked(false);
        auto* grid = controls::grid(group);
        depthEdit = new QLabel(group);
        depthEdit->setObjectName("boxTreeDepth");
        depthEdit->setText(tree ? QString::number(tree->depth()) : QString("—"));
        controls::field(grid, 0, "Tree depth", depthEdit,
            "Uniform subdivisions of the mesh bounding box: level 0 is one box; each level splits every box into eight. "
            "Configured by fmmTreeDepth in inputFiles/DD.txt. Change the input and reload the simulation to change the tree depth.");
        levelEdit = new QSpinBox(group);
        levelEdit->setObjectName("boxTreeLevel");
        levelEdit->setRange(0, tree ? tree->depth() : 0);
        levelEdit->setValue(tree ? tree->depth() : 0);
        controls::field(grid, 1, "Display level", levelEdit,
            "Show one existing level of the hierarchy. This changes visibility only and does not rebuild boxes.");
        opacityEdit = new QSlider(Qt::Horizontal, group);
        opacityEdit->setObjectName("boxTreeOpacity");
        opacityEdit->setRange(0, 100);
        opacityEdit->setValue(60);
        controls::field(grid, 2, "Opacity", opacityEdit, "Opacity of the box edges: 0 is transparent, 100 is opaque.");
        widthEdit = new QDoubleSpinBox(group);
        widthEdit->setObjectName("boxTreeLineWidth");
        widthEdit->setRange(0.5, 5.0);
        widthEdit->setSingleStep(0.5);
        widthEdit->setValue(1.0);
        controls::field(grid, 3, "Line width", widthEdit, "Width of box edges in screen pixels.");
        cutsEdit = new QCheckBox("Glide-plane cuts", group);
        cutsEdit->setObjectName("boxTreeCuts");
        cutsEdit->setToolTip("Display the infinite glide plane cut into closed polygons by boxes at the display level. These polygons are not clipped to grains or loops.");
        grid->addWidget(cutsEdit, 4, 0, 1, 2);
        planeEdit = new QComboBox(group);
        planeEdit->setObjectName("boxTreePlane");
        controls::field(grid, 5, "Glide plane", planeEdit,
            "Select an existing glide plane from the current configuration. Cut geometry is cached on that plane and shared by its loops. Orange edges are shared between neighboring polygons.");
        summary = new QLabel(group);
        summary->setObjectName("boxTreeSummary");
        summary->setWordWrap(true);
        summary->setToolTip("Counts for the displayed level. Shared faces and edges are represented once.");
        grid->addWidget(summary, 6, 0, 1, 2);
        layout->addWidget(group);
        status = new QLabel(this);
        status->setObjectName("boxTreeStatus");
        status->setWordWrap(true);
        layout->addWidget(status);
        layout->addStretch(1);

        mapper->SetInputData(polydata);
        mapper->ScalarVisibilityOff();
        actor->SetMapper(mapper);
        actor->GetProperty()->SetColor(0.25, 0.35, 0.65);
        actor->GetProperty()->LightingOff();
        actor->SetVisibility(false);
        actor->PickableOff();
        renderer->AddActor(actor);
        cutMapper->SetInputData(cutPolydata);
        cutMapper->ScalarVisibilityOff();
        cutActor->SetMapper(cutMapper);
        cutActor->GetProperty()->SetColor(1.0, 0.45, 0.1);
        cutActor->GetProperty()->SetOpacity(0.35);
        cutActor->GetProperty()->SetLineWidth(2.0);
        cutActor->GetProperty()->LightingOff();
        cutActor->SetVisibility(false);
        cutActor->PickableOff();
        renderer->AddActor(cutActor);
        if(!tree)
        {
            group->setEnabled(false);
            status->setText("Box topology is available when DislocationDynamics is enabled.");
        }
        else
        {
            status->setText("Boxes enclose the mesh; they are not clipped to grains.");
        }
        connect(group, &QGroupBox::toggled, this, [this]() { modify(); });
        connect(cutsEdit, &QCheckBox::toggled, this, [this]() { modify(); });
        connect(planeEdit, &QComboBox::currentIndexChanged, this, [this]() { modify(); });
        connect(levelEdit, &QSpinBox::valueChanged, this, [this]() { modify(); });
        connect(opacityEdit, &QSlider::valueChanged, this, [this]() { modify(); });
        connect(widthEdit, &QDoubleSpinBox::valueChanged, this, [this]() { modify(); });
    }

    BoxTreeActor::~BoxTreeActor()
    {
        renderer->RemoveActor(actor);
        renderer->RemoveActor(cutActor);
    }

    void BoxTreeActor::refreshGeometry()
    {
        if(!tree)
        {
            return;
        }
        const unsigned level = std::min(static_cast<unsigned>(levelEdit->value()), tree->depth());
        if(displayedRevision == tree->revision() && displayedLevel == static_cast<int>(level))
        {
            return;
        }
        const QSignalBlocker levelBlocker(levelEdit);
        depthEdit->setText(QString::number(tree->depth()));
        levelEdit->setRange(0, tree->depth());
        levelEdit->setValue(level);
        auto points = vtkSmartPointer<vtkPoints>::New();
        points->SetDataTypeToDouble();
        auto lines = vtkSmartPointer<vtkCellArray>::New();
        std::map<BoxEntity<3>::Key, vtkIdType> vertices;
        const auto edges = tree->edges(level);
        for(const auto& edge : edges)
        {
            vtkIdType ids[2];
            for(size_t side = 0; side < 2; ++side)
            {
                const auto& vertex = edge->boundary().at(side).entity;
                auto found = vertices.find(vertex->key());
                if(found == vertices.end())
                {
                    found = vertices.emplace(vertex->key(), points->InsertNextPoint(vertex->lower().data())).first;
                }
                ids[side] = found->second;
            }
            lines->InsertNextCell(2, ids);
        }
        polydata->SetPoints(points);
        polydata->SetLines(lines);
        polydata->Modified();
        displayedRevision = tree->revision();
        displayedLevel = level;
        const unsigned n = 1u << level;
        const Eigen::Vector3d size = (tree->upper() - tree->lower()) / n;
        summary->setText(QString("%1 boxes, %2 faces, %3 edges, %4 vertices\nBox size: %5 × %6 × %7 b")
            .arg(static_cast<qulonglong>(tree->boxes(level).size())).arg(3u * n * n * (n + 1))
            .arg(static_cast<qulonglong>(edges.size())).arg(static_cast<qulonglong>(vertices.size()))
            .arg(size[0], 0, 'g', 5).arg(size[1], 0, 'g', 5).arg(size[2], 0, 'g', 5));
    }

    void BoxTreeActor::modify()
    {
        if(group->isChecked())
        {
            refreshGeometry();
            refreshCuts();
        }
        actor->SetVisibility(tree && group->isChecked());
        cutActor->SetVisibility(tree && group->isChecked() && cutsEdit->isChecked() && bool(displayedCut));
        actor->GetProperty()->SetOpacity(opacityEdit->value() / 100.0);
        actor->GetProperty()->SetLineWidth(widthEdit->value());
        renderWindow->Render();
    }

    void BoxTreeActor::updateConfiguration()
    {
        // Fixed geometry is not rebuilt on frame changes. Future users can
        // replace the hierarchy via setDepth(), detected through revision().
        if(group->isChecked())
        {
            refreshGeometry();
            refreshCuts();
        }
    }

    void BoxTreeActor::refreshCuts()
    {
        const auto selected = planeEdit->currentData();
        const QSignalBlocker blocker(planeEdit);
        availablePlanes.clear();
        planeEdit->clear();
        if(planeProvider)
        {
            for(const auto& plane : planeProvider())
            {
                if(plane)
                {
                    availablePlanes.emplace(plane->sID, plane);
                    planeEdit->addItem(QString("%1 (grain %2)").arg(plane->sID).arg(plane->regionIDs.first),
                                       QVariant::fromValue<qulonglong>(plane->sID));
                }
            }
        }
        const int oldIndex = planeEdit->findData(selected);
        if(oldIndex >= 0)
        {
            planeEdit->setCurrentIndex(oldIndex);
        }
        planeEdit->setEnabled(cutsEdit->isChecked() && !availablePlanes.empty());
        std::shared_ptr<const PlaneBoxIntersection<3>> cut;
        QString cutError;
        if(cutsEdit->isChecked() && tree && planeEdit->currentIndex() >= 0)
        {
            const auto plane = availablePlanes.at(planeEdit->currentData().toULongLong()).lock();
            if(plane)
            {
                try
                {
                    cut = plane->boxIntersection(tree, levelEdit->value());
                }
                catch(const std::exception& error)
                {
                    cutError = QString("Could not compute plane cuts: %1").arg(error.what());
                }
            }
        }
        if(cut != displayedCut)
        {
            auto points = vtkSmartPointer<vtkPoints>::New();
            points->SetDataTypeToDouble();
            auto polygons = vtkSmartPointer<vtkCellArray>::New();
            auto lines = vtkSmartPointer<vtkCellArray>::New();
            std::map<PlaneBoxIntersection<3>::Key, vtkIdType> pointIDs;
            std::map<PlaneBoxIntersection<3>::SegmentKey, PlaneBoxIntersection<3>::SegmentPtr> edges;
            if(cut)
            {
                for(const auto& [key, polygon] : cut->polygons())
                {
                    std::vector<vtkIdType> ids;
                    for(const auto& point : polygon.vertices)
                    {
                        auto found = pointIDs.find(point->support->key());
                        if(found == pointIDs.end())
                        {
                            found = pointIDs.emplace(point->support->key(), points->InsertNextPoint(point->position.data())).first;
                        }
                        ids.push_back(found->second);
                    }
                    polygons->InsertNextCell(ids.size(), ids.data());
                    for(const auto& side : polygon.boundary)
                    {
                        const auto& edge = side.segment;
                        edges.emplace(std::make_pair(edge->first->support->key(), edge->second->support->key()), edge);
                    }
                }
                for(const auto& [key, edge] : edges)
                {
                    const vtkIdType ids[] = {pointIDs.at(key.first), pointIDs.at(key.second)};
                    lines->InsertNextCell(2, ids);
                }
            }
            cutPolydata->SetPoints(points);
            cutPolydata->SetPolys(polygons);
            cutPolydata->SetLines(lines);
            cutPolydata->Modified();
            displayedCut = cut;
        }
        cutActor->SetVisibility(group->isChecked() && cutsEdit->isChecked() && bool(displayedCut));
        if(!cutError.isEmpty())
        {
            status->setText(cutError);
        }
        else if(cutsEdit->isChecked())
        {
            status->setText(displayedCut ? QString("%1 plane-cut polygons; not clipped to grains or loops.").arg(displayedCut->polygons().size())
                                        : QString("No glide planes in the current configuration."));
        }
        else
        {
            status->setText("Boxes enclose the mesh; they are not clipped to grains.");
        }
    }

    vtkPolyData* BoxTreeActor::cutData() const
    {
        return cutPolydata;
    }

    vtkPolyData* BoxTreeActor::data() const
    {
        return polydata;
    }
}
