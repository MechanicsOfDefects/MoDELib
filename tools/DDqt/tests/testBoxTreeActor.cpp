#include <BoxTreeActor.h>
#include <DefectiveCrystal/DefectiveCrystalBase.h>
#include <DislocationDynamics/DislocationDynamicsBase.h>
#include <GlidePlane.h>
#include <QCheckBox>
#include <QComboBox>
#include <filesystem>
#include <fstream>
#include <chrono>
#include <QApplication>
#include <QDoubleSpinBox>
#include <QGroupBox>
#include <QLabel>
#include <QSpinBox>
#include <QSlider>
#include <vtkActor.h>
#include <vtkActorCollection.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkPolyData.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <iostream>
#include <stdexcept>

void require(bool condition, const char* message)
{
    if(!condition)
    {
        throw std::runtime_error(message);
    }
}
int main(int argc, char** argv)
{
    QApplication app(argc, argv);
    try
    {
        auto window = vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
        window->SetReadyForRendering(false);
        auto renderer = vtkSmartPointer<vtkRenderer>::New();
        window->AddRenderer(renderer);
        auto tree = std::make_shared<model::BoxTree<3>>(Eigen::Vector3d(-1, 2, 3), Eigen::Vector3d(3, 6, 11));
        {
            model::BoxTreeActor actor(window, renderer, tree);
            auto* group = actor.findChild<QGroupBox*>("boxTreeGroup");
            auto* depth = actor.findChild<QLabel*>("boxTreeDepth");
            auto* level = actor.findChild<QSpinBox*>("boxTreeLevel");
            auto* opacity = actor.findChild<QSlider*>("boxTreeOpacity");
            auto* width = actor.findChild<QDoubleSpinBox*>("boxTreeLineWidth");
            require(group && depth && level && opacity && width, "Missing box controls");
            auto* actors = renderer->GetActors();
            actors->InitTraversal();
            auto* prop = actors->GetNextActor();
            require(prop, "Missing box wireframe actor");
            require(!group->isChecked() && !prop->GetVisibility(), "Boxes should start hidden");
            group->setChecked(true);
            require(actor.data()->GetNumberOfLines() == 12 && actor.data()->GetNumberOfPoints() == 8,
                    "Root wireframe is incorrect");
            require(!actor.findChild<QSpinBox*>("boxTreeDepth") && depth->text() == "0",
                    "Tree depth must be displayed without an editable control");
            tree->setDepth(2);
            actor.updateConfiguration();
            level->setValue(2);
            require(depth->text() == "2" && level->maximum() == 2,
                    "Viewer did not synchronize the existing hierarchy");
            require(actor.data()->GetNumberOfLines() == 300 && actor.data()->GetNumberOfPoints() == 125,
                    "Level 2 wireframe has missing or duplicate edges");
            const auto root = tree->root();
            const auto revision = tree->revision();
            const auto modified = actor.data()->GetMTime();
            actor.updateConfiguration();
            require(tree->root() == root && actor.data()->GetMTime() == modified, "Frame update rebuilt static boxes");
            level->setValue(1);
            require(actor.data()->GetNumberOfLines() == 54 && tree->revision() == revision,
                    "Display selection rebuilt tree or produced incorrect edges");
            double bounds[6];
            actor.data()->GetBounds(bounds);
            require(bounds[0] == -1 && bounds[1] == 3 && bounds[2] == 2 && bounds[3] == 6 &&
                    bounds[4] == 3 && bounds[5] == 11, "Box wireframe is positioned incorrectly");
            opacity->setValue(25);
            width->setValue(2.5);
            require(prop->GetProperty()->GetOpacity() == .25 && prop->GetProperty()->GetLineWidth() == 2.5,
                    "Style controls do not reach the actor");
            group->setChecked(false);
            require(!prop->GetVisibility(), "Disabling boxes leaves actor visible");
            group->setChecked(true);
            require(level->value() == 1 && opacity->value() == 25 && tree->revision() == revision,
                    "Toggling visibility lost selections or rebuilt boxes");
            tree->setDepth(0);
            actor.updateConfiguration();
            require(actor.data()->GetNumberOfLines() == 12 && level->value() == 0 && depth->text() == "0",
                    "Actor did not synchronize an external topology revision");
        }
        require(renderer->GetActors()->GetNumberOfItems() == 0, "Box actor leaked into renderer");
        tree->setDepth(2);
        {
            model::BoxTreeActor actor(window, renderer, tree);
            require(actor.findChild<QLabel*>("boxTreeDepth")->text() == "2" &&
                    actor.findChild<QSpinBox*>("boxTreeLevel")->value() == 2,
                    "Viewer did not start at the configured tree depth");
            actor.findChild<QGroupBox*>("boxTreeGroup")->setChecked(true);
            require(actor.data()->GetNumberOfLines() == 300,
                    "Configured tree depth was not visualized");
        }
        {
            model::BoxTreeActor actor(window, renderer, nullptr);
            auto* group = actor.findChild<QGroupBox*>("boxTreeGroup");
            require(!group->isEnabled(), "Box controls enabled without DD geometry");
            actor.updateConfiguration();
        }
        require(argc == 2, "Expected repository root for plane-cache viewer test");
        const auto folder = std::filesystem::temp_directory_path() /
            ("box-plane-viewer-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        const auto input = folder / "inputFiles";
        std::filesystem::create_directories(input);
        std::filesystem::create_directories(folder / "F");
        const std::filesystem::path root(argv[1]);
        std::filesystem::copy_file(root / "Library/Materials/W.txt", input / "W.txt");
        std::filesystem::copy_file(root / "Library/Meshes/unitCube24.msh", input / "mesh.msh");
        std::ofstream(input / "DefectiveCrystal.txt")
            << "physics=DislocationDynamics;\nfmmTreeDepth=2;\nuseFEM=0;\nNsteps=0;\nmaxResolveSteps=0;\ndtMax=1;\n"
            << "outputFrequency=1;\noutputBinary=0;\nstartAtTimeStep=0;\nperiodicImageSize=0 0 0;\nEwaldLengthFactor=0;\n";
        std::ofstream(input / "polycrystal.txt")
            << "materialFile=W.txt;\nmeshFile=mesh.msh;\nabsoluteTemperature=300;\nC2G1=1 0 0 0 1 0 0 0 1;\n"
            << "F=100 0 0 0 100 0 0 0 100;\nX0=0 0 0;\nperiodicFaceIDs=;\n";
        {
            model::DefectiveCrystalBase<3> base(folder.string());
            auto& dd = base.dislocationDynamics();
            const auto grain = base.poly.grain(1);
            auto plane = dd.glidePlaneFactory.getFromKey(model::GlidePlaneKey<3>(
                (base.mesh.xMin() + base.mesh.xMax()) * .5,
                grain->reciprocalLatticeDirection(Eigen::Vector3d::UnitZ())));
            std::vector<std::shared_ptr<model::GlidePlane<3>>> planes{plane};
            model::BoxTreeActor actor(window, renderer, base.fmmContext->boxTree, [&planes]() { return planes; });
            actor.findChild<QGroupBox*>("boxTreeGroup")->setChecked(true);
            auto* cuts = actor.findChild<QCheckBox*>("boxTreeCuts");
            auto* selector = actor.findChild<QComboBox*>("boxTreePlane");
            auto* level = actor.findChild<QSpinBox*>("boxTreeLevel");
            require(selector->count() == 1 && actor.cutData()->GetNumberOfPolys() == 0,
                    "Plane cuts should be opt-in");
            cuts->setChecked(true);
            const auto expected = plane->boxIntersection(base.fmmContext->boxTree, 2);
            require(actor.cutData()->GetNumberOfPolys() > 0 &&
                    size_t(actor.cutData()->GetNumberOfPolys()) == expected->polygons().size(),
                    "Viewer did not display the cached plane polygons");
            const auto modified = actor.cutData()->GetMTime();
            actor.updateConfiguration();
            require(actor.cutData()->GetMTime() == modified, "Frame update rebuilt unchanged plane cuts");
            const auto revision = base.fmmContext->boxTree->revision();
            level->setValue(0);
            require(actor.cutData()->GetNumberOfPolys() == 1 && base.fmmContext->boxTree->revision() == revision,
                    "Display level did not select existing cut level");
            planes.clear();
            plane.reset();
            actor.updateConfiguration();
            require(selector->count() == 0 && actor.cutData()->GetNumberOfPolys() == 0,
                    "Disappearing glide plane left stale cut polygons");
        }
        std::filesystem::remove_all(folder);
        require(renderer->GetActors()->GetNumberOfItems() == 0, "Plane-cut actor leaked into renderer");
        std::cout << "Box viewer controls, plane cuts, geometry, revision handling and cleanup passed\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
