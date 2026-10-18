#include <DefectiveCrystalActor.h>
#include <DDqtMainWindow.h>
#include <QApplication>
#include <QPointer>
#include <QEventLoop>
#include <vtkRenderer.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor2D.h>
#include <vtkActorCollection.h>
#include <vtkDataSetMapper.h>
#include <vtkCallbackCommand.h>
#include <vtkExecutive.h>
#include <vtkProperty.h>
#include <filesystem>
#include <fstream>
#include <chrono>
#include <iostream>
#include <atomic>
#include <limits>

void check(bool condition, const char *message)
{
    if(!condition)
    {
        throw std::runtime_error(message);
    }
}

void checkDataSetMapperInputs(vtkRenderer* renderer)
{
    size_t checked=0;
    auto* actors=renderer->GetActors();
    actors->InitTraversal();
    while(auto* actor=actors->GetNextActor())
    {
        auto* mapper=vtkDataSetMapper::SafeDownCast(actor->GetMapper());
        if(mapper)
        {
            ++checked;
            check(mapper->GetNumberOfInputConnections(0)==1,
                  "Dataset mapper has no input when inclusion physics is absent");
            bool error=false;
            auto callback=vtkSmartPointer<vtkCallbackCommand>::New();
            callback->SetClientData(&error);
            callback->SetCallback([](vtkObject*,unsigned long,void* data,void*)
            {
                *static_cast<bool*>(data)=true;
            });
            const auto observer=mapper->GetExecutive()->AddObserver(vtkCommand::ErrorEvent,callback);
            mapper->Update();
            mapper->GetExecutive()->RemoveObserver(observer);
            check(!error,"Dataset mapper reported a pipeline error");
        }
    }
    check(checked>0,"Expected at least the polyhedral inclusion dataset mapper");
}

void checkPeriodicImages()
{
    auto window = vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
    window->SetReadyForRendering(false);
    auto renderer = vtkSmartPointer<vtkRenderer>::New();
    window->AddRenderer(renderer);
    auto source = vtkSmartPointer<vtkActor>::New();
    auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    auto data = vtkSmartPointer<vtkPolyData>::New();
    mapper->SetInputData(data);
    source->SetMapper(mapper);
    source->SetPosition(1, 2, 3);
    source->RotateZ(30);
    source->SetScale(2, 3, 4);
    renderer->AddActor(source);
    auto overlay = vtkSmartPointer<vtkActor2D>::New();
    renderer->AddActor2D(overlay);
    {
        model::PeriodicImages images(renderer);
        images.setTranslations({Eigen::Vector3d(10, 20, 30)});
        check(renderer->GetViewProps()->GetNumberOfItems() == 3, "Incorrect periodic copy count");
        auto* copy = vtkActor::SafeDownCast(renderer->GetViewProps()->GetLastProp());
        check(copy && copy != source && copy->GetMapper() == mapper, "Copy does not share source geometry");
        for(int row = 0; row < 4; ++row)
        {
            for(int col = 0; col < 4; ++col)
            {
                const double shift = col == 3 && row < 3 ? 10.0 * (row + 1) : 0.0;
                check(std::abs(copy->GetMatrix()->GetElement(row, col) -
                      source->GetMatrix()->GetElement(row, col) - shift) < 1.e-12,
                      "Periodic translation corrupted source transformation");
            }
        }
        source->VisibilityOff();
        source->GetProperty()->SetOpacity(0.25);
        auto replacement = vtkSmartPointer<vtkPolyDataMapper>::New();
        replacement->SetInputData(data);
        source->SetMapper(replacement);
        renderer->InvokeEvent(vtkCommand::StartEvent);
        check(!copy->GetVisibility() && copy->GetProperty()->GetOpacity() == 0.25 &&
              copy->GetMapper() == replacement, "Render callback did not synchronize actor settings");
        check(renderer->GetViewProps()->GetNumberOfItems() == 3, "Copies recursively replicated");
        renderer->RemoveActor(source);
        renderer->InvokeEvent(vtkCommand::StartEvent);
        check(renderer->GetViewProps()->GetNumberOfItems() == 1, "Removed source retained its copies");
        renderer->AddActor(source);
        renderer->InvokeEvent(vtkCommand::StartEvent);
        check(renderer->GetViewProps()->GetNumberOfItems() == 3, "New actor was not replicated");
    }
    check(renderer->GetViewProps()->GetNumberOfItems() == 2, "Periodic images leaked on destruction");
    Eigen::Matrix<double, 3, Eigen::Dynamic> basis(3, 2);
    basis << 10, 3, 0, 20, 0, 4;
    {
        model::PeriodicImagesWidget controls(window, renderer, basis, nullptr);
        auto* group = controls.findChild<QGroupBox*>("periodicImagesGroup");
        auto* first = controls.findChild<QSpinBox*>("periodicImageCount1");
        auto* second = controls.findChild<QSpinBox*>("periodicImageCount2");
        check(group && first && second && !group->isChecked(), "Missing periodic controls or wrong default");
        first->setValue(1);
        second->setValue(1);
        group->setChecked(true);
        check(renderer->GetViewProps()->GetNumberOfItems() == 10, "Two directions did not create nine cells");
        std::set<std::array<int, 3>> actual;
        vtkCollectionSimpleIterator iterator;
        renderer->GetViewProps()->InitTraversal(iterator);
        while(auto* prop = renderer->GetViewProps()->GetNextProp(iterator))
        {
            if(auto* actor = vtkActor::SafeDownCast(prop))
            {
                std::array<int, 3> shift;
                for(int row = 0; row < 3; ++row)
                {
                    shift[row] = std::lround(actor->GetMatrix()->GetElement(row, 3) - source->GetMatrix()->GetElement(row, 3));
                }
                actual.insert(shift);
            }
        }
        for(int i = -1; i <= 1; ++i)
        {
            for(int j = -1; j <= 1; ++j)
            {
                check(actual.count({10*i + 3*j, 20*j, 4*j}) == 1, "Missing oblique lattice image");
            }
        }
        group->setChecked(false);
        check(renderer->GetViewProps()->GetNumberOfItems() == 2, "Disabling images did not restore original scene");
    }
}

struct FailingElastic : model::ElasticDeformation<3>
{
    using model::ElasticDeformation<3>::ElasticDeformation;
    mutable std::atomic<bool> failNext{false};
    VectorDim displacement(const VectorDim& point, const NodeType* node, const ElementType* element,
                           const SimplexDim* guess) const override
    {
        if(failNext.exchange(false))
        {
            throw std::runtime_error("Injected field failure");
        }
        return model::ElasticDeformation<3>::displacement(point, node, element, guess);
    }
};

int main(int argc, char **argv)
{
    QApplication app(argc, argv);
    try
    {
        checkPeriodicImages();
        check(argc == 2, "Expected repository root");
        const std::filesystem::path root(argv[1]);
        const auto dir =
            std::filesystem::temp_directory_path() /
            ("model-surface-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        for(const auto *name : {"inputFiles", "evl", "F"})
        {
            std::filesystem::create_directories(dir / name);
        }
        const auto input = dir / "inputFiles";
        std::filesystem::copy_file(root / "Library/Materials/W.txt", input / "W.txt");
        std::filesystem::copy_file(root / "Library/Meshes/unitCube24.msh", input / "mesh.msh");
        {
            std::ofstream f(input / "DefectiveCrystal.txt");
            f << "physics=ElasticDeformation;\nuseFEM=0;\nNsteps=1;\nmaxResolveSteps=0;\ndtMax=1;\n"
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
        // Reopening a visualization-only tutorial leaves empty history files.
        // Explicit startAtTimeStep=0 must not attempt to dereference a missing row.
        std::ofstream(dir / "F/F_0.txt").close();
        std::ofstream(dir / "F/F_labels.txt").close();
        {
            model::DislocationDynamicsBase<3> base(dir.string());
            check(base.simulationParameters.runID == 0 && base.simulationParameters.totalTime == 0.0,
                  "Empty history did not retain the default restart state");

            model::DDconfigIO<3> config((dir / "evl").string());
            model::DDauxIO<3> aux((dir / "evl").string());
            for(size_t frame : {0, 10, 20})
            {
                config.write(frame, false);
                aux.write(frame, false);
            }
            auto window = vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
            window->SetReadyForRendering(false);
            auto renderer = vtkSmartPointer<vtkRenderer>::New();
            window->AddRenderer(renderer);
            auto *viewer = new model::DefectiveCrystalActor(window, renderer, nullptr, base);
            QPointer<model::DefectiveCrystalActor> guardedViewer(viewer);
            checkDataSetMapperInputs(renderer);
            check(viewer->updateConfiguration(0), "Initial frame load failed");
            checkDataSetMapperInputs(renderer);
            // Add a sampling plane so each successful load must refresh Fields.
            viewer->ddField->findChildren<QSpinBox*>().front()->setValue(1);
            auto *plane = viewer->ddField->findChild<model::DDPlaneField*>();
            check(plane && !plane->dataPnts().empty(), "Missing field sampling plane");
            int notifications = 0;
            QObject::connect(viewer, &model::DefectiveCrystalActor::frameChanged, viewer, [&](qulonglong frame)
            {
                ++notifications;
                check(viewer->frameIDedit->text().toULongLong() == frame, "Frame ID not committed before notification");
                check(plane->dataPnts().front().displacement.front().allFinite(), "Fields not ready at frame notification");
            });
            auto poisonField = [&]()
            {
                plane->dataPnts().front().displacement.front().setConstant(std::numeric_limits<double>::quiet_NaN());
            };
            viewer->frameIncrementEdit->setText("10");
            poisonField();
            viewer->plusFrameButton->click();
            check(notifications == 1 && plane->dataPnts().front().displacement.front().allFinite(), "Next did not update Fields");
            poisonField();
            viewer->frameIDedit->setText("0");
            viewer->frameIDedit->setModified(true);
            QMetaObject::invokeMethod(viewer->frameIDedit, "editingFinished", Qt::DirectConnection);
            check(notifications == 2 && plane->dataPnts().front().displacement.front().allFinite(), "Frame edit did not update Fields");
            poisonField();
            viewer->stopFrameButton->click();
            QMetaObject::invokeMethod(viewer, "stopConfigurations", Qt::DirectConnection);
            check(!plane->dataPnts().front().displacement.front().allFinite(), "Stop recomputed Fields");
            QMetaObject::invokeMethod(viewer, "playConfigurations", Qt::DirectConnection);
            QMetaObject::invokeMethod(viewer, "playConfigurations", Qt::DirectConnection);
            auto timers = viewer->findChildren<QTimer*>("framePlaybackTimer");
            check(timers.size() == 1 && timers.front()->isActive(), "Play created duplicate timers");
            QMetaObject::invokeMethod(timers.front(), "timeout", Qt::DirectConnection);
            check(viewer->frameIDedit->text() == "10" && plane->dataPnts().front().displacement.front().allFinite(), "Playback did not update Fields");
            checkDataSetMapperInputs(renderer);
            viewer->stopConfigurations();
            check(!timers.front()->isActive(), "Stop did not stop the timer immediately");
            const int committed = notifications;
            const auto pointsTime = viewer->nodes->nodePolyData->GetMTime();
            check(!viewer->updateConfiguration(999), "Missing frame unexpectedly loaded");
            check(notifications == committed && viewer->frameIDedit->text() == "10" &&
                  viewer->nodes->nodePolyData->GetMTime() == pointsTime, "Failed load changed the scene or frame ID");
            config.write(30, false); // Missing auxiliary data must fail before touching actors.
            check(!viewer->updateConfiguration(30) && viewer->nodes->nodePolyData->GetMTime() == pointsTime,
                  "Missing auxiliary file partially updated the scene");
            auto failingPhysics = std::make_shared<FailingElastic>(viewer->defectiveCrystal);
            viewer->defectiveCrystal.microstructures().front() = failingPhysics;
            const auto *previousConfig = &viewer->configIO();
            const auto previousCamera = Eigen::Vector3d(renderer->GetActiveCamera()->GetPosition());
            failingPhysics->failNext = true;
            check(!viewer->updateConfiguration(20), "Field failure did not abort the frame");
            check(&viewer->configIO() == previousConfig && viewer->frameIDedit->text() == "10" &&
                  notifications == committed, "Field failure did not roll back configuration and frame ID");
            check(plane->dataPnts().front().displacement.front().allFinite(), "Rollback did not restore field samples");
            check(Eigen::Vector3d(renderer->GetActiveCamera()->GetPosition()).isApprox(previousCamera),
                  "Failed frame changed the camera");
            viewer->saveImage->setChecked(true);
            check(viewer->updateConfiguration(20) && viewer->frameStatus->text().contains("capture failed"),
                  "Capture failure incorrectly undid a valid frame");
            check(plane->dataPnts().front().displacement.front().allFinite(), "Capture preceded field refresh");
            viewer->saveImage->setChecked(false);
            const Eigen::Vector3d center=0.5*(base.mesh.xMin()+base.mesh.xMax());
            auto* camera=renderer->GetActiveCamera();
            const Eigen::Vector3d initialPosition=center+Eigen::Vector3d(10,0,0);
            camera->SetPosition(initialPosition.data());
            camera->SetFocalPoint(center.data());
            camera->SetViewUp(0,0,1);
            viewer->cameraAdvanceEdit->setText("1");
            viewer->rotationAngleEdit->setText("0");
            auto manualFrame=[&](const char* text)
            {
                viewer->frameIDedit->setText(text);
                viewer->frameIDedit->setModified(true);
                QMetaObject::invokeMethod(viewer->frameIDedit,"editingFinished",Qt::DirectConnection);
            };
            manualFrame("10");
            check(Eigen::Vector3d(camera->GetPosition()).isApprox(center),
                  "Manual frame selection did not advance camera to center");
            manualFrame("20");
            check(Eigen::Vector3d(camera->GetPosition()).isApprox(center-Eigen::Vector3d(10,0,0)),
                  "Manual frame selection did not advance camera past center");
            const int beforeRepeatedFinish=notifications;
            QMetaObject::invokeMethod(viewer->frameIDedit,"editingFinished",Qt::DirectConnection);
            check(notifications==beforeRepeatedFinish,"Focus loss reloaded an already committed frame");
            manualFrame("0");
            check(Eigen::Vector3d(camera->GetPosition()).isApprox(initialPosition),
                  "Manual backward selection did not restore camera position");
            // The main window must actually destroy removed pages and their timers.
            model::DDqtMainWindow mainWindow;
            auto *tabs = qobject_cast<QTabWidget*>(mainWindow.centralWidget());
            const int index = tabs->addTab(viewer, "Test viewer");
            QMetaObject::invokeMethod(viewer, "playConfigurations", Qt::DirectConnection);
            QPointer<QTimer> timer(timers.front());
            QMetaObject::invokeMethod(&mainWindow, "closeTab", Qt::DirectConnection, Q_ARG(int, index));
            QCoreApplication::sendPostedEvents(nullptr, QEvent::DeferredDelete);
            check(guardedViewer.isNull() && timer.isNull(), "Closing the tab retained the viewer or timer");
        }
        std::filesystem::remove_all(dir);
        std::cout << "PASS: frame navigation, Fields ordering, playback, failed loads, capture failure and tab cleanup\n";
        return 0;
    }
    catch(const std::exception &error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
