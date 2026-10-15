#include <DefectiveCrystalActor.h>
#include <DDqtMainWindow.h>
#include <QApplication>
#include <QPointer>
#include <QEventLoop>
#include <vtkRenderer.h>
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
        {
            model::DislocationDynamicsBase<3> base(dir.string());

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
            check(viewer->updateConfiguration(0), "Initial frame load failed");
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
            QMetaObject::invokeMethod(viewer->frameIDedit, "returnPressed", Qt::DirectConnection);
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
