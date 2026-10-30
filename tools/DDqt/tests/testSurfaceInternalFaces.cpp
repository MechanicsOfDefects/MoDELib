#include <SurfaceMeshActor.h>
#include <QApplication>
#include <filesystem>
#include <fstream>
#include <chrono>
#include <iostream>
void check(bool condition,const char* message)
{
    if(!condition)
    {
        throw std::runtime_error(message);
    }
}
int main(int argc,char** argv)
{
    QApplication app(argc,argv);
    try
    {
        check(argc == 2, "Expected repository root");
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
        std::filesystem::copy_file(root / "Library/Meshes/bicrystal_48.msh", input / "mesh.msh");
        {
            std::ofstream f(input / "DefectiveCrystal.txt");
            f << "physics=ElasticDeformation;\nfmmTreeDepth=0;\nuseFEM=0;\nNsteps=1;\nmaxResolveSteps=0;\ndtMax=1;\n"
              << "outputFrequency=1;\noutputBinary=0;\nstartAtTimeStep=0;\n";
            std::ofstream p(input / "polycrystal.txt");
            p << "materialFile=W.txt;\nmeshFile=mesh.msh;\nabsoluteTemperature=300;\n"
              << "C2G1=1 0 0 0 1 0 0 0 1;\nC2G2=0.8 0 0.6 0 1 0 -0.6 0 0.8;\nF=100 0 0 0 100 0 0 0 100;\nX0=0 0 0;\nperiodicFaceIDs=;\n";
            std::ofstream e(input / "ElasticDeformation.txt");
            e << "useElasticDeformationFEM=0;\ninertiaReliefPenaltyFactor=1e-3;\n"
              << "ExternalStress0=0 0 .01 0 0 0;\nExternalStressRate=0 0 .01 0 0 0;\n"
              << "ExternalStrain0=0 0 0 0 0 0;\nExternalStrainRate=0 0 0 0 0 0;\n"
              << "stiffnessRatio=0 0 0 0 0 0;\n";
        }

        {
            model::DefectiveCrystalBase<3> base(dir.string());
            std::set<size_t> externalIDs,internalIDs;
            std::set<const model::Simplex<3,2>*> internalTriangles;
            for(const auto& region : base.mesh.regions())
            {
                for(const auto& entry : region.second->faces())
                {
                    const auto& face=*entry.second;
                    if(face.isExternal())
                    {
                        externalIDs.insert(face.sID);
                    }
                    else
                    {
                        internalIDs.insert(face.sID);
                        internalTriangles.insert(face.simplices().begin(),face.simplices().end());
                    }
                }
            }
            check(!internalIDs.empty(),"Bicrystal has no internal faces");
            const auto inside=model::SurfaceMeshGeometry::build(base.mesh,625,internalIDs);
            double area=0;
            for(size_t i=0;i<inside.triangles.size();++i)
            {
                const auto& t=inside.triangles[i];
                const Eigen::Vector3d cross=(inside.referencePoints[t[1]]-inside.referencePoints[t[0]])
                    .cross(inside.referencePoints[t[2]]-inside.referencePoints[t[0]]);
                const auto* face=inside.boundaryFaces[i];
                check(internalTriangles.count(face),"Internal selection included external face");
                check(cross.dot(face->outNormal(model::SurfaceMeshGeometry::samplingParent(face)->region->regionID))>0,"Internal orientation invalid");
                area+=0.5*cross.norm();
            }
            double expected=0;
            for(const auto* face : internalTriangles)
            {
                expected+=face->vol0;
            }
            check(std::abs(area-expected)<1.e-8,"Internal faces duplicated or lost during refinement");
            model::SurfaceMeshGeometry::Trace trace{Eigen::Vector3d::Zero(),Eigen::Vector3d::UnitY(),internalTriangles};
            const auto cut=inside.split({trace});
            check(!cut.walls.empty(),"Internal interface did not acquire slip-step walls");
            std::set<size_t> all=externalIDs;
            all.insert(internalIDs.begin(),internalIDs.end());
            const auto both=model::SurfaceMeshGeometry::build(base.mesh,1.e20,all);
            const auto outside=model::SurfaceMeshGeometry::build(base.mesh,1.e20,externalIDs);
            check(both.triangles.size()==outside.triangles.size()+internalTriangles.size(),"Combined faces duplicated interfaces");
            check(model::SurfaceMeshGeometry::build(base.mesh,625,{}).triangles.empty(),"Empty selection rendered faces");
            model::DefectiveCrystal<3> crystal(base);
            model::DefectiveCrystalConfigIO<3> config((dir/"config").string());
            crystal.initializeConfiguration(config);
            auto window=vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
            window->SetReadyForRendering(false);
            auto renderer=vtkSmartPointer<vtkRenderer>::New();
            window->AddRenderer(renderer);
            model::SurfaceMeshActor surface(window,renderer,crystal);
            check(!surface.externalFaces->isChecked() && !surface.internalFaces->isChecked(),"Incorrect group defaults");
            for(const auto id : externalIDs)
            {
                auto* item=surface.findChild<QCheckBox*>(QString("surfaceExternalFace_%1").arg(id));
                check(item && !item->isChecked(),"External face selected by default");
            }
            for(const auto id : internalIDs)
            {
                auto* item=surface.findChild<QCheckBox*>(QString("surfaceInternalFace_%1").arg(id));
                check(item && !item->isChecked(),"Internal face selected by default");
            }
            check(surface.findChildren<QCheckBox*>(QRegularExpression("surfaceInternalFace_.*")).size()==static_cast<int>(internalIDs.size()),
                  "Internal face appears twice in UI");
            surface.meshSizeEdit->setText("1e20");
            surface.externalFaces->setChecked(false);
            surface.updateConfiguration();
            surface.showSurface->setChecked(true);
            check(surface.polyData->GetNumberOfCells()==0,"Disabled groups rendered selected faces");
            surface.internalFaces->setChecked(true);
            check(surface.polyData->GetNumberOfCells()==static_cast<vtkIdType>(internalTriangles.size()),"Internal face selection did not render");
            for(const auto id : internalIDs)
            {
                auto* item=surface.findChild<QCheckBox*>(QString("surfaceInternalFace_%1").arg(id));
                check(item->isChecked(),"Enabling internal group did not select all faces");
                item->setChecked(false);
            }
            check(surface.polyData->GetNumberOfCells()==0,"Individual deselection did not hide internal faces");
            surface.internalFaces->setChecked(false);
            surface.internalFaces->setChecked(true);
            check(surface.polyData->GetNumberOfCells()==0,"Re-enabling internal group lost the empty selection");
            for(const auto id : internalIDs)
            {
                auto* item=surface.findChild<QCheckBox*>(QString("surfaceInternalFace_%1").arg(id));
                check(!item->isChecked(),"Internal checkbox selection was reset");
                item->setChecked(true);
            }
            double oldPoint[3];
            surface.polyData->GetPoint(0,oldPoint);
            base.simulationParameters.totalTime+=1.0;
            surface.updateConfiguration();
            double newPoint[3];
            surface.polyData->GetPoint(0,newPoint);
            check((Eigen::Map<Eigen::Vector3d>(oldPoint)-Eigen::Map<Eigen::Vector3d>(newPoint)).norm()>1.e-10,
                  "Internal displacement did not update with frame physics");
            surface.internalFaces->setChecked(false);
            check(surface.polyData->GetNumberOfCells()==0,"Hidden internal faces left stale geometry");
            surface.internalFaces->setChecked(true);
            check(surface.polyData->GetNumberOfCells()==static_cast<vtkIdType>(internalTriangles.size()),"Internal group did not restore selected faces");
        }
        std::filesystem::remove_all(dir);
        return 0;
    }
    catch(const std::exception& e)
    {
        std::cerr<<e.what()<<std::endl;
        return 1;
    }
}
