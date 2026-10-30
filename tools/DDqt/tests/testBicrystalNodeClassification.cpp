#include <NetworkNodeActor.h>
#include <NetworkLinkActor.h>
#include <vtkActorCollection.h>
#include <vtkTubeFilter.h>
#include <DislocationNetwork.h>
#include <DefectiveCrystal/MicrostructureContainer.h>
#include <DefectiveCrystal/DefectiveCrystalBase.h>
#include <QApplication>
#include <vtkRenderer.h>
#include <vtkPointData.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <iostream>
#include <iterator>

void require(bool condition,const char* message)
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
        require(argc == 2, "Expected repository root");
        const std::filesystem::path root(argv[1]);
        const auto folder = std::filesystem::temp_directory_path() /
            ("bicrystal-node-classification-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        const auto input = folder / "inputFiles";
        std::filesystem::create_directories(input);
        std::filesystem::create_directories(folder / "F");
        std::filesystem::copy_file(root / "Library/Materials/W.txt", input / "W.txt");

        std::filesystem::copy_file(root / "Library/Meshes/bicrystal_48.msh", input / "mesh.msh");
        std::filesystem::copy_file(root / "Library/DislocationDynamics/DD.txt", input / "DD.txt");
        std::filesystem::create_directories(folder / "config");
        const double length=3671.9477120460206;
        {
            std::ofstream f(input / "DefectiveCrystal.txt");
            f << "physics=DislocationDynamics;\nfmmTreeDepth=0;\nuseFEM=0;\nNsteps=0;\nmaxResolveSteps=0;\ndtMax=1;\noutputFrequency=1;\noutputBinary=0;\nstartAtTimeStep=0;\nperiodicImageSize=0 0 0;\nEwaldLengthFactor=0;\n";
            std::ofstream p(input / "polycrystal.txt");
            p << std::setprecision(17) << "materialFile=W.txt;\nmeshFile=mesh.msh;\nabsoluteTemperature=300;\nC2G1=1 0 0 0 1 0 0 0 1;\n"
              << "C2G2=0.8 0 0.6 0 1 0 -0.6 0 0.8;\nF=" << length << " 0 0 0 " << length << " 0 0 0 " << length
              << ";\nX0=0 0 0;\nperiodicFaceIDs=-1;\n";
        }
        // Match the tutorial remeshing length, without depending on its inputs.
        {
            const auto dd=input/"DD.txt";
            std::ifstream source(dd);
            std::string contents((std::istreambuf_iterator<char>(source)),std::istreambuf_iterator<char>());
            const auto location=contents.find("Lmax=100;");
            require(location!=std::string::npos,"Expected library remeshing default");
            contents.replace(location,std::string("Lmax=100;").size(),"Lmax=150;");
            std::ofstream output(dd);
            output<<contents;
        }

        auto window=vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
        window->SetReadyForRendering(false);
        auto renderer=vtkSmartPointer<vtkRenderer>::New();
        window->AddRenderer(renderer);
        model::NetworkNodeActor actor(window,renderer);
        for(int frame : {430,500})
        {
            model::DefectiveCrystalBase<3> base(folder.string());
            model::DefectiveCrystal<3> crystal(base);
            auto& network=*crystal.getUniqueTypedMicrostructure<model::DislocationNetwork<3>>();
            model::DefectiveCrystalConfigIO<3> saved(folder.string());
            std::ifstream stream(root/"tests/data"/("bicrystalBoundary"+std::to_string(frame)+".txt"));
            require(bool(stream),"Missing bicrystal fixture");
            saved.readTxtStream(stream);
            network.initializeConfiguration(saved,std::ofstream(),std::ofstream());
            size_t internal=0,periodicGB=0,sameGrain=0,bulk=0,corners=0,stale=0;
            for(const auto& entry : network.networkNodes())
            {
                const auto node=entry.second.lock();
                bool internalFace=false,crossGrainFace=false,sameGrainFace=false,external=false;
                for(const auto* face : node->meshFaces())
                {
                    if(face->regionIDs.first!=face->regionIDs.second)
                    {
                        internalFace=true;
                    }
                    else
                    {
                        external=true;
                        const auto* twin=face->periodicFacePair.second;
                        if(twin)
                        {
                            crossGrainFace|=face->regionIDs.first!=twin->regionIDs.first;
                            sameGrainFace|=face->regionIDs.first==twin->regionIDs.first;
                        }
                    }
                }
                const bool gb=internalFace || crossGrainFace;
                const int expected=gb?2:(external?1:0);
                require(node->isGrainBoundaryNode()==gb,"Wrong physical GB classification");
                require(node->isOnInternalBoundary()==internalFace,"Geometric internal-face classification changed");
                require(node->isBoundaryNode()==external,"Geometric external-face classification changed");
                require(node->meshLocation()==expected,"GB does not take priority at periodic boundary");
                const model::DislocationNodeIO<3> output(*node);
                require(output.meshLocation==expected,"New output does not serialize updated classification");
                internal+=internalFace;
                periodicGB+=crossGrainFace;
                sameGrain+=sameGrainFace && !gb;
                bulk+=!gb && !external;
                corners+=gb && sameGrainFace;
            }
            require(internal && periodicGB && sameGrain && bulk && corners,
                    "Bicrystal fixture lacks required node categories");
            {
                size_t internalSegments=0,periodicSegments=0,boundarySegments=0;
                for(const auto& entry : network.networkLinks())
                {
                    const auto link=entry.second.lock();
                    bool gb=false;
                    for(const auto* face : link->meshFaces())
                    {
                        gb|=face->regionIDs.first!=face->regionIDs.second ||
                            (face->periodicFacePair.second && face->regionIDs!=face->periodicFacePair.second->regionIDs);
                    }
                    require(link->isGrainBoundarySegment()==gb,"Wrong GB segment classification");
                    if(gb)
                    {
                        require(link->meshLocation()==2,"GB segment has boundary-only output label");
                    }
                    if(!link->hasZeroBurgers())
                    {
                        internalSegments+=link->isOnInternalBoundary();
                        periodicSegments+=gb && !link->isOnInternalBoundary();
                        boundarySegments+=gb || link->isBoundarySegment();
                    }
                }
                require(internalSegments && periodicSegments,"Fixture lacks both GB segment types");
                actor.updateConfiguration(saved,&network);
                const int actorCount=renderer->GetActors()->GetNumberOfItems();
                model::NetworkLinkActor links(window,renderer,crystal);
                links.updateConfiguration(actor.nodePolyData,{});
                auto* boundaryActor=vtkActor::SafeDownCast(renderer->GetActors()->GetItemAsObject(actorCount+1));
                require(boundaryActor && boundaryActor->GetMapper(),"Missing boundary segment actor");
                auto* tube=vtkTubeFilter::SafeDownCast(boundaryActor->GetMapper()->GetInputAlgorithm());
                auto* data=vtkPolyData::SafeDownCast(tube->GetInputDataObject(0,0));
                require(data && data->GetNumberOfLines()==static_cast<vtkIdType>(boundarySegments),"DDqt omitted GB segments");
                auto* show=links.findChild<QCheckBox*>("showSegments");
                require(bool(show),"Missing segment visibility control");
                auto* opacity=links.findChild<QSlider*>("boundarySegmentOpacity");
                require(opacity && opacity->value()==0 && !boundaryActor->GetVisibility(),
                        "Boundary segments should default to transparent");
                opacity->setValue(50);
                require(boundaryActor->GetVisibility() && std::abs(boundaryActor->GetProperty()->GetOpacity()-0.5)<1.e-12,
                        "Boundary opacity control did not update the actor");
                show->setChecked(false);
                require(!boundaryActor->GetVisibility(),"Segments control did not hide GB segments");
                show->setChecked(true);
                require(boundaryActor->GetVisibility(),"Segments control did not show GB segments");
            }
            actor.updateConfiguration(saved,&network);
            auto* colors=actor.nodePolyData->GetPointData()->GetScalars();
            require(colors && colors->GetNumberOfTuples()==static_cast<vtkIdType>(saved.nodes().size()),
                    "Wrong color tuple count");
            size_t i=0;
            for(const auto& savedNode : saved.nodes())
            {
                const auto node=network.networkNodes().at(savedNode.sID).lock();
                const int location=node->meshLocation();
                for(int component=0;component<3;++component)
                {
                    require(colors->GetComponent(i,component)==actor.nodeClr[location][component],
                            "DDqt uses stale serialized classification");
                }
                stale+=savedNode.meshLocation!=location;
                ++i;
            }
            require(stale>0,"Fixture does not exercise legacy boundary labels");
            // With no DD network available, retain support for serialized-only data.
            actor.updateConfiguration(saved);
            i=0;
            for(const auto& node : saved.nodes())
            {
                const int index=node.meshLocation>=0 && node.meshLocation<=2?node.meshLocation:3;
                for(int component=0;component<3;++component)
                {
                    require(actor.nodePolyData->GetPointData()->GetScalars()->GetComponent(i,component)==actor.nodeClr[index][component],
                            "Serialized-only node coloring changed");
                }
                ++i;
            }
            std::cout<<"Frame "<<frame<<": internal GB "<<internal<<", periodic GB "<<periodicGB
                     <<", same-grain seam "<<sameGrain<<", corners "<<corners<<", corrected labels "<<stale<<std::endl;
        }
        std::filesystem::remove_all(folder);
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr<<error.what()<<std::endl;
        return 1;
    }
}
