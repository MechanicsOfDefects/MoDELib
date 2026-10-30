#include <DislocationNetwork.h>
#include <NetworkLinkActor.h>
#include <NetworkLoopActor.h>
#include <QApplication>
#include <vtkActorCollection.h>
#include <vtkCellData.h>
#include <vtkRenderer.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <iostream>

void require(bool condition,const char* message)
{
    if(!condition)
    {
        throw std::runtime_error(message);
    }
}

void checkLine(vtkPolyData* data,vtkIdType index,vtkIdType source,vtkIdType sink)
{
    vtkIdType count=0;
    const vtkIdType* ids=nullptr;
    data->GetLines()->GetNextCell(count,ids);
    require(index<data->GetNumberOfLines() && count==2 && ids[0]==source && ids[1]==sink,
            "Cached connectivity differs from legacy networkID ordering");
}

void checkActors(model::DefectiveCrystal<3>& crystal,bool changeTopology)
{
    auto network=crystal.getUniqueTypedMicrostructure<model::DislocationNetwork<3>>();
    auto window=vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
    window->SetReadyForRendering(false);
    auto renderer=vtkSmartPointer<vtkRenderer>::New();
    window->AddRenderer(renderer);
    model::NetworkLinkActor links(window,renderer,crystal);
    const int loopStart=renderer->GetActors()->GetNumberOfItems();
    model::NetworkLoopActor loops(window,renderer,crystal);
    const auto verify=[&]()
    {
        const auto nodeCount=network->networkNodes().size();
        const auto loopNodeCount=network->loopNodes().size();
        std::map<std::pair<size_t,size_t>,std::set<size_t>> sharing;
        auto points=vtkSmartPointer<vtkPoints>::New();
        for(const auto& entry : network->networkNodes())
        {
            points->InsertNextPoint(entry.second.lock()->get_P().data());
        }
        auto input=vtkSmartPointer<vtkPolyData>::New();
        input->SetPoints(points);
        for(const auto& entry : network->networkLinks())
        {
            for(const auto* ll : entry.second.lock()->loopLinks())
            {
                sharing[entry.first].insert(ll->loop->sID);
            }
        }
        links.updateConfiguration(input,{});
        loops.updateConfiguration({});
        vtkPolyData* data[3];
        vtkIdType counts[3]={0,0,0};
        for(int group=0;group<3;++group)
        {
            auto* actor=vtkActor::SafeDownCast(renderer->GetActors()->GetItemAsObject(group));
            auto* tube=vtkTubeFilter::SafeDownCast(actor->GetMapper()->GetInputAlgorithm());
            require(tube,"Missing segment tube filter");
            data[group]=vtkPolyData::SafeDownCast(tube->GetInputDataObject(0,0));
            require(data[group] && data[group]->GetPoints()==points,"Segment points changed");
            data[group]->GetLines()->InitTraversal();
        }
        for(const auto& entry : network->networkLinks())
        {
            const auto link=entry.second.lock();
            const int group=link->hasZeroBurgers()?2:((link->isBoundarySegment() || link->isGrainBoundarySegment())?1:0);
            checkLine(data[group],counts[group],link->source->networkID(),link->sink->networkID());
            if(group!=2)
            {
                const auto expected=links.computeColor(link->burgers(),link->chord(),link->slipSystem(),{});
                auto* colors=group==0?data[group]->GetCellData()->GetArray("Colors"):data[group]->GetCellData()->GetScalars();
                require(colors,"Missing segment colors");
                for(int k=0;k<3;++k)
                {
                    require(colors->GetComponent(counts[group],k)==static_cast<unsigned char>(expected(k)),"Segment color changed");
                }
            }
            std::set<size_t> memberships;
            for(const auto* ll : link->loopLinks())
            {
                memberships.insert(ll->loop->sID);
            }
            require(memberships==sharing.at(entry.first),"Actor update changed junction sharing");
            ++counts[group];
        }
        for(int group=0;group<3;++group)
        {
            require(counts[group]==data[group]->GetNumberOfLines(),"Segment count changed");
        }
        auto* loopActor=vtkActor::SafeDownCast(renderer->GetActors()->GetItemAsObject(loopStart+1));
        auto* loopData=vtkPolyData::SafeDownCast(loopActor->GetMapper()->GetInputDataObject(0,0));
        require(loopData && loopData->GetNumberOfPoints()==static_cast<vtkIdType>(loopNodeCount),"Loop occurrences were merged");
        require(loopData->GetNumberOfLines()==static_cast<vtkIdType>(network->loopLinks().size()),"Loop line count changed");
        loopData->GetLines()->InitTraversal();
        vtkIdType i=0;
        for(const auto& entry : network->loopLinks())
        {
            checkLine(loopData,i++,entry.second.source->networkID(),entry.second.sink->networkID());
        }
        for(const auto& entry : network->loopNodes())
        {
            const auto node=entry.second.lock();
            double xyz[3];
            loopData->GetPoint(node->networkID(),xyz);
            for(int k=0;k<3;++k)
            {
                require(xyz[k]==static_cast<float>(node->get_P()(k)),"Loop point position changed");
            }
        }
        require(network->networkNodes().size()==nodeCount && network->loopNodes().size()==loopNodeCount,
                "Actor update changed topology");
    };
    verify();
    if(changeTopology)
    {
        std::shared_ptr<model::DislocationNetwork<3>::NetworkLinkType> shared;
        for(const auto& entry : network->networkLinks())
        {
            const auto link=entry.second.lock();
            if(link->loopLinks().size()>1 && link->chordLength()>1 && !link->isBoundarySegment())
            {
                shared=link;
                break;
            }
        }
        require(bool(shared),"Fixture has no shared junction");
        const auto memberships=shared->loopLinks().size();
        auto inserted=network->networkNodes().create(shared->get_r(0.5));
        network->expandNetworkLink(shared,inserted);
        shared.reset();
        require(inserted->loopNodes().size()==memberships,"Expansion lost junction sharing");
        network->updateGeometry();
        verify();
        network->removeNetworkNode(inserted->sID);
        inserted.reset();
        network->updateGeometry();
        verify();
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
            ("actor-index-caches-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
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

        for(const std::string fixture : {"bicrystalBoundary430.txt","bicrystalBoundary500.txt","bicrystalSharedConstraints1915.txt"})
        {
            model::DefectiveCrystalBase<3> base(folder.string());
            model::DefectiveCrystal<3> crystal(base);
            auto network=crystal.getUniqueTypedMicrostructure<model::DislocationNetwork<3>>();
            model::DefectiveCrystalConfigIO<3> saved(folder.string());
            std::ifstream stream(root/"tests/data"/fixture);
            require(bool(stream),"Missing bicrystal fixture");
            saved.readTxtStream(stream);
            network->initializeConfiguration(saved,std::ofstream(),std::ofstream());
            checkActors(crystal,fixture=="bicrystalSharedConstraints1915.txt");
        }
        std::filesystem::remove_all(folder);
        std::cout<<"PASS: actor points, connectivity, colors and junction sharing match legacy indexing\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr<<error.what()<<std::endl;
        return 1;
    }
}
