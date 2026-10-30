#include <DislocationDynamics/DislocationDynamicsBase.h>
#include <DislocationNetwork.h>
#include <DefectiveCrystal/MicrostructureContainer.h>
#include <DefectiveCrystal/DefectiveCrystalBase.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <iostream>
#include <set>

void require(bool condition, const char* message)
{
    if(!condition)
    {
        throw std::runtime_error(message);
    }
}

int main(int argc, char** argv)
{
    try
    {
        require(argc == 2, "Expected repository root");
        const std::filesystem::path root(argv[1]);
        const auto folder = std::filesystem::temp_directory_path() /
            ("periodic-fields-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        const auto input = folder / "inputFiles";
        std::filesystem::create_directories(input);
        std::filesystem::create_directories(folder / "F");
        std::filesystem::copy_file(root / "Library/Materials/W.txt", input / "W.txt");
        const auto meshFile = root / "Library/Meshes/unitCube24.msh";
        std::filesystem::copy_file(meshFile, input / "mesh.msh");
        model::SimplicialMesh<3> mesh(meshFile.string(), Eigen::Matrix3d::Identity(), Eigen::Vector3d::Zero(), {-1});
        std::vector<int> seeds;
        for(int axis = 0; axis < 3; ++axis)
        {
            for(const auto& region : mesh.regions())
            {
                for(const auto& face : region.second->faces())
                {
                    const auto& pair = face.second->periodicFacePair;
                    if(pair.second && pair.first[axis] > 0.9 &&
                       pair.first.normalized().cross(Eigen::Vector3d::Unit(axis)).norm() < 1.e-8)
                    {
                        seeds.push_back(face.first);
                    }
                }
            }
        }
        require(seeds.size() == 3, "Expected three cube periodic directions");

        std::filesystem::copy_file(root / "Library/DislocationDynamics/DD.txt", input / "DD.txt");
        const double length=200.0/std::sqrt(3.0);
        for(int activeAxis : {-1,0,1})
        {
            const bool periodic=activeAxis==0;
            {
                std::ofstream f(input / "DefectiveCrystal.txt");
                f << "physics=DislocationDynamics;\nfmmTreeDepth=0;\nuseFEM=0;\nNsteps=0;\nmaxResolveSteps=0;\ndtMax=1;\n"
                  << "outputFrequency=1;\noutputBinary=0;\nstartAtTimeStep=0;\nperiodicImageSize=0;\nEwaldLengthFactor=0;\n";
                std::ofstream p(input / "polycrystal.txt");
                p << std::setprecision(17) << "materialFile=W.txt;\nmeshFile=mesh.msh;\nabsoluteTemperature=300;\n"
                  << "C2G1=1 0 0 0 1 0 0 0 1;\nF=" << length << " 0 0 0 " << length << " 0 0 0 " << length
                  << ";\nX0=0 0 0;\nperiodicFaceIDs=";
                if(activeAxis>=0)
                {
                    p << seeds[activeAxis];
                }
                p << ";\n";
            }
            model::DefectiveCrystalBase<3> base(folder.string());
            model::MicrostructureContainer<3> container(base);
            using Network=model::DislocationNetwork<3>;
            using V=Eigen::Vector3d;
            Network network(container);
            const auto grain=base.poly.grain(1);
            const auto plane=base.dislocationDynamics().glidePlaneFactory.getFromKey(model::GlidePlaneKey<3>(V::Zero(),grain->reciprocalLatticeDirection(V::UnitZ())));
            const V burgers=2.0/std::sqrt(3.0)*V::UnitX();
            const auto first=network.loops().create(burgers,plane);
            const auto second=network.loops().create(2.0*burgers,plane);
            // Two oppositely oriented loops share AB. B is beyond the face in
            // loop space, but its network representative is inside the RVE.
            const V a=length*V(0.2,0,0);
            const V b=length*V(0.7,0,0);
            const auto na=network.networkNodes().create(a);
            const auto nb=network.networkNodes().create(length*V(periodic ? -0.3 : 0.49,0,0));
            const auto makeNode=[&](const auto& loop,const auto& net,const V& position)
            {
                return network.loopNodes().create(loop,net,position,loop->periodicGlidePlane->getPatch(V::Zero()),
                    std::make_pair(std::shared_ptr<model::PeriodicPlaneEdge<3>>(),std::shared_ptr<model::PeriodicPlaneEdge<3>>()));
            };
            auto a1=makeNode(first,na,a);
            auto b1=makeNode(first,nb,b);
            auto a2=makeNode(second,na,a);
            auto b2=makeNode(second,nb,b);
            const auto c1=makeNode(first,network.networkNodes().create(length*V(0.1,0.2,0)),length*V(0.1,0.2,0));
            const auto c2=makeNode(second,network.networkNodes().create(length*V(0.1,-0.2,0)),length*V(0.1,-0.2,0));
            network.insertLoop(first,std::vector<std::shared_ptr<Network::LoopNodeType>>{a1,b1,c1});
            network.insertLoop(second,std::vector<std::shared_ptr<Network::LoopNodeType>>{b2,a2,c2});
            // An unrelated coincident polygon must not acquire junction sharing.
            const auto unrelated=network.loops().create(burgers,plane);
            network.insertLoop(unrelated,std::vector<std::shared_ptr<Network::LoopNodeType>>{
                makeNode(unrelated,network.networkNodes().create(a),a),
                makeNode(unrelated,network.networkNodes().create(nb->get_P()),b),
                makeNode(unrelated,network.networkNodes().create(c1->get_P()),c1->get_P())});
            // A trimmed network node may still belong to an unaffected loop.
            const auto survivor=network.loops().create(burgers,plane);
            network.insertLoop(survivor,std::vector<std::shared_ptr<Network::LoopNodeType>>{
                makeNode(survivor,nb,nb->get_P()),
                makeNode(survivor,network.networkNodes().create(length*V(0,0.1,0)),length*V(0,0.1,0)),
                makeNode(survivor,network.networkNodes().create(length*V(0,0.3,0)),length*V(0,0.3,0))});
            // A real mesh-search failure must not change a network position
            // before throwing (set_P is also used by the commit transaction).
            const V originalPosition=na->get_P();
            const auto originalSimplex=na->includingSimplex();
            bool rejected=false;
            try
            {
                na->set_P(length*V(10,10,10));
            }
            catch(const std::runtime_error&)
            {
                rejected=true;
            }
            require(rejected && (na->get_P()-originalPosition).norm()==0 && na->includingSimplex()==originalSimplex,
                    "Failed position validation changed a committed network node");
            std::map<size_t,std::pair<const Network::LoopNodeType*,V>> oldLoopNodes;
            std::map<size_t,std::pair<const Network::NetworkNodeType*,V>> oldNetworkNodes;
            std::map<Network::LoopLinkType::KeyType,const Network::LoopLinkType*> oldLinks;
            std::map<const Network::NetworkLinkType*,V> oldBurgers;
            std::map<size_t,V> oldDisplacements;
            for(const auto& entry : network.loopNodes())
            {
                const auto node=entry.second.lock();
                oldLoopNodes.emplace(entry.first,std::make_pair(node.get(),node->get_P()));
            }
            for(const auto& entry : network.networkNodes())
            {
                const auto node=entry.second.lock();
                oldNetworkNodes.emplace(entry.first,std::make_pair(node.get(),node->get_P()));
                oldDisplacements.emplace(entry.first,node->get_deltaP());
            }
            for(const auto& entry : network.loopLinks())
            {
                oldLinks.emplace(entry.first,&entry.second);
            }
            for(const auto& entry : network.networkLinks())
            {
                const auto link=entry.second.lock();
                oldBurgers.emplace(link.get(),link->burgers());
            }
            size_t failAt=1;
            for(;failAt<500;++failAt)
            {
                size_t operations=0;
                try
                {
                    network.updateBoundaryNodes([&]()
                    {
                        if(++operations==failAt)
                        {
                            throw std::runtime_error("Injected boundary commit failure");
                        }
                    });
                    break;
                }
                catch(const std::runtime_error& error)
                {
                    require(std::string(error.what())=="Injected boundary commit failure","Unexpected update failure");
                }
                require(network.loopNodes().size()==oldLoopNodes.size(),"Failed update leaked loop nodes");
                require(network.networkNodes().size()==oldNetworkNodes.size(),"Failed update leaked network nodes");
                require(network.loopLinks().size()==oldLinks.size(),"Failed update changed connectivity");
                require(network.networkLinks().size()==oldBurgers.size(),"Failed update changed junction count");
                for(const auto& entry : oldLoopNodes)
                {
                    const auto node=network.loopNodes().get(entry.first);
                    require(node.get()==entry.second.first && (node->get_P()-entry.second.second).norm()==0,
                            "Failed update changed loop-node identity or position");
                    require(!node->periodicContourVertex,"Failed update changed contour flags");
                }
                for(const auto& entry : oldNetworkNodes)
                {
                    const auto node=network.networkNodes().get(entry.first);
                    require(node.get()==entry.second.first && (node->get_P()-entry.second.second).norm()==0,
                            "Failed update changed network-node identity or position");
                    require((node->get_deltaP()-oldDisplacements.at(entry.first)).norm()==0,
                            "Failed update changed the stored motion increment");
                }
                for(const auto& entry : oldLinks)
                {
                    require(&network.loopLinks().at(entry.first)==entry.second,"Rollback replaced an original link");
                    require(entry.second->source->next.second==entry.second && entry.second->sink->prev.second==entry.second,
                            "Rollback broke directed incidence");
                }
                for(const auto& entry : oldBurgers)
                {
                    require((entry.first->burgers()-entry.second).norm()==0,"Rollback changed junction Burgers vector");
                }
            }
            require(failAt>5 && failAt<500,"Commit fault injection did not complete");

            require(first->isLoop() && second->isLoop(),"Clipped topology is not closed");
            size_t shared=0;
            size_t junctionSegments=0;
            std::set<size_t> ids;
            for(const auto* link : first->linkSequence())
            {
                ids.insert(link->source->sID);
                const auto n=link->source->networkNode;
                if(std::abs(n->get_P().x()-0.5*length)<1.e-6 && std::abs(n->get_P().y())<1.e-6)
                {
                    require(n->loopNodes().size()==2,"Junction intersection is not shared");
                    ++shared;
                }
                if(!periodic)
                {
                    require(link->source!=b1,"Trimmed original node remains in loop");
                    require(link->source->get_P().x()<=0.5*length+1.e-6,"Contour outside physical face");
                }
                if(const auto segment=link->networkLink())
                {
                    if(segment->loopLinks().size()==2)
                    {
                        ++junctionSegments;
                        V expected=V::Zero();
                        for(const auto* participant : segment->loopLinks())
                        {
                            const V flow=participant->loop==first ? burgers : 2.0*burgers;
                            expected+=(participant->source->networkNode==segment->source ? flow : -flow);
                        }
                        require((segment->burgers()-expected).norm()<1.e-8,"Junction Burgers orientation changed");
                    }
                }
            }
            require(junctionSegments>0,"Shared junction segments disappeared");
            require(shared==1,"Expected one shared intersection on positive face");
            if(!periodic)
            {
                require(b1->prev.second==nullptr && b1->next.second==nullptr,"Removed node retains links");
                require(b2->prev.second==nullptr && b2->next.second==nullptr,"Removed partner retains links");
            }
            b1.reset();
            b2.reset();
            if(!periodic)
            {
                require(nb->loopNodes().size()==1,"Trimming removed unrelated membership or retained obsolete nodes");
            }
            network.updateBoundaryNodes();
            require(survivor->isLoop(),"Trimming damaged unaffected loop");
            std::set<size_t> after;
            for(const auto* link : first->linkSequence())
            {
                after.insert(link->source->sID);
            }
            require(ids==after,"Repeated update replaced stable contour nodes");
            if(periodic)
            {
                // A stale loop subdivision is replaced without breaking the
                // other loop's ownership of its shared network representative.
                std::shared_ptr<Network::LoopNodeType> stale;
                for(const auto* link : first->linkSequence())
                {
                    if(link->source->periodicPlaneEdge.first &&
                       std::abs(link->source->networkNode->get_P().x()-0.5*length)<1.e-6 &&
                       std::abs(link->source->get_P().y())<1.e-6)
                    {
                        stale=link->source;
                        break;
                    }
                }
                require(bool(stale),"Missing periodic subdivision");
                const auto sharedNetwork=stale->networkNode;
                network.danglingBoundaryLoopNodes.insert(stale.get());
                network.updateBoundaryNodes();
                require(stale->prev.second==nullptr && stale->next.second==nullptr,"Stale subdivision remains connected");
                stale.reset();
                require(sharedNetwork->loopNodes().size()==2,"Replacing stale subdivision broke junction sharing");
                require(network.danglingBoundaryLoopNodes.empty(),"Dangling node registry not cleared");
            }
        }
        std::filesystem::remove_all(folder);
        std::cout << "PASS: shared junctions across periodic and physical faces, trimming, repeated updates\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
