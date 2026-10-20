#include <DislocationDynamics/DislocationDynamicsBase.h>
#include <DislocationNetwork.h>
#include <DefectiveCrystal/MicrostructureContainer.h>
#include <DefectiveCrystal/DefectiveCrystalBase.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <iostream>

namespace model
{
    struct JunctionFormationTestAccess
    {
        static void glissile(DislocationNetwork<3>& network,double minimumLength)
        {
            network.junctionsMaker.glissileJunctions(minimumLength);
        }
    };
}

void require(bool value,const char* message)
{
    if(!value)
    {
        throw std::runtime_error(message);
    }
}

int main(int argc,char** argv)
{
    try
    {
        require(argc>=2,"Expected repository root");
#ifdef NDEBUG
        require(argc<3 || std::string(argv[2])!="--require-assertions","Assertions are disabled in this test build");
#else
        std::cout<<"Assertions enabled in periodic glissile-junction regression\n";
#endif
        const std::filesystem::path root(argv[1]);
        const auto folder=std::filesystem::temp_directory_path()/
            ("periodic-glissile-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        const auto input=folder/"inputFiles";
        std::filesystem::create_directories(input);
        std::filesystem::create_directories(folder/"F");
        std::filesystem::create_directories(folder/"config");
        std::filesystem::copy_file(root/"Library/Materials/W.txt",input/"W.txt");
        std::filesystem::copy_file(root/"Library/DislocationDynamics/DD.txt",input/"DD.txt");
        std::filesystem::copy_file(root/"Library/Meshes/unitCube24.msh",input/"mesh.msh");
        using Network=model::DislocationNetwork<3>;
        using V=Eigen::Vector3d;
        model::SimplicialMesh<3> mesh((input/"mesh.msh").string(),Eigen::Matrix3d::Identity(),V::Zero(),{-1});
        std::array<int,3> seeds{-1,-1,-1};
        for(const auto& region : mesh.regions())
        {
            for(const auto& face : region.second->faces())
            {
                const auto& pair=face.second->periodicFacePair;
                if(pair.second)
                {
                    for(int axis=0;axis<3;++axis)
                    {
                        if(pair.first[axis]>0.9 && pair.first.normalized().cross(V::Unit(axis)).norm()<1.e-8)
                        {
                            seeds[axis]=face.first;
                        }
                    }
                }
            }
        }
        require(seeds[0]>=0 && seeds[1]>=0,"Missing periodic face seeds");
        model::StaticID<model::PlanarMeshFace<3>>::force_count(0);
        const double length=200.0/std::sqrt(3.0);
        for(int scenario=0;scenario<4;++scenario)
        {
            const bool corner=scenario==1 || scenario==2;
            const bool twoEnds=scenario!=1;
            const int rank=corner? 3 : (scenario==3? 2 : 1);
            {
                std::ofstream f(input/"DefectiveCrystal.txt");
                f<<"physics=DislocationDynamics;\nuseFEM=0;\nNsteps=0;\nmaxResolveSteps=0;\ndtMax=1;\noutputFrequency=1;\noutputBinary=0;\nstartAtTimeStep=0;\nperiodicImageSize=";
                for(int k=0;k<rank;++k)
                {
                    f<<"0 ";
                }
                f<<";\nEwaldLengthFactor=0;\n";
                std::ofstream p(input/"polycrystal.txt");
                p<<std::setprecision(17)<<"materialFile=W.txt;\nmeshFile=mesh.msh;\nabsoluteTemperature=300;\nC2G1=1 0 0 0 1 0 0 0 1;\nF="<<length<<" 0 0 0 "<<length<<" 0 0 0 "<<length<<";\nX0=0 0 0;\nperiodicFaceIDs=";
                for(int k=0;k<rank;++k)
                {
                    p<<seeds[k]<<' ';
                }
                p<<";\n";
            }
            model::DefectiveCrystalBase<3> base(folder.string());
            model::MicrostructureContainer<3> container(base);
            Network network(container);
            const auto grain=base.poly.grain(1);
            const V direction=corner? V(1,1,1) : (scenario==3? V(1,1,0) : V(1,0,0));
            const V parentNormal=corner? V(1,-1,0).normalized().eval() : V::UnitZ().eval();
            const V side=parentNormal.cross(direction).normalized();
            const auto planeAt=[&](const V& point)
            {
                return base.dislocationDynamics().glidePlaneFactory.getFromKey(model::GlidePlaneKey<3>(point,grain->reciprocalLatticeDirection(parentNormal)));
            };
            const auto plane=planeAt(V::Zero());
            V burgers=V::Zero();
            for(const auto& [id,slip] : grain->slipSystems())
            {
                if(std::abs(slip->unitNormal.dot(direction))<1.e-8 &&
                   std::abs(slip->s.cartesian().dot(parentNormal))>1.e-8)
                {
                    burgers=slip->s.cartesian();
                    break;
                }
            }
            require(burgers.norm()>0.0,"No glissile Burgers vector for test direction");
            const auto insert=[&](const V& b,const auto& gp,const auto& netNodes)
            {
                const auto loop=network.loops().create(b,gp);
                std::vector<std::shared_ptr<Network::LoopNodeType>> nodes;
                for(const auto& node : netNodes)
                {
                    nodes.push_back(network.loopNodes().create(loop,node,node->get_P(),
                        loop->periodicGlidePlane->getPatch(V::Zero()),std::make_pair(nullptr,nullptr)));
                }
                network.insertLoop(loop,nodes);
                return loop;
            };
            const auto a=network.networkNodes().create(V((twoEnds? -0.3 : 0.0)*length*direction));
            const auto b=network.networkNodes().create(V(0.3*length*direction));
            const auto c=network.networkNodes().create(V(0.12*length*side));
            const auto d=network.networkNodes().create(V(-0.12*length*side));
            const auto first=insert(2.0*burgers,plane,std::vector{a,b,c});
            const auto second=insert(-burgers,plane,std::vector{a,b,d});
            network.updateBoundaryNodes();
            network.updateGeometry();
            network.trialNodeMotion=true;
            if(twoEnds)
            {
                a->trySet_P(-0.9*length*direction);
            }
            b->trySet_P(0.9*length*direction);
            network.updateBoundaryNodes();
            network.trialNodeMotion=false;
            network.updateGeometry();
            const double threshold=twoEnds? 0.8*length*direction.norm() : 0.0;
            std::shared_ptr<Network::NetworkLinkType> selected;
            std::vector<std::shared_ptr<Network::NetworkLinkType>> shared;
            for(const auto& entry : network.networkLinks())
            {
                const auto link=entry.second.lock();
                if(link->loopLinks().size()==2)
                {
                    shared.push_back(link);
                    if(link->chordLength()>threshold && !selected)
                    {
                        selected=link;
                    }
                }
            }
            require(selected && shared.size()>=2,"Missing shared periodic continuation");
            if(twoEnds)
            {
                for(const auto* occurrence : selected->loopLinks())
                {
                    require(occurrence->source->isPeriodicSubdivision() && occurrence->sink->isPeriodicSubdivision(),
                        "Test did not isolate a segment with two periodic endpoints");
                }
            }
            // An unrelated loop visits a periodic endpoint. It must not be
            // selected as a continuation, modified, or detached from that node.
            auto endpoint=selected->source;
            if(!(*selected->loopLinks().begin())->source->isPeriodicSubdivision())
            {
                endpoint=selected->sink;
            }
            const auto unrelated=network.loops().create(burgers,plane);
            std::map<size_t,std::shared_ptr<Network::NetworkNodeType>> copies;
            std::vector<std::shared_ptr<Network::LoopNodeType>> unrelatedNodes;
            for(const auto* link : first->linkSequence())
            {
                const auto& node=link->source;
                auto& copy=copies[node->networkNode->sID];
                if(!copy)
                {
                    copy=node->networkNode==endpoint? endpoint : network.networkNodes().clone(node->networkNode->sID);
                }
                unrelatedNodes.push_back(network.loopNodes().create(unrelated,copy,node->get_P(),
                    node->periodicPlanePatch(),node->periodicPlaneEdge));
            }
            network.insertLoop(unrelated,unrelatedNodes);
            unrelatedNodes.clear();
            copies.clear();
            network.updateGeometry();
            const V unrelatedArea=unrelated->currentAreaVector();
            const auto memberships=endpoint->loopNodes();
            const Eigen::Matrix3d stressBefore=network.stress(V(7,11,19),nullptr,nullptr,nullptr);
            model::JunctionFormationTestAccess::glissile(network,threshold);
            network.updateGeometry();
            require(network.loops().size()==4,"Expected exactly one glissile seed for the whole continuation");
            for(auto* occurrence : memberships)
            {
                require(endpoint->loopNodes().count(occurrence)==1,"Shared periodic endpoint was detached");
            }
            require(endpoint->loopNodes().size()>memberships.size(),"Seed did not share the periodic endpoint");
            require((unrelated->currentAreaVector()-unrelatedArea).norm()<1.e-8,"Unrelated loop at periodic endpoint changed");
            require((network.stress(V(7,11,19),nullptr,nullptr,nullptr)-stressBefore).norm()<1.e-9*std::max(1.e-12,stressBefore.norm()),
                "Zero-area junction seed changed the stress field");
            for(const auto& link : shared)
            {
                require(link->loopLinks().size()==3 && link->hasZeroBurgers(),"Seed lost junction sharing or Burgers cancellation");
            }
            std::shared_ptr<Network::LoopType> seed;
            for(const auto& entry : network.loops())
            {
                const auto loop=entry.second.lock();
                if(loop!=first && loop!=second && loop!=unrelated)
                {
                    seed=loop;
                }
            }
            require(seed && seed->isLoop() && seed->currentAreaVector().norm()<1.e-8,"Invalid glissile seed");
            size_t corners=0;
            for(const auto* link : seed->linkSequence())
            {
                const auto node=link->source;
                corners+=bool(node->periodicPlaneEdge.second);
                require((node->get_P()+node->periodicPlanePatch()->shift-node->networkNode->get_P()).norm()<1.e-7,"Incorrect seed patch translation");
                for(const auto& edge : {node->periodicPlaneEdge.first,node->periodicPlaneEdge.second})
                {
                    if(edge)
                    {
                        require(edge->patch==node->periodicPlanePatch().get() && edge->meshIntersection->contains(node->networkNode->get_P()),
                            "Invalid corner edge metadata");
                    }
                }
            }
            require(!corner || corners>0,"Corner crossing lost its second edge");
            // Move the free interior seed node and let the ordinary boundary
            // walk rebuild both sides. This exercises the stored corner metadata.
            for(const auto& entry : network.networkNodes())
            {
                entry.second.lock()->set_V(V::Zero(),false);
            }
            for(const auto* link : seed->linkSequence())
            {
                const auto node=link->source;
                if(!node->isPeriodicSubdivision() && node->networkNode->loopNodes().size()==1)
                {
                    node->networkNode->set_V(seed->glidePlane->unitNormal.cross(direction).normalized(),false);
                }
            }
            double sharedLengthBefore=0.0;
            for(const auto& link : shared)
            {
                sharedLengthBefore+=link->chordLength();
            }
            shared.clear();
            selected.reset();
            network.moveNodes(0.1);
            require(seed->currentAreaVector().norm()>1.e-8,"Glissile seed did not expand");
            double sharedLengthAfter=0.0;
            for(const auto& entry : network.networkLinks())
            {
                const auto segment=entry.second.lock();
                bool hasFirst=false,hasSecond=false,hasSeed=false;
                for(const auto* occurrence : segment->loopLinks())
                {
                    hasFirst=hasFirst || occurrence->loop==first;
                    hasSecond=hasSecond || occurrence->loop==second;
                    hasSeed=hasSeed || occurrence->loop==seed;
                }
                if(hasFirst && hasSecond && hasSeed)
                {
                    require(segment->hasZeroBurgers(),"Motion lost Burgers cancellation on shared continuation");
                    sharedLengthAfter+=segment->chordLength();
                }
            }
            require(std::abs(sharedLengthAfter-sharedLengthBefore)<1.e-7,"Motion lost sharing along the periodic continuation");
            for(const auto& entry : network.loopNodes())
            {
                const auto node=entry.second.lock();
                require((node->get_P()+node->periodicPlanePatch()->shift-node->networkNode->get_P()).norm()<1.e-7,
                    "Motion lost periodic mapping");
            }
            require((unrelated->currentAreaVector()-unrelatedArea).norm()<1.e-7,"Seed motion altered unrelated loop");
            require(first->isLoop() && second->isLoop() && unrelated->isLoop() && seed->isLoop(),"Motion broke loop closure");
            std::cout<<"PASS periodic glissile scenario "<<scenario<<'\n';
        }
        std::filesystem::remove_all(folder);
    }
    catch(const std::exception& error)
    {
        std::cerr<<error.what()<<'\n';
        return 1;
    }
}
