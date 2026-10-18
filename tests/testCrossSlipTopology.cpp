#include <MicrostructureGenerator.h>
#include <DislocationNetwork.h>
#include <MicrostructureContainer.h>
#include <DislocationDynamicsBase.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <iostream>
#include <set>
#include <array>
#include <sstream>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <DislocationCrossSlip.h>
namespace model
{
struct CrossSlipTestAccess
{
    template<class N> static void queue(N& network,const std::deque<std::shared_ptr<typename N::LoopNodeType>>& nodes,int slip)
    {
        network.crossSlipMaker.csNodes.clear();
        network.crossSlipMaker.csNodes.emplace_back(nodes,slip);
    }
};
}


void require(bool condition, const char* message)
{
    if(!condition)
    {
        throw std::runtime_error(message);
    }
}


// Compare the production parallel post-event refresh with the original serial
// refresh on the identical topology, avoiding restart and allocation differences.
void checkSerialRefresh(model::DislocationNetwork<3>& network)
{
    const auto capture=[&]()
    {
        std::vector<Eigen::Matrix<double,22,1>> values;
        for(const auto& entry : network.networkLinks())
        {
            for(const auto& point : entry.second.lock()->quadraturePoints())
            {
                Eigen::Matrix<double,22,1> value;
                value.head<9>()=Eigen::Map<const Eigen::Matrix<double,9,1>>(point.stress.data());
                value.segment<3>(9)=point.pkForce;
                value.segment<3>(12)=point.lineTensionForce;
                value.segment<3>(15)=point.stackingFaultForce;
                value.segment<3>(18)=point.velocity;
                value(21)=point.coreEnergyPerLength;
                values.push_back(value);
            }
        }
        return values;
    };
    const auto parallel=capture();
    require(!parallel.empty(),"Cross-slip refresh fixture has no quadrature points");
    const auto velocities=network.glideSolver->getNodeVelocities().eval();
    for(const auto& entry : network.networkLinks())
    {
        entry.second.lock()->createQuadraturePoints(false);
    }
    for(const auto& entry : network.networkLinks())
    {
        entry.second.lock()->updateQuadraturePoints(false);
    }
    const auto serial=capture();
    require(serial.size()==parallel.size(),"Parallel refresh changed quadrature count");
    for(size_t k=0;k<serial.size();++k)
    {
        require((serial[k].array()==parallel[k].array()).all(),"Parallel refresh changed quadrature fields or forces");
    }
    require((velocities.array()==network.glideSolver->getNodeVelocities().array()).all(),
            "Parallel refresh changed nodal velocities");
}

int main(int argc, char** argv)
{
    try
    {
#ifdef _OPENMP
        omp_set_num_threads(4);
#endif
        require(argc == 2, "Expected repository root");
        const std::filesystem::path root(argv[1]);
        const auto folder = std::filesystem::temp_directory_path() /
            ("cross-slip-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        const auto input = folder / "inputFiles";
        std::filesystem::create_directories(input);
        std::filesystem::create_directories(folder / "F");
        std::filesystem::copy_file(root / "Library/Materials/W.txt", input / "W.txt");

        std::filesystem::copy_file(root / "Library/Meshes/bicrystal_48.msh", input / "mesh.msh");
        std::filesystem::copy_file(root / "Library/DislocationDynamics/DD.txt", input / "DD.txt");
        {
            std::ifstream inputDD(input/"DD.txt");
            std::string contents((std::istreambuf_iterator<char>(inputDD)),std::istreambuf_iterator<char>());
            const auto pos=contents.find("crossSlipModel=0;");
            require(pos!=std::string::npos,"Missing cross-slip setting");
            contents.replace(pos,17,"crossSlipModel=1;");
            std::ofstream outputDD(input/"DD.txt");
            outputDD<<contents;
        }
        std::filesystem::create_directories(folder / "evl");
        model::SimplicialMesh<3> mesh((input/"mesh.msh").string(),Eigen::Matrix3d::Identity(),Eigen::Vector3d::Zero(),{-1});
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
                        if(pair.first.normalized().cross(Eigen::Vector3d::Unit(axis)).norm()<1.e-8)
                        {
                            seeds[axis]=face.first;
                        }
                    }
                }
            }
        }
        require(seeds[0]>=0 && seeds[1]>=0,"Missing periodic seeds");
        model::StaticID<model::PlanarMeshFace<3>>::force_count(0);
        const double length=3671.9477120460206;
        for(int scenario=0;scenario<6;++scenario)
        {
            const bool full=scenario==1;
            const bool periodic=scenario!=2;
            {
                std::ofstream f(input/"DefectiveCrystal.txt");
                f<<"physics=DislocationDynamics;\nuseFEM=0;\nNsteps=0;\nmaxResolveSteps=0;\ndtMax=1;\noutputFrequency=1;\noutputBinary=0;\nstartAtTimeStep=0;\nperiodicImageSize="<<(full ? "0 0 0" : "0")<<";\nEwaldLengthFactor=0;\n";
                std::ofstream p(input/"polycrystal.txt");
                p<<std::setprecision(17)<<"materialFile=W.txt;\nmeshFile=mesh.msh;\nabsoluteTemperature=300;\nC2G1=1 0 0 0 1 0 0 0 1;\nC2G2=0.8 0 0.6 0 1 0 -0.6 0 0.8;\nF="<<length<<" 0 0 0 "<<length<<" 0 0 0 "<<length<<";\nX0=0 0 0;\nperiodicFaceIDs=";
                if(periodic)
                {
                    p<<(full ? -1 : seeds[1]);
                }
                p<<";\n";
            }
            model::DislocationDynamicsBase<3> base(folder.string());
            model::MicrostructureContainer<3> container(base);
            using Network=model::DislocationNetwork<3>;
            using V=Eigen::Vector3d;
            Network network(container);
            const auto grain=base.poly.grain(1);
            const auto slip=grain->slipSystems().begin()->second;
            int other=-1;
            for(size_t k=0;k<grain->slipSystems().size();++k)
            {
                const auto candidate=grain->slipSystems().at(k);
                if((candidate->s-slip->s).squaredNorm()==0 && candidate->unitNormal.cross(slip->unitNormal).norm()>0.1)
                {
                    other=k;
                    break;
                }
            }
            require(other>=0,"No cross-slip system");
            const auto css=grain->slipSystems().at(other);
            const V b=slip->s.cartesian();
            const V center=scenario==4 ? V::Zero().eval() : V(length*V(-0.5,periodic ? 0.49 : 0.0,0));
            const auto plane=base.glidePlaneFactory.getFromKey(model::GlidePlaneKey<3>(center,slip->n));
            const long spacing=std::max(1L,static_cast<long>(20.0/css->n.planeSpacing()));
            const long height=std::lround(double(css->n.closestPlaneIndexOfPoint(center))/spacing)*spacing;
            model::PlanePlaneIntersection<3> intersection(plane->P,plane->unitNormal,height*css->n.planeSpacing()*css->unitNormal,css->unitNormal);
            V origin=intersection.P+(center-intersection.P).dot(intersection.d)*intersection.d;
            V along=intersection.d;
            if(scenario==4 && along.x()<0.0)
            {
                along=-along;
            }
            V side=plane->unitNormal.cross(along).normalized();
            if(scenario==4)
            {
                origin=-0.15*length*along;
                if(side.x()>0.0)
                {
                    side=-side;
                }
            }
            if(scenario==5)
            {
                origin+=(0.5*length+2.0-origin.y())/along.y()*along;
            }
            const V a=origin-0.15*length*along;
            const V middle=origin;
            const V end=origin+0.15*length*along;
            std::vector<V> polygon{a,middle,end,end+100.0*side,a+100.0*side};
            if(scenario==5)
            {
                polygon[1]-=4.0/side.y()*side;
            }
            if(scenario==3)
            {
                polygon[0]+=0.2*side;
                polygon[1]+=0.2*side;
                polygon[2]+=2.0*css->n.planeSpacing()/css->unitNormal.dot(side)*side;
            }
            model::MicrostructureGenerator generator(base);
            require(generator.insertJunctionLoop(polygon,base.periodicGlidePlaneFactory.get(plane->key),b,plane->unitNormal,plane->P,1,model::DislocationLoopIO<3>::GLISSILELOOP),"Parent insertion failed");
            network.initializeConfiguration(generator.config(),std::ofstream(),std::ofstream());

            const auto parent=network.loops().begin()->second.lock();
            std::deque<std::shared_ptr<Network::LoopNodeType>> branch;
            for(const auto& point : std::vector<V>{polygon[0],polygon[1],polygon[2]})
            {
                for(const auto* link : parent->linkSequence())
                {
                    if((link->source->get_P()-point).norm()<1.e-6 && !link->source->isPeriodicSubdivision())
                    {
                        branch.push_back(link->source);
                        break;
                    }
                }
            }
            require(branch.size()==3,"Missing branch anchors");
            // A second loop shares the final endpoint. Its plane either permits
            // alignment, or locks the displaced endpoint against that motion.
            const auto endpoint=branch.back()->networkNode;
            const V normal=scenario==3 ? css->unitNormal : plane->unitNormal;
            const auto extraPlane=base.glidePlaneFactory.getFromKey(model::GlidePlaneKey<3>(endpoint->get_P(),grain->reciprocalLatticeDirection(normal)));
            const auto extra=network.loops().create(b,extraPlane);
            V u=normal.cross(along).normalized();
            if(scenario==4 && u.x()>0.0)
            {
                u=-u;
            }
            const auto e1=network.networkNodes().create(V(endpoint->get_P()-30.0*along));
            const auto e2=network.networkNodes().create(V(endpoint->get_P()-30.0*along+20.0*u));
            std::vector<std::shared_ptr<Network::LoopNodeType>> ring;
            for(const auto& node : std::vector{endpoint,e1,e2})
            {
                ring.push_back(network.loopNodes().create(extra,node,node->get_P(),extra->periodicGlidePlane->getPatch(V::Zero()),std::make_pair(nullptr,nullptr)));
            }
            network.insertLoop(extra,ring);
            network.updateGeometry();
            const auto snapshot=[&]()
            {
                std::ostringstream out;
                out<<std::setprecision(17);
                for(const auto& entry : network.networkNodes())
                {
                    const auto node=entry.second.lock();
                    out<<node->sID<<' '<<node->get_P().transpose()<<' '<<node->get_V().transpose()<<'\n';
                }
                for(const auto& entry : network.loops())
                {
                    out<<entry.first<<'\n';
                    for(const auto* link : entry.second.lock()->linkSequence())
                    {
                        out<<link->source->sID<<' '<<link->source->networkNode->sID<<' '<<link->source->get_P().transpose()<<' '<<link->sink->sID<<'\n';
                    }
                }
                return out.str();
            };
            if(scenario!=3 && scenario!=5)
            {
                for(const auto& entry : network.networkLinks())
                {
                    const auto segment=entry.second.lock();
                    segment->createQuadraturePoints(false);
                    for(auto& point : segment->quadraturePoints())
                    {
                        point.pkForce=css->unitNormal.cross(b.normalized());
                    }
                }
                std::deque<model::CrossSlipBranch<Network::LoopNodeType>> discovered;
                parent->crossSlipBranches(discovered);
                bool found=false;
                for(const auto& candidate : discovered)
                {
                    bool hasFirst=false,hasLast=false;
                    for(const auto& node : candidate.first)
                    {
                        hasFirst=hasFirst || node==branch.front();
                        hasLast=hasLast || node==branch.back();
                    }
                    found=found || (hasFirst && hasLast);
                }
                if(!found)
                {
                    std::cerr<<"scenario "<<scenario<<" desired slip "<<other<<" discovered "<<discovered.size()<<"\n";
                    for(const auto* link : parent->linkSequence())
                    {
                        if(link->networkLink())
                        {
                            const auto selected=network.crossSlipModel->isCrossSlipLink(*link->networkLink());
                            std::cerr<<link->source->get_P().transpose()<<" -> "<<link->sink->get_P().transpose()<<" length "<<link->networkLink()->chordLength()<<" selection "<<selected.first<<" "<<selected.second.second<<"\n";
                        }
                    }
                }
                require(found,"Branch discovery lost periodic or physical endpoint");
            }
            const auto before=snapshot();
            const auto endpointID=endpoint->sID;
            const auto middleNode=branch[1]->networkNode;
            const V oldMiddle=middleNode->get_P();
            model::CrossSlipTestAccess::queue(network,branch,other);
            network.crossSlipMaker.execute();
            if(scenario==3)
            {
                require(snapshot()==before,"Rejected alignment changed network");
                continue;
            }
            checkSerialRefresh(network);
            if(scenario==5)
            {
                require((middleNode->get_P()-oldMiddle).norm()>0.5*length,"Alignment did not cross the periodic face");
                require(network.networkNodes().get(middleNode->sID)==middleNode,"Alignment lost network identity at seam");
            }
            branch.clear();
            ring.clear();
            require(network.loops().size()==3,"Cross-slip seed was discarded");
            require(network.networkNodes().get(endpointID)==endpoint && endpoint->loopNodes().size()==3,"Shared endpoint was detached");
            std::shared_ptr<Network::LoopType> seed;
            for(const auto& entry : network.loops())
            {
                const auto loop=entry.second.lock();
                if(loop!=parent && loop!=extra)
                {
                    seed=loop;
                }
            }
            require(seed && seed->loopLinks().size()>=4 && seed->slipSystem(),"Invalid zero-area seed");
            if(scenario==4)
            {
                bool physical=false;
                for(const auto* node : endpoint->loopNodes())
                {
                    if(node->loop()==seed)
                    {
                        physical=node->isPhysicalBoundaryNode();
                    }
                }
                require(physical,"Cross-slip lost physical endpoint metadata");
            }
            require(seed->currentAreaVector().norm()<1.e-6,"Seed added slipped area before motion");
            size_t seams=0,cancelled=0;
            for(const auto* link : seed->linkSequence())
            {
                seams+=link->source->isPeriodicSubdivision();
                if(link->networkLink() && link->networkLink()->loopLinks().size()==2)
                {
                    require(link->networkLink()->hasZeroBurgers(),"Parent branch Burgers vector not cancelled");
                    ++cancelled;
                }
                require((link->source->get_P()+link->source->periodicPlanePatch()->shift-link->source->networkNode->get_P()).norm()<1.e-6,"Incorrect periodic mapping");
            }
            require(cancelled>=2,"No cancellation on original branch");
            require((seams>0)==(periodic && scenario!=4),"Periodic crossing was not preserved");
            const auto expected=network.glideSolver->getNodeVelocities();
            size_t k=0;
            for(const auto& entry : network.networkNodes())
            {
                const auto node=entry.second.lock();
                const V stored=node->get_V();
                node->set_V(expected.segment<3>(3*k),false);
                require((stored-node->get_V()).norm()<1.e-7,"Stale nodal velocity after cross-slip");
                ++k;
            }
            // Drive only the free side of the seed on its new plane. This
            // checks the ordinary motion/boundary walk after seed creation.
            for(const auto& entry : network.networkNodes())
            {
                entry.second.lock()->set_V(V::Zero(),false);
            }
            for(const auto* link : seed->linkSequence())
            {
                const auto node=link->source->networkNode;
                if(node->loopNodes().size()==1)
                {
                    node->set_V(css->unitNormal.cross(along).normalized(),false);
                }
            }
            network.moveNodes(0.1);
            require(seed->currentAreaVector().norm()>1.e-6,"Seed failed to expand under motion");
            require(endpoint->loopNodes().size()==3,"Motion lost shared endpoint");
        }
        std::filesystem::remove_all(folder);
        std::cout<<"PASS: cross-slip seed, periodic mapping, shared endpoint, fresh velocities and atomic rejection\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr<<error.what()<<'\n';
        return 1;
    }
}
