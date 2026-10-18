#include <MicrostructureGenerator.h>
#include <DislocationNetwork.h>
#include <MicrostructureContainer.h>
#include <DislocationDynamicsBase.h>
#include <DislocationCrossSlip.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <iostream>

namespace model
{
    struct CrossSlipTestAccess
    {
        template<class N> static void queue(N& network,const std::deque<CrossSlipBranch<typename N::LoopNodeType>>& branches)
        {
            network.crossSlipMaker.csNodes=branches;
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
        require(argc==2,"Expected repository root");
        const std::filesystem::path root(argv[1]);
        const auto folder=std::filesystem::temp_directory_path()/
            ("thermal-branches-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        const auto input=folder/"inputFiles";
        std::filesystem::create_directories(input);
        std::filesystem::create_directories(folder/"F");
        std::filesystem::create_directories(folder/"evl");
        std::filesystem::copy_file(root/"Library/Materials/Zr.txt",input/"Zr.txt");
        std::filesystem::copy_file(root/"Library/Meshes/bicrystal_48.msh",input/"mesh.msh");
        {
            std::ofstream material(input/"Zr.txt",std::ios::app);
            // Zero barriers and volumes give a known, spatially uniform rate.
            material<<"\ncsw_SI=1e15;\ncsLb_SI=9.699e-8;\ncsLp_SI=9.699e-8;\ncsb2pE_eV=0;\ncsb2pVeg_SI=0;\ncsb2pVsc_SI=0;\ncsb2pVec_SI=0;\ncsp2bE_eV=0;\ncsp2bVeg_SI=0;\ncsp2bVsc_SI=0;\ncsp2bVec_SI=0;\n";
            std::ifstream dd(root/"Library/DislocationDynamics/DD.txt");
            std::string text((std::istreambuf_iterator<char>(dd)),std::istreambuf_iterator<char>());
            const auto pos=text.find("crossSlipModel=0;");
            require(pos!=std::string::npos,"Missing cross-slip setting");
            text.replace(pos,17,"crossSlipModel=2;");
            std::ofstream(input/"DD.txt")<<text;
            std::ofstream(input/"DefectiveCrystal.txt")<<"physics=DislocationDynamics;\nuseFEM=0;\nNsteps=0;\nmaxResolveSteps=0;\ndtMax=1;\noutputFrequency=1;\noutputBinary=0;\nstartAtTimeStep=0;\nperiodicImageSize=;\nEwaldLengthFactor=0;\n";
        }
        const double box=3671.9477120460206;
        std::vector<double> referenceLocations;
        for(const int refinement : {1,2,12,200})
        {
            {
                std::ofstream p(input/"polycrystal.txt");
                p<<std::setprecision(17)<<"materialFile=Zr.txt;\nmeshFile=mesh.msh;\nabsoluteTemperature=300;\nC2G1="<<(refinement==200? "0 -1 0 1 0 0 0 0 1" : "1 0 0 0 1 0 0 0 1")<<";\nC2G2=1 0 0 0 1 0 0 0 1;\nF="<<box<<" 0 0 0 "<<box<<" 0 0 0 "<<box<<";\nX0=0 0 0;\nperiodicFaceIDs=;\n";
            }
            using Network=model::DislocationNetwork<3>;
            using V=Eigen::Vector3d;
            model::DislocationDynamicsBase<3> base(folder.string());
            model::MicrostructureContainer<3> container(base);
            Network network(container);
            const auto grain=base.poly.grain(1);
            const auto slip=grain->slipSystems().begin()->second;
            const auto plane=base.glidePlaneFactory.getFromKey(model::GlidePlaneKey<3>(V(-0.5*box,0,0),slip->n));
            const V center=plane->snapToPlane(V(-0.5*box,0,0));
            const V along=slip->s.cartesian().normalized();
            const V side=plane->unitNormal.cross(along);
            const V a=center-0.15*box*along,end=center+0.15*box*along;
            std::vector<V> polygon;
            for(int k=0;k<=refinement;++k)
            {
                polygon.push_back(a+(end-a)*(double(k)/refinement));
            }
            polygon.push_back(end+100.0*side);
            polygon.push_back(a+100.0*side);
            model::MicrostructureGenerator generator(base);
            require(generator.insertJunctionLoop(polygon,base.periodicGlidePlaneFactory.get(plane->key),slip->s.cartesian(),plane->unitNormal,plane->P,1,model::DislocationLoopIO<3>::GLISSILELOOP),"Parent insertion failed");
            network.initializeConfiguration(generator.config(),std::ofstream(),std::ofstream());
            const auto parent=network.loops().begin()->second.lock();
            base.simulationParameters.dt=1.0;
            auto* thermal=dynamic_cast<model::EscaigCrossSlipModelHEX<Network>*>(network.crossSlipModel.get());
            require(thermal,"Missing production HEX thermal model");
            double totalRate=0.0,totalLength=0.0;
            for(const auto* link : parent->linkSequence())
            {
                if(link->networkLink())
                {
                    link->networkLink()->createQuadraturePoints(false);
                    for(auto& q : link->networkLink()->quadraturePoints())
                    {
                        q.stress.setZero();
                    }
                    if(thermal->isCrossSlipLink(*link->networkLink()).first)
                    {
                        totalLength+=link->networkLink()->chordLength();
                        for(const auto& [target,rate] : thermal->targetRates(*link->networkLink()))
                        {
                            totalRate+=rate;
                        }
                    }
                }
            }
            require(std::abs(totalLength-0.6*box)<1.e-6,"Thermal eligibility depends on remeshing length");
            require(totalRate>0.0 && std::abs(totalRate-thermal->csw*totalLength/thermal->csLb)<1.e-8*totalRate,"Integrated branch hazard depends on segmentation");
            std::deque<model::CrossSlipBranch<Network::LoopNodeType>> candidates;
            parent->crossSlipBranches(candidates,0.0);
            require(candidates.empty(),"Thermal sampler ignored the supplied timestep");
            thermal->gen.seed(77);
            parent->crossSlipBranches(candidates,1.0);
            require(candidates.size()==2,"Branch discovery depends on segmentation");
            std::vector<double> locations;
            for(const auto& candidate : candidates)
            {
                require(candidate.thermal && std::abs(candidate.thermal->length-300.0)<1.e-9,"Wrong physical nucleus length");
                const auto target=grain->slipSystems().at(candidate.second);
                require(target->n.base()(2)==0,"Wrong thermal target plane family");
                locations.push_back(candidate.thermal->centerFraction);
            }
            std::sort(locations.begin(),locations.end());
            if(referenceLocations.empty())
            {
                referenceLocations=locations;
            }
            else
            {
                for(size_t k=0;k<locations.size();++k)
                {
                    require(std::abs(locations[k]-referenceLocations[k])<1.e-9,"Nucleation position depends on segmentation");
                }
            }
            if(refinement==2)
            {
                // A point junction divides the screw branch even though its
                // adjacent segments each still belong to only one loop.
                std::shared_ptr<Network::LoopNodeType> junction;
                for(const auto* link : parent->linkSequence())
                {
                    if((link->source->get_P()-center).norm()<1.e-6)
                    {
                        junction=link->source;
                    }
                }
                require(bool(junction),"Missing junction anchor");
                const auto other=network.loops().create(slip->s.cartesian(),plane);
                std::vector<std::shared_ptr<Network::LoopNodeType>> nodes;
                for(const auto& net : std::vector{junction->networkNode,
                    network.networkNodes().create(V(center+30.0*side)),
                    network.networkNodes().create(V(center+30.0*side+30.0*along))})
                {
                    nodes.push_back(network.loopNodes().create(other,net,net->get_P(),
                        other->periodicGlidePlane->getPatch(V::Zero()),std::make_pair(nullptr,nullptr)));
                }
                network.insertLoop(other,nodes);
                network.updateGeometry();
                candidates.clear();
                parent->crossSlipBranches(candidates,1.0);
                require(candidates.size()==3,"Thermal branch did not stop at a shared junction");
                for(const auto& candidate : candidates)
                {
                    for(size_t k=1;k+1<candidate.first.size();++k)
                    {
                        require(candidate.first[k]->networkNode!=junction->networkNode,"Junction remained inside thermal branch");
                    }
                }
                continue;
            }
            model::CrossSlipTestAccess::queue(network,candidates);
            candidates.clear();
            network.crossSlipMaker.execute();
            require(network.loops().size()==3,"Production thermal seeds were not inserted");
            for(const auto& entry : network.loops())
            {
                const auto loop=entry.second.lock();
                if(loop!=parent)
                {
                    double perimeter=0.0;
                    for(const auto* link : loop->linkSequence())
                    {
                        perimeter+=(link->sink->get_P()-link->source->get_P()).norm();
                    }
                    require(std::abs(perimeter-600.0)<1.e-6,"Production seed length depends on segmentation");
                }
            }
        }
        std::filesystem::remove_all(folder);
        std::cout<<"PASS: production HEX branch discovery, rates and nucleus geometry invariant under refinement\n";
    }
    catch(const std::exception& error)
    {
        std::cerr<<error.what()<<'\n';
        return 1;
    }
}
