#include <DefectiveCrystal/DefectiveCrystalBase.h>
#include <DislocationNetwork.h>
#include <DefectiveCrystal/MicrostructureContainer.h>
#include <filesystem>
#include <fstream>
#include <chrono>
#include <iostream>
#include <numbers>

void require(bool condition,const char* message)
{
    if(!condition)
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
        const auto fixture=root/"tests/data/polycrystal1103";
        const auto folder=std::filesystem::temp_directory_path()/
            ("polycrystal-enclosed-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        const auto input=folder/"inputFiles";
        std::filesystem::create_directories(input);
        std::filesystem::create_directories(folder/"F");
        std::filesystem::create_directories(folder/"config");
        for(const auto* name : {"polycrystal.txt","DD.txt"})
        {
            std::filesystem::copy_file(fixture/name,input/name);
        }
        std::filesystem::copy_file(root/"Library/Materials/W.txt",input/"W.txt");
        std::filesystem::copy_file(root/"Library/Meshes/poly30_periodicXYZ.msh",input/"poly30_periodicXYZ.msh");
        {
            std::ofstream out(input/"DefectiveCrystal.txt");
            out<<"physics=DislocationDynamics;\nfmmTreeDepth=0;\nuseFEM=0;\nNsteps=0;\nmaxResolveSteps=0;\ndtMax=1;\noutputFrequency=1;\noutputBinary=0;\nstartAtTimeStep=0;\nperiodicImageSize=0 0 0;\nEwaldLengthFactor=0;\n";
        }
        model::DefectiveCrystalBase<3> base(folder.string());
        model::MicrostructureContainer<3> container(base);
        model::DislocationNetwork<3> network(container);
        model::DefectiveCrystalConfigIO<3> saved(folder.string());
        std::ifstream stream(fixture/"evl_1103.txt");
        require(bool(stream),"Missing frame 1103 fixture");
        saved.readTxtStream(stream);
        network.initializeConfiguration(saved,std::ofstream(),std::ofstream());
        // Collision selection must compare planes from distinct grains in the
        // global frame, not cross their unrelated integer lattice coordinates.
        size_t crossGrainChecks=0;
        for(const auto& entry : network.networkLinks())
        {
            const auto segment=entry.second.lock();
            if(segment->glidePlanes().empty())
            {
                continue;
            }
            for(const auto& other : network.loops())
            {
                const auto target=other.second.lock()->glidePlane;
                if(!target)
                {
                    continue;
                }
                const auto pairs=segment->parallelAndCoincidentGlidePlanes({target.get()});
                for(const auto* source : segment->glidePlanes())
                {
                    if(&source->n.lattice!=&target->n.lattice)
                    {
                        bool found=false;
                        for(const auto& pair : pairs)
                        {
                            found=found || (pair.first==source && pair.second==target.get());
                        }
                        require(found==(source->unitNormal.cross(target->unitNormal).norm()<FLT_EPSILON),
                            "Cross-grain plane comparison used crystal coordinates");
                        ++crossGrainChecks;
                    }
                }
            }
            if(crossGrainChecks)
            {
                break;
            }
        }
        require(crossGrainChecks>0,"Missing cross-grain plane comparison coverage");
        using V=Eigen::Vector3d;
        auto loop=network.loops().at(87).lock();
        const auto area=[&]()->V
        {
            V result=V::Zero();
            for(const auto* link : loop->linkSequence())
            {
                result+=0.5*link->source->get_P().cross(link->sink->get_P());
            }
            return result;
        };
        auto patch=loop->periodicGlidePlane->getPatch(V::Zero());
        V expected=V::Zero();
        for(const auto& edge : patch->edges())
        {
            const V a=loop->periodicGlidePlane->referencePlane->globalPosition(*edge->source);
            const V b=loop->periodicGlidePlane->referencePlane->globalPosition(*edge->sink);
            expected+=0.5*a.cross(b);
        }
        if(expected.dot(area())<0)
        {
            expected=-expected;
        }
        require(expected.norm()>1,"Expected a resolved grain cross-section");
        network.moveNodes(168.9662550850772);
        const auto check=[&]()
        {
            require((area()-expected).norm()<1.e-8*expected.norm(),"Fully enclosed grain lost slipped area");
            require(loop->linkSequence().size()==patch->edges().size(),"Retained contour is not the grain boundary");
            for(const auto* link : loop->linkSequence())
            {
                require(link->source->isPhysicalBoundaryNode(),"Enclosed contour is not confined to physical faces");
            }
            for(const auto& entry : network.loopNodes())
            {
                const auto node=entry.second.lock();
                require((node->get_P()+node->periodicPlanePatch()->shift-node->networkNode->get_P()).norm()<1.e-7,
                        "Boundary update broke loop/network periodic mapping");
            }
        };
        check();
        const auto count=network.loopNodes().size();
        network.updateBoundaryNodes();
        check();
        require(network.loopNodes().size()==count,"Repeated boundary update changed topology");
        network.executeSingleGlideStepDiscreteEvents(1103);
        check();
        // Drive every eligible segment, so a missing GB exclusion cannot pass
        // merely because this snapshot happens to have zero resolved force.
        size_t drivenPoints=0;
        for(const auto& entry : network.networkLinks())
        {
            const auto segment=entry.second.lock();
            segment->createQuadraturePoints(false);
            if(segment->isGrainBoundarySegment())
            {
                require(segment->quadraturePoints().empty(),"GB segment generated glide quadrature points");
            }
            for(auto& point : segment->quadraturePoints())
            {
                point.velocity=V(0.1,0.2,0.3);
                ++drivenPoints;
            }
        }
        require(drivenPoints>0,"No mobile segments exercise the solver");
        for(const auto* link : loop->linkSequence())
        {
            require(bool(link->networkLink()),"Boundary loop lost its network segment");
            require(link->networkLink()->isGrainBoundarySegment(),"Expected loop entirely on GBs");
            require(link->source->networkNode->loopNodes().size()==1,"Stationary-loop fixture shares a mobile node");
        }
        const Eigen::VectorXd velocities=network.glideSolver->getNodeVelocities();
        require(velocities.allFinite() && velocities.norm()>0,"No finite mobile response");
        std::map<size_t,V> positions;
        for(const auto* link : loop->linkSequence())
        {
            const auto node=link->source->networkNode;
            const V velocity=velocities.segment<3>(3*node->networkID());
            require(velocity.norm()==0,"Fully GB-confined loop has nonzero solver velocity");
            // Verify that committing the solution clears any previous motion.
            node->set_V(V(0.1,0.2,0.3),false);
            node->set_V(velocity,false);
            require(node->get_V().norm()==0,"GB node retained its previous velocity");
            positions.emplace(node->sID,node->get_P());
        }
        const Eigen::Matrix3d expectedPlastic=-loop->burgers()*expected.transpose()/base.mesh.volume();
        const Eigen::Matrix3d fullPlastic=network.averagePlasticDistortion();
        require(expectedPlastic.allFinite() && expectedPlastic.norm()>1.e-12,
                "GB loop has no resolved plastic distortion");
        const auto checkPlastic=[&]()
        {
            require((loop->averagePlasticDistortion()-expectedPlastic).norm()<1.e-8*expectedPlastic.norm(),
                    "Stationary GB loop lost its plastic distortion");
            require((network.averagePlasticDistortion()-fullPlastic).norm()<1.e-8*expectedPlastic.norm(),
                    "Network plastic distortion changed without motion");
            loop->updateRates();
            require(loop->averagePlasticDistortionRate().allFinite() &&
                    loop->averagePlasticDistortionRate().norm()==0,
                    "Stationary GB loop has a nonzero plastic distortion rate");
        };
        checkPlastic();
        // Probe off the loop plane to avoid dislocation cores and symmetry zeros.
        V center=V::Zero();
        for(const auto& entry : positions)
        {
            center+=entry.second;
        }
        center/=positions.size();
        const double scale=std::sqrt(expected.norm());
        const V normal=expected.normalized();
        std::vector<V> probes;
        std::vector<Eigen::Matrix3d> fullStress;
        std::vector<Eigen::Matrix3d> boundaryStress;
        std::vector<V> fullDisplacement;
        std::vector<V> boundaryDisplacement;
        for(const double offset : {0.13,0.31,0.57})
        {
            const V point=center+offset*scale*normal;
            probes.push_back(point);
            fullStress.push_back(network.stress(point,nullptr,nullptr,nullptr));
            fullDisplacement.push_back(network.displacement(point,nullptr,nullptr,nullptr));
            V displacement=V::Zero();
            for(const auto& shift : base.periodicShifts)
            {
                displacement-=loop->solidAngle(point+shift)*loop->burgers()/(4.0*std::numbers::pi);
            }
            Eigen::Matrix3d contribution=Eigen::Matrix3d::Zero();
            for(const auto* link : loop->linkSequence())
            {
                const auto segment=link->networkLink();
                require(segment->quadraturePoints().empty(),"Stress-source GB segment has quadrature points");
                for(const auto& shift : base.periodicShifts)
                {
                    contribution+=segment->straight.stress(point+shift);
                    displacement+=segment->straight.displacement(point+shift);
                }
            }
            require(contribution.allFinite() && contribution.norm()>1.e-12,
                    "Stationary GB loop produces no resolved stress");
            boundaryStress.push_back(contribution);
            require(displacement.allFinite() && displacement.norm()>1.e-12,
                    "Stationary GB loop produces no resolved displacement");
            boundaryDisplacement.push_back(displacement);
        }
        // Hold unrelated loops fixed to isolate persistence of the enclosed loop.
        for(const auto& entry : network.networkNodes())
        {
            entry.second.lock()->set_V(V::Zero(),false);
        }
        for(int step=0;step<3;++step)
        {
            network.moveNodes(100.0);
            check();
            checkPlastic();
            for(const auto& entry : positions)
            {
                const auto node=network.networkNodes().at(entry.first).lock();
                require((node->get_P()-entry.second).norm()<1.e-10,"GB-confined loop moved");
            }
        }
        for(size_t k=0;k<probes.size();++k)
        {
            const Eigen::Matrix3d current=network.stress(probes[k],nullptr,nullptr,nullptr);
            require((current-fullStress[k]).norm()<1.e-8*boundaryStress[k].norm(),
                    "Stress changed while the GB loop remained stationary");
            const V displacement=network.displacement(probes[k],nullptr,nullptr,nullptr);
            require((displacement-fullDisplacement[k]).norm()<1.e-8*boundaryDisplacement[k].norm(),
                    "Displacement changed while the GB loop remained stationary");
        }
        // Test the public network field, not just the segment kernel: removing
        // this isolated loop must remove exactly its nonzero field contributions.
        network.deleteLoop(loop->sID);
        require(loop->linkSequence().empty(),"GB loop was not removed for the field comparison");
        // Release the test's owning reference so the erased loop is absent from
        // the network's weak-pointer loop container and its field sums.
        patch.reset();
        loop.reset();
        require((fullPlastic-network.averagePlasticDistortion()-expectedPlastic).norm()<1.e-8*expectedPlastic.norm(),
                "Network plastic distortion omitted or miscounted the GB loop");
        for(size_t k=0;k<probes.size();++k)
        {
            const Eigen::Matrix3d withoutLoop=network.stress(probes[k],nullptr,nullptr,nullptr);
            require((fullStress[k]-withoutLoop-boundaryStress[k]).norm()<1.e-8*boundaryStress[k].norm(),
                    "Network stress omitted or miscounted the stationary GB loop");
            const V displacementWithoutLoop=network.displacement(probes[k],nullptr,nullptr,nullptr);
            require((fullDisplacement[k]-displacementWithoutLoop-boundaryDisplacement[k]).norm()<1.e-8*boundaryDisplacement[k].norm(),
                    "Network displacement omitted or miscounted the stationary GB loop");
        }
        std::filesystem::remove_all(folder);
        std::cout<<"PASS: frame 1103 stationary grain-boundary loop retains stress, displacement and plastic distortion with zero distortion rate\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr<<error.what()<<'\n';
        return 1;
    }
}
