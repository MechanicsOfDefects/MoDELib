#include <MicrostructureGenerator.h>
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
        require(argc == 2 || argc == 3, "Expected repository root");
        const std::filesystem::path root(argv[1]);
        const auto folder = std::filesystem::temp_directory_path() /
            ("bicrystal-boundary-solver-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        const auto input = folder / "inputFiles";
        std::filesystem::create_directories(input);
        std::filesystem::create_directories(folder / "F");
        std::filesystem::copy_file(root / "Library/Materials/W.txt", input / "W.txt");

        std::filesystem::copy_file(root / "Library/Meshes/bicrystal_48.msh", input / "mesh.msh");
        std::filesystem::copy_file(root / "Library/DislocationDynamics/DD.txt", input / "DD.txt");
        std::filesystem::create_directories(folder / "config");
        const std::string mode=argc==3 ? argv[2] : "full";
        const bool partial=mode=="partial";
        const bool physicalOnly=mode=="physical";
        std::string faceIDs="-1";
        if(partial)
        {
            model::SimplicialMesh<3> mesh((input/"mesh.msh").string(),Eigen::Matrix3d::Identity(),Eigen::Vector3d::Zero(),{-1});
            bool found=false;
            for(const auto& region : mesh.regions())
            {
                for(const auto& face : region.second->faces())
                {
                    const auto& pair=face.second->periodicFacePair;
                    if(face.second->isExternal() && pair.second &&
                       pair.first.normalized().cross(Eigen::Vector3d::UnitY()).norm()<1.e-8)
                    {
                        faceIDs=std::to_string(face.first);
                        found=true;
                    }
                }
            }
            require(found,"No Y periodic seed found");
            model::StaticID<model::Polyhedron<3,3-1>>::force_count(0);
        }
        if(physicalOnly)
        {
            faceIDs="";
        }
        const double length=3671.9477120460206;
        {
            std::ofstream f(input / "DefectiveCrystal.txt");
            f << "physics=DislocationDynamics;\nuseFEM=0;\nNsteps=0;\nmaxResolveSteps=0;\ndtMax=1;\noutputFrequency=1;\noutputBinary=0;\nstartAtTimeStep=0;\nperiodicImageSize=" << (partial ? "0" : physicalOnly ? "" : "0 0 0") << ";\nEwaldLengthFactor=0;\n";
            std::ofstream p(input / "polycrystal.txt");
            p << std::setprecision(17) << "materialFile=W.txt;\nmeshFile=mesh.msh;\nabsoluteTemperature=300;\nC2G1=1 0 0 0 1 0 0 0 1;\n"
              << "C2G2=0.8 0 0.6 0 1 0 -0.6 0 0.8;\nF=" << length << " 0 0 0 " << length << " 0 0 0 " << length
              << ";\nX0=0 0 0;\nperiodicFaceIDs=" << faceIDs << ";\n";
        }
        model::DefectiveCrystalBase<3> base(folder.string());
        model::MicrostructureGenerator generator(base);
        model::PrismaticLoopIndividualSpecification spec;
        spec.slipSystemIDs={0};
        spec.loopRadii={1.e-7};
        spec.glideSteps={50};
        spec.loopCenters.resize(1,3);
        spec.loopCenters << -0.5*length,0,0;
        generator.addPrismaticLoopIndividual(spec);
        require(generator.config().loops().size()==7,"Expected complete prismatic construction");
        model::MicrostructureContainer<3> container(base);
        model::DislocationNetwork<3> network(container);
        network.initializeConfiguration(generator.config(),std::ofstream(),std::ofstream());
        // Duplicate the glide loops using the same network nodes. This makes
        // the surviving mobile branches genuine shared junction segments.
        std::vector<std::shared_ptr<model::DislocationNetwork<3>::LoopType>> originals;
        for(const auto& entry : network.loops())
        {
            const auto loop=entry.second.lock();
            if(loop->loopType==model::DislocationLoopIO<3>::GLISSILELOOP)
            {
                originals.push_back(loop);
            }
        }
        for(const auto& original : originals)
        {
            const auto duplicate=network.loops().create(original->burgers(),original->glidePlane);
            std::vector<std::shared_ptr<model::DislocationNetwork<3>::LoopNodeType>> nodes;
            for(const auto* link : original->linkSequence())
            {
                const auto* node=link->source.get();
                nodes.push_back(network.loopNodes().create(duplicate,node->networkNode,node->get_P(),
                    node->periodicPlanePatch(),node->periodicPlaneEdge));
            }
            network.insertLoop(duplicate,nodes);
        }
        network.alphaLineTension=1.0;
        require(!network.stochasticForceGenerator,"Regression requires deterministic mobility");
        using V=Eigen::Vector3d;
        if(argc==3 && !partial && !physicalOnly)
        {
            model::DefectiveCrystalConfigIO<3> saved(folder.string());
            std::ifstream stream(argv[2]);
            require(bool(stream),"Cannot read saved configuration");
            saved.readTxtStream(stream);
            network.setConfiguration(saved);
            network.moveNodes(10.0);
        }
        else
        {
            const V center(-0.5*length,0,0);
            network.trialNodeMotion=true;
            for(const auto& entry : network.networkNodes())
            {
                const auto node=entry.second.lock();
                // Retain bulk branches while reaching both physical and periodic faces.
                node->trySet_P(center+(physicalOnly ? 1.0 : 60.0)*(node->get_P()-center));
            }
            network.updateBoundaryNodes();
            network.trialNodeMotion=false;
            network.updateGeometry();
        }
        size_t sharedPhysical=0;
        size_t periodic=0;
        for(const auto& entry : network.networkNodes())
        {
            const auto node=entry.second.lock();
            for(const auto* ln : node->loopNodes())
            {
                if(ln->isPhysicalBoundaryNode() && (physicalOnly || node->isBoundaryNode()) && node->loopNodes().size()>1)
                {
                    ++sharedPhysical;
                }
                if(ln->isPeriodicSubdivision())
                {
                    ++periodic;
                }
            }
        }
        require(physicalOnly || sharedPhysical>0,"No shared physical boundary node");
        require(physicalOnly ? periodic==0 : periodic>0,"Unexpected same-grain periodic crossing count");
        require(base.periodicLatticeBasis.cols()==(physicalOnly ? 0 : (partial ? 1 : 3)),"Incorrect active periodic rank");
        for(const auto& entry : network.networkLinks())
        {
            const auto segment=entry.second.lock();
            segment->createQuadraturePoints(false);
            // A spatially varying prescribed velocity exposes missing interpolation constraints.
            for(auto& point : segment->quadraturePoints())
            {
                point.velocity=(1.0+0.15*std::sin(point.r.sum()/length))*V(0.1,0.2,0.3);
            }
        }
        const Eigen::VectorXd velocities=network.glideSolver->getNodeVelocities();
        require(velocities.size()==static_cast<Eigen::Index>(3*network.networkNodes().size()),"Incorrect solver output size");
        require(velocities.allFinite(),"Nonfinite constrained velocity");
        require(velocities.norm()>0,"Solver test has no driven nodes");
        if(physicalOnly)
        {
            // Independent reference: original unconstrained lumping must be
            // unchanged when no periodic subdivisions exist.
            std::deque<Eigen::Triplet<double>> triplets;
            Eigen::VectorXd force=Eigen::VectorXd::Zero(velocities.size());
            Eigen::VectorXd mass=Eigen::VectorXd::Zero(velocities.size());
            for(const auto& entry : network.networkLinks())
            {
                const auto segment=entry.second.lock();
                if(!segment->isBoundarySegment())
                {
                    segment->addToGlobalAssembly(triplets,force);
                }
            }
            for(const auto& entry : triplets)
            {
                mass(entry.row())+=0.5*entry.value();
                mass(entry.col())+=0.5*entry.value();
            }
            Eigen::VectorXd expected=Eigen::VectorXd::Zero(velocities.size());
            for(Eigen::Index k=0;k<mass.size();++k)
            {
                if(std::abs(mass(k))>FLT_EPSILON)
                {
                    expected(k)=force(k)/mass(k);
                }
            }
            require((velocities-expected).norm()<1.e-10,"Nonperiodic solve changed unconstrained velocities");
        }

        size_t movingPhysical=0;
        size_t movingSharedPhysical=0;
        for(const auto& entry : network.networkNodes())
        {
            const auto node=entry.second.lock();
            const V velocity=velocities.segment<3>(3*node->networkID());
            bool physical=false;
            bool subdivision=!node->loopNodes().empty();
            for(const auto* ln : node->loopNodes())
            {
                physical=physical || ln->isPhysicalBoundaryNode();
                subdivision=subdivision && ln->isPeriodicSubdivision();
            }
            if(physical && velocity.norm()>1.e-8)
            {
                ++movingPhysical;
                if(node->loopNodes().size()>1)
                {
                    ++movingSharedPhysical;
                }
                require(velocity.dot(V(0.1,0.2,0.3))>0,
                        "Physical endpoint lost response to mobile segments");
            }
            if(subdivision)
            {
                for(const auto* ln : node->loopNodes())
                {
                    const auto prev=ln->periodicPrev();
                    const auto next=ln->periodicNext();
                    if(prev && next)
                    {
                        const double a=(ln->get_P()-prev->get_P()).norm();
                        const double b=(ln->get_P()-next->get_P()).norm();
                        const V expected=(b*velocities.segment<3>(3*prev->networkNode->networkID())
                                         +a*velocities.segment<3>(3*next->networkNode->networkID()))/(a+b);
                        require((velocity-expected).norm()<1.e-8,"Periodic endpoint interpolation changed");
                    }
                }
            }
            node->set_V(velocity,true);
            if(physical)
            {
                for(const auto& normal : node->physicalBoundaryNormals())
                {
                    require(node->get_V().dot(normal)>=-1.e-8,
                            "Trial velocity allows inward detachment from a physical face");
                }
            }
        }
        require(physicalOnly || movingPhysical>0,"Physical nodes were all pinned by the solver");
        require(physicalOnly || movingSharedPhysical>0,"Shared physical junctions lost motion from mobile segments");
        std::filesystem::remove_all(folder);
        std::cout << "PASS: bicrystal physical junctions and periodic crossings in glide solver\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
