#ifdef _OPENMP
#include <omp.h>
#endif
#include <chrono>
#include "ClusterPeriodicFixture.h"
#include <ElasticDeformation.h>
#include <DefectiveCrystal/DefectiveCrystal.h>

// Exercise periodic sinks with an intentionally nonperiodic source stress.
// Line-tension-driven climb and the actual DD concentration kernels remain active.
struct PeriodicReferenceStressNetwork : DislocationNetwork<3>
{
    using DislocationNetwork<3>::DislocationNetwork;
    MatrixDim stress(const VectorDim& point,const NodeType* const,const ElementType* const,const SimplexDim* const) const override
    {
        return (1.e-8*point(0))*MatrixDim::Identity();
    }
};

int main(int argc,char** argv)
{
    try
    {
        require(argc==2,"Expected repository root");
        const fs::path root(argv[1]);
        for(int rank : {0,1,3})
        {
            const auto folder=fixture(root,rank,0,false);
            auto text=read(folder/"inputFiles/DefectiveCrystal.txt");
            parameter(text,"physics","ElasticDeformation ClusterDynamics DislocationDynamics");
            text+="\nsolvedQuantityMaxIterations=150;\nsolvedQuantityRelativeTolerance=1e-6;\nsolvedQuantityRelaxation=0.8;\n";
            std::ofstream(folder/"inputFiles/DefectiveCrystal.txt")<<text;
            auto dd=read(root/"Library/DislocationDynamics/DD.txt");
            parameter(dd,"glideSolverType","none");
            parameter(dd,"computeDDinteractions","0");
            std::ofstream(folder/"inputFiles/DD.txt")<<dd;
            auto ed=read(root/"Library/ElasticDeformation/ElasticDeformation.txt");
            parameter(ed,"useElasticDeformationFEM","1");
            std::ofstream(folder/"inputFiles/ElasticDeformation.txt")<<ed;
            DefectiveCrystalBase<3> base(folder.string());
            MicrostructureGenerator generator(base);
            const auto grain=base.poly.grain(1);
            const Eigen::Vector3d normal=Eigen::Vector3d::UnitZ();
            const auto key=GlidePlaneKey<3>(Eigen::Vector3d::Zero(),grain->reciprocalLatticeDirection(normal));
            const auto plane=base.dislocationDynamics().periodicGlidePlaneFactory.getFromKey(key);
            for(int sign : {-1,1})
            {
                const Eigen::Vector3d center(15.0*sign,0,0);
                std::vector<Eigen::Vector3d> polygon;
                for(int i=0;i<8;++i)
                {
                    const double angle=2.0*M_PI*i/8;
                    polygon.emplace_back(center+Eigen::Vector3d(6*std::cos(angle),6*std::sin(angle),0));
                }
                require(generator.insertJunctionLoop(polygon,plane,sign*std::sqrt(8.0/3.0)*normal,normal,center,1,DislocationLoopIO<3>::SESSILELOOP),"Failed loop insertion");
            }
            std::unique_ptr<MicrostructureContainer<3>> owner;
            if(rank==0)
            {
                owner=std::make_unique<DefectiveCrystal<3>>(base);
            }
            else
            {
                owner=std::make_unique<MicrostructureContainer<3>>(base);
                owner->push_back(std::make_shared<ElasticDeformation<3>>(*owner));
                owner->push_back(std::make_shared<ClusterDynamics<3>>(*owner));
                owner->push_back(std::make_shared<PeriodicReferenceStressNetwork>(*owner));
                owner->finalizeSolvedQuantities();
            }
            auto& crystal=*owner;
            crystal.initializeConfiguration(generator.config(),std::ofstream(),std::ofstream());
            auto network=crystal.getUniqueTypedMicrostructure<DislocationNetwork<3>>();
            auto cd=crystal.getUniqueTypedMicrostructure<ClusterDynamics<3>>();
            auto elastic=crystal.getUniqueTypedMicrostructure<ElasticDeformation<3>>();
            std::map<size_t,Eigen::Vector3d> positions;
            for(const auto& node : network->networkNodes())
            {
                positions.emplace(node.first,node.second.lock()->get_P());
            }
            const auto loops=network->loops().size();
            const double time=base.simulationParameters.totalTime;
            const Eigen::VectorXd initialVelocity=network->solvedQuantityState();
            const Eigen::VectorXd initialConcentration=cd->clusterDynamicsFEM->mobileClusters->dofVector();
            crystal.solvedQuantities.options.maxIterations=1;
            bool failed=false;
            try
            {
                crystal.solve();
            }
            catch(const std::runtime_error& error)
            {
                failed=std::string(error.what()).find("did not converge")!=std::string::npos;
                if(!failed)
                {
                    throw;
                }
            }
            require(failed,"Expected a deliberately unconverged solve to fail");
            require((network->solvedQuantityState()-initialVelocity).norm()==0.0,"Failure changed climb rates or velocities");
            require((cd->clusterDynamicsFEM->mobileClusters->dofVector()-initialConcentration).norm()<1.e-14*initialConcentration.norm(),"Failure changed concentration");
            crystal.solvedQuantities.options.maxIterations=150;
            crystal.solve();
            require(base.simulationParameters.totalTime==time && network->loops().size()==loops,"Coupled solve changed time or topology");
            for(const auto& node : network->networkNodes())
            {
                require(node.second.lock()->get_P()==positions.at(node.first),"Coupled solve moved a shared node");
            }
            bool coupled=false;
            for(const auto& group : crystal.solvedQuantities.groups)
            {
                if(group.cyclic)
                {
                    require(group.iterations>1 && group.error<=1.0,"Coupled group did not converge");
                    coupled=true;
                }
            }
            require(coupled,"DD-CD cycle was not detected");
            // Each periodic sink class receives one value, even though source
            // stresses differ at its copies. The total concentration, not the
            // FEM correction alone, must agree at those copies.
            const auto& sinkFEM=*cd->clusterDynamicsFEM;
            std::map<size_t,std::vector<size_t>> sinkClasses;
            for(const auto& sink : sinkFEM.sinkNormals)
            {
                sinkClasses[base.fe->periodicRepresentative(sink.first)].push_back(sink.first);
            }
            for(const auto& entry : sinkClasses)
            {
                Eigen::Matrix3d meanStress=Eigen::Matrix3d::Zero();
                for(const size_t id : entry.second)
                {
                    const auto& node=base.fe->node(id);
                    meanStress+=owner->stress(node.P0,&node,nullptr,nullptr)/double(entry.second.size());
                }
                const auto& normals=sinkFEM.sinkNormals.at(entry.second.front());
                double traction=0.0;
                for(const auto& normal : normals)
                {
                    traction+=normal.dot(meanStress*normal)/double(normals.size());
                }
                const auto expected=cd->cdp.boundaryMobileConcentration(meanStress.trace(),traction).matrix().transpose().eval();
                for(const size_t id : entry.second)
                {
                    require((cd->solvedQuantitySinkConcentrations.at(id)-expected).norm()<1.e-12*std::max(1.e-50,expected.norm()),"Periodic sink stress averaging mismatch");
                }
            }

#ifdef _OPENMP
            const int previousThreads=omp_get_max_threads();
            omp_set_num_threads(1);
            elastic->solveStress();
            const Eigen::VectorXd serialDisplacement=elastic->elasticDeformationFEM->u->dofVector();
            omp_set_num_threads(4);
            elastic->solveStress();
            require((serialDisplacement-elastic->elasticDeformationFEM->u->dofVector()).norm()==0.0,"Parallel ED traction assembly changed the solution");
            omp_set_num_threads(previousThreads);
#endif

            require(elastic->weakProblemSolver->factorizationCount()==1,"ED factorization repeated inside coupling");
            require(cd->clusterDynamicsFEM->mSolver->factorizationCount()==1,"CD factorization repeated inside coupling");
            require(cd->solvedQuantityInfluenceBuilds==2,"H must rebuild once per solve, including retries");
#ifdef _OPENMP
            const auto rebuildCD=[&](int threads)
            {
                omp_set_num_threads(threads);
                cd->solvedQuantityInfluence.clear();
                cd->solvedQuantitySourceIDs.clear();
                cd->solvedQuantityTargetRows.clear();
                cd->solvedQuantitySinkConcentrations.clear();
                cd->solvedQuantityBoundaryCache=true;
                cd->applyBoundaryConditions();
                cd->solvedQuantityBoundaryCache=false;
            };
            rebuildCD(1);
            const auto serialH=cd->solvedQuantityInfluence;
            const auto serialSinks=cd->solvedQuantitySinkConcentrations;
            const auto serialBCs=cd->clusterDynamicsFEM->mobileClusters->dirichletConditions();
            rebuildCD(4);
            for(size_t species=0;species<serialH.size();++species)
            {
                require((serialH[species]-cd->solvedQuantityInfluence[species]).norm()==0.0,"Parallel H assembly differs from serial");
            }
            for(const auto& entry : serialSinks)
            {
                require((entry.second-cd->solvedQuantitySinkConcentrations.at(entry.first)).norm()==0.0,"Parallel sink concentrations differ from serial");
            }
            require(serialBCs==cd->clusterDynamicsFEM->mobileClusters->dirichletConditions(),"Parallel CD boundary values differ from serial");
            omp_set_num_threads(previousThreads);
#endif
            // Compare H against independent segment evaluation with nonuniform source rates.
            const auto originalState=network->solvedQuantityState();
            Eigen::VectorXd probe=originalState;
            for(Eigen::Index n=0;n<probe.size()/(ClusterDynamics<3>::mSize+3);++n)
            {
                for(int species=0;species<ClusterDynamics<3>::mSize;++species)
                {
                    probe(n*(ClusterDynamics<3>::mSize+3)+species)=1.e-10*(1+species)*std::sin(double(n+1));
                }
            }
            network->restoreSolvedQuantityState(probe);
            for(const auto& node : cd->clusterDynamicsFEM->mobileClusters->fe().nodes())
            {
                const auto row=cd->solvedQuantityTargetRows.find(node.gID);
                if(row==cd->solvedQuantityTargetRows.end())
                {
                    continue;
                }
                const auto direct=network->mobileConcentration(node.P0,&node,nullptr,nullptr);
                for(int species=0;species<ClusterDynamics<3>::mSize;++species)
                {
                    double cached=0.0;
                    for(size_t column=0;column<cd->solvedQuantitySourceIDs.size();++column)
                    {
                        cached+=cd->solvedQuantityInfluence[species](row->second,column)*network->networkNodes().get(cd->solvedQuantitySourceIDs[column])->climbVelocityScalar(species);
                    }
                    require(std::abs(cached-direct(species))<1.e-12*std::max(std::abs(direct(species)),1.e-100),"H disagrees with segment concentration evaluation");
                }
            }
            cd->solvedQuantityBoundaryCache=true;
            const auto cachedStart=std::chrono::steady_clock::now();
            cd->applyBoundaryConditions();
            const double cachedSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-cachedStart).count();
            const auto cachedBCs=cd->clusterDynamicsFEM->mobileClusters->dirichletConditions();
            cd->solvedQuantityBoundaryCache=false;
            const auto directStart=std::chrono::steady_clock::now();
            cd->applyBoundaryConditions();
            const double directSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-directStart).count();
            const auto& directBCs=cd->clusterDynamicsFEM->mobileClusters->dirichletConditions();
            for(const auto& bc : cachedBCs)
            {
                require(std::abs(bc.second-directBCs.at(bc.first))<1.e-12*std::max(std::abs(bc.second),1.e-100),"Cached boundary values differ from direct evaluation");
            }
            std::cout<<"H cache boundary update: "<<cachedSeconds<<" s; direct: "<<directSeconds<<" s\n";
            network->restoreSolvedQuantityState(originalState);
            const auto state=network->solvedQuantityState();
            const auto concentration=cd->clusterDynamicsFEM->mobileClusters->dofVector().eval();
            // Independent full solve rebuilds quadrature and climb operators, including stresses.
            network->solve();
            const auto refreshed=network->solvedQuantityState();
            require((state-refreshed).norm()<1.e-4*std::max(state.norm(),1.e-30),"Cached stresses changed the converged climb result");
            cd->solve();
            require((concentration-cd->clusterDynamicsFEM->mobileClusters->dofVector()).norm()<1.e-5*std::max(concentration.norm(),1.e-100),"CD fixed point is inconsistent");
            fs::remove_all(folder);
        }
        std::cout<<"PASS: coupled DD/CD/ED, frozen geometry, cache equivalence, partial/full periodicity\n";
        return 0;
    }
    catch(const std::exception& e)
    {
        std::cerr<<e.what()<<'\n';
        return 1;
    }
}
