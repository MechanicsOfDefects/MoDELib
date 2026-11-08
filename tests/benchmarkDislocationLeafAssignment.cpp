#ifdef MODELIB_BENCHMARK_KOKKOS
#include <FMM/KokkosBoxParticleIndex.h>
#include <DislocationDynamics/KokkosDislocationFieldBuffers.h>
#endif
#include <DefectiveCrystal/DefectiveCrystal.h>
#include <DislocationDynamics/DislocationFieldSampling.h>
#include <FMM/BoxParticleIndex.h>
#include <chrono>
#include <iomanip>
using Clock=std::chrono::steady_clock;
template<class F> double timing(F&& f,int repeats)
{
    const auto start=Clock::now();
    for(int i=0;i<repeats;++i) { f(); }
    return std::chrono::duration<double,std::milli>(Clock::now()-start).count()/repeats;
}
int main(int argc,char** argv)
{
    try
    {
        if(argc<3 || argc>4) { throw std::invalid_argument("Usage: benchmarkDislocationLeafAssignment tutorial frame [repeats]"); }
        const int repeats=argc==4?std::stoi(argv[3]):20;
        if(repeats<1) { throw std::invalid_argument("Positive repeat count required"); }
#ifdef MODELIB_BENCHMARK_KOKKOS
        Kokkos::ScopeGuard guard;
#endif
        model::DefectiveCrystalBase<3> base(argv[1]);
        model::DefectiveCrystalConfigIO<3> config(base.simulationParameters.traitsIO.configFolder);
        config.read(std::stol(argv[2]));
        model::MicrostructureContainer<3> container(base);
        model::DislocationNetwork<3> network(container);
        network.initializeConfiguration(config,std::ofstream(),std::ofstream());
        network.updateGeometry();
        const double density=model::TextFileParser(base.simulationParameters.traitsIO.inputFilesFolder+"/DD.txt").readScalar<double>("quadPerLength",true);
        const auto legacyCreate=[&]
        {
            for(const auto& entry:network.networkLinks())
            {
                entry.second.lock()->createQuadraturePoints(false);
            }
        };
        const double legacyFirst=timing(legacyCreate,1);
        const double legacyAverage=timing(legacyCreate,repeats);
        std::size_t legacyCount=0;
        for(const auto& entry:network.networkLinks())
        {
            legacyCount+=entry.second.lock()->quadraturePoints().size();
        }
        std::cout<<"LEGACY_QUADRATURE frame="<<argv[2]<<" points="<<legacyCount
                 <<" first_ms="<<legacyFirst<<" creation_ms="<<legacyAverage
                 <<" elastic_energy="<<network.computeElasticEnergyPerLength<<std::endl;
        auto& evaluator=network.fieldEvaluator;
        double sampling=0,assignment=0;
        const auto firstTimes=evaluator.rebuild(density,density,model::DislocationTargetPurpose::glide);
        const double first=firstTimes.assignmentMilliseconds;
        for(int i=0;i<repeats;++i)
        {
            const auto times=evaluator.rebuild(density,density,model::DislocationTargetPurpose::glide);
            sampling+=times.samplingMilliseconds/repeats;
            assignment+=times.assignmentMilliseconds/repeats;
        }
        const auto& buffers=evaluator.buffers();
        const auto& index=evaluator.leafIndex();
        const auto& tree=evaluator.tree();
        // Verify every point against the existing topology lookup, outside timing.
        const auto verify=[&](const auto& positions,const auto& indices)
        {
            for(std::size_t p=0;p<positions[0].size();++p)
            {
                Eigen::Vector3d x(positions[0][p],positions[1][p],positions[2][p]);
                if(tree.locate(x,tree.depth())!=tree.boxes(tree.depth())[indices.leafOfPoint[p]])
                { throw std::runtime_error("Leaf verification failed"); }
            }
        };
        verify(buffers.sources().positions.values,index.sources());
        verify(buffers.targets().positions.values,index.targets());
        std::cout<<std::setprecision(6)<<"LEAF_BENCHMARK frame="<<argv[2]<<" links="<<network.networkLinks().size()
                 <<" sources="<<buffers.sources().size()<<" targets="<<buffers.targets().size()
                 <<" depth="<<tree.depth()<<" leaves="<<tree.boxes(tree.depth()).size()<<" repeats="<<repeats
                 <<" sampling_ms="<<sampling<<" first_assignment_ms="<<first<<" assignment_ms="<<assignment<<std::endl;
#ifdef MODELIB_BENCHMARK_KOKKOS
        const auto upload=[&]{evaluator.upload();};
        const double firstUpload=timing(upload,1);
        const double uploadMs=timing(upload,repeats);
        std::cout<<"LEAF_BENCHMARK backend="<<Kokkos::DefaultExecutionSpace::name()<<" first_upload_ms="<<firstUpload<<" upload_ms="<<uploadMs<<std::endl;
#endif
    }
    catch(const std::exception& e) { std::cerr<<e.what()<<std::endl; return 1; }
}
