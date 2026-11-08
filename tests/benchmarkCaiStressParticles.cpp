#include <DislocationDynamics/CaiStressParticleKernel.h>
#include <DislocationDynamics/KokkosDislocationFieldBuffers.h>
#include <DefectiveCrystal/DefectiveCrystal.h>
#include <DislocationDynamics/DislocationNetwork.h>
#include <DislocationDynamics/DislocationFieldBase.h>
#include <chrono>
#include <algorithm>
#include <iomanip>
using Clock=std::chrono::steady_clock;
int main(int argc,char** argv)
{
    try
    {
        if(argc<3 || argc>5) { throw std::invalid_argument("Usage: benchmarkCaiStressParticles tutorial frame [targets=1024] [repeats=9]"); }
        const std::size_t requested=argc>3?std::stoul(argv[3]):1024;
        const int repeats=argc>4?std::stoi(argv[4]):9;
        if(!requested || repeats<3) { throw std::invalid_argument("Positive target count and at least 3 repeats required"); }
        Kokkos::ScopeGuard guard;
        model::DefectiveCrystalBase<3> base(argv[1]);
        model::DefectiveCrystalConfigIO<3> config(base.simulationParameters.traitsIO.configFolder);
        config.read(std::stol(argv[2]));
        model::MicrostructureContainer<3> container(base);
        model::DislocationNetwork<3> network(container);
        network.initializeConfiguration(config,std::ofstream(),std::ofstream());
        network.updateGeometry();
        const double density=model::TextFileParser(base.simulationParameters.traitsIO.inputFilesFolder+"/DD.txt").readScalar<double>("quadPerLength",true);
        network.fieldEvaluator.rebuild(density,density,model::DislocationTargetPurpose::glide);
        const auto& buffers=network.fieldEvaluator.buffers();
        model::DislocationTargetArrays<3> targets;
        const auto count=std::min(requested,buffers.targets().size());
        if(!count || !buffers.sources().size()) { throw std::runtime_error("Benchmark requires sources and targets"); }
        for(std::size_t t=0;t<count;++t)
        {
            const auto p=t*buffers.targets().size()/count;
            targets.positions.append({buffers.targets().positions.values[0][p],buffers.targets().positions.values[1][p],buffers.targets().positions.values[2][p]});
        }
        model::KokkosDislocationFieldBuffers<3> device;
        device.uploadSources(buffers.sources()); device.uploadTargets(targets);
        const model::CaiStressParticleKernel kernel(base.poly.mu,base.poly.nu,model::DislocationFieldBase<3>::a);
        const auto run=[&]()
        {
            device.executionSpace().fence();
            const auto start=Clock::now();
            model::evaluateCaiParticles(device.executionSpace(),kernel,device.sourcePositions(),device.sourceTensors(),device.targetPositions(),device.stresses(),device.sourceCount(),device.targetCount());
            device.executionSpace().fence();
            return std::chrono::duration<double,std::milli>(Clock::now()-start).count();
        };
        model::DislocationStressArrays<3> result;
        run(); device.downloadStresses(result);
        double norm=0;
        for(const auto& component:result.stresses.values)
        {
            for(double value:component)
            {
                if(!std::isfinite(value)) { throw std::runtime_error("Nonfinite stress"); }
                norm+=value*value;
            }
        }
        std::vector<double> times;
        for(int i=0;i<repeats;++i) { times.push_back(run()); }
        std::sort(times.begin(),times.end());
        std::cout<<std::setprecision(8)<<"CAI_BENCHMARK frame="<<argv[2]<<" backend="<<Kokkos::DefaultExecutionSpace::name()
                 <<" concurrency="<<device.executionSpace().concurrency()<<" sources="<<device.sourceCount()<<" targets="<<count<<" available_targets="<<buffers.targets().size()
                 <<" repeats="<<repeats<<" six_median_ms="<<times[repeats/2]
                 <<" six_min_ms="<<times.front()<<" six_max_ms="<<times.back()
                 <<" stress_norm="<<std::sqrt(norm)<<std::endl;
    }
    catch(const std::exception& error) { std::cerr<<error.what()<<std::endl; return 1; }
}
