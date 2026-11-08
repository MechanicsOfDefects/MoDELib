#include <DislocationDynamics/LazarStressParticleKernel.h>
#include <DislocationDynamics/KokkosDislocationFieldBuffers.h>
#include <DefectiveCrystal/DefectiveCrystal.h>
#include <DislocationDynamics/DislocationNetwork.h>
#include <DislocationDynamics/DislocationFieldBase.h>
#include <FMM/KokkosStressFMM.h>
#include <FMM/KokkosShortRangeParticleNeighbors.h>
#include <chrono>
#include <algorithm>
#include <iomanip>
using Clock=std::chrono::steady_clock;
int main(int argc,char** argv)
{
    try
    {
        if(argc<3 || argc>5) { throw std::invalid_argument("Usage: benchmarkLazarFMM tutorial frame [targets=1024] [repeats=9]"); }
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

        const double ell=model::DislocationFieldBase<3>::a,h=4*ell,rc=12*h;
        const model::LazarStressParticleKernel<model::LazarPart::longRange> kernel(base.poly.mu,base.poly.nu,ell,h);
        const auto directStart=Clock::now();
        model::evaluateStressParticles(device.executionSpace(),kernel,device.sourcePositions(),device.sourceTensors(),device.targetPositions(),device.stresses(),device.sourceCount(),device.targetCount());
        device.executionSpace().fence();
        const double directMs=std::chrono::duration<double,std::milli>(Clock::now()-directStart).count();
        model::DislocationStressArrays<3> reference; device.downloadStresses(reference);
        const model::LazarStressParticleKernel<model::LazarPart::full> fullKernel(base.poly.mu,base.poly.nu,ell,h);
        model::evaluateStressParticles(device.executionSpace(),fullKernel,device.sourcePositions(),device.sourceTensors(),device.targetPositions(),device.stresses(),device.sourceCount(),device.targetCount());
        model::DislocationStressArrays<3> full; device.downloadStresses(full);
        const model::ShortRangeParticleNeighbors list(buffers.sources().positions.values,targets.positions.values,rc);
        model::KokkosShortRangeParticleNeighbors<> neighbors(device.executionSpace(),list);
        const model::LazarStressParticleKernel<model::LazarPart::shortRange> shortKernel(base.poly.mu,base.poly.nu,ell,h);
        neighbors.evaluate(device.executionSpace(),shortKernel,device.sourcePositions(),device.sourceTensors(),device.targetPositions(),device.stresses(),device.sourceCount(),device.targetCount());
        model::DislocationStressArrays<3> shortResult; device.downloadStresses(shortResult);
        for(int order:{3,4})
        {
        for(int mode:{0,1,2})
        {
            const auto start=Clock::now();
            model::KokkosStressFMM<> fmm(device.executionSpace(),network.fieldEvaluator.tree(),buffers.sources().positions.values,targets.positions.values,order);
            fmm.setTranslationCacheBudget(mode==0?0:256*1024*1024);
            fmm.setInterpolateNear(mode==2);
            fmm.setProfiling(true);
            const double setup=std::chrono::duration<double,std::milli>(Clock::now()-start).count();
            const auto run=[&]()
            {
                const auto t=Clock::now();
                fmm.evaluate(device.executionSpace(),kernel,device.sourcePositions(),device.sourceTensors(),device.targetPositions(),device.stresses(),device.sourceCount(),device.targetCount());
                return std::chrono::duration<double,std::milli>(Clock::now()-t).count();
            };
            const double firstMs=run(); std::vector<double> times;
            for(int i=0;i<repeats;++i) { times.push_back(run()); }
            std::sort(times.begin(),times.end());
            model::DislocationStressArrays<3> result; device.downloadStresses(result);
            double err=0,norm=0,totalError=0,totalNorm=0,maxError=0;
            for(std::size_t t=0;t<count;++t)
            {
                double pt=0;
                for(int c=0;c<9;++c)
                {
                    const double r=reference.stresses.values[c][t],v=result.stresses.values[c][t],f=full.stresses.values[c][t];
                    const double d=v-r,dt=v+shortResult.stresses.values[c][t]-f;
                    pt+=d*d; norm+=r*r; totalError+=dt*dt; totalNorm+=f*f;
                }
                err+=pt; maxError=std::max(maxError,std::sqrt(pt));
            }
            if(!std::isfinite(err)) { throw std::runtime_error("Nonfinite FMM error"); }
            std::cout<<std::setprecision(10)<<"FMM_BENCH frame="<<argv[2]<<" depth="<<network.fieldEvaluator.tree().depth()<<" order="<<order
                     <<" sources="<<device.sourceCount()<<" targets="<<count<<" far_box_pairs="<<fmm.farBoxPairs()<<" near_box_pairs="<<fmm.nearBoxPairs()
                     <<" mode="<<mode<<" cache_bytes="<<fmm.translationCacheBytes()<<" first_ms="<<firstMs<<" setup_ms="<<setup<<" median_ms="<<times[repeats/2]<<" direct_single_ms="<<directMs
                     <<" upward_last_ms="<<fmm.lastTimings()[1]<<" downward_last_ms="<<fmm.lastTimings()[2]<<" near_last_ms="<<fmm.lastTimings()[3]
                     <<" relative_long="<<std::sqrt(err/norm)<<" relative_total="<<std::sqrt(totalError/totalNorm)
                     <<" max_error_over_long_rms="<<maxError/std::sqrt(norm/count)<<std::endl;
        }
        }
    }
    catch(const std::exception& error) { std::cerr<<error.what()<<std::endl; return 1; }
}
