#include <DislocationDynamics/LazarStressParticleKernel.h>
#include <DislocationDynamics/KokkosDislocationFieldBuffers.h>
#include <DefectiveCrystal/DefectiveCrystal.h>
#include <DislocationDynamics/DislocationNetwork.h>
#include <DislocationDynamics/DislocationFieldBase.h>
#include <FMM/KokkosShortRangeParticleNeighbors.h>
#include <chrono>
#include <algorithm>
#include <iomanip>
using Clock=std::chrono::steady_clock;
int main(int argc,char** argv)
{
    try
    {
        if(argc<3 || argc>5) { throw std::invalid_argument("Usage: benchmarkLazarNeighbors tutorial frame [targets=1024] [repeats=9]"); }
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
        const model::LazarStressParticleKernel<model::LazarPart::shortRange> kernel(base.poly.mu,base.poly.nu,ell,h);
        for(int periodic=0;periodic<2;++periodic)
        {
            const model::ShortRangeParticleNeighbors::Basis basis=periodic?base.mesh.periodicBasis():model::ShortRangeParticleNeighbors::Basis(3,0);
            if(periodic && basis.cols()==0) { continue; }
            const auto start=Clock::now();
            const model::ShortRangeParticleNeighbors list(buffers.sources().positions.values,targets.positions.values,rc,basis);
            const double build=std::chrono::duration<double,std::milli>(Clock::now()-start).count();
            const auto uploadStart=Clock::now();
            model::KokkosShortRangeParticleNeighbors<> neighbors(device.executionSpace(),list);
            const double upload=std::chrono::duration<double,std::milli>(Clock::now()-uploadStart).count();
            const auto run=[&]()
            {
                const auto t=Clock::now();
                neighbors.evaluate(device.executionSpace(),kernel,device.sourcePositions(),device.sourceTensors(),device.targetPositions(),device.stresses(),device.sourceCount(),device.targetCount());
                device.executionSpace().fence();
                return std::chrono::duration<double,std::milli>(Clock::now()-t).count();
            };
            run(); std::vector<double> times;
            for(int i=0;i<repeats;++i) { times.push_back(run()); }
            std::sort(times.begin(),times.end());
            model::DislocationStressArrays<3> result; device.downloadStresses(result);
            // Brute-force Kokkos reference uses the same translations, without cell culling.
            using ShiftView=Kokkos::View<double**,Kokkos::LayoutLeft>;
            ShiftView shifts("reference shifts",list.shifts.size(),3); auto host=Kokkos::create_mirror(shifts);
            for(std::size_t i=0;i<list.shifts.size();++i) { for(int j=0;j<3;++j) { host(i,j)=list.shifts[i][j]; } }
            Kokkos::deep_copy(shifts,host);
            const auto src=device.sourcePositions(),tensor=device.sourceTensors(),tgt=device.targetPositions();
            const auto output=device.stresses(); const std::size_t ns=device.sourceCount(),ni=list.shifts.size();
            const auto bruteStart=Clock::now();
            Kokkos::parallel_for("Brute periodic cutoff",Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(device.executionSpace(),0,count),
                KOKKOS_LAMBDA(const int t)
                {
                    double sum[6]={};
                    for(std::size_t im=0;im<ni;++im) { for(std::size_t s=0;s<ns;++s)
                    {
                        double r[3],q[9],d2=0;
                        for(int j=0;j<3;++j) { r[j]=(tgt(t,j)-src(s,j))-shifts(im,j); d2+=r[j]*r[j]; }
                        if(d2<=rc*rc)
                        {
                            for(int j=0;j<9;++j) { q[j]=tensor(s,j); }
                            kernel.addSymmetric(r,q,sum);
                        }
                    } }
                    int k=0; for(int i=0;i<3;++i) { for(int j=i;j<3;++j) { output(t,3*i+j)=output(t,3*j+i)=sum[k++]; } }
                });
            device.executionSpace().fence();
            const double brute=std::chrono::duration<double,std::milli>(Clock::now()-bruteStart).count();
            model::DislocationStressArrays<3> ref; device.downloadStresses(ref);
            double err=0,norm=0;
            for(int j=0;j<9;++j) { for(std::size_t t=0;t<count;++t)
            {
                const double d=result.stresses.values[j][t]-ref.stresses.values[j][t]; err+=d*d;
                norm+=ref.stresses.values[j][t]*ref.stresses.values[j][t];
            } }
            if(!std::isfinite(err)||std::sqrt(err/std::max(norm,1e-300))>1e-11) { throw std::runtime_error("Neighbor versus brute force stress mismatch"); }
            std::cout<<std::setprecision(10)<<"NEIGHBORS frame="<<argv[2]<<" periodic_rank="<<basis.cols()<<" sources="<<ns<<" targets="<<count
                     <<" images="<<ni<<" candidates="<<list.candidates<<" pairs="<<list.pairs.size()
                     <<" build_ms="<<build<<" upload_ms="<<upload<<" eval_median_ms="<<times[repeats/2]
                     <<" total_ms="<<build+upload+times[repeats/2]<<" brute_single_ms="<<brute
                     <<" relative_error="<<std::sqrt(err/std::max(norm,1e-300))<<std::endl;
        }
    }
    catch(const std::exception& error) { std::cerr<<error.what()<<std::endl; return 1; }
}
