#include <DislocationDynamics/LazarStressParticleKernel.h>
#include <DislocationDynamics/KokkosDislocationFieldBuffers.h>
#include <DefectiveCrystal/DefectiveCrystal.h>
#include <DislocationDynamics/DislocationNetwork.h>
#include <DislocationDynamics/DislocationFieldBase.h>
#include <chrono>
#include <algorithm>
#include <iomanip>

// Benchmark-only cutoff wrapper; production force evaluation is unchanged.
struct CutoffKernel
{
    model::LazarStressParticleKernel<model::LazarPart::shortRange> kernel;
    double radiusSquared;
    KOKKOS_INLINE_FUNCTION void addSymmetric(const double* r,const double* q,double* out) const
    {
        if(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]<=radiusSquared) { kernel.addSymmetric(r,q,out); }
    }
};
int main(int argc,char** argv)
{
    try
    {
        if(argc<3 || argc>4) { throw std::invalid_argument("Usage: benchmarkLazarCutoff tutorial frame [targets=1024]"); }
        const std::size_t requested=argc>3?std::stoul(argv[3]):1024;
        if(!requested) { throw std::invalid_argument("Positive target count required"); }
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

        const double ell=model::DislocationFieldBase<3>::a,h=4*ell;
        using Short=model::LazarStressParticleKernel<model::LazarPart::shortRange>;
        const Short shortKernel(base.poly.mu,base.poly.nu,ell,h);
        const auto evaluate=[&](const auto& kernel)
        {
            model::evaluateStressParticles(device.executionSpace(),kernel,device.sourcePositions(),device.sourceTensors(),device.targetPositions(),device.stresses(),device.sourceCount(),device.targetCount());
            device.executionSpace().fence();
            model::DislocationStressArrays<3> result;
            device.downloadStresses(result);
            return result;
        };
        const auto reference=evaluate(shortKernel);
        const auto full=evaluate(model::LazarStressParticleKernel<model::LazarPart::full>(base.poly.mu,base.poly.nu,ell,h));
        const auto longField=evaluate(model::LazarStressParticleKernel<model::LazarPart::longRange>(base.poly.mu,base.poly.nu,ell,h));
        double shortNorm2=0,fullNorm2=0,splitError2=0;
        for(int c=0;c<9;++c) { for(std::size_t t=0;t<count;++t)
        {
            const double a=reference.stresses.values[c][t],b=full.stresses.values[c][t];
            shortNorm2+=a*a; fullNorm2+=b*b;
            const double d=a+longField.stresses.values[c][t]-b; splitError2+=d*d;
        } }
        if(!(shortNorm2>0 && fullNorm2>0) || std::sqrt(splitError2/fullNorm2)>1e-9)
        { throw std::runtime_error("Invalid reference or full/short/long reconstruction"); }
        std::cout<<std::setprecision(10)<<"CUTOFF_SETUP frame="<<argv[2]<<" sources="<<device.sourceCount()<<" targets="<<count
                 <<" ell="<<ell<<" h="<<h<<" full_rms="<<std::sqrt(fullNorm2/count)
                 <<" split_relative_error="<<std::sqrt(splitError2/fullNorm2)<<std::endl;
        for(double multiple:{2.,3.,4.,5.,6.,8.,10.,12.,16.,20.})
        {
            const double cutoff=multiple*h;
            const auto truncated=evaluate(CutoffKernel{shortKernel,cutoff*cutoff});
            const auto src=device.sourcePositions(),tgt=device.targetPositions();
            const std::size_t ns=device.sourceCount();
            unsigned long long retained=0;
            Kokkos::parallel_reduce("Count retained cutoff pairs",Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(device.executionSpace(),0,count),
                KOKKOS_LAMBDA(const int t,unsigned long long& sum)
                {
                    for(std::size_t p=0;p<ns;++p)
                    {
                        double d2=0; for(int j=0;j<3;++j) { const double d=tgt(t,j)-src(p,j); d2+=d*d; }
                        if(d2<=cutoff*cutoff) { ++sum; }
                    }
                },retained);
            device.executionSpace().fence();
            double error2=0,maxError=0;
            for(std::size_t t=0;t<count;++t)
            {
                double pointError2=0;
                for(int c=0;c<9;++c)
                {
                    const double d=truncated.stresses.values[c][t]-reference.stresses.values[c][t];
                    if(!std::isfinite(d)) { throw std::runtime_error("Nonfinite cutoff error"); }
                    pointError2+=d*d;
                }
                error2+=pointError2; maxError=std::max(maxError,std::sqrt(pointError2));
            }
            std::cout<<"CUTOFF frame="<<argv[2]<<" rc_over_h="<<multiple<<" retained="<<retained
                     <<" retained_fraction="<<double(retained)/(double(ns)*count)
                     <<" relative_short="<<std::sqrt(error2/shortNorm2)<<" relative_full="<<std::sqrt(error2/fullNorm2)
                     <<" rms_absolute="<<std::sqrt(error2/count)<<" max_absolute="<<maxError
                     <<" max_over_full_rms="<<maxError/std::sqrt(fullNorm2/count)<<std::endl;
        }
    }
    catch(const std::exception& error) { std::cerr<<error.what()<<std::endl; return 1; }
}
