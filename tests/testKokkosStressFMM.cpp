#include <DislocationDynamics/LazarStressParticleKernel.h>
#include <DislocationDynamics/KokkosDislocationFieldBuffers.h>
#include <FMM/KokkosStressFMM.h>
#include <iostream>
#include <random>
void require(bool x,const char* m) { if(!x) { throw std::runtime_error(m); } }
struct ConstantKernel
{
    KOKKOS_INLINE_FUNCTION void addSymmetric(const double*,const double* q,double* out) const
    { for(int i=0;i<6;++i) { out[i]+=q[i]; } }
};
struct PolynomialKernel
{
    KOKKOS_INLINE_FUNCTION void addSymmetric(const double* r,const double* q,double* out) const
    {
        const double w=1+r[0]+r[1]*r[1]+r[2]*r[2]*r[2];
        for(int i=0;i<6;++i) { out[i]+=w*q[i]; }
    }
};
int main()
{
    try
    {
        Kokkos::ScopeGuard guard;
        std::mt19937 rng(831); std::uniform_real_distribution<double> u(0.,1.);
        model::DislocationSourceArrays<3> src; model::DislocationTargetArrays<3> tgt;
        for(int i=0;i<96;++i)
        {
            src.positions.append({u(rng),u(rng),u(rng)});
            std::array<double,9> q; for(auto& x:q) { x=u(rng)-.5; } src.weightedBurgersTensors.append(q);
        }
        // Include domain endpoints and a coincident source/target.
        tgt.positions.append({0,0,0}); tgt.positions.append({1,1,1});
        tgt.positions.append({src.positions.values[0][0],src.positions.values[1][0],src.positions.values[2][0]});
        for(int i=0;i<21;++i) { tgt.positions.append({u(rng),u(rng),u(rng)}); }
        model::KokkosDislocationFieldBuffers<3> data; data.uploadSources(src); data.uploadTargets(tgt);
        const auto ex=data.executionSpace();
        const model::LazarStressParticleKernel<model::LazarPart::longRange> kernel(1,.3,.01,.08);
        const auto direct=[&](const auto& k)
        {
            model::evaluateStressParticles(ex,k,data.sourcePositions(),data.sourceTensors(),data.targetPositions(),data.stresses(),src.size(),tgt.size());
            model::DislocationStressArrays<3> r; data.downloadStresses(r); return r;
        };
        const auto ref=direct(kernel),constantRef=direct(ConstantKernel{});
        const auto error=[&](const auto& reference)
        {
            model::DislocationStressArrays<3> actual; data.downloadStresses(actual); double err=0,norm=0;
            for(int c=0;c<9;++c) { for(std::size_t t=0;t<tgt.size();++t)
            { const double d=actual.stresses.values[c][t]-reference.stresses.values[c][t]; err+=d*d; norm+=reference.stresses.values[c][t]*reference.stresses.values[c][t]; } }
            require(std::isfinite(err),"Nonfinite FMM stress"); return std::sqrt(err/norm);
        };
        double first=0,last=0;
        for(int depth:{0,2,3})
        {
            model::BoxTree<3> tree(Eigen::Vector3d::Zero(),Eigen::Vector3d::Ones(),depth);
            for(int p:{2,4,6})
            {
                if(depth==0 && p!=2) { continue; }
                if(depth==3 && p!=4) { continue; }
                model::KokkosStressFMM<> fmm(ex,tree,src.positions.values,tgt.positions.values,p);
                fmm.evaluate(ex,ConstantKernel{},data.sourcePositions(),data.sourceTensors(),data.targetPositions(),data.stresses(),src.size(),tgt.size());
                require(error(constantRef)<1e-12,"FMM partition or transfer conservation");
                fmm.evaluate(ex,kernel,data.sourcePositions(),data.sourceTensors(),data.targetPositions(),data.stresses(),src.size(),tgt.size());
                const double e=error(ref);
                std::cout<<"FMM depth="<<depth<<" order="<<p<<" relative="<<e<<" far_pairs="<<fmm.farBoxPairs()<<std::endl;
                if(depth==0) { require(e<1e-13,"Depth zero must be direct"); }
                if(depth==2 && p==2) { first=e; }
                if(depth==2 && p==6) { last=e; require(e<2e-3,"FMM order six accuracy"); }
                if(depth==3) { require(e<1e-2,"FMM multilevel accuracy"); }
                // Reusing the plan must reset moments and locals, not accumulate results.
                fmm.evaluate(ex,kernel,data.sourcePositions(),data.sourceTensors(),data.targetPositions(),data.stresses(),src.size(),tgt.size());
                require(std::abs(error(ref)-e)<1e-14,"FMM repeated evaluation reset");
            }
        }
        require(last<first*.1,"FMM order convergence");
        // Cached coefficients must agree with on-the-fly translations, including changed material/core.
        {
            model::BoxTree<3> cacheTree(Eigen::Vector3d::Zero(),Eigen::Vector3d::Ones(),2);
            model::KokkosStressFMM<> cached(ex,cacheTree,src.positions.values,tgt.positions.values,3);
            for(double screen:{.08,.16})
            {
                const model::LazarStressParticleKernel<model::LazarPart::longRange> changed(2,.27,.01,screen);
                cached.setTranslationCacheBudget(0);
                cached.evaluate(ex,changed,data.sourcePositions(),data.sourceTensors(),data.targetPositions(),data.stresses(),src.size(),tgt.size());
                model::DislocationStressArrays<3> baseline; data.downloadStresses(baseline);
                cached.setTranslationCacheBudget(256*1024*1024);
                cached.evaluate(ex,changed,data.sourcePositions(),data.sourceTensors(),data.targetPositions(),data.stresses(),src.size(),tgt.size());
                require(cached.translationCacheBytes()>0,"Translation cache populated");
                require(error(baseline)<1e-12,"Cached and uncached translations differ");
                cached.setTranslationCacheBudget(27*27*7*sizeof(double));
                cached.evaluate(ex,changed,data.sourcePositions(),data.sourceTensors(),data.targetPositions(),data.stresses(),src.size(),tgt.size());
                require(cached.translationCacheBytes()==27*27*7*sizeof(double),"Cache budget is respected");
                require(error(baseline)<1e-12,"Partially cached translations differ");
            }
            // Change kernel parameters without resetting the cache explicitly.
            cached.evaluate(ex,kernel,data.sourcePositions(),data.sourceTensors(),data.targetPositions(),data.stresses(),src.size(),tgt.size());
            model::DislocationStressArrays<3> changedCached; data.downloadStresses(changedCached);
            cached.setTranslationCacheBudget(0);
            cached.evaluate(ex,kernel,data.sourcePositions(),data.sourceTensors(),data.targetPositions(),data.stresses(),src.size(),tgt.size());
            require(error(changedCached)<1e-12,"Kernel changes invalidate translation cache");
        }
        // Depth zero isolates near interpolation. Check exact polynomial reproduction and
        // Lazar convergence; convergence alone is not an accuracy guarantee.
        {
            const model::LazarStressParticleKernel<model::LazarPart::longRange> smooth(1,.3,.01,2.);
            const auto smoothRef=direct(smooth),polynomialRef=direct(PolynomialKernel{});
            model::BoxTree<3> nearTree(Eigen::Vector3d::Zero(),Eigen::Vector3d::Ones(),0);
            double coarse=0;
            for(int p:{2,3,4})
            {
                model::KokkosStressFMM<> near(ex,nearTree,src.positions.values,tgt.positions.values,p);
                near.setInterpolateNear(true);
                near.evaluate(ex,smooth,data.sourcePositions(),data.sourceTensors(),data.targetPositions(),data.stresses(),src.size(),tgt.size());
                const double e=error(smoothRef);
                std::cout<<"Near interpolation order="<<p<<" relative="<<e<<std::endl;
                if(p==2) { coarse=e; }
                if(p==4) { require(e<coarse*.5,"Near interpolation order convergence");
                    near.evaluate(ex,PolynomialKernel{},data.sourcePositions(),data.sourceTensors(),data.targetPositions(),data.stresses(),src.size(),tgt.size());
                    require(error(polynomialRef)<1e-12,"Near interpolation cubic reproduction"); }
                near.evaluate(ex,ConstantKernel{},data.sourcePositions(),data.sourceTensors(),data.targetPositions(),data.stresses(),src.size(),tgt.size());
                require(error(constantRef)<1e-12,"Near interpolation conserves constant field");
                near.setInterpolateNear(false);
                near.evaluate(ex,smooth,data.sourcePositions(),data.sourceTensors(),data.targetPositions(),data.stresses(),src.size(),tgt.size());
                require(error(smoothRef)<1e-13,"Near direct fallback");
            }
        }

        model::BoxTree<3> tree(Eigen::Vector3d::Zero(),Eigen::Vector3d::Ones(),1);
        model::KokkosStressFMM<> stale(ex,tree,src.positions.values,tgt.positions.values,2); tree.setDepth(2);
        bool rejected=false;
        try { stale.evaluate(ex,kernel,data.sourcePositions(),data.sourceTensors(),data.targetPositions(),data.stresses(),src.size(),tgt.size()); }
        catch(const std::logic_error&) { rejected=true; } require(rejected,"FMM rejects stale tree");
        // Empty source set produces zero stress, rather than retaining old local fields.
        model::DislocationSourceArrays<3> empty;
        model::KokkosDislocationFieldBuffers<3> zero; zero.uploadSources(empty); zero.uploadTargets(tgt);
        model::KokkosStressFMM<> zeroPlan(ex,tree,empty.positions.values,tgt.positions.values,2);
        zeroPlan.evaluate(ex,kernel,zero.sourcePositions(),zero.sourceTensors(),zero.targetPositions(),zero.stresses(),0,tgt.size());
        model::DislocationStressArrays<3> zeros; zero.downloadStresses(zeros);
        for(const auto& component:zeros.stresses.values) { for(double v:component) { require(v==0,"Empty source result"); } }
        model::KokkosStressFMM<> reuse(ex,tree,src.positions.values,tgt.positions.values,3);
        reuse.evaluate(ex,kernel,data.sourcePositions(),data.sourceTensors(),data.targetPositions(),data.stresses(),src.size(),tgt.size());
        model::DislocationStressArrays<3> before; data.downloadStresses(before);
        for(auto& component:src.weightedBurgersTensors.values) { for(auto& v:component) { v*=2; } }
        data.uploadSources(src);
        reuse.evaluate(ex,kernel,data.sourcePositions(),data.sourceTensors(),data.targetPositions(),data.stresses(),src.size(),tgt.size());
        model::DislocationStressArrays<3> after; data.downloadStresses(after);
        for(int c=0;c<9;++c) { for(std::size_t t=0;t<tgt.size();++t)
        { require(std::abs(after.stresses.values[c][t]-2*before.stresses.values[c][t])<1e-12,"Updated source strengths"); } }
        rejected=false;
        try { model::KokkosStressFMM<> invalid(ex,tree,src.positions.values,tgt.positions.values,1000000); }
        catch(const std::invalid_argument&) { rejected=true; } require(rejected,"Invalid interpolation order");
        std::cout<<"FMM tests passed"<<std::endl;
    }
    catch(const std::exception& e) { std::cerr<<e.what()<<std::endl; return 1; }
}
