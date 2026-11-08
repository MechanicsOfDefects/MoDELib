#include <FMM/KokkosFMMInteraction.h>
#include <FMM/HelmholtzDiffusionKernel.h>
#include <FMM/HelmholtzStressKernelFamily.h>
#include <iostream>
#include <random>
#include <functional>
using namespace model;
void require(bool v,const char* m) { if(!v) throw std::runtime_error(m); }
template<class F> void rejects(F f) { bool caught=false; try { f(); } catch(const std::exception&) { caught=true; } require(caught,"Expected rejection"); }
struct CustomKernel
{
    KOKKOS_INLINE_FUNCTION void add(const double* r,const double* q,double* out) const
    { out[0]+=q[0];out[1]+=q[1]*r[0];out[2]+=(q[0]+q[1])*r[1]*r[1]; }
};
int main()
{
    try
    {
        Kokkos::ScopeGuard guard;
        Kokkos::DefaultExecutionSpace ex;
        using Scalar=KokkosFMMInteraction<HelmholtzDiffusionFamily>;
        using View=Scalar::View;
        Eigen::Matrix3d R=Eigen::AngleAxisd(.73,Eigen::Vector3d(1,2,3).normalized()).toRotationMatrix();
        Eigen::Matrix3d D=R*Eigen::Vector3d(.3,1.,3.).asDiagonal()*R.transpose();
        const double a=.02,h=.12;
        HelmholtzDiffusionFamily family(D,a,h);
        auto full=family.fullKernel(); auto sh=family.shortKernel(); auto lo=family.longKernel();
        const double d=std::cbrt(D.determinant()),pi=std::acos(-1.);
        const double zero[3]={};
        require(std::abs(lo.value(zero)*4*pi*d*(a+h)-1)<1e-13,"Long self value");
        require(std::abs(full.value(zero)*4*pi*d*a-1)<1e-13,"Full self value");
        for(int k=-12;k<=3;++k)
        {
            const double scale=std::pow(10.,k),r[3]={scale,.2*scale,-.1*scale};
            require(std::abs((sh.value(r)+lo.value(r))/full.value(r)-1)<2e-13,"Split reconstruction");
            double opposite[3]={-r[0],-r[1],-r[2]};
            require(lo.value(r)==lo.value(opposite),"Kernel reciprocity");
        }
        rejects([&]{ HelmholtzDiffusionFamily bad(-D,a,h); });
        rejects([&]{ HelmholtzDiffusionFamily bad(D,a,a); });
        // Covariance under rotation with matching anisotropic metric.
        HelmholtzDiffusionFamily diagonal(Eigen::Vector3d(.3,1.,3.).asDiagonal(),a,h);
        Eigen::Vector3d x(.3,.4,.5),rx=R*x;
        require(std::abs(lo.value(rx.data())/diagonal.longKernel().value(x.data())-1)<1e-13,"Rotation covariance");
        FMMManager registry;
        auto scalar=registry.registerKernel("CD/vacancies",family);
        auto stress=registry.registerKernel("DD/stress",HelmholtzStressKernelFamily(1,.3,a,h));
        require(registry.kernelCount()==2,"Registry count");
        rejects([&]{registry.registerKernel("CD/vacancies",family);});
        rejects([&]{registry.kernel<HelmholtzStressKernelFamily>("CD/vacancies");});
        Scalar::Positions src,tgt;
        std::mt19937 gen(3401); std::uniform_real_distribution<double> u(0,1);
        for(int n=0;n<48;++n) for(int c=0;c<3;++c) src[c].push_back(u(gen));
        for(int n=0;n<12;++n) for(int c=0;c<3;++c) tgt[c].push_back(u(gen));
        for(int c=0;c<3;++c) tgt[c].push_back(src[c][0]);
        View q("strengths",48,1),out("concentration",13,1);
        auto qh=Kokkos::create_mirror(q);
        for(int n=0;n<48;++n) qh(n,0)=.5+u(gen);
        Kokkos::deep_copy(q,qh);
        auto direct=[&](auto kernel,double cutoff)
        {
            std::vector<double> result(13);
            for(int t=0;t<13;++t) for(int s=0;s<48;++s)
            {
                double r[3];for(int c=0;c<3;++c) r[c]=tgt[c][t]-src[c][s];
                if(kernel.distance(r)<=cutoff) result[t]+=kernel.value(r)*qh(s,0);
            }
            return result;
        };
        auto error=[&](const std::vector<double>& ref)
        {
            auto oh=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),out);
            double e=0,n=0;for(int t=0;t<13;++t){e+=std::pow(oh(t,0)-ref[t],2);n+=ref[t]*ref[t];}
            return std::sqrt(e/n);
        };
        auto tree0=std::make_shared<BoxTree<3>>(Eigen::Vector3d::Zero(),Eigen::Vector3d::Ones(),0);
        FMMInteractionOptions opt; opt.boxShortRange=false;opt.interpolateScreenedNear=false; opt.shortCutoff=.2; opt.interpolationOrder=4;
        auto interaction=registry.registerInteraction<Scalar>("climb",scalar,ex,tree0,src,tgt,opt);
        interaction->evaluateShort(q,out);
        require(error(direct(sh,.2))<1e-13,"Metric cutoff short evaluation");
        interaction->evaluateLong(q,out);
        require(error(direct(lo,1e100))<1e-13,"Depth zero long evaluation");
        auto expected=direct(lo,1e100),sr=direct(sh,.2);
        for(int i=0;i<13;++i) expected[i]+=sr[i];
        interaction->evaluate(q,out); require(error(expected)<1e-13,"Combined split evaluation");
        auto tree=std::make_shared<BoxTree<3>>(Eigen::Vector3d::Zero(),Eigen::Vector3d::Ones(),2);
        interaction->rebuildGeometry(tree,src,tgt);
        interaction->evaluateLong(q,out);
        const double fmmError=error(direct(lo,1e100));
        std::cout<<"Scalar FMM relative error "<<fmmError<<std::endl;
        require(fmmError<.003,"Scalar FMM accuracy");
        const auto cache=interaction->translationCacheBytes(); require(cache>0,"Scalar translation cache unused");
        for(int i=0;i<48;++i) qh(i,0)*=-2;
        Kokkos::deep_copy(q,qh);interaction->evaluateLong(q,out);
        require(error(direct(lo,1e100))<.003,"Updated strengths not used");
        require(cache==interaction->translationCacheBytes(),"Cache not reused");
        rejects([&]{auto bad=opt;bad.boundary=FMMBoundaryMode::periodic;registry.registerInteraction<Scalar>("periodic",scalar,ex,tree,src,tgt,bad);});
        FMMManager another;rejects([&]{another.registerInteraction<Scalar>("foreign",scalar,ex,tree,src,tgt,opt);});
        rejects([&]{View wrong("bad shape",48,2);interaction->evaluate(wrong,out);});
        // Stress family shares traversal but retains its six independent tensor outputs.
        using Stress=KokkosFMMInteraction<HelmholtzStressKernelFamily>;
        auto si=registry.registerInteraction<Stress>("stress",stress,ex,tree0,src,tgt,opt);
        View tensor("Q",48,9),sigma("stress6",13,6);
        auto th=Kokkos::create_mirror(tensor);for(int i=0;i<48;++i)for(int c=0;c<9;++c)th(i,c)=u(gen)-.5;
        Kokkos::deep_copy(tensor,th); si->evaluateLong(tensor,sigma);
        auto result=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),sigma);
        auto sk=stress.family().longKernel();
        for(int t=0;t<13;++t)
        {
            double sum[6]={};for(int s=0;s<48;++s){double r[3],qq[9];for(int c=0;c<3;++c)r[c]=tgt[c][t]-src[c][s];for(int c=0;c<9;++c)qq[c]=th(s,c);sk.add(r,qq,sum);}
            for(int c=0;c<6;++c)require(std::abs(result(t,c)-sum[c])<1e-11*(1+std::abs(sum[c])),"Registered stress differs from direct");
        }
        // A third-party 2-input/3-output kernel needs neither stress methods nor a cache.
        View customQ("custom Q",48,2),customOut("custom output",13,3);
        auto ch=Kokkos::create_mirror(customQ);
        for(int i=0;i<48;++i){ch(i,0)=qh(i,0);ch(i,1)=.5;}
        Kokkos::deep_copy(customQ,ch);
        View sp("source xyz",48,3),tp("target xyz",13,3);
        auto sph=Kokkos::create_mirror(sp),tph=Kokkos::create_mirror(tp);
        for(int c=0;c<3;++c){for(int i=0;i<48;++i)sph(i,c)=src[c][i];for(int i=0;i<13;++i)tph(i,c)=tgt[c][i];}
        Kokkos::deep_copy(sp,sph);Kokkos::deep_copy(tp,tph);
        KokkosParticleFMM<2,3> custom(ex,*tree,src,tgt,4);
        custom.evaluate(ex,CustomKernel{},sp,customQ,tp,customOut,48,13);
        auto customResult=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),customOut);
        for(int t=0;t<13;++t){double sum[3]={};for(int i=0;i<48;++i){double r[3],qq[2]={ch(i,0),ch(i,1)};for(int c=0;c<3;++c)r[c]=tgt[c][t]-src[c][i];CustomKernel{}.add(r,qq,sum);}for(int c=0;c<3;++c)require(std::abs(customResult(t,c)-sum[c])<1e-11*(1+std::abs(sum[c])),"Generic component dimensions");}
        // All-box interpolation: finite self interactions and neighboring boxes,
        // with no target-side particle sum (even when source count is below p^3).
        custom.setInterpolateNearBoxes(true);
        custom.evaluate(ex,CustomKernel{},sp,customQ,tp,customOut,48,13);
        auto smoothCustom=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),customOut);
        for(int t=0;t<13;++t)for(int c=0;c<3;++c)
            require(std::abs(smoothCustom(t,c)-customResult(t,c))<1.e-10,"Smooth neighbor translations changed polynomial field");
        auto smoothFamily=std::make_shared<HelmholtzDiffusionFamily>(D,a,4.0);
        // Cache self/neighbor translations under the same budget as ordinary M2L.
        // Depth zero has no far boxes: any allocated operator must be the self operator.
        for(auto geometry:{tree0,tree})
        {
            KokkosParticleFMM<1,1> cached(ex,*geometry,src,tgt,4);
            cached.setInterpolateNearBoxes(true);
            cached.setTranslationCacheBudget(0);
            cached.evaluate(ex,smoothFamily->longKernel(),sp,q,tp,out,48,13);
            const auto reference=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),out);
            for(std::size_t budget:{std::size_t(64*64*sizeof(double)),std::size_t(256*1024*1024)})
            {
                cached.setTranslationCacheBudget(budget);
                cached.evaluate(ex,smoothFamily->longKernel(),sp,q,tp,out,48,13);
                require(cached.translationCacheBytes()>0 && cached.translationCacheBytes()<=budget,"Near translation cache budget");
                const auto values=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),out);
                for(int i=0;i<13;++i)
                {
                    require(std::abs(values(i,0)-reference(i,0))<1e-12*(1+std::abs(reference(i,0))),"Cached neighbor translations changed result");
                }
                const auto bytes=cached.translationCacheBytes();
                Kokkos::deep_copy(q,qh);
                cached.evaluate(ex,smoothFamily->longKernel(),sp,q,tp,out,48,13);
                require(cached.translationCacheBytes()==bytes,"Near translation cache reuse");
            }
            // Changing the kernel must invalidate the cached coefficients.
            const HelmholtzDiffusionFamily changed(D,a,2.0);
            cached.evaluate(ex,changed.longKernel(),sp,q,tp,out,48,13);
            const auto changedCached=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),out);
            cached.setTranslationCacheBudget(0);
            cached.evaluate(ex,changed.longKernel(),sp,q,tp,out,48,13);
            const auto changedDirect=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),out);
            for(int i=0;i<13;++i)
            {
                require(std::abs(changedCached(i,0)-changedDirect(i,0))<1e-12*(1+std::abs(changedDirect(i,0))),"Near translation cache invalidation");
            }
        }
        auto smoothOptions=opt;smoothOptions.interpolateScreenedNear=true;smoothOptions.interpolationOrder=6;
        for(auto geometry : {tree0,tree})
        {
            Scalar smooth(smoothFamily,ex,geometry,src,tgt,smoothOptions);
            smooth.evaluateLong(q,out);
            const double smoothError=error(direct(smoothFamily->longKernel(),1.e100));
            std::cout<<"Smooth diffusion all-box relative error "<<smoothError<<std::endl;
            require(smoothError<.003,"Smooth diffusion all-box accuracy");
            smooth.evaluate(q,out);
            auto combined=direct(smoothFamily->longKernel(),1.e100),local=direct(smoothFamily->shortKernel(),opt.shortCutoff);
            for(int t=0;t<13;++t)combined[t]+=local[t];
            require(error(combined)<.003,"Smooth split double counted near interactions");
        }
        // Production box stencil: compare to an independent all-pairs cell test.
        // A tiny legacy cutoff deliberately proves no radial filter remains.
        auto boxOpt=opt;boxOpt.boxShortRange=true;boxOpt.interpolateScreenedNear=true;boxOpt.shortCutoff=1.e-30;
        Scalar boxes(smoothFamily,ex,tree,src,tgt,boxOpt);
        std::vector<double> boxRef(13);std::size_t pairCount=0;
        for(int t=0;t<13;++t)for(int s=0;s<48;++s)
        {
            bool adjacent=true;double r[3];
            for(int c=0;c<3;++c){adjacent=adjacent && std::abs(int(tgt[c][t]*4)-int(src[c][s]*4))<=1;r[c]=tgt[c][t]-src[c][s];}
            if(adjacent){boxRef[t]+=smoothFamily->shortKernel().value(r)*qh(s,0);++pairCount;}
        }
        boxes.evaluateShort(q,out);
        require(error(boxRef)<1.e-13,"Box short field differs from independent cell stencil");
        require(boxes.shortPairCount()==pairCount && pairCount<48*13,"Box pair count or locality");
        // Corner and outside-domain edge leaves: exactly one visit per source.
        Scalar::Positions edgeSrc={{{0.,.49,.51,1.},{0.,.49,0.,1.},{0.,.49,0.,1.}}};
        Scalar::Positions edgeTgt={{{0.,1.},{0.,1.},{0.,1.}}};
        LeafBoxNeighbors edgeList(*tree,edgeSrc,edgeTgt);
        std::vector<int> seen(4);edgeList.forEach(0,[&](std::size_t s){++seen[s];});
        require(seen==std::vector<int>({1,1,0,0}),"Corner stencil must include diagonal neighbor and exclude second neighbor");
        seen.assign(4,0);edgeList.forEach(1,[&](std::size_t s){++seen[s];});
        require(seen==std::vector<int>({0,0,0,1}),"Upper root boundary stencil");
        // Split kernels: shifted direct short work and shifted long P2P/M2L.
        auto shiftedTree=std::make_shared<BoxTree<3>>(Eigen::Vector3d::Zero(),Eigen::Vector3d::Ones(),2);
        FMMInteractionOptions shiftedOptions;shiftedOptions.interpolateScreenedNear=false;
        shiftedOptions.localImageShifts={Eigen::Vector3d::Zero(),Eigen::Vector3d::UnitX(),-Eigen::Vector3d::UnitX(),Eigen::Vector3d::UnitY()};
        std::vector<double> shiftedShort(13,0),shiftedLong(13,0);
        for(int t=0;t<13;++t) { for(int s=0;s<48;++s)
        {
            const Eigen::Vector3d x(tgt[0][t],tgt[1][t],tgt[2][t]),y(src[0][s],src[1][s],src[2][s]);
            const auto tb=shiftedTree->locate(x,2),sb=shiftedTree->locate(y,2);
            for(const auto& shift:shiftedOptions.localImageShifts)
            {
                const Eigen::Vector3d r=x-y-shift;
                shiftedLong[t]+=lo.value(r.data())*qh(s,0);
                if(!shiftedTree->admissible(tb->center()-sb->center()-shift,2))
                { shiftedShort[t]+=sh.value(r.data())*qh(s,0); }
            }
        } }
        Scalar shifted(std::make_shared<const HelmholtzDiffusionFamily>(family),ex,shiftedTree,src,tgt,shiftedOptions);
        shifted.evaluateShort(q,out);require(error(shiftedShort)<1.e-12,"Shifted short kernel mismatch");
        shifted.evaluateLong(q,out);require(error(shiftedLong)<.003,"Shifted long kernel mismatch");
        shiftedOptions.interpolateScreenedNear=true;
        Scalar shiftedCached(std::make_shared<const HelmholtzDiffusionFamily>(family),ex,shiftedTree,src,tgt,shiftedOptions);
        shiftedCached.evaluateLong(q,out);
        const auto savedShifted=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),out);
        std::vector<double> shiftedCachedValues(13);
        for(int i=0;i<13;++i) { shiftedCachedValues[i]=savedShifted(i,0); }
        shiftedOptions.translationCacheBytes=0;
        Scalar shiftedUncached(std::make_shared<const HelmholtzDiffusionFamily>(family),ex,shiftedTree,src,tgt,shiftedOptions);
        shiftedUncached.evaluateLong(q,out);
        require(error(shiftedCachedValues)<1.e-12,"Shifted self/neighbor cache differs from uncached translations");
        tree->setDepth(1);
        rejects([&]{interaction->evaluateLong(q,out);});
        // Empty snapshots and lifetime independent of the manager's name entry.
        Scalar::Positions empty;View eq("empty q",0,1),eo("empty out",0,1);
        interaction->rebuildGeometry(tree0,empty,empty);interaction->evaluate(eq,eo);
        registry.unregisterInteraction("climb");interaction->evaluate(eq,eo);
        std::cout<<"FMM registration tests passed\n";
    }
    catch(const std::exception& e) {std::cerr<<e.what()<<std::endl;return 1;}
}
