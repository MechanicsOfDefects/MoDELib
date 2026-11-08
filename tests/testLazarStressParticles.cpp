#include <DislocationDynamics/LazarStressParticleKernel.h>
#include "LazarStressParticleReference.h"
#include <random>
#include <iostream>
#include <algorithm>
#include <type_traits>
void check(bool b,const char* msg) { if(!b) { throw std::runtime_error(msg); } }
using model::LazarPart;
using K=model::LazarStressParticleKernel<LazarPart::full>;
struct Der { long double d1=0,d2=0,d3=0; };
Der expTerm(long double r,long double t,bool overR)
{
    const long double e=std::exp(-t*r);
    if(overR) { return {-e*(t/r+1/(r*r)),e*(t*t/r+2*t/(r*r)+2/(r*r*r)),-e*(t*t*t/r+3*t*t/(r*r)+6*t/(r*r*r)+6/(r*r*r*r))}; }
    return {-t*e,t*t*e,-t*t*t*e};
}
void add(Der& a,Der b,long double w) { a.d1+=w*b.d1; a.d2+=w*b.d2; a.d3+=w*b.d3; }
int main()
{
    try {
    Kokkos::ScopeGuard guard;
    double equivalenceError=0;
    const auto equivalence=[&](auto part)
    {
        constexpr auto which=decltype(part)::value;
        std::mt19937 gen(731);
        std::uniform_real_distribution<double> uniform(-1,1);
        for(double ell:{0.2,1.0,3.0}) { for(double scale:{4.0,8.0,20.0})
        {
            model::LazarStressParticleKernel<which> kernel(2.3,.28,ell,ell*scale);
            model::LazarStressParticleReference<which> reference(2.3,.28,ell,ell*scale);
            for(int sample=0;sample<500;++sample)
            {
                const double radius=ell*std::pow(10.0,-10+14.0*sample/499);
                double r[3]={uniform(gen),uniform(gen),uniform(gen)},q[9];
                const double norm=std::sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
                for(double& value:r) { value*=radius/norm; }
                for(double& value:q) { value=uniform(gen); }
                double old[6]={},now[6]={}; reference.addSymmetric(r,q,old); kernel.addSymmetric(r,q,now);
                double err=0,ref=0;
                for(int j=0;j<6;++j) { err+=(old[j]-now[j])*(old[j]-now[j]); ref+=old[j]*old[j]; }
                const double scaled=std::sqrt(err)/std::max(1/(ell*ell),std::sqrt(ref));
                equivalenceError=std::max(equivalenceError,scaled);
                check(scaled<2e-9,"Pre-optimization stress equivalence");
            }
        } }
    };
    equivalence(std::integral_constant<LazarPart,LazarPart::full>{});
    equivalence(std::integral_constant<LazarPart,LazarPart::shortRange>{});
    equivalence(std::integral_constant<LazarPart,LazarPart::longRange>{});
    std::cout<<"Maximum scaled pre-optimization difference="<<equivalenceError<<std::endl;
    double maxError=0;
    for(double ell:{0.2,1.0,3.0}) { for(double ratio:{4.0,8.0,20.0})
    {
        K full(1,0.3,ell,ell*ratio);
        model::LazarStressParticleKernel<LazarPart::shortRange> sh(1,0.3,ell,ell*ratio);
        model::LazarStressParticleKernel<LazarPart::longRange> lo(1,0.3,ell,ell*ratio);
        for(double y:{1e-10,0.02,0.099999,0.100001,0.3,1.0,4.0,20.0,100.0})
        {
            const double r=y*ell; double B,C,b,c,bl,cl;
            full.radial(r,B,C); sh.radial(r,b,c); lo.radial(r,bl,cl);
            check(std::abs(B-b-bl)*ell*ell<2e-9 && std::abs(C-c-cl)*ell*ell<2e-9,"Split reconstruction");
            if(y==1e-10) { check(std::abs(B*ell*ell+0.25)<1e-9 && std::abs(C*ell*ell-0.25)<1e-9,"Core limit"); continue; }
            const long double length=ell;
            Der f; f.d1=1; add(f,expTerm(r,0,true),2*length*length); add(f,expTerm(r,1/length,true),-2*length*length);
            const long double rb=f.d2/r-f.d1/(r*r),rc=f.d3-3*rb;
            maxError=std::max(maxError,double(std::max(std::abs(B-rb),std::abs(C-rc))*ell*ell));
            check(std::abs(B-rb)*ell*ell<2e-7 && std::abs(C-rc)*ell*ell<2e-7,"Full radial reference");
            if(ratio>=2)
            {
                const long double l=ell,h=l*ratio,D=h*h-l*l;
                Der s; add(s,expTerm(r,1/h,true),2*h*h*h*h*l*l/(D*D));
                add(s,expTerm(r,1/l,true),-2*h*h*h*h*l*l/(D*D)); add(s,expTerm(r,1/h,false),-h*h*h/D);
                const long double sb=s.d2/r-s.d1/(r*r),sc=s.d3-3*sb;
                check(std::abs(b-sb)*ell*ell<2e-7 && std::abs(c-sc)*ell*ell<2e-7,"Short radial reference");
            }
        }
    } }
    // Check oddness and the classical far-field limit of the six-component stress.
    const double r[3]={0.7,-0.4,1.1},Q[9]={0.2,0.4,-0.1,0.3,-0.8,0.7,0.9,0.2,-0.3};
    const double radius=std::sqrt(1.86),mu=1,nu=.3;
    double B,C,out[6]={}; K full(mu,nu,1,4); full.radial(radius,B,C); full.addSymmetric(r,Q,out);
    // Independent classical limit: Cai with vanishing core at large radius.
    const double far[3]={7000,-4000,11000}; double ca[6]={},la[6]={},neg[6]={};
    model::CaiStressParticleKernel cai(mu,nu,1e-5);
    cai.addSymmetric(far,Q,ca); full.addSymmetric(far,Q,la);
    const double minus[3]={-r[0],-r[1],-r[2]}; full.addSymmetric(minus,Q,neg);
    for(int i=0;i<6;++i)
    {
        check(std::abs(ca[i]-la[i])<1e-15,"Classical far limit");
        check(std::abs(neg[i]+out[i])<1e-14,"Odd stress kernel");
    }
    // Exercise the same device evaluator as the benchmark and the exact-coincidence convention.
    Kokkos::View<double**> src("src",1,3),tgt("tgt",2,3),tensor("Q",1,9),stress("stress",2,9);
    auto ts=Kokkos::create_mirror_view(tgt); auto qs=Kokkos::create_mirror_view(tensor);
    for(int i=0;i<3;++i) { ts(0,i)=r[i]; ts(1,i)=0; } for(int i=0;i<9;++i) { qs(0,i)=Q[i]; }
    Kokkos::deep_copy(src,0.0); Kokkos::deep_copy(tgt,ts); Kokkos::deep_copy(tensor,qs);
    const auto deviceCheck=[&](const auto& kernel)
    {
        double expected[6]={}; kernel.addSymmetric(r,Q,expected);
        model::evaluateStressParticles(Kokkos::DefaultExecutionSpace(),kernel,src,tensor,tgt,stress,1,2);
        auto hs=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),stress);
        int z=0; for(int i=0;i<3;++i) { for(int j=i;j<3;++j) { check(std::abs(hs(0,3*i+j)-expected[z++])<1e-14,"Device agreement"); } }
        for(int j=0;j<9;++j) { check(hs(1,j)==0,"Coincidence convention"); }
    };
    deviceCheck(full);
    deviceCheck(model::LazarStressParticleKernel<LazarPart::shortRange>(mu,nu,1,4));
    deviceCheck(model::LazarStressParticleKernel<LazarPart::longRange>(mu,nu,1,4));
    for(double h:{-1.0,0.0,0.5,1.0,1.000000001,2.0,std::nextafter(4.0,0.0)})
    {
        const auto rejects=[&](auto part)
        {
            bool rejected=false;
            try { model::LazarStressParticleKernel<decltype(part)::value> bad(mu,nu,1,h); }
            catch(const std::invalid_argument&) { rejected=true; }
            check(rejected,"Reject insufficient scale separation at construction");
        };
        rejects(std::integral_constant<LazarPart,LazarPart::full>{});
        rejects(std::integral_constant<LazarPart,LazarPart::shortRange>{});
        rejects(std::integral_constant<LazarPart,LazarPart::longRange>{});
    }
    std::cout<<"Lazar checks passed; maximum scaled radial error="<<maxError<<std::endl;
    } catch(const std::exception& e) { std::cerr<<e.what()<<std::endl; return 1; }
}
