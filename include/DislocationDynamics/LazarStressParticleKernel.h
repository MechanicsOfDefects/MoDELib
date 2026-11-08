/* This file is part of MoDELib, distributed under GNU GPL v2. */
#ifndef model_LazarStressParticleKernel_H_
#define model_LazarStressParticleKernel_H_
#include <DislocationDynamics/CaiStressParticleKernel.h>
#include <array>
namespace model
{
    enum class LazarPart { full, shortRange, longRange };
    template<LazarPart Part> struct LazarStressParticleKernel
    {
        // Validate scale separation once; no h/ell branches in particle evaluation.
        double ell,h,ratio,first,second,invEll,invEll2,q2,q3,screenA,screenB,seriesB[19],seriesC[19];
        LazarStressParticleKernel(double mu,double nu,double length,double screen)
            :ell(length),h(screen),ratio(length/screen),first(-mu/(8*3.14159265358979323846)),second(mu/(4*3.14159265358979323846*(1-nu)))
        {
            if(!std::isfinite(mu)||mu<=0||!std::isfinite(nu)||nu<=-1||nu>=0.5||
               !std::isfinite(ell)||ell<=0||!std::isfinite(h)||h<=0||ratio>0.25||!std::isfinite(1/(ell*ell))||ratio==0)
            { throw std::invalid_argument("Lazar kernel requires finite mu>0, -1<nu<0.5 and h>=4*ell>0"); }
            invEll=1/ell; invEll2=invEll*invEll;
            q2=ratio*ratio; q3=q2*ratio;
            const double d=1-q2;
            screenA=2/(d*d); screenB=-1/(ratio*d);
            // Dimensionless potential A/ell in y=r/ell. Precompute its series.
            double ex[22],psi[22]; ex[0]=1; psi[0]=0.5;
            for(int n=1;n<22;++n) { ex[n]=ex[n-1]*(-ratio)/n; psi[n]=psi[n-1]*(-(1-ratio))/(n+2); }
            double full=2;
            for(int n=0;n<22;++n)
            {
                if(n) { full*=-1.0/(n+1); }
                double sh=ex[n];
                for(int j=0;j<n;++j) { sh+=2*ratio*ex[n-1-j]*psi[j]; }
                sh*=-1/(ratio*(1+ratio)*(1+ratio));
                const double f=full+(n==1?1:0);
                if(n>=3)
                {
                    const double value=Part==LazarPart::full?f:(Part==LazarPart::shortRange?sh:f-sh);
                    seriesB[n-3]=n*(n-2)*value*invEll2;
                    seriesC[n-3]=n*(n-2)*(n-4)*value*invEll2;
                }
            }
        }
        // Radial third-derivative tensor coefficients.
        KOKKOS_INLINE_FUNCTION void radial(double r,double& B,double& C) const
        {
            const double y=r*invEll;
            if(y<0.1)
            {
                B=0; C=0;
                for(int n=18;n>=0;--n)
                {
                    B=B*y+seriesB[n];
                    C=C*y+seriesC[n];
                }
            }
            else
            {
                const double inv=1/y,i2=inv*inv,i4=i2*i2,e=Kokkos::exp(-y);
                const double u=i2*(1+inv*(3+3*inv));
                const double v=inv*(1+inv*(6+inv*(15+15*inv)));
                double bf=0,cf=0;
                if constexpr(Part!=LazarPart::shortRange)
                {
                    bf=-i2+6*i4-2*e*u;
                    cf=3*i2-30*i4+2*e*v;
                }
                if constexpr(Part==LazarPart::full) { B=bf; C=cf; }
                else
                {
                    const double eh=Kokkos::exp(-ratio*y);
                    const double uh=i2*(q2+inv*(3*ratio+3*inv));
                    const double vh=inv*(q3+inv*(6*q2+inv*(15*ratio+15*inv)));
                    B=screenA*(eh*uh-e*u)+screenB*eh*inv*(q2+ratio*inv);
                    C=screenA*(e*v-eh*vh)-screenB*eh*(q3+inv*(3*q2+3*ratio*inv));
                    if constexpr(Part==LazarPart::longRange) { B=bf-B; C=cf-C; }
                }
                B*=invEll2; C*=invEll2;
            }
        }

        // Geometry/material-only translation data; reusable when source strengths change.
        std::array<double,5> translationKey() const
        { return {ell,h,first,second,static_cast<double>(Part)}; }
        KOKKOS_INLINE_FUNCTION void prepareTranslation(const double* r,double* data) const
        {
            const double radius=Kokkos::sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
            for(int i=0;i<7;++i) { data[i]=0; }
            if(radius==0) { return; }
            const double invRadius=1/radius;
            for(int i=0;i<3;++i) { data[i]=r[i]*invRadius; }
            double B,C; radial(radius,B,C);
            data[3]=-first*(C+5*B); data[4]=second*B;
            data[5]=second*C; data[6]=-second*(4*B+C);
        }
        KOKKOS_INLINE_FUNCTION void addSymmetric(const double* r,const double* Q,double* out) const
        {
            double data[7]; prepareTranslation(r,data); applyTranslation(data,Q,out);
        }
        KOKKOS_INLINE_FUNCTION void applyTranslation(const double* data,const double* Q,double* out) const
        {
            const double x=data[0],y=data[1],z=data[2];
            const double cx=Q[5]-Q[7],cy=Q[6]-Q[2],cz=Q[1]-Q[3];
            const double dot=cx*x+cy*y+cz*z;
            // Four shared radial/material combinations, then six scalar contractions.
            const double f=data[3],g=data[4];
            const double e=data[5]*dot,d=data[6]*dot;
            out[0]+=2*f*(y*Q[6]-z*Q[3])+2*g*cx*x+e*x*x+d;
            out[1]+=f*(y*Q[7]-z*Q[4]+z*Q[0]-x*Q[6])+g*(cx*y+cy*x)+e*x*y;
            out[2]+=f*(y*Q[8]-z*Q[5]+x*Q[3]-y*Q[0])+g*(cx*z+cz*x)+e*x*z;
            out[3]+=2*f*(z*Q[1]-x*Q[7])+2*g*cy*y+e*y*y+d;
            out[4]+=f*(z*Q[2]-x*Q[8]+x*Q[4]-y*Q[1])+g*(cy*z+cz*y)+e*y*z;
            out[5]+=2*f*(x*Q[5]-y*Q[2])+2*g*cz*z+e*z*z+d;
        }
    };
}
#endif
