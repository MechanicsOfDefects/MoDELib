// Frozen pre-optimization implementation for numerical equivalence tests.
/* This file is part of MoDELib, distributed under GNU GPL v2. */
#ifndef model_LazarStressParticleReference_H_
#define model_LazarStressParticleReference_H_
#include <DislocationDynamics/CaiStressParticleKernel.h>
namespace model
{
    template<LazarPart Part> struct LazarStressParticleReference
    {
        // Validate scale separation once; no h/ell branches in particle evaluation.
        double ell,h,ratio,first,second,coeff[22];
        LazarStressParticleReference(double mu,double nu,double length,double screen)
            :ell(length),h(screen),ratio(length/screen),first(-mu/(8*3.14159265358979323846)),second(mu/(4*3.14159265358979323846*(1-nu)))
        {
            if(!std::isfinite(mu)||mu<=0||!std::isfinite(nu)||nu<=-1||nu>=0.5||
               !std::isfinite(ell)||ell<=0||!std::isfinite(h)||h<=0||ratio>0.25||!std::isfinite(1/(ell*ell))||ratio==0)
            { throw std::invalid_argument("Lazar kernel requires finite mu>0, -1<nu<0.5 and h>=4*ell>0"); }
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
                coeff[n]=Part==LazarPart::full?f:(Part==LazarPart::shortRange?sh:f-sh);
            }
        }
        // Radial third-derivative tensor coefficients.
        KOKKOS_INLINE_FUNCTION void radial(double r,double& B,double& C) const
        {
            const double y=r/ell;
            if(y<0.1)
            {
                B=0; C=0;
                for(int n=21;n>=3;--n)
                {
                    B=B*y+n*(n-2)*coeff[n];
                    C=C*y+n*(n-2)*(n-4)*coeff[n];
                }
            }
            else
            {
                const double inv=1/y,i2=inv*inv,i3=i2*inv,i4=i2*i2,e=Kokkos::exp(-y);
                const double bf=-i2+6*i4-2*e*(i2+3*i3+3*i4);
                const double cf=3*i2-30*i4+2*e*(inv+6*i2+15*i3+15*i4);
                if constexpr(Part==LazarPart::full) { B=bf; C=cf; }
                else
                {
                    const double q=ratio,q2=q*q,q3=q2*q,eh=Kokkos::exp(-q*y),d=1-q2;
                    const double a=2/(d*d),b=-1/(q*d);
                    B=a*(eh*(q2*i2+3*q*i3+3*i4)-e*(i2+3*i3+3*i4))+b*eh*(q2*inv+q*i2);
                    C=a*(-eh*(q3*inv+6*q2*i2+15*q*i3+15*i4)+e*(inv+6*i2+15*i3+15*i4))-b*eh*(q3+3*q2*inv+3*q*i2);
                    if constexpr(Part==LazarPart::longRange) { B=bf-B; C=cf-C; }
                }
            }
            B/=ell*ell; C/=ell*ell;
        }
        KOKKOS_INLINE_FUNCTION void addSymmetric(const double* r,const double* Q,double* out) const
        {
            const double radius=Kokkos::sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
            // Odd-kernel principal-value convention at exact coincidence.
            // Helmholtz stress has a direction-dependent finite limit, not a continuous zero limit.
            if(radius==0) { return; }
            double n[3]={r[0]/radius,r[1]/radius,r[2]/radius},B,C;
            radial(radius,B,C);
            const double c[3]={Q[5]-Q[7],Q[6]-Q[2],Q[1]-Q[3]};
            const double s=c[0]*n[0]+c[1]*n[1]+c[2]*n[2];
            double v[9];
            for(int j=0;j<3;++j)
            {
                v[j]=n[1]*Q[6+j]-n[2]*Q[3+j];
                v[3+j]=n[2]*Q[j]-n[0]*Q[6+j];
                v[6+j]=n[0]*Q[3+j]-n[1]*Q[j];
            }
            int k=0;
            for(int i=0;i<3;++i) { for(int j=i;j<3;++j)
            {
                out[k++]+=-first*(C+5*B)*(v[3*i+j]+v[3*j+i])+second*(B*(c[i]*n[j]+c[j]*n[i])+C*s*n[i]*n[j]-(i==j?(4*B+C)*s:0));
            } }
        }
    };
}
#endif
