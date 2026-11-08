/* This file is part of MoDELib, distributed under GNU GPL v2. */
#ifndef model_CaiStressParticleKernel_H_
#define model_CaiStressParticleKernel_H_
#include <cfloat>
#include <Kokkos_Core.hpp>
#include <cmath>
#include <stdexcept>
namespace model
{
    // Unscreened isotropic Cai kernel, x-y convention, row-major Q=b*t*dL.
    // Core radius is positive. No cutoff, images, or special self subtraction.
    struct CaiStressParticleKernel
    {
        double a2, first, second;
        CaiStressParticleKernel(double mu,double nu,double a)
            :a2(a*a),first(-mu/(8.0*3.14159265358979323846)),
             second(mu/(4.0*3.14159265358979323846*(1.0-nu)))
        {
            if(!std::isfinite(mu) || mu<=0 || !std::isfinite(nu) || nu<=-1 || nu>=0.5 ||
               !std::isfinite(a) || a<=0 || !std::isfinite(a2) || a2==0)
            { throw std::invalid_argument("Invalid isotropic Cai kernel parameters"); }
        }
        struct Terms
        {
            double c[3], v[9], f,g,d,e;
        };
        KOKKOS_INLINE_FUNCTION Terms terms(const double* r,const double* q) const
        {
            Terms z;
            const double inv=1.0/Kokkos::sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]+a2);
            const double inv3=inv*inv*inv,inv5=inv3*inv*inv;
            z.c[0]=q[5]-q[7]; z.c[1]=q[6]-q[2]; z.c[2]=q[1]-q[3];
            const double s=z.c[0]*r[0]+z.c[1]*r[1]+z.c[2]*r[2];
            z.f=first*(2*inv3+3*a2*inv5);
            z.g=-second*inv3;
            z.d=second*s*(inv3+3*a2*inv5);
            z.e=3*second*s*inv5;
            for(int j=0;j<3;++j)
            {
                z.v[j]=r[1]*q[6+j]-r[2]*q[3+j];
                z.v[3+j]=r[2]*q[j]-r[0]*q[6+j];
                z.v[6+j]=r[0]*q[3+j]-r[1]*q[j];
            }
            return z;
        }
        // Six independent sigma entries: xx, xy, xz, yy, yz, zz.
        KOKKOS_INLINE_FUNCTION void addSymmetric(const double* r,const double* q,double* s) const
        {
            const auto z=terms(r,q);
            int k=0;
            for(int i=0;i<3;++i)
            {
                for(int j=i;j<3;++j)
                {
                    s[k++]+=z.f*(z.v[3*i+j]+z.v[3*j+i])+z.g*(z.c[i]*r[j]+z.c[j]*r[i])+
                            z.e*r[i]*r[j]+(i==j?z.d:0.0);
                }
            }
        }
    };
    // Six-component accumulation with full symmetric tensor output.
    template<class ExecutionSpace,class Kernel,class PositionView,class TensorView,class OutputView>
    void evaluateStressParticles(const ExecutionSpace& execution,const Kernel& kernel,
                              PositionView sources,TensorView tensors,PositionView targets,
                              OutputView stresses,std::size_t sourceCount,std::size_t targetCount)
    {
        if(sourceCount>sources.extent(0) || sourceCount>tensors.extent(0) || targetCount>targets.extent(0) ||
           targetCount>stresses.extent(0) || sources.extent(1)!=3 || targets.extent(1)!=3 ||
           tensors.extent(1)!=9 || stresses.extent(1)!=9)
        { throw std::invalid_argument("Invalid stress particle view extents"); }
        Kokkos::parallel_for("Stress six components",
            Kokkos::RangePolicy<ExecutionSpace,Kokkos::IndexType<std::size_t>>(execution,0,targetCount),
            KOKKOS_LAMBDA(const std::size_t t)
        {
            double sum[6]={};
            const double x[3]={targets(t,0),targets(t,1),targets(t,2)};
            for(std::size_t p=0;p<sourceCount;++p)
            {
                double r[3],q[9];
                for(int a=0;a<3;++a) { r[a]=x[a]-sources(p,a); }
                for(int a=0;a<9;++a) { q[a]=tensors(p,a); }
                kernel.addSymmetric(r,q,sum);
            }
            int k=0;
            for(int i=0;i<3;++i)
            {
                for(int j=i;j<3;++j)
                {
                    const double value=sum[k++];
                    stresses(t,3*i+j)=value;
                    stresses(t,3*j+i)=value;
                }
            }
        });
    }
    template<class ExecutionSpace,class PositionView,class TensorView,class OutputView>
    void evaluateCaiParticles(const ExecutionSpace& execution,const CaiStressParticleKernel& kernel,
                             PositionView sources,TensorView tensors,PositionView targets,OutputView stresses,
                             std::size_t sourceCount,std::size_t targetCount)
    {
        evaluateStressParticles(execution,kernel,sources,tensors,targets,stresses,sourceCount,targetCount);
    }

}
#endif
