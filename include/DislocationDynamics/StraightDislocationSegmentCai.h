/* MoDELib: analytical unscreened Cai segment stress, displacement, energy, concentration and flux. */
#ifndef model_StraightDislocationSegmentCai_H_
#define model_StraightDislocationSegmentCai_H_
#include <Eigen/Dense>
#include <DislocationDynamics/StraightDislocationSegmentBase.h>
#include <cmath>
#include <numbers>
namespace model
{
    template<typename Scalar>
    Eigen::Matrix<Scalar,6,1> caiSegmentEndpointStress(
        const Eigen::Matrix<Scalar,3,1>& r,const Eigen::Matrix<Scalar,3,1>& b,
        const Eigen::Matrix<Scalar,3,1>& t,const Eigen::Matrix<Scalar,3,1>& bCrossT,
        const Scalar a2,const Scalar oneMinusNu)
    {
        using Vector=Eigen::Matrix<Scalar,3,1>;
        const Scalar R2=r.squaredNorm()+a2;
        const Scalar R=std::sqrt(R2);
        const Vector Y=r+R*t;
        const Scalar Yt=Y.dot(t),Y2=Y.squaredNorm()+a2,f=2/Y2;
        const Vector bY=b.cross(Y);
        const Scalar bYt=bY.dot(t),g=f*bYt;
        const Vector v=(f*oneMinusNu*(1+a2/Y2))*bY
                      +(f*oneMinusNu*Scalar(0.5)*a2/R2)*b.cross(r)-(g/Yt)*r;
        const Vector w=-f*bCrossT;
        const Scalar diagonal=-g*(1+2*a2/Y2+a2/R2);
        const Scalar rr=-2*g*(R+Yt)/(R*Y2),tt=-g*R/Yt;
        Eigen::Matrix<Scalar,6,1> result;
        int k=0;
        for(int i=0;i<3;++i)
        {
            for(int j=i;j<3;++j)
            {
                result[k++]=t[i]*v[j]+t[j]*v[i]+Y[i]*w[j]+Y[j]*w[i]
                    +rr*r[i]*r[j]+tt*t[i]*t[j]+(i==j?diagonal:Scalar(0));
            }
        }
        return result;
    }
    template<typename Scalar>
    Eigen::Matrix<Scalar,3,3> caiSegmentStress(
        const Eigen::Matrix<Scalar,3,1>& p0,const Eigen::Matrix<Scalar,3,1>& p1,
        const Eigen::Matrix<Scalar,3,1>& x,const Eigen::Matrix<Scalar,3,1>& b,
        const Eigen::Matrix<Scalar,3,1>& t,const Eigen::Matrix<Scalar,3,1>& bCrossT,
        const Scalar a2,const Scalar oneMinusNu,const Scalar prefactor)
    {
        const Eigen::Matrix<Scalar,6,1> value=prefactor*(
            caiSegmentEndpointStress<Scalar>(p1-x,b,t,bCrossT,a2,oneMinusNu)
           -caiSegmentEndpointStress<Scalar>(p0-x,b,t,bCrossT,a2,oneMinusNu));
        Eigen::Matrix<Scalar,3,3> result;
        int k=0;
        for(int i=0;i<3;++i)
        {
            for(int j=i;j<3;++j)
            {
                result(i,j)=result(j,i)=value[k++];
            }
        }
        return result;
    }
    template<typename Scalar>
    Eigen::Matrix<Scalar,3,1> caiSegmentEndpointDisplacement(
        const Eigen::Matrix<Scalar,3,1>& r,const Eigen::Matrix<Scalar,3,1>& t,
        const Eigen::Matrix<Scalar,3,1>& bCrossT,const Scalar a2,const Scalar oneMinusNu)
    {
        using Vector=Eigen::Matrix<Scalar,3,1>;
        const Scalar Ra(std::sqrt(r.squaredNorm()+a2));
        const Vector Ya(r+Ra*t);
        const Scalar Yat(Ya.dot(t));
        return -(2.0-0.5/oneMinusNu+(2.0-1.0/oneMinusNu)*std::log(Yat)-a2/Ra/Yat)/8.0/std::numbers::pi_v<Scalar>*bCrossT
        /*  */ +bCrossT.dot(r)/Yat/Ra/oneMinusNu/8.0/std::numbers::pi_v<Scalar>*Ya;
    }

    // Line-integral contribution only; excludes the solid-angle displacement.
    template<typename Scalar>
    Eigen::Matrix<Scalar,3,1> caiSegmentDisplacement(
        const Eigen::Matrix<Scalar,3,1>& p0,const Eigen::Matrix<Scalar,3,1>& p1,
        const Eigen::Matrix<Scalar,3,1>& x,const Eigen::Matrix<Scalar,3,1>& t,
        const Eigen::Matrix<Scalar,3,1>& bCrossT,const Scalar a2,const Scalar oneMinusNu)
    {
        return caiSegmentEndpointDisplacement<Scalar>(p1-x,t,bCrossT,a2,oneMinusNu)
             - caiSegmentEndpointDisplacement<Scalar>(p0-x,t,bCrossT,a2,oneMinusNu);
    }

    template<typename Scalar>
    Scalar caiSegmentEndpointInteractionEnergy(
        const Eigen::Matrix<Scalar,3,1>& z,const Eigen::Matrix<Scalar,3,1>& t,
        const Eigen::Matrix<Scalar,3,1>& b,const Eigen::Matrix<Scalar,3,1>& tA,
        const Eigen::Matrix<Scalar,3,1>& bA,const Scalar a2,const Scalar oneMinusNu,const Scalar nu)
    {
        using Vector=Eigen::Matrix<Scalar,3,1>;
        const Scalar za(std::sqrt(z.squaredNorm()+a2));
        const Vector Ya(z+za*t);
        const Scalar Yat(Ya.dot(t));
        const Scalar logYat(std::log(Yat));
        const Scalar zaYat(za*Yat);

        const Scalar bAtA(bA.dot(tA));
        const Scalar bt(b.dot(t));
        const Scalar bAt(bA.dot(t));
        const Scalar btA(b.dot(tA));
        const Scalar bAb(bA.dot(b));
        const Scalar tAt(tA.dot(t));
        const Scalar zt(z.dot(t));
        const Scalar zb(z.dot(b));
        const Scalar zbA(z.dot(bA));
        const Scalar energyInf((0.5*oneMinusNu*bAtA*bt+nu*bAt*btA)*(2.0+2.0*logYat-a2/zaYat)
                           /*  */-bAb*tAt*(1.5+logYat-a2/zaYat)
                           /*  */-bAt*bt*tAt*(0.5+logYat+zt/Yat)
                           /*  */+tAt/Yat*(zbA*zb/za+bAt*zb+bt*zbA));
        return energyInf;
    }

    // Interaction with the local line element (xA,tA,bA), integrated over this segment.
    template<typename Scalar>
    Scalar caiSegmentInteractionEnergy(
        const Eigen::Matrix<Scalar,3,1>& p0,const Eigen::Matrix<Scalar,3,1>& p1,
        const Eigen::Matrix<Scalar,3,1>& xA,const Eigen::Matrix<Scalar,3,1>& t,
        const Eigen::Matrix<Scalar,3,1>& b,const Eigen::Matrix<Scalar,3,1>& tA,
        const Eigen::Matrix<Scalar,3,1>& bA,const Scalar a2,const Scalar oneMinusNu,
        const Scalar nu,const Scalar prefactor)
    {
        return -Scalar(0.5)*prefactor*(
            caiSegmentEndpointInteractionEnergy<Scalar>(p1-xA,t,b,tA,bA,a2,oneMinusNu,nu)
           -caiSegmentEndpointInteractionEnergy<Scalar>(p0-xA,t,b,tA,bA,a2,oneMinusNu,nu));
    }

    // Endpoint integrals of (1-u,u)/sqrt(u*u+(b/a)*u+c/a), u in [0,1].
    // The physical Green-function normalization includes the remaining 1/sqrt(a).
    template<typename Scalar>
    Eigen::Matrix<Scalar,2,1> caiSegmentConcentrationShapeIntegrals(Scalar a,Scalar b,Scalar c)
    {
        const Scalar ba=b/a,ca=c/a;
        const Scalar R1=std::sqrt(1+ba+ca),R0=std::sqrt(ca);
        const Scalar logarithm=std::log((2*R1+2+ba)/(2*R0+ba));
        return {(1+Scalar(.5)*ba)*logarithm-R1+R0,-Scalar(.5)*ba*logarithm+R1-R0};
    }

    // -D grad of the line integral with source(u)=source0+(source1-source0)*u.
    // Material/species normalization, line length, atomic volume and image sums
    // are applied by the caller. a=chord dot inverse(D) chord; core is in that metric.
    template<typename Scalar>
    Eigen::Matrix<Scalar,3,1> caiSegmentFluxIntegral(
        const Eigen::Matrix<Scalar,3,1>& r,const Eigen::Matrix<Scalar,3,1>& chord,
        const Eigen::Matrix<Scalar,3,3>& inverse,Scalar a,Scalar core,Scalar source0,Scalar source1)
    {
        using Vector=Eigen::Matrix<Scalar,3,1>;
        const Scalar u0=chord.dot(inverse*r)/a;
        const Vector perpendicular=r-u0*chord;
        const Scalar h2=(perpendicular.dot(inverse*perpendicular)+core)/a;
        const Scalar h=std::sqrt(h2);
        const Scalar t0=-u0;
        const Scalar t1=1-u0;
        const Scalar R0=std::hypot(t0,h);
        const Scalar R1=std::hypot(t1,h);
        Scalar j0;
        if(t0*t1>0)
        {
            j0=(t1+t0)/(R0*R1*(t1*R0+t0*R1));
        }
        else
        {
            j0=(t1/R1-t0/R0)/h2;
        }
        const Scalar j1=(t1+t0)/((R1+R0)*R0*R1);
        const Scalar j2=std::asinh(t1/h)-std::asinh(t0/h)-h2*j0;
        const Scalar delta=source1-source0;
        const Scalar center=source0+delta*u0;
        const Vector integral=(perpendicular*(center*j0+delta*j1)
                                  -chord*(center*j1+delta*j2))/(a*std::sqrt(a));
        return integral;
    }

    inline Eigen::Vector2d caiSegmentConcentrationIntegrals(const Eigen::Matrix3d& diffusion,double core,
        const Eigen::Vector3d& p0,const Eigen::Vector3d& p1,const Eigen::Vector3d& x,double begin=0,double end=1)
    {
        const Eigen::Vector3d chord=p1-p0,delta=(end-begin)*chord,r=x-p0-begin*chord;
        if(delta.squaredNorm()==0)
        {
            return Eigen::Vector2d::Zero();
        }
        const Eigen::Matrix3d inv=diffusion.inverse();
        const double det=diffusion.determinant(),a=delta.dot(inv*delta);
        const auto v=caiSegmentConcentrationShapeIntegrals(a,-2*r.dot(inv*delta),r.dot(inv*r)+core*core/std::cbrt(det));
        const double factor=delta.norm()/(4*std::numbers::pi*std::sqrt(a*det));
        return factor*Eigen::Vector2d((1-begin)*v[0]+(1-end)*v[1],begin*v[0]+end*v[1]);
    }

    class CaiSegmentConcentration final : public SegmentConcentrationKernelBase
    {
        Eigen::Matrix3d diffusion_;
        Eigen::Matrix<long double,3,3> inverse_;
        double core_;
        long double determinant_;
    public:
        CaiSegmentConcentration(const Eigen::Matrix3d& diffusion,double core)
            :diffusion_(diffusion),inverse_(diffusion.inverse().cast<long double>()),core_(core),determinant_(diffusion.determinant()) {}
        Eigen::Vector2d concentrationIntegrals(const Eigen::Vector3d& a,const Eigen::Vector3d& b,const Eigen::Vector3d& x,double begin=0,double end=1) const override
        {
            return caiSegmentConcentrationIntegrals(diffusion_,core_,a,b,x,begin,end);
        }
        Eigen::Vector3d flux(const Eigen::Vector3d& p0,const Eigen::Vector3d& p1,const Eigen::Vector3d& x,double q0,double q1,double begin=0,double end=1) const override
        {
            using Vec=Eigen::Matrix<long double,3,1>;
            const Vec chord=(p1-p0).cast<long double>();
            const Vec delta=(end-begin)*chord,r=(x-p0).cast<long double>()-begin*chord;
            if(delta.squaredNorm()==0)
            {
                return Eigen::Vector3d::Zero();
            }
            const long double a=delta.dot(inverse_*delta),core=core_*core_/std::cbrt(determinant_);
            const auto integral=caiSegmentFluxIntegral<long double>(r,delta,inverse_,a,core,q0+(q1-q0)*begin,q0+(q1-q0)*end);
            return (delta.norm()/(4*std::numbers::pi_v<long double>*std::sqrt(determinant_))*integral).cast<double>();
        }
    };
    class StraightDislocationSegmentCai final : public StraightDislocationSegmentBase
    {
        SegmentCoreParameters parameters_;
    public:
        explicit StraightDislocationSegmentCai(const SegmentCoreParameters& parameters):parameters_(parameters) {}
        const char* name() const override { return "Cai"; }
        Eigen::Matrix3d stress(const Eigen::Vector3d& x,const Segment& s) const override
        {
            const Eigen::Vector3d chord=s.p1-s.p0;
            if(chord.squaredNorm()==0)
            {
                return Eigen::Matrix3d::Zero();
            }
            const Eigen::Vector3d t=chord.normalized(),bc=s.b.cross(t);
            return caiSegmentStress<double>(s.p0,s.p1,x,s.b,t,bc,parameters_.coreWidth*parameters_.coreWidth,
                1-parameters_.nu,parameters_.mu/(4*std::numbers::pi*(1-parameters_.nu)));
        }
        Eigen::Matrix3d nearStress(const Eigen::Vector3d& x,const Segment& s) const override { return stress(x,s); }
        Eigen::Matrix3d shortStress(const Eigen::Vector3d&,const Segment&) const override { throw std::logic_error("Cai Ewald split is not implemented"); }
        Eigen::Matrix3d longStress(const Eigen::Vector3d&,const Segment&) const override { throw std::logic_error("Cai Ewald split is not implemented"); }
        Eigen::Vector3d displacement(const Eigen::Vector3d& x,const Segment& s) const override
        {
            const Eigen::Vector3d t=(s.p1-s.p0).normalized(),bc=s.b.cross(t);
            return caiSegmentDisplacement<double>(s.p0,s.p1,x,t,bc,parameters_.coreWidth*parameters_.coreWidth,1-parameters_.nu);
        }
        double elasticInteractionEnergy(const Eigen::Vector3d& x,const Eigen::Vector3d& tA,const Eigen::Vector3d& bA,const Segment& s) const override
        {
            const Eigen::Vector3d t=(s.p1-s.p0).normalized();
            return caiSegmentInteractionEnergy<double>(s.p0,s.p1,x,t,s.b,tA,bA,parameters_.coreWidth*parameters_.coreWidth,
                1-parameters_.nu,parameters_.nu,parameters_.mu/(4*std::numbers::pi*(1-parameters_.nu)));
        }
        ConcentrationPtr concentration(const Eigen::Matrix3d& D,bool=false,SegmentKernelPart part=SegmentKernelPart::full) const override
        {
            if(part!=SegmentKernelPart::full)
            {
                throw std::logic_error("Cai Ewald split is not implemented");
            }
            return std::make_shared<CaiSegmentConcentration>(D,parameters_.coreWidth);
        }
        SegmentMomentIntegrator::StressMoments stressMoments(const SegmentMomentIntegrator& helper,const Segment& segment,const Eigen::Vector3d& center,const Eigen::Vector3d& width) const override
        {
            return helper.stress(segment,center,width);
        }
        SegmentMomentIntegrator::ConcentrationMoments concentrationMoments(const SegmentMomentIntegrator& helper,const Eigen::Vector3d& p0,const Eigen::Vector3d& p1,const Eigen::Vector3d& center,const Eigen::Vector3d& width,double begin,double end,const std::array<double,2>& weights) const override
        {
            return helper.concentration(p0,p1,center,width,begin,end,weights);
        }
    };

}
#endif
