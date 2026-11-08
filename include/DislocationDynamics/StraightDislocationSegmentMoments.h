#ifndef model_StraightDislocationSegmentMoments_H_
#define model_StraightDislocationSegmentMoments_H_
#include <DislocationDynamics/HelmholtzSegmentGeometry.h>
#include <vector>
#include <array>
#include <cmath>
#include <numbers>
#include <stdexcept>
namespace model
{
struct SegmentGaussRule {
 std::vector<double>x,w;
 explicit SegmentGaussRule(int n):x(n),w(n){
  for(int i=0;i<(n+1)/2;++i){
   double z=std::cos(std::numbers::pi*(i+.75)/(n+.5)),dp=0;
   for(int it=0;it<30;++it){double p0=1,p1=z;for(int k=2;k<=n;++k){double p=((2*k-1)*z*p1-(k-1)*p0)/k;p0=p1;p1=p;}
    dp=n*(z*p1-p0)/(z*z-1);double d=p1/dp;z-=d;if(std::abs(d)<2e-16) {
          break;
        }
      }
   // Recompute derivative at the final root.
   double p0=1,p1=z;for(int k=2;k<=n;++k){double p=((2*k-1)*z*p1-(k-1)*p0)/k;p0=p1;p1=p;}
   dp=n*(z*p1-p0)/(z*z-1);
   x[i]=-z;x[n-1-i]=z;w[i]=w[n-1-i]=2/((1-z*z)*dp*dp);
  }
 }
};

    // Core-independent, polynomial-exact integration of tensor interpolation bases.
    class SegmentMomentIntegrator
    {
        static int checkedOrder(int p)
        {
            if(p<2||p>6)
            {
                throw std::invalid_argument("Segment moment order must be 2..6");
            }
            return p;
        }
        int p_;
        Eigen::Matrix3d coordinateTransform_;
        SegmentGaussRule stressRule_,concentrationRule_;
        std::vector<double> nodes_;
        void basis(const Eigen::Vector3d& x,const Eigen::Vector3d& center,const Eigen::Vector3d& width,double w[3][6]) const
        {
            const Eigen::Vector3d delta=coordinateTransform_*(x-center);
            for(int c=0;c<3;++c)
            {
                const double u=2*delta[c]/width[c];
                for(int j=0;j<p_;++j)
                {
                    w[c][j]=1;
                    for(int a=0;a<p_;++a)
                    {
                        if(a!=j)
                        {
                            w[c][j]*=(u-nodes_[a])/(nodes_[j]-nodes_[a]);
                        }
                    }
                }
            }
        }
    public:
        using StressMoments=Eigen::Matrix<double,Eigen::Dynamic,9>;
        using ConcentrationMoments=Eigen::Matrix<double,Eigen::Dynamic,2>;
        explicit SegmentMomentIntegrator(int p,const Eigen::Matrix3d& transform=Eigen::Matrix3d::Identity()):p_(checkedOrder(p)),coordinateTransform_(transform),stressRule_((3*p_-1)/2),concentrationRule_(3*p_/2),nodes_(p_)
        {
            if(p<2||p>6)
            {
                throw std::invalid_argument("Segment moment order must be 2..6");
            }
            for(int j=0;j<p;++j)
            {
                nodes_[j]=std::cos(std::numbers::pi*j/(p-1));
            }
        }
        std::size_t stressRuleSize() const { return stressRule_.x.size(); }
        StressMoments stress(const helmholtzSegment::Segment& s,const Eigen::Vector3d& center,const Eigen::Vector3d& width) const
        {
            const int P=p_*p_*p_;
            StressMoments result=StressMoments::Zero(P,9);
            const Eigen::Vector3d chord=s.p1-s.p0;
            for(std::size_t k=0;k<stressRule_.x.size();++k)
            {
                const Eigen::Vector3d x=s.p0+(.5+.5*stressRule_.x[k])*chord;
                double w[3][6];basis(x,center,width,w);
                for(int j=0;j<P;++j)
                {
                    const double weight=.5*stressRule_.w[k]*w[0][j%p_]*w[1][(j/p_)%p_]*w[2][j/(p_*p_)];
                    for(int a=0;a<3;++a)
                    {
                        for(int c=0;c<3;++c)
                        {
                            result(j,3*a+c)+=weight*s.b[a]*chord[c];
                        }
                    }
                }
            }
            return result;
        }
        ConcentrationMoments concentration(const Eigen::Vector3d& p0,const Eigen::Vector3d& p1,
            const Eigen::Vector3d& center,const Eigen::Vector3d& width,double begin,double end,const std::array<double,2>& sourceWeights) const
        {
            const int P=p_*p_*p_;
            ConcentrationMoments result=ConcentrationMoments::Zero(P,2);
            const Eigen::Vector3d chord=p1-p0;
            for(std::size_t k=0;k<concentrationRule_.x.size();++k)
            {
                const double u=(begin+end)*.5+(end-begin)*.5*concentrationRule_.x[k];
                const Eigen::Vector3d x=p0+u*chord;
                double w[3][6];basis(x,center,width,w);
                for(int j=0;j<P;++j)
                {
                    const double weight=chord.norm()*(end-begin)*.5*concentrationRule_.w[k]*w[0][j%p_]*w[1][(j/p_)%p_]*w[2][j/(p_*p_)];
                    result(j,0)+=weight*(1-u)*sourceWeights[0];
                    result(j,1)+=weight*u*sourceWeights[1];
                }
            }
            return result;
        }
    };
}
#endif
