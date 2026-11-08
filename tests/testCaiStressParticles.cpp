#include <DislocationDynamics/CaiStressParticleKernel.h>
#include <DislocationDynamics/KokkosDislocationFieldBuffers.h>
#include <DislocationDynamics/StraightDislocationSegment.h>
#include <Eigen/Geometry>
#include <DislocationDynamics/StraightDislocationSegmentCai.h>
#include <DislocationDynamics/StraightDislocationSegmentHelmholtz.h>
#include <iostream>
#include <random>
void require(bool ok,const char* message)
{
    if(!ok) { throw std::runtime_error(message); }
}
Eigen::Matrix3d stress(const model::CaiStressParticleKernel& kernel,const Eigen::Vector3d& r,const Eigen::Matrix3d& Q)
{
    double q[9],sum[6]={};
    for(int i=0;i<3;++i) { for(int j=0;j<3;++j) { q[3*i+j]=Q(i,j); } }
    kernel.addSymmetric(r.data(),q,sum);
    Eigen::Matrix3d result;
    int k=0;
    for(int i=0;i<3;++i)
    {
        for(int j=i;j<3;++j)
        {
            result(i,j)=result(j,i)=sum[k++];
        }
    }
    return result;
}
Eigen::Matrix3d integrate(const model::CaiStressParticleKernel& kernel,const Eigen::Vector3d& p0,
                        const Eigen::Vector3d& p1,const Eigen::Vector3d& b,const Eigen::Vector3d& x,int intervals)
{
    Eigen::Matrix3d result=Eigen::Matrix3d::Zero();
    for(int n=0;n<=intervals;++n)
    {
        const double w=(n==0 || n==intervals)?1.0:(n%2?4.0:2.0);
        const Eigen::Vector3d y=p0+(p1-p0)*(double(n)/intervals);
        result+=stress(kernel,x-y,b*(p1-p0).transpose()*(w/(3.0*intervals)));
    }
    return result;
}
void checkCaiConcentrationFlux()
{
    const Eigen::Matrix3d inverse=(Eigen::Vector3d(.7,1.2,2.3)).asDiagonal();
    const Eigen::Vector3d chord(1.3,-.7,.5);
    const double a=chord.dot(inverse*chord),core=.2;
    double concentrationError=0,fluxError=0;
    for(double longitudinal:{-4.,-.1,0.,.5,1.,1.1,4.})
    {
        for(double transverse:{0.,.01,1.})
        {
            const Eigen::Vector3d r=longitudinal*chord+transverse*chord.unitOrthogonal();
            const double b=-2*r.dot(inverse*chord),c=r.dot(inverse*r)+core;
            const Eigen::Vector2d concentration=model::caiSegmentConcentrationShapeIntegrals(a,b,c)/std::sqrt(a);
            for(const auto& strengths:{Eigen::Vector2d(1,1),Eigen::Vector2d(.3,-.7)})
            {
                const Eigen::Vector3d actual=model::caiSegmentFluxIntegral<double>(r,chord,inverse,a,core,strengths[0],strengths[1]);
                Eigen::Vector2d concentrationReference=Eigen::Vector2d::Zero();
                Eigen::Vector3d fluxReference=Eigen::Vector3d::Zero();
                constexpr int panels=8192;
                for(int i=0;i<=panels;++i)
                {
                    const double u=double(i)/panels,weight=(i==0||i==panels)?1:(i%2?4:2);
                    const Eigen::Vector3d separation=r-u*chord;
                    const double R=std::sqrt(separation.dot(inverse*separation)+core);
                    concentrationReference+=weight/(3*panels*R)*Eigen::Vector2d(1-u,u);
                    fluxReference+=weight/(3*panels*R*R*R)*((1-u)*strengths[0]+u*strengths[1])*separation;
                }
                concentrationError=std::max(concentrationError,(concentration-concentrationReference).norm()/concentrationReference.norm());
                fluxError=std::max(fluxError,(actual-fluxReference).norm()/std::max(1.0,fluxReference.norm()));
            }
        }
    }
    require(concentrationError<1e-11,"Cai concentration endpoint integral differs from line quadrature");
    require(fluxError<1e-10,"Cai flux endpoint integral differs from line quadrature");
    std::cout<<"Cai concentration relative / flux scaled errors: "<<concentrationError<<", "<<fluxError<<'\n';
}

int main(int argc,char** argv)
{
    try
    {
        require(argc==2,"Expected repository root");
        Kokkos::ScopeGuard guard(argc,argv);
        checkCaiConcentrationFlux();
        model::PolycrystallineMaterialBase material(std::string(argv[1])+"/Library/Materials/Al.txt",300);
        const Eigen::Matrix3d rotation=Eigen::AngleAxisd(0.71,Eigen::Vector3d(1,2,3).normalized()).toRotationMatrix();
        const Eigen::Vector3d start(-1,0,0),end(1,0,0),burgers(0,1,0),target(0,2,1);
        const auto snapshot=[&]()
        {
            model::StraightDislocationSegment<3> original(material,start,end,burgers);
            return model::StraightDislocationSegment<3>(original);
        }();
        require((snapshot.P0-start).norm()==0 && (snapshot.P1-end).norm()==0,"Snapshot copy lost geometry");
        Eigen::Vector3d liveEnd=end,tangent=(end-start).normalized();
        double length=(end-start).norm();
        model::StraightDislocationSegment<3> live(material,start,liveEnd,burgers,length,tangent);
        liveEnd+=Eigen::Vector3d(0,1,0);
        length=(liveEnd-start).norm();
        tangent=(liveEnd-start)/length;
        live.updateGeometry();
        model::StraightDislocationSegment<3> updated(material,start,liveEnd,burgers);
        require((live.stress(target)-updated.stress(target)).norm()<1e-14,"Live segment geometry was not updated");
        for(bool full:{false,true})
        {
            model::helmholtzSegment::Evaluator helmholtz(1,8,material.mu,material.nu,8,64,32,full);
            require((snapshot.stress(target,helmholtz)-helmholtz.stress(target,{start,end,burgers})).norm()<1e-14,
                    "Segment Helmholtz adapter changed stress");
            model::StraightDislocationSegment<3> zero(material,start,start,burgers);
            require(zero.stress(target).norm()==0 && zero.stress(target,helmholtz).norm()==0,"Zero-length segment stress");
        }
        double maxError=0;
        for(double a:{0.2,1.0,3.0})
        {
            model::DislocationFieldBase<3>::a=a; model::DislocationFieldBase<3>::a2=a*a;
            model::CaiStressParticleKernel kernel(material.mu,material.nu,a);
            const Eigen::Vector3d p0(-2,0.4,-0.7),p1(1.5,-0.1,0.8);
            const auto t=(p1-p0).normalized().eval();
            const Eigen::Vector3d edge=t.unitOrthogonal();
            for(const Eigen::Vector3d b:{t,edge,(t+edge).eval()})
            {
                model::StraightDislocationSegment<3> reference(material,p0,p1,b);
                for(const Eigen::Vector3d x:{Eigen::Vector3d(0.1,0.2,0.3),(0.5*(p0+p1)).eval(),(p0+0.01*edge).eval(),Eigen::Vector3d(9,-7,5)})
                {
                    const auto expected=reference.stress(x);
                    const auto numerical=integrate(kernel,p0,p1,b,x,8192);
                    const double error=(numerical-expected).norm()/std::max(1.0,expected.norm());
                    maxError=std::max(maxError,error);
                    require(error<1.e-9,"Particle integral differs from analytical segment stress");
                }
            }
            std::mt19937 generator(71);
            std::uniform_real_distribution<double> distribution(-3,3);
            for(int sample=0;sample<100;++sample)
            {
                Eigen::Vector3d r;
                Eigen::Matrix3d q;
                for(int i=0;i<3;++i)
                {
                    r[i]=distribution(generator);
                    for(int j=0;j<3;++j) { q(i,j)=distribution(generator); }
                }
                const auto full=stress(kernel,r,q);
                require((stress(kernel,-r,q)+full).norm()<1.e-13,"Separation sign");
                require((stress(kernel,r,-q)+full).norm()<1.e-13,"Burgers sign");
                require((stress(kernel,rotation*r,rotation*q*rotation.transpose())-rotation*full*rotation.transpose()).norm()<1.e-12,"Rotation covariance");
            }
            require(stress(kernel,Eigen::Vector3d::Zero(),Eigen::Matrix3d::Identity()).norm()==0,"Coincident pair not finite zero");
        }
        model::DislocationSourceArrays<3> sources;
        model::DislocationTargetArrays<3> targets;
        for(int p=0;p<17;++p)
        {
            sources.positions.append({double(p),0.2*p,-0.1*p});
            sources.weightedBurgersTensors.append({1,2,3,4,5,6,7,8,double(p)});
            targets.positions.append({double(p)+0.5,1,2});
        }
        model::KokkosDislocationFieldBuffers<3> device;
        device.uploadSources(sources); device.uploadTargets(targets);
        model::CaiStressParticleKernel kernel(1,0.3,1);
        model::DislocationStressArrays<3> symmetric;
        model::evaluateCaiParticles(device.executionSpace(),kernel,device.sourcePositions(),device.sourceTensors(),device.targetPositions(),device.stresses(),17,17);
        device.downloadStresses(symmetric);
        for(int p=0;p<17;++p)
        {
            Eigen::Matrix3d expected=Eigen::Matrix3d::Zero();
            for(int s=0;s<17;++s)
            {
                Eigen::Vector3d r; Eigen::Matrix3d q;
                for(int i=0;i<3;++i)
                {
                    r[i]=targets.positions.values[i][p]-sources.positions.values[i][s];
                    for(int j=0;j<3;++j) { q(i,j)=sources.weightedBurgersTensors.values[3*i+j][s]; }
                }
                expected+=stress(kernel,r,q);
            }
            for(int i=0;i<3;++i)
            {
                for(int j=0;j<3;++j)
                {
                    require(std::abs(symmetric.stresses.values[3*i+j][p]-expected(i,j))<1.e-12,"Six-component device differs from CPU");
                }
            }
        }
        model::evaluateCaiParticles(device.executionSpace(),kernel,device.sourcePositions(),device.sourceTensors(),device.targetPositions(),device.stresses(),0,17);
        device.downloadStresses(symmetric);
        for(const auto& v:symmetric.stresses.values) { for(double value:v) { require(value==0,"Empty source set not zero"); } }
        std::cout<<"Cai particle tests passed; max segment scaled error="<<maxError<<std::endl;
    }
    catch(const std::exception& error) { std::cerr<<error.what()<<std::endl; return 1; }
}
