#include <cfloat>
#include <DislocationDynamics/StraightDislocationSegmentFactory.h>
#include <type_traits>
#include <FMM/HelmholtzDiffusionKernel.h>
#include <FMM/HelmholtzStressKernelFamily.h>
#include <FMM/CaiStressKernelFamily.h>
#include <DislocationDynamics/ScreenedDislocationAssembly.h>
#include <DislocationDynamics/StraightDislocationSegmentHelmholtz.h>
#include <chrono>
#include <DislocationDynamics/StraightDislocationSegmentCai.h>
#include <iostream>
#include <random>
using namespace model;
void require(bool v,const char* msg){if(!v) {
    throw std::runtime_error(msg);
  }
}
void checkDiffusionEndpointTables()
{
    Eigen::Matrix3d D=(Eigen::Vector3d(.7,1.3,2.1)).asDiagonal();
    const double d=std::cbrt(D.determinant());const Eigen::Matrix3d M=d*D.inverse();
    model::HelmholtzDiffusionFamily family(D,1,8);
    model::helmholtzSegment::Rule rule(64);
    for(bool screened:{false,true})
    {
        const auto start=std::chrono::steady_clock::now();
        model::helmholtzSegment::ConcentrationEvaluator kernel(D,1,8,screened,8,64,2000);
        double error2=0,ref2=0,maxRelative=0,maxScaled=0;
        for(double rho:{0.,1e-10,1e-4,.01,.1,1.,10.,100.})
        {
            for(double za:{-100.,-1.,0.,1.,100.})
            {
                for(double length:{1e-8,.01,1.,20.})
                {
                    Eigen::Vector3d p0(za,0,0),p1(za+length,0,0),x(0,rho,0);
                    const auto actual=kernel.concentrationIntegrals(p0,p1,x);
                    const Eigen::Vector3d chord=p1-p0,r=x-p0;const double lm=std::sqrt(chord.dot(M*chord)),center=r.dot(M*chord)/(lm*lm);
                    const Eigen::Vector3d perp=r-center*chord;
                    const double scale=std::sqrt(perp.dot(M*perp)+1);
                    const double lo=std::asinh(-center*lm/scale),hi=std::asinh((1-center)*lm/scale);
                    Eigen::Vector2d reference=Eigen::Vector2d::Zero();
                    std::vector<double> cuts{lo};
                    if(lo<0&&hi>0)
                    {
                        cuts.push_back(0);
                    }
                    cuts.push_back(hi);
                    for(std::size_t c=1;c<cuts.size();++c)
                    {
                        int panels=std::max(1,int(std::ceil(2*(cuts[c]-cuts[c-1]))));
                        for(int j=0;j<panels;++j)
                        {
                            double a=cuts[c-1]+(cuts[c]-cuts[c-1])*j/panels,b=cuts[c-1]+(cuts[c]-cuts[c-1])*(j+1)/panels;
                            for(std::size_t k=0;k<rule.x.size();++k)
                            {
                                double q=(a+b)*.5+(b-a)*.5*rule.x[k],u=center+scale*std::sinh(q)/lm;
                                Eigen::Vector3d sep=x-(p0+u*chord);
                                double g=screened?family.shortKernel().value(sep.data()):family.fullKernel().value(sep.data());
                                double w=g*chord.norm()*scale*std::cosh(q)/lm*(b-a)*.5*rule.w[k];
                                reference+=w*Eigen::Vector2d(1-u,u);
                            }
                        }
                    }
                    double err=(actual-reference).norm();error2+=err*err;ref2+=reference.squaredNorm();
                    double relative=err/std::max(reference.norm(),1e-100);maxRelative=std::max(relative,maxRelative);
                    maxScaled=std::max(maxScaled,err/length);

                }
            }
        }
        require(std::sqrt(error2/ref2)<1e-7,"Diffusion endpoint table aggregate error");
        require(maxScaled<2e-8,"Diffusion endpoint table absolute error per unit length");
        std::cout<<"screened="<<screened<<" relative L2="<<std::sqrt(error2/ref2)<<" maxRelative="<<maxRelative<<" maxScaled="<<maxScaled<<" seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<'\n';
    }
}

void checkCoreSelection()
{
    const Eigen::Matrix3d D=Eigen::Vector3d(.7,1.3,2.1).asDiagonal();
    const Eigen::Vector3d p0(.1,.2,.3),p1(.8,.7,.6),x(.4,.8,.9);
    const double core=.1,d=std::cbrt(D.determinant());
    const Eigen::Matrix3d metric=d*D.inverse();
    helmholtzSegment::Rule rule(100);
    Eigen::Vector2d cai=Eigen::Vector2d::Zero();
    HelmholtzDiffusionFamily family(D,core,4*core,CoreType::Cai);
    for(std::size_t j=0;j<rule.x.size();++j)
    {
        const double u=(rule.x[j]+1)/2;
        const Eigen::Vector3d r=x-p0-u*(p1-p0);
        cai+=rule.w[j]*(p1-p0).norm()/2*family.fullKernel().value(r.data())*Eigen::Vector2d(1-u,u);
    }
    const Eigen::Vector2d endpoint=caiSegmentConcentrationIntegrals(D,core,p0,p1,x);
    require((endpoint-cai).norm()<1e-12,"Cai diffusion endpoints and FMM kernel disagree");
    const Eigen::Vector2d pieces=caiSegmentConcentrationIntegrals(D,core,p0,p1,x,0,.37)+caiSegmentConcentrationIntegrals(D,core,p0,p1,x,.37,1);
    require((pieces-endpoint).norm()<1e-12,"Cai split pieces lost original shape functions");
    helmholtzSegment::ConcentrationEvaluator full(D,core,4*core,false,8,64,4,true);
    const Eigen::Vector2d weights(.7,-.2);
    Eigen::Vector3d gradient;
    const double step=2e-5;
    for(int k=0;k<3;++k)
    {
        Eigen::Vector3d plus=x,minus=x;plus[k]+=step;minus[k]-=step;
        gradient[k]=(full.concentrationIntegrals(p0,p1,plus).dot(weights)-full.concentrationIntegrals(p0,p1,minus).dot(weights))/(2*step);
    }
    const Eigen::Vector3d flux=full.flux(p0,p1,x,weights[0],weights[1]);
    require((flux+D*gradient).norm()<1e-7*flux.norm(),"Helmholtz flux is inconsistent with concentration gradient");
    auto manager=std::make_shared<FMMManager>();
    auto tree=std::make_shared<BoxTree<3>>(Eigen::Vector3d::Zero(),Eigen::Vector3d::Ones(),0);
    ClimbKernelChannel channel;channel.diffusion=D;channel.segments.push_back({p0,p1,{0,1},{1,1}});
    ConcentrationTarget target;target.position={x[0],x[1],x[2]};
    Eigen::VectorXd velocity=weights;
    for(CoreType type:{CoreType::Cai,CoreType::Helmholtz})
    {
        ScreenedDislocationOptions options;options.coreType=type;options.screened=false;options.coreWidth=core;options.screenWidth=4*core;
        const double expected=(type==CoreType::Cai?endpoint:full.concentrationIntegrals(p0,p1,x)).dot(weights);
        for(unsigned depth:{0u,2u})
        {
            auto grid=std::make_shared<BoxTree<3>>(Eigen::Vector3d::Zero(),Eigen::Vector3d::Ones(),depth);
            options.interpolationOrder=6;
            DislocationConcentrationEvaluator evaluator(manager,grid,2,1,{channel},{target},options);
            require(std::abs(evaluator.evaluate(0,velocity)[0]-expected)<(depth?1e-4*std::abs(expected):1e-9),"Concentration evaluator ignored coreType");
        }
    }
    std::cout<<"Core selection, Cai split endpoints and Helmholtz flux gradient passed\n";
}
void checkSharedSplits()
{
    BoxTree<3> tree(Eigen::Vector3d::Zero(),Eigen::Vector3d::Ones(),2);
    SegmentSplitCache cache;
    std::vector<SegmentSplitSnapshot::Endpoints> segments={
        {Eigen::Vector3d(.02,.11,.19),Eigen::Vector3d(.97,.81,.71)},
        {Eigen::Vector3d(.25,.25,.25),Eigen::Vector3d(.75,.75,.75)},
        {Eigen::Vector3d(.4,.4,.4),Eigen::Vector3d(.4,.4,.4)}};
    const auto first=cache.get(tree,segments);
    require(cache.buildCount()==1,"Initial segment cache not built once");
    auto duplicate=segments;duplicate.insert(duplicate.end(),segments.begin(),segments.end());
    require(cache.get(tree,duplicate)==first&&cache.buildCount()==1,"Repeated channels rebuilt shared cuts");
    for(const auto& segment:segments)
    {
        double end=0;
        for(const auto& piece:first->pieces(segment.first,segment.second))
        {
            require(piece.begin==end&&piece.end>piece.begin,"Missing or overlapping split interval");
            require((piece.p0-(segment.first+piece.begin*(segment.second-segment.first))).norm()<1e-14,"Split start inconsistent");
            require((piece.p1-(segment.first+piece.end*(segment.second-segment.first))).norm()<1e-14,"Split end inconsistent");
            end=piece.end;
        }
        require(end==((segment.first==segment.second)?0:1),"Split intervals do not cover segment");
    }
    auto moved=segments;moved[0].second[0]-=.01;
    const auto second=cache.get(tree,moved);
    require(second!=first&&cache.buildCount()==2,"Motion failed to invalidate cuts");
    require(cache.get(tree,moved)==second&&cache.buildCount()==2,"Unchanged geometry rebuilt cuts");
    moved.pop_back();
    cache.get(tree,moved);require(cache.buildCount()==3,"Topology change failed to invalidate cuts");
    tree.setDepth(3);
    require(!second->matches(tree),"Old split snapshot accepted new tree");
    cache.get(tree,moved);require(cache.buildCount()==4,"Tree change failed to invalidate cuts");
    std::cout<<"Shared split cache reuse, motion, topology and tree invalidation passed\n";
}
void checkFluxTables()
{
    using namespace model::helmholtzSegment;
    const Eigen::Matrix3d D=Eigen::Vector3d(.7,1.3,2.1).asDiagonal();
    const double d=std::cbrt(D.determinant());
    const Eigen::Matrix3d M=d*D.inverse();
    ConcentrationEvaluator table(D,1,8,false,8,64,2000,true);
    Rule rule(64);
    double error2=0,reference2=0,maxScaled=0,maxRelative=0,splitError=0;
    for(double rho:{0.,1e-10,1e-4,.01,.1,1.,10.,100.})
    {
        for(double za:{-100.,-1.,0.,1.,100.})
        {
            for(double length:{1e-8,.01,1.,20.})
            {
                const Eigen::Vector3d p0(za,0,0),p1(za+length,0,0),x(0,rho,0),chord=p1-p0;
                const Eigen::Vector3d actual=table.flux(p0,p1,x,.7,-.2);
                const Eigen::Vector3d pieces=table.flux(p0,p1,x,.7,-.2,0,.37)+table.flux(p0,p1,x,.7,-.2,.37,1);
                splitError=std::max(splitError,(actual-pieces).norm()/length);
                const double lm=std::sqrt(chord.dot(M*chord)),center=(x-p0).dot(M*chord)/(lm*lm);
                const double scale=std::max(1e-12,rho*std::sqrt(M(1,1)));
                const double lo=std::asinh(-center*lm/scale),hi=std::asinh((1-center)*lm/scale);
                const int panels=std::max(1,int(std::ceil(2*(hi-lo))));
                Eigen::Vector3d reference=Eigen::Vector3d::Zero();
                for(int panel=0;panel<panels;++panel)
                {
                    const double a=lo+(hi-lo)*panel/panels,b=lo+(hi-lo)*(panel+1)/panels;
                    for(std::size_t k=0;k<rule.x.size();++k)
                    {
                        const double t=(a+b)*.5+(b-a)*.5*rule.x[k],u=center+scale*std::sinh(t)/lm;
                        const Eigen::Vector3d r=x-p0-u*chord;
                        const double R=std::sqrt(r.dot(M*r));
                        if(R>0)
                        {
                            const double f=R<.01?R*R*(.5+R*(-1./3+R*(1./8+R*(-1./30+R*(1./144+R*(-1./840+R/5760)))))):-std::expm1(-R)-R*std::exp(-R);
                            reference+=length/(4*pi)*(.7-.9*u)*f/(R*R*R)*r*scale*std::cosh(t)/lm*(b-a)*.5*rule.w[k];
                        }
                    }
                }
                if(length<1e-4)
                {
                    reference.setZero();
                    for(std::size_t k=0;k<rule.x.size();++k)
                    {
                        const double u=(rule.x[k]+1)*.5;
                        const Eigen::Vector3d r=x-p0-u*chord;
                        const double R=std::sqrt(r.dot(M*r));
                        if(R>0)
                        {
                            const double f=R<.01?R*R*(.5+R*(-1./3+R*(1./8+R*(-1./30+R*(1./144+R*(-1./840+R/5760)))))):-std::expm1(-R)-R*std::exp(-R);
                            reference+=chord.norm()/(4*pi)*(.7-.9*u)*f/(R*R*R)*r*.5*rule.w[k];
                        }
                    }
                }
                const double err=(actual-reference).norm();error2+=err*err;reference2+=reference.squaredNorm();
                const double scaled=err/length;
                if(scaled>maxScaled)
                {
                    maxScaled=scaled;

                }
                maxRelative=std::max(maxRelative,err/std::max(reference.norm(),1e-100));
            }
        }
    }
    std::cout<<"relative L2="<<std::sqrt(error2/reference2)<<" maxScaled="<<maxScaled<<" maxRelative="<<maxRelative<<" splitError="<<splitError<<'\n';
    require(std::sqrt(error2/reference2)<1e-7,"Helmholtz flux table aggregate error");
    require(maxScaled<2e-8,"Helmholtz flux table absolute error per segment length");
    require(splitError<1e-12,"Helmholtz flux split shapes are inconsistent");

}


void checkSegmentCoreDispatch()
{
    static_assert(std::is_abstract_v<StraightDislocationSegmentBase>);
    static_assert(std::is_abstract_v<SegmentConcentrationKernelBase>);
    SegmentCoreParameters parameters;
    parameters.coreWidth=.1;parameters.ewaldWidth=.8;parameters.maxDistance=4;
    parameters.mu=1.7;parameters.nu=.27;
    const helmholtzSegment::Segment segment{{.1,.2,.3},{.8,.7,.6},{.2,-.3,.7}};
    const Eigen::Vector3d x(.4,.8,.9),chord=segment.p1-segment.p0,t=chord.normalized(),bc=segment.b.cross(t);
    const Eigen::Matrix3d D=Eigen::Vector3d(.7,1.3,2.1).asDiagonal();
    const auto cai=makeStraightDislocationSegmentCore(CoreType::Cai,false,parameters);
    const auto full=makeStraightDislocationSegmentCore(CoreType::Helmholtz,false,parameters);
    const auto ewald=makeStraightDislocationSegmentCore(CoreType::Helmholtz,true,parameters);
    const auto reference=caiSegmentStress<double>(segment.p0,segment.p1,x,segment.b,t,bc,.01,1-parameters.nu,parameters.mu/(4*std::numbers::pi*(1-parameters.nu)));
    require((cai->stress(x,segment)-reference).norm()<1e-13,"Virtual Cai stress changed");
    helmholtzSegment::Evaluator fullReference(.1,.8,1.7,.27,8,64,4,true),shortReference(.1,.8,1.7,.27,8,64,4,false);
    require((full->stress(x,segment)-fullReference.stress(x,segment)).norm()==0,"Virtual Helmholtz stress changed");
    require((ewald->nearStress(x,segment)-shortReference.stress(x,segment)).norm()==0,"Virtual Ewald near stress changed");
    require((ewald->stress(x,segment)-ewald->shortStress(x,segment)-ewald->longStress(x,segment)).norm()<1e-14,"Ewald stress split changed");
    for(const auto& core:{cai,full,ewald})
    {
        auto concentration=core->concentration(D,true);
        Eigen::Vector2d expected;
        if(core==cai)
        {
            expected=caiSegmentConcentrationIntegrals(D,.1,segment.p0,segment.p1,x);
        }
        else
        {
            helmholtzSegment::ConcentrationEvaluator table(D,.1,.8,false,8,64,4,true);
            expected=table.concentrationIntegrals(segment.p0,segment.p1,x);
            require((concentration->flux(segment.p0,segment.p1,x,.7,-.2)-table.flux(segment.p0,segment.p1,x,.7,-.2)).norm()==0,"Virtual Helmholtz flux changed");
        }
        require((concentration->concentrationIntegrals(segment.p0,segment.p1,x)-expected).norm()<1e-14,"Virtual concentration changed");
        for(int order=2;order<=6;++order)
        {
            SegmentMomentIntegrator helper(order);
            const auto moments=core->stressMoments(helper,segment,Eigen::Vector3d::Constant(.5),Eigen::Vector3d::Ones());
            const auto weights=core->concentrationMoments(helper,segment.p0,segment.p1,Eigen::Vector3d::Constant(.5),Eigen::Vector3d::Ones(),.2,.9,{.7,-.2});
            for(int a=0;a<3;++a)
            {
                for(int b=0;b<3;++b)
                {
                    require(std::abs(moments.col(3*a+b).sum()-segment.b[a]*chord[b])<1e-13,"Stress moment partition of unity");
                }
            }
            require(std::abs(weights.col(0).sum()-chord.norm()*.7*.7*(1-.55))<1e-13,"Source endpoint zero moment");
            require(std::abs(weights.col(1).sum()-chord.norm()*.7*(-.2)*.55)<1e-13,"Source endpoint one moment");
        }
    }
    const auto cf=ewald->concentration(D,true),cs=ewald->concentration(D,true,SegmentKernelPart::shortRange),cl=ewald->concentration(D,true,SegmentKernelPart::longRange);
    require((cf->concentrationIntegrals(segment.p0,segment.p1,x)-cs->concentrationIntegrals(segment.p0,segment.p1,x)-cl->concentrationIntegrals(segment.p0,segment.p1,x)).norm()<1e-14,"Concentration split");
    require((cf->flux(segment.p0,segment.p1,x,.7,-.2)-cs->flux(segment.p0,segment.p1,x,.7,-.2)-cl->flux(segment.p0,segment.p1,x,.7,-.2)).norm()<1e-14,"Flux split");
    bool rejected=false;
    try
    {
        makeStraightDislocationSegmentCore(CoreType::Cai,true,parameters);
    }
    catch(const std::invalid_argument&)
    {
        rejected=true;
    }
    require(rejected,"Cai Ewald must be rejected");
    std::cout<<"Segment core polymorphism and shared moments passed\n";
}

int main()
{
    try
    {
        Kokkos::ScopeGuard guard;
        checkSegmentCoreDispatch();
        checkFluxTables();
        checkSharedSplits();
        checkCoreSelection();
        checkDiffusionEndpointTables();
        auto manager=std::make_shared<FMMManager>();
        auto tree=std::make_shared<BoxTree<3>>(Eigen::Vector3d::Zero(),Eigen::Vector3d::Ones(),0);
        ScreenedDislocationOptions opt;opt.boxShortRange=false;opt.interpolateScreenedNear=false;opt.coreWidth=.02;opt.screenWidth=.15;opt.shortCutoff=.3;
        std::mt19937 random(301);std::uniform_real_distribution<double> u(.01,.99);
        ClimbAssemblySnapshot snapshot;snapshot.nodeCount=7;snapshot.speciesCount=2;
        snapshot.rhs={Eigen::VectorXd::Ones(7),2*Eigen::VectorXd::Ones(7)};
        for(int t=0;t<19;++t)
        {
            ClimbQuadratureMap p;p.position={u(random),u(random),u(random)};p.nodes={std::size_t(t%7),std::size_t((t+1)%7)};
            p.weights={.03*(1-u(random)),-.01*u(random)};snapshot.targets.push_back(p);
        }
        for(int c=0;c<3;++c)
        {
            ClimbKernelChannel channel;channel.species=c%2;
            const auto R=Eigen::AngleAxisd(.4*c,Eigen::Vector3d(1,2,3).normalized()).toRotationMatrix();
            channel.diffusion=R*Eigen::Vector3d(.5+c,1+c,2+c).asDiagonal()*R.transpose();
            for(int s=0;s<27;++s)
            {
                ClimbQuadratureMap p;p.position={u(random),u(random),u(random)};p.nodes={std::size_t(s%7),std::size_t((s+2)%7)};
                p.weights={.04*(1-u(random)),.02*u(random)};channel.sources.push_back(p);
            }
            channel.sources[0].position=snapshot.targets[0].position; // Finite self interaction.
            snapshot.channels.push_back(channel);
        }
        std::vector<Eigen::MatrixXd> shortRef(2,Eigen::MatrixXd::Zero(7,7)),longRef=shortRef;
        for(const auto& ch:snapshot.channels)
        {
            HelmholtzDiffusionFamily family(ch.diffusion,opt.coreWidth,opt.screenWidth);
            for(const auto& target:snapshot.targets) {
        for(const auto& source:ch.sources)
            {
                double r[3];for(int c=0;c<3;++c) {
            r[c]=target.position[c]-source.position[c];
          }
          const double sh=family.shortKernel().distance(r)<=opt.shortCutoff?family.shortKernel().value(r):0;
                const double lo=family.longKernel().value(r);
                for(int i=0;i<2;++i) {
            for(int j=0;j<2;++j)
                {shortRef[ch.species](target.nodes[i],source.nodes[j])+=target.weights[i]*sh*source.weights[j];longRef[ch.species](target.nodes[i],source.nodes[j])+=target.weights[i]*lo*source.weights[j];}
          }
        }
      }
    }
        Eigen::VectorXd v(7);for(int i=0;i<7;++i) {
      v(i)=u(random)-.5;
    }
    {
            ScreenedClimbAssembly assembly(manager,tree,snapshot,opt);
            require(manager->interactionCount()==3,"Missing material/species registration");
            for(int s=0;s<2;++s)
            {
                require((Eigen::MatrixXd(assembly.shortMatrix(s))-shortRef[s]).norm()<1e-13,"Sparse local assembly mismatch");
                const auto ref=(shortRef[s]+longRef[s])*v;
                require((assembly.apply(s,v)-ref).norm()<1e-13,"Projected FMM operator mismatch");
                require((assembly.applyLong(s,2*v)-2*assembly.applyLong(s,v)).norm()<1e-13,"Long operator linearity");
                require((assembly.rhs(s)-snapshot.rhs[s]).norm()==0,"RHS changed");
            }
        }
        require(manager->interactionCount()==0&&manager->kernelCount()==0,"Assembly registration leak");
        // A nontrivial hierarchy uses the same assembled maps and sparse matrix.
        tree->setDepth(2);
        {
            ScreenedClimbAssembly assembly(manager,tree,snapshot,opt);
            for(int s=0;s<2;++s)
            {
                const Eigen::VectorXd ref=longRef[s]*v;
                const double error=(assembly.applyLong(s,v)-ref).norm()/ref.norm();
                std::cout<<"Climb long relative error species "<<s<<": "<<error<<'\n';
                require(error<.01,"Hierarchical Galerkin accuracy");
            }
        }
        // The sparse climb operator and concentration field must use identical
        // source/target box stencils, also for rotated anisotropic diffusion.
        {
            auto settings=opt;settings.boxShortRange=true;settings.interpolateScreenedNear=true;settings.screenWidth=4.;settings.interpolationOrder=6;
            settings.shortCutoff=1.e-30; // Unused by the production stencil.
            ScreenedClimbAssembly assembly(manager,tree,snapshot,settings);
            std::vector<ConcentrationTarget> targets;
            for(const auto& t:snapshot.targets){ConcentrationTarget p;p.position=t.position;targets.push_back(p);}
            DislocationConcentrationEvaluator field(manager,tree,7,2,snapshot.channels,targets,settings,false);
            for(int species=0;species<2;++species)
            {
                Eigen::MatrixXd reference=Eigen::MatrixXd::Zero(7,7);
                for(const auto& channel:snapshot.channels) {
          if(channel.species==species)
                {
                    const auto kernel=HelmholtzDiffusionFamily(channel.diffusion,settings.coreWidth,settings.screenWidth).shortKernel();
                    for(const auto& t:snapshot.targets) {
              for(const auto& s:channel.sources)
                    {
                        bool adjacent=true;double r[3];
                        for(int c=0;c<3;++c){adjacent=adjacent && std::abs(int(t.position[c]*4)-int(s.position[c]*4))<=1;r[c]=t.position[c]-s.position[c];}
                        if(adjacent) {
                  for(int i=0;i<2;++i) {
                    for(int j=0;j<2;++j) {
                      reference(t.nodes[i],s.nodes[j])+=t.weights[i]*kernel.value(r)*s.weights[j];
                    }
                  }
                }
              }
            }
          }
        }
        require((Eigen::MatrixXd(assembly.shortMatrix(species))-reference).norm()<1.e-13,"Climb box sparse assembly mismatch");
                const auto values=field.evaluate(species,v);Eigen::VectorXd projected=Eigen::VectorXd::Zero(7);
                for(std::size_t i=0;i<snapshot.targets.size();++i) {
          for(int j=0;j<2;++j) {
            projected(snapshot.targets[i].nodes[j])+=snapshot.targets[i].weights[j]*values(i);
          }
        }
        require((assembly.apply(species,v)-projected).norm()<1.e-12,"Climb and concentration use different box operators");
            }
        }
        // Unscreened mode applies the full kernel once; its local matrix is only a preconditioner.
        std::vector<Eigen::MatrixXd> fullRef(2,Eigen::MatrixXd::Zero(7,7));
        for(const auto& ch:snapshot.channels)
        {
            HelmholtzDiffusionFamily family(ch.diffusion,opt.coreWidth,opt.screenWidth);
            for(const auto& target:snapshot.targets) {
        for(const auto& source:ch.sources)
            {
                double r[3];for(int c=0;c<3;++c) {
            r[c]=target.position[c]-source.position[c];
          }
          const double g=family.fullKernel().value(r);
                for(int i=0;i<2;++i) {
            for(int j=0;j<2;++j) {
              fullRef[ch.species](target.nodes[i],source.nodes[j])+=target.weights[i]*g*source.weights[j];
            }
          }
        }
      }
    }
        for(int depth:{0,2})
        {
            tree->setDepth(depth);auto fullOpt=opt;fullOpt.screened=false;
            ScreenedClimbAssembly assembly(manager,tree,snapshot,fullOpt);
            for(int species=0;species<2;++species)
            {
                const Eigen::VectorXd ref=fullRef[species]*v;
                require((assembly.apply(species,v)-ref).norm()/ref.norm()<(depth?0.01:1e-12),"Full FMM climb mismatch or double counting");
                require((assembly.apply(species,v)-assembly.applyLong(species,v)).norm()<1e-13,"Unscreened local matrix added twice");
            }
        }
        // Field evaluation filters each grain channel, including shared-boundary targets.
        for(std::size_t i=0;i<snapshot.channels.size();++i) {
      snapshot.channels[i].grain=i+1;
    }
    std::vector<ConcentrationTarget> fieldTargets;
        for(std::size_t i=0;i<snapshot.targets.size();++i)
        {
            ConcentrationTarget target;target.position=snapshot.targets[i].position;
            if(i%3==0) {
        target.grains={1,2,3}; // Shared grain boundary.
      } else if(i%3==1) {
        target.grains={2}; // One grain.
      } else {
        target.grains={99}; // No eligible source.
      }
      fieldTargets.push_back(target);
        }
        tree->setDepth(0);
        for(bool screened:{true,false})
        {
            auto settings=opt;settings.screened=screened;
            DislocationConcentrationEvaluator field(manager,tree,7,2,snapshot.channels,fieldTargets,settings,true);
            for(std::size_t species=0;species<2;++species)
            {
                Eigen::VectorXd reference=Eigen::VectorXd::Zero(fieldTargets.size());
                for(const auto& ch:snapshot.channels) {
          if(ch.species==species)
                {
                    HelmholtzDiffusionFamily family(ch.diffusion,opt.coreWidth,opt.screenWidth);
                    for(std::size_t t=0;t<fieldTargets.size();++t)
                    {
                        const auto& target=fieldTargets[t];
                        if(std::find(target.grains.begin(),target.grains.end(),ch.grain)==target.grains.end()) {
                continue;
              }
              for(const auto& source:ch.sources)
                        {
                            double r[3];for(int c=0;c<3;++c) {
                  r[c]=target.position[c]-source.position[c];
                }
                double g=family.fullKernel().value(r);
                            if(screened) {
                  g=family.longKernel().value(r)+(family.shortKernel().distance(r)<=opt.shortCutoff?family.shortKernel().value(r):0);
                }
                reference(t)+=g*(source.weights[0]*v(source.nodes[0])+source.weights[1]*v(source.nodes[1]));
                        }
                    }
                }
        }
        require((field.evaluate(species,v)-reference).norm()<1e-12*reference.norm(),"Concentration grain filter/direct kernel mismatch");
                require((field.evaluate(species,2*v)-2*reference).norm()<2e-12*reference.norm(),"Concentration velocity update mismatch");
            }
        }
        // Stress assembly returns the full symmetric tensor and sums the same split.
        DislocationSourceArrays<3> src;DislocationTargetArrays<3> tgt;
        for(int s=0;s<31;++s){src.positions.append({u(random),u(random),u(random)});std::array<double,9> q;for(auto& x:q) {
        x=.01*(u(random)-.5);
      }
      src.weightedBurgersTensors.append(q);}
        for(int t=0;t<11;++t) {
      tgt.positions.append({u(random),u(random),u(random)});
    }
    tree->setDepth(2);
        for(bool screened:{true,false})
        {
            opt.screened=screened;
            FMMStressEvaluator stress(manager,tree,src,tgt,1,.3,opt);
            const auto actual=stress.evaluate();HelmholtzStressKernelFamily family(1,.3,opt.coreWidth,opt.screenWidth);
            for(int t=0;t<11;++t)
            {
                double sum[6]={};
                for(int s=0;s<31;++s){double r[3],q[9];for(int c=0;c<3;++c) {
            r[c]=tgt.positions.values[c][t]-src.positions.values[c][s];
          }
          for(int c=0;c<9;++c) {
            q[c]=src.weightedBurgersTensors.values[c][s];
          }
          if(screened){family.longKernel().add(r,q,sum);if(family.shortKernel().distance(r)<=opt.shortCutoff) {
              family.shortKernel().add(r,q,sum);
            }
          }else {
            family.fullKernel().add(r,q,sum);
          }
        }
                Eigen::Vector3d x;for(int axis=0;axis<3;++axis) {
          x[axis]=tgt.positions.values[axis][t];
        }
        const auto queried=stress.stressAt(x);
                for(int i=0;i<3;++i) {
          for(int j=0;j<3;++j) {
            require(std::abs(queried(i,j)-actual.stresses.values[3*i+j][t])<1e-11,"Cached point query differs from FMM batch");
          }
        }
        int c=0;for(int i=0;i<3;++i) {
          for(int j=i;j<3;++j){require(std::abs(actual.stresses.values[3*i+j][t]-sum[c++])<.005,"Stress split mismatch");require(actual.stresses.values[3*i+j][t]==actual.stresses.values[3*j+i][t],"Stress symmetry");}
        }
      }
        }
        // Production smooth-near path must agree between cached scalar and batch queries.
        opt.boxShortRange=true;opt.interpolateScreenedNear=true;opt.screened=true;opt.screenWidth=4.0;opt.interpolationOrder=6;
        for(int depth : {0,2})
        {
            tree->setDepth(depth);
            FMMStressEvaluator smooth(manager,tree,src,tgt,1,.3,opt);
            const auto batch=smooth.evaluate();
            require(smooth.stressAt(Eigen::Vector3d(-std::numeric_limits<double>::epsilon(),.5,.5)).allFinite(),"Boundary roundoff rejected");
            bool outsideRejected=false;try{smooth.stressAt(Eigen::Vector3d(-.01,.5,.5));}catch(const std::out_of_range&){outsideRejected=true;}
            require(outsideRejected,"Genuinely outside stress target accepted");
            opt.interpolateScreenedNear=false;
            FMMStressEvaluator direct(manager,tree,src,tgt,1,.3,opt);
            const auto reference=direct.evaluate();
            opt.interpolateScreenedNear=true;
            double e=0,n=0;
            for(int t=0;t<11;++t)
            {
                Eigen::Vector3d x;for(int c=0;c<3;++c) {
          x[c]=tgt.positions.values[c][t];
        }
        const auto point=smooth.stressAt(x);
                for(int c=0;c<9;++c)
                {
                    const double value=batch.stresses.values[c][t],ref=reference.stresses.values[c][t];
                    require(std::abs(point(c/3,c%3)-value)<1.e-10,"Smooth cached stress double counts neighbors");
                    e+=(value-ref)*(value-ref);n+=ref*ref;
                }
            }
            std::cout<<"Smooth stress relative error "<<std::sqrt(e/n)<<std::endl;
            require(std::sqrt(e/n)<.003,"Smooth stress all-box accuracy");
        }
        opt.screenWidth=.15;opt.boxShortRange=false;opt.interpolateScreenedNear=false;
        opt.coreType=CoreType::Cai;opt.screened=false;
        for(int depth : {0,2})
        {
            tree->setDepth(depth);
            FMMStressEvaluator stress(manager,tree,src,tgt,1,.3,opt);
            const auto actual=stress.evaluate();
            CaiStressParticleKernel direct(1,.3,opt.coreWidth);
            double error2=0,reference2=0;
            for(int t=0;t<11;++t)
            {
                double sum[6]={};
                for(int s=0;s<31;++s)
                {
                    double r[3],q[9];
                    for(int c=0;c<3;++c) {
            r[c]=tgt.positions.values[c][t]-src.positions.values[c][s];
          }
          for(int c=0;c<9;++c) {
            q[c]=src.weightedBurgersTensors.values[c][s];
          }
          direct.addSymmetric(r,q,sum);
                }
                Eigen::Vector3d x;for(int c=0;c<3;++c) {
          x[c]=tgt.positions.values[c][t];
        }
        const auto queried=stress.stressAt(x);
                int c=0;
                for(int i=0;i<3;++i) {
          for(int j=i;j<3;++j)
                {
                    const double value=actual.stresses.values[3*i+j][t];
                    require(std::abs(queried(i,j)-value)<1.e-11,"Cai cached stress query mismatch");
                    error2+=(value-sum[c])*(value-sum[c]);reference2+=sum[c]*sum[c];++c;
                }
        }
      }
            require(std::sqrt(error2/reference2)<(depth==0?1.e-12:.01),"Cai FMM disagrees with direct particles");
        }
        opt.screened=true;
        bool rejected=false;
        try { FMMStressEvaluator invalid(manager,tree,src,tgt,1,.3,opt); }
        catch(const std::invalid_argument&) { rejected=true; }
        require(rejected,"Cai must reject unsupported Ewald screening");
        require(manager->interactionCount()==0&&manager->kernelCount()==0,"Stress registration leak");
        // Segment sources retain linear endpoint velocities across leaf boundaries.
        {
            ClimbKernelChannel channel;
            channel.diffusion=(Eigen::Vector3d(.7,1.3,2.1)).asDiagonal();
            channel.segments.push_back({Eigen::Vector3d(.03,.13,.21),Eigen::Vector3d(.97,.83,.71),{0,1},{.8,-.3}});
            channel.segments.push_back({Eigen::Vector3d(.12,.91,.87),Eigen::Vector3d(.81,.04,.31),{2,3},{-.2,1.1}});
            std::vector<ConcentrationTarget> points;
            for(const Eigen::Vector3d x:{Eigen::Vector3d(.5,.48,.46),Eigen::Vector3d(.06,.05,.96),Eigen::Vector3d(.9,.9,.05)})
            {
                ConcentrationTarget q;
                for(int c=0;c<3;++c)
                {
                    q.position[c]=x[c];
                }
                q.grains={0};points.push_back(q);
            }
            Eigen::VectorXd velocity(4);velocity<<.7,-.4,.2,.9;
            ScreenedDislocationOptions options;
            options.coreWidth=.02;options.screenWidth=.12;options.interpolationOrder=6;
            for(bool screened:{false,true})
            {
                options.screened=screened;
                for(int depth:{0,2})
                {
                    // Keep the smoothing width resolved by the leaf interpolation.
                    options.screenWidth=.5/(1<<depth);
                    auto segmentTree=std::make_shared<BoxTree<3>>(Eigen::Vector3d::Zero(),Eigen::Vector3d::Ones(),depth);
                    DislocationConcentrationEvaluator evaluator(manager,segmentTree,4,1,{channel},points,options);
                    const Eigen::VectorXd actual=evaluator.evaluate(0,velocity);
                    require((evaluator.evaluate(0,2*velocity)-2*actual).norm()<1e-12,"Segment source moments were not refreshed for new velocities");
                    HelmholtzDiffusionFamily family(channel.diffusion,options.coreWidth,options.screenWidth);
                    Eigen::VectorXd reference=Eigen::VectorXd::Zero(points.size());
                    for(std::size_t i=0;i<points.size();++i)
                    {
                        const Eigen::Vector3d x=Eigen::Map<const Eigen::Vector3d>(points[i].position.data());
                        for(const auto& segment:channel.segments)
                        {
                            const Eigen::Vector3d chord=segment.p1-segment.p0;
                            constexpr int count=32768;
                            for(int j=0;j<=count;++j)
                            {
                                const double u=double(j)/count;
                                const Eigen::Vector3d y=segment.p0+u*chord,separation=x-y;
                                double g=family.fullKernel().value(separation.data());
                                if(screened)
                                {
                                    bool near=true;
                                    for(int c=0;c<3;++c)
                                    {
                                        near=near&&std::abs(int(x[c]*(1<<depth))-int(y[c]*(1<<depth)))<=1;
                                    }
                                    g=family.longKernel().value(separation.data())+(near?family.shortKernel().value(separation.data()):0);
                                }
                                const double weight=(j==0||j==count)?1:(j%2?4:2);
                                reference[i]+=g*chord.norm()*weight/(3*count)*((1-u)*segment.weights[0]*velocity[segment.nodes[0]]+u*segment.weights[1]*velocity[segment.nodes[1]]);
                            }
                        }
                    }
                    const double error=(actual-reference).norm()/reference.norm();
                    std::cout<<"Segment concentration screened="<<screened<<", depth="<<depth<<", relative error="<<error<<'\n';
                    require(error<(!screened?(depth==0?1e-9:1e-5):(depth==0?.004:.0002)),"Segment concentration disagrees with independent line quadrature");
                    ClimbAssemblySnapshot segmentSnapshot;segmentSnapshot.nodeCount=4;segmentSnapshot.speciesCount=1;
                    segmentSnapshot.channels={channel};segmentSnapshot.rhs={Eigen::VectorXd::Zero(4)};
                    Eigen::VectorXd expected=Eigen::VectorXd::Zero(4);
                    for(std::size_t i=0;i<points.size();++i)
                    {
                        ClimbQuadratureMap q;q.position=points[i].position;q.nodes={i%4,(i+1)%4};q.weights={.3,-.2};
                        segmentSnapshot.targets.push_back(q);
                        for(int k=0;k<2;++k)
                        {
                            expected[q.nodes[k]]+=q.weights[k]*actual[i];
                        }
                    }
                    ScreenedClimbAssembly assembly(manager,segmentTree,segmentSnapshot,options);
                    require((assembly.apply(0,velocity)-expected).norm()<1e-12,"Climb assembly and CD segment evaluation differ");
                    auto filtered=points;filtered[0].grains={1};
                    DislocationConcentrationEvaluator filteredEvaluator(manager,segmentTree,4,1,{channel},filtered,options,true);
                    const auto filteredValues=filteredEvaluator.evaluate(0,velocity);
                    require(filteredValues[0]==0 && (filteredValues.tail(2)-actual.tail(2)).norm()<1e-12,"Segment concentration grain filter");
                }
            }
        }
        std::cout<<"Screened stress and climb assembly tests passed\n";
    }
    catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}
}
