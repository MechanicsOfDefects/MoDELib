#include <FiniteElement.h>
#include <WeakProblem.h>
#include <Eigen/SparseLU>
#include <Eigen/IterativeLinearSolvers>
#include <filesystem>
#include <numeric>
#include <iostream>

using namespace model;
using Sparse=Eigen::SparseMatrix<double>;
using Direct=Eigen::SparseLU<Sparse>;
using Iterative=Eigen::BiCGSTAB<Sparse>;
using Mesh=SimplicialMesh<3>;

void require(bool value,const char* message)
{
    if(!value)
    {
        throw std::runtime_error(message);
    }
}
template<class F>
void rejects(F&& operation,const char* message)
{
    bool rejected=false;
    try
    {
        operation();
    }
    catch(const std::exception&)
    {
        rejected=true;
    }
    require(rejected,message);
}

// Independent oracle: construct equivalence classes from the pair edges,
// offsets from coordinates, and perform the explicit dense T^T*A*T product.
template<class Field,class BWF,class LWF,class Problem>
void compare(const std::shared_ptr<Field>& u,const std::shared_ptr<BWF>& bwf,
             const std::shared_ptr<LWF>& lwf,const std::shared_ptr<Problem>& problem)
{
    const auto& fe=u->fe();
    constexpr int components=Field::dofPerNode;
    std::vector<size_t> roots(fe.nodeSize());
    std::iota(roots.begin(),roots.end(),0);
    for(const auto& pair : fe.periodicNodePairs())
    {
        const auto a=roots[pair.source];
        const auto b=roots[pair.sink];
        for(auto& root : roots)
        {
            if(root==a || root==b)
            {
                root=std::min(a,b);
            }
        }
    }
    const Eigen::Index n=u->gSize();
    Eigen::VectorXd g(n);
    std::map<size_t,double> fixed;
    for(size_t node=0;node<fe.nodeSize();++node)
    {
        g.segment(node*components,components)=u->periodicGradient()*(fe.node(node).P0-fe.node(roots[node]).P0);
    }
    for(const auto& bc : u->dirichletConditions())
    {
        fixed[roots[bc.first/components]*components+bc.first%components]=bc.second-g(bc.first);
    }
    std::map<size_t,size_t> columns;
    for(size_t i=0;i<size_t(n);++i)
    {
        const auto root=roots[i/components]*components+i%components;
        if(!fixed.count(root))
        {
            columns.emplace(root,columns.size());
        }
    }
    Eigen::MatrixXd T=Eigen::MatrixXd::Zero(n,columns.size());
    for(size_t i=0;i<size_t(n);++i)
    {
        const auto root=roots[i/components]*components+i%components;
        if(fixed.count(root))
        {
            g(i)+=fixed.at(root);
        }
        else
        {
            T(i,columns.at(root))=1.0;
        }
    }
    Eigen::MatrixXd A=Eigen::MatrixXd::Zero(n,n);
    for(const auto& t : bwf->globalTriplets())
    {
        A(t.row(),t.col())+=t.value();
    }
    const Eigen::MatrixXd reduced=T.transpose()*A*T;
    const Eigen::VectorXd b=T.transpose()*(lwf->globalVector()-A*g);
    const Eigen::VectorXd expected=T*reduced.fullPivLu().solve(b)+g;
    require((Eigen::MatrixXd(problem->reducedMatrix())-reduced).norm()<1.e-10*std::max(1.0,reduced.norm()),"Constrained matrix differs from T^T*A*T");
    require((problem->reducedForceVector()-b).norm()<1.e-10*std::max(1.0,b.norm()),"Constrained force differs from T^T*(b-A*g)");
    require((u->dofVector()-expected).norm()<2.e-8*std::max(1.0,expected.norm()),"Reconstructed solution differs from reference");
}

template<int order,int components>
void exercise(const std::filesystem::path& root,bool periodic,bool iterative)
{
    auto mesh=std::make_shared<Mesh>((root/"Library/Meshes/bicrystal_48.msh").string(),Eigen::Matrix3d::Identity(),Eigen::Vector3d::Zero(),periodic ? std::set<int>{-1} : std::set<int>{});
    using FE=FiniteElement<LagrangeElement<3,order>>;
    auto fe=std::make_shared<FE>(mesh);
    auto u=fe->template trial<components>();
    auto other=fe->template trial<components>();
    const auto dV=fe->template domain<EntireDomain,4,GaussLegendre>();
    const auto mass=(test(u),u)*dV;
    const auto bwf=mass+(test(grad(u)),grad(u))*dV;
    const Eigen::Matrix<double,components,1> ones=Eigen::Matrix<double,components,1>::Ones();
    auto lwf=(test(u),make_constant(ones))*dV;
    WeakProblemOptions options;
    options.useDirectSolver=!iterative;
    options.tolerance=1.e-9;
    options.maxIterations=4000;
    auto problem=weakProblem<Direct,Iterative>(bwf,lwf,options);
    require(!problem->solveWithMasterSlaveConstraints().reusedFactorization,"First solve claimed reuse");
    compare(u,bwf,lwf,problem);
    require((u->dofVector().array()-1.0).matrix().norm()<1.e-7,"Constant solution changed");
    require(problem->factorizationCount()==1,"Initial setup count incorrect");

    lwf=(test(u),make_constant(Eigen::Matrix<double,components,1>(3.0*ones)))*dV;
    problem->setRightHandSide(lwf);
    require(problem->solveWithMasterSlaveConstraints().reusedFactorization,"New RHS lost factorization");
    compare(u,bwf,lwf,problem);
    require(problem->factorizationCount()==1,"New RHS refactorized matrix");
    const auto previous=u->dofVector().eval();
    const auto foreign=(test(other),make_constant(ones))*dV;
    rejects([&]() { problem->setRightHandSide(foreign); },"Foreign RHS accepted");
    rejects([&]() { auto invalid=weakProblem<Direct,Iterative>(bwf,foreign); },"Mixed fields accepted at construction");
    require((u->dofVector()-previous).norm()==0.0,"Failed RHS replacement modified the field");

    u->dirichletConditions()[0]=2.0;
    problem->solveWithMasterSlaveConstraints();
    compare(u,bwf,lwf,problem);
    require(problem->factorizationCount()==2,"New constraint mapping did not rebuild");
    u->dirichletConditions()[0]=-0.4;
    require(problem->solveWithMasterSlaveConstraints().reusedFactorization,"Prescribed-value change refactorized matrix");
    compare(u,bwf,lwf,problem);

    Eigen::Matrix<double,components,3> gradient;
    for(int c=0;c<components;++c)
    {
        gradient.row(c)<<0.2*(c+1),-0.3*(c+1),0.4*(c+1);
    }
    u->setPeriodicGradient(gradient);
    require(problem->solveWithMasterSlaveConstraints().reusedFactorization,"Affine offset change refactorized matrix");
    compare(u,bwf,lwf,problem);

    if(periodic)
    {
        const auto pair=fe->periodicNodePairs().front();
        const size_t source=pair.source*components;
        const size_t sink=pair.sink*components;
        u->dirichletConditions().clear();
        u->dirichletConditions()[source]=0.7;
        u->dirichletConditions()[sink]=0.7+(gradient*pair.shift)(0);
        problem->solveWithMasterSlaveConstraints();
        compare(u,bwf,lwf,problem);
        const auto saved=u->dofVector().eval();
        u->dirichletConditions()[sink]+=1.0;
        rejects([&]() { problem->solveWithMasterSlaveConstraints(); },"Conflicting periodic Dirichlet values accepted");
        require((u->dofVector()-saved).norm()==0.0,"Failed constraint validation modified field");
        u->dirichletConditions().clear();
        // Select one active lattice direction, then remove periodicity entirely.
        int seed=-1;
        for(const auto& region : mesh->regions())
        {
            for(const auto& face : region.second->faces())
            {
                const auto& p=face.second->periodicFacePair;
                if(p.second && p.first.cross(mesh->periodicBasis().col(0)).norm()<1.e-8)
                {
                    seed=face.first;
                }
            }
        }
        require(seed>=0,"No partial-periodic seed face");
        mesh->identifyParallelFaces({seed});
        fe->updatePeriodicNodePairs();
        require(!problem->solveWithMasterSlaveConstraints().reusedFactorization,"Partial-periodic mapping did not rebuild");
        compare(u,bwf,lwf,problem);
        mesh->identifyParallelFaces({});
        fe->updatePeriodicNodePairs();
        require(!problem->solveWithMasterSlaveConstraints().reusedFactorization,"Removed periodicity did not rebuild");
        compare(u,bwf,lwf,problem);
    }
    problem->reassembleMatrix();
    require(!problem->solveWithMasterSlaveConstraints().reusedFactorization,"Explicit reassembly reused stale setup");
    compare(u,bwf,lwf,problem);
}

void failuresAndReuse(const std::filesystem::path& root)
{
    auto mesh=std::make_shared<Mesh>((root/"Library/Meshes/bicrystal_48.msh").string(),Eigen::Matrix3d::Identity(),Eigen::Vector3d::Zero(),std::set<int>{});
    using FE=FiniteElement<LagrangeElement<3,1>>;
    auto fe=std::make_shared<FE>(mesh);
    auto u=fe->trial<1>();
    auto coefficient=fe->trial<1>();
    *coefficient=1.0;
    auto dV=fe->domain<EntireDomain,4,GaussLegendre>();
    auto bwf=(test(u),eval(coefficient)*u)*dV;
    auto lwf=(test(u),make_constant(1.0))*dV;
    auto problem=weakProblem<Direct,Iterative>(bwf,lwf);
    problem->solveWithMasterSlaveConstraints();
    *coefficient=2.0;
    require(problem->solveWithMasterSlaveConstraints().reusedFactorization,"Coefficient snapshot unexpectedly rebuilt");
    require((u->dofVector().array()-1.0).matrix().norm()<1.e-9,"Cached coefficient snapshot changed implicitly");
    problem->reassembleMatrix();
    problem->solveWithMasterSlaveConstraints();
    require((u->dofVector().array()-0.5).matrix().norm()<1.e-9,"Changed coefficient not reassembled");
    require(problem->factorizationCount()==2,"Coefficient change setup count incorrect");
    compare(u,bwf,lwf,problem);
    const auto saved=u->dofVector().eval();
    *coefficient=0.0;
    problem->reassembleMatrix();
    rejects([&]() { problem->solveWithMasterSlaveConstraints(); },"Singular system accepted");
    require((u->dofVector()-saved).norm()==0.0,"Failed solve modified field");
    *coefficient=4.0;
    problem->reassembleMatrix();
    problem->solveWithMasterSlaveConstraints();
    require((u->dofVector().array()-0.25).matrix().norm()<1.e-9,"Recovery after failed factorization failed");
    auto rhsField=fe->trial<1>();
    *rhsField=2.0;
    auto liveRhs=(test(u),eval(rhsField))*dV;
    problem->setRightHandSide(liveRhs);
    require(problem->solveWithMasterSlaveConstraints().reusedFactorization,"Live RHS lost setup");
    require((u->dofVector().array()-0.5).matrix().norm()<1.e-9,"Live RHS not evaluated");
    *rhsField=3.0;
    require(problem->solveWithMasterSlaveConstraints().reusedFactorization,"Changed live RHS lost setup");
    require((u->dofVector().array()-0.75).matrix().norm()<1.e-9,"Changed live RHS was stale");
    const auto beforeInvalid=u->dofVector().eval();
    *rhsField=std::numeric_limits<double>::quiet_NaN();
    rejects([&]() { problem->solveWithMasterSlaveConstraints(); },"Nonfinite RHS accepted");
    require((u->dofVector()-beforeInvalid).norm()==0.0,"Invalid RHS changed trial");
    problem->setRightHandSide(lwf);
    for(size_t i=0;i<rhsField->gSize();++i)
    {
        rhsField->dofVector()(i)=std::sin(0.37*double(i));
    }
    WeakProblemOptions limited;
    limited.useDirectSolver=false;
    limited.maxIterations=1;
    limited.tolerance=1.e-14;
    auto nonconverging=weakProblem<Direct,Iterative>(bwf,liveRhs,limited);
    rejects([&]() { nonconverging->solveWithMasterSlaveConstraints(); },"Insufficient iterative solve accepted");
    require((u->dofVector()-beforeInvalid).norm()==0.0,"Nonconverged iteration changed trial");
    auto spd=weakProblem<Eigen::SimplicialLDLT<Sparse>,Eigen::ConjugateGradient<Sparse>>(bwf,lwf);
    spd->solveWithMasterSlaveConstraints();
    require((u->dofVector().array()-0.25).matrix().norm()<1.e-9,"Alternative solver types failed");
    for(size_t i=0;i<u->gSize();++i)
    {
        u->dirichletConditions()[i]=double(i)/10;
    }
    const size_t count=problem->factorizationCount();
    problem->solveWithMasterSlaveConstraints();
    require(problem->reducedMatrix().rows()==0,"Fully prescribed problem has unknowns");
    require(problem->factorizationCount()==count,"Empty system called solver");
    for(size_t i=0;i<u->gSize();++i)
    {
        require(u->dofVector()(i)==double(i)/10,"Prescribed-only reconstruction failed");
    }
    const auto allFixed=u->dofVector().eval();
    u->dirichletConditions()[u->gSize()]=0.0;
    rejects([&]() { problem->solveWithMasterSlaveConstraints(); },"Out-of-range constraint accepted");
    require((u->dofVector()-allFixed).norm()==0.0,"Invalid constraint modified field");
}
void ownership(const std::filesystem::path& root)
{
    using FE=FiniteElement<LagrangeElement<3,1>>;
    using Field=TrialFunction<'u',1,FE>;
    std::weak_ptr<Field> weakField;
    std::weak_ptr<FE> weakSpace;
    auto problem=[&]()
    {
        auto mesh=std::make_shared<Mesh>((root/"Library/Meshes/bicrystal_48.msh").string(),Eigen::Matrix3d::Identity(),Eigen::Vector3d::Zero(),std::set<int>{-1});
        auto fe=std::make_shared<FE>(mesh);
        auto u=fe->trial<1>();
        weakField=u;
        weakSpace=fe;
        auto dV=fe->domain<EntireDomain,4,GaussLegendre>();
        return weakProblem((test(u),u)*dV,(test(u),make_constant(2.0))*dV);
    }();
    problem->solveWithMasterSlaveConstraints();
    require(!weakField.expired() && !weakSpace.expired(),"Problem did not retain trial and FEM space");
    require((weakField.lock()->dofVector().array()-2.0).matrix().norm()<1.e-9,"Detached problem solution invalid");
    problem.reset();
    require(weakField.expired() && weakSpace.expired(),"WeakProblem leaked an ownership cycle");
}

int main(int argc,char** argv)
{
    try
    {
        require(argc==2,"Expected repository root");
        exercise<1,1>(argv[1],false,false);
        exercise<1,1>(argv[1],true,false);
        exercise<2,1>(argv[1],true,false);
        exercise<1,3>(argv[1],true,false);
        exercise<1,1>(argv[1],true,true);
        failuresAndReuse(argv[1]);
        ownership(argv[1]);
        std::cout<<"PASS WeakProblem constrained triplets, affine periodicity, repeated RHS, reuse and atomic failure"<<std::endl;
    }
    catch(const std::exception& e)
    {
        std::cerr<<e.what()<<std::endl;
        return 1;
    }
}
