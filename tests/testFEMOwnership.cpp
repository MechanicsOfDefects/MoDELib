#include <FiniteElement.h>
#include <FixedDirichletSolver.h>
#include <Eigen/SparseCholesky>
#include <filesystem>
#include <iostream>
#include <optional>

using namespace model;
using Mesh=SimplicialMesh<3>;
using Space=FiniteElement<LagrangeElement<3,1>>;
using Field=TrialFunction<'c',1,Space>;

void require(bool ok,const std::string& message)
{
    if(!ok)
    {
        throw std::runtime_error(message);
    }
}
template<class F>
void rejects(F&& f,const std::string& message)
{
    bool threw=false;
    try
    {
        f();
    }
    catch(const std::exception&)
    {
        threw=true;
    }
    require(threw,message);
}
template<class Form>
Eigen::MatrixXd matrix(const Form& form)
{
    Eigen::SparseMatrix<double> result(form->field().gSize(),form->field().gSize());
    const auto entries=form->globalTriplets();
    result.setFromTriplets(entries.begin(),entries.end());
    return result;
}
std::shared_ptr<const Mesh> makeMesh(const std::filesystem::path& root,double scale=1.0)
{
    return std::make_shared<const Mesh>((root/"Library/Meshes/bicrystal_48.msh").string(),scale*Eigen::Matrix3d::Identity(),Eigen::Vector3d::Zero(),std::set<int>{-1});
}

void independenceAndAssembly(const std::filesystem::path& root)
{
    auto mesh=makeMesh(root);
    auto fe=std::make_shared<Space>(mesh);
    auto other=std::make_shared<Space>(makeMesh(root,2.0));
    auto u=std::make_shared<Field>(*fe);
    *u=2.0;
    auto v=std::make_shared<Field>(*fe);
    auto w=std::make_shared<Field>(*other);
    *v=7.0;
    *w=11.0;
    const auto& ele=fe->elements().begin()->second;
    const Eigen::Vector4d bary=Eigen::Vector4d::Constant(.25);
    require(std::abs((*eval(u))(ele,bary)(0)-2.0)<1.e-12,"Second same-type field overwrote first field");
    require(std::abs((*eval(v))(ele,bary)(0)-7.0)<1.e-12,"Independent field values lost");
    require(&u->fe()==fe.get() && &w->fe()==other.get(),"FEM pointer shared across instances");
    u->dirichletConditions().emplace(0,2.0);
    u->dirichletNodeMap()[0][0]=true;
    require(v->dirichletConditions().empty() && w->dirichletNodeMap().empty(),"Constraints leaked across fields");
    auto alias=u;
    *alias=3.0;
    require(std::abs((*eval(u))(ele,bary)(0)-3.0)<1.e-12,"Copy handle lost field identity");
    auto beforeAssignment=eval(u);
    *u=*v;
    *v=9.0;
    require(std::abs((*beforeAssignment)(ele,bary)(0)-7.0)<1.e-12,"Assignment rebound an existing expression");
    auto savedGradient=[&]()
    {
        auto local=2.0*grad(u);
        return eval(local);
    }();
    require((*savedGradient)(ele,bary).norm()<1.e-11,"Saved intermediate expression lost ownership");
    std::weak_ptr<Field> coefficientField;
    auto coefficient=[&]()
    {
        auto c=std::make_shared<Field>(*fe);
        *c=4.0;
        coefficientField=c;
        return eval(c)*u;
    }();
    require(std::abs((*eval(coefficient))(ele,bary)(0)-28.0)<1.e-11,"Evaluated coefficient did not retain its field");

    require(!coefficientField.expired(),"Coefficient field was not retained");
    auto sum=u+u;
    static_assert(std::is_same_v<decltype(sum),std::shared_ptr<TrialSum<Field,Field>>>);
    static_assert(!std::is_copy_constructible_v<Field>);
    require(sum->op1.get()==u.get() && sum->op2.get()==u.get(),"Sum copied its operands");
    require(std::abs((*eval(sum))(ele,bary)(0)-14.0)<1.e-12,"Trial sum evaluation changed");
    rejects([&]() { auto invalid=u+v; },"Sum accepted independent fields");
    rejects([&]() { auto invalid=grad(std::shared_ptr<Field>{}); },"Null operand accepted");
    require((*eval(grad(2.0*u)))(ele,bary).norm()<1.e-11,"Gradient of scalar product changed");

    auto dV=fe->domain<EntireDomain,4,GaussLegendre>();
    auto mass=(test(u),u)*dV;
    auto stiffness=(test(grad(u)),grad(u))*dV;
    auto load=(test(u),make_constant(1.0))*dV;
    auto M=matrix(mass);
    auto K=matrix(stiffness);
    Eigen::MatrixXd expectedM=Eigen::MatrixXd::Zero(u->gSize(),u->gSize());
    Eigen::MatrixXd expectedK=expectedM;
    Eigen::VectorXd expectedLoad=Eigen::VectorXd::Zero(u->gSize());
    for(const auto& entry : fe->elements())
    {
        const auto& e=entry.second;
        Eigen::Matrix4d coordinates;
        for(int i=0;i<4;++i)
        {
            coordinates(i,0)=1.0;
            coordinates.row(i).tail<3>()=e.node(i).P0.transpose();
        }
        const Eigen::Matrix<double,3,4> gradients=coordinates.inverse().bottomRows<3>();
        for(int i=0;i<4;++i)
        {
            const size_t I=e.node(i).gID;
            expectedLoad(I)+=e.simplex.vol0/4.0;
            for(int j=0;j<4;++j)
            {
                const size_t J=e.node(j).gID;
                expectedM(I,J)+=e.simplex.vol0*(i==j ? 2.0 : 1.0)/20.0;
                expectedK(I,J)+=e.simplex.vol0*gradients.col(i).dot(gradients.col(j));
            }
        }
    }
    require((M-expectedM).norm()<1.e-11*expectedM.norm(),"Mass assembly changed");
    require((K-expectedK).norm()<1.e-11*expectedK.norm(),"Diffusion assembly changed");
    require((load->globalVector()-expectedLoad).norm()<1.e-11*expectedLoad.norm(),"Load assembly changed");
    require((matrix(mass+mass)-2*M).norm()<1.e-11*M.norm(),"Bilinear sum changed");
    require(((load+load-load)->globalVector()-expectedLoad).norm()<1.e-11*expectedLoad.norm(),"Linear sum/difference changed");
    rejects([&]() { matrix((test(u),v)*dV); },"Mixed field bilinear form accepted");
    rejects([&]() { matrix(mass+(test(v),v)*dV); },"Mixed field bilinear sum accepted");
    rejects([&]() { (load+(test(v),make_constant(1.0))*dV)->globalVector(); },"Mixed field load sum accepted");

    rejects([&]() { matrix((test(w),w)*dV); },"Foreign FEM integration domain accepted");

    auto retainedMass=[&]()
    {
        auto localDomain=fe->domain<EntireDomain,4,GaussLegendre>();
        auto localForm=(test(u),u);
        return localForm*localDomain;
    }();
    require((matrix(retainedMass)-M).norm()<1.e-12,"Temporary integration domain was borrowed");
    auto boundary=[&]()
    {
        auto localDomain=fe->boundary<ExternalBoundary,4,GaussLegendre>();
        return (test(u),make_constant(1.0))*localDomain;
    }();
    double area=0.0;
    for(const auto& e : fe->boundary<ExternalBoundary,4,GaussLegendre>())
    {
        area+=e.first->simplex.child(e.second).vol0;
    }
    require(std::abs(boundary->globalVector().sum()-area)<1.e-10*area,"Boundary assembly changed");

    FixedDirichletSolver<Eigen::SimplicialLLT<Eigen::SparseMatrix<double>>,Eigen::ConjugateGradient<Eigen::SparseMatrix<double>>> solver(true,1.e-12);
    u->dirichletConditions()[0]=2.0;
    solver.compute(mass);
    Eigen::VectorXd desired=Eigen::VectorXd::Constant(u->gSize(),2.0);
    require((solver.solve(M*desired)-desired).norm()<1.e-10,"Fixed Dirichlet solve changed");
    // Constraints are intentionally live after compute; their set stays fixed.
    u->dirichletConditions()[0]=5.0;
    desired.setConstant(5.0);
    require((solver.solve(M*desired)-desired).norm()<1.e-10,"Solver lost instance-owned constraints");
    require(!fe->periodicNodePairs().empty(),"Periodic node pairing lost");
}

void lifetime(const std::filesystem::path& root)
{
    std::weak_ptr<const Mesh> weakMesh;
    std::weak_ptr<Space> weakSpace;
    std::weak_ptr<const void> weakField;
    {
        auto expression=[&]()
        {
            auto mesh=makeMesh(root);
            auto fe=std::make_shared<Space>(mesh);
            auto u=std::make_shared<Field>(*fe);
            *u=6.0;
            weakMesh=mesh;
            weakSpace=fe;
            weakField=u;
            auto intermediate=2.0*u;
            return eval(intermediate);
        }();
        require(!weakField.expired() && !weakSpace.expired() && !weakMesh.expired(),"Expression failed to retain field/space/mesh");
        const auto& e=expression->field().fe().elements().begin()->second;
        require(std::abs((*expression)(e,Eigen::Vector4d::Constant(.25))(0)-12.0)<1.e-11,"Detached expression evaluates incorrectly");
        auto copy=expression;
        require(std::abs((*copy)(e,Eigen::Vector4d::Constant(.25))(0)-12.0)<1.e-11,"Copied evaluator evaluates incorrectly");
    }
    require(weakField.expired() && weakSpace.expired() && weakMesh.expired(),"Shared ownership cycle leaked FEM objects");
    {
        auto retainedDomain=[&]()
        {
            auto fe=std::make_shared<Space>(makeMesh(root));
            weakSpace=fe;
            return fe->domain<EntireDomain,4,GaussLegendre>();
        }();
        require(!weakSpace.expired(),"Domain did not retain its FEM space");
        require(retainedDomain.element(0).simplex.vol0>0,"Retained domain is invalid");
    }
    require(weakSpace.expired(),"Domain leaked its FEM space");
    {
        FixedDirichletSolver<Eigen::SimplicialLLT<Eigen::SparseMatrix<double>>,Eigen::ConjugateGradient<Eigen::SparseMatrix<double>>> solver(true,1.e-12);
        Eigen::VectorXd rhs;
        {
            auto fe=std::make_shared<Space>(makeMesh(root));
            auto u=std::make_shared<Field>(*fe);
            *u=3.0;
            u->dirichletConditions()[0]=3.0;
            weakField=u->ownership();
            weakSpace=fe;
            const auto form=(test(u),u)*fe->domain<EntireDomain,4,GaussLegendre>();
            rhs=matrix(form)*u->dofVector();
            solver.compute(form);
        }
        require(!weakField.expired() && !weakSpace.expired(),"Solver failed to retain field state");
        require((solver.solve(rhs).array()-3.0).matrix().norm()<1.e-10,"Detached solver lost constraints or DOFs");
    }
    require(weakField.expired() && weakSpace.expired(),"Solver leaked its field state");
    auto mesh=makeMesh(root);
    auto expired=[&]()
    {
        Space stackSpace(*mesh);
        auto u=std::make_shared<Field>(stackSpace);
        return eval(u);
    }();
    rejects([&]() { (*expired)(Eigen::Vector3d::Zero()); },"Expired borrowed FEM space was dereferenced");
    auto borrowedMesh=[&]()
    {
        auto temporaryMesh=makeMesh(root);
        auto fe=std::make_shared<Space>(*temporaryMesh);
        auto u=std::make_shared<Field>(*fe);
        return eval(u);
    }();
    rejects([&]() { (*borrowedMesh)(Eigen::Vector3d::Zero()); },"Expired borrowed mesh was dereferenced");
}

template<int order>
void affine(const std::filesystem::path& root)
{
    using FE=FiniteElement<LagrangeElement<3,order>>;
    auto fe=std::make_shared<FE>(makeMesh(root));
    auto u=fe->template trial<3>();
    const Eigen::Matrix3d A=(Eigen::Matrix3d()<<1,2,3,4,5,6,7,8,9).finished();
    for(size_t n=0;n<fe->nodeSize();++n)
    {
        u->dofVector().template segment<3>(3*n)=A*fe->node(n).P0+Eigen::Vector3d(2,3,4);
    }
    auto gradient=eval(grad(u));
    auto sumGradient=eval(grad(u+u));
    auto scaledGradient=eval(grad(2.0*u));
    auto sumStrain=eval(def(u+u));
    auto strain=eval(def(u));
    Eigen::Matrix<double,6,1> expectedStrain;
    expectedStrain<<1,5,9,6,14,10;
    for(const auto& entry : fe->elements())
    {
        const auto& e=entry.second;
        const Eigen::Vector4d bary=Eigen::Vector4d::Constant(.25);
        require(((*eval(u))(e,bary)-A*e.position(bary)-Eigen::Vector3d(2,3,4)).norm()<1.e-10,"Affine interpolation changed");
        const Eigen::Matrix<double,9,1> expected=A.transpose().reshaped();
        require(((*gradient)(e,bary)-expected).norm()<1.e-10,"Affine gradient changed");
        require(((*sumGradient)(e,bary)-2*expected).norm()<1.e-10,"Gradient of trial sum changed");
        require(((*scaledGradient)(e,bary)-2*expected).norm()<1.e-10,"Gradient of scaled trial changed");
        require(((*sumStrain)(e,bary)-2*expectedStrain).norm()<1.e-10,"Strain of trial sum changed");
        require(((*strain)(e,bary)-expectedStrain).norm()<1.e-10,"Engineering strain changed");
    }
    auto dV=fe->template domain<EntireDomain,4,GaussLegendre>();
    const auto energy=matrix((test(def(u)),def(u))*dV);
    const double integrated=u->dofVector().dot(energy*u->dofVector());
    const double exact=fe->mesh.volume()*expectedStrain.squaredNorm();
    require(std::abs(integrated-exact)<1.e-10*exact,"Elastic strain energy changed");
}
int main(int argc,char** argv)
{
    try
    {
        require(argc==2,"Expected repository root");
        independenceAndAssembly(argv[1]);
        lifetime(argv[1]);
        affine<1>(argv[1]);
        affine<2>(argv[1]);
        affine<3>(argv[1]);
        std::cout<<"PASS FEM field ownership, lifetimes, mass/diffusion/load/elastic assembly and constraints"<<std::endl;
    }
    catch(const std::exception& e)
    {
        std::cerr<<e.what()<<std::endl;
        return 1;
    }
}
