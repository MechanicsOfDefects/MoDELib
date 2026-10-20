#include <FiniteElement.h>
#include <WeakProblem.h>
#include <Eigen/QR>
#include <filesystem>
#include <iostream>
using namespace model;
void require(bool value,const char* message)
{
    if(!value)
    {
        throw std::runtime_error(message);
    }
}
template<class Trial>
struct Load
{
    std::shared_ptr<Trial> trial;
    Eigen::VectorXd values;
    const auto& field() const
    {
        return trial->field();
    }
    Eigen::VectorXd globalVector() const
    {
        return values;
    }
};
template<int components,int order>
void exercise(const std::filesystem::path& root,int rank,bool clamp)
{
    auto mesh=std::make_shared<SimplicialMesh<3>>((root/"Library/Meshes/bicrystal_48.msh").string(),Eigen::Matrix3d::Identity(),Eigen::Vector3d::Zero(),std::set<int>{-1});
    std::set<int> seeds;
    for(const auto& region : mesh->regions())
    {
        for(const auto& face : region.second->faces())
        {
            const auto& pair=face.second->periodicFacePair;
            if(pair.second && (pair.first.normalized().cross(Eigen::Vector3d::UnitY()).norm()<1.e-8 ||
                (rank==2 && pair.first.normalized().cross(Eigen::Vector3d::UnitX()).norm()<1.e-8)))
            {
                seeds.insert(face.first);
            }
        }
    }
    mesh->identifyParallelFaces(seeds);
    using FE=FiniteElement<LagrangeElement<3,order>>;
    auto fe=std::make_shared<FE>(mesh);
    auto u=fe->template trial<components>();
    using Trial=typename decltype(u)::element_type;
    const auto dV=fe->template domain<EntireDomain,4,GaussLegendre>();
    auto bwf=(test(u),u)*dV+(test(grad(u)),grad(u))*dV;
    Eigen::SparseMatrix<double> K(u->gSize(),u->gSize());
    const auto triplets=bwf->globalTriplets();
    K.setFromTriplets(triplets.begin(),triplets.end());
    const auto& B=mesh->periodicBasis();
    const Eigen::MatrixXd Q=B.householderQr().householderQ()*Eigen::MatrixXd::Identity(3,rank);
    Eigen::Matrix<double,components,3> H,G;
    for(int c=0;c<components;++c)
    {
        H.row(c)<<0.13*(c+1),-0.17*(c+1),0.03*(c+1);
        G.row(c)<<0.2*(c+1),0.4*(c+1),0.3*(c+1);
        if(clamp)
        {
            G.row(c)<<0,0,0.3*(c+1);
        }
    }
    u->setPeriodicGradient(H);
    Eigen::MatrixXd modes=Eigen::MatrixXd::Zero(u->gSize(),components*rank);
    Eigen::VectorXd expected(u->gSize());
    for(size_t n=0;n<fe->nodeSize();++n)
    {
        const Eigen::VectorXd shift=Q.transpose()*fe->periodicShift(n);
        for(int j=0;j<rank;++j)
        {
            modes.block(n*components,j*components,components,components)=shift(j)*Eigen::Matrix<double,components,components>::Identity();
        }
        const Eigen::Vector3d x=fe->node(n).P0-Eigen::Vector3d(0,0,mesh->xMin()(2));
        expected.segment(n*components,components)=G*x;
        if(!clamp)
        {
            expected.segment(n*components,components).array()+=0.7;
        }
        if(clamp && std::abs(fe->node(n).P0(2)-mesh->xMin()(2))<1.e-10)
        {
            for(int c=0;c<components;++c)
            {
                u->dirichletConditions()[n*components+c]=0.0;
            }
        }
    }
    const auto pair=fe->periodicNodePairs().front();
    if(!clamp)
    {
        for(int c=0;c<components;++c)
        {
            u->dirichletConditions()[pair.sink*components+c]=expected(pair.sink*components+c);
        }
    }
    auto rhs=std::make_shared<Load<Trial>>(Load<Trial>{u,K*expected});
    auto problem=weakProblem(bwf,rhs);
    problem->setAffineModes(modes);
    const auto verify=[&]()
    {
        require((u->dofVector()-expected).norm()<1.e-8*std::max(1.0,expected.norm()),"Coupled affine/Dirichlet solution is wrong");
        Eigen::Matrix<double,components,3> gradient=H;
        for(int j=0;j<rank;++j)
        {
            gradient+=problem->affineCoefficients().segment(j*components,components)*Q.col(j).transpose();
        }
        require(((gradient-G)*Q).norm()<1.e-8,"Original affine coefficients were not reconstructed");
        for(const auto& bc : u->dirichletConditions())
        {
            require(std::abs(u->dofVector()(bc.first)-bc.second)<1.e-9,"Prescribed node moved");
        }
    };
    problem->solveWithMasterSlaveConstraints();
    verify();
    if(!clamp)
    {
        // Move the anchor to its periodic partner: the nodal mapping is unchanged,
        // but the affine part changes and must invalidate the factorization.
        u->dirichletConditions().clear();
        for(int c=0;c<components;++c)
        {
            u->dirichletConditions()[pair.source*components+c]=expected(pair.source*components+c);
        }
        require(!problem->solveWithMasterSlaveConstraints().reusedFactorization,"Changed affine anchor reused a stale matrix");
        verify();
        for(int c=0;c<components;++c)
        {
            u->dirichletConditions()[pair.sink*components+c]=expected(pair.sink*components+c);
        }
        problem->solveWithMasterSlaveConstraints();
        verify();
    }
    // A prescribed periodic jump is redundant with the clamps or paired pins.
    // Duplicate it deliberately; consistent redundancy must not cause singularity.
    Eigen::MatrixXd rows=Eigen::MatrixXd::Zero(2,u->gSize());
    size_t source=pair.source,sink=pair.sink;
    if(clamp)
    {
        for(const auto& p : fe->periodicNodePairs())
        {
            if(u->dirichletConditions().count(p.source*components) && u->dirichletConditions().count(p.sink*components))
            {
                source=p.source;
                sink=p.sink;
                break;
            }
        }
    }
    rows(0,source*components)=-1.0;
    rows(0,sink*components)=1.0;
    rows.row(1)=2.0*rows.row(0);
    problem->setLinearConstraints(rows);
    problem->setLinearConstraintValues(rows*expected);
    problem->solveWithMasterSlaveConstraints();
    verify();
    require(problem->solveWithMasterSlaveConstraints().reusedFactorization,"Unchanged coupled constraints did not reuse factorization");
    // New RHS and compatible prescribed values change the particular solution,
    // including fixed affine coefficients, without changing the free subspace.
    expected*=1.4;
    rhs->values=K*expected;
    for(auto& bc : u->dirichletConditions())
    {
        bc.second=expected(bc.first);
    }
    G*=1.4;
    problem->setLinearConstraintValues(rows*expected);
    require(problem->solveWithMasterSlaveConstraints().reusedFactorization,"Compatible value changes refactorized the matrix");
    verify();
    const Eigen::VectorXd saved=u->dofVector();
    const Eigen::VectorXd coefficients=problem->affineCoefficients();
    auto incompatible=(rows*expected).eval();
    incompatible(0)+=1.0;
    incompatible(1)+=2.0;
    problem->setLinearConstraintValues(incompatible);
    bool rejected=false;
    try
    {
        problem->solveWithMasterSlaveConstraints();
    }
    catch(const std::invalid_argument& e)
    {
        rejected=std::string(e.what()).find("Incompatible")!=std::string::npos;
    }
    require(rejected,"Conflicting affine loading and clamps were not diagnosed");
    require((u->dofVector()-saved).norm()==0.0 && (problem->affineCoefficients()-coefficients).norm()==0.0,"Rejected update changed the solution");
    problem->setLinearConstraintValues(rows*expected);
    problem->solveWithMasterSlaveConstraints();
    verify();
    if(clamp)
    {
        // Contradictory jumps along different edges cannot be fitted by one affine mode.
        u->dirichletConditions().at(sink*components)+=1.0;
        rejected=false;
        try
        {
            problem->solveWithMasterSlaveConstraints();
        }
        catch(const std::invalid_argument& e)
        {
            rejected=std::string(e.what()).find("Incompatible")!=std::string::npos;
        }
        require(rejected,"Inconsistent prescribed periodic copies were accepted");
        require((u->dofVector()-saved).norm()<1.e-12,"Incompatible nodal conditions changed field");
    }
    // All nodal and affine unknowns may be prescribed; no factorization is needed.
    u->dirichletConditions().clear();
    for(size_t i=0;i<size_t(expected.size());++i)
    {
        u->dirichletConditions()[i]=expected(i);
    }
    problem->solveWithMasterSlaveConstraints();
    verify();
}
int main(int argc,char** argv)
{
    try
    {
        require(argc==2,"Expected repository root");
        for(int rank : {1,2})
        {
            for(bool clamp : {false,true})
            {
                exercise<1,1>(argv[1],rank,clamp);
                exercise<3,2>(argv[1],rank,clamp);
            }
        }
        std::cout<<"PASS affine modes with pins, clamps, compatibility checks and factorization reuse\n";
    }
    catch(const std::exception& e)
    {
        std::cerr<<e.what()<<'\n';
        return 1;
    }
}
