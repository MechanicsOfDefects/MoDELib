#include <cfloat>
#include <FiniteElement.h>
#include <WeakProblem.h>
#include <cmath>
using namespace model;

Eigen::VectorXd solveElasticity(std::shared_ptr<SimplicialMesh<3>> mesh)
{
    using FE = FiniteElement<LagrangeElement<3,2>>;
    auto fe = std::make_shared<FE>(mesh);
    auto u = fe->trial<3>();                 // displacement
    const double mu = 1.0, nu = 0.3;
    const double lambda = 2*mu*nu/(1-2*nu);
    Eigen::Matrix<double,6,6> C = Eigen::Matrix<double,6,6>::Zero();
    C.topLeftCorner<3,3>().setConstant(lambda);
    C.diagonal().head<3>().array() += 2*mu;
    C.diagonal().tail<3>().setConstant(mu); // engineering shear strain

    for (size_t n=0; n<fe->nodeSize(); ++n)
        if (std::abs(fe->node(n).P0(0)-mesh->xMin()(0)) < 1.e-10)
            for (int c=0; c<3; ++c)
                u->dirichletConditions()[3*n+c] = 0.0; // clamp x-min face

    auto dV = fe->domain<EntireDomain,4,GaussLegendre>();
    auto e = def(u);
    auto stiffness = (test(e), C*e)*dV;
    const Eigen::Vector3d force(0.0, 0.0, -1.e-3); // body force per volume
    auto load = (test(u), make_constant(force))*dV;
    auto problem = weakProblem(stiffness, load);
    problem->solveWithMasterSlaveConstraints(); // assemble, constrain, solve
    return u->dofVector();                     // three entries per FE node
}
