#include <DislocationDynamics/DislocationDynamicsBase.h>
#include <CTM.h>
namespace model
{
    template<int dim>
    DislocationDynamicsBase<dim>::DislocationDynamicsBase(
        std::shared_ptr<const SimplicialMesh<dim>> mesh,
        std::shared_ptr<const Polycrystal<dim>> poly,
        const Eigen::Matrix<double,dim,Eigen::Dynamic>& basis,
        const std::string& inputFile) :
        meshOwner(std::move(mesh)),
        polyOwner(std::move(poly)),
        glidePlaneFactory(*polyOwner),
        periodicGlidePlaneFactory(*polyOwner,glidePlaneFactory),
        EwaldLength(basis.cols()>0 ? TextFileParser(inputFile).readScalar<double>("EwaldLengthFactor",true)
            *std::pow(std::sqrt((basis.transpose()*basis).determinant())/CTM::factorial(basis.cols()),1.0/basis.cols()) : 0.0)
    {
    }
    template struct DislocationDynamicsBase<3>;
}
