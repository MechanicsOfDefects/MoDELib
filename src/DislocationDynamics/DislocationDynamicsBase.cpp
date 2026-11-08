#include <DislocationDynamics/DislocationDynamicsBase.h>
namespace model
{
    template<int dim>
    DislocationDynamicsBase<dim>::DislocationDynamicsBase(
        std::shared_ptr<const SimplicialMesh<dim>> mesh,
        std::shared_ptr<const Polycrystal<dim>> poly) :
        meshOwner(std::move(mesh)),
        polyOwner(std::move(poly)),
        glidePlaneFactory(*polyOwner),
        periodicGlidePlaneFactory(*polyOwner,glidePlaneFactory)
    {
    }
    template struct DislocationDynamicsBase<3>;
}
