#ifndef model_DislocationDynamicsBase_H_
#define model_DislocationDynamicsBase_H_
#include <GlidePlaneFactory.h>
#include <PeriodicGlidePlaneFactory.h>
namespace model
{
    template<int dim>
    struct DislocationDynamicsBase
    {
        const std::shared_ptr<const SimplicialMesh<dim>> meshOwner;
        const std::shared_ptr<const Polycrystal<dim>> polyOwner;
        GlidePlaneFactory<dim> glidePlaneFactory;
        PeriodicGlidePlaneFactory<dim> periodicGlidePlaneFactory;
        const double EwaldLength;
        DislocationDynamicsBase(std::shared_ptr<const SimplicialMesh<dim>> mesh,
                                std::shared_ptr<const Polycrystal<dim>> poly,
                                const Eigen::Matrix<double,dim,Eigen::Dynamic>& basis,
                                const std::string& inputFile);
        DislocationDynamicsBase(const DislocationDynamicsBase&)=delete;
        DislocationDynamicsBase& operator=(const DislocationDynamicsBase&)=delete;
    };
}
#endif
