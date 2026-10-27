#ifndef model_DislocationDynamicsBase_H_
#define model_DislocationDynamicsBase_H_
#include <GlidePlaneFactory.h>
#include <FMM/BoxTree.h>
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
        const unsigned fmmTreeDepth;
        // Geometry only; retained across frames and shared with the viewer.
        const std::shared_ptr<BoxTree<dim>> boxTree;
        DislocationDynamicsBase(std::shared_ptr<const SimplicialMesh<dim>> mesh,
                                std::shared_ptr<const Polycrystal<dim>> poly,
                                const Eigen::Matrix<double,dim,Eigen::Dynamic>& basis,
                                const std::string& inputFile,
                                const std::string& ddInputFile);
        DislocationDynamicsBase(const DislocationDynamicsBase&)=delete;
        DislocationDynamicsBase& operator=(const DislocationDynamicsBase&)=delete;
    };
}
#endif
