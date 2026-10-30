/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * Additional contributors:
 *   Nicholas Huebner Julian <njulian@lanl.gov>
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_DefectiveCrystalBase_H_
#define model_DefectiveCrystalBase_H_

#include <FMM/FMMContext.h>
#include <string>
#include <memory>
#include <vector>
#include <Eigen/Dense>

#include <DefectiveCrystal/DefectiveCrystalParameters.h>
#include <SimplicialMesh.h>
#include <Polycrystal.h>
#include <FiniteElement.h>
#include <Math/VoigtTraits.h>

namespace model
{
    template<int dim> struct DislocationDynamicsBase;
    template<int dim> struct ClusterDynamicsBase;
    template<int dim> struct ElasticDeformationBase;

    template <int _dim>
    struct DefectiveCrystalBase
    {
        // Shared data for MicrostructureGenerator and DefectiveCrystal.
        static constexpr int dim=_dim;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef DefectiveCrystalBase<dim> DefectiveCrystalBaseType;
        typedef LagrangeElement<3,2> ElementType;
        typedef typename ElementType::BaryType BaryType;

        DefectiveCrystalBase( const std::string& folderName);
        DefectiveCrystalParameters simulationParameters;
        const SymmetricVoigtTraits<dim> voigtTraits;
        const std::shared_ptr<const SimplicialMesh<dim>> meshOwner;
        const SimplicialMesh<dim>& mesh;

        const std::vector<int> periodicImageSize;
        const Eigen::Matrix<double,dim,Eigen::Dynamic> periodicLatticeBasis;
        const Eigen::Matrix<double,dim,Eigen::Dynamic> periodicLatticeReciprocalBasis;
        const std::vector<VectorDim> periodicShifts;
        const std::shared_ptr<const Polycrystal<dim>> polyOwner;
        const Polycrystal<dim>& poly;
        std::shared_ptr<FiniteElement<ElementType>> fe;
        const std::vector<std::string> activePhysics;
        const std::shared_ptr<FMMContext<dim>> fmmContext;
        const std::shared_ptr<DislocationDynamicsBase<dim>> dislocationDynamicsBase;
        const std::shared_ptr<ClusterDynamicsBase<dim>> clusterDynamicsBase;
        const std::shared_ptr<ElasticDeformationBase<dim>> elasticDeformationBase;

        bool hasPhysics(const std::string& name) const;
        DislocationDynamicsBase<dim>& dislocationDynamics() const;
        ClusterDynamicsBase<dim>& clusterDynamics() const;
        ElasticDeformationBase<dim>& elasticDeformation() const;
        DefectiveCrystalBase(const DefectiveCrystalBase&)=delete;
        DefectiveCrystalBase& operator=(const DefectiveCrystalBase&)=delete;

        Eigen::VectorXd periodicCoordinates(const VectorDim& x) const;

    };
} // namespace model
#endif
