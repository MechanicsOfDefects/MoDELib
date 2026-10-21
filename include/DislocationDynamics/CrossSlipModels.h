/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_CrossSlipModels_H
#define model_CrossSlipModels_H

#include <deque>
#include <CrossSlipBranch.h>
#include <map>
#include <tuple>

#include <Eigen/Dense>

#include <TypeTraits.h>
#include <LatticeModule.h>
#include <BCClattice.h>
#include <FCClattice.h>
#include <Grain.h>
#include <TerminalColors.h>
#include <DefectiveCrystal/DefectiveCrystalTraitsIO.h>

namespace model
{
    template<int dim> class DislocationNetwork;
    template <int dim>
    struct BaseCrossSlipModel
    {
        typedef TypeTraits<DislocationNetwork<dim>> TraitsType;
        typedef typename TraitsType::NetworkLinkType NetworkLinkType;
        typedef typename TraitsType::NetworkNodeType NetworkNodeType;
        typedef std::tuple<std::shared_ptr<NetworkNodeType>,std::shared_ptr<NetworkNodeType>,size_t,size_t> CrossSlipTupleType;
        typedef std::deque<CrossSlipTupleType> CrossSlipContainerType;
        
        
        BaseCrossSlipModel(const DefectiveCrystalTraitsIO& traitsIO);
        virtual ~BaseCrossSlipModel() = default;

        virtual void addToCrossSlip(const NetworkLinkType& link,CrossSlipContainerType& crossSlipDeq) = 0;
        virtual std::pair<bool,std::pair<int,int>> isCrossSlipLink(const NetworkLinkType& link) = 0;
        std::pair<bool,int> isBaseCrossSlipLink(const NetworkLinkType& link, const bool enforceRemeshLength=true) const;
        virtual bool isThermal() const { return false; }
        virtual std::optional<ThermalCrossSlipEvent> sampleBranch(const std::vector<const NetworkLinkType*>&,const double) { return std::nullopt; }
        
        const double crossSlipDeg;
        const double sinCrossSlip;
    };

    template <int dim>
    struct AthermalCrossSlipModel : public BaseCrossSlipModel<dim>
    {
        typedef TypeTraits<DislocationNetwork<dim>> TraitsType;
        typedef typename TraitsType::NetworkLinkType NetworkLinkType;
        typedef typename TraitsType::NetworkNodeType NetworkNodeType;
        typedef std::tuple<std::shared_ptr<NetworkNodeType>,std::shared_ptr<NetworkNodeType>,size_t,size_t> CrossSlipTupleType;
        typedef std::deque<CrossSlipTupleType> CrossSlipContainerType;
        
        AthermalCrossSlipModel(const DefectiveCrystalTraitsIO& traitsIO);
        void addToCrossSlip(const NetworkLinkType& link,CrossSlipContainerType& crossSlipDeq) override;
        std::pair<bool,std::pair<int,int>> isCrossSlipLink(const NetworkLinkType& link) override;
    };

//    template <int dim>   //For FCC lattice the same is BCC
//    struct EscaigCrossSlipModelFCC : public BaseCrossSlipModel<dim>
//    {
//        typedef TypeTraits<DislocationNetwork<dim>> TraitsType;
//        static constexpr int dim=TraitsType::dim;
//        typedef typename TraitsType::NetworkLinkType NetworkLinkType;
//        typedef typename TraitsType::NetworkNodeType NetworkNodeType;
//        typedef std::tuple<std::shared_ptr<NetworkNodeType>,std::shared_ptr<NetworkNodeType>,size_t,size_t> CrossSlipTupleType;
//        typedef std::deque<CrossSlipTupleType> CrossSlipContainerType;
//
//        std::random_device rd;
//        std::mt19937 gen;
//        std::uniform_real_distribution<> dist;
//        // const bool enableBasalSlipSystems;
//        const double kB_DD;
//        const double lengthC;
//        const double activeVol;
//        const double Tcs;
//        const double omegaA;
//        const double Ea;
//        const double Tau_3;
//
//       
//        EscaigCrossSlipModelFCC(const PolycrystallineMaterialBase& material,const DefectiveCrystalTraitsIO& traitsIO);
//        double crossSlipRateKernel(const int& k,const NetworkLinkType& link,const SlipSystem& css,
//                                   const double& w,
//                                   const double& L,
//                                   const double& dE,
//                                   const double& Veff
//                                   ) const;
//        void addToCrossSlip(const NetworkLinkType& link,CrossSlipContainerType& crossSlipDeq) override;
//        std::pair<bool,std::pair<int,int>> isCrossSlipLink(const NetworkLinkType& link) override;
//    };

    template <int dim>
    struct EscaigCrossSlipModelHEX : public BaseCrossSlipModel<dim>
    {
        typedef TypeTraits<DislocationNetwork<dim>> TraitsType;
        typedef typename TraitsType::NetworkLinkType NetworkLinkType;
        typedef typename TraitsType::NetworkNodeType NetworkNodeType;
        typedef std::tuple<std::shared_ptr<NetworkNodeType>,std::shared_ptr<NetworkNodeType>,size_t,size_t> CrossSlipTupleType;
        typedef std::deque<CrossSlipTupleType> CrossSlipContainerType;
        
        std::random_device rd;
        std::mt19937 gen;
        std::uniform_real_distribution<> dist;
        const bool enableBasalSlipSystems;
        const bool enablePrismaticSlipSystems;
        const bool enablePyramidalSlipSystems;
        const double csw;
        const double csLb;
        const double csLp;
        const double csb2pE;
        const double csb2pVeg;
        const double csb2pVsc;
        const double csb2pVec;
        const double csp2bE;
        const double csp2bVeg;
        const double csp2bVsc;
        const double csp2bVec;
        
        EscaigCrossSlipModelHEX(const PolycrystallineMaterialBase& material,const DefectiveCrystalTraitsIO& traitsIO);
        bool isThermal() const override { return true; }
        std::optional<ThermalCrossSlipEvent> sampleBranch(const std::vector<const NetworkLinkType*>& links,const double dt) override;
        std::map<int,double> targetRates(const NetworkLinkType& link) const;
        double nucleationLength(const NetworkLinkType& link) const;
        double crossSlipRateKernel(const int& k,const NetworkLinkType& link,const SlipSystem& css,
                                   const double& w,
                                   const double& L,
                                   const double& dE,
                                   const double& Veg,
                                   const double& Vsc,
                                   const double& Vec) const;
        void addToCrossSlip(const NetworkLinkType& link,CrossSlipContainerType& crossSlipDeq) override;
        std::pair<bool,std::pair<int,int>> isCrossSlipLink(const NetworkLinkType& link) override;
    };
}
#endif
