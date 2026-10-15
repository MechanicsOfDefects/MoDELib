/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LoopSinkReaction_H_
#define model_LoopSinkReaction_H_

#include <Eigen/Dense>
#include <EvalFunction.h>
#include <EvalExpression.h>
#include <ImmobileRateEquations.h>

namespace model
{

    /*!\brief What the immobile clusters contribute to the equations of the mobile
     * species, evaluated from ImmobileRateEquations so that the two sets of
     * equations exchange exactly the same defects.
     *
     * LoopSinkReaction is the first-order reaction matrix of the absorption at
     * the immobile clusters. It is diagonal: the entry of mobile species m is
     * minus the sum over the families of the coefficients A_km of
     * ImmobileRateEquations::absorptionCoefficients.
     */
    template <typename MobileTrialFunctionType, typename ImmobileTrialFunctionType>
    struct LoopSinkReaction : public EvalFunction<LoopSinkReaction<MobileTrialFunctionType,ImmobileTrialFunctionType>>
    {
        constexpr static int rows=MobileTrialFunctionType::rows;
        constexpr static int cols=rows;
        constexpr static int iSize=ImmobileTrialFunctionType::rows;
        constexpr static int dim=MobileTrialFunctionType::dim;
        typedef ImmobileRateEquations<dim> RateEquationsType;

        const RateEquationsType& rateEquations;
        const EvalExpression<ImmobileTrialFunctionType> immobileEval;

        LoopSinkReaction(const ImmobileTrialFunctionType& immobile,const RateEquationsType& rateEquations_in) :
        /* init */ rateEquations(rateEquations_in)
        /* init */,immobileEval(immobile)
        {

        }

        template<typename ElementType, typename BaryType>
        const Eigen::Matrix<double,rows,cols> operator() (const ElementType& ele, const BaryType& bary) const
        {
            const Eigen::Matrix<double,iSize,1> value(immobileEval(ele,bary));
            // the interpolated fields can be slightly negative between nodes
            const typename RateEquationsType::ArrayF n(value.template segment<iSize/2>(0).transpose().array().max(0.0)); // clusters per atom
            const typename RateEquationsType::ArrayF c(value.template segment<iSize/2>(iSize/2).transpose().array().max(0.0));
            Eigen::Matrix<double,rows,cols> temp(Eigen::Matrix<double,rows,cols>::Zero());
//          temp.diagonal()=-rateEquations.absorptionCoefficients(n,c,rateEquations.cdp.effectiveClusterRadius(n,c)).colwise().sum().matrix().transpose();
            temp.diagonal()=-rateEquations.absorptionCoefficients(n,c,rateEquations.cdp.captureLength(n,c)).colwise().sum().matrix().transpose();
            return temp;
        }
    };

    /*!\brief The defects that the immobile clusters return to the mobile species (Eq. 55):
     * vacancies from the dissolution and thermal emission of the vacancy clusters,
     * interstitials from the dissolution of the interstitial clusters. A column with
     * one entry per mobile species; only the entries of the single vacancy and of
     * the single interstitial are not zero.
     */
    template <typename MobileTrialFunctionType, typename ImmobileTrialFunctionType>
    struct ImmobileClusterSource : public EvalFunction<ImmobileClusterSource<MobileTrialFunctionType,ImmobileTrialFunctionType>>
    {
        constexpr static int rows=MobileTrialFunctionType::rows;
        constexpr static int cols=1;
        constexpr static int iSize=ImmobileTrialFunctionType::rows;
        constexpr static int dim=MobileTrialFunctionType::dim;
        typedef ImmobileRateEquations<dim> RateEquationsType;

        const RateEquationsType& rateEquations;
        const EvalExpression<ImmobileTrialFunctionType> immobileEval;

        ImmobileClusterSource(const ImmobileTrialFunctionType& immobile,const RateEquationsType& rateEquations_in) :
        /* init */ rateEquations(rateEquations_in)
        /* init */,immobileEval(immobile)
        {

        }

        template<typename ElementType, typename BaryType>
        const Eigen::Matrix<double,rows,cols> operator() (const ElementType& ele, const BaryType& bary) const
        {
            Eigen::Matrix<double,rows,cols> temp(Eigen::Matrix<double,rows,cols>::Zero());

            const Eigen::Matrix<double,iSize,1> value(immobileEval(ele,bary));
            const typename RateEquationsType::ArrayF n(value.template segment<iSize/2>(0).transpose().array().max(0.0)); // clusters per atom
            const typename RateEquationsType::ArrayF c(value.template segment<iSize/2>(iSize/2).transpose().array().max(0.0));

            // Eq. (55): the content removed by annealing and thermal emission returns to the single vacancy and the single interstitial
//          const typename RateEquationsType::ArrayF released(rateEquations.Ac(n,c)+rateEquations.Aemit(n,c,rateEquations.cdp.effectiveClusterRadius(n,c)));
            const typename RateEquationsType::ArrayF released(rateEquations.Ac(n,c)+rateEquations.Aemit(n,c,rateEquations.cdp.captureLength(n,c)));
            if(rateEquations.cdp.vacancyIndex>=0)
            {
                temp(rateEquations.cdp.vacancyIndex)=released.matrix().dot(rateEquations.cdp.Iv.matrix());
            }
            if(rateEquations.cdp.interstitialIndex>=0)
            {
                temp(rateEquations.cdp.interstitialIndex)=released.matrix().dot(rateEquations.cdp.Ii.matrix());
            }

            return temp;
        }
    };

}
#endif
