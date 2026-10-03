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
            const typename RateEquationsType::ArrayF n(value.template segment<iSize/2>(0).transpose().array().max(0.0)*rateEquations.cdp.omega);
            const typename RateEquationsType::ArrayF c(value.template segment<iSize/2>(iSize/2).transpose().array().max(0.0));
            Eigen::Matrix<double,rows,cols> temp(Eigen::Matrix<double,rows,cols>::Zero());
            temp.diagonal()=-rateEquations.absorptionCoefficients(n,c).colwise().sum().matrix().transpose();
            return temp;
        }
    };

    /*!\brief The vacancies that the vacancy clusters return to the mobile
     * species by dissolution and by thermal emission. A column with one entry per
     * mobile species; only the entry of the single vacancy is not zero.
     */
    template <typename MobileTrialFunctionType, typename ImmobileTrialFunctionType>
    struct VacancyLoopSource : public EvalFunction<VacancyLoopSource<MobileTrialFunctionType,ImmobileTrialFunctionType>>
    {
        constexpr static int rows=MobileTrialFunctionType::rows;
        constexpr static int cols=1;
        constexpr static int iSize=ImmobileTrialFunctionType::rows;
        constexpr static int dim=MobileTrialFunctionType::dim;
        typedef ImmobileRateEquations<dim> RateEquationsType;

        const RateEquationsType& rateEquations;
        const EvalExpression<ImmobileTrialFunctionType> immobileEval;

        VacancyLoopSource(const ImmobileTrialFunctionType& immobile,const RateEquationsType& rateEquations_in) :
        /* init */ rateEquations(rateEquations_in)
        /* init */,immobileEval(immobile)
        {

        }

        template<typename ElementType, typename BaryType>
        const Eigen::Matrix<double,rows,cols> operator() (const ElementType& ele, const BaryType& bary) const
        {
            Eigen::Matrix<double,rows,cols> temp(Eigen::Matrix<double,rows,cols>::Zero());
            if(rateEquations.vacancyIndex>=0)
            {
                const Eigen::Matrix<double,iSize,1> value(immobileEval(ele,bary));
                const typename RateEquationsType::ArrayF n(value.template segment<iSize/2>(0).transpose().array().max(0.0)*rateEquations.cdp.omega);
                const typename RateEquationsType::ArrayF c(value.template segment<iSize/2>(iSize/2).transpose().array().max(0.0));
                temp(rateEquations.vacancyIndex)=rateEquations.vacancySource(n,c);
            }
            return temp;
        }
    };

}
#endif
