/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_TrialDef_H_
#define model_TrialDef_H_

#include <utility>
#include <type_traits>
#include <TypeTraits.h>
#include <Simplex.h>
#include <TrialBase.h>
#include <TrialExpressionBase.h>
#include <memory>
#include <stdexcept>

namespace model
{

    template <typename T>
    struct TrialDef : public TrialExpressionBase<TrialDef<T> >
    {

        typedef typename T::TrialFunctionType TrialFunctionType;

        typedef typename TypeTraits<TrialFunctionType>::ElementType ElementType;
        typedef typename TypeTraits<TrialFunctionType>::BaryType BaryType;
        constexpr static int dim=TypeTraits<TrialFunctionType>::dim;
        constexpr static int dofPerElement=TypeTraits<TrialFunctionType>::dofPerElement;
        constexpr static int rows=(TrialFunctionType::dim*(TrialFunctionType::dim+1))/2;
        constexpr static int cols=T::cols;
        typedef Eigen::Matrix<double,rows,dofPerElement> ShapeFunctionMatrixType;
        typedef Eigen::Matrix<double,rows,1> EvalMatrixType;

        static_assert(TypeTraits<TrialFunctionType>::dim==TypeTraits<TrialFunctionType>::nComponents,"SYMMETRIC GRADIENT (DEF) CAN ONLY BE COMPUTED IF nComponents==dim.");

        const std::shared_ptr<T> op;

        explicit TrialDef(std::shared_ptr<T> x) : op(std::move(x))
        {
            if(!op)
            {
                throw std::invalid_argument("Null FEM expression operand.");
            }
        }

        Eigen::Matrix<double,rows,dofPerElement> sfm(const ElementType& ele,
                                                     const BaryType& bary) const
        {
            return op->sfmDef(ele,bary);
        }

        const auto& field() const
        {
            return op->field();
        }

    };

    template<class T, std::enable_if_t<std::is_base_of_v<TrialExpressionBase<T>,T>,int> = 0>
    auto def(const std::shared_ptr<T>& x)
    {
        return std::make_shared<TrialDef<T>>(x);
    }
}
#endif
