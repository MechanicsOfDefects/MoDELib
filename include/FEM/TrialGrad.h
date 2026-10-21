/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_TrialGrad_H_
#define model_TrialGrad_H_

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
    class TrialGrad :
     public TrialExpressionBase<TrialGrad<T> >
    {

    public:

        typedef typename T::TrialFunctionType TrialFunctionType;

        typedef typename TypeTraits<TrialFunctionType>::ElementType ElementType;
        typedef typename TypeTraits<TrialFunctionType>::BaryType BaryType;
        constexpr static int dim=TypeTraits<TrialFunctionType>::dim;
        constexpr static int dofPerElement=TypeTraits<TrialFunctionType>::dofPerElement;
        constexpr static int rows=TrialFunctionType::dim*T::rows;
        constexpr static int cols=T::cols;
        typedef Eigen::Matrix<double,rows,dofPerElement> ShapeFunctionMatrixType;
        typedef Eigen::Matrix<double,rows,1> EvalMatrixType;

        const std::shared_ptr<T> op;

        explicit TrialGrad(std::shared_ptr<T> exp) : op(std::move(exp))
        {
            if(!op)
            {
                throw std::invalid_argument("Null FEM expression operand.");
            }
        }

        ShapeFunctionMatrixType sfm(const ElementType& ele,
                                    const BaryType& bary) const
        {
            return op->sfmGrad(ele,bary);
        }

        const auto& field() const
        {
            return op->field();
        }

    };

    template<class T, std::enable_if_t<std::is_base_of_v<TrialExpressionBase<T>,T>,int> = 0>
    auto grad(const std::shared_ptr<T>& x)
    {
        return std::make_shared<TrialGrad<T>>(x);
    }
}
#endif
