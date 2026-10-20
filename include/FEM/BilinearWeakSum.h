/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_BilinearWeakSum_H_
#define model_BilinearWeakSum_H_

#include <LowRankBilinear.h>
#include <iterator>
#include <type_traits>

#include <BilinearWeakExpression.h>
#include <memory>
#include <stdexcept>

namespace model
{

    template <typename T1,typename T2>
	struct  BilinearWeakSum : public BilinearWeakExpression<BilinearWeakSum<T1,T2> >
    {

        static_assert(std::is_same<typename T1::TrialFunctionType,typename T2::TrialFunctionType>::value,"YOU ARE SUMMING BilinearWeakForms OF DIFFERENT TRIALFUNCTIONS.");
        typedef typename T1::TrialFunctionType TrialFunctionType;

        const std::shared_ptr<T1> op1;

        const std::shared_ptr<T2> op2;
        BilinearWeakSum(std::shared_ptr<T1> x, std::shared_ptr<T2> y) :
            op1(std::move(x)), op2(std::move(y))
        {
            if(!op1 || !op2)
            {
                throw std::invalid_argument("Null FEM expression operand.");
            }
        }

        std::vector<Eigen::Triplet<double> >  globalTriplets() const
        {
            if(!(op2->field().sameField(field())))
            {
                throw std::runtime_error("Cannot assemble expressions from different FEM field instances.");
            }
            std::vector<Eigen::Triplet<double> > gt1=op1->globalTriplets();
            std::vector<Eigen::Triplet<double> > gt2=op2->globalTriplets();
            gt1.insert(
                        gt1.end(),
                        std::make_move_iterator(gt2.begin()),
                        std::make_move_iterator(gt2.end())
                        );
            return gt1;
        }

        Eigen::MatrixXd lowRankPenaltyRows() const
        {
            if(!op2->field().sameField(field()))
            {
                throw std::runtime_error("Cannot combine different FEM fields.");
            }
            const auto a=lowRankBilinearRows(*op1);
            const auto b=lowRankBilinearRows(*op2);
            Eigen::MatrixXd result(a.rows()+b.rows(),field().gSize());
            result.topRows(a.rows())=a;
            result.bottomRows(b.rows())=b;
            return result;
        }

        const auto& field() const
        {
            return op1->field();
        }

    };

    template<class T1,class T2, std::enable_if_t<std::is_base_of_v<BilinearWeakExpression<T1>,T1> && std::is_base_of_v<BilinearWeakExpression<T2>,T2>,int> = 0>
    auto operator+(const std::shared_ptr<T1>& x,const std::shared_ptr<T2>& y)
    {
        return std::make_shared<BilinearWeakSum<T1,T2>>(x,y);
    }

}
#endif
