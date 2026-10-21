/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_EvalExpression_H_
#define model_EvalExpression_H_

#include <utility>
#include <EvalFunction.h>
#include <memory>
#include <stdexcept>

namespace model
{

    template<typename TrialExpressionType>
    struct EvalExpression : public EvalFunction<EvalExpression<TrialExpressionType>>
    {
        typedef typename TrialExpressionType::TrialFunctionType TrialFunctionType;
        constexpr static int rows=TrialExpressionType::rows;
        constexpr static int cols=TrialExpressionType::cols;

        typedef Eigen::Matrix<double,rows,1> EvalMatrixType;

        typedef typename TypeTraits<TrialFunctionType>::ElementType ElementType;
        typedef typename TypeTraits<TrialFunctionType>::BaryType BaryType;
        typedef typename TypeTraits<TrialFunctionType>::NodeType NodeType;
        constexpr static int dim=TypeTraits<TrialFunctionType>::dim;

        typedef TrialBase<TrialFunctionType> TrialBaseType;

        const std::shared_ptr<TrialExpressionType> trialExp;

        explicit EvalExpression(std::shared_ptr<TrialExpressionType> te) : trialExp(std::move(te))
        {
            if(!trialExp)
            {
                throw std::invalid_argument("Null FEM expression operand.");
            }
        }

        EvalMatrixType operator()(const ElementType& ele,
                                  const BaryType& bary) const
        {
            return trialExp->sfm(ele,bary)*field().dofs(ele);
        }

        EvalMatrixType operator()(const Eigen::Matrix<double,dim,1>& P,
                                  const Simplex<dim,dim>* const guess) const
        {
            const std::pair<bool,const ElementType*> temp(field().fe().searchWithGuess(P,guess));
            if(temp.first)
            {
                auto bary(temp.second->simplex.pos2bary(P));
                for(int d=0; d<dim+1; d++)
                {
                    if(bary(d)<0.0)
                    {
                        bary(d)=0.0;
                    }
                    if(bary(d)>1.0)
                    {
                        bary(d)=1.0;
                    }
                }

                return trialExp->sfm(*(temp.second),bary)*field().dofs(*(temp.second));
            }
            else
            {
                std::cout<<"WARNING: EVALUATING EXPRESSION OUTSIDE DOMAIN"<<std::endl;
                return EvalMatrixType::Zero();
            }
        }

        EvalMatrixType operator()(const Eigen::Matrix<double,dim,1>& P) const
        {
            return this->operator()(P,field().mesh().simplices().begin()->second.get());
        }

        EvalMatrixType operator()(const NodeType& node) const
        {
            EvalMatrixType temp(EvalMatrixType::Zero());
            for(const auto& ele : node)
            {
                const Eigen::Matrix<double,dim+1,1> bary(ele->simplex.pos2bary(node.P0));
                temp+=trialExp->sfm(*ele,bary)*field().dofs(*ele);
            }
            return temp/node.size();
        }

        const auto& field() const
        {
            return trialExp->field();
        }

    };

}
#endif
