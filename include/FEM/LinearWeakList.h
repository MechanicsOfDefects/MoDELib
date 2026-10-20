/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LinearWeakList_H_
#define model_LinearWeakList_H_

#include <TypeTraits.h>
#include <TerminalColors.h>
#include <TestExpression.h>
#include <LinearWeakExpression.h>
#include <JGNselector.h>
#include <IntegrationList.h>
#include <memory>
#include <stdexcept>

namespace model
{

    template <typename TestType, typename _IntegrationListType>
    class LinearWeakList : public LinearWeakExpression<LinearWeakList<TestType,_IntegrationListType> >
    {

    public:

        typedef typename TestType::TrialFunctionType TrialFunctionType;

        constexpr static int dim=TrialFunctionType::dim;

    private:
        typedef typename TypeTraits<TrialFunctionType>::ElementType ElementType;
        constexpr static int dofPerNode=TypeTraits<TrialFunctionType>::dofPerNode;
        constexpr static int dofPerElement=TypeTraits<TrialFunctionType>::dofPerElement;
        typedef Eigen::Matrix<double,dofPerElement,1> ElementVectorType;

    public:

        const std::shared_ptr<TestExpression<TestType>> testExp;
        const _IntegrationListType list;
        LinearWeakList(std::shared_ptr<TestExpression<TestType>> f,const _IntegrationListType& dom) : testExp(std::move(f)), list(dom)
        {
            if(!testExp)
            {
                throw std::invalid_argument("Null FEM form.");
            }
        }

        Eigen::Matrix<double,Eigen::Dynamic,1> globalVector() const
        {
            Eigen::Matrix<double,Eigen::Dynamic,1> _globalVector(Eigen::Matrix<double,Eigen::Dynamic,1>::Zero(field().gSize()));

            for (size_t k=0;k<list.size();++k)
            {
                const ElementType& ele(list[k].ele);
                const int f(list[k].boundaryFace);
                const Eigen::Matrix<double,dim+1,1>& bary(list[k].domainBary);

                const ElementVectorType ve(testExp->sfm(ele,bary).transpose()
                                           *list[k]
                                           *JGNselector<3>::jGN(ele.jGN(bary,f))
                                           *list[k].weight
                                           );

                for (int i=0;i<dofPerElement;++i)
                {
                    const size_t  nodeID_I(i/dofPerNode);
                    const size_t nodeDof_I(i%dofPerNode);
                    const size_t gI= ele.node(nodeID_I).gID*dofPerNode+nodeDof_I;
                    _globalVector(gI) += ve(i);
                }
            }

            return _globalVector;
        }

        const auto& field() const
        {
            return testExp->field();
        }

    };

    template<class T,int D,class P>
    auto operator,(const std::shared_ptr<TestExpression<T>>& t,const IntegrationList<D,P>& list)
    {
        return std::make_shared<LinearWeakList<T,IntegrationList<D,P>>>(t,list);
    }

}
#endif
