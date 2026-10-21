/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_BilinearWeakForm_H_
#define model_BilinearWeakForm_H_

#include <chrono>
#include <vector>

#include <Eigen/Sparse>

#include <TrialBase.h>
#include <LinearWeakForm.h>
#include <TerminalColors.h>
#include <BilinearForm.h>

#include <BilinearWeakExpression.h>
#include <BilinearWeakSum.h>

namespace model
{

    template <typename _BilinearFormType,typename _IntegrationDomainType>
    struct BilinearWeakForm : public BilinearWeakExpression<BilinearWeakForm<_BilinearFormType,_IntegrationDomainType> >
    {

        typedef _BilinearFormType BilinearFormType;
        typedef _IntegrationDomainType IntegrationDomainType;
        typedef BilinearWeakForm<BilinearFormType,IntegrationDomainType> BilinearWeakFormType;
        typedef typename BilinearFormType::TrialFunctionType TrialFunctionType;
        typedef typename BilinearFormType::TestExpressionType TestExpressionType;
        typedef typename BilinearFormType::TrialExpressionType TrialExpressionType;
        typedef typename TrialFunctionType::ElementType ElementType;
        typedef typename IntegrationDomainType::QuadratureType QuadratureType;
        typedef typename QuadratureType::VectorDim AbscissaType;

        constexpr static int dim=TypeTraits<TrialFunctionType>::dim;
        constexpr static int nodesPerElement=TypeTraits<TrialFunctionType>::nodesPerElement;
        constexpr static int dofPerNode=TypeTraits<TrialFunctionType>::dofPerNode;
        constexpr static int dofPerElement=TypeTraits<TrialFunctionType>::dofPerElement;

        typedef Eigen::Matrix<double,dofPerElement,dofPerElement> ElementMatrixType;

        const std::shared_ptr<BilinearFormType> bilinearForm;

        const IntegrationDomainType domain;
        BilinearWeakForm(std::shared_ptr<BilinearFormType> f,const IntegrationDomainType& dom) : bilinearForm(std::move(f)), domain(dom)
        {
            if(!bilinearForm)
            {
                throw std::invalid_argument("Null FEM form.");
            }
            domain.validateSpace(field().fe());
        }




        std::vector<Eigen::Triplet<double> >  globalTriplets() const
        {
            domain.validateSpace(field().fe());
            std::vector<Eigen::Triplet<double> > glbtrip;
            glbtrip.reserve(dofPerElement*dofPerElement*field().elementSize());

            for (size_t k=0;k<domain.size();++k)
            {
                const ElementType& ele(domain.element(k));
                ElementMatrixType ke(ElementMatrixType::Zero());
                QuadratureType::integrate(this,ke,&BilinearWeakFormType::elementAssemblyKernel<AbscissaType>,ele);

                for (int i=0;i<dofPerElement;++i)
                {
                    for (int j=0;j<dofPerElement;++j)
                    {
                        if(ke(i,j)!=0.0)
                        {
                            const size_t  nodeID_I(i/dofPerNode);
                            const size_t nodeDof_I(i%dofPerNode);
                            const size_t gI= ele.node(nodeID_I).gID*dofPerNode+nodeDof_I;

                            const size_t  nodeID_J(j/dofPerNode);
                            const size_t nodeDof_J(j%dofPerNode);
                            const size_t gJ=ele.node(nodeID_J).gID*dofPerNode+nodeDof_J;

                            glbtrip.emplace_back(gI,gJ,ke(i,j));
                        }
                    }
                }
            }
            return glbtrip;
        }

        template <typename AbscissaType>
        ElementMatrixType elementAssemblyKernel(const AbscissaType& a, const ElementType& ele) const
        {
            const Eigen::Matrix<double,dim+1,1> bary(BarycentricTraits<dim>::x2l(a));

            return (*bilinearForm)(ele,bary)*ele.absJ(bary);
        }

        const auto& field() const
        {
            return bilinearForm->field();
        }

    };

    template<class T1,class T2,class FE,int D,int Q,template<short unsigned int,size_t> class Rule>
    auto operator*(const std::shared_ptr<BilinearForm<T1,T2>>& f,const IntegrationDomain<FE,D,Q,Rule>& dom)
    {
        return std::make_shared<BilinearWeakForm<BilinearForm<T1,T2>,IntegrationDomain<FE,D,Q,Rule>>>(f,dom);
    }

}
#endif
