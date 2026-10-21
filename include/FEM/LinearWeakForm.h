/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LinearWeakForm_H_
#define model_LinearWeakForm_H_

#include <chrono>
#include <vector>

#include <Eigen/Dense>

#include <TerminalColors.h>
#include <IntegrationDomain.h>
#include <LinearWeakExpression.h>
#include <LinearForm.h>
#include <LinearWeakSum.h>
#include <LinearWeakDiff.h>
#include <LinearWeakAssembler.h>

#include <memory>
#include <stdexcept>

namespace model
{

    template <typename _LinearFormType, typename _IntegrationDomainType>
	struct LinearWeakForm : public LinearWeakExpression<LinearWeakForm<_LinearFormType,_IntegrationDomainType> >
    {

        typedef _LinearFormType LinearFormType;
        typedef _IntegrationDomainType IntegrationDomainType;
        typedef LinearWeakForm<LinearFormType,IntegrationDomainType> LinearWeakFormType;

        typedef typename LinearFormType::TrialFunctionType TrialFunctionType;
        constexpr static int dim=IntegrationDomainType::dim;
        static_assert(dim==TypeTraits<TrialFunctionType>::dim, "dim MUST BE THE SAME");
        constexpr static int domainDim=IntegrationDomainType::domainDim;
        static_assert(domainDim==dim || domainDim==dim-1, "DOMAIN DIMENSIONALITY MUST BE EITHER dim (volume integration) OR dim-1 (boundary integration)");

        const std::shared_ptr<LinearFormType> linearForm;

        const IntegrationDomainType domain;
        LinearWeakForm(std::shared_ptr<LinearFormType> f,const IntegrationDomainType& dom) : linearForm(std::move(f)), domain(dom)
        {
            if(!linearForm)
            {
                throw std::invalid_argument("Null FEM form.");
            }
            domain.validateSpace(field().fe());
        }




        Eigen::Matrix<double,Eigen::Dynamic,1> globalVector() const
        {
            domain.validateSpace(field().fe());
            return LinearWeakAssembler<LinearWeakFormType,dim-domainDim>(*this).assemble();
        }

        const auto& field() const
        {
            return linearForm->field();
        }

    };

    template<class T1,class T2,class FE,int D,int Q,template<short unsigned int,size_t> class Rule>
    auto operator*(const std::shared_ptr<LinearForm<T1,T2>>& f,const IntegrationDomain<FE,D,Q,Rule>& dom)
    {
        return std::make_shared<LinearWeakForm<LinearForm<T1,T2>,IntegrationDomain<FE,D,Q,Rule>>>(f,dom);
    }

}
#endif
