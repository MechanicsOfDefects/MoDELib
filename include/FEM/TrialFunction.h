/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_TrialFunction_H_
#define model_TrialFunction_H_

#include <utility>
#include <deque>
#include <map>
#include <list>
#include <assert.h>

#include <TerminalColors.h>

#include <TrialBase.h>

#include <TrialFunctionTraits.h>
#include <Constant.h>

#include <TrialProd.h>
#include <TrialGrad.h>
#include <TrialDef.h>

namespace model
{

    template<char _name, int _nComponents, typename _FiniteElementType>
    class TrialFunction : public TrialExpressionBase<TrialFunction<_name,_nComponents,_FiniteElementType>>,
                          public TrialBase<TrialFunction<_name,_nComponents,_FiniteElementType>>
    {

    public:

        typedef double Scalar;
        typedef _FiniteElementType FiniteElementType;
        typedef typename FiniteElementType::ElementType ElementType;
        typedef typename FiniteElementType::NodeType NodeType;

        typedef TrialFunction<_name,_nComponents,FiniteElementType> TrialFunctionType;

        constexpr static int rows=_nComponents;
        constexpr static int cols=1;

        constexpr static int dim=ElementType::dim;
        constexpr static int nodesPerElement=ElementType::nodesPerElement;
        constexpr static int dofPerNode=TypeTraits<TrialFunctionType>::dofPerNode;
        constexpr static int dofPerElement=TypeTraits<TrialFunctionType>::dofPerElement;

        typedef std::map<size_t,double> DirichletConditionContainerType;
        typedef std::map<size_t,std::array<bool,dofPerNode>>  DirichletNodeMapType;

        typedef typename TypeTraits<TrialFunctionType>::BaryType BaryType;
        typedef typename TypeTraits<TrialFunctionType>::ShapeFunctionMatrixType     ShapeFunctionMatrixType;
        typedef typename TypeTraits<TrialFunctionType>::ShapeFunctionGradMatrixType ShapeFunctionGradMatrixType;
        typedef typename TypeTraits<TrialFunctionType>::ShapeFunctionDefMatrixType  ShapeFunctionDefMatrixType;


        explicit TrialFunction(FiniteElementType& fe): TrialBase<TrialFunctionType>(fe) {}
        TrialFunction(const TrialFunctionType&)=delete;
        TrialFunction(TrialFunctionType&&)=delete;
        // Assignment changes values, never the identity captured by expressions.
        TrialFunctionType& operator=(const TrialFunctionType& other)
        {
            return operator=(other.dofVector());
        }

        FiniteElementType& fe() const
        {
            return TrialBase<TrialFunctionType>::fe();
        }

        TrialFunctionType& operator=(const double& c)
        {
            TrialBase<TrialFunctionType>::dofVector().setConstant(c);
            return *this;
        }

        TrialFunctionType& setConstant(const Eigen::Matrix<double,_nComponents,1>& c)
        {
            for(size_t n=0; n< nodeSize();++n)
            {
                TrialBase<TrialFunctionType>::dofVector().template segment<_nComponents>(n*_nComponents)=c;
            }
            return *this;
        }

        TrialFunctionType& operator=(const Eigen::VectorXd& temp)
        {

            assert(size_t(temp.size())==TrialBase<TrialFunctionType>::gSize() && "DOF SIZE MISMATCH");
            TrialBase<TrialFunctionType>::dofVector()=temp;
            return *this;
        }

        TrialFunctionType& operator+=(const Eigen::VectorXd& temp)
        {

            assert(size_t(temp.size())==TrialBase<TrialFunctionType>::gSize() && "DOF SIZE MISMATCH");
            TrialBase<TrialFunctionType>::dofVector()+=temp;
            return *this;
        }

        TrialFunctionType& operator-=(const Eigen::VectorXd& temp)
        {

            assert(size_t(temp.size())==TrialBase<TrialFunctionType>::gSize() && "DOF SIZE MISMATCH");
            TrialBase<TrialFunctionType>::dofVector()-=temp;
            return *this;
        }

        size_t gSize() const
        {
            return TrialBase<TrialFunctionType>::gSize();
        }

        size_t elementSize() const
        {
            return TrialBase<TrialFunctionType>::elementSize();
        }

        size_t nodeSize() const
        {
            return TrialBase<TrialFunctionType>::nodeSize();
        }

        const NodeType& node(const size_t& n) const
        {
            return TrialBase<TrialFunctionType>::node(n);
        }

        static ShapeFunctionMatrixType sfm(const ElementType& ele, const BaryType& bary)
        {
            return ele.template sfm<TrialFunctionType>(bary);
        }

        static ShapeFunctionGradMatrixType sfmGrad(const ElementType& ele, const BaryType& bary)
        {
            return ele.template sfmGrad<TrialFunctionType>(bary);
        }

        static ShapeFunctionDefMatrixType sfmDef(const ElementType& ele, const BaryType& bary)
        {
            return ele.template sfmDef<TrialFunctionType>(bary);
        }

        template<typename Condition>
        void addDirichletCondition(const size_t& nodeListID,
                                   const Condition& cond,
                                   const std::array<bool,dofPerNode>& constrainDof)
        {

            TrialBase<TrialFunctionType>::addDirichletCondition(nodeListID,cond,constrainDof);
        }

        void clearDirichletConditions()
        {
            TrialBase<TrialFunctionType>::clearDirichletConditions();

        }

        DirichletConditionContainerType& dirichletConditions() const
        {
            return TrialBase<TrialFunctionType>::dirichletConditions();
        }

        DirichletNodeMapType& dirichletNodeMap() const
        {
            return TrialBase<TrialFunctionType>::dirichletNodeMap();
        }

        Eigen::VectorXd& dofVector()
        {
            return TrialBase<TrialFunctionType>::dofVector();
        }

        const Eigen::Matrix<double,Eigen::Dynamic,1>& dofVector() const
        {
            return TrialBase<TrialFunctionType>::dofVector();
        }

        Eigen::Matrix<double,dofPerNode,1> dofs(const NodeType& node) const
        {

            return TrialBase<TrialFunctionType>::dofs(node);
        }

        Eigen::Matrix<double,dofPerNode,1> dofs(const size_t& gID) const
        {

            return TrialBase<TrialFunctionType>::dofs(gID);
        }

        Eigen::Matrix<double,dofPerElement,1> dofs(const ElementType& ele) const
        {

            return TrialBase<TrialFunctionType>::dofs(ele);

        }

    };

}
#endif
