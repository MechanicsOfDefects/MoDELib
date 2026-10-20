/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_TrialBase_H_
#define model_TrialBase_H_

#include <iostream>

#include <type_traits>
#include <map>
#include <cfloat>
#include <Eigen/Dense>
#include <TypeTraits.h>
#include <TerminalColors.h>

namespace model
{

    template<typename _TrialFunctionType>
    class TrialBase
    {

    public:
        typedef double Scalar;
        typedef typename TypeTraits<_TrialFunctionType>::FiniteElementType FiniteElementType;
        typedef typename FiniteElementType::MeshType MeshType;
        typedef typename TypeTraits<_TrialFunctionType>::ElementType ElementType;
        typedef typename TypeTraits<_TrialFunctionType>::BaryType BaryType;

        typedef typename TypeTraits<_TrialFunctionType>::NodeType NodeType;
        typedef _TrialFunctionType TrialFunctionType;
        static constexpr int dim=ElementType::dim;
        static constexpr int nodesPerElement=ElementType::nodesPerElement;
        static constexpr int dofPerNode=TypeTraits<TrialFunctionType>::dofPerNode;
        static constexpr int dofPerElement=TypeTraits<TrialFunctionType>::dofPerElement;

        typedef std::map<size_t,std::array<bool,TypeTraits<TrialFunctionType>::dofPerNode>> DirichletNodeMapType;
        typedef std::map<size_t,double> DirichletConditionContainerType;

    private:
        struct State
        {
            std::shared_ptr<FiniteElementType> owner;
            FiniteElementType* space;
            std::weak_ptr<const int> lifetime;
            Eigen::VectorXd dofvector;
            Eigen::VectorXd periodicNodalOffset;
            DirichletNodeMapType nodes;
            DirichletConditionContainerType conditions;
            Eigen::Matrix<double,dofPerNode,dim> periodicGradient=Eigen::Matrix<double,dofPerNode,dim>::Zero();
            explicit State(FiniteElementType& fe):
                owner(fe.weak_from_this().lock()), space(&fe), lifetime(fe.lifetimeToken()),
                dofvector(Eigen::VectorXd::Zero(fe.nodeSize()*dofPerNode)) {}
        };
        std::shared_ptr<State> state;
    public:
        explicit TrialBase(FiniteElementType& fe): state(std::make_shared<State>(fe)) {}
        // Each field owns independent state; only shared pointers alias the field.
        TrialBase(const TrialBase&)=delete;
        TrialBase(TrialBase&&)=delete;
        const TrialBase& field() const
        {
            return *this;
        }
        std::shared_ptr<const void> ownership() const
        {
            return state;
        }
        bool sameField(const TrialBase& other) const
        {
            return state==other.state;
        }
        FiniteElementType& fe() const
        {
            if(state->lifetime.expired())
            {
                throw std::runtime_error("FEM field outlived its borrowed space. Use a shared FEM space for retained expressions.");
            }
            state->space->validateMeshLifetime();
            return *state->space;
        }
        size_t gSize() const
        {
            return nodeSize()*dofPerNode;
        }

        Eigen::Matrix<double,Eigen::Dynamic,1>& dofVector() const
        {
            return state->dofvector;
        }

        // Affine jump across a periodic translation p: x(X+p)-x(X)=gradient*p.
        void setPeriodicGradient(const Eigen::Matrix<double,dofPerNode,dim>& gradient) const
        {
            if(!gradient.allFinite())
            {
                throw std::invalid_argument("Nonfinite FEM periodic gradient.");
            }
            state->periodicGradient=gradient;
        }

        const Eigen::Matrix<double,dofPerNode,dim>& periodicGradient() const
        {
            return state->periodicGradient;
        }

        // Optional known field whose differences supplement affine periodic jumps.
        void setPeriodicNodalOffset(const Eigen::VectorXd& values) const
        {
            if((values.size()!=0 && values.size()!=Eigen::Index(gSize())) || !values.allFinite())
            {
                throw std::invalid_argument("Invalid FEM periodic nodal offset.");
            }
            state->periodicNodalOffset=values;
        }
        const Eigen::VectorXd& periodicNodalOffset() const
        {
            return state->periodicNodalOffset;
        }

        void setFullDofVector(const Eigen::VectorXd& values) const
        {
            if(values.size()!=Eigen::Index(gSize()) || !values.allFinite())
            {
                throw std::invalid_argument("Invalid full FEM DOF vector.");
            }
            state->dofvector=values;
        }

        Eigen::Matrix<double,dofPerNode,1> dofs(const NodeType& node) const
        {
            return state->dofvector.template segment<dofPerNode>(node.gID*dofPerNode);
        }

        Eigen::Matrix<double,dofPerNode,1> dofs(const size_t& gID) const
        {
            return state->dofvector.template segment<dofPerNode>(gID*dofPerNode);
        }

        Eigen::Matrix<double,dofPerElement,1> dofs(const ElementType& ele) const
        {
            Eigen::Matrix<double,dofPerElement,1> temp;
            for(int n=0;n<nodesPerElement;++n)
            {
                temp.template segment<dofPerNode>(n*dofPerNode)=dofs(ele.node(n));
            }
            return temp;
        }

        const ElementType& element(const size_t& n) const
        {
            return fe().element(n);
        }

        const MeshType& mesh() const
        {
            return fe().mesh;
        }

        size_t elementSize() const
        {
            return fe().elementSize();
        }

        size_t nodeSize() const
        {
            return fe().nodeSize();
        }

        const NodeType& node(const size_t& n) const
        {
            return fe().node(n);
        }

        void clearDirichletConditions() const
        {
            state->nodes.clear();
            state->conditions.clear();
        }

        DirichletNodeMapType& dirichletNodeMap() const
        {
            return state->nodes;
        }

        DirichletConditionContainerType& dirichletConditions() const
        {
            return state->conditions;
        }

        template<typename Condition>
        void addDirichletCondition(const size_t& nodeListID,
                                          const Condition& cond,
                                          const std::array<bool,dofPerNode>& constrainDof) const
        {

            bool isConstrained(false);
            for(int dof=0;dof<dofPerNode;++dof)
            {
                isConstrained+=constrainDof[dof];
            }

            if(isConstrained)
            {

                for(const auto& pNode :  TrialBase<TrialFunctionType>::fe().nodeList(nodeListID))
                {

                    Eigen::Matrix<double,dofPerNode,1> value(Eigen::Matrix<double,dofPerNode,1>::Zero());
                    cond(*pNode,value);

                    for(int dof=0;dof<dofPerNode;++dof)
                    {
                        state->nodes[pNode->gID][dof]+=constrainDof[dof];

                        if(constrainDof[dof])
                        {
                            const auto temp=state->conditions.emplace(dofPerNode*(pNode->gID)+dof,value(dof));

                            if(!(temp.second || std::abs(temp.first->second-value(dof))<FLT_EPSILON ))
                            {
                                std::cout<<"FEM node "<<pNode->gID<<", dof "<<dof<<std::endl;
                                std::cout<<"exising dirichlet-condition= "<<temp.first->second<<std::endl;
                                std::cout<<"    new dirichlet-condition= "<<value(dof)<<std::endl;
                                throw std::runtime_error("Conflicting FEM Dirichlet conditions.");
                            }
                        }
                    }
                }
            }

        }

    };

}
#endif
