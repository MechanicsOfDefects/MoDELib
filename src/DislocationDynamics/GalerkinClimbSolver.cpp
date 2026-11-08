/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GalerkinClimbSolver_cpp_
#define model_GalerkinClimbSolver_cpp_

#include <DislocationDynamics/DislocationFieldBase.h>
#include <limits>
#include <deque>
#include <Utilities/SolverStageTimer.h>


#include <ClusterDynamics/ClusterDynamicsParameters.h>
#include <GalerkinClimbSolver.h>
#include <ConstrainedVelocitySolve.h>
#include <TerminalColors.h>
#include <EqualIteratorRange.h>

#include <DislocationNetwork.h>
#include <DislocationDynamics/DislocationConcentrationSources.h>

namespace model
{

    template <int dim>
    GalerkinClimbSolver<dim>::GalerkinClimbSolver(const DislocationNetwork<dim>& DN_in,const ClusterDynamics<dim>* const CD_in) :
    /* init */ DislocationClimbSolverBase<dim>(DN_in,CD_in)
    {
        std::cout<<greenBoldColor<<"Creating GalerkinClimbSolver"<<defaultColor<<std::endl;
    }

    template <int dim>
    typename GalerkinClimbSolver<dim>::ForceVectorMatrixType GalerkinClimbSolver<dim>::clusterForceKernel(const int& k,const NetworkLinkType& fieldSegment) const
    {
        const Eigen::Matrix<double,1,mSize> deltaC(fieldSegment.quadraturePoint(k).cDD-fieldSegment.quadraturePoint(k).cCD);
        const double u(fieldSegment.quadraturePoint(k).abscissa);
        return  (ForceVectorMatrixType()<<(1.0-u)*deltaC,u*deltaC).finished();
    }

    template <int dim>
    typename GalerkinClimbSolver<dim>::ForceVectorMatrixType GalerkinClimbSolver<dim>::clusterForceVector(const NetworkLinkType& fieldSegment) const
    {
        ForceVectorMatrixType F(ForceVectorMatrixType::Zero());
        const auto bxt(fieldSegment.burgers().cross(fieldSegment.chord()));
        NetworkLinkType::QuadratureDynamicType::integrate(fieldSegment.quadraturePoints().size(),this,F,&GalerkinClimbSolver<dim>::clusterForceKernel,fieldSegment);
        F.row(0)*=bxt.dot(fieldSegment.source->climbDirection());
        F.row(1)*=bxt.dot(fieldSegment.sink->climbDirection());
        return F;
    }

    template <int dim>
    typename GalerkinClimbSolver<dim>::StiffnessMatrixType GalerkinClimbSolver<dim>::clusterStiffnessKernel(const int& k,const NetworkLinkType& fieldSegment,const NetworkLinkType& sourceSegment) const
    {
        StiffnessMatrixType temp(StiffnessMatrixType::Zero());
        if(sourceSegment.grains().size() == 1)
        {
            const auto bxt(fieldSegment.burgers().cross(fieldSegment.chord()));
            const double u(fieldSegment.quadraturePoint(k).abscissa);
            const Eigen::Matrix<double,2,1> m(( Eigen::Matrix<double,2,1>()<<(1.0-u)*bxt.dot(fieldSegment.source->climbDirection()),
                                               /*                               */ u*bxt.dot(fieldSegment.sink->climbDirection())).finished());
            const auto sourceConcentratoinMatrices(sourceSegment.concentrationMatrices(fieldSegment.quadraturePoint(k).r,this->CD->cdp));
            for(int kc=0; kc < mSize; kc++)
            {
                temp.template block<2,2>(kc*2,0)+=m*sourceConcentratoinMatrices.row(kc);
            }
        }
        return temp;
    }

    template <int dim>
    typename GalerkinClimbSolver<dim>::StiffnessMatrixType GalerkinClimbSolver<dim>::clusterStiffnessMatrix(const NetworkLinkType& fieldSegment,const NetworkLinkType& sourceSegment) const
    {
        StiffnessMatrixType K(StiffnessMatrixType::Zero());
        NetworkLinkType::QuadratureDynamicType::integrate(fieldSegment.quadraturePoints().size(),this,K,&GalerkinClimbSolver<dim>::clusterStiffnessKernel,fieldSegment,sourceSegment);
        return K;
    }

    template<int dim>
    ClimbAssemblySnapshot GalerkinClimbSolver<dim>::screenedAssemblySnapshot() const
    {
        if(!this->CD)throw std::logic_error("Screened climb requires ClusterDynamics");
        if(this->DN.ddBase.periodicLatticeBasis.cols()!=0)
            throw std::logic_error("Screened climb periodic correction is not implemented");
        ClimbAssemblySnapshot snapshot;
        snapshot.nodeCount=this->DN.networkNodes().size();snapshot.speciesCount=mSize;
        snapshot.rhs.assign(mSize,Eigen::VectorXd::Zero(snapshot.nodeCount));
        const auto& cdp=this->CD->cdp;
        for(const auto& item:this->DN.networkLinks())
        {
            const auto link=item.second.lock();
            if(!link)throw std::logic_error("Expired link during climb assembly");
            if(link->hasZeroBurgers()||link->isBoundarySegment()||link->isGrainBoundarySegment()||
               link->chordLength()<=std::numeric_limits<float>::epsilon())continue;
            const std::array<std::size_t,2> nodes{link->source->gID(),link->sink->gID()};
            if(nodes[0]>=snapshot.nodeCount||nodes[1]>=snapshot.nodeCount)throw std::logic_error("Noncontiguous climb node indices");
            if(link->isSessile())
            {
                if(link->quadraturePoints().empty())throw std::logic_error("Update segment quadrature fields before assembling climb RHS");
                const auto force=clusterForceVector(*link);
                for(int species=0;species<mSize;++species){snapshot.rhs[species](nodes[0])+=force(0,species);snapshot.rhs[species](nodes[1])+=force(1,species);}
                const auto bxt=link->burgers().cross(link->chord());
                for(const auto& qp:link->quadraturePoints())
                {
                    ClimbQuadratureMap map;map.nodes=nodes;
                    for(int c=0;c<3;++c)map.position[c]=qp.r[c];
                    map.weights={qp.weight*(1-qp.abscissa)*bxt.dot(link->source->climbDirection()),
                                 qp.weight*qp.abscissa*bxt.dot(link->sink->climbDirection())};
                    snapshot.targets.push_back(map);
                }
            }
        }
        snapshot.channels=concentrationSources(this->DN,cdp,true);
        return snapshot;
    }
    template<int dim>
    std::unique_ptr<ScreenedClimbAssembly> GalerkinClimbSolver<dim>::assembleScreened(
        const ScreenedDislocationOptions& options) const
    {
        auto snapshot=screenedAssemblySnapshot();
        const auto context=this->DN.ddBase.fmmContext;
        auto sharedOptions=options;
        sharedOptions.segmentSplits=this->DN.segmentSplits();
        return std::make_unique<ScreenedClimbAssembly>(context->manager,context->boxTree,snapshot,sharedOptions);
    }

    template <int dim>
    Eigen::VectorXd GalerkinClimbSolver<dim>::getNodeVelocitiesPipe() const
    {
        Eigen::VectorXd nodeVelocities(Eigen::VectorXd::Zero(this->DN.networkNodes().size()*dim));
        return nodeVelocities;
    }

    template <int dim>
    void GalerkinClimbSolver<dim>::computeClimbScalarVelocities(const bool& updateStiffness)
    {
        std::cout<<", climbSolver "<<std::flush;
        if(this->DN.ddBase.fmmContext->usesFMM())computeClimbScalarVelocitiesFMM(updateStiffness);
        else computeClimbScalarVelocitiesBulk(updateStiffness);
    }

    template<int dim>
    void GalerkinClimbSolver<dim>::computeClimbScalarVelocitiesFMM(bool updateStiffness)
    {
        const auto context=this->DN.ddBase.fmmContext;
        if(this->DN.ddBase.periodicLatticeBasis.cols()!=0)throw std::runtime_error("Periodic FMM climb is not implemented; set multipoleMode=none");
        const auto count=this->DN.networkNodes().size();
        if(updateStiffness||!screenedAssembly||cachedNodeCount!=count)
        {
            SolverStageTimer stage("FMM climb sparse assembly and long-range plans");
            ScreenedDislocationOptions options;options.segmentSplits=this->DN.segmentSplits();options.coreType=this->DN.multipoleSettings.coreType;options.segmentChebyshevOrder=this->DN.multipoleSettings.ewaldSegmentChebyshevOrder;options.segmentTableResolution=this->DN.multipoleSettings.ewaldSegmentTableResolution;options.screened=context->usesEwaldFMM();options.coreWidth=DislocationFieldBase<dim>::a;options.screenWidth=this->DN.multipoleSettings.ewaldWidth;
            options.shortCutoff=this->DN.multipoleSettings.preconditionerCutoff;options.interpolationOrder=context->interpolationOrder;
            std::shared_ptr<ScreenedClimbAssembly> next=assembleScreened(options);
            const auto Z=periodicClimbVelocityMap<dim>(this->DN);
            std::vector<std::unique_ptr<CoupledClimbSolve>> solvers;
            for(int species=0;species<mSize;++species)
            {
                auto solver=std::make_unique<CoupledClimbSolve>();
                solver->computeOperator([next,species](const Eigen::VectorXd& v){return next->apply(species,v);},next->shortMatrix(species),Z);
                solvers.push_back(std::move(solver));
            }
            screenedAssembly=std::move(next);coupledSolvers=std::move(solvers);cachedNodeCount=count;
        }
        // Concentrations change within the coupled iteration; geometry/maps do not.
        Fc.assign(mSize,Eigen::VectorXd::Zero(count));
        for(const auto& entry:this->DN.networkLinks())
        {
            const auto link=entry.second.lock();
            if(!link||link->hasZeroBurgers()||!link->isSessile()||link->isBoundarySegment()||link->isGrainBoundarySegment()||link->chordLength()<=std::numeric_limits<float>::epsilon())continue;
            const auto force=clusterForceVector(*link);
            for(int s=0;s<mSize;++s){Fc[s](link->source->gID())+=force(0,s);Fc[s](link->sink->gID())+=force(1,s);}
        }
        this->scalarVelocities().assign(count,Eigen::Array<double,1,mSize>::Zero());
        for(int s=0;s<mSize;++s)
        {
            SolverStageTimer stage("FMM climb species "+std::to_string(s)+" GMRES");
            const auto v=coupledSolvers[s]->solve(Fc[s]);
            for(std::size_t i=0;i<count;++i)this->scalarVelocities()[i](s)=v(i);
        }
    }

    template <int dim>
    void GalerkinClimbSolver<dim>::computeClimbScalarVelocitiesBulk(const bool& updateStiffness)
    {
        const size_t nodeCount(this->DN.networkNodes().size());
        const bool rebuild=updateStiffness || coupledSolvers.size()!=mSize || cachedNodeCount!=nodeCount;
        std::cout<<"\n  Climb system: nodes="<<nodeCount<<", species="<<mSize
                 <<", matrices="<<(rebuild?"rebuild":"reuse")<<std::flush;
        SolverStageTimer stage("Climb RHS/dense matrix allocation");
        Fc.assign(mSize,Eigen::VectorXd::Zero(nodeCount));
        std::vector<Eigen::MatrixXd> matrices;
        if(rebuild)
        {
            matrices.assign(mSize,Eigen::MatrixXd::Zero(nodeCount,nodeCount));
        }

        stage.start("Climb field/source preparation and RHS assembly");
        using NetworkLinkPtrType=std::shared_ptr<NetworkLinkType>;
        struct SourceData
        {
            NetworkLinkPtrType link;
            size_t j0;
            size_t j1;
        };
        struct FieldQuadratureData
        {
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW
            VectorDim r;
            double weightedM0;
            double weightedM1;
        };
        using FieldQuadratureContainer=std::vector<FieldQuadratureData,Eigen::aligned_allocator<FieldQuadratureData>>;
        struct FieldData
        {
            NetworkLinkPtrType link;
            size_t i0;
            size_t i1;
            ForceVectorMatrixType force;
            FieldQuadratureContainer quadrature;
        };

        // Lock and filter each network link once, outside the O(N^2) loop.
        std::vector<SourceData> sources;
        std::vector<FieldData> fields;
        sources.reserve(this->DN.networkLinks().size());
        fields.reserve(this->DN.networkLinks().size());
        for(const auto& networkLink : this->DN.networkLinks())
        {
            const auto link(networkLink.second.lock());
            if(!link)
            {
                continue;
            }

            const bool isActive(!link->hasZeroBurgers()
                                && !link->isBoundarySegment()
                                && !link->isGrainBoundarySegment()
                                && link->chordLength()>std::numeric_limits<float>::epsilon());
            // clusterStiffnessKernel returned zero for all other sources.
            if(isActive && link->grains().size()==1)
            {
                sources.push_back({link,link->source->gID(),link->sink->gID()});
            }

            if(isActive && link->isSessile())
            {
                FieldData field{link,link->source->gID(),link->sink->gID(),clusterForceVector(*link),{}};
                field.quadrature.reserve(link->quadraturePoints().size());

                // These quantities depend only on the field segment and its
                // quadrature point, so compute them once rather than per source.
                const auto bxt(link->burgers().cross(link->chord()));
                const double sourceProjection(bxt.dot(link->source->climbDirection()));
                const double sinkProjection(bxt.dot(link->sink->climbDirection()));
                for(const auto& qp : link->quadraturePoints())
                {
                    field.quadrature.push_back({qp.r,
                                                qp.weight*(1.0-qp.abscissa)*sourceProjection,
                                                qp.weight*qp.abscissa*sinkProjection});
                }
                fields.emplace_back(std::move(field));
            }
        }

        // Each field contributes to two shared rows. Atomics protect only the
        // accumulated pair blocks, avoiding a dense matrix copy per thread.
        for(const auto& field : fields)
        {
            for(int kc=0;kc<mSize;++kc)
            {
                Fc[kc](field.i0)+=field.force(0,kc);
                Fc[kc](field.i1)+=field.force(1,kc);
            }
        }
        std::cout<<"\n  Climb interactions: field segments="<<fields.size()<<", source segments="<<sources.size()<<std::flush;
        if(rebuild)
        {
            stage.start("Climb dense interaction assembly (OpenMP)");
#ifdef _OPENMP
#pragma omp parallel for schedule(guided)
#endif
            for(size_t fi=0;fi<fields.size();++fi)
            {
                const auto& field=fields[fi];
                for(const auto& source : sources)
                {
                    StiffnessMatrixType block=StiffnessMatrixType::Zero();
                    for(const auto& qp : field.quadrature)
                    {
                        const auto c=source.link->concentrationMatrices(qp.r,this->CD->cdp);
                        for(int kc=0;kc<mSize;++kc)
                        {
                            block.row(2*kc)+=qp.weightedM0*c.row(kc);
                            block.row(2*kc+1)+=qp.weightedM1*c.row(kc);
                        }
                    }
                    const size_t rows[2]={field.i0,field.i1};
                    const size_t cols[2]={source.j0,source.j1};
                    for(int kc=0;kc<mSize;++kc)
                    {
                        for(int i=0;i<2;++i)
                        {
                            for(int j=0;j<2;++j)
                            {
#ifdef _OPENMP
#pragma omp atomic update
#endif
                                matrices[kc](rows[i],cols[j])+=block(2*kc+i,j);
                            }
                        }
                    }
                }
            }
            stage.start("Climb constraint map and local sparsity pattern");
            const auto Z=periodicClimbVelocityMap<dim>(this->DN);
            // Topologically local preconditioner: self and immediate segment
            // neighbors. Every other interaction remains in the dense operator.
            std::set<std::pair<size_t,size_t>> pattern;
            for(size_t i=0;i<nodeCount;++i)
            {
                pattern.emplace(i,i);
            }
            for(const auto& source : sources)
            {
                pattern.emplace(source.j0,source.j1);
                pattern.emplace(source.j1,source.j0);
            }
            std::vector<std::unique_ptr<CoupledClimbSolve>> next;
            for(int kc=0;kc<mSize;++kc)
            {
                stage.start("Climb species "+std::to_string(kc)+" local matrix construction");
                std::vector<Eigen::Triplet<double>> entries;
                for(const auto& ij : pattern)
                {
                    entries.emplace_back(ij.first,ij.second,matrices[kc](ij.first,ij.second));
                }
                SparseMatrixType local(nodeCount,nodeCount);
                local.setFromTriplets(entries.begin(),entries.end());
                auto solver=std::make_unique<CoupledClimbSolve>();
                stage.finish();
                std::cout<<"\n  Climb species "<<kc<<" operator setup"<<std::flush;
                solver->compute(matrices[kc],local,Z);
                next.emplace_back(std::move(solver));
            }
            coupledSolvers=std::move(next);
            cachedNodeCount=nodeCount;
        }
        stage.finish();
        std::vector<Eigen::Array<double,1,mSize>> nodeV(nodeCount,Eigen::Array<double,1,mSize>::Zero());
        for(int kc=0;kc<mSize;++kc)
        {
            SolverStageTimer speciesTimer("Climb species "+std::to_string(kc)+" linear solve");
            const Eigen::VectorXd constrained=coupledSolvers[kc]->solve(Fc[kc]);
            for(size_t n=0;n<nodeCount;++n)
            {
                nodeV[n](kc)=constrained(n);
            }
        }

        this->scalarVelocities().resize(nodeCount,Eigen::Array<double,1,mSize>::Zero());
        for(size_t k=0;k<nodeCount;++k)
        {
            this->scalarVelocities()[k]=nodeV[k];
        }
    }

    template <int dim>
    Eigen::VectorXd GalerkinClimbSolver<dim>::getNodeVelocitiesBulk() const
    {
        Eigen::VectorXd nodeVelocities(Eigen::VectorXd::Zero(this->DN.networkNodes().size()*dim));
        size_t k=0;
        for(auto& node: this->DN.networkNodes())
        {
            const double vScalarTotal(-1.0*(this->scalarVelocities()[k]*this->CD->cdp.msVector/this->CD->cdp.msVector.abs()).matrix().sum());
            nodeVelocities.template segment<dim>(k*dim)= vScalarTotal*node.second.lock()->climbDirection();
            k++;
        }
        return interpolateClimbSubdivisionVelocities<dim>(this->DN,nodeVelocities);
    }

    template <int dim>
    Eigen::VectorXd GalerkinClimbSolver<dim>::getNodeVelocities() const
    {
        return getNodeVelocitiesBulk()+getNodeVelocitiesPipe();
    }

    template class GalerkinClimbSolver<3>;

}
#endif
