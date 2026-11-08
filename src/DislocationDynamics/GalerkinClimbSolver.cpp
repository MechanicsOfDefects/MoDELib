#include <Utilities/Profiling.h>
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
#include <DislocationDynamics/StraightDislocationSegmentFactory.h>
#include <Utilities/ParallelFor.h>
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
        if(sourceSegment.grains().size()==1 && fieldSegment.grains().size()==1
           && *sourceSegment.grains().begin()==*fieldSegment.grains().begin())
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
        if(!this->CD) { throw std::logic_error("Climb requires ClusterDynamics"); }
        ClimbAssemblySnapshot snapshot;
        snapshot.nodeCount=this->DN.networkNodes().size();snapshot.speciesCount=mSize;
        snapshot.rhs.assign(mSize,Eigen::VectorXd::Zero(snapshot.nodeCount));
        const auto& cdp=this->CD->cdp;
        for(const auto& item:this->DN.networkLinks())
        {
            const auto link=item.second.lock();
            if(!link) { throw std::logic_error("Expired link during climb assembly"); }
            if(link->hasZeroBurgers()||link->isBoundarySegment()||link->isGrainBoundarySegment()||
               link->chordLength()<=std::numeric_limits<float>::epsilon()||link->grains().size()!=1) { continue; }
            const std::array<std::size_t,2> nodes{link->source->gID(),link->sink->gID()};
            if(nodes[0]>=snapshot.nodeCount||nodes[1]>=snapshot.nodeCount) { throw std::logic_error("Noncontiguous climb node indices"); }
            if(link->isSessile())
            {
                if(link->quadraturePoints().empty()) { throw std::logic_error("Update segment quadrature fields before assembling climb RHS"); }
                const auto force=clusterForceVector(*link);
                for(int species=0;species<mSize;++species){snapshot.rhs[species](nodes[0])+=force(0,species);snapshot.rhs[species](nodes[1])+=force(1,species);}
                const auto bxt=link->burgers().cross(link->chord());
                for(const auto& qp:link->quadraturePoints())
                {
                    ClimbQuadratureMap map;map.nodes=nodes;map.grain=(*link->grains().begin())->grainID;
                    for(int c=0;c<3;++c) { map.position[c]=qp.r[c]; }
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
        if(this->DN.ddBase.periodicLatticeBasis.cols()!=0)
        {
            throw std::logic_error("Periodic FMM climb is not implemented");
        }
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
        if(profilingEnabled())
        {
            std::cout<<", climbSolver "<<std::flush;
        }
        if(this->DN.ddBase.fmmContext->usesFMM()) { computeClimbScalarVelocitiesFMM(updateStiffness); }
        else { computeClimbScalarVelocitiesBulk(updateStiffness); }
    }

    template<int dim>
    void GalerkinClimbSolver<dim>::computeClimbScalarVelocitiesFMM(bool rebuild)
    {
        computeGrainClimbVelocities(rebuild,true);
    }

    template<int dim>
    void GalerkinClimbSolver<dim>::computeClimbScalarVelocitiesBulk(const bool& rebuild)
    {
        computeGrainClimbVelocities(rebuild,false);
    }

    template<int dim>
    void GalerkinClimbSolver<dim>::computeGrainClimbVelocities(bool updateStiffness,bool useFMM)
    {
        const auto context=this->DN.ddBase.fmmContext;
        if(useFMM && this->DN.ddBase.periodicLatticeBasis.cols()!=0)
        {
            throw std::runtime_error("Periodic FMM climb is not implemented; set multipoleMode=none");
        }
        const auto count=this->DN.networkNodes().size();
        const bool rebuild=updateStiffness||!systemsReady||cachedNodeCount!=count||cachedFMM!=useFMM;
        if(rebuild)
        {
            SolverStageTimer stage("Climb velocity assembly: independent grain systems",true);
            auto groups=partitionClimbByGrain(screenedAssemblySnapshot());
            const auto maps=grainClimbConstraints(periodicClimbVelocityMap<dim>(this->DN),groups);
            ScreenedDislocationOptions options;
            options.coreType=this->DN.multipoleSettings.coreType;
            options.segmentChebyshevOrder=this->DN.multipoleSettings.ewaldSegmentChebyshevOrder;
            options.segmentTableResolution=this->DN.multipoleSettings.ewaldSegmentTableResolution;
            options.screened=context->usesEwaldFMM();
            options.coreWidth=DislocationFieldBase<dim>::a;
            options.screenWidth=this->DN.multipoleSettings.ewaldWidth;
            options.shortCutoff=this->DN.multipoleSettings.preconditionerCutoff;
            options.interpolationOrder=context->interpolationOrder;
            if(useFMM) { options.segmentSplits=this->DN.segmentSplits(); }
            std::vector<GrainSystem> next;
            for(std::size_t g=0;g<groups.size();++g)
            {
                GrainSystem system;
                system.partition=std::move(groups[g]);
                const auto& snapshot=system.partition.snapshot;
                if(useFMM)
                {
                    system.assembly=std::make_shared<ScreenedClimbAssembly>(context->manager,context->boxTree,snapshot,options);
                    for(int species=0;species<mSize;++species)
                    {
                        auto solver=std::make_unique<CoupledClimbSolve>();
                        const auto assembly=system.assembly;
                        solver->computeOperator([assembly,species](const Eigen::VectorXd& v) { return assembly->apply(species,v); },assembly->shortMatrix(species),maps[g]);
                        system.solvers.push_back(std::move(solver));
                    }
                }
                else
                {
                    SegmentCoreParameters parameters;
                    parameters.coreWidth=options.coreWidth;
                    parameters.order=options.segmentChebyshevOrder;
                    parameters.resolution=options.segmentTableResolution;
                    const auto core=makeStraightDislocationSegmentCore(options.coreType,false,parameters);
                    for(int species=0;species<mSize;++species)
                    {
                        SolverStageTimer timer("Climb grain "+std::to_string(system.partition.grain)+": direct matrix assembly");
                        Eigen::MatrixXd matrix=Eigen::MatrixXd::Zero(snapshot.nodeCount,snapshot.nodeCount);
                        std::vector<std::vector<std::pair<std::size_t,double>>> rows(snapshot.nodeCount);
                        for(std::size_t q=0;q<snapshot.targets.size();++q)
                        {
                            for(int k=0;k<2;++k) { rows[snapshot.targets[q].nodes[k]].emplace_back(q,snapshot.targets[q].weights[k]); }
                        }
                        std::set<std::pair<std::size_t,std::size_t>> pattern;
                        for(std::size_t i=0;i<snapshot.nodeCount;++i) { pattern.emplace(i,i); }
                        for(const auto& channel:snapshot.channels)
                        {
                            if(channel.species!=std::size_t(species)) { continue; }
                            const auto kernel=core->concentration(channel.diffusion);
                            independentParallelFor(snapshot.nodeCount,[&](std::size_t row)
                            {
                                for(const auto& term:rows[row])
                                {
                                    if(term.second==0.0) { continue; }
                                    const auto& target=snapshot.targets[term.first];
                                    const Eigen::Vector3d x(target.position.data());
                                    for(const auto& source:channel.segments)
                                    {
                                        Eigen::Vector2d value=Eigen::Vector2d::Zero();
                                        for(const auto& shift:this->DN.ddBase.periodicShifts)
                                        {
                                            value+=kernel->concentrationIntegrals(source.p0,source.p1,x+shift);
                                        }
                                        for(int k=0;k<2;++k) { matrix(row,source.nodes[k])+=term.second*value[k]*source.weights[k]; }
                                    }
                                }
                            });
                            for(const auto& source:channel.segments)
                            {
                                pattern.emplace(source.nodes[0],source.nodes[1]);
                                pattern.emplace(source.nodes[1],source.nodes[0]);
                            }
                        }
                        std::vector<Eigen::Triplet<double>> entries;
                        for(const auto& ij:pattern) { entries.emplace_back(ij.first,ij.second,matrix(ij.first,ij.second)); }
                        SparseMatrixType local(snapshot.nodeCount,snapshot.nodeCount);
                        local.setFromTriplets(entries.begin(),entries.end());
                        auto solver=std::make_unique<CoupledClimbSolve>();
                        solver->compute(matrix,local,maps[g]);
                        system.solvers.push_back(std::move(solver));
                    }
                }
                next.push_back(std::move(system));
            }
            grainSystems=std::move(next);
            cachedNodeCount=count;cachedFMM=useFMM;systemsReady=true;
        }
        SolverStageTimer rhsStage("Climb velocity assembly: RHS",true);
        Fc.assign(mSize,Eigen::VectorXd::Zero(count));
        for(const auto& entry:this->DN.networkLinks())
        {
            const auto link=entry.second.lock();
            if(!link||link->hasZeroBurgers()||!link->isSessile()||link->isBoundarySegment()||link->isGrainBoundarySegment()||link->grains().size()!=1||link->chordLength()<=std::numeric_limits<float>::epsilon()) { continue; }
            const auto force=clusterForceVector(*link);
            for(int s=0;s<mSize;++s)
            {
                Fc[s](link->source->gID())+=force(0,s);
                Fc[s](link->sink->gID())+=force(1,s);
            }
        }
        rhsStage.finish();
        this->scalarVelocities().assign(count,Eigen::Array<double,1,mSize>::Zero());
        // Each operator parallelizes internally. Keep grain solves sequential to
        // avoid nested Kokkos/OpenMP execution; their states and histories are independent.
        for(auto& system:grainSystems)
        {
            for(int s=0;s<mSize;++s)
            {
                SolverStageTimer stage("Climb grain "+std::to_string(system.partition.grain)+", species "+std::to_string(s)+" GMRES",true);
                const Eigen::VectorXd velocity=system.solvers[s]->solve(system.partition.gather(Fc[s]));
                for(std::size_t i=0;i<system.partition.nodes.size();++i)
                {
                    if(system.partition.active[i]) { this->scalarVelocities()[system.partition.nodes[i]](s)=velocity[i]; }
                }
            }
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
