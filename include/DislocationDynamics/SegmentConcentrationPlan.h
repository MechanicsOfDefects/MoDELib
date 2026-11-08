#ifndef model_SegmentConcentrationPlan_H_
#define model_SegmentConcentrationPlan_H_
#include <FMM/BoxBoundaryTolerance.h>
#include <FMM/ShiftedBoxNeighbors.h>
#include <DislocationDynamics/ScreenedDislocationAssembly.h>
#include <DislocationDynamics/StraightDislocationSegmentHelmholtz.h>
#include <Eigen/Eigenvalues>
#include <DislocationDynamics/StraightDislocationSegmentFactory.h>
#include <exception>
#include <mutex>
#include <Utilities/ParallelFor.h>
namespace model
{
    // Geometry-only maps: nodal velocities -> leaf coefficients and local concentrations.
    // Endpoint weights include b cross tangent dot climbDirection / species bias, not dL.
    class SegmentConcentrationPlan
    {
        struct Piece
        {
            ClimbSourceSegment segment;
            double begin,end;
        };
        using Matrix=Eigen::SparseMatrix<double,Eigen::RowMajor>;
        std::vector<std::vector<Piece>> bins;

        Eigen::Vector3d lower,width;
        int n,p;
    public:
        std::array<std::vector<double>,3> markers;
        Matrix moments,local;
        SegmentConcentrationPlan(const BoxTree<3>& tree,const ClimbKernelChannel& channel,
            std::size_t nodeCount,const std::array<std::vector<double>,3>& targets,
            const ScreenedDislocationOptions& options,bool buildLocal,bool buildMoments=true)
            :lower(tree.lower()),width(Eigen::Vector3d::Constant(1.0/(1<<tree.depth()))),
             n(1<<tree.depth()),p(options.interpolationOrder)
        {
            if(p<2||p>6)
            {
                throw std::invalid_argument("Segment concentration interpolation order must be 2..6");
            }
            if(!options.boxShortRange || (options.screened&&!options.interpolateScreenedNear))
            {
                throw std::invalid_argument("Segment concentration requires box-local FMM and screened near interpolation");
            }
            const ShiftedBoxNeighbors neighborhood(tree,options.localImageShifts);
            bins.resize(std::size_t(n)*n*n);
            const double tolerance=boxBoundaryTolerance(lower,tree.upper());
            if(!options.segmentSplits || !options.segmentSplits->matches(tree))
            {
                throw std::logic_error("Concentration plan requires a current shared segment split snapshot");
            }
            for(const auto& s:channel.segments)
            {
                if(s.nodes[0]>=nodeCount||s.nodes[1]>=nodeCount||!std::isfinite(s.weights[0])||!std::isfinite(s.weights[1]))
                {
                    throw std::invalid_argument("Invalid climb source segment");
                }
                for(const auto& piece:options.segmentSplits->pieces(s.p0,s.p1))
                {
                    bins[piece.leaf].push_back({s,piece.begin,piece.end});
                }
            }
            SegmentCoreParameters parameters;
            parameters.coreWidth=options.coreWidth;parameters.ewaldWidth=options.screenWidth;
            parameters.order=options.segmentChebyshevOrder;parameters.resolution=options.segmentTableResolution;
            parameters.maxDistance=2*tree.cellDiameter(tree.depth());
            const auto core=makeStraightDislocationSegmentCore(options.coreType,options.screened,parameters);
            const SegmentMomentIntegrator helper(p,tree.inverseBasis());
            std::vector<Eigen::Triplet<double>> entries;
            if(buildMoments)
            {
            const int P=p*p*p;
            std::vector<std::size_t> occupied;
            for(std::size_t b=0;b<bins.size();++b)
            {
                if(!bins[b].empty()) { occupied.push_back(b); }
            }
            for(auto& component:markers) { component.resize(occupied.size()); }
            std::vector<std::vector<Eigen::Triplet<double>>> boxEntries(occupied.size());
            independentParallelFor(occupied.size(),[&](std::size_t i)
            {
                const auto b=occupied[i];
                const Eigen::Vector3d center=tree.cellCenter(tree.depth(),b);
                for(int c=0;c<3;++c) { markers[c][i]=center[c]; }
                auto& output=boxEntries[i];
                output.reserve(bins[b].size()*P*2);
                for(const auto& piece:bins[b])
                {
                    const auto& s=piece.segment;
                    const auto weights=core->concentrationMoments(helper,s.p0,s.p1,center,width,piece.begin,piece.end,s.weights);
                    for(int j=0;j<P;++j)
                    {
                        for(int k=0;k<2;++k) { output.emplace_back(b*P+j,s.nodes[k],weights(j,k)); }
                    }
                }
            },32);
            appendOrderedRows(entries,boxEntries);
            moments.resize(bins.size()*P,nodeCount);
            moments.setFromTriplets(entries.begin(),entries.end());
            moments.makeCompressed();
            }
            local.resize(targets[0].size(),nodeCount);
            if(!buildLocal)
            {
                return;
            }
            const auto kernel=core->concentration(channel.diffusion);
            entries.clear();
            std::vector<std::vector<Eigen::Triplet<double>>> rowEntries(targets[0].size());
            std::exception_ptr error;
            std::mutex errorMutex;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
            for(std::size_t i=0;i<targets[0].size();++i)
            {
                try
                {
                const Eigen::Vector3d x(targets[0][i],targets[1][i],targets[2][i]);
                if(!tree.contains(x,tolerance))
                {
                    throw std::out_of_range("Concentration target outside FMM tree");
                }
                const Eigen::Vector3d fractional=tree.fractional(x);
                Eigen::Vector3i cell;
                for(int c=0;c<3;++c)
                {
                    cell[c]=std::clamp(int(fractional[c]/width[c]),0,n-1);
                }
                const auto leaf=cell[0]+n*(cell[1]+n*cell[2]);
                for(const auto& entry:neighborhood.rows[leaf])
                {
                    const Eigen::Vector3d shiftedTarget=x-options.localImageShifts[entry.image];
                            for(const auto& piece:bins[entry.source])
                            {
                                const Eigen::Vector2d values=kernel->concentrationIntegralsNear(piece.segment.p0,piece.segment.p1,shiftedTarget,piece.begin,piece.end);
                                for(int k=0;k<2;++k)
                                {
                                    rowEntries[i].emplace_back(i,piece.segment.nodes[k],values[k]*piece.segment.weights[k]);
                                }
                            }
                    }
                }
                catch(...)
                {
                    std::lock_guard<std::mutex> lock(errorMutex);
                    if(!error)
                    {
                        error=std::current_exception();
                    }
                }
            }
            if(error)
            {
                std::rethrow_exception(error);
            }
            appendOrderedRows(entries,rowEntries);
            local.setFromTriplets(entries.begin(),entries.end());
            local.makeCompressed();
        }
    };
}
#endif
