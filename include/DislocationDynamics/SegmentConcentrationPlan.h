#ifndef model_SegmentConcentrationPlan_H_
#define model_SegmentConcentrationPlan_H_
#include <DislocationDynamics/ScreenedDislocationAssembly.h>
#include <DislocationDynamics/StraightDislocationSegmentHelmholtz.h>
#include <Eigen/Eigenvalues>
#include <DislocationDynamics/StraightDislocationSegmentFactory.h>
#include <exception>
#include <mutex>
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
            :lower(tree.lower()),width((tree.upper()-tree.lower())/(1<<tree.depth())),
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
            bins.resize(std::size_t(n)*n*n);
            const double tolerance=64*std::numeric_limits<double>::epsilon()*std::max({1.0,lower.norm(),tree.upper().norm()});
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
            parameters.maxDistance=2*width.norm();
            const auto core=makeStraightDislocationSegmentCore(options.coreType,options.screened,parameters);
            const SegmentMomentIntegrator helper(p);
            std::vector<Eigen::Triplet<double>> entries;
            if(buildMoments)
            {
            const int P=p*p*p;
            for(std::size_t b=0;b<bins.size();++b)
            {
                if(bins[b].empty())
                {
                    continue;
                }
                const Eigen::Vector3d center=lower+width.cwiseProduct(Eigen::Vector3d(b%n+.5,(b/n)%n+.5,b/(n*n)+.5));
                for(int c=0;c<3;++c)
                {
                    markers[c].push_back(center[c]);
                }
                for(const auto& piece:bins[b])
                {
                    const auto& s=piece.segment;
                    const auto weights=core->concentrationMoments(helper,s.p0,s.p1,center,width,piece.begin,piece.end,s.weights);
                    for(int j=0;j<P;++j)
                    {
                        for(int k=0;k<2;++k)
                        {
                            entries.emplace_back(b*P+j,s.nodes[k],weights(j,k));
                        }
                    }
                }
            }
            moments.resize(bins.size()*P,nodeCount);
            moments.setFromTriplets(entries.begin(),entries.end());
            moments.makeCompressed();
            }
            local.resize(targets[0].size(),nodeCount);
            if(!buildLocal)
            {
                return;
            }
            const auto kernel=core->concentration(channel.diffusion,false,options.screened?SegmentKernelPart::shortRange:SegmentKernelPart::full);
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
                if(!x.allFinite()||(x.array()<lower.array()-tolerance).any()||(x.array()>tree.upper().array()+tolerance).any())
                {
                    throw std::out_of_range("Concentration target outside FMM tree");
                }
                Eigen::Vector3i cell;
                for(int c=0;c<3;++c)
                {
                    cell[c]=std::clamp(int((x[c]-lower[c])/width[c]),0,n-1);
                }
                for(int z=std::max(0,cell[2]-1);z<std::min(n,cell[2]+2);++z)
                {
                    for(int y=std::max(0,cell[1]-1);y<std::min(n,cell[1]+2);++y)
                    {
                        for(int bx=std::max(0,cell[0]-1);bx<std::min(n,cell[0]+2);++bx)
                        {
                            for(const auto& piece:bins[bx+n*(y+n*z)])
                            {
                                const Eigen::Vector2d values=kernel->concentrationIntegrals(piece.segment.p0,piece.segment.p1,x,piece.begin,piece.end);
                                for(int k=0;k<2;++k)
                                {
                                    rowEntries[i].emplace_back(i,piece.segment.nodes[k],values[k]*piece.segment.weights[k]);
                                }
                            }
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
            for(auto& row:rowEntries)
            {
                entries.insert(entries.end(),row.begin(),row.end());
            }
            local.setFromTriplets(entries.begin(),entries.end());
            local.makeCompressed();
        }
    };
}
#endif
