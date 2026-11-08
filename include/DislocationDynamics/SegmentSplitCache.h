#ifndef model_SegmentSplitCache_H_
#define model_SegmentSplitCache_H_
#include <FMM/BoxBoundaryTolerance.h>
#include <FMM/BoxTree.h>
#include <array>
#include <map>
#include <mutex>
#include <algorithm>
#include <limits>
#include <exception>
#ifdef _OPENMP
#include <omp.h>
#endif
namespace model
{
    // Only geometry lives here. Coincident oriented segments may share their cuts;
    // each physics retains its own identities, nodal maps and source weights.
    class SegmentSplitSnapshot
    {
    public:
        using Endpoints=std::pair<Eigen::Vector3d,Eigen::Vector3d>;
        using Key=std::array<double,6>;
        struct Piece
        {
            std::size_t leaf;
            double begin,end;
            Eigen::Vector3d p0,p1;
        };
        using Geometry=std::map<Key,Endpoints>;
    private:
        const BoxTree<3>* tree_;
        std::uint64_t revision_;
        Geometry geometry_;
        std::map<Key,std::vector<Piece>> pieces_;
        std::size_t leafCount_;
    public:
        static Key key(const Eigen::Vector3d& a,const Eigen::Vector3d& b)
        {
            if(!a.allFinite()||!b.allFinite())
            {
                throw std::invalid_argument("Nonfinite segment split geometry");
            }
            return {a[0],a[1],a[2],b[0],b[1],b[2]};
        }
        static Geometry geometry(const std::vector<Endpoints>& segments)
        {
            Geometry result;
            for(const auto& segment:segments)
            {
                result.emplace(key(segment.first,segment.second),segment);
            }
            return result;
        }
        SegmentSplitSnapshot(const BoxTree<3>& tree,Geometry geometry)
            :tree_(&tree),revision_(tree.revision()),geometry_(std::move(geometry))
        {
            const int n=1<<tree.depth();
            const Eigen::Vector3d lo=Eigen::Vector3d::Zero(),width=Eigen::Vector3d::Constant(1.0/n);
            leafCount_=std::size_t(n)*n*n;
            const double tolerance=boxBoundaryTolerance(tree.lower(),tree.upper());
            struct Work
            {
                const Endpoints* endpoints;
                std::vector<Piece>* pieces;
            };
            std::vector<Work> work;
            work.reserve(geometry_.size());
            // Create map entries before launching workers: each worker owns only
            // its segment's output vector; no concurrent map insertion occurs.
            for(const auto& entry:geometry_)
            {
                auto inserted=pieces_.emplace(entry.first,std::vector<Piece>{});
                work.push_back({&entry.second,&inserted.first->second});
            }
            std::vector<std::exception_ptr> failures(work.size());
            const auto split=[&](std::size_t i)
            {
                const auto& a=work[i].endpoints->first;
                const auto& b=work[i].endpoints->second;
                auto& pieces=*work[i].pieces;
                if(!tree.contains(a,tolerance)||!tree.contains(b,tolerance))
                { throw std::out_of_range("Segment outside split-cache box tree"); }
                const Eigen::Vector3d qa=tree.fractional(a),dq=tree.inverseBasis()*(b-a);
                const Eigen::Vector3d chord=b-a;
                if(chord.squaredNorm()==0)
                {
                    return;
                }
                std::vector<double> cuts{0,1};
                for(int c=0;c<3;++c)
                {
                    if(dq[c]!=0)
                    {
                        for(int j=1;j<n;++j)
                        {
                            const double u=(j*width[c]-qa[c])/dq[c];
                            if(u>32*std::numeric_limits<double>::epsilon()&&u<1-32*std::numeric_limits<double>::epsilon())
                            {
                                cuts.push_back(u);
                            }
                        }
                    }
                }
                std::sort(cuts.begin(),cuts.end());
                // Crossings of multiple faces at one edge/vertex can differ by
                // roundoff after the affine transform. Do not create zero-scale pieces.
                cuts.erase(std::unique(cuts.begin(),cuts.end(),[](double a,double b)
                { return std::abs(a-b)<=32*std::numeric_limits<double>::epsilon(); }),cuts.end());
                for(std::size_t j=1;j<cuts.size();++j)
                {
                    const double begin=cuts[j-1],end=cuts[j];
                    const Eigen::Vector3d mid=qa+(begin+(end-begin)*.5)*dq;
                    Eigen::Vector3i cell;
                    for(int c=0;c<3;++c)
                    {
                        cell[c]=std::clamp(int(std::floor((mid[c]-lo[c])/width[c])),0,n-1);
                    }
                    const Eigen::Vector3d p0=begin==0?a:Eigen::Vector3d(a+begin*chord);
                    const Eigen::Vector3d p1=end==1?b:Eigen::Vector3d(a+end*chord);
                    if((p1-p0).squaredNorm()>0)
                    {
                        pieces.push_back({std::size_t(cell[0]+n*(cell[1]+n*cell[2])),begin,end,p0,p1});
                    }
                }
            };
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(work.size()>=256 && !omp_in_parallel())
#endif
            for(std::ptrdiff_t i=0;i<static_cast<std::ptrdiff_t>(work.size());++i)
            {
                try
                {
                    split(static_cast<std::size_t>(i));
                }
                catch(...)
                {
                    failures[i]=std::current_exception();
                }
            }
            // Exceptions cannot cross an OpenMP region. Rethrow in deterministic
            // segment order before the cache publishes a replacement snapshot.
            for(const auto& failure:failures)
            {
                if(failure)
                {
                    std::rethrow_exception(failure);
                }
            }
        }
        bool matches(const BoxTree<3>& tree) const
        {
            return tree_==&tree && revision_==tree.revision();
        }
        bool matches(const Geometry& geometry) const
        {
            if(geometry.size()!=geometry_.size())
            {
                return false;
            }
            auto old=geometry_.begin();
            for(const auto& entry:geometry)
            {
                if(entry.first!=(old++)->first)
                {
                    return false;
                }
            }
            return true;
        }
        std::size_t leafCount() const { return leafCount_; }
        const std::vector<Piece>& pieces(const Eigen::Vector3d& a,const Eigen::Vector3d& b) const
        {
            const auto found=pieces_.find(key(a,b));
            if(found==pieces_.end())
            {
                throw std::logic_error("Segment absent from shared split snapshot; refresh geometry");
            }
            return found->second;
        }
    };
    class SegmentSplitCache
    {
        mutable std::mutex mutex_;
        std::shared_ptr<const SegmentSplitSnapshot> snapshot_;
        std::size_t builds_=0;
    public:
        std::shared_ptr<const SegmentSplitSnapshot> get(const BoxTree<3>& tree,
            const std::vector<SegmentSplitSnapshot::Endpoints>& segments)
        {
            auto geometry=SegmentSplitSnapshot::geometry(segments);
            std::lock_guard<std::mutex> lock(mutex_);
            if(!snapshot_||!snapshot_->matches(tree)||!snapshot_->matches(geometry))
            {
                auto next=std::make_shared<SegmentSplitSnapshot>(tree,std::move(geometry));
                snapshot_=std::move(next);
                ++builds_;
            }
            return snapshot_;
        }
        std::size_t buildCount() const
        {
            std::lock_guard<std::mutex> lock(mutex_);
            return builds_;
        }
    };
}
#endif
