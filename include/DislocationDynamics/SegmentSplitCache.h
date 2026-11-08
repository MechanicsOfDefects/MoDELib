#ifndef model_SegmentSplitCache_H_
#define model_SegmentSplitCache_H_
#include <FMM/BoxTree.h>
#include <array>
#include <map>
#include <mutex>
#include <algorithm>
#include <limits>
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
            const Eigen::Vector3d lo=tree.lower(),hi=tree.upper(),width=(hi-lo)/n;
            leafCount_=std::size_t(n)*n*n;
            const double tolerance=64*std::numeric_limits<double>::epsilon()*std::max({1.,lo.norm(),hi.norm()});
            for(const auto& entry:geometry_)
            {
                const auto& a=entry.second.first;
                const auto& b=entry.second.second;
                auto& pieces=pieces_[entry.first];
                if((a.array()<lo.array()-tolerance).any()||(b.array()<lo.array()-tolerance).any()||
                   (a.array()>hi.array()+tolerance).any()||(b.array()>hi.array()+tolerance).any())
                {
                    throw std::out_of_range("Segment outside split-cache box tree");
                }
                const Eigen::Vector3d chord=b-a;
                if(chord.squaredNorm()==0)
                {
                    continue;
                }
                std::vector<double> cuts{0,1};
                for(int c=0;c<3;++c)
                {
                    if(chord[c]!=0)
                    {
                        for(int j=1;j<n;++j)
                        {
                            const double u=(lo[c]+j*width[c]-a[c])/chord[c];
                            if(u>0&&u<1)
                            {
                                cuts.push_back(u);
                            }
                        }
                    }
                }
                std::sort(cuts.begin(),cuts.end());
                cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
                for(std::size_t j=1;j<cuts.size();++j)
                {
                    const double begin=cuts[j-1],end=cuts[j];
                    const Eigen::Vector3d mid=a+(begin+(end-begin)*.5)*chord;
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
