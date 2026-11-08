/* This file is part of MoDELib, distributed under GNU GPL v2. */
#ifndef model_ShortRangeParticleNeighbors_H_
#define model_ShortRangeParticleNeighbors_H_
#include <Eigen/Dense>
#include <array>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace model
{
    // Immutable geometric snapshot. Rebuild after ANY source/target motion or reorder.
    // Positions remain in their original arrays; periodic images carry a translation.
    // Cartesian hash cells of side cutoff are independent of the FMM depth.
    class ShortRangeParticleNeighbors
    {
    public:
        using Positions=std::array<std::vector<double>,3>;
        using Basis=Eigen::Matrix<double,3,Eigen::Dynamic>;
        struct Pair { std::size_t source,image; };
        std::vector<std::size_t> offsets;
        std::vector<Pair> pairs;
        std::vector<Eigen::Vector3d> shifts;
        std::size_t sourceCount=0,targetCount=0,candidates=0;
        double cutoff=0;
    private:
        using Cell=std::array<long long,3>;
        struct Hash
        {
            std::size_t operator()(const Cell& c) const
            {
                std::size_t h=0;
                for(auto v:c) { h^=std::hash<long long>{}(v)+0x9e3779b9+(h<<6)+(h>>2); }
                return h;
            }
        };
        static std::pair<Eigen::Vector3d,Eigen::Vector3d> bounds(const Positions& p)
        {
            Eigen::Vector3d lo,hi;
            for(int j=0;j<3;++j)
            {
                if(p[j].size()!=p[0].size()) { throw std::invalid_argument("Inconsistent neighbor positions"); }
                lo[j]=std::numeric_limits<double>::infinity(); hi[j]=-lo[j];
                for(double x:p[j])
                {
                    if(!std::isfinite(x)) { throw std::invalid_argument("Nonfinite neighbor position"); }
                    lo[j]=std::min(lo[j],x); hi[j]=std::max(hi[j],x);
                }
            }
            return {lo,hi};
        }
    public:
        ShortRangeParticleNeighbors(const Positions& sources,const Positions& targets,double radius,
                                    const Basis& periodicBasis=Basis(3,0))
            :sourceCount(sources[0].size()),targetCount(targets[0].size()),cutoff(radius)
        {
            if(!std::isfinite(radius)||radius<=0||!std::isfinite(radius*radius)||radius*radius==0)
            { throw std::invalid_argument("Invalid neighbor cutoff"); }
            if(periodicBasis.cols()>3||!periodicBasis.allFinite()) { throw std::invalid_argument("Invalid periodic basis"); }
            const auto sb=bounds(sources),tb=bounds(targets);
            const int rank=periodicBasis.cols();
            Eigen::MatrixXd reciprocal(rank,3);
            if(rank)
            {
                Eigen::ColPivHouseholderQR<Basis> qr(periodicBasis);
                if(qr.rank()!=rank) { throw std::invalid_argument("Dependent periodic generators"); }
                reciprocal=qr.solve(Eigen::Matrix3d::Identity());
            }
            offsets.assign(targetCount+1,0);
            if(!sourceCount||!targetCount) { shifts.push_back(Eigen::Vector3d::Zero()); return; }
            // A pair inside the cutoff satisfies B^+(target-source)=n+B^+r.
            // Project the complete displacement bounding box; this covers skew cells,
            // partial periodicity, and cutoffs spanning several replicas.
            std::array<int,3> low{0,0,0},high{0,0,0};
            std::size_t imageCount=1;
            for(int k=0;k<rank;++k)
            {
                double lower=0,upper=0;
                for(int j=0;j<3;++j)
                {
                    const double a=reciprocal(k,j)*(tb.first[j]-sb.second[j]);
                    const double b=reciprocal(k,j)*(tb.second[j]-sb.first[j]);
                    lower+=std::min(a,b); upper+=std::max(a,b);
                }
                const double pad=radius*reciprocal.row(k).norm();
                lower=std::floor(lower-pad)-1; upper=std::ceil(upper+pad)+1;
                if(!std::isfinite(lower)||!std::isfinite(upper)||lower<-1000000||upper>1000000)
                { throw std::length_error("Periodic neighbor image range too large"); }
                low[k]=int(lower); high[k]=int(upper);
                const auto width=std::size_t(high[k]-low[k]+1);
                if(width>1000000/imageCount) { throw std::length_error("Too many periodic neighbor images"); }
                imageCount*=width;
            }
            for(int i=low[0];i<=high[0];++i) { for(int j=low[1];j<=high[1];++j) { for(int k=low[2];k<=high[2];++k)
            {
                Eigen::Vector3d shift=Eigen::Vector3d::Zero(); const int n[3]={i,j,k};
                for(int a=0;a<rank;++a) { shift+=periodicBasis.col(a)*n[a]; }
                // Discard translated source bounding boxes that cannot reach any target.
                double d2=0;
                for(int a=0;a<3;++a)
                {
                    const double gap=std::max({0.0,sb.first[a]+shift[a]-tb.second[a],tb.first[a]-sb.second[a]-shift[a]});
                    d2+=gap*gap;
                }
                if(d2<=radius*radius*(1+32*std::numeric_limits<double>::epsilon())) { shifts.push_back(shift); }
            } } }
            const auto cell=[&](const Eigen::Vector3d& x)
            {
                Cell c;
                for(int j=0;j<3;++j)
                {
                    const double value=std::floor((x[j]-sb.first[j])/radius);
                    if(!std::isfinite(value)||std::abs(value)>1e15) { throw std::out_of_range("Neighbor cell index too large"); }
                    c[j]=static_cast<long long>(value);
                }
                return c;
            };
            std::unordered_map<Cell,std::vector<std::size_t>,Hash> grid;
            grid.reserve(sourceCount);
            for(std::size_t s=0;s<sourceCount;++s)
            { grid[cell(Eigen::Vector3d(sources[0][s],sources[1][s],sources[2][s]))].push_back(s); }
            for(std::size_t t=0;t<targetCount;++t)
            {
                const Eigen::Vector3d x(targets[0][t],targets[1][t],targets[2][t]);
                for(std::size_t im=0;im<shifts.size();++im)
                {
                    const Eigen::Vector3d y=x-shifts[im];
                    if((y.array()<sb.first.array()-radius).any()||(y.array()>sb.second.array()+radius).any()) { continue; }
                    const auto c=cell(y);
                    for(int i=-1;i<=1;++i) { for(int j=-1;j<=1;++j) { for(int k=-1;k<=1;++k)
                    {
                        const auto found=grid.find(Cell{c[0]+i,c[1]+j,c[2]+k});
                        if(found==grid.end()) { continue; }
                        for(const auto s:found->second)
                        {
                            ++candidates;
                            double d2=0;
                            for(int a=0;a<3;++a) { const double d=(x[a]-sources[a][s])-shifts[im][a]; d2+=d*d; }
                            if(d2<=radius*radius) { pairs.push_back({s,im}); }
                        }
                    } } }
                }
                // Stable image/source order matches a brute-force image-first sum.
                std::sort(pairs.begin()+offsets[t],pairs.end(),[](const Pair& a,const Pair& b)
                { return a.image<b.image || (a.image==b.image && a.source<b.source); });
                offsets[t+1]=pairs.size();
            }
        }
    };
}
#endif
