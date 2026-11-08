#ifndef model_ShiftedBoxNeighbors_H_
#define model_ShiftedBoxNeighbors_H_
#include <FMM/BoxTree.h>
#include <set>
#include <functional>
#include <cmath>
#include <stdexcept>
namespace model
{
    // Explicit image neighborhood only, NOT a complete periodic Green operator.
    // Source geometry and coefficients are stored once in the central tree.
    inline void validateImageShifts(const std::vector<Eigen::Vector3d>& shifts)
    {
        std::set<std::array<double,3>> seen;
        if(shifts.empty()) { throw std::invalid_argument("Empty FMM image neighborhood"); }
        for(const auto& s:shifts)
        {
            if(!s.allFinite() || !seen.insert({s[0],s[1],s[2]}).second)
            { throw std::invalid_argument("Nonfinite or duplicate FMM image shift"); }
        }
    }
    // Minimal explicit root neighborhood under the same sphere criterion.
    // All remaining images belong to the deferred exterior periodic operator.
    inline std::vector<Eigen::Vector3d> explicitFMMImageNeighborhood(
        const BoxTree<3>& tree,const Eigen::Matrix<double,3,Eigen::Dynamic>& basis)
    {
        const int rank=basis.cols();
        if(rank==0) { return {Eigen::Vector3d::Zero()}; }
        if(rank>3 || !basis.allFinite()) { throw std::invalid_argument("Invalid periodic FMM basis"); }
        const auto qr=basis.colPivHouseholderQr();
        if(qr.rank()!=rank) { throw std::invalid_argument("Dependent periodic FMM basis"); }
        const Eigen::MatrixXd inverse=qr.solve(Eigen::Matrix3d::Identity());
        const double radius=1.01*tree.cellDiameter(0);
        Eigen::VectorXi limit(rank),index=Eigen::VectorXi::Zero(rank);
        std::size_t candidates=1;
        for(int c=0;c<rank;++c)
        {
            const double bound=std::ceil(radius*inverse.row(c).norm());
            if(!std::isfinite(bound)||bound>1000000 || candidates>1000000/(2*std::size_t(bound)+1))
            { throw std::invalid_argument("Periodic FMM neighborhood too large; reduce the skewed lattice basis"); }
            limit[c]=int(bound);candidates*=2*std::size_t(bound)+1;
        }
        std::vector<Eigen::Vector3d> shifts;
        const std::function<void(int)> visit=[&](int c)
        {
            if(c==rank)
            {
                const Eigen::Vector3d shift=basis*index.cast<double>();
                if(!tree.admissible(shift,0)) { shifts.push_back(shift); }
                return;
            }
            for(index[c]=-limit[c];index[c]<=limit[c];++index[c]) { visit(c+1); }
        };
        visit(0);
        validateImageShifts(shifts);
        return shifts;
    }
    struct ShiftedBoxNeighbors
    {
        struct Entry { std::size_t source,image; };
        std::vector<std::vector<Entry>> rows;
        ShiftedBoxNeighbors(const BoxTree<3>& tree,const std::vector<Eigen::Vector3d>& shifts)
        {
            validateImageShifts(shifts);
            const int n=1<<tree.depth();
            rows.resize(std::size_t(n)*n*n);
            Eigen::Vector3d radius;
            for(int c=0;c<3;++c)
            { radius[c]=1.01*tree.cellDiameter(tree.depth())*tree.inverseBasis().row(c).norm()*n; }
            for(std::size_t b=0;b<rows.size();++b)
            {
                const Eigen::Vector3d center=tree.cellCenter(tree.depth(),b);
                for(std::size_t image=0;image<shifts.size();++image)
                {
                    const Eigen::Vector3d q=n*tree.fractional(center-shifts[image])-Eigen::Vector3d::Constant(.5);
                    Eigen::Vector3i lo,hi;
                    for(int c=0;c<3;++c)
                    {
                        lo[c]=int(std::clamp(std::ceil(q[c]-radius[c]),0.0,double(n)));
                        hi[c]=int(std::clamp(std::floor(q[c]+radius[c]),-1.0,double(n-1)));
                    }
                    for(int z=lo[2];z<=hi[2];++z) { for(int y=lo[1];y<=hi[1];++y) { for(int x=lo[0];x<=hi[0];++x)
                    {
                        const std::size_t source=x+n*(y+n*z);
                        if(!tree.admissible(center-tree.cellCenter(tree.depth(),source)-shifts[image],tree.depth()))
                        { rows[b].push_back({source,image}); }
                    } } }
                }
            }
        }
    };
}
#endif
