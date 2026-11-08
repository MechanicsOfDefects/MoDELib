/* This file is part of MoDELib, distributed under GNU GPL v2. */
#ifndef model_BoxParticleIndex_H_
#define model_BoxParticleIndex_H_
#include <FMM/BoxBoundaryTolerance.h>
#include <FMM/BoxTree.h>
#include <array>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <limits>
#include <sstream>
#include <iomanip>
namespace model
{
    // Stable counting-sort index; particle arrays and their ownership never move.
    // Leaf IDs index BoxTree::boxes(depth). No periodic wrapping is implicit:
    // points must be expressed in the root cell before assignment.
    template<int dim>
    class BoxParticleIndex
    {
    public:
        struct Points
        {
            std::vector<std::size_t> offsets, permutation, leafOfPoint;
        };
    private:
        const BoxTree<dim>* tree_ = nullptr;
        std::uint64_t revision_ = 0;
        Points sources_, targets_;
        static Points assign(const BoxTree<dim>& tree,
                             const std::vector<std::size_t>& leafIDs,
                             const std::array<std::vector<double>,dim>& positions,const char* role)
        {
            Points result;
            const auto count=positions[0].size();
            for(const auto& component : positions)
            {
                if(component.size()!=count)
                {
                    throw std::invalid_argument("Inconsistent particle position components");
                }
            }
            result.offsets.assign(leafIDs.size()+1,0);
            result.leafOfPoint.resize(count);
            result.permutation.resize(count);
            const unsigned n=1u<<tree.depth();
            for(std::size_t p=0;p<count;++p)
            {
                typename BoxTree<dim>::Vector physicalPoint;
                for(int c=0;c<dim;++c) { physicalPoint[c]=positions[c][p]; }
                const auto fractional=tree.fractional(physicalPoint);
                std::size_t cell=0, stride=1;
                for(int a=0;a<dim;++a)
                {
                    const double x=fractional[a];
                    const double lo=0,hi=1;
                    // Interpolating points on a boundary can round a few ULPs
                    // outside it. Clamp only the indexing coordinate, never the
                    // physical source/target geometry. Reject real excursions.
                    const double tolerance=boxBoundaryTolerance(tree.lower(),tree.upper())*tree.inverseBasis().row(a).norm();
                    if(!std::isfinite(x) || x<lo-tolerance || x>hi+tolerance)
                    {
                        std::ostringstream message;
                        message<<std::setprecision(17)<<"Particle outside BoxTree root: "<<role
                               <<" index="<<p<<", axis="<<a<<", coordinate="<<x
                               <<", bounds=["<<lo<<", "<<hi<<"], tolerance="<<tolerance;
                        throw std::out_of_range(message.str());
                    }
                    const double coordinate=std::clamp(x,lo,hi);
                    const auto i=std::min(n-1,static_cast<unsigned>(std::floor(
                        (coordinate-lo)/(hi-lo)*n)));
                    cell+=stride*i;
                    stride*=n;
                }
                const auto leaf=leafIDs[cell];
                result.leafOfPoint[p]=leaf;
                ++result.offsets[leaf+1];
            }
            for(std::size_t i=1;i<result.offsets.size();++i)
            {
                result.offsets[i]+=result.offsets[i-1];
            }
            auto cursor=result.offsets;
            for(std::size_t p=0;p<count;++p)
            {
                result.permutation[cursor[result.leafOfPoint[p]]++]=p;
            }
            return result;
        }
    public:
        const Points& sources() const { return sources_; }
        const Points& targets() const { return targets_; }
        bool matches(const BoxTree<dim>& tree) const
        {
            return tree_==&tree && revision_==tree.revision();
        }
        // Transactional: malformed/outside targets cannot leave sources updated.
        void rebuild(const BoxTree<dim>& tree,
                     const std::array<std::vector<double>,dim>& sources,
                     const std::array<std::vector<double>,dim>& targets)
        {
            const auto& leaves=tree.boxes(tree.depth());
            std::vector<std::size_t> leafIDs(leaves.size());
            const auto n=1u<<tree.depth();
            const auto width=BoxTree<dim>::latticeSize/n;
            for(std::size_t id=0;id<leaves.size();++id)
            {
                std::size_t cell=0,stride=1;
                for(int a=0;a<dim;++a)
                {
                    cell+=stride*(leaves[id]->key()[a]/width);
                    stride*=n;
                }
                leafIDs[cell]=id;
            }
            auto nextSources=assign(tree,leafIDs,sources,"source");
            auto nextTargets=assign(tree,leafIDs,targets,"target");
            sources_=std::move(nextSources);
            targets_=std::move(nextTargets);
            tree_=&tree;
            revision_=tree.revision();
        }
    };
}
#endif
