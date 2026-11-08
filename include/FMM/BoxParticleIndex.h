/* This file is part of MoDELib, distributed under GNU GPL v2. */
#ifndef model_BoxParticleIndex_H_
#define model_BoxParticleIndex_H_
#include <FMM/BoxTree.h>
#include <array>
#include <cmath>
#include <stdexcept>
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
                             const std::array<std::vector<double>,dim>& positions)
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
                std::size_t cell=0, stride=1;
                for(int a=0;a<dim;++a)
                {
                    const double x=positions[a][p];
                    if(!std::isfinite(x) || x<tree.lower()[a] || x>tree.upper()[a])
                    {
                        throw std::out_of_range("Particle outside BoxTree root");
                    }
                    const auto i=std::min(n-1,static_cast<unsigned>(std::floor(
                        (x-tree.lower()[a])/(tree.upper()[a]-tree.lower()[a])*n)));
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
            auto nextSources=assign(tree,leafIDs,sources);
            auto nextTargets=assign(tree,leafIDs,targets);
            sources_=std::move(nextSources);
            targets_=std::move(nextTargets);
            tree_=&tree;
            revision_=tree.revision();
        }
    };
}
#endif
