/* This file is part of MoDELib, distributed under GNU GPL v2. */
#ifndef model_LeafBoxNeighbors_H_
#define model_LeafBoxNeighbors_H_
#include <FMM/BoxParticleIndex.h>
#include <FMM/ShiftedBoxNeighbors.h>
namespace model
{
    // Physical near stencil on the same hierarchy as the long-field plan.
    // Cubes have 27 neighbors; skewed/elongated cells may require more.
    // Stores box adjacency, never an all-particle pair list. Free-space bounds;
    // periodic images must be supplied explicitly by the caller.
    struct LeafBoxNeighbors
    {
        BoxParticleIndex<3> index;
        std::vector<std::size_t> offsets,boxes,images;
        using Positions=std::array<std::vector<double>,3>;
        LeafBoxNeighbors(const BoxTree<3>& tree,const Positions& src,const Positions& tgt,const std::vector<Eigen::Vector3d>& shifts={Eigen::Vector3d::Zero()})
        {
            index.rebuild(tree,src,tgt);
            const int n=1<<tree.depth();
            const auto step=BoxTree<3>::latticeSize/n;
            const auto& leaves=tree.boxes(tree.depth());
            std::vector<std::size_t> lex(leaves.size());
            for(std::size_t b=0;b<leaves.size();++b)
            {
                const auto k=leaves[b]->key();
                lex[k[0]/step+n*(k[1]/step+n*(k[2]/step))]=b;
            }
            const ShiftedBoxNeighbors neighborhood(tree,shifts);
            offsets.push_back(0);
            for(std::size_t b=0;b<leaves.size();++b)
            {
                const auto k=leaves[b]->key();
                const int x=k[0]/step,y=k[1]/step,z=k[2]/step;
                for(const auto& entry:neighborhood.rows[x+n*(y+n*z)])
                {
                    const auto source=lex[entry.source];
                    if(index.sources().offsets[source+1]>index.sources().offsets[source])
                    { boxes.push_back(source);images.push_back(entry.image); }
                }
                offsets.push_back(boxes.size());
            }
        }
        template<class F> void forEachImage(std::size_t target,F f) const
        {
            const auto b=index.targets().leafOfPoint[target];
            const auto& s=index.sources();
            for(auto k=offsets[b];k<offsets[b+1];++k)
            { for(auto j=s.offsets[boxes[k]];j<s.offsets[boxes[k]+1];++j) { f(s.permutation[j],images[k]); } }
        }
        template<class F> void forEach(std::size_t target,F f) const
        {
            const auto b=index.targets().leafOfPoint[target];
            const auto& s=index.sources();
            for(auto k=offsets[b];k<offsets[b+1];++k) {
      for(auto j=s.offsets[boxes[k]];j<s.offsets[boxes[k]+1];++j) {
        f(s.permutation[j]);
      }
    }
  }
    };
}
#endif
