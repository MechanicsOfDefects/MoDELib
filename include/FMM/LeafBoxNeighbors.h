/* This file is part of MoDELib, distributed under GNU GPL v2. */
#ifndef model_LeafBoxNeighbors_H_
#define model_LeafBoxNeighbors_H_
#include <FMM/BoxParticleIndex.h>
namespace model
{
    // Compact 27-box stencil on the same hierarchy as the long-field plan.
    // Stores box adjacency, never an all-particle pair list. Free-space bounds;
    // periodic images must be supplied explicitly by the caller.
    struct LeafBoxNeighbors
    {
        BoxParticleIndex<3> index;
        std::vector<std::size_t> offsets,boxes;
        using Positions=std::array<std::vector<double>,3>;
        LeafBoxNeighbors(const BoxTree<3>& tree,const Positions& src,const Positions& tgt)
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
            offsets.push_back(0);
            for(std::size_t b=0;b<leaves.size();++b)
            {
                const auto k=leaves[b]->key();
                const int x=k[0]/step,y=k[1]/step,z=k[2]/step;
                for(int dz=-1;dz<=1;++dz)
                {
                    for(int dy=-1;dy<=1;++dy)
                    {
                        for(int dx=-1;dx<=1;++dx)
                        {
                            const int i=x+dx,j=y+dy,l=z+dz;
                            if(i>=0&&j>=0&&l>=0&&i<n&&j<n&&l<n)
                            {
                                const auto s=lex[i+n*(j+n*l)];
                                if(index.sources().offsets[s+1]>index.sources().offsets[s]) {
                boxes.push_back(s);
              }
            }
                        }
                    }
                }
                offsets.push_back(boxes.size());
            }
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
