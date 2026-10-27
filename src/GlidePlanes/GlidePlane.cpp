/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GLIDEPLANE_cpp
#define model_GLIDEPLANE_cpp

#include <GlidePlane.h>

namespace model
{

    template <int dim>
    GlidePlane<dim>::GlidePlane(const GlidePlaneFactoryType* const gpF,
                                const GlidePlaneKeyType& key_in) :
    /* init */ LatticePlane(key_in.planeIndex(),ReciprocalLatticeDirection<dim>(key_in.reciprocalDirectionComponents(),*gpF->poly.grain(key_in.latticeID()))) // BETTER TO CONSTRUCT N WITH PRIMITIVE VECTORS ON THE PLANE
    /* init */,MeshPlane<dim>(gpF->poly.mesh,key_in.latticeID(),this->planeOrigin(),this->n.cartesian())
    /* init */,glidePlaneFactory(*gpF)
    /* init */,grain(*gpF->poly.grain(key_in.latticeID()))
    /* init */,key(key_in)
    {
        VerboseGlidePlane(1,"Creating GlidePlane "<<this->sID<<std::endl;);
    }

    template <int dim>
    GlidePlane<dim>::~GlidePlane()
    {
        VerboseGlidePlane(1,"Destroying GlidePlane "<<this->sID<<std::endl;);
    }

    template <int dim>
    std::set<std::shared_ptr<SlipSystem>> GlidePlane<dim>::slipSystems() const
    {
        std::set<std::shared_ptr<SlipSystem>> temp;
        for(const auto& ss : grain.slipSystems())
        {
            if(this->n.cross(ss.second->n).base().squaredNorm()==0)
            {
                temp.emplace(ss.second);
            }
        }
        return temp;
    }

template <int dim>
std::shared_ptr<GlidePlaneBase> GlidePlane<dim>::planeBase() const
{
    for(const auto& ss : grain.planeNormals())
    {
        if(this->n.cross(*ss.second).base().squaredNorm()==0)
        {
            return ss.second;
        }
    }
    return std::shared_ptr<GlidePlaneBase>(nullptr);
}

template <int dim>
int GlidePlane<dim>::planeBaseID() const
{
    for(const auto& ss : grain.planeNormals())
    {
        if(this->n.cross(*ss.second).base().squaredNorm()==0)
        {
            return int(ss.first);
        }
    }
    return -1;
}


template<int dim>
typename GlidePlane<dim>::BoxCutPtr GlidePlane<dim>::boxIntersection(
    const std::shared_ptr<const BoxTree<dim>>& tree, unsigned level) const
{
    if(!tree)
    {
        throw std::invalid_argument("GlidePlane box intersection requires a tree");
    }
    const std::lock_guard<std::mutex> lock(boxCutMutex);
    for(auto it = boxCuts.begin(); it != boxCuts.end();)
    {
        if(it->first.expired())
        {
            it = boxCuts.erase(it);
        }
        else
        {
            ++it;
        }
    }
    auto& levels = boxCuts[TreeOwner(tree)];
    if(!levels.empty() && levels.begin()->second->revision() != tree->revision())
    {
        levels.clear();
    }
    const auto found = levels.find(level);
    if(found != levels.end())
    {
        return found->second;
    }
    const auto cut = std::make_shared<const PlaneBoxIntersection<dim>>(*tree, level, this->P, this->unitNormal);
    levels.emplace(level, cut);
    return cut;
}

template struct GlidePlane<3>;
}
#endif

