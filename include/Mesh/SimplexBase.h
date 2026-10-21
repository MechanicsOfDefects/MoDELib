/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SimplexBase_H_
#define model_SimplexBase_H_

#include <StaticID.h>
//#include <NonCopyable.h>
#include <SimplexTraits.h>


namespace model
{
    
	/**************************************************************************/
	/**************************************************************************/	
	template<short int _dim, short int order>
	struct SimplexBase : public StaticID<SimplexBase<_dim,order> >
    {
     
        static constexpr short int dim=_dim;
        typedef SimplexTraits<dim,order> SimplexTraitsType;
        typedef typename SimplexTraitsType::SimplexIDType SimplexIDType;
        
        SimplicialMesh<dim>* const mesh;
        using KeyType=SimplexIDType;
        const SimplexIDType xID;
        const SimplexIDType& key;
        
		/**********************************************************************/
        SimplexBase(SimplicialMesh<dim>* const m,
                    const SimplexIDType& vIN) :
        /* init */ mesh(m)
        /* init */,xID(SimplexTraitsType::sortID(vIN))
        /* init */,key(xID)
        {/*!
          */

        }
        
        SimplexBase(const SimplexBase&)=delete;
        SimplexBase& operator=(const SimplexBase&)=delete;

        std::array<const Simplex<dim,0>*, SimplexTraits<dim,order>::nVertices> vertices() const
        {
            std::array<const Simplex<dim,0>*, SimplexTraits<dim,order>::nVertices> temp{};
            const auto* simplex=static_cast<const Simplex<dim,order>*>(this);
            if constexpr (order==0)
            {
                temp[0]=simplex;
            }
            else
            {
                // Follow owned children, so retained simplices keep their original geometry
                // even when the mesh is reloaded or destroyed.
                for (const auto& child : simplex->children())
                {
                    for (const auto* vertex : child->vertices())
                    {
                        for (size_t v=0;v<temp.size();++v)
                        {
                            if (xID[v]==vertex->xID[0])
                            {
                                temp[v]=vertex;
                            }
                        }
                    }
                }
            }
            return temp;
        }
        
	};
    
}	// close namespace
#endif
