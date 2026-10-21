/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SimplicialMesh_H_
#define model_SimplicialMesh_H_

#include <chrono>
#include <sstream>
#include <fstream>
#include <assert.h>
#include <utility>      // std::pair, std::make_pair
#include <set>

#include <map>
//#include <VertexReader.h>
#include <TerminalColors.h>
#include <SimplexTraits.h>
#include <Simplex.h>
#include <SimplexReader.h>
#include <IDreader.h>
//#include <MeshStats.h>
 // defines mode::cout
#include <MeshRegionBoundary.h>
#include <SimplexFactory.h>
//#include <SimplicialMeshFace.h>
#include <GmshReader.h>
#include <PeriodicSimplexPair.h>
#include <PolyhedralMesh.h>


namespace model
{
    template<int dim> class PolyhedralMesh;
    
    
    /**************************************************************************/
    /**************************************************************************/
    template<int _dim>
    class SimplicialMesh : public SimplexReader<_dim>
    /*                  */,public std::map<typename SimplexTraits<_dim,_dim>::SimplexIDType, // key
    /*                                */ std::shared_ptr<Simplex<_dim,_dim>>>
    /*                  */,public std::map<std::pair<size_t,size_t>,MeshRegionBoundary<_dim>> // MeshRegionBoundary container
    //    /*                  */,public std::deque<SimplicialMeshFace<_dim>> // MeshRegionBoundary container
    {
     
        typedef Eigen::Matrix<double,_dim,1> VectorDim;
        
        SimplexFactories<_dim> _simplexFactories{this};
        VectorDim _xMin;
        VectorDim _xMax;
        VectorDim _xC;

        double vol0;
        std::shared_ptr<const int> lifetime=std::make_shared<const int>(0);
        PeriodicSimplexPairMap<_dim> _periodicSimplexPairs;
        Eigen::Matrix<double,_dim,Eigen::Dynamic> _periodicBasis = Eigen::Matrix<double,_dim,Eigen::Dynamic>(_dim,0);
        
        std::shared_ptr<PolyhedralMesh<_dim>> _polyhedralMesh = std::make_shared<PolyhedralMesh<_dim>>();
        void createMesh(const std::set<int>&);
        
        
    public:
        
        
        static constexpr int dim=_dim;
        
        
        typedef std::map<typename SimplexTraits<dim,dim>::SimplexIDType, // key
        /*            */ std::shared_ptr<Simplex<dim,dim>>>  SimplexMapType;
        typedef IDreader<1,dim+2,size_t> ElementReaderType;
        typedef Polyhedron<dim,dim> MeshRegionType;

        typedef MeshRegionBoundary<dim> MeshRegionBoundaryType;
        typedef std::pair<size_t,size_t> MeshRegionIDType;
        typedef std::map<MeshRegionIDType,MeshRegionBoundaryType> MeshRegionBoundaryContainerType;
        typedef Eigen::Matrix<double,dim,Eigen::Dynamic> PeriodicBasisType;
        
        SimplicialMesh();
        
        SimplicialMesh(const std::string& meshFileName,const Eigen::Matrix<double,dim,dim>& A,const Eigen::Matrix<double,dim,1>& x0,const std::set<int>&);
        
        void readMesh(const std::string& meshFileName,const Eigen::Matrix<double,dim,dim>& A,const Eigen::Matrix<double,dim,1>& x0,const std::set<int>&);
        
        const SimplexMapType& simplices() const;
        
        SimplexMapType& simplices();
        
        template<short int order> auto& simplexFactory()
        {
            return _simplexFactories.template factory<order>();
        }
        template<short int order> const auto& simplexFactory() const
        {
            return _simplexFactories.template factory<order>();
        }
        const PolyhedralMesh<_dim>& polyhedralMesh() const;
        const auto& regions() const { return _polyhedralMesh->template factory<_dim>(); }
        Polyhedron<_dim,_dim>* region(const size_t& id) { return regions().at(id).get(); }
        const Polyhedron<_dim,_dim>* region(const size_t& id) const { return regions().at(id).get(); }
        std::shared_ptr<Polyhedron<_dim,_dim>> getSharedRegion(size_t id) { return _polyhedralMesh->getRegion(id); }


        void updateRegions();
        
        void updateRegionBoundaries();
        
        // Discover the full valid lattice first, then select active translations using
        // external face IDs. Only integer-combination pairs are activated; unmatched
        // faces remain physical boundaries. -1 enables all; an empty set disables pairing.
        void identifyParallelFaces(const std::set<int>&);
        
        Eigen::Matrix<double,dim,Eigen::Dynamic> periodicBasis() const;
        std::weak_ptr<const int> lifetimeToken() const
        {
            return lifetime;
        }
        const PeriodicSimplexPairMap<dim>& periodicSimplexPairs() const
        {
            return _periodicSimplexPairs;
        }
        

        void insertSimplex(const typename SimplexTraits<dim,dim>::SimplexIDType& xIN,const int& regionID);
        
        std::pair<bool,const Simplex<_dim,_dim>*> search(const Eigen::Matrix<double,dim,1>& P) const;
        
        std::pair<bool,const Simplex<_dim,_dim>*> searchWithGuess(const Eigen::Matrix<double,dim,1>& P,
                                                                const Simplex<dim,dim>* const guess) const;
        
        std::pair<bool,const Simplex<_dim,_dim>*> searchRegion(const int& regionID,
                                                             const Eigen::Matrix<double,dim,1>& P) const;
        
        std::pair<bool,const Simplex<_dim,_dim>*> searchRegionWithGuess(const Eigen::Matrix<double,dim,1>& P,
                                                                      const Simplex<dim,dim>* const guess) const;
        
        std::pair<bool,const Simplex<_dim,_dim>*> searchWithGuess(const bool& searchAllRegions,
                                                                const Eigen::Matrix<double,dim,1>& P,
                                                                const Simplex<dim,dim>* const guess,
                                                                std::set<const Simplex<dim,dim>*>& searchSet) const;
        void checkSearch(const bool& searchAllRegions,
                         const Eigen::Matrix<double,dim,1>& P,
                         const Simplex<dim,dim>* const guess,
                         const std::pair<bool,const Simplex<dim,dim>*>& lastSearched) const;
        std::pair<bool,const Simplex<_dim,_dim>*> isStrictlyInsideMesh(const Eigen::Matrix<double,dim,1>& P,
                                                                     const Simplex<dim,dim>* const guess,
                                                                     const double& tol) const;
        std::pair<bool,const Simplex<_dim,_dim>*> isOnMeshBoundary(const Eigen::Matrix<double,dim,1>& P, const Simplex<dim,dim>* const guess, const double& tol) const;
        const Eigen::Matrix<double,_dim,1>& xMin() const;
        const double& xMin(const int& k) const;
        const Eigen::Matrix<double,_dim,1>& xMax() const;
        const double& xMax(const int& k) const;
        const Eigen::Matrix<double,_dim,1>& xCenter() const;
        const double& volume() const;
        const MeshRegionBoundaryContainerType& regionBoundaries() const;
        MeshRegionBoundaryContainerType& regionBoundaries();
        const MeshRegionBoundaryType& regionBoundary(const int& i,const int& j) const;
        
    };
    
}
#endif
