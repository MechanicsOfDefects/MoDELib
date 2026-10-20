/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SimplicialMesh_cpp_
#define model_SimplicialMesh_cpp_

#include <chrono>
#include <sstream>
#include <fstream>
#include <assert.h>
#include <utility>      // std::pair, std::make_pair
#include <set>
#include <map>

#include <Eigen/LU>

//#include <VertexReader.h>
#include <TerminalColors.h>
#include <SimplexTraits.h>
#include <Simplex.h>
#include <SimplexReader.h>
//#include <MeshStats.h>
#include <MeshRegionObserver.h>
// defines mode::cout
#include <MeshRegionBoundary.h>
#include <SimplexObserver.h>
//#include <SimplicialMeshFace.h>
#include <GmshReader.h>
#include <SimplicialMesh.h>
#include <PeriodicMeshFacePairing.h>


namespace model
{

    template<int dim>
    void SimplicialMesh<dim>::createMesh(const std::set<int>& periodicFaceIDs)
    {/*!
      */
        
        const auto t0= std::chrono::system_clock::now();
        vol0=0.0;
        size_t eleConter=0;
        for (const auto& eIter : this->simplexReader().elements())
        {
            insertSimplex(eIter.second.first,eIter.second.second);
            eleConter++;
            std::cout<<greenBoldColor<<"\r"<<"Creating mesh "<<eleConter*100/this->simplexReader().elements().size()<<"%"<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::flush;
        }
        std::cout<<defaultColor<<std::endl;
        
        this->info(); // print mesh info
        
        if(simplices().size())
        {
            _xMin=this->template observer<0>().begin()->second->P0;
            _xMax=this->template observer<0>().begin()->second->P0;
            
            for (const auto& nIter : this->template observer<0>())
            {
                for(int d=0;d<dim;++d)
                {
                    if (nIter.second->P0(d)<_xMin(d))
                    {
                        _xMin(d)=nIter.second->P0(d);
                    }
                    if (nIter.second->P0(d)>_xMax(d))
                    {
                        _xMax(d)=nIter.second->P0(d);
                    }
                }
            }
            _xC=0.5*(_xMax+_xMin);
        }
        else
        {
            std::cout<<"Mesh is empty."<<std::endl;
        }
        std::cout<<"  xMin="<< std::setprecision(15)<<std::scientific<<_xMin.transpose()<<std::endl;
        std::cout<<"  xMax="<< std::setprecision(15)<<std::scientific<<_xMax.transpose()<<std::endl;
        std::cout<<"  xCenter="<< std::setprecision(15)<<std::scientific<<_xC.transpose()<<std::endl;
        
        // Populate typename SimplicialMesh<dim>::MeshRegionBoundaryContainerType
        regionBoundaries().clear();
        
        size_t bndSimplexCount=0;
        size_t rgnBndSimplexCount=0;
        for (const auto& simpl : this->template observer<dim-1>())
        {
            
            if(simpl.second->isBoundarySimplex())
            {// count number of bonudary simplices for later check
                bndSimplexCount++;
            }
            
            if(simpl.second->isRegionBoundarySimplex())
            {
                const auto regionIDset=simpl.second->regionIDs();
                std::pair<size_t,size_t> regionIDs(std::make_pair(*regionIDset.begin(),*regionIDset.rbegin()));
                const auto regionBndIter=regionBoundaries().find(regionIDs);
                if(regionBndIter!=regionBoundaries().end())
                {
                    regionBndIter->second.simplices().insert(simpl.second);
                }
                else
                {
                    regionBoundaries().emplace(regionIDs,regionIDs).first->second.simplices().insert(simpl.second);
                }
                rgnBndSimplexCount++;
            }
        }
        
        updateRegions();
        updateRegionBoundaries();
        identifyParallelFaces(periodicFaceIDs);
        
        // Check that all simplices have been accounted for
        size_t bndFaceSimplexSum=0;
        for(auto region : MeshRegionObserverType::regions())
        {// Sum number of external faces for final check
            std::cout<<magentaColor<<"MeshRegion "<<region.second->regionID<<defaultColor<<std::endl;
            std::cout<<"    simplices: "<<region.second->simplices().size()<<" Simplex<"<<dim<<","<<dim<<">"<<std::endl;
            for(auto& face : region.second->faces())
            {
                std::cout<<"    face "<<face.second->sID<<": size="<<face.second->size()<<",hullPts="<<face.second->convexHull().size()<<", outNormal "<<face.second->outNormal().transpose()<<std::endl;
                bndFaceSimplexSum+=face.second->size();
            }
        }
        
        size_t rgnBndFaceSimplexSum=0;
        for(const auto& rgnBnd : regionBoundaries())
        {// Sum number of internal faces for final check
            std::cout<<magentaColor<<"MeshRegionBoundary ("<<rgnBnd.second.regionBndID.first<<","<<rgnBnd.second.regionBndID.second<<")"<<defaultColor<<std::endl;
            std::cout<<"    simplices: "<<rgnBnd.second.simplices().size()<<" Simplex<"<<dim<<","<<dim-1<<">"<<std::endl;
            for(auto& face : rgnBnd.second.faces())
            {
                std::cout<<"    face "<<face.second->sID<<": hullPts="<<face.second->convexHull().size()<<", outNormal "<<face.second->outNormal().transpose()<<std::endl;
                rgnBndFaceSimplexSum+=face.second->size();
                bndFaceSimplexSum-=2*face.second->size(); // each region boundary face was added to two regions
            }
        }
        
        if(bndFaceSimplexSum!=bndSimplexCount)
        {
            std::cout<<"boundary simplices="<<bndSimplexCount<<std::endl;
            std::cout<<"simplices in external faces="<<bndFaceSimplexSum<<std::endl;
            throw std::runtime_error("WRONG NUMBER OF BOUNDAY FACE SIMPLICES");
        }
        
        if(rgnBndFaceSimplexSum!=rgnBndSimplexCount)
        {
            std::cout<<"region-boundary simplices="<<rgnBndSimplexCount<<std::endl;
            std::cout<<"simplices in internal faces="<<rgnBndFaceSimplexSum<<std::endl;
            throw std::runtime_error("WRONG NUMBER OF REGION-BOUNDARY FACE SIMPLICES");
        }
    }

    template<int dim>
    SimplicialMesh<dim>::SimplicialMesh() :
    /* init */ _xMin(Eigen::Matrix<double,dim,1>::Zero())
    /* init */,_xMax(Eigen::Matrix<double,dim,1>::Zero())
    /* init */,_xC(Eigen::Matrix<double,dim,1>::Zero())
    /* init */,vol0(0.0)
    {
    }

    template<int dim>
    SimplicialMesh<dim>::SimplicialMesh(const std::string& meshFileName,const Eigen::Matrix<double,dim,dim>& A,const Eigen::Matrix<double,dim,1>& x0,const std::set<int>& periodicFaceIDs) :
    /* init */ _xMin(Eigen::Matrix<double,dim,1>::Zero())
    /* init */,_xMax(Eigen::Matrix<double,dim,1>::Zero())
    /* init */,_xC(Eigen::Matrix<double,dim,1>::Zero())
    /* init */,vol0(0.0)
    {
        this->read(meshFileName,A,x0);
        createMesh(periodicFaceIDs);
        this->simplexReader().clear();
    }

    template<int dim>
    void SimplicialMesh<dim>::readMesh(const std::string& meshFileName,const Eigen::Matrix<double,dim,dim>& A,const Eigen::Matrix<double,dim,1>& x0,const std::set<int>& periodicFaceIDs)
    {
        _periodicSimplexPairs.clear();
        simplices().clear();
        this->read(meshFileName,A,x0);
        createMesh(periodicFaceIDs);
        this->simplexReader().clear();
    }

    template<int dim>
    const typename SimplicialMesh<dim>::SimplexMapType& SimplicialMesh<dim>::simplices() const
    {
        return *this;
    }

    template<int dim>
    typename SimplicialMesh<dim>::SimplexMapType& SimplicialMesh<dim>::simplices()
    {
        return *this;
    }

    template<int dim>
    void SimplicialMesh<dim>::updateRegions()
    {
        for(auto region : MeshRegionObserverType::regions())
        {
            region.second->update();
        }
    }

    template<int dim>
    void SimplicialMesh<dim>::updateRegionBoundaries()
    {
        for(auto& rgnBnd : regionBoundaries())
        {
            rgnBnd.second.update();
            MeshRegionType* const region1(this->region(rgnBnd.second.regionBndID.first));
            MeshRegionType* const region2(this->region(rgnBnd.second.regionBndID.second));
            for(const auto& face : rgnBnd.second.faces())
            {// add each face to the region boundary to the corresponding two regions
                region1->faces().emplace(face.second->sID,face.second);
                region2->faces().emplace(face.second->sID,face.second);
                
            }
        }
    }

    template<int dim>
    void SimplicialMesh<dim>::identifyParallelFaces(const std::set<int>& periodicFaceIDs)
    {
        std::vector<PlanarMeshFace<dim>*> faces;
        std::set<int> available;
        for(const auto& region : this->regions())
        {
            for(const auto& face : region.second->faces())
            {
                if(face.second->isExternal())
                {
                    faces.push_back(face.second.get());
                    available.insert(face.second->sID);
                }
            }
        }
        for(const auto id : periodicFaceIDs)
        {
            if(id != -1 && !available.count(id))
            {
                throw std::runtime_error("Unknown external periodic face ID " + std::to_string(id));
            }
        }
        const double tolerance = 1.e-8*(xMax()-xMin()).norm()
            + 64.0*std::numeric_limits<double>::epsilon()*std::max(xMin().norm(),xMax().norm());
        PeriodicMeshLattice<dim> lattice;
        if(!periodicFaceIDs.empty())
        {
            const auto full = discoverPeriodicMeshLattice(faces,tolerance);
            lattice = selectPeriodicMeshLattice(full,periodicFaceIDs,tolerance);
        }
        const auto& selected = lattice.pairs;
        PeriodicSimplexPairMap<dim> simplexPairs;
        for(const auto& entry : selected)
        {
            PeriodicSimplexPairMap<dim> matches;
            VectorDim shift;
            if(!isPeriodicFacePair<dim>(*entry.first,*entry.second.second,shift,tolerance,&matches) ||
               (shift-entry.second.first).norm()>tolerance)
            {
                throw std::runtime_error("Active periodic face has inconsistent simplex correspondence.");
            }
            simplexPairs.insert(matches.begin(),matches.end());
        }
        for(const auto& entry : simplexPairs)
        {
            const auto reverse=simplexPairs.find(entry.second.partner);
            if(reverse==simplexPairs.end() || reverse->second.partner!=entry.first ||
               (entry.second.shift+reverse->second.shift).norm()>tolerance)
            {
                throw std::runtime_error("Nonreciprocal periodic simplex correspondence.");
            }
            for(size_t i=0;i<dim;++i)
            {
                if(reverse->second.vertexPermutation[entry.second.vertexPermutation[i]]!=i)
                {
                    throw std::runtime_error("Nonreciprocal periodic vertex permutation.");
                }
            }
        }
        // Commit only after face and simplex validation have both succeeded.
        _periodicBasis = lattice.basis;
        _periodicSimplexPairs.swap(simplexPairs);
        for(auto* face : faces)
        {
            const auto found = selected.find(face);
            face->periodicFacePair = found == selected.end()
                ? std::make_pair(VectorDim(VectorDim::Zero()), static_cast<const PlanarMeshFace<dim>*>(nullptr))
                : std::make_pair(found->second.first, static_cast<const PlanarMeshFace<dim>*>(found->second.second));
        }
    }

    template<int dim>
    typename SimplicialMesh<dim>::PeriodicBasisType SimplicialMesh<dim>::periodicBasis() const
    {
        return _periodicBasis;
    }

    template<int dim>
    void SimplicialMesh<dim>::insertSimplex(const typename SimplexTraits<dim,dim>::SimplexIDType& xIN,const int& regionID)
    {/*!@param[in] xIN the (unsorted) array of mesh node IDs defining a simplex
      *\brief Inserts a Simplex<dim> into the mesh, with ID equal to the
      * sorted array xIN.
      */
        const typename SimplexTraits<dim,dim>::SimplexIDType xID(SimplexTraits<dim,dim>::sortID(xIN));
        const auto pair=simplices().emplace(std::piecewise_construct,
                                            std::make_tuple(xID),
                                            std::make_tuple(this,xID, regionID)
                                            );
        assert(pair.second);
        vol0+=pair.first->second.vol0;
    }

    template<int dim>
    std::pair<bool,const Simplex<dim,dim>*> SimplicialMesh<dim>::search(const Eigen::Matrix<double,dim,1>& P) const
    {/*!@param[in] P position to search for
      *\returns a pair, where:
      * -pair.first is a boolean indicating whether the
      * search succesfully found a Simplex<dim,dim> which includes P.
      * -pair.second is a pointer to the last Simplex<dim,dim> searched.
      *
      * By default the search starts at this->begin()->second
      */
        return searchWithGuess(P,&(simplices().begin()->second));
    }

    template<int dim>
    std::pair<bool,const Simplex<dim,dim>*> SimplicialMesh<dim>::searchWithGuess(const Eigen::Matrix<double,dim,1>& P,
                                                                                 const Simplex<dim,dim>* const guess) const
    {/*!@param[in] P position to search for
      * @param[in] guess Simplex* where the search starts
      *\returns a pair, where:
      * -pair.first is a boolean indicating whether the
      * search succesfully found a Simplex<dim,dim> which includes P.
      * -pair.second is a pointer to the last Simplex<dim,dim> searched.
      */
        
        std::set<const Simplex<dim,dim>*> searchSet;
        if(guess)
        {
            return searchWithGuess(true,P,guess,searchSet);
        }
        else
        {
            return searchWithGuess(true,P,&(simplices().begin()->second),searchSet);
        }
    }

    template<int dim>
    std::pair<bool,const Simplex<dim,dim>*> SimplicialMesh<dim>::searchRegion(const int& regionID,
                                                                              const Eigen::Matrix<double,dim,1>& P) const
    {/*!@param[in] P position to search for
      *\returns a pair, where:
      * -pair.first is a boolean indicating whether the
      * search succesfully found a Simplex<dim,dim> which includes P.
      * -pair.second is a pointer to the last Simplex<dim,dim> searched.
      *
      * By default the search starts at this->begin()->second
      */
        return searchRegionWithGuess(P,*this->region(regionID)->simplices().begin());
    }

    template<int dim>
    std::pair<bool,const Simplex<dim,dim>*> SimplicialMesh<dim>::searchRegionWithGuess(const Eigen::Matrix<double,dim,1>& P,
                                                                                       const Simplex<dim,dim>* const guess) const
    {/*!@param[in] P position to search for
      * @param[in] guess Simplex* where the search starts
      *\returns a pair, where:
      * -pair.first is a boolean indicating whether the
      * search succesfully found a Simplex<dim,dim> which includes P.
      * -pair.second is a pointer to the last Simplex<dim,dim> searched.
      */
        std::set<const Simplex<dim,dim>*> searchSet;
        return searchWithGuess(false,P,guess,searchSet);
    }

    template<int dim>
    std::pair<bool,const Simplex<dim,dim>*> SimplicialMesh<dim>::searchWithGuess(const bool& searchAllRegions,
                                                                                 const Eigen::Matrix<double,dim,1>& P,
                                                                                 const Simplex<dim,dim>* const guess,
                                                                                 std::set<const Simplex<dim,dim>*>& searchSet) const
    {/*!@param[in] P position to search for
      * @param[in] guess Simplex* where the search starts. If searchAllRegions=false, only the region of guess is searched
      *\returns a pair, where:
      * -pair.first is a boolean indicating whether the
      * search succesfully found a Simplex<dim,dim> which includes P.
      * -pair.second is a pointer to the last Simplex<dim,dim> searched.
      */
        
        searchSet.clear();
        std::pair<bool,const Simplex<dim,dim>*> lastSearched(false,NULL);
        guess->convexDelaunaynSearch(searchAllRegions,P,lastSearched,searchSet);
        checkSearch(searchAllRegions,P,guess,lastSearched);
        
        if(!lastSearched.first)
        {// if search was not successful, try to search neighbors of last searched Simplex
            const std::pair<bool,const Simplex<dim,dim>*> temp(lastSearched);
            const Eigen::Matrix<double,dim+1,1> bary=temp.second->pos2bary(P);
            // Do not overwrite a successful result while checking another face.
            for(int k=0;k<dim+1 && !lastSearched.first;++k)
            {
                if(bary(k)<=FLT_EPSILON)
                {
                    //                        for(typename Simplex<dim,dim-1>::ParentContainerType::const_iterator pIter=temp.second->child(k).parentBegin();
                    //                            /*                                                            */ pIter!=temp.second->child(k).parentEnd();++pIter)
                    //                        {
                    //                            if((*pIter)->region->regionID==temp.second->region->regionID || searchAllRegions)
                    //                            {
                    //                                (*pIter)->convexDelaunaynSearch(searchAllRegions,P,lastSearched,searchSet);
                    //                                if (lastSearched.first)
                    //                                {
                    //                                    break;
                    //                                }
                    //                            }
                    //                        }
                    for(const auto& pIter : temp.second->child(k).parents())
                    {
                        if(pIter.second->region->regionID==temp.second->region->regionID || searchAllRegions)
                        {
                            pIter.second->convexDelaunaynSearch(searchAllRegions,P,lastSearched,searchSet);
                            if (lastSearched.first)
                            {
                                break;
                            }
                        }
                    }
                }
            }
            
            if(searchAllRegions)
            {// if search is still unsuccessful, reset lastSearched to temp
                if(!lastSearched.first && !lastSearched.second->isBoundarySimplex())
                {
                    lastSearched=temp;
                }
            }
            else
            {// if search is still unsuccessful, reset lastSearched to temp
                if(!lastSearched.first && !lastSearched.second->isBoundarySimplex() && !lastSearched.second->isRegionBoundarySimplex())
                {
                    lastSearched=temp;
                }
            }
        }
        
        checkSearch(searchAllRegions,P,guess,lastSearched);

        if(searchAllRegions && !lastSearched.first)
        {
            // The global grain assembly can be non-convex: a directed walk
            // may hit an external face before reaching an interior point in
            // another grain. Search each convex grain before declaring failure.
            // Keep the original boundary simplex if the point is truly outside.
            for(const auto& regionEntry : this->regions())
            {
                std::set<const Simplex<dim,dim>*> regionSearchSet;
                const auto regional=searchWithGuess(false,P,
                    *regionEntry.second->simplices().begin(),regionSearchSet);
                searchSet.insert(regionSearchSet.begin(),regionSearchSet.end());
                if(regional.first)
                {
                    return regional;
                }
            }
        }
        return lastSearched;
    }

    template<int dim>
    void SimplicialMesh<dim>::checkSearch(const bool& searchAllRegions,
                                          const Eigen::Matrix<double,dim,1>& P,
                                          const Simplex<dim,dim>* const guess,
                                          const std::pair<bool,const Simplex<dim,dim>*>& lastSearched) const
    {
        if(searchAllRegions)
        {// Check that search was successful, or that it ended on a boundary (point outside)
            if(!(lastSearched.first || lastSearched.second->isBoundarySimplex()))
            {
                std::cout<<"P="<<std::setprecision(15)<<std::scientific<<P.transpose()<<std::endl;
                std::cout<<"guess="<<guess->xID<<std::endl;
                std::cout<<"lastSearched="<<lastSearched.second->xID<<std::endl;
                assert(0 && "SEARCH DID NOT END ON BOUNDARY SIMPLEX");
            }
        }
        else
        {// Check that search was successful, or that it ended on a boundary or region-boundary
            if(!(lastSearched.first || lastSearched.second->isBoundarySimplex() || lastSearched.second->isRegionBoundarySimplex()))
            {
                std::cout<<"P="<<std::setprecision(15)<<std::scientific<<P.transpose()<<std::endl;
                std::cout<<"guess="<<guess->xID<<std::endl;
                std::cout<<"lastSearched="<<lastSearched.second->xID<<std::endl;
                assert(0 && "SEARCH DID NOT END ON BOUNDARY/REGION-BOUNDARY SIMPLEX");
            }
        }
    }

    template<int dim>
    std::pair<bool,const Simplex<dim,dim>*> SimplicialMesh<dim>::isStrictlyInsideMesh(const Eigen::Matrix<double,dim,1>& P,
                                                                                      const Simplex<dim,dim>* const guess,
                                                                                      const double& tol) const
    {
        std::pair<bool,const Simplex<dim,dim>*> temp(searchWithGuess(P,guess));
        if(temp.first)
        {
            const Eigen::Matrix<double,dim+1,1> bary(temp.second->pos2bary(P));
            int kMin;
            const double baryMin(bary.minCoeff(&kMin));
            if (std::fabs(baryMin)<tol && guess->child(kMin).isBoundarySimplex()) // on a boundary face
            {
                temp.first=false;
            }
        }
        return temp;
    }

    template<int dim>
    std::pair<bool,const Simplex<dim,dim>*> SimplicialMesh<dim>::isOnMeshBoundary(const Eigen::Matrix<double,dim,1>& P, const Simplex<dim,dim>* const guess, const double& tol) const
    {
        std::pair<bool,const Simplex<dim,dim>*> temp(searchWithGuess(P,guess));
        if(temp.first)
        {
            const Eigen::Matrix<double,dim+1,1> bary(temp.second->pos2bary(P));
            int kMin;
            const double baryMin(bary.minCoeff(&kMin));
            if (!(std::fabs(baryMin)<tol && guess->child(kMin).isBoundarySimplex())) // not on a boundary face
            {
                temp.first=false;
            }
        }
        return temp;
    }

    template<int dim>
    const Eigen::Matrix<double,dim,1>& SimplicialMesh<dim>::xMin() const
    {
        return _xMin;
    }

    template<int dim>
    const double& SimplicialMesh<dim>::xMin(const int& k) const
    {
        return _xMin(k);
    }

    template<int dim>
    const Eigen::Matrix<double,dim,1>& SimplicialMesh<dim>::xMax() const
    {
        return _xMax;
    }

    template<int dim>
    const double& SimplicialMesh<dim>::xMax(const int& k) const
    {
        return _xMax(k);
    }

    template<int dim>
    const Eigen::Matrix<double,dim,1>& SimplicialMesh<dim>::xCenter() const
    {
        return _xC;
    }

    template<int dim>
    const double& SimplicialMesh<dim>::volume() const
    {
        return vol0;
    }

    template<int dim>
    const typename SimplicialMesh<dim>::MeshRegionBoundaryContainerType& SimplicialMesh<dim>::regionBoundaries() const
    {
        return *this;
    }

    template<int dim>
    typename SimplicialMesh<dim>::MeshRegionBoundaryContainerType& SimplicialMesh<dim>::regionBoundaries()
    {
        return *this;
    }

    template<int dim>
    const typename SimplicialMesh<dim>::MeshRegionBoundaryType& SimplicialMesh<dim>::regionBoundary(const int& i,const int& j) const
    {
        return regionBoundaries().at(std::make_pair(std::min(i,j),std::max(i,j)));
    }

    template class SimplicialMesh<3>;
}
#endif
