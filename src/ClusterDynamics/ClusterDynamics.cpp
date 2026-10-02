/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ClusterDynamics_cpp_
#define model_ClusterDynamics_cpp_

#ifdef _OPENMP
#include <omp.h>
#endif

#include <random>
#include <array>
#include <numeric>
#include <ClusterDynamics.h>
//#include <ExternalAndInternalBoundary.h>
//#include <Fix.h>
#include <DislocationNetwork.h>
#include <MicrostructureGenerator.h>

namespace model
{

    //template <int dim>
    //struct BoundaryConcentration
    //{
    //    const IrradiationClimbParameters<dim>& icp;
    //    const Eigen::Matrix<double,dim,dim>& stress;
    //
    //    /**********************************************************************/
    //    BoundaryConcentration(const IrradiationClimbParameters<dim>& icp_in,const Eigen::Matrix<double,dim,dim>& stress_in) :
    //    /* */ icp(icp_in)
    //    /* */,stress(stress_in)
    //    {
    //
    //    }
    //
    //    /**********************************************************************/
    //    template <typename NodeType,int dofPerNode>
    //    Eigen::Matrix<double,dofPerNode,1> operator()(const NodeType& node,
    //                                                  Eigen::Matrix<double,dofPerNode,1>& val) const
    //    {
    //
    //        const auto v1(icp.equilibriumClusterConcentration(stress.trace()));
    //
    //        const auto outNormal(node.outNormal()); // used to compute traction
    //        const double tn = outNormal.dot(stress*outNormal);
    //
    //        val=(v1.array()*exp(-tn*icp.omega/icp.kB/icp.T*icp.speciesVector)).matrix().transpose();
    //
    //        return val;
    //
    //    }
    //
    //};

template <int dim>
typename ClusterDynamics<dim>::UniformControllerContainerType ClusterDynamics<dim>::getUniformControllers(const DislocationDynamicsBase<dim>& ddBase,const ClusterDynamicsParameters<dim>& cdp)
{
    typedef typename UniformControllerType::MatrixVoigt MatrixVoigt;
    typedef typename UniformControllerType::VectorVoigt VectorVoigt;
    
    UniformControllerContainerType temp;
    for(const auto& pair : cdp.D)
    {
        MatrixVoigt C(MatrixVoigt::Zero());
        const auto& grainID(pair.first);
        for(size_t k=0;k<pair.second.size();++k)
        {
            const auto& Dm(pair.second[k]);
            C.template block<dim,dim>(k*dim,k*dim)=Dm;
        }
        
        const auto f0(VectorVoigt::Zero()); // flux
        const auto f0Dot(VectorVoigt::Zero()); // flux rate
        const auto g0(VectorVoigt::Zero()); // Concentration gradient
        const auto g0Dot(VectorVoigt::Zero()); // Concentration gradient rate
        const auto stiffnessRatio(VectorVoigt::Ones()*1.0e20); // impose concentration gradient
        const double t0(0.0);

        temp.emplace(grainID,new UniformControllerType(t0,C,stiffnessRatio,g0,g0Dot,f0,f0Dot));
    }
    return temp;
}


    template<int dim>
    ClusterDynamics<dim>::ClusterDynamics(MicrostructureContainerType& mc) :
    /* init */ MicrostructureBase<dim>("ClusterDynamics",mc)
    /* init */,cdp(this->microstructures.ddBase)
//    /* init */,ClusterDynamicsBase<dim>(this->microstructures.ddBase)
    // /* init */,useClusterDynamicsFEM(this->microstructures.ddBase.fe? bool(TextFileParser(this->microstructures.ddBase.simulationParameters.traitsIO.ddFile).readScalar<int>("useClusterDynamicsFEM",true)) : false )
    /* init */,useClusterDynamicsFEM(this->microstructures.ddBase.fe? bool(TextFileParser(this->microstructures.ddBase.simulationParameters.traitsIO.inputFilesFolder+"/ClusterDynamics.txt").readScalar<int>("useClusterDynamicsFEM",true)) : false )
//    /* init */,useClusterDynamicsFEM(bool(TextFileParser(this->microstructures.ddBase.simulationParameters.traitsIO.ddFile).readScalar<int>("useClusterDynamicsFEM",true)))
    /* init */,clusterDynamicsFEM(useClusterDynamicsFEM?  new ClusterDynamicsFEM<dim>(this->microstructures.ddBase,cdp) : nullptr)
    /* init */,uniformControllers(useClusterDynamicsFEM? UniformControllerContainerType() : getUniformControllers(this->microstructures.ddBase,this->cdp))
    /* init */,discreteLoopsInitialized(false)
//    /* init */,nodeListInternalExternal(this->microstructures.ddBase.isPeriodicDomain ? -1 : this->microstructures.ddBase.fe->template createNodeList<ExternalAndInternalBoundary>())
//    /* init */,mobileClustersIncrement(this->microstructures.ddBase.fe->template trial<'d',mSize>())
//    /* init */,dV(this->microstructures.ddBase.fe->template domain<EntireDomain,dVorder,GaussLegendre>())
//    /* init */,mBWF((test(this->mobileGrad),-this->microstructures.ddBase.poly.Omega*this->mobileFlux)*dV)
//    /* init */,dmBWF((test(grad(mobileClustersIncrement)),-this->microstructures.ddBase.poly.Omega*(FluxMatrix<dim>(this->cdp)*grad(mobileClustersIncrement)))*dV)
//    /* init */,mSolver(true,FLT_EPSILON)
////    /* init */,mSolver(mBWF,true,FLT_EPSILON)
//    /* init */,cascadeGlobalProduction(((test(clusterDynamicsFEM->mobileClusters),make_constant(this->cdp.G))*dV).globalVector())
    ///* init */,cascadeGlobalProduction(((test(clusterDynamicsFEM->mobileClusters),make_constant(this->cdp.G.transpose().eval()))*dV).globalVector())
    {
//        mobileClustersIncrement.setConstant(Eigen::Matrix<double,mSize,1>::Zero());    
//        std::cout<<cascadeGlobalProduction<<std::endl;
    }

    

    template<int dim>
    void ClusterDynamics<dim>::initializeConfiguration(const DDconfigIO<dim>& configIO,const std::ofstream& f_file,const std::ofstream& F_labels)
    {
        this->lastUpdateTime=this->microstructures.ddBase.simulationParameters.totalTime;
        
        if(clusterDynamicsFEM)
        {
            clusterDynamicsFEM->initializeConfiguration(configIO,f_file,F_labels);
        }
        else
        {// nothing to initialize for uniformControllers
            
        }
    }

template<int dim>
void ClusterDynamics<dim>::applyBoundaryConditions()
{
 
    const auto& nodesInternalExternal(clusterDynamicsFEM->mobileClusters.fe().nodeList(clusterDynamicsFEM->nodeListInternalExternal));
    
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(size_t k=0;k<nodesInternalExternal.size();++k)
    {
        const auto& node(nodesInternalExternal[k]);
        const auto outNormal(node->outNormal()); // used to compute traction
        const MatrixDim sigma(this->microstructures.stress(node->P0,node,nullptr,nullptr));
        const double normalTraction(outNormal.dot(sigma*outNormal));
        const auto bndConcentration(cdp.boundaryMobileConcentration(sigma.trace(),normalTraction));
        
        VectorMSize otherConcentration(VectorMSize::Zero());
        for(const auto& microstructure : this->microstructures)
        {
            if(microstructure.get()!=static_cast<const MicrostructureBase<dim>* const>(this))
            {// not the ClusterDynamics physics
                otherConcentration += microstructure->mobileConcentration(node->P0,node,nullptr,nullptr);
            }
        }
        for(int k=0;k<mSize;++k)
        {
            clusterDynamicsFEM->mobileClusters.dirichletConditions().at(mSize*node->gID+k) = bndConcentration(k) - otherConcentration(k);
        }
    }
}

    template<int dim>
    void ClusterDynamics<dim>::solve()
    {
        if(clusterDynamicsFEM)
        {
            if(!clusterDynamicsFEM->solverInitialized)
            {
                clusterDynamicsFEM->initializeSolver();
            }
            std::cout<<", mobile BCs"<<std::flush;
            applyBoundaryConditions();

            auto DN(this->microstructures.template getUniqueTypedMicrostructure<DislocationNetwork<dim>>());
            const bool hasDiscreteLoops(DN? (DN->loops().size()>0? true : false) : false);
            clusterDynamicsFEM->solve(hasDiscreteLoops); 
        }
        else
        {

        }
    }

    template<int dim>
    void ClusterDynamics<dim>::reSolve()
    {
        solve();
    }

    template<int dim>
    void ClusterDynamics<dim>::updateConfiguration()
    {
        if(clusterDynamicsFEM)
        {
            auto DN(this->microstructures.template getUniqueTypedMicrostructure<DislocationNetwork<dim>>());
            const bool hasDiscreteLoops(DN? (DN->loops().size()>0 ? true : false) : false);
            const double dt(this->microstructures.getDt());

            if constexpr (iSize > 0)
            {
                if(!hasDiscreteLoops)
                {
                    clusterDynamicsFEM->updateImmobileClusters(dt);
                }

                if(   DN
                   && !hasDiscreteLoops
                   && !discreteLoopsInitialized
                   && cdp.discretizationTime <= this->microstructures.ddBase.simulationParameters.totalTime)
                {// replace the immobile fields by discrete loops. Needs DislocationDynamics in the physics
                    initializeDiscreteClimbLoops();
                }
            }
        }
        this->lastUpdateTime=this->microstructures.ddBase.simulationParameters.totalTime;
    }

    template<int dim>
    double ClusterDynamics<dim>::getDt() const
    {
        return this->microstructures.ddBase.simulationParameters.dtMax;
    }

    template<int dim>
    void ClusterDynamics<dim>::output(DDconfigIO<dim>& configIO,DDauxIO<dim>& auxIO,std::ofstream& f_file,std::ofstream& F_labels) const
    {
        if(clusterDynamicsFEM)
        {
            const size_t nNodes(clusterDynamicsFEM->mobileClusters.fe().nodes().size());
            configIO.cdMatrix().resize(nNodes,mSize+iSize);
            configIO.cdMatrix().block(0,0,nNodes,mSize)=clusterDynamicsFEM->mobileClusters.dofVector().reshaped(mSize,nNodes).transpose();
            configIO.cdMatrix().block(0,mSize,nNodes,iSize)=clusterDynamicsFEM->immobileClusters.dofVector().reshaped(iSize,nNodes).transpose();
        }
        else
        {
            
        }
    }

    template<int dim>
    typename ClusterDynamics<dim>::VectorDim ClusterDynamics<dim>::inelasticDisplacementRate(const VectorDim& x, const NodeType* const node, const ElementType* const ele,const SimplexDim* const guess) const
    {
        if(clusterDynamicsFEM)
        {
            return clusterDynamicsFEM->inelasticDisplacementRate(x,node,ele,guess);
        }
        else
        {
            return VectorDim::Zero();
        }
    }

    template<int dim>
    typename ClusterDynamics<dim>::MatrixDim ClusterDynamics<dim>::averagePlasticDistortion() const
    {
        if(clusterDynamicsFEM)
        {
            return clusterDynamicsFEM->averageBetaP()+clusterDynamicsFEM->averageBetaV();
        }
        else
        {
            return MatrixDim::Zero();
        }
    }

    template<int dim>
    typename ClusterDynamics<dim>::MatrixDim ClusterDynamics<dim>::averagePlasticDistortionRate() const
    {
        return MatrixDim::Zero();
    }

    template<int dim>
    typename ClusterDynamics<dim>::VectorDim ClusterDynamics<dim>::displacement(const VectorDim&,const NodeType* const,const ElementType* const,const SimplexDim* const) const
    {
        return VectorDim::Zero();
    }

    template<int dim>
    typename ClusterDynamics<dim>::MatrixDim ClusterDynamics<dim>::stress(const VectorDim&,const NodeType* const,const ElementType* const,const SimplexDim* const) const
    {
        return MatrixDim::Zero();
    }

    template<int dim>
    typename ClusterDynamics<dim>::VectorMSize ClusterDynamics<dim>::mobileConcentration(const VectorDim& x,const NodeType* const node,const ElementType* const ele,const SimplexDim* const guess) const
    {
        if(uniformControllers.size())
        {
            const auto pointGrains(this->pointGrains(x,node,ele,guess));
            VectorMSize averageVal(VectorMSize::Zero());
            for(const auto& grain : pointGrains)
            {
                const auto grainID(grain->grainID);
                const auto& controller(uniformControllers.at(grainID));
                averageVal += controller->grad(this->microstructures.ddBase.simulationParameters.totalTime).reshaped(dim,mSize).transpose()*x;
            }
            return averageVal/pointGrains.size();
        }
        else
        {
            if(node)
            {
                return eval(clusterDynamicsFEM->mobileClusters)(*node);
            }
            else
            {
                if(ele)
                {
                    return eval(clusterDynamicsFEM->mobileClusters)(*ele,ele->simplex.pos2bary(x));
                }
                else
                {
                    return eval(clusterDynamicsFEM->mobileClusters)(x,guess);
                }
            }
        }
    }

    template<int dim>
    typename ClusterDynamics<dim>::VectorISize ClusterDynamics<dim>::immobileClusters(const VectorDim& x,const NodeType* const node,const ElementType* const ele,const SimplexDim* const guess) const
    {
        if(node)
        {
            return eval(clusterDynamicsFEM->immobileClusters)(*node);
        }
        else
        {
            if(ele)
            {
                return eval(clusterDynamicsFEM->immobileClusters)(*ele,ele->simplex.pos2bary(x));
            }
            else
            {
                return eval(clusterDynamicsFEM->immobileClusters)(x,guess);
            }
        }
    }


    template<int dim>
    typename ClusterDynamics<dim>::MatrixDim ClusterDynamics<dim>::averageStress() const
    {
        return MatrixDim::Zero();
    }


    /*!\brief A discrete loop drawn from the immobile fields, before it is inserted in the network
     */
    template<int dim>
    struct ClusterDynamicsDiscreteLoop
    {
        Eigen::Matrix<double,dim,1> P; // center
        double r;                      // radius of the disc that stores its defects
        int merged;                    // number of drawn loops that it holds
    };

    /*!\brief Merges the overlapping loops of one family in one grain, conserving their total area.
     *
     * Two loops merge when the distance of their centers is less than the sum
     * of their radii and the centers lie within one plane spacing of a common
     * habit plane: loops on different parallel planes cannot become one loop
     * without climbing. The merged loop has the area of the two, and its center
     * is their area-weighted centroid.
     *
     * Merging at fixed total area lowers the number of loops N, and the ratio of
     * the loop diameter to the loop spacing grows as N^(-1/6). A pass that would
     * make the mean diameter larger than the mean spacing is therefore rejected,
     * and the population is reported as saturated.
     *
     * \returns true if the population is saturated
     */
    template<int dim>
    bool coalesceDiscreteLoops(std::vector<ClusterDynamicsDiscreteLoop<dim>>& loops,
                               const Eigen::Matrix<double,dim,1>& unitNormal,
                               const double& planeSpacing,
                               const double& volume)
    {
        const int maxPasses(20);
        for(int pass=0;pass<maxPasses;++pass)
        {
            if(loops.size()<2)
            {
                break;
            }

            double rMax(0.0);
            for(const auto& loop : loops)
            {
                rMax=std::max(rMax,loop.r);
            }
            if(rMax<=0.0)
            {
                break;
            }

            // Bin the centers in cells of size 2*rMax, so that overlapping loops are in neighbor cells
            const double cellSize(2.0*rMax);
            std::map<std::array<long int,3>,std::vector<size_t>> cells;
            const auto cellOf=[&](const Eigen::Matrix<double,dim,1>& P)
            {
                return std::array<long int,3>{(long int)std::floor(P(0)/cellSize),(long int)std::floor(P(1)/cellSize),(long int)std::floor(P(2)/cellSize)};
            };
            for(size_t i=0;i<loops.size();++i)
            {
                cells[cellOf(loops[i].P)].push_back(i);
            }

            // Union-find over the overlapping pairs
            std::vector<size_t> parent(loops.size());
            std::iota(parent.begin(),parent.end(),0);
            const auto find=[&](size_t i)
            {
                while(parent[i]!=i)
                {
                    parent[i]=parent[parent[i]];
                    i=parent[i];
                }
                return i;
            };

            size_t nPairs(0);
            for(size_t i=0;i<loops.size();++i)
            {
                const auto ci(cellOf(loops[i].P));
                for(long int dx=-1;dx<=1;++dx)
                {
                    for(long int dy=-1;dy<=1;++dy)
                    {
                        for(long int dz=-1;dz<=1;++dz)
                        {
                            const auto cellIter(cells.find(std::array<long int,3>{ci[0]+dx,ci[1]+dy,ci[2]+dz}));
                            if(cellIter!=cells.end())
                            {
                                for(const size_t& j : cellIter->second)
                                {
                                    if(j>i)
                                    {
                                        const Eigen::Matrix<double,dim,1> dP(loops[i].P-loops[j].P);
                                        if(   dP.norm()<loops[i].r+loops[j].r
                                           && std::fabs(unitNormal.dot(dP))<=planeSpacing)
                                        {
                                            const size_t ri(find(i));
                                            const size_t rj(find(j));
                                            if(ri!=rj)
                                            {
                                                parent[std::max(ri,rj)]=std::min(ri,rj);
                                                nPairs++;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if(nPairs==0)
            {
                break;
            }

            // Merge each group into one loop
            std::map<size_t,ClusterDynamicsDiscreteLoop<dim>> groups;
            for(size_t i=0;i<loops.size();++i)
            {
                const size_t root(find(i));
                const double area(loops[i].r*loops[i].r); // area/pi
                auto groupIter(groups.find(root));
                if(groupIter==groups.end())
                {
                    groups.emplace(root,ClusterDynamicsDiscreteLoop<dim>{loops[i].P*area,area,loops[i].merged});
                }
                else
                {
                    groupIter->second.P+=loops[i].P*area;
                    groupIter->second.r+=area;
                    groupIter->second.merged+=loops[i].merged;
                }
            }
            std::vector<ClusterDynamicsDiscreteLoop<dim>> candidates;
            double rSum(0.0);
            for(auto& pair : groups)
            {
                const double area(pair.second.r);
                pair.second.P/=area;
                pair.second.r=std::sqrt(area);
                rSum+=pair.second.r;
                candidates.push_back(pair.second);
            }

            const double spacing(std::pow(volume/candidates.size(),1.0/3.0));
            if(2.0*rSum/candidates.size()>spacing)
            {// the merged loops would be wider than their own spacing
                return true;
            }
            loops=candidates;
        }
        return false;
    }

    template<int dim>
    void ClusterDynamics<dim>::initializeDiscreteClimbLoops()
    {/*! Replaces the immobile cluster fields by discrete prismatic loops in
      * the DislocationNetwork, and removes from the fields what the loops hold.
      *
      * For each family and each grain:
      * - Number. The number of loops is the integral of the number density,
      *   divided by clusterDiscretizationFactor (the number of clusters that one
      *   discrete loop stands for; 1 gives the population of the fields).
      * - Position. The loops are placed in the elements with a probability
      *   proportional to the number of clusters in each element, so that the
      *   spatial variation of the field is kept.
      * - Size. A loop holds the defects of the clusters it stands for. Its
      *   radius r is that of the disc that stores them, pi*r^2*|b|/Omega. The radii
      *   are then scaled, so that the loops hold exactly the defects of the field.
      * - Coalescence. Overlapping loops in the same habit plane are merged at
      *   constant area (see coalesceDiscreteLoops).
      * - Shape. Each loop is a polygon whose area is that of the disc.
      * - Loops smaller than minimumLoopSize, and loops that do not fit in their
      *   grain, are not created. Their share of the defects stays in the fields:
      *   number density and content are scaled by the same factor, which leaves
      *   the cluster size unchanged.
      *
      * The number of loops, the stored defects and the sink strength before and
      * after the conversion are printed for each family.
      */
        if constexpr (iSize == 0)
        {
            return;
        }
        else
        {
            discreteLoopsInitialized=true;
            auto DN(this->microstructures.template getUniqueTypedMicrostructure<DislocationNetwork<dim>>());
            if(!DN || !clusterDynamicsFEM)
            {
                return;
            }
            std::cout<<greenBoldColor<<"\nClusterDynamics: converting the immobile fields to discrete loops"<<defaultColor<<std::endl;

            constexpr int nF(iSize/2);
            const int loopSides(12);
            const double polygonFactor(std::sqrt(2.0*M_PI/(loopSides*std::sin(2.0*M_PI/loopSides)))); // circumradius of the polygon/radius of the disc of equal area
            const double lumping(std::max(cdp.discretizationFactor,1.0));
            const int maxFitAttempts(30);
            const auto& ddBase(this->microstructures.ddBase);
            typedef ClusterDynamicsDiscreteLoop<dim> DiscreteLoopType;
            typedef Eigen::Matrix<double,dim+1,1> BaryType;

            // Quadrature points of the tetrahedron, exact for the quadratic fields
            const double qa(0.5854101966249685);
            const double qb(0.1381966011250105);

            MicrostructureGenerator mg(this->microstructures.ddBase);
            std::mt19937 generator(0);
            std::exponential_distribution<double> expDistribution(1.0);

            for(int k=0;k<nF;++k)
            {
                const double storedPerArea(M_PI*cdp.b*cdp.immobileSpeciesBurgersMagnitude(k)/cdp.omega); // defects stored per r^2

                // Clusters, defects and sink strength held by the field in each element
                std::vector<const ElementType*> elements;
                std::vector<double> elementNumber;
                std::vector<double> elementDefects;
                std::map<size_t,double> grainVolume;
                std::map<size_t,VectorDim> grainCenter;
                double fieldNumber(0.0);
                double fieldDefects(0.0);
                double fieldSink(0.0);
                for(const auto& elePair : ddBase.fe->elements())
                {
                    const auto& ele(elePair.second);
                    double Ne(0.0);
                    double Ce(0.0);
                    double Se(0.0);
                    for(int q=0;q<dim+1;++q)
                    {
                        BaryType bary(BaryType::Constant(qb));
                        bary(q)=qa;
                        const Eigen::Matrix<double,iSize,1> value(eval(clusterDynamicsFEM->immobileClusters)(ele,bary));
                        const double Nq(std::max(value(k),0.0));
                        const double cq(std::max(value(nF+k),cdp.n_min(k)*cdp.omega*Nq)); // same lower bound as in the sink strength
                        Ne+=Nq;
                        Ce+=cq/cdp.omega;
                        if(Nq>0.0)
                        {
                            Se+=2.0*M_PI*std::sqrt(cq/cdp.omega/Nq/storedPerArea)*Nq;
                        }
                    }
                    const double weight(ele.simplex.vol0/(dim+1));
                    elements.push_back(&ele);
                    elementNumber.push_back(Ne*weight);
                    elementDefects.push_back(Ce*weight);
                    fieldNumber+=Ne*weight;
                    fieldDefects+=Ce*weight;
                    fieldSink+=Se*weight;

                    const size_t grainID(ele.simplex.region->regionID);
                    grainVolume[grainID]+=ele.simplex.vol0;
                    auto centerIter(grainCenter.find(grainID));
                    if(centerIter==grainCenter.end())
                    {
                        grainCenter.emplace(grainID,ele.simplex.center()*ele.simplex.vol0);
                    }
                    else
                    {
                        centerIter->second+=ele.simplex.center()*ele.simplex.vol0;
                    }
                }
                for(auto& pair : grainCenter)
                {
                    pair.second/=grainVolume.at(pair.first);
                }

                // Draw the loops
                const long int nDraw(std::llround(fieldNumber/lumping));
                std::map<size_t,std::vector<DiscreteLoopType>> grainLoops;
                if(nDraw>0 && fieldDefects>0.0)
                {
                    std::discrete_distribution<size_t> elementDistribution(elementNumber.begin(),elementNumber.end());
                    for(long int i=0;i<nDraw;++i)
                    {
                        const size_t e(elementDistribution(generator));
                        const auto& ele(*elements[e]);
                        BaryType bary;
                        for(int q=0;q<dim+1;++q)
                        {// uniform point in the element
                            bary(q)=expDistribution(generator);
                        }
                        bary/=bary.sum();
                        const double defectsPerLoop(elementDefects[e]/elementNumber[e]*lumping);
                        grainLoops[ele.simplex.region->regionID].push_back(DiscreteLoopType{ele.position(bary),std::sqrt(defectsPerLoop/storedPerArea),1});
                    }
                }

                // Scale the radii, so that the loops hold the defects of the field
                double drawnArea(0.0);
                double drawnRadius(0.0);
                for(const auto& pair : grainLoops)
                {
                    for(const auto& loop : pair.second)
                    {
                        drawnArea+=loop.r*loop.r;
                        drawnRadius+=loop.r;
                    }
                }
                const double radiusScale(drawnArea>0.0? std::sqrt(fieldDefects/storedPerArea/drawnArea) : 1.0);

                // Coalesce and insert
                size_t coalescedNumber(0);
                size_t insertedNumber(0);
                size_t refusedSmall(0);
                size_t refusedFit(0);
                double coalescedRadius(0.0);
                double insertedDefects(0.0);
                double insertedSink(0.0);
                bool saturated(false);
                for(auto& pair : grainLoops)
                {
                    const size_t& grainID(pair.first);
                    auto& loops(pair.second);
                    for(auto& loop : loops)
                    {
                        loop.r*=radiusScale;
                    }

                    const auto& grain(ddBase.poly.grain(grainID));
                    const VectorDim b(grain->latticeBasis*cdp.immobileSpeciesBurgers.col(k));
                    const ReciprocalLatticeDirection<dim> n(grain->reciprocalLatticeDirection(b));

                    saturated=coalesceDiscreteLoops(loops,b.normalized(),n.planeSpacing(),grainVolume.at(grainID)) || saturated;
                    coalescedNumber+=loops.size();

                    for(const auto& loop : loops)
                    {
                        coalescedRadius+=loop.r;
                        if(loop.r<cdp.minimumLoopSize)
                        {
                            refusedSmall++;
                            continue;
                        }

                        VectorDim center(loop.P);
                        bool inserted(false);
                        for(int attempt=0;attempt<maxFitAttempts && !inserted;++attempt)
                        {
                            const long int planeIndex(n.closestPlaneIndexOfPoint(center));
                            GlidePlaneKey<dim> glidePlaneKey(planeIndex,n);
                            std::shared_ptr<PeriodicGlidePlane<dim>> periodicGlidePlane(this->microstructures.ddBase.periodicGlidePlaneFactory.get(glidePlaneKey));
                            const VectorDim P0(periodicGlidePlane->referencePlane->snapToPlane(center));
                            std::vector<VectorDim> pts;
                            for(int j=0;j<loopSides;++j)
                            {
                                const Eigen::Matrix<double,dim-1,1> pLocal(std::cos(2.0*M_PI/loopSides*j),std::sin(2.0*M_PI/loopSides*j));
                                pts.push_back(periodicGlidePlane->referencePlane->globalPosition(polygonFactor*loop.r*pLocal)-periodicGlidePlane->referencePlane->P+P0);
                            }
                            VectorDim nA(VectorDim::Zero());
                            for(size_t j=0;j<pts.size();++j)
                            {
                                nA+=0.5*(pts[j]-pts.front()).cross(pts[(j+1)%pts.size()]-pts[j]);
                            }
                            // b along the area normal for vacancy loops, opposite to it for interstitial loops
                            const VectorDim loopBurgers((b.dot(nA)*cdp.immobileSpeciesVector(k)<0.0)? b : (-b).eval());
                            inserted=mg.insertJunctionLoop(pts,periodicGlidePlane,loopBurgers,periodicGlidePlane->referencePlane->unitNormal,P0,grainID,DislocationLoopIO<dim>::SESSILELOOP);
                            if(!inserted)
                            {// some nodes are outside the grain: move the loop towards the center of the grain
                                center+=0.1*(grainCenter.at(grainID)-center);
                            }
                        }
                        if(inserted)
                        {
                            insertedNumber++;
                            insertedDefects+=storedPerArea*loop.r*loop.r;
                            insertedSink+=2.0*M_PI*loop.r;
                        }
                        else
                        {
                            refusedFit++;
                        }
                    }
                }

                // Remove from the field what the loops hold
                const double fractionInserted(fieldDefects>0.0? std::min(insertedDefects/fieldDefects,1.0) : 0.0);
                if(insertedNumber>0)
                {
                    clusterDynamicsFEM->scaleImmobileFamily(k,1.0-fractionInserted);
                }

                const double meshVolume(ddBase.mesh.volume());
                const auto diameterOverSpacing=[&](const double& radiusSum,const double& number)
                {
                    return number>0.0? 2.0*radiusSum/number/std::pow(meshVolume/number,1.0/3.0) : 0.0;
                };
                std::cout<<"  family "<<k<<" (b="<<cdp.immobileSpeciesBurgers.col(k).transpose()<<", "<<(cdp.immobileSpeciesVector(k)<0.0? "vacancy" : "interstitial")<<")\n"
                /*     */<<"    loops:         field "<<fieldNumber<<", drawn "<<nDraw<<", after coalescence "<<coalescedNumber<<", inserted "<<insertedNumber
                /*     */<<" (below minimumLoopSize "<<refusedSmall<<", not fitting in the grain "<<refusedFit<<")\n"
                /*     */<<"    defects:       field "<<fieldDefects<<", inserted "<<insertedDefects<<", left in the field "<<fieldDefects-insertedDefects<<"\n"
                /*     */<<"    sink strength: field "<<fieldSink<<", inserted "<<insertedSink<<" [b]\n"
                /*     */<<"    diameter/spacing: drawn "<<diameterOverSpacing(drawnRadius*radiusScale,nDraw)<<", after coalescence "<<diameterOverSpacing(coalescedRadius,coalescedNumber)
                /*     */<<(saturated? " (saturated: the loops are wider than their spacing)" : "")<<std::endl;
            }

            if(mg.config().loops().size())
            {
                DN->addConfiguration(mg.config());
            }
            std::cout<<"  "<<mg.config().loops().size()<<" loops added to the DislocationNetwork"<<std::endl;
        }
    }



    template struct ClusterDynamics<3>;
}
#endif


//template<int dim>
//void ClusterDynamics<dim>::solveDiffusiveDisplacement()
//{
//    std::cout<<"Solving diffusiveDisplacementRate..."<<std::flush;
//    const auto t0= std::chrono::system_clock::now();
//    diffusiveDisplacementRate=Eigen::VectorXd::Zero(this->diffusiveDisplacement.gSize());
//    for(const auto& node: this->microstructures.ddBase.fe->nodes())
//    {
//        const Eigen::Matrix<double,dim*mSize,1> speciesFlux(eval(this->mobileFlux)(node));
//        Eigen::Matrix<double,dim,1> netFlux(Eigen::Matrix<double,dim,1>::Zero());
//        for(int i=0; i<mSize; i++)
//        {
//            const int mSgn(this->cdp.msVector(i)/std::abs(this->cdp.msVector(i)));
//            netFlux+= speciesFlux.template block<dim,1>(i*dim,0)*mSgn;
//        }
//        diffusiveDisplacementRate.template segment<dim>(dim*node.gID)=netFlux;
//    }
//    std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
//
//}
