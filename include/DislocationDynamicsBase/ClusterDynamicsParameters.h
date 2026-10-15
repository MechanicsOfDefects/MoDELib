/* This file is part of GreatWhite, a library for integrating
 * MOOSE and MoDeLib.
 *
 * (c) 2018 Multiscale Materials Solutions, LLC
 * ALL RIGTHS REVERVED
 *   Prepared by 2018 Multiscale Materials Solutions
 *     under contract N. 126180
 *   with the NAVAL NUCLEAR LABORATORY
 *
 *  See COPYRIGHT for full restriction
*/


#ifndef model_ClusterDynamicsParameters_H_
#define model_ClusterDynamicsParameters_H_

#include <TerminalColors.h>
#include <tuple>
#include <vector>
//#include <DislocatedMaterialBase.h>
#include <Polycrystal.h>
#include <DislocationDynamicsBase.h>
//#include <TrialBase.h>
//#include <EvalFunction.h>
//#include <EvalExpression.h>
//#include <DislocationStress.h>

#ifndef MODELIB_CD_MSIZE
#define MODELIB_CD_MSIZE 4
#endif
#ifndef MODELIB_CD_ISIZE
#define MODELIB_CD_ISIZE 8
#endif

namespace model
{

template<int dim>
struct ClusterDynamicsParameters
{
    // Numbers of mobile and immobile cluster-dynamics species. They are compile-time constants,
    // set by the CMake options MODELIB_CD_MSIZE and MODELIB_CD_ISIZE, and must match the material file.
    static constexpr int mSize=MODELIB_CD_MSIZE;     // e.g. Cv, Ci, C2i, C3i
    static constexpr int iSize=MODELIB_CD_ISIZE;  // e.g. Nc, Na1, Na2, Na3, cv, ca1, ca2, ca3
    static_assert(iSize>=0 && iSize%2==0,"MODELIB_CD_ISIZE must be even: each immobile family carries a number density and a defect content");

    typedef Eigen::Matrix<double,dim,1> VectorDim;
    typedef Eigen::Matrix<double,dim,dim> MatrixDim;
    typedef typename Polycrystal<dim>::GrainContainerType GrainContainerType;

    // Materials parameters
    const double kB;
    const double T; // Studied temperature in K
    const double omega; // Atomic volume
    const double b; // Burgers vector
    const double G0; // dpa/s, Dose rate

    // Mobile Species
    const Eigen::Array<double,1,mSize> msVector;
    const Eigen::Array<double,1,mSize> msRelRelaxVol;
    const Eigen::Array<double,1,mSize> msEf;
    const Eigen::Matrix<double,mSize,dim*(dim+1)/2> msEm; // for each species (row) em11 em12 em13 em22 em23 em33
    const Eigen::Matrix<double,mSize,dim*(dim+1)/2> msD0; // diffusion pre-exponential coeff. for each species (row) D011 D012 D013 D022 D023 D033
    const std::map<size_t,std::vector<Eigen::Matrix<double,dim,dim>>> D;   // map<grain_ID, vector of diffusion coeff for each species>
    const std::map<size_t,std::vector<Eigen::Matrix<double,dim,dim>>> invD;
    const std::map<size_t,Eigen::Array<double,mSize,1>> detD;
    const Eigen::Array<double,1,mSize> Eb; // binding energy of each mobile cluster, for its dissociation. Not used for the single v and i
    const Eigen::Array<double,1,mSize> mobileMeanDiffusivity; // orientation-averaged diffusion coefficient of each mobile species, det(D)^(1/3) (first grain)
    const Eigen::Array<double,1,mSize> mobileCaptureRadius; // capture radius of each mobile species: sphere of one atomic volume for v and i, loop of |m| defects otherwise
    const int vacancyIndex; // index of the single vacancy among the mobile species, or -1
    const int interstitialIndex; // index of the single interstitial among the mobile species, or -1
    const Eigen::Array<double,1,mSize> polarityM; // p_m: +1 interstitial, -1 vacancy

    // Irradiation Production
    const Eigen::Matrix<double,1,mSize> msCascadeFractions;
    const double vacancyLoopCascadeFraction; // epsilon_vL = 1 - f.M_v (Eq. 42): fraction of the cascade vacancies born in immobile clusters
    const double interstitialLoopCascadeFraction; // epsilon_iL = 1 - f.M_i (Eq. 42): fraction of the cascade interstitials born in immobile clusters
    const Eigen::Matrix<double,mSize,1> G; // dpa/s, Effective dose rate

    // First-order reaction
    const double dislocationSinkDensity; // density of the dislocation network: sink of the mobile species and target of loop-network coalescence
    const Eigen::Array<double,1,mSize> mobileDislocationSinkEfficiency; // sink efficiency (bias) of the network for each mobile species
    const Eigen::Array<double,1,mSize> otherSinks; // sink strength of the network for each mobile species, efficiency*density
    const Eigen::Array<double,1,iSize> initialImmobileClusters;
    const std::map<std::pair<int,int>,double> reactionMap;
    const Eigen::Matrix<double,mSize,mSize> mobileDissociationRate; // E: dissociation of the mobile clusters, content space
    const Eigen::Matrix<double,mSize,mSize> mobileOtherSinkRate;    // diagonal -otherSinks*Dbar: background part of D, content space
    // Second-order reaction
    const std::vector<Eigen::Matrix<double,mSize,mSize>> mobileReactionRate; // R_2: second-order reactions between the mobile species, content space

    // Bias factors
    const Eigen::Array<double,iSize/2,mSize> immobileSpeciesBias; // capture bias Z0_km of each immobile family (row) for each mobile species (column), fitted factor included; the diffusion-anisotropy factor is applied in getLoopBias

    // Immobile Species
    const Eigen::Array<double,1,iSize/2> immobileSpeciesVector;
    const Eigen::Array<double,1,iSize/2> immobileSpeciesRelRelaxVol;
    const Eigen::Matrix<double,dim,iSize/2> immobileSpeciesBurgers;
    const Eigen::Array<double,1,iSize/2> immobileSpeciesBurgersMagnitude;
    const Eigen::Array<double,1,iSize/2> polarityF; // p_k of each immobile family: +1 interstitial, -1 vacancy
    const Eigen::Array<double,iSize/2,mSize> polarityFM; // p_k*p_m: +1 when family k and mobile species m have the same sign
    const Eigen::Array<double,1,iSize/2> Iv; // vacancy selector of the families (report I_v)
    const Eigen::Array<double,1,iSize/2> Ii; // interstitial selector of the families (report I_i)
//  const double pyramidRelaxationVolume; // volume change of one bi-pyramid cluster (strain only). Not needed: the volume change is -Omega per defect in both shapes
    const double shapeTransitionWidth; // width of the shape sigmoid, as a fraction of shapeTransitionMaxSize-shapeTransitionMinSize (w0, Eq. 33)
    const Eigen::Array<double,1,iSize/2> shapeTransitionMinSize; // defects per cluster below which a cluster is a bi-pyramid (m_min, Eq. 33)
    const Eigen::Array<double,1,iSize/2> shapeTransitionMaxSize; // defects per cluster above which a cluster is a loop (m_max, Eq. 33)
    const double emptyFamilyThreshold; // clusters or defects per atom below which a family is treated as absent

    // Nucleation of immobile clusters (optional keys of the material file; without them the number densities are constant)
    const Eigen::Array<double,1,iSize/2> nucleatedClusterSize; // defects in a newly nucleated cluster

    // Kinetics of the immobile clusters integrated by ImmobileSpeciesSolver (optional keys of the material file)
    const Eigen::Array<double,1,iSize/2> clusterBindingEnergy; // binding energy of a vacancy to the clusters of each vacancy family, for thermal emission (Eq. 53). Not used for the interstitial families
    const Eigen::Array<double,1,iSize/2> clusterDissolutionRate; // nu_a=nu_a0*exp(-E_a/kB/T) of each immobile family (report Eq. 40). Zero: the family does not dissolve
    const Eigen::Array<double,1,iSize/2> cLL; // rate coefficient of like-loop coalescence
    const Eigen::Array<double,1,iSize/2> cLN; // rate coefficient of loop-network coalescence
    const double kappaLL; // overlap coefficient of like-loop coalescence
    const double kappaLN; // overlap coefficient of loop-network coalescence
    const Eigen::Array<double,1,iSize/2> minimumImmobileClusterSize; // smallest immobile cluster of each family: largest mobile cluster of the same sign plus one. The shrink gate closes at this size
    const std::vector<std::tuple<int,int,double,int>> nucleationChannels; // reactions between mobile species whose product is not mobile: (species a, species b, rate coefficient in content space, 0 vacancy/1 interstitial)
    const Eigen::Array<double,iSize/2,mSize> loopBias; // capture bias Z_km of each family (row) for each mobile species (column): immobileSpeciesBias times the diffusion-anisotropy factor


    // Reaction map (types: parameters)
    const bool computeReactions;


    // Discrete Loop Generation
    const double discretizationFactor;
    const double discretizationTime;
    const double minimumLoopSize;
    const Eigen::Array<double,2,mSize> discreteDislocationBias; // bias of the discrete dislocation loops for each mobile species (first row vacancy loops, second row interstitial loops), used only by the discrete dislocation segments
    
    ClusterDynamicsParameters(const DislocationDynamicsBase<dim>& ddBase /*const Polycrystal<dim>& poly*/);
    std::map<std::pair<int,int>,double> getMap(const Eigen::Array<double,mSize*(mSize+1)/2,3> matrix_in) const;
    Eigen::Matrix<double,mSize,mSize> getMobileDissociationRate() const;
    Eigen::Matrix<double,mSize,mSize> getMobileOtherSinkRate() const;
    std::vector<Eigen::Matrix<double,mSize,mSize>> getMobileReactionRate() const;
    static Eigen::Array<double,1,iSize/2> getOptionalFamilyArray(const std::string& materialFile,const std::string& key,const double& defaultValue);
    static int getOptionalInt(const std::string& materialFile,const std::string& key,const int& defaultValue);
    static double getOptionalDouble(const std::string& materialFile,const std::string& key,const double& defaultValue);
    Eigen::Array<double,iSize/2,mSize> getLoopBias() const;
    Eigen::Array<double,1,mSize> getMobileMeanDiffusivity() const;
    Eigen::Array<double,1,mSize> getMobileCaptureRadius() const;
    int getMobileSpeciesIndex(const double& m) const;
    std::vector<std::tuple<int,int,double,int>> getNucleationChannels() const;
    bool isBasalFamily(const int& k) const;
    Eigen::Array<double,1,iSize/2> getImmobileSpeciesBurgersMagnitude(const GrainContainerType& grains) const;
    std::map<size_t,std::vector<Eigen::Matrix<double,dim,dim>>> getD(const GrainContainerType& grains) const;
    std::vector<Eigen::Matrix<double,dim,dim>> getDlocal() const;
    std::map<size_t,std::vector<Eigen::Matrix<double,dim,dim>>> getInvD() const;
    std::map<size_t,Eigen::Array<double,mSize,1>> getDetD() const;
    Eigen::Array<double,1,iSize> getInitialImmobileClusters(const Eigen::Array<double,1,iSize> initialImmobileClusters_SI, const double b_SI) const;
    Eigen::Array<double,1,mSize> equilibriumMobileConcentration(const double& stressTrace) const;
    Eigen::Array<double,1,mSize> dislocationMobileConcentration(const VectorDim& b,const VectorDim& t,const VectorDim& fPK,const MatrixDim& stress) const;
    Eigen::Array<double,1,mSize> boundaryMobileConcentration(const double& stressTrace,const double& normalTraction) const;
    // Shape of the clusters (Eqs. 30-33), as functions of the number of defects per cluster m of each family
    Eigen::Array<double,1,iSize/2> shapeSigmoid(const Eigen::Array<double,1,iSize/2>& m) const;
    Eigen::Array<double,1,iSize/2> pyramidRadius(const Eigen::Array<double,1,iSize/2>& m) const;
    Eigen::Array<double,1,iSize/2> loopRadius(const Eigen::Array<double,1,iSize/2>& m) const;
    Eigen::Array<double,1,iSize/2> effectiveRadiusDerivative(const Eigen::Array<double,1,iSize/2>& m) const;
    Eigen::Array<double,1,iSize/2> shapeInterpolation(const Eigen::Array<double,1,iSize/2>& m, const Eigen::Array<double,1,iSize/2>& lowValue, const Eigen::Array<double,1,iSize/2>& highValue) const;
    Eigen::Array<double,dim,dim> shapeInterpolation(const Eigen::Array<double,1,iSize/2>& m, const Eigen::Array<double,dim,dim>& lowValue, const Eigen::Array<double,dim,dim>& highValue, const int& k) const;
    Eigen::Array<bool,1,iSize/2> populated(const Eigen::Array<double,1,iSize/2>& n, const Eigen::Array<double,1,iSize/2>& c) const;
    Eigen::Array<double,1,iSize/2> effectiveClusterRadius(const Eigen::Array<double,1,iSize/2>& n, const Eigen::Array<double,1,iSize/2>& c) const;
    Eigen::Array<double,1,iSize/2> captureLength(const Eigen::Array<double,1,iSize/2>& n, const Eigen::Array<double,1,iSize/2>& c) const;

    // Former sigmoidal interpolation and plotting helpers (not compiled), kept for reference
//  Eigen::Array<double,1,iSize/2> sigmoidalVectorInterpolation(const Eigen::Array<double,1,iSize/2>& CI, const Eigen::Array<double,1,iSize/2>& N, const Eigen::Array<double,1,iSize/2>& lowValue, const Eigen::Array<double,1,iSize/2>& highValue) const;
//  Eigen::Array<double,dim,dim> sigmoidalMatrixInterpolation(const Eigen::Array<double,1,iSize/2>& CI, const Eigen::Array<double,1,iSize/2>& N, const Eigen::Array<double,dim,dim>& lowValue, const Eigen::Array<double,dim,dim>& highValue, const int& index) const;
//  Eigen::Array<double,1,iSize/2> sigmoidalPlotVectorInterpolation(const Eigen::Array<double,1,iSize/2>& CI, const Eigen::Array<double,1,iSize/2>& N, const Eigen::Array<double,1,iSize/2>& lowValue, const Eigen::Array<double,1,iSize/2>& highValue) const;
//  Eigen::Array<double,1,iSize/2> clusterPlotRadius(const Eigen::Array<double,1,iSize/2>& CI, const Eigen::Array<double,1,iSize/2>& N) const;
};

}
#endif
