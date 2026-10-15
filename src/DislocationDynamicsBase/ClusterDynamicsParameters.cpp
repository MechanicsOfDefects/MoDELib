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


#ifndef model_ClusterDynamicsParameters_cpp_
#define model_ClusterDynamicsParameters_cpp_

#include <ClusterDynamicsParameters.h>

namespace model
{

    template<int dim>
    ClusterDynamicsParameters<dim>::ClusterDynamicsParameters(const DislocationDynamicsBase<dim>& ddBase) :
    /* init */ kB(ddBase.poly.kB),
    /* init */ T(ddBase.poly.T),
    /* init */ omega(ddBase.poly.Omega), // atomic volume of the lattice (Polycrystal::getAtomicVolume), in b^3
    /* init */ b(ddBase.poly.b),
    /* init */ G0(true? TextFileParser(ddBase.poly.materialFile).readScalar<double>("doseRate_dpaPerSec",true)*(ddBase.poly.b_SI/ddBase.poly.cs_SI) : 0.0),
    /* MOBILE SPECIES */
    /* init */ msVector(true? TextFileParser(ddBase.poly.materialFile).readMatrix<int>("mobileSpeciesVector",1,mSize,true).array().template cast<double>().eval() : Eigen::Array<double,1,mSize>::Zero()),
    /* init */ msRelRelaxVol(true? TextFileParser(ddBase.poly.materialFile).readMatrix<double>("mobileSpeciesRelRelaxVol",1,mSize,true).array().eval() : Eigen::Array<double,1,mSize>::Zero()),
    /* init */ msEf(true? (TextFileParser(ddBase.poly.materialFile).readMatrix<double>("mobileSpeciesEnergyFormation_eV",1,mSize,true).array()*ddBase.poly.eV2J/ddBase.poly.mu_SI/pow(ddBase.poly.b_SI,3)).eval() : Eigen::Array<double,1,mSize>::Zero()),
    /* init */ msEm(true? (TextFileParser(ddBase.poly.materialFile).readMatrix<double>("mobileSpeciesEnergyMigration_eV",mSize,dim*(dim+1)/2,true)*ddBase.poly.eV2J/ddBase.poly.mu_SI/pow(ddBase.poly.b_SI,3)).eval() : Eigen::Matrix<double,mSize,dim*(dim+1)/2>::Zero()),
    /* init */ msD0(true? (TextFileParser(ddBase.poly.materialFile).readMatrix<double>("mobileSpeciesD0_SI",mSize,dim*(dim+1)/2,true)/ddBase.poly.b_SI/ddBase.poly.cs_SI).eval() : Eigen::Matrix<double,mSize,dim*(dim+1)/2>::Zero()),
    /* init */ D(getD(ddBase.poly.grains)),
    /* init */ invD(getInvD()),
    /* init */ detD(getDetD()),
    /* init */ Eb((true && iSize>0) ? TextFileParser(ddBase.poly.materialFile).readMatrix<double,1,mSize>("mobileSpeciesBindingEnergy_eV",true).array()*ddBase.poly.eV2J/ddBase.poly.mu_SI/pow(ddBase.poly.b_SI,3) : Eigen::Array<double,1,mSize>::Zero().eval()),
    /* init */ mobileMeanDiffusivity(getMobileMeanDiffusivity()),
    /* init */ mobileCaptureRadius(getMobileCaptureRadius()),
    /* init */ vacancyIndex(getMobileSpeciesIndex(-1.0)),
    /* init */ interstitialIndex(getMobileSpeciesIndex(1.0)),
    /* init */ polarityM(msVector.sign()),
    /* init */ msCascadeFractions((mSize>1 && G0>0.0) ? TextFileParser(ddBase.poly.materialFile).readMatrix<double>("mobileSpeciesCascadeFractions",1,mSize,true) : Eigen::Matrix<double,1,mSize>::Ones()),
    /* init */ vacancyLoopCascadeFraction(std::max(1.0-(msCascadeFractions.array()*(polarityM<0.0).template cast<double>()).sum(),0.0)),
    /* init */ interstitialLoopCascadeFraction(std::max(1.0-(msCascadeFractions.array()*(polarityM>0.0).template cast<double>()).sum(),0.0)),
    /* init */ G(G0*msCascadeFractions), // P_M = G f (Eq. 3): G the nominal dpa rate, f the fractions of the cascade defects born as each mobile species
    /* init */ dislocationSinkDensity(TextFileParser(ddBase.poly.materialFile).readScalar<double>("dislocationSinkDensity_SI",true)*ddBase.poly.b_SI*ddBase.poly.b_SI),
    /* init */ mobileDislocationSinkEfficiency(TextFileParser(ddBase.poly.materialFile).readMatrix<double,1,mSize>("mobileDislocationSinkEfficiency",true).array()),
    /* init */ otherSinks((TextFileParser(ddBase.poly.materialFile).readScalar<double>("dislocationSinkDensity_SI",true)*mobileDislocationSinkEfficiency*ddBase.poly.b_SI*ddBase.poly.b_SI).eval()),
    /* init */ initialImmobileClusters(iSize>0 ? getInitialImmobileClusters(TextFileParser(ddBase.poly.materialFile).readMatrix<double,1,iSize>("initialImmobileClusters_SI",true),ddBase.poly.b_SI) : Eigen::Array<double,1,iSize>::Zero().eval()), // read only with immobile species
    /* init */ reactionMap((true && mSize>1) ? getMap(TextFileParser(ddBase.poly.materialFile).readMatrix<double,mSize*(mSize+1)/2,3>("reactionPrefactorMap",true)) : std::map<std::pair<int,int>,double>()),
    /* init */ mobileDissociationRate(mSize>1 ? getMobileDissociationRate() : Eigen::Matrix<double,mSize,mSize>::Zero().eval()),
    /* init */ mobileOtherSinkRate(mSize>1 ? getMobileOtherSinkRate() : Eigen::Matrix<double,mSize,mSize>::Zero().eval()),
    /* init */ mobileReactionRate(mSize>1 ? getMobileReactionRate() : std::vector<Eigen::Matrix<double,mSize,mSize>>(2,Eigen::Matrix<double,mSize,mSize>::Zero())),
    /* init */ immobileSpeciesBias(iSize>0 ? TextFileParser(ddBase.poly.materialFile).readMatrix<double,iSize/2,mSize>("immobileSpeciesBias",true).eval() : Eigen::Array<double,iSize/2,mSize>::Zero()),
    /* IMMOBILE SPECIES */
    /* init */ immobileSpeciesVector((true && iSize>0) ? TextFileParser(ddBase.poly.materialFile).readMatrix<int>("immobileSpeciesVector",1,iSize/2,true).array().template cast<double>() : Eigen::Array<double,1,iSize/2>::Zero().eval()),
    /* init */ immobileSpeciesRelRelaxVol((true && iSize>0) ? TextFileParser(ddBase.poly.materialFile).readMatrix<double,1,iSize/2>("immobileSpeciesRelRelaxVol",true).array() : Eigen::Array<double,1,iSize/2>::Zero().eval()),
    /* init */ immobileSpeciesBurgers((true && iSize>0) ? TextFileParser(ddBase.poly.materialFile).readMatrix<double,dim,iSize/2>("immobileSpeciesBurgers",true) : Eigen::Matrix<double,dim,iSize/2>::Zero()),
    /* init */ immobileSpeciesBurgersMagnitude((true && iSize>0) ? getImmobileSpeciesBurgersMagnitude(ddBase.poly.grains) : Eigen::Array<double,1,iSize/2>::Zero().eval()),
    /* init */ polarityF(immobileSpeciesVector.sign()),
    /* init */ polarityFM((polarityF.transpose().matrix()*polarityM.matrix()).array()),
    /* init */ Iv((polarityF<0.0).template cast<double>()),
    /* init */ Ii((polarityF>0.0).template cast<double>()),
//  /* init */ pyramidRelaxationVolume((true && iSize>0) ? TextFileParser(ddBase.poly.materialFile).readScalar<double>("pyramidRelaxationVolume_SI",true)/pow(ddBase.poly.b_SI,3) : 0.0),
    /* init */ shapeTransitionWidth((true && iSize>0) ? TextFileParser(ddBase.poly.materialFile).readScalar<double>("shapeTransitionWidth",true) : 0.0),
    /* init */ shapeTransitionMinSize((true && iSize>0) ? TextFileParser(ddBase.poly.materialFile).readMatrix<double,1,iSize/2>("shapeTransitionMinSize",true) : Eigen::Array<double,1,iSize/2>::Zero()),
    /* init */ shapeTransitionMaxSize((true && iSize>0) ? TextFileParser(ddBase.poly.materialFile).readMatrix<double,1,iSize/2>("shapeTransitionMaxSize",true) : Eigen::Array<double,1,iSize/2>::Zero()),
    /* init */ emptyFamilyThreshold(1.0e-30),
    /* init */ nucleatedClusterSize(getOptionalFamilyArray(ddBase.poly.materialFile,"nucleatedClusterSize",1.0)),
    /* init */ clusterBindingEnergy(iSize>0 ? (TextFileParser(ddBase.poly.materialFile).readMatrix<double,1,iSize/2>("clusterBindingEnergy_eV",true).array()*ddBase.poly.eV2J/ddBase.poly.mu_SI/pow(ddBase.poly.b_SI,3)).eval() : Eigen::Array<double,1,iSize/2>::Zero().eval()),
    /* init */ clusterDissolutionRate(getOptionalFamilyArray(ddBase.poly.materialFile,"clusterDissolutionRatePrefactor",0.0)*(ddBase.poly.b_SI/ddBase.poly.cs_SI)
    /*      */        *(-getOptionalFamilyArray(ddBase.poly.materialFile,"clusterDissolutionEnergy",0.0)*ddBase.poly.eV2J/ddBase.poly.mu_SI/pow(ddBase.poly.b_SI,3)/kB/T).exp()),
    /* init */ cLL(getOptionalFamilyArray(ddBase.poly.materialFile,"cLL",0.0)),
    /* init */ cLN(getOptionalFamilyArray(ddBase.poly.materialFile,"cLN",0.0)),
    /* init */ kappaLL(getOptionalDouble(ddBase.poly.materialFile,"kappaLL",0.0)),
    /* init */ kappaLN(getOptionalDouble(ddBase.poly.materialFile,"kappaLN",0.0)),
    /* init */ minimumImmobileClusterSize((polarityF<0.0).select(Eigen::Array<double,1,iSize/2>::Constant(std::max(-msVector.minCoeff(),0.0)+1.0),Eigen::Array<double,1,iSize/2>::Constant(std::max(msVector.maxCoeff(),0.0)+1.0))),
    /* init */ nucleationChannels(getNucleationChannels()), // clusters are always born from the reactions whose product is not mobile; set their reactionPrefactorMap entries to zero to turn it off
    /* init */ loopBias(getLoopBias()),
    /* init */ computeReactions((true && mSize>0 && mSize+iSize>1)? TextFileParser(ddBase.poly.materialFile).readScalar<int>("computeReactions",true) : 0),
    /* DISCRETIZATION */
    /* init */ discretizationFactor((true && iSize>0) ? TextFileParser(ddBase.poly.materialFile).readScalar<double>("clusterDiscretizationFactor",true) : 0.0),
    /* init */ discretizationTime((true && iSize>0) ? TextFileParser(ddBase.poly.materialFile).readScalar<double>("clusterDiscretizationTime",true)*(ddBase.poly.cs_SI/ddBase.poly.b_SI) : 0.0),
    /* init */ minimumLoopSize((true && iSize>0) ? TextFileParser(ddBase.poly.materialFile).readScalar<double>("minimumLoopSize",true)/ddBase.poly.b_SI : 0.0 ),
    /* init */ discreteDislocationBias(true? TextFileParser(ddBase.poly.materialFile).readMatrix<double,2,mSize>("discreteDislocationBias",true).eval() : Eigen::Array<double,2,mSize>::Zero())

    {
        
        double vacancySum(0.0);
        double interstitialSum(0.0);
        for(int k=0;k<mSize;++k)
        {
            if(msVector(k)<0)
            {// vacancy cluster
                vacancySum+=msCascadeFractions(k);
            }
            else if(msVector(k)>0)
            {// interstitial cluster
                interstitialSum+=msCascadeFractions(k);
            }
        }

        if(vacancySum>1.0+FLT_EPSILON)
        {
            std::cout<<redBoldColor<<"Warning: vacancy cluster fractions sum to more than one."<<defaultColor<<std::endl;
        }
        else if(vacancySum<1.0-FLT_EPSILON)
        {// the remainder is born in immobile clusters
            std::cout<<"vacancy cluster fractions sum to "<<vacancySum<<": the remainder "<<1.0-vacancySum<<" is born in immobile clusters (epsilon_vL, Eq. 42)."<<std::endl;
        }
        if(interstitialSum>1.0+FLT_EPSILON)
        {
            std::cout<<redBoldColor<<"Warning: interstitial cluster fractions sum to more than one."<<defaultColor<<std::endl;
        }
        else if(interstitialSum<1.0-FLT_EPSILON)
        {// the remainder is born in immobile clusters
            std::cout<<"interstitial cluster fractions sum to "<<interstitialSum<<": the remainder "<<1.0-interstitialSum<<" is born in immobile clusters (epsilon_iL, Eq. 42)."<<std::endl;
        }
        for(const auto& pair: reactionMap)
        {
            const auto key=pair.first;
            std::cout<<"find the interaction between "<<key.first<<" - "<<key.second<<" and parameter: "<<pair.second<<std::endl;
        }

        std::cout<<"mobile dissociation rate, content space (in 1/s): "<<std::endl;
        std::cout<<mobileDissociationRate*ddBase.poly.cs_SI/ddBase.poly.b_SI<<std::endl;
        std::cout<<"mobile other-sink rate, content space (in 1/s): "<<std::endl;
        std::cout<<mobileOtherSinkRate*ddBase.poly.cs_SI/ddBase.poly.b_SI<<std::endl;
        
        for(size_t k=0; k<mSize; k++)
        {
            std::cout<<"second-order interaction matrix, content space (in 1/s) for "<<static_cast<int>(msVector(k))<<"-species is: "<<std::endl;
            std::cout<<mobileReactionRate[k]*ddBase.poly.cs_SI/ddBase.poly.b_SI<<std::endl;
        }

        if(iSize>0)
        {
            std::cout<<"immobile cluster nucleation: epsilon_vL="<<vacancyLoopCascadeFraction<<", epsilon_iL="<<interstitialLoopCascadeFraction<<", nucleatedClusterSize="<<nucleatedClusterSize<<std::endl;
            for(const auto& channel : nucleationChannels)
            {
                std::cout<<"  clustering channel "<<static_cast<int>(msVector(std::get<0>(channel)))<<" + "<<static_cast<int>(msVector(std::get<1>(channel)))<<", rate coefficient, content space (in 1/s) "<<std::get<2>(channel)*ddBase.poly.cs_SI/ddBase.poly.b_SI<<std::endl;
            }
        }
    }

    template<int dim>
    std::map<std::pair<int,int>,double> ClusterDynamicsParameters<dim>::getMap(const Eigen::Array<double,mSize*(mSize+1)/2,3> matrix_in) const
    {
        std::map<std::pair<int,int>,double> tempMap;
        for(int k=0; k<mSize*(mSize+1)/2; k++)
        {
            const int key0= matrix_in(k,0)<0? static_cast<int>(matrix_in(k,0)-msVector(0)): static_cast<int>(matrix_in(k,0)-msVector(0)-1);
            const int key1= matrix_in(k,1)<0? static_cast<int>(matrix_in(k,1)-msVector(0)): static_cast<int>(matrix_in(k,1)-msVector(0)-1);
            if(key0>key1)
            {
                throw std::runtime_error("Please keep the reaction map orders, key0<=key1.");
            }
            const std::pair<int,int> key(key0,key1);
            tempMap.insert(std::pair<std::pair<int,int>,double>(key,matrix_in(k,2)));
        }
        return tempMap;
    }

    template<int dim>
    Eigen::Matrix<double,ClusterDynamicsParameters<dim>::mSize,ClusterDynamicsParameters<dim>::mSize> ClusterDynamicsParameters<dim>::getMobileDissociationRate() const
    {/*!\returns the dissociation (emission) matrix E of the mobile clusters, in content space:
      * the entry (i,j) is the rate at which the content of species j goes to species i.
      */
        const int vIndex(vacancyIndex);
        const int iIndex(interstitialIndex);
        if(iIndex-vIndex!=1 && msVector.matrix().squaredNorm()>0)
        { // Expected sequence should look like ... -2,-1,+1,+2,+3 ...
            throw std::runtime_error("intersitital Index must be after vacancy Index.");
        }
        
        const Eigen::Array<double,1,mSize>& aveD(mobileMeanDiffusivity);
        
        const Eigen::Array<double,1,mSize>& rn(mobileCaptureRadius);
                
        Eigen::Array<double,1,mSize> alpha(Eigen::Array<double,1,mSize>::Zero());
        if(vIndex>0)
        { // Dissociation rate for v-clusters [Christien and Barbu 2005 JNM 346]
            for(int k=0;k<vIndex;k++)
            {
                auto it=reactionMap.find(std::pair<int,int>(k+1,vIndex));
                if(it!=reactionMap.end())
                {
                    alpha(k)=it->second*2.0*M_PI*rn(k)*aveD(vIndex)/omega*exp(-Eb(k)/kB/T);
                }
            }
        }
                
        if(mSize-1-iIndex>1)
        { // Dissociation rate for i-clusters
            for(int k=iIndex+1;k<mSize;k++)
            {
                auto it=reactionMap.find(std::pair<int,int>(iIndex,k-1));
                if(it!=reactionMap.end())
                {
                    alpha(k)=it->second*8.0*M_PI*rn(iIndex)*aveD(iIndex)/omega*exp(-Eb(k)/kB/T);
                }
            }
        }
        
        Eigen::Matrix<double,mSize,mSize> tempDR(Eigen::Matrix<double,mSize,mSize>::Zero());
        for(int k=0;k<mSize;k++)
        {
            tempDR(k,k) += (-alpha(k));
            if(k>0 && k<=vIndex)
            {
                tempDR(k,k-1) += alpha(k-1)*std::fabs(msVector(k))/std::fabs(msVector(k-1)); // content space
                tempDR(vIndex,k-1) += alpha(k-1)/std::fabs(msVector(k-1)); // content space
            }
            else if(k>=iIndex && k<mSize-1)
            {
                tempDR(k,k+1) += alpha(k+1)*std::fabs(msVector(k))/std::fabs(msVector(k+1)); // content space
                tempDR(iIndex,k+1) += alpha(k+1)/std::fabs(msVector(k+1)); // content space
            }
            
        }
        
        return tempDR;
    }

    template<int dim>
    Eigen::Matrix<double,ClusterDynamicsParameters<dim>::mSize,ClusterDynamicsParameters<dim>::mSize> ClusterDynamicsParameters<dim>::getMobileOtherSinkRate() const
    {/*!\returns the diagonal matrix of the absorption of the mobile species at the
      * sinks other than the immobile clusters, -otherSinks_m*Dbar_m. Diagonal, so the
      * same in content and cluster space.
      */
        Eigen::Matrix<double,mSize,mSize> temp(Eigen::Matrix<double,mSize,mSize>::Zero());
        for(int k=0;k<mSize;k++)
        {
            temp(k,k)=-otherSinks(k)*mobileMeanDiffusivity(k);
        }
        return temp;
    }

    template<int dim>
    std::vector<Eigen::Matrix<double,ClusterDynamicsParameters<dim>::mSize,ClusterDynamicsParameters<dim>::mSize>> ClusterDynamicsParameters<dim>::getMobileReactionRate() const
    {

        const Eigen::Array<double,1,mSize>& aveD(mobileMeanDiffusivity);
        
        const Eigen::Array<double,1,mSize>& rn(mobileCaptureRadius);
        
        std::vector<Eigen::Matrix<double,mSize,mSize>> tempR2;
        Eigen::Matrix<double,mSize,mSize> R2sum(Eigen::Matrix<double,mSize,mSize>::Zero());
        double maxCoeff(0.0);
        
        for(size_t k=0;k<mSize;++k)
        {
            Eigen::Matrix<double,mSize,mSize> tempR2_k(Eigen::Matrix<double,mSize,mSize>::Zero());
            for(const auto& pair : reactionMap)
            {
                const auto key=pair.first;
                const double contentFactor(std::fabs(msVector(k))/std::fabs(msVector(key.first))/std::fabs(msVector(key.second))); // content space

                if(int(k)==key.first || int(k)==key.second)
                {// (key.first, key.second) consume k_specie
                    tempR2_k(key.first ,key.second) -= contentFactor*pair.second*4.0*M_PI*(rn(key.first)+rn(key.second))*(aveD(key.first)+aveD(key.second))/omega;
                    tempR2_k(key.second ,key.first) -= contentFactor*pair.second*4.0*M_PI*(rn(key.first)+rn(key.second))*(aveD(key.first)+aveD(key.second))/omega;
                }
                else if( fabs(msVector(key.first)+msVector(key.second)-msVector(k))<FLT_EPSILON )
                {// (key.first, key.second) generate k_specie
                    const double same= (key.first==key.second) ? 0.5 : 1.0;
                    tempR2_k(key.first ,key.second) += contentFactor*same*pair.second*4.0*M_PI*(rn(key.first)+rn(key.second))*(aveD(key.first)+aveD(key.second))/omega;
                    tempR2_k(key.second ,key.first) += contentFactor*same*pair.second*4.0*M_PI*(rn(key.first)+rn(key.second))*(aveD(key.first)+aveD(key.second))/omega;
                }
            }
            
            tempR2.push_back(tempR2_k);
            const double maxCoeffLocal(std::max(std::abs(tempR2.back().maxCoeff()),std::abs(tempR2.back().minCoeff())));
            
            if((tempR2.back()-tempR2.back().transpose()).norm()/maxCoeffLocal>FLT_EPSILON)
            {
                throw std::runtime_error("ClusterDynamicsParameters: R2_"+std::to_string(k)+" is not symmetric.");
            }
            maxCoeff=std::max(maxCoeff,maxCoeffLocal);
            R2sum+=tempR2.back()*(msVector(k)>0.0? 1.0 : -1.0); // defect conservation in content space
        }
        
        if(R2sum.norm()/maxCoeff>FLT_EPSILON)
        {
            std::cout<<redBoldColor<<"Warning: Sum of R2 is not zero."<<defaultColor<<std::endl;
        }
        
        return tempR2;
    }

    template<int dim>
    Eigen::Array<double,1,ClusterDynamicsParameters<dim>::iSize/2> ClusterDynamicsParameters<dim>::getOptionalFamilyArray(const std::string& materialFile,const std::string& key,const double& defaultValue)
    {/*!\returns the per-family values of an optional key of the material file,
      * or defaultValue for every family when the key is absent.
      */
        if(iSize>0)
        {
            Eigen::Matrix<double,1,iSize/2> temp;
            try
            {
                temp=TextFileParser(materialFile).readMatrix<double,1,iSize/2>(key,true);
            }
            catch(const std::runtime_error& e)
            {
                if(std::string(e.what()).find("does not cointain line with format")==std::string::npos)
                {// the key is present, and wrong
                    throw;
                }
                return Eigen::Array<double,1,iSize/2>::Constant(defaultValue);
            }
            return temp.array();
        }
        else
        {
            return Eigen::Array<double,1,iSize/2>::Constant(defaultValue);
        }
    }

    template<int dim>
    int ClusterDynamicsParameters<dim>::getOptionalInt(const std::string& materialFile,const std::string& key,const int& defaultValue)
    {
        try
        {
            return TextFileParser(materialFile).readScalar<int>(key,true);
        }
        catch(const std::runtime_error& e)
        {
            if(std::string(e.what()).find("does not cointain line with format")==std::string::npos)
            {// the key is present, and wrong
                throw;
            }
            return defaultValue;
        }
    }

    template<int dim>
    double ClusterDynamicsParameters<dim>::getOptionalDouble(const std::string& materialFile,const std::string& key,const double& defaultValue)
    {
        try
        {
            return TextFileParser(materialFile).readScalar<double>(key,true);
        }
        catch(const std::runtime_error& e)
        {
            if(std::string(e.what()).find("does not cointain line with format")==std::string::npos)
            {// the key is present, and wrong
                throw;
            }
            return defaultValue;
        }
    }

    template<int dim>
    Eigen::Array<double,ClusterDynamicsParameters<dim>::iSize/2,ClusterDynamicsParameters<dim>::mSize> ClusterDynamicsParameters<dim>::getLoopBias() const
    {/*!\returns the capture bias Z_km = Z0_km*ZDAD_k(p_m) of each immobile family k
      * for each mobile species m.
      * Z0 is the row of immobileSpeciesBias of the family (fitted factor included).
      * ZDAD is the diffusion-anisotropy factor, with p_m=(D_c/D_a)^(1/6):
      * p_m for clusters in the basal plane, (p_m+p_m^-2)/2 for clusters in a
      * plane that contains the c axis. It is 1 for isotropic diffusion.
      */
        Eigen::Array<double,iSize/2,mSize> temp(Eigen::Array<double,iSize/2,mSize>::Ones());
        if(iSize>0)
        {
            const std::vector<Eigen::Matrix<double,dim,dim>> Dlocal(getDlocal());
            for(int k=0;k<iSize/2;++k)
            {
                for(int m=0;m<mSize;++m)
                {
                    const double p(pow(Dlocal[m](2,2)/Dlocal[m](0,0),1.0/6.0));
                    const double ZDAD(isBasalFamily(k)? p : 0.5*(p+1.0/(p*p)));
                    temp(k,m)=immobileSpeciesBias(k,m)*ZDAD;
                }
            }
        }
        return temp;
    }

    template<int dim>
    bool ClusterDynamicsParameters<dim>::isBasalFamily(const int& k) const
    {/*!\returns true if the Burgers vector of family k is along the third lattice vector
      * (the c axis of a HEX crystal), so that the clusters lie in the basal plane.
      */
        return fabs(immobileSpeciesBurgers(2,k))>FLT_EPSILON
        /*  */ && fabs(immobileSpeciesBurgers(0,k))<FLT_EPSILON
        /*  */ && fabs(immobileSpeciesBurgers(1,k))<FLT_EPSILON;
    }

    template<int dim>
    Eigen::Array<double,1,ClusterDynamicsParameters<dim>::mSize> ClusterDynamicsParameters<dim>::getMobileMeanDiffusivity() const
    {/*!\returns the orientation-averaged diffusion coefficient of each mobile
      * species, det(D)^(1/3), in the first grain. Zero without diffusion coefficients.
      */
        Eigen::Array<double,1,mSize> temp(Eigen::Array<double,1,mSize>::Zero());
        if(!detD.empty())
        {
            for(int k=0;k<mSize;k++)
            {
                temp(k)=pow(detD.begin()->second(k),1.0/3.0);
            }
        }
        return temp;
    }

    template<int dim>
    Eigen::Array<double,1,ClusterDynamicsParameters<dim>::mSize> ClusterDynamicsParameters<dim>::getMobileCaptureRadius() const
    {/*!\returns the capture radius of each mobile species: a sphere of one atomic
      * volume for the single vacancy and interstitial, a loop of |m| defects for a cluster.
      */
        Eigen::Array<double,1,mSize> temp(Eigen::Array<double,1,mSize>::Zero());
        for(int k=0;k<mSize;k++)
        {
            if(fabs(msVector(k)-1)<FLT_EPSILON || fabs(msVector(k)+1)<FLT_EPSILON)
            {
                temp(k)=pow(3.0*this->omega/4.0/M_PI,1.0/3.0);
            }
            else
            {
                temp(k)=pow(fabs(msVector(k)*this->omega/this->b/M_PI),1.0/2.0);
            }
        }
        return temp;
    }

    template<int dim>
    int ClusterDynamicsParameters<dim>::getMobileSpeciesIndex(const double& m) const
    {/*!\returns the index of the mobile species of size m (-1 the single vacancy,
      * +1 the single interstitial), or -1 if there is none.
      */
        for(int k=0;k<mSize;k++)
        {
            if(fabs(msVector(k)-m)<FLT_EPSILON)
            {
                return k;
            }
        }
        return -1;
    }

    template<int dim>
    std::vector<std::tuple<int,int,double,int>> ClusterDynamicsParameters<dim>::getNucleationChannels() const
    {/*!\returns the reactions between mobile species a and b whose product is
      * not mobile: m_a+m_b larger than the largest mobile cluster (the product
      * joins the interstitial families) or smaller than the smallest one (the
      * vacancy families). The rate coefficient is that of the mobile equations,
      * K_ab=p_ab*4*pi*(r_a+r_b)*(D_a+D_b)/Omega, divided here by |m_a||m_b| so that
      * it multiplies the contents c_a c_b (content space). The last entry is 0 for a
      * vacancy product and 1 for an interstitial product.
      */
        std::vector<std::tuple<int,int,double,int>> temp;
        if(detD.empty() || iSize==0)
        {
            return temp;
        }
        const double mMax(msVector.maxCoeff());
        const double mMin(msVector.minCoeff());
        for(const auto& pair : reactionMap)
        {
            const int a(pair.first.first);
            const int c(pair.first.second);
            const double product(msVector(a)+msVector(c));
            if(pair.second>0.0 && (product>mMax+FLT_EPSILON || product<mMin-FLT_EPSILON))
            {
                temp.emplace_back(a,c,pair.second*4.0*M_PI*(mobileCaptureRadius(a)+mobileCaptureRadius(c))*(mobileMeanDiffusivity(a)+mobileMeanDiffusivity(c))/omega/std::fabs(msVector(a))/std::fabs(msVector(c)),product>0.0? 1 : 0);
            }
        }
        return temp;
    }

    template<int dim>
    Eigen::Array<double,1,ClusterDynamicsParameters<dim>::iSize/2> ClusterDynamicsParameters<dim>::getImmobileSpeciesBurgersMagnitude(const GrainContainerType& grains) const
    {
        Eigen::Array<double,1,iSize/2> temp(Eigen::Array<double,1,iSize/2>::Zero());
        const Eigen::Matrix<double,dim,dim> lat(grains.begin()->second->latticeBasis);
        const Eigen::Matrix<double,dim,iSize/2> localBurgers(lat*immobileSpeciesBurgers);
        
        for(size_t k=0; k<iSize/2; k++)
        {
            temp(k) = (localBurgers.col(k)).norm();
        }
        
        std::cout<< "immobileSpeciesBurgersMagnitude: "<<temp<<std::endl;
        return temp;
    }


    template<int dim>
    std::map<size_t,std::vector<Eigen::Matrix<double,dim,dim>>> ClusterDynamicsParameters<dim>::getD(const GrainContainerType& grains) const
    {
        std::map<size_t,std::vector<Eigen::Matrix<double,dim,dim>>> temp;
        
        for(const auto& grainPair : grains)
        {
            auto grainIter(temp.emplace(grainPair.first,std::vector<Eigen::Matrix<double,dim,dim>>()));
            
            for(size_t k=0; k<mSize; k++)
            {
                Eigen::Matrix<double,dim,dim> Dlocal;
                Dlocal<<   msD0(k,0)*exp(-msEm(k,0)/kB/T), msD0(k,1)*exp(-msEm(k,1)/kB/T), msD0(k,2)*exp(-msEm(k,2)/kB/T),
                msD0(k,1)*exp(-msEm(k,1)/kB/T), msD0(k,3)*exp(-msEm(k,3)/kB/T), msD0(k,4)*exp(-msEm(k,4)/kB/T),
                msD0(k,2)*exp(-msEm(k,2)/kB/T), msD0(k,4)*exp(-msEm(k,4)/kB/T), msD0(k,5)*exp(-msEm(k,5)/kB/T);
                
                // g^T inv(Dg)*g = (C2G*c)^T*inv(Dg)*C2G*c = c^T*C2G^T*inv(Dg)*C2G*c= c^T*inv(inv(C2G)*Dg*inv(C2G^T))*c
                // c^T*inv(C2G^T*Dg*C2G)*c
                // Dc=C2G^T*Dg*C2G
                // Dg = C2G*Dc*C2G^T
                const Eigen::Matrix<double,dim,dim> Dglobal(grainPair.second->C2G*Dlocal*grainPair.second->C2G.transpose());
                grainIter.first->second.emplace_back(Dglobal);
                
                if( fabs(Dglobal.eigenvalues().imag()(0)*Dglobal.eigenvalues().imag()(1)*Dglobal.eigenvalues().imag()(2)) > FLT_EPSILON)
                {
                    throw std::runtime_error("Diffusion tensor does not have proper eigen values.");
                }
                std::cout<<greenColor<<"  Diffusion coefficient for "<<k<<":\n"<<Dglobal<<std::endl;
            }
        }
        
        return temp;
    }


    template<int dim>
    std::vector<Eigen::Matrix<double,dim,dim>> ClusterDynamicsParameters<dim>::getDlocal() const
    {
        std::vector<Eigen::Matrix<double,dim,dim>> temp;
        
        for(size_t k=0; k<mSize; k++)
        {
            Eigen::Matrix<double,dim,dim> Dlocal;
            Dlocal<<   msD0(k,0)*exp(-msEm(k,0)/kB/T), msD0(k,1)*exp(-msEm(k,1)/kB/T), msD0(k,2)*exp(-msEm(k,2)/kB/T),
            msD0(k,1)*exp(-msEm(k,1)/kB/T), msD0(k,3)*exp(-msEm(k,3)/kB/T), msD0(k,4)*exp(-msEm(k,4)/kB/T),
            msD0(k,2)*exp(-msEm(k,2)/kB/T), msD0(k,4)*exp(-msEm(k,4)/kB/T), msD0(k,5)*exp(-msEm(k,5)/kB/T);
            
            temp.emplace_back(Dlocal);
        }
        
        return temp;
    }

    template<int dim>
    std::map<size_t,std::vector<Eigen::Matrix<double,dim,dim>>> ClusterDynamicsParameters<dim>::getInvD() const
    {
        std::map<size_t,std::vector<Eigen::Matrix<double,dim,dim>>> temp;
        for(const auto& pair : D)
        {
            auto mapIter(temp.emplace(pair.first,std::vector<Eigen::Matrix<double,dim,dim>>()));
            for(const auto& d : pair.second)
            {
                mapIter.first->second.push_back(d.inverse());
            }
        }
        return temp;
    }

    template<int dim>
    std::map<size_t,Eigen::Array<double,ClusterDynamicsParameters<dim>::mSize,1>> ClusterDynamicsParameters<dim>::getDetD() const
    {
        std::map<size_t,Eigen::Array<double,mSize,1>> temp;// (1,mSize);
        for(const auto& pair : D)
        {
            auto mapIter(temp.emplace(pair.first,Eigen::Array<double,mSize,1>()));
            for(int k=0;k<mSize;++k)
            {
                mapIter.first->second(k)=pair.second[k].determinant();
            }
        }
        return temp;
    }

    template<int dim>
    Eigen::Array<double,1,ClusterDynamicsParameters<dim>::iSize> ClusterDynamicsParameters<dim>::getInitialImmobileClusters(const Eigen::Array<double,1,iSize> initialImmobileClusters_SI, const double b_SI) const
    {
        const Eigen::Array<double,1,iSize/2> initialDensity = initialImmobileClusters_SI.template block<1,iSize/2>(0,0)*b_SI*b_SI*b_SI;
        const Eigen::Array<double,1,iSize/2> initialContent = initialImmobileClusters_SI.template block<1,iSize/2>(0,iSize/2);
        
        Eigen::Array<double,1,iSize> temp;
        temp<< initialDensity,initialContent;
        
        return temp;
    }

    template<int dim>
    Eigen::Array<double,1,ClusterDynamicsParameters<dim>::mSize> ClusterDynamicsParameters<dim>::equilibriumMobileConcentration(const double& stressTrace) const
    {

        return msVector.abs()*exp(-(msEf-stressTrace*msVector*msRelRelaxVol*omega/3.0)/kB/T);
    }

    template<int dim>
    Eigen::Array<double,1,ClusterDynamicsParameters<dim>::mSize> ClusterDynamicsParameters<dim>::boundaryMobileConcentration(const double& stressTrace,const double& normalTraction) const
    {
        return equilibriumMobileConcentration(stressTrace)*exp(-normalTraction*msVector*omega/kB/T);
    }

    template<int dim>
    Eigen::Array<double,1,ClusterDynamicsParameters<dim>::mSize> ClusterDynamicsParameters<dim>::dislocationMobileConcentration(const VectorDim& b,
                                                                                                                                const VectorDim& t,
                                                                                                                                const VectorDim& fPK,
                                                                                                                                const MatrixDim& stress) const
    {
        const VectorDim bxt(b.cross(t));
        const double bxtNorm2(bxt.squaredNorm());
        if(bxtNorm2>FLT_EPSILON)
        {
            const double fc(fPK.dot(bxt));
            return equilibriumMobileConcentration(stress.trace())*exp(-(msVector*omega*fc)/(kB*T*(bxtNorm2+0.05*b.squaredNorm())));
        }
        else
        {// screw direction
            return Eigen::Array<double,1,ClusterDynamicsParameters<dim>::mSize>::Zero();
        }
    }


    template<int dim>
    Eigen::Array<double,1,ClusterDynamicsParameters<dim>::iSize/2> ClusterDynamicsParameters<dim>::shapeSigmoid(const Eigen::Array<double,1,iSize/2>& m) const
    {/*!\returns S_k(m) (Eqs. 32, 33), the switch from bi-pyramid (0) to loop (1) of each
      * family for clusters of m defects. A step at the midpoint when the width is zero.
      */
        const Eigen::Array<double,1,iSize/2> m0 = (shapeTransitionMinSize+shapeTransitionMaxSize)*(0.5);
        const Eigen::Array<double,1,iSize/2> w = shapeTransitionWidth*(shapeTransitionMaxSize-shapeTransitionMinSize);
        return (w>0.0).select(1.0/(1.0+exp(-(m - m0)/w)),(m>=m0).select(Eigen::Array<double,1,iSize/2>::Ones(),Eigen::Array<double,1,iSize/2>::Zero()));
    }

    template<int dim>
    Eigen::Array<double,1,ClusterDynamicsParameters<dim>::iSize/2> ClusterDynamicsParameters<dim>::pyramidRadius(const Eigen::Array<double,1,iSize/2>& m) const
    {/*!\returns R_pyr (Eq. 31), the radius of a bi-pyramid of m defects: volume sqrt(8)*R^3 = m*Omega
      */
        return pow(m*omega/sqrt(8),1.0/3.0);
    }

    template<int dim>
    Eigen::Array<double,1,ClusterDynamicsParameters<dim>::iSize/2> ClusterDynamicsParameters<dim>::loopRadius(const Eigen::Array<double,1,iSize/2>& m) const
    {/*!\returns R_loop (Eq. 31), the radius of a planar loop of m defects: pi*R^2*b_k = m*Omega
      */
        return sqrt(m*omega/(M_PI*b*immobileSpeciesBurgersMagnitude));
    }

    template<int dim>
    Eigen::Array<double,1,ClusterDynamicsParameters<dim>::iSize/2> ClusterDynamicsParameters<dim>::effectiveRadiusDerivative(const Eigen::Array<double,1,iSize/2>& m) const
    {/*!\returns dR_k/dm, the derivative of the effective radius (Eq. 30) with respect to the defects per cluster,
      * including the change of shape: S*r_loop/(2m) + (1-S)*R_pyr/(3m) + dS/dm*(r_loop-R_pyr), with dS/dm = S*(1-S)/w.
      * The last term is the faster growth of the radius while the clusters turn from bi-pyramids into loops.
      */
        const Eigen::Array<double,1,iSize/2> S(shapeSigmoid(m));
        const Eigen::Array<double,1,iSize/2> w(shapeTransitionWidth*(shapeTransitionMaxSize-shapeTransitionMinSize));
        const Eigen::Array<double,1,iSize/2> dS((w>0.0).select(S*(1.0-S)/w,0.0)); // a step when the width is zero: no derivative
        const Eigen::Array<double,1,iSize/2> rl(loopRadius(m));
        const Eigen::Array<double,1,iSize/2> rp(pyramidRadius(m));
        return S*rl/(2.0*m)+(1.0-S)*rp/(3.0*m)+dS*(rl-rp);
    }

    template<int dim>
    Eigen::Array<double,1,ClusterDynamicsParameters<dim>::iSize/2> ClusterDynamicsParameters<dim>::shapeInterpolation(const Eigen::Array<double,1,iSize/2>& m, const Eigen::Array<double,1,iSize/2>& lowValue, const Eigen::Array<double,1,iSize/2>& highValue) const
    {/*!\returns (1-S)*lowValue + S*highValue for each family: lowValue for the bi-pyramid, highValue for the loop
      */
        const Eigen::Array<double,1,iSize/2> S(shapeSigmoid(m));
        return highValue*S + lowValue*(1.0 - S);
    }

    template<int dim>
    Eigen::Array<double,dim,dim> ClusterDynamicsParameters<dim>::shapeInterpolation(const Eigen::Array<double,1,iSize/2>& m, const Eigen::Array<double,dim,dim>& lowValue, const Eigen::Array<double,dim,dim>& highValue, const int& k) const
    {/*!\returns (1-S_k)*lowValue + S_k*highValue for the tensors of family k (Eq. 35)
      */
        const double S(shapeSigmoid(m)(k));
        return highValue*S + lowValue*(1.0 - S);
    }

    template<int dim>
    Eigen::Array<bool,1,ClusterDynamicsParameters<dim>::iSize/2> ClusterDynamicsParameters<dim>::populated(const Eigen::Array<double,1,iSize/2>& n, const Eigen::Array<double,1,iSize/2>& c) const
    {/*!@param[in] n clusters per atom, c defects per atom, of each family
      *\returns true for the families that hold clusters. The terms that involve the clusters are zero for the others.
      */
        return (n>emptyFamilyThreshold) && (c>emptyFamilyThreshold);
    }

    template<int dim>
    Eigen::Array<double,1,ClusterDynamicsParameters<dim>::iSize/2> ClusterDynamicsParameters<dim>::effectiveClusterRadius(const Eigen::Array<double,1,iSize/2>& n, const Eigen::Array<double,1,iSize/2>& c) const
    {/*!@param[in] n clusters per atom, c defects per atom, of each family
      *\returns R_k (Eq. 30), the effective radius of the clusters of each family, of m=c/n defects. Zero for an absent family.
      */
        const Eigen::Array<bool,1,iSize/2> pop(populated(n,c));
        const Eigen::Array<double,1,iSize/2> m(pop.select(c/n,1.0)); // defects per cluster (1 for an absent family, not used)
        return pop.select(shapeInterpolation(m,pyramidRadius(m),loopRadius(m)),0.0);
    }

    template<int dim>
    Eigen::Array<double,1,ClusterDynamicsParameters<dim>::iSize/2> ClusterDynamicsParameters<dim>::captureLength(const Eigen::Array<double,1,iSize/2>& n, const Eigen::Array<double,1,iSize/2>& c) const
    {/*!@param[in] n clusters per atom, c defects per atom, of each family
      *\returns L_k, the capture length of the clusters of each family, of m=c/n defects. Zero for an absent family.
      * The sink strength of one cluster for species m is Z_km*L_k. L_k blends the sink strength of each shape
      * with the shape sigmoid (as in MoDELib3), rather than using the blended radius in the loop formula:
      *
      *   L_k = (1-S_k)*4*pi*R^pyr_k + S_k*2*pi*R^loop_k
      *
      * 4*pi*R for the 3D bi-pyramid (sphere), 2*pi*R for the planar loop.
      */
        const Eigen::Array<bool,1,iSize/2> pop(populated(n,c));
        const Eigen::Array<double,1,iSize/2> m(pop.select(c/n,1.0)); // defects per cluster (1 for an absent family, not used)
        return pop.select(shapeInterpolation(m,4.0*M_PI*pyramidRadius(m),2.0*M_PI*loopRadius(m)),0.0);
    }

    template struct ClusterDynamicsParameters<3>;

    // Former sigmoidal interpolation and plotting helpers (not compiled), kept for reference

//  template<int dim>
//  Eigen::Array<double,1,ClusterDynamicsParameters<dim>::iSize/2> ClusterDynamicsParameters<dim>::sigmoidalVectorInterpolation(const Eigen::Array<double,1,iSize/2>& CI, const Eigen::Array<double,1,iSize/2>& N, const Eigen::Array<double,1,iSize/2>& lowValue, const Eigen::Array<double,1,iSize/2>& highValue) const
//  {
//      const Eigen::Array<double,1,iSize/2> n = CI/N/omega;
        
//      return highValue*sigmoid(n) + lowValue*(1.0 - sigmoid(n));
//  }

//  template<int dim>
//  Eigen::Array<double,dim,dim> ClusterDynamicsParameters<dim>::sigmoidalMatrixInterpolation(const Eigen::Array<double,1,iSize/2>& CI, const Eigen::Array<double,1,iSize/2>& N, const Eigen::Array<double,dim,dim>& lowValue, const Eigen::Array<double,dim,dim>& highValue, const int& index) const
//  {
//      const Eigen::Array<double,1,iSize/2> n = CI/N/omega;
        
//      return highValue*sigmoid(n)(index) + lowValue*(1.0 - sigmoid(n)(index));
//  }

//  template<int dim>
//  Eigen::Array<double,1,ClusterDynamicsParameters<dim>::iSize/2> ClusterDynamicsParameters<dim>::sigmoidalPlotVectorInterpolation(const Eigen::Array<double,1,iSize/2>& CI, const Eigen::Array<double,1,iSize/2>& N, const Eigen::Array<double,1,iSize/2>& lowValue, const Eigen::Array<double,1,iSize/2>& highValue) const
//  {
//      const Eigen::Array<double,1,iSize/2> n = CI/N/omega;
        
//      return lowValue*sigmoid(n) + highValue*(1.0 - sigmoid(n));
//  }

//  template<int dim>
//  Eigen::Array<double,1,ClusterDynamicsParameters<dim>::iSize/2> ClusterDynamicsParameters<dim>::clusterPlotRadius(const Eigen::Array<double,1,iSize/2>& CI, const Eigen::Array<double,1,iSize/2>& N) const
//  {
//      const Eigen::Array<double,1,iSize/2> n = CI/N/omega;
        
//      return sigmoidalPlotVectorInterpolation(CI,N,rpyr(n),rloop(n));
//  }


}
#endif
