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
    /* init */ omega(getAtomicVolume(ddBase)),
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
    /* init */ msCascadeFractions((mSize>1 && G0>0.0) ? TextFileParser(ddBase.poly.materialFile).readMatrix<double>("mobileSpeciesCascadeFractions",1,mSize,true) : Eigen::Matrix<double,1,mSize>::Ones()),
    /* init */ msSurvivingEfficiency(G0>0.0 ? TextFileParser(ddBase.poly.materialFile).readScalar<double>("mobileSpeciesSurvivingEfficiency",true) : 1.0),
    /* init */ G(G0*msSurvivingEfficiency*msCascadeFractions),
    /* init */ otherSinks(true? (TextFileParser(ddBase.poly.materialFile).readMatrix<double,1,mSize>("otherSinks_SI",true)*ddBase.poly.b_SI*ddBase.poly.b_SI).eval() : Eigen::Array<double,1,mSize>::Zero()),
    /* init */ dislocationSinks(iSize>0 ? (TextFileParser(ddBase.poly.materialFile).readMatrix<double,1,iSize/2>("dislocationSinks_SI",true)*ddBase.poly.b_SI*ddBase.poly.b_SI).array().eval() : Eigen::Array<double,1,iSize/2>::Zero().eval()), // read only with immobile species
    /* init */ initLoopSinks(iSize>0 ? getInitLoopSinks(TextFileParser(ddBase.poly.materialFile).readMatrix<double,1,iSize>("initLoopSinks_SI",true),ddBase.poly.b_SI) : Eigen::Array<double,1,iSize>::Zero().eval()), // read only with immobile species
    /* init */ reactionMap((true && mSize>1) ? getMap(TextFileParser(ddBase.poly.materialFile).readMatrix<double,mSize*(mSize+1)/2,3>("reactionPrefactorMap",true)) : std::map<std::pair<int,int>,double>()),
    /* init */ R1(mSize>1 ? getR1() : Eigen::Matrix<double,mSize,mSize>::Zero()),
    /* init */ R1cd(iSize>0 ? msVector.abs().matrix().asDiagonal()*R1*(1.0/msVector.abs()).matrix().asDiagonal() : Eigen::Matrix<double,mSize,mSize>::Zero().eval()),
    /* init */ R2(mSize>1 ? getR2() : std::vector<Eigen::Matrix<double,mSize,mSize>>(2,Eigen::Matrix<double,mSize,mSize>::Zero())),
    /* init */ discreteDislocationBias(true? TextFileParser(ddBase.poly.materialFile).readMatrix<double,2,mSize>("discreteDislocationBias",true).eval() : Eigen::Array<double,2,mSize>::Zero()),
    /* IMMOBILE SPECIES */
    /* init */ immobileSpeciesVector((true && iSize>0) ? TextFileParser(ddBase.poly.materialFile).readMatrix<int>("immobileSpeciesVector",1,iSize/2,true).array().template cast<double>() : Eigen::Array<double,1,iSize/2>::Zero().eval()),
    /* init */ immobileSpeciesRelRelaxVol((true && iSize>0) ? TextFileParser(ddBase.poly.materialFile).readMatrix<double,1,iSize/2>("immobileSpeciesRelRelaxVol",true).array() : Eigen::Array<double,1,iSize/2>::Zero().eval()),
    /* init */ immobileSpeciesBurgers((true && iSize>0) ? TextFileParser(ddBase.poly.materialFile).readMatrix<double,dim,iSize/2>("immobileSpeciesBurgers",true) : Eigen::Matrix<double,dim,iSize/2>::Zero()),
    /* init */ immobileSpeciesBurgersMagnitude((true && iSize>0) ? getImmobileSpeciesBurgersMagnitude(ddBase.poly.grains) : Eigen::Array<double,1,iSize/2>::Zero().eval()),
    /* init */ immobileBias((true && iSize>0) ? TextFileParser(ddBase.poly.materialFile).readMatrix<double,1,mSize>("immobileBias",true) : Eigen::Array<double,1,mSize>::Zero()),
    /* init */ a_bp((true && iSize>0) ? TextFileParser(ddBase.poly.materialFile).readScalar<double>("alpha_bp",true) : 0.0),
    /* init */ delVPyramid((true && iSize>0) ? TextFileParser(ddBase.poly.materialFile).readScalar<double>("delVPyramid",true)/pow(ddBase.poly.b_SI,3) : 0.0),
    /* init */ w0((true && iSize>0) ? TextFileParser(ddBase.poly.materialFile).readScalar<double>("w0",true) : 0.0),
    /* init */ n_s((true && iSize>0) ? TextFileParser(ddBase.poly.materialFile).readScalar<double>("n_s",true) : 0.0 ),
//    /* init */ evc((true && iSize>0) ? TextFileParser(ddBase.poly.materialFile).readScalar<double>("evc",true) : 0.0 ),
//    /* init */ Nvmax((true && iSize>0) ? TextFileParser(ddBase.poly.materialFile).readScalar<double>("Nvmax",true)*pow(ddBase.poly.b_SI,3) : 0.0 ),
    /* init */ nmin((true && iSize>0) ? TextFileParser(ddBase.poly.materialFile).readMatrix<double,1,iSize/2>("nmin",true) : Eigen::Array<double,1,iSize/2>::Zero()),
    /* init */ nmax((true && iSize>0) ? TextFileParser(ddBase.poly.materialFile).readMatrix<double,1,iSize/2>("nmax",true) : Eigen::Array<double,1,iSize/2>::Zero()),
    /* init */ n_min((true && iSize>0) ? TextFileParser(ddBase.poly.materialFile).readMatrix<double,1,iSize/2>("n_min",true).array() : Eigen::Array<double,1,iSize/2>::Zero().eval()),
    /* init */ loopCascadeFractions(getOptionalFamilyArray(ddBase.poly.materialFile,"loopCascadeFractions",0.0)),
    /* init */ nNuc(getOptionalFamilyArray(ddBase.poly.materialFile,"nNuc",1.0)),
    /* init */ loopG(G0*msSurvivingEfficiency*loopCascadeFractions),
    /* init */ useClusteringNucleation(iSize>0 && mSize>1 ? getOptionalInt(ddBase.poly.materialFile,"loopClusteringNucleation",0)>0 : false),
    /* init */ loopNucChannels(useClusteringNucleation ? getLoopNucChannels() : std::map<std::pair<int,int>,double>()),
    /* init */ clusteringShare(getClusteringShare()),
    /* init */ tauVac(getOptionalDouble(ddBase.poly.materialFile,"tau0_vLoop_SI",0.0)*(ddBase.poly.cs_SI/ddBase.poly.b_SI)
    /*      */        *exp(getOptionalDouble(ddBase.poly.materialFile,"Ea_vLoop_eV",0.0)*ddBase.poly.eV2J/ddBase.poly.mu_SI/pow(ddBase.poly.b_SI,3)/kB/T)),
    /* init */ dissolveEmbryosOnly(iSize>0 ? getOptionalInt(ddBase.poly.materialFile,"dissolveEmbryosOnly",0)>0 : false),
    /* init */ cLL(getOptionalFamilyArray(ddBase.poly.materialFile,"cLL",0.0)),
    /* init */ cLN(getOptionalFamilyArray(ddBase.poly.materialFile,"cLN",0.0)),
    /* init */ kappaLL(getOptionalDouble(ddBase.poly.materialFile,"kappaLL",0.0)),
    /* init */ kappaLN(getOptionalDouble(ddBase.poly.materialFile,"kappaLN",0.0)),
    /* init */ rhoNetwork(getOptionalDouble(ddBase.poly.materialFile,"rhoNetwork_SI",0.0)*ddBase.poly.b_SI*ddBase.poly.b_SI),
    /* init */ r_min(getOptionalFamilyArray(ddBase.poly.materialFile,"r_min",0.0)/ddBase.poly.b_SI),
    /* init */ loopSinkScale(getOptionalFamilyArray(ddBase.poly.materialFile,"loopSinkScale",1.0)),
    /* init */ clusteringNucleationInODE(getOptionalInt(ddBase.poly.materialFile,"loopClusteringNucleation",1)>0),
    /* init */ nucleationPerReaction(getOptionalInt(ddBase.poly.materialFile,"loopNucleationPerReaction",0)>0),
    /* init */ coalescenceClimbBurgers(getOptionalInt(ddBase.poly.materialFile,"coalescenceClimbBurgers",1)>0),
    /* init */ loopBias(getLoopBias()),
    /* init */ computeReactions((true && mSize>0 && mSize+iSize>1)? TextFileParser(ddBase.poly.materialFile).readScalar<int>("computeReactions",true) : 0),
//    /* init */ use0DsinkStrength((true && iSize>0) ? TextFileParser(ddBase.poly.materialFile).readScalar<int>("use0DsinkStrength",true) : 0),
//    /* init */ Zv((true && iSize>0) ? TextFileParser(ddBase.poly.materialFile).readMatrix<double,1,dim>("Zv",true) : Eigen::Array<double,1,dim>::Zero()),
//    /* init */ Zi((true && iSize>0) ? TextFileParser(ddBase.poly.materialFile).readMatrix<double,1,dim>("Zi",true) : Eigen::Array<double,1,dim>::Zero()),
//    /* init */ ZVec((true && iSize>0) ? TextFileParser(ddBase.poly.materialFile).readMatrix<double,1,mSize>("ZVec",true) : Eigen::Array<double,1,mSize>::Zero()),
//    /* init */ rc_il((true && iSize>0) ? TextFileParser(ddBase.poly.materialFile).readMatrix<double,1,dim>("rc_il",true)/ddBase.poly.b_SI : Eigen::Array<double,1,dim>::Zero().eval()),
    // /* init */ discreteDistanceFactor((true && iSize>0) ? TextFileParser(ddBase.poly.materialFile).readScalar<double>("distanceFactor",true) : 0.0),
    /* init */ discretizationFactor((true && iSize>0) ? TextFileParser(ddBase.poly.materialFile).readScalar<double>("clusterDiscretizationFactor",true) : 0.0),
    /* init */ discretizationTime((true && iSize>0) ? TextFileParser(ddBase.poly.materialFile).readScalar<double>("clusterDiscretizationTime",true)*(ddBase.poly.cs_SI/ddBase.poly.b_SI) : 0.0),
    /* init */ minimumLoopSize((true && iSize>0) ? TextFileParser(ddBase.poly.materialFile).readScalar<double>("minimumLoopSize",true)/ddBase.poly.b_SI : 0.0 )
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
        for(int k=0;k<iSize/2;++k)
        {// defects born in cascades as immobile clusters are part of the same balance
            if(loopCascadeFractions(k)<0.0)
            {
                throw std::runtime_error("ClusterDynamicsParameters: loopCascadeFractions must not be negative.");
            }
            if(loopCascadeFractions(k)>0.0 && nNuc(k)<1.0)
            {
                throw std::runtime_error("ClusterDynamicsParameters: nNuc must be at least 1 for each family with a non-zero loopCascadeFractions.");
            }
            if(immobileSpeciesVector(k)<0)
            {
                vacancySum+=loopCascadeFractions(k);
            }
            else if(immobileSpeciesVector(k)>0)
            {
                interstitialSum+=loopCascadeFractions(k);
            }
        }
        
        if(vacancySum>1.0+FLT_EPSILON)
        {
            std::cout<<redBoldColor<<"Warning: vacancy cluster fractions sum to more than one."<<defaultColor<<std::endl;
        }
        else if(vacancySum<1.0-FLT_EPSILON)
        {// with immobileIntegrator=cvode the remainder is born in immobile clusters. With euler it is not produced
            std::cout<<"vacancy cluster fractions sum to "<<vacancySum<<": the remainder "<<1.0-vacancySum<<" is not born as mobile or cascade-born clusters of loopCascadeFractions."<<std::endl;
        }
        if(interstitialSum>1.0+FLT_EPSILON)
        {
            std::cout<<redBoldColor<<"Warning: interstitial cluster fractions sum to more than one."<<defaultColor<<std::endl;
        }
        else if(interstitialSum<1.0-FLT_EPSILON)
        {// with immobileIntegrator=cvode the remainder is born in immobile clusters. With euler it is not produced
            std::cout<<"interstitial cluster fractions sum to "<<interstitialSum<<": the remainder "<<1.0-interstitialSum<<" is not born as mobile or cascade-born clusters of loopCascadeFractions."<<std::endl;
        }
        
        for(const auto& pair: reactionMap)
        {
            const auto key=pair.first;
            std::cout<<"find the interaction between "<<key.first<<" - "<<key.second<<" and parameter: "<<pair.second<<std::endl;
        }
        
        std::cout<<"first-order interaction matrix (in 1/s) is: "<<std::endl;
        std::cout<<R1*ddBase.poly.cs_SI/ddBase.poly.b_SI<<std::endl;
        
        for(size_t k=0; k<mSize; k++)
        {
            std::cout<<"second-order interaction matrix (in 1/s) for "<<static_cast<int>(msVector(k))<<"-species is: "<<std::endl;
            std::cout<<R2[k]*ddBase.poly.cs_SI/ddBase.poly.b_SI<<std::endl;
        }

        if(hasNucleation())
        {
            std::cout<<"immobile cluster nucleation: loopCascadeFractions="<<loopCascadeFractions<<", nNuc="<<nNuc<<std::endl;
            for(const auto& pair : loopNucChannels)
            {
                std::cout<<"  clustering channel "<<static_cast<int>(msVector(pair.first.first))<<" + "<<static_cast<int>(msVector(pair.first.second))<<", rate coefficient (in 1/s) "<<pair.second*ddBase.poly.cs_SI/ddBase.poly.b_SI<<std::endl;
            }
            if(loopNucChannels.size())
            {
                std::cout<<"  share of each family in the clustering nucleation: "<<clusteringShare<<std::endl;
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
    Eigen::Matrix<double,ClusterDynamicsParameters<dim>::mSize,ClusterDynamicsParameters<dim>::mSize> ClusterDynamicsParameters<dim>::getR1() const
    {
        int vIndex, iIndex;
        for(int k=0;k<mSize;k++)
        {
            if(fabs(msVector(k)+1)<FLT_EPSILON)
            {
                vIndex=k;
            }
            if(fabs(msVector(k)-1)<FLT_EPSILON)
            {
                iIndex=k;
            }
        }
        if(iIndex-vIndex!=1 && msVector.matrix().squaredNorm()>0)
        { // Expected sequence should look like ... -2,-1,+1,+2,+3 ...
            throw std::runtime_error("intersitital Index must be after vacancy Index.");
        }
        
        Eigen::Array<double,1,mSize> aveD(Eigen::Array<double,1,mSize>::Zero());
        for(int k=0;k<mSize;k++)
        {
            aveD(k)=pow(detD.begin()->second(k),1.0/3.0);
        }
        
        Eigen::Array<double,1,mSize> rn(Eigen::Array<double,1,mSize>::Zero());
        for(int k=0;k<mSize;k++)
        {
            if(fabs(msVector(k)-1)<FLT_EPSILON || fabs(msVector(k)+1)<FLT_EPSILON)
            {
                rn(k)=pow(3.0*this->omega/4.0/M_PI,1.0/3.0);
            }
            else
            {
                rn(k)=pow(fabs(msVector(k)*this->omega/this->b/M_PI),1.0/2.0);
            }
        }
                
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
                    // alpha(k)=it->second*2.0*M_PI*rn(k)*aveD(iIndex)/omega*exp(-Eb(k)/kB/T);
                    alpha(k)=it->second*8.0*M_PI*rn(iIndex)*aveD(iIndex)/omega*exp(-Eb(k)/kB/T);
                }
            }
        }
        
        Eigen::Matrix<double,mSize,mSize> tempR1(Eigen::Matrix<double,mSize,mSize>::Zero());
        for(int k=0;k<mSize;k++)
        {
            // Add other sinks
            tempR1(k,k) += (-otherSinks(k)*aveD(k));
            // Add dissociations
            tempR1(k,k) += (-alpha(k));
            if(k>0 && k<=vIndex)
            {
                tempR1(k,k-1) += alpha(k-1);
                tempR1(vIndex,k-1) += alpha(k-1);
            }
            else if(k>=iIndex && k<mSize-1)
            {
                tempR1(k,k+1) += alpha(k+1);
                tempR1(iIndex,k+1) += alpha(k+1);
            }
            
        }
        
        return tempR1;
    }

    template<int dim>
    std::vector<Eigen::Matrix<double,ClusterDynamicsParameters<dim>::mSize,ClusterDynamicsParameters<dim>::mSize>> ClusterDynamicsParameters<dim>::getR2() const
    {

        Eigen::Array<double,1,mSize> aveD(Eigen::Array<double,1,mSize>::Zero());
        for(int k=0;k<mSize;k++)
        {
            aveD(k)=pow(detD.begin()->second(k),1.0/3.0);
        }
        
        Eigen::Array<double,1,mSize> rn(Eigen::Array<double,1,mSize>::Zero());
        for(int k=0;k<mSize;k++)
        {
            if(fabs(msVector(k)-1)<FLT_EPSILON || fabs(msVector(k)+1)<FLT_EPSILON)
            {
                rn(k)=pow(3.0*this->omega/4.0/M_PI,1.0/3.0);
            }
            else
            {
                rn(k)=pow(fabs(msVector(k)*this->omega/this->b/M_PI),1.0/2.0);
            }
        }
        
        std::vector<Eigen::Matrix<double,mSize,mSize>> tempR2;
        Eigen::Matrix<double,mSize,mSize> R2sum(Eigen::Matrix<double,mSize,mSize>::Zero());
        double maxCoeff(0.0);
        
        for(size_t k=0;k<mSize;++k)
        {
            Eigen::Matrix<double,mSize,mSize> tempR2_k(Eigen::Matrix<double,mSize,mSize>::Zero());
            for(const auto& pair : reactionMap)
            {
                const auto key=pair.first;

                if(int(k)==key.first || int(k)==key.second)
                {// (key.first, key.second) consume k_specie
//                    if(msVector(key.second)>1.0 && use0DsinkStrength)
//                    { // 3D mobile + 1D immobile
//                        tempR2_k(key.first ,key.second) -= pair.second*2.0*M_PI*(rn(key.second))*(aveD(key.first))/omega;
//                        tempR2_k(key.second,key.first)  -= pair.second*2.0*M_PI*(rn(key.second))*(aveD(key.first))/omega;
//                    }
//                    else
//                    { // 3D mobile + 3D mobile
//                        tempR2_k(key.first ,key.second) -= pair.second*4.0*M_PI*(rn(key.first)+rn(key.second))*(aveD(key.first)+aveD(key.second))/omega;
//                        tempR2_k(key.second ,key.first) -= pair.second*4.0*M_PI*(rn(key.first)+rn(key.second))*(aveD(key.first)+aveD(key.second))/omega;
//                    }
                    tempR2_k(key.first ,key.second) -= pair.second*4.0*M_PI*(rn(key.first)+rn(key.second))*(aveD(key.first)+aveD(key.second))/omega;
                    tempR2_k(key.second ,key.first) -= pair.second*4.0*M_PI*(rn(key.first)+rn(key.second))*(aveD(key.first)+aveD(key.second))/omega;
                }
                else if( fabs(msVector(key.first)+msVector(key.second)-msVector(k))<FLT_EPSILON )
                {// (key.first, key.second) generate k_specie
                    const double same= (key.first==key.second) ? 0.5 : 1.0;
//                    if(msVector(key.second)>1.0 && use0DsinkStrength)
//                    { // 3D mobile + 1D immobile
//                        tempR2_k(key.first ,key.second) += same*pair.second*2.0*M_PI*(rn(key.second))*(aveD(key.first))/omega;
//                        tempR2_k(key.second ,key.first) += same*pair.second*2.0*M_PI*(rn(key.second))*(aveD(key.first))/omega;
//                    }
//                    else
//                    { // 3D mobile + 3D mobile
//                        tempR2_k(key.first ,key.second) += same*pair.second*4.0*M_PI*(rn(key.first)+rn(key.second))*(aveD(key.first)+aveD(key.second))/omega;
//                        tempR2_k(key.second ,key.first) += same*pair.second*4.0*M_PI*(rn(key.first)+rn(key.second))*(aveD(key.first)+aveD(key.second))/omega;
//                    }
                    tempR2_k(key.first ,key.second) += same*pair.second*4.0*M_PI*(rn(key.first)+rn(key.second))*(aveD(key.first)+aveD(key.second))/omega;
                    tempR2_k(key.second ,key.first) += same*pair.second*4.0*M_PI*(rn(key.first)+rn(key.second))*(aveD(key.first)+aveD(key.second))/omega;
                }
            }
            
            tempR2.push_back(tempR2_k);
            const double maxCoeffLocal(std::max(std::abs(tempR2.back().maxCoeff()),std::abs(tempR2.back().minCoeff())));
            
            if((tempR2.back()-tempR2.back().transpose()).norm()/maxCoeffLocal>FLT_EPSILON)
            {
                throw std::runtime_error("ClusterDynamicsParameters: R2_"+std::to_string(k)+" is not symmetric.");
            }
            maxCoeff=std::max(maxCoeff,maxCoeffLocal);
            R2sum+=tempR2.back()*msVector(k);
        }
        
        if(R2sum.norm()/maxCoeff>FLT_EPSILON)
        {
            //            throw std::runtime_error("Sum of R2's must be zero.");
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
    double ClusterDynamicsParameters<dim>::getAtomicVolume(const DislocationDynamicsBase<dim>& ddBase)
    {/*!\returns the atomic volume used by cluster dynamics, in units of b^3.
      * It is that of the lattice, unless the material file gives atomicVolume_SI:
      * a parameter set fitted with another value keeps its radii and rate coefficients.
      */
        const double atomicVolume_SI(getOptionalDouble(ddBase.poly.materialFile,"atomicVolume_SI",0.0));
        return atomicVolume_SI>0.0? atomicVolume_SI/pow(ddBase.poly.b_SI,3) : ddBase.poly.Omega;
    }

    template<int dim>
    Eigen::Array<double,ClusterDynamicsParameters<dim>::iSize/2,ClusterDynamicsParameters<dim>::mSize> ClusterDynamicsParameters<dim>::getLoopBias() const
    {/*!\returns the capture bias Z_km = Z0_km*ZDAD_k(p_m) of each immobile family k
      * for each mobile species m.
      * Z0 is the row of discreteDislocationBias for the sign of the family
      * (first row vacancy clusters, second row interstitial clusters).
      * ZDAD is the diffusion-anisotropy factor, with p_m=(D_c/D_a)^(1/6):
      * p_m for clusters in the basal plane, (p_m+p_m^-2)/2 for clusters in a
      * plane that contains the c axis. It is 1 for isotropic diffusion.
      * The fitted factor loopSinkScale of the family multiplies the result, so that
      * it acts on the absorption, on the loss of the mobile species and on the emission alike.
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
                    temp(k,m)=discreteDislocationBias(immobileSpeciesVector(k)<0.0? 0 : 1,m)*ZDAD*loopSinkScale(k);
                }
            }
        }
        return temp;
    }

    template<int dim>
    bool ClusterDynamicsParameters<dim>::hasNucleation() const
    {
        return (iSize>0) && ((loopG>0.0).any() || loopNucChannels.size()>0);
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
    std::map<std::pair<int,int>,double> ClusterDynamicsParameters<dim>::getLoopNucChannels() const
    {/*!\returns the reactions between mobile species that nucleate immobile clusters,
      * with the rate coefficient that getR2() assembles,
      *
      *      K_ab = p_ab*4*pi*(r_a+r_b)*(D_a+D_b)/Omega.
      *
      * A pair (a,b) of the reaction map is a nucleating channel when
      *  (i)  the two species have the same sign (a mixed pair is recombination), and
      *  (ii) no mobile species has the size of the product, so that the product
      *       leaves the mobile ladder and becomes an immobile cluster.
      * For the species {v,i,2i,3i} the channels are i+3i, 2i+2i and 2i+3i.
      * getR2() debits the reactants of these channels and credits no product,
      * which is the origin of the warning "Sum of R2 is not zero".
      */
        std::map<std::pair<int,int>,double> temp;
        if(detD.empty())
        {
            return temp;
        }

        Eigen::Array<double,1,mSize> aveD(Eigen::Array<double,1,mSize>::Zero());
        Eigen::Array<double,1,mSize> rn(Eigen::Array<double,1,mSize>::Zero());
        for(int k=0;k<mSize;k++)
        {// same radii and average diffusion coefficients as in getR2()
            aveD(k)=pow(detD.begin()->second(k),1.0/3.0);
            if(fabs(msVector(k)-1)<FLT_EPSILON || fabs(msVector(k)+1)<FLT_EPSILON)
            {
                rn(k)=pow(3.0*this->omega/4.0/M_PI,1.0/3.0);
            }
            else
            {
                rn(k)=pow(fabs(msVector(k)*this->omega/this->b/M_PI),1.0/2.0);
            }
        }

        for(const auto& pair : reactionMap)
        {
            const int a(pair.first.first);
            const int c(pair.first.second);
            if(pair.second<=0.0)
            {// channel switched off in the material file
                continue;
            }
            if(msVector(a)*msVector(c)<=0.0)
            {// opposite signs: not clustering
                continue;
            }
            bool productIsMobile(false);
            for(int k=0;k<mSize;k++)
            {
                if(fabs(msVector(a)+msVector(c)-msVector(k))<FLT_EPSILON)
                {
                    productIsMobile=true;
                    break;
                }
            }
            if(!productIsMobile)
            {
                temp.emplace(pair.first,pair.second*4.0*M_PI*(rn(a)+rn(c))*(aveD(a)+aveD(c))/omega);
            }
        }
        return temp;
    }

    template<int dim>
    Eigen::Array<double,1,ClusterDynamicsParameters<dim>::iSize/2> ClusterDynamicsParameters<dim>::getClusteringShare() const
    {/*!\returns the share of the clustering nucleation taken by each family.
      * The shares are those of loopCascadeFractions within the families of
      * the same sign, or equal shares when all those fractions are zero.
      */
        Eigen::Array<double,1,iSize/2> temp(Eigen::Array<double,1,iSize/2>::Zero());
        Eigen::Array<double,1,2> fracTot(Eigen::Array<double,1,2>::Zero());
        Eigen::Array<double,1,2> cntTot(Eigen::Array<double,1,2>::Zero());
        for(int k=0;k<iSize/2;++k)
        {
            const int p(immobileSpeciesVector(k)<0.0? 0 : 1);
            fracTot(p)+=loopCascadeFractions(k);
            cntTot(p)+=1.0;
        }
        for(int k=0;k<iSize/2;++k)
        {
            const int p(immobileSpeciesVector(k)<0.0? 0 : 1);
            temp(k)= fracTot(p)>0.0 ? loopCascadeFractions(k)/fracTot(p) : 1.0/cntTot(p);
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
    Eigen::Array<double,1,ClusterDynamicsParameters<dim>::iSize> ClusterDynamicsParameters<dim>::getInitLoopSinks(const Eigen::Array<double,1,iSize> initloopSinks_SI, const double b_SI) const
    {
        const Eigen::Array<double,1,iSize/2> initDen = initloopSinks_SI.template block<1,iSize/2>(0,0)*b_SI*b_SI*b_SI;
        // const Eigen::Array<double,1,iSize/2> initRad = initloopSinks_SI.template block<1,iSize/2>(0,iSize/2)/b_SI;
        const Eigen::Array<double,1,iSize/2> initRad = initloopSinks_SI.template block<1,iSize/2>(0,iSize/2);
        
        Eigen::Array<double,1,iSize> temp;
        temp<< initDen,initRad;
        
        return temp;
    }

    template<int dim>
    Eigen::Array<double,1,ClusterDynamicsParameters<dim>::mSize> ClusterDynamicsParameters<dim>::equilibriumMobileConcentration(const double& stressTrace) const
    {

        return msVector.abs()*exp(-(msEf-stressTrace*msVector*msRelRelaxVol*omega/3.0)/kB/T);
        // return msVector.abs()*exp(-msEf/kB/T);
    }

    template<int dim>
    Eigen::Array<double,1,ClusterDynamicsParameters<dim>::mSize> ClusterDynamicsParameters<dim>::boundaryMobileConcentration(const double& stressTrace,const double& normalTraction) const
    {
        return equilibriumMobileConcentration(stressTrace)*exp(-normalTraction*msVector*omega/kB/T);
        // return msVector.abs()*exp(-msEf/kB/T);
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


    // sigmoid function of number of vacancies FOR ALL SPECIES
    template<int dim>
    Eigen::Array<double,1,ClusterDynamicsParameters<dim>::iSize/2> ClusterDynamicsParameters<dim>::sigmoid(const Eigen::Array<double,1,iSize/2>& n) const
    {
        // parameters
        const Eigen::Array<double,1,iSize/2> n0 = ((nmin+nmax)*(0.5) - n_s);
        const Eigen::Array<double,1,iSize/2> w = w0*(nmax-nmin);
        
        return 1.0/(1.0+exp(-(n - n0)/w));
    }

    // pyramid radius function of V
    template<int dim>
    Eigen::Array<double,1,ClusterDynamicsParameters<dim>::iSize/2> ClusterDynamicsParameters<dim>::rpyr(const Eigen::Array<double,1,iSize/2>& n) const
    {
        return pow(n*omega/sqrt(8),1.0/3.0); // Vp = sqrt(8)*r^3;
    }

    // loop radius function of V
    template<int dim>
    Eigen::Array<double,1,ClusterDynamicsParameters<dim>::iSize/2> ClusterDynamicsParameters<dim>::rloop(const Eigen::Array<double,1,iSize/2>& n) const
    {
        return sqrt(n*omega/(M_PI*b*immobileSpeciesBurgersMagnitude)); // Vl = pi*b*r^2;
    }

    template<int dim>
    Eigen::Array<double,1,ClusterDynamicsParameters<dim>::iSize/2> ClusterDynamicsParameters<dim>::sigmoidalVectorInterpolation(const Eigen::Array<double,1,iSize/2>& CI, const Eigen::Array<double,1,iSize/2>& N, const Eigen::Array<double,1,iSize/2>& lowValue, const Eigen::Array<double,1,iSize/2>& highValue) const
    {
        const Eigen::Array<double,1,iSize/2> n = CI/N/omega;
        
        return highValue*sigmoid(n) + lowValue*(1.0 - sigmoid(n));
    }

    template<int dim>
    Eigen::Array<double,1,ClusterDynamicsParameters<dim>::iSize/2> ClusterDynamicsParameters<dim>::clusterRadius(const Eigen::Array<double,1,iSize/2>& CI, const Eigen::Array<double,1,iSize/2>& N) const
    {
        const Eigen::Array<double,1,iSize/2> n = CI/N/omega;
        
        return sigmoidalVectorInterpolation(CI,N,rpyr(n),rloop(n));
    }

    template<int dim>
    Eigen::Array<double,1,ClusterDynamicsParameters<dim>::iSize/2> ClusterDynamicsParameters<dim>::clusterDensity(const Eigen::Array<double,1,iSize/2>& CI, const Eigen::Array<double,1,iSize/2>& N) const
    {
        const Eigen::Array<double,1,iSize/2> n = CI/N/omega;
        const Eigen::Array<double,1,iSize/2> LoopS = 2.0*M_PI*rloop(n)*N;
        const Eigen::Array<double,1,iSize/2> PyrS = a_bp*4.0*M_PI*rpyr(n)*N;
        
        return sigmoidalVectorInterpolation(CI,N,PyrS,LoopS);
    }

    template<int dim>
    Eigen::Array<double,dim,dim> ClusterDynamicsParameters<dim>::sigmoidalMatrixInterpolation(const Eigen::Array<double,1,iSize/2>& CI, const Eigen::Array<double,1,iSize/2>& N, const Eigen::Array<double,dim,dim>& lowValue, const Eigen::Array<double,dim,dim>& highValue, const int& index) const
    {
        const Eigen::Array<double,1,iSize/2> n = CI/N/omega;
        
        return highValue*sigmoid(n)(index) + lowValue*(1.0 - sigmoid(n)(index));
    }

    template<int dim>
    Eigen::Array<double,1,ClusterDynamicsParameters<dim>::iSize/2> ClusterDynamicsParameters<dim>::sigmoidalPlotVectorInterpolation(const Eigen::Array<double,1,iSize/2>& CI, const Eigen::Array<double,1,iSize/2>& N, const Eigen::Array<double,1,iSize/2>& lowValue, const Eigen::Array<double,1,iSize/2>& highValue) const
    {
        const Eigen::Array<double,1,iSize/2> n = CI/N/omega;
        
        return lowValue*sigmoid(n) + highValue*(1.0 - sigmoid(n));
    }

    template<int dim>
    Eigen::Array<double,1,ClusterDynamicsParameters<dim>::iSize/2> ClusterDynamicsParameters<dim>::clusterPlotRadius(const Eigen::Array<double,1,iSize/2>& CI, const Eigen::Array<double,1,iSize/2>& N) const
    {
        const Eigen::Array<double,1,iSize/2> n = CI/N/omega;
        
        return sigmoidalPlotVectorInterpolation(CI,N,rpyr(n),rloop(n));
    }

    template struct ClusterDynamicsParameters<3>;

}
#endif
