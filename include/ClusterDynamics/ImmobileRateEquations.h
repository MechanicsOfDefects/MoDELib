/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ImmobileRateEquations_H_
#define model_ImmobileRateEquations_H_

#include <cmath>
#include <cfloat>
#include <algorithm>
#include <vector>
#include <Eigen/Dense>
#include <ClusterDynamicsParameters.h>

namespace model
{

    /*!\brief The rate equations of the immobile clusters at one point of the crystal.
     *
     * Each immobile family k is described by its number of clusters per atom n_k
     * and by its defect content per atom c_k. With the mobile concentrations cM
     * held fixed, the two fields obey the ordinary differential equations
     *
     *   dn/dt = (s_nuc + G_L)/n_nuc - a_n - DeltaN
     *   dc/dt = G_c + G_L + s_nuc - A_c - A_emit - DeltaC
     *
     * where
     *   G_c      net absorption of mobile defects, from the rates Gamma_km
     *   G_L      content of the clusters born in cascades
     *   s_nuc    content of the clusters born from reactions between mobile species
     *   a_n, A_c dissolution of the vacancy clusters, with the lifetime tauVac
     *   A_emit   thermal emission of vacancies from the vacancy clusters
     *   DeltaN   coalescence, like-loop and loop-network
     *   DeltaC   content removed by loop-network coalescence
     *
     * The class holds no finite-element type: the same functions give the rates
     * integrated by ImmobileODESolver, and the sink and source terms that the
     * immobile clusters contribute to the equations of the mobile species
     * (LoopSinkReaction, VacancyLoopSource).
     *
     * Units: MoDELib units (lengths in b, times in b/cs). The state passed to the
     * functions is [n_k, c_k]: both per atom. The finite-element fields store the
     * number density per unit volume, N_k=n_k/Omega.
     */
    template <int dim>
    struct ImmobileRateEquations
    {
        static constexpr int mSize=ClusterDynamicsParameters<dim>::mSize;
        static constexpr int iSize=ClusterDynamicsParameters<dim>::iSize;
        static constexpr int nF=iSize/2;

        typedef Eigen::Array<double,1,nF> ArrayF;
        typedef Eigen::Array<double,1,mSize> ArrayM;
        typedef Eigen::Array<double,nF,mSize> ArrayFM;
        typedef Eigen::Array<double,2,nF> ArrayChi; // row 0: vacancy families, row 1: interstitial families
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;

        const ClusterDynamicsParameters<dim>& cdp;
        const double floorValue; // below this number per atom a family is treated as absent
        const ArrayM Dbar; // orientation-averaged diffusion coefficient of each mobile species
        const ArrayM polarityM; // +1 interstitial, -1 vacancy
        const ArrayF polarityF;
        const int vacancyIndex; // index of the single vacancy in the mobile species, or -1
        const double epsilonV; // fraction of the cascade vacancies born in immobile clusters
        const double epsilonI; // fraction of the cascade interstitials born in immobile clusters
        const double G; // defect production rate per atom
        const std::vector<std::tuple<int,int,double,int>> nucleationChannels; // (species a, species b, rate coefficient, 0 vacancy/1 interstitial)

        /**********************************************************************/
        ImmobileRateEquations(const ClusterDynamicsParameters<dim>& cdp_in) :
        /* init */ cdp(cdp_in)
        /* init */,floorValue(1.0e-30)
        /* init */,Dbar(getDbar())
        /* init */,polarityM(cdp.msVector.sign())
        /* init */,polarityF(cdp.immobileSpeciesVector.sign())
        /* init */,vacancyIndex(getVacancyIndex())
        /* init */,epsilonV(std::max(1.0-(cdp.msCascadeFractions.array()*(polarityM<0.0).template cast<double>()).sum(),0.0))
        /* init */,epsilonI(std::max(1.0-(cdp.msCascadeFractions.array()*(polarityM>0.0).template cast<double>()).sum(),0.0))
        /* init */,G(cdp.G0*cdp.msSurvivingEfficiency)
        /* init */,nucleationChannels(cdp.clusteringNucleationInODE? getNucleationChannels() : std::vector<std::tuple<int,int,double,int>>())
        {

        }

        /**********************************************************************/
        ArrayM getDbar() const
        {
            ArrayM temp(ArrayM::Zero());
            if(!cdp.detD.empty())
            {
                for(int m=0;m<mSize;++m)
                {
                    temp(m)=std::pow(std::max(cdp.detD.begin()->second(m),0.0),1.0/3.0);
                }
            }
            return temp;
        }

        /**********************************************************************/
        int getVacancyIndex() const
        {
            for(int m=0;m<mSize;++m)
            {
                if(std::fabs(cdp.msVector(m)+1.0)<FLT_EPSILON)
                {
                    return m;
                }
            }
            return -1;
        }

        /**********************************************************************/
        std::vector<std::tuple<int,int,double,int>> getNucleationChannels() const
        {/*!\returns the reactions between mobile species a and b whose product is
          * not mobile: m_a+m_b larger than the largest mobile cluster (the product
          * joins the interstitial families) or smaller than the smallest one (the
          * vacancy families). The rate coefficient is that of the mobile equations,
          * K_ab=p_ab*4*pi*(r_a+r_b)*(D_a+D_b)/Omega.
          */
            std::vector<std::tuple<int,int,double,int>> temp;
            if(cdp.detD.empty() || nF==0)
            {
                return temp;
            }
            const double mMax(cdp.msVector.maxCoeff());
            const double mMin(cdp.msVector.minCoeff());
            ArrayM rn(ArrayM::Zero());
            for(int k=0;k<mSize;k++)
            {// capture radii, as in ClusterDynamicsParameters::getR2
                rn(k)= std::fabs(std::fabs(cdp.msVector(k))-1.0)<FLT_EPSILON ? std::pow(3.0*cdp.omega/4.0/M_PI,1.0/3.0)
                /*                                                        */ : std::sqrt(std::fabs(cdp.msVector(k))*cdp.omega/cdp.b/M_PI);
            }
            for(const auto& pair : cdp.reactionMap)
            {
                const int a(pair.first.first);
                const int b(pair.first.second);
                const double product(cdp.msVector(a)+cdp.msVector(b));
                if(pair.second>0.0 && (product>mMax+FLT_EPSILON || product<mMin-FLT_EPSILON))
                {
                    temp.emplace_back(a,b,pair.second*4.0*M_PI*(rn(a)+rn(b))*(Dbar(a)+Dbar(b))/cdp.omega,product>0.0? 1 : 0);
                }
            }
            return temp;
        }

        /**********************************************************************/
        bool populated(const double& n,const double& c) const
        {
            return n>floorValue && c>floorValue;
        }

        /**********************************************************************/
        ArrayF effectiveRadius(const ArrayF& n,const ArrayF& c) const
        {/*!\returns the radius R_k of the clusters of each family: that of a
          * bi-pyramid for small clusters and of a loop for large ones, with the
          * sigmoidal interpolation of ClusterDynamicsParameters. Zero for an absent family.
          */
            ArrayF temp(ArrayF::Zero());
            for(int k=0;k<nF;++k)
            {
                if(populated(n(k),c(k)))
                {
                    const double m(c(k)/n(k)); // defects per cluster
                    const double m0(0.5*(cdp.nmin(k)+cdp.nmax(k))-cdp.n_s);
                    const double w(cdp.w0*(cdp.nmax(k)-cdp.nmin(k)));
                    const double S(w>0.0? 1.0/(1.0+std::exp(-(m-m0)/w)) : (m>=m0? 1.0 : 0.0));
                    temp(k)=(1.0-S)*std::pow(m*cdp.omega/std::sqrt(8.0),1.0/3.0)+S*loopRadius(k,m);
                }
            }
            return temp;
        }

        /**********************************************************************/
        double loopRadius(const int& k,const double& m) const
        {/*!\returns the radius of a planar loop of family k that stores m defects
          */
            return std::sqrt(m*cdp.omega/(M_PI*cdp.b*cdp.immobileSpeciesBurgersMagnitude(k)));
        }

        /**********************************************************************/
        ArrayF shrinkGate(const ArrayF& R) const
        {/*!\returns the factor in [0,1] applied to the absorption of defects of
          * opposite sign, which makes a cluster shrink. It closes smoothly as the
          * radius approaches r_min, so that a family is not driven to negative
          * content. One when r_min=0. The same factor multiplies the loss of the
          * mobile species, so that the exchange conserves the defects.
          */
            ArrayF temp(ArrayF::Ones());
            for(int k=0;k<nF;++k)
            {
                if(cdp.r_min(k)>0.0)
                {
                    const double x((R(k)/cdp.r_min(k)-1.0)/0.3);
                    temp(k)= x<=0.0? 0.0 : (x>=1.0? 1.0 : x*x*(3.0-2.0*x));
                }
            }
            return temp;
        }

        /**********************************************************************/
        ArrayFM absorptionCoefficients(const ArrayF& n,const ArrayF& c) const
        {/*!\returns the coefficients A_km such that Gamma_km=A_km*cM_m is the
          * number of defects of mobile species m absorbed by family k, per atom
          * and per unit time:
          *
          *   A_km = gate_km * D_m * 2*pi*Z_km*R_k*n_k/Omega
          *
          * gate_km is shrinkGate for species of sign opposite to the family, and 1 otherwise.
          */
            const ArrayF R(effectiveRadius(n,c));
            const ArrayF gate(shrinkGate(R));
            ArrayFM temp(ArrayFM::Zero());
            for(int k=0;k<nF;++k)
            {
                if(R(k)>0.0)
                {
                    for(int m=0;m<mSize;++m)
                    {
                        const double g(polarityF(k)*polarityM(m)<0.0? gate(k) : 1.0);
                        temp(k,m)=g*Dbar(m)*2.0*M_PI*cdp.loopBias(k,m)*R(k)*n(k)/cdp.omega;
                    }
                }
            }
            return temp;
        }

        /**********************************************************************/
        double emissionCoefficient(const int& k,const double& n,const double& c) const
        {/*!\returns the rate alpha_k*c_k at which the vacancy family k loses content
          * by thermal emission of single vacancies, per atom and per unit time.
          * Detailed balance with the absorption of vacancies gives
          *
          *   alpha_k*c_k = D_v * 2*pi*Z_kv*R_k*n_k/Omega * exp(-Eb/kB/T)
          *
          * with Eb the binding energy of a vacancy to the cluster
          * (first entry of mobileSpeciesBindingEnergy_eV for the vacancy).
          */
            if(vacancyIndex<0 || polarityF(k)>0.0 || !populated(n,c))
            {
                return 0.0;
            }
            ArrayF nk(ArrayF::Zero());
            ArrayF ck(ArrayF::Zero());
            nk(k)=n;
            ck(k)=c;
            const double R(effectiveRadius(nk,ck)(k));
            return Dbar(vacancyIndex)*2.0*M_PI*cdp.loopBias(k,vacancyIndex)*R*n/cdp.omega*std::exp(-cdp.Eb(vacancyIndex)/cdp.kB/cdp.T);
        }

        /**********************************************************************/
        double dissolvedContent(const int& k,const double& n,const double& c) const
        {/*!\returns the content that the vacancy family k loses by the dissolution
          * of its clusters, per atom and per unit time. Clusters dissolve at the rate
          * n/tau. By default each takes the mean content of the family, c/n, so that
          * the content is lost at the rate c/tau. With dissolveEmbryosOnly the clusters
          * that dissolve are embryos, which hold nNuc defects, and the loops that
          * have grown keep their content: the rate is nNuc*n/tau.
          */
            if(cdp.tauVac<=0.0 || polarityF(k)>0.0)
            {
                return 0.0;
            }
            return cdp.dissolveEmbryosOnly ? std::min(std::max(cdp.nNuc(k),1.0)*n,c)/cdp.tauVac : c/cdp.tauVac;
        }

        /**********************************************************************/
        double vacancySource(const ArrayF& n,const ArrayF& c) const
        {/*!\returns the vacancies returned to the mobile species by the vacancy
          * families, per atom and per unit time: dissolution c_k/tau plus thermal
          * emission. It is, term by term, what the content equation removes.
          */
            double temp(0.0);
            for(int k=0;k<nF;++k)
            {
                if(polarityF(k)<0.0 && populated(n(k),c(k)))
                {
                    temp+=dissolvedContent(k,n(k),c(k));
                    temp+=emissionCoefficient(k,n(k),c(k));
                }
            }
            return temp;
        }

        /**********************************************************************/
        ArrayChi nucleationWeights(const MatrixDim& stress,const Eigen::Matrix<double,dim,nF>& burgers) const
        {/*!@param[in] stress the local stress
          * @param[in] burgers the Burgers vectors of the families, in the global frame
          *\returns the fractions chi_v (row 0) and chi_i (row 1) of the nucleating
          * vacancy and interstitial content taken by each family. The weight of
          * a family is exp(stress:Omega_k/kB/T), with Omega_k the formation-volume
          * tensor of its nucleus of n_nuc defects, normalized over the families of
          * the same sign. Without stress the families of a sign have equal fractions.
          */
            ArrayChi temp(ArrayChi::Zero());
            Eigen::Array<double,1,2> sum(Eigen::Array<double,1,2>::Zero());
            for(int k=0;k<nF;++k)
            {
                const double m(std::max(cdp.nNuc(k),1.0));
                const double m0(0.5*(cdp.nmin(k)+cdp.nmax(k))-cdp.n_s);
                const double w(cdp.w0*(cdp.nmax(k)-cdp.nmin(k)));
                const double S(w>0.0? 1.0/(1.0+std::exp(-(m-m0)/w)) : (m>=m0? 1.0 : 0.0));
                const double b2(burgers.col(k).squaredNorm());
                const MatrixDim loopTensor(b2>0.0? (burgers.col(k)*burgers.col(k).transpose()/b2).eval() : MatrixDim::Zero());
                const MatrixDim formationVolume((S*loopTensor+(1.0-S)*MatrixDim::Identity()/3.0)*m*cdp.omega);
                const int row(polarityF(k)<0.0? 0 : 1);
                temp(row,k)=std::exp((stress.array()*formationVolume.array()).sum()/cdp.kB/cdp.T);
                sum(row)+=temp(row,k);
            }
            for(int k=0;k<nF;++k)
            {
                for(int row=0;row<2;++row)
                {
                    if(sum(row)>0.0)
                    {
                        temp(row,k)/=sum(row);
                    }
                }
            }
            return temp;
        }

        /**********************************************************************/
        ArrayChi nucleationWeights() const
        {// without stress
            return nucleationWeights(MatrixDim::Zero(),Eigen::Matrix<double,dim,nF>::Zero());
        }

        /**********************************************************************/
        Eigen::Array<double,1,2> clusteringContent(const ArrayM& cM,const bool& events=false) const
        {/*!@param[in] cM the mobile concentrations, in defects per atom
          * @param[in] events return the number of reactions instead of their content
          *\returns the content transferred from the mobile species to the immobile
          * vacancy (0) and interstitial (1) families by the reactions whose product
          * is not mobile, in defects per atom and per unit time; with events, the
          * number of these reactions per atom and per unit time.
          */
            Eigen::Array<double,1,2> temp(Eigen::Array<double,1,2>::Zero());
            for(const auto& channel : nucleationChannels)
            {
                const int a(std::get<0>(channel));
                const int b(std::get<1>(channel));
                const double Ca(std::max(cM(a),0.0)/std::fabs(cdp.msVector(a))); // clusters per atom
                const double Cb(std::max(cM(b),0.0)/std::fabs(cdp.msVector(b)));
                const double loss(std::get<2>(channel)*Ca*Cb); // clusters of each reactant lost per atom and unit time
                temp(std::get<3>(channel))+= events ? loss
                /*                               */ : a==b ? std::fabs(cdp.msVector(a))*loss
                /*                               */ : (std::fabs(cdp.msVector(a))+std::fabs(cdp.msVector(b)))*loss;
            }
            return temp;
        }

        /**********************************************************************/
        void rates(const ArrayM& cM,const ArrayChi& chi,const ArrayF& nIn,const ArrayF& cIn,ArrayF& nDot,ArrayF& cDot) const
        {/*!@param[in] cM the mobile concentrations, in defects per atom
          * @param[in] chi the nucleation fractions of the families
          * @param[in] nIn the clusters per atom of each family
          * @param[in] cIn the defects per atom stored in each family
          * @param[out] nDot,cDot their rates
          */
            const ArrayF n(nIn.max(0.0));
            const ArrayF c(cIn.max(0.0));
            const ArrayM cMp(cM.max(0.0));

            const ArrayFM A(absorptionCoefficients(n,c));
            const Eigen::Array<double,1,2> sNuc(clusteringContent(cMp));
            const Eigen::Array<double,1,2> eventsNuc(cdp.nucleationPerReaction? clusteringContent(cMp,true) : Eigen::Array<double,1,2>::Zero());

            for(int k=0;k<nF;++k)
            {
                const int row(polarityF(k)<0.0? 0 : 1);

                // Nucleation: cascades and reactions between mobile species
                const double source(chi(row,k)*(G*(row==0? epsilonV : epsilonI)+sNuc(row)));
                nDot(k)=cdp.nucleationPerReaction? chi(row,k)*(G*(row==0? epsilonV : epsilonI)/std::max(cdp.nNuc(k),1.0)+eventsNuc(row))
                /*                             */ : source/std::max(cdp.nNuc(k),1.0);
                cDot(k)=source;

                if(populated(n(k),c(k)))
                {
                    // Absorption: like defects add to the content, opposite defects remove from it
                    double gain(0.0);
                    for(int m=0;m<mSize;++m)
                    {
                        const double Gamma(A(k,m)*cMp(m));
                        cDot(k)+=polarityF(k)*polarityM(m)*Gamma;
                        if(polarityF(k)*polarityM(m)>0.0)
                        {
                            gain+=Gamma;
                        }
                    }

                    // Dissolution and thermal emission of the vacancy clusters
                    if(row==0)
                    {
                        if(cdp.tauVac>0.0)
                        {
                            nDot(k)-=n(k)/cdp.tauVac;
                        }
                        cDot(k)-=dissolvedContent(k,n(k),c(k));
                        cDot(k)-=emissionCoefficient(k,n(k),c(k));
                    }

                    // Coalescence. The rate scale is the climb speed of the loop from the absorbed like defects
                    if(cdp.cLL(k)>0.0 || cdp.cLN(k)>0.0)
                    {
                        const double m(c(k)/n(k));
                        const double r(loopRadius(k,m));
                        const double N(n(k)/cdp.omega); // loops per unit volume
                        const double vAbs(0.5*r*gain/c(k)*(cdp.coalescenceClimbBurgers? 1.0 : cdp.immobileSpeciesBurgersMagnitude(k)));
                        const double rLL(std::min(r,std::pow(N,-1.0/3.0)));
                        const double nuLL(cdp.cLL(k)*vAbs*std::pow(N,1.0/3.0));
                        const double phiLL(1.0-std::exp(-cdp.kappaLL*4.0/3.0*M_PI*rLL*rLL*rLL*N));
                        double nuLNphiLN(0.0);
                        if(cdp.rhoNetwork>0.0)
                        {
                            const double rLN(std::min(r,1.0/std::sqrt(cdp.rhoNetwork)));
                            nuLNphiLN=cdp.cLN(k)*vAbs*std::sqrt(cdp.rhoNetwork)*(1.0-std::exp(-cdp.kappaLN*M_PI*rLN*rLN*cdp.rhoNetwork));
                        }
                        nDot(k)-=(nuLL*phiLL+nuLNphiLN)*n(k);
                        cDot(k)-=nuLNphiLN*c(k); // like-loop coalescence conserves the content
                    }
                }
            }
        }

    };

}
#endif
