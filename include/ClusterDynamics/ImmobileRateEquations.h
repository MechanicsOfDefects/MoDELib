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
     * (Eqs. 36 and 41 of the report)
     *
     *   dn/dt = Nnuc - an - DeltaN
     *   dc/dt = PI + Gc + snuc - Ac - Aemit - DeltaC
     *
     * Each term is a member function that returns the vector over the families:
     *   Nnuc     nucleation rate, (snuc + PI)/n_nuc                     (Eq. 38)
     *   an       annealing (dissolution) of the clusters, n*nu_a       (Eq. 40)
     *   DeltaN   coalescence, like-loop and loop-network              (Eq. 37)
     *   PI       content of the clusters born in cascades             (Eq. 42)
     *   Gc       net absorption of mobile defects, from Gamma_km       (Eqs. 44-46)
     *   snuc     content of the clusters born from reactions between mobile species (Eq. 47)
     *   Ac       content removed by annealing, c*nu_a                  (Eq. 51)
     *   Aemit    thermal emission of vacancies from the vacancy clusters (Eq. 54)
     *   DeltaC   content removed by loop-network coalescence          (Eq. 43)
     * rates() is their sum, and the defects returned to the mobile species are
     * (Ac+Aemit).Iv and (Ac+Aemit).Ii (Eq. 55, ImmobileClusterSource).
     *
     * The class holds no finite-element type: the same functions give the rates
     * integrated by ImmobileSpeciesSolver, and the sink and source terms that the
     * immobile clusters contribute to the equations of the mobile species
     * (LoopSinkReaction, ImmobileClusterSource).
     *
     * Units: MoDELib units (lengths in b, times in b/cs). The state passed to the
     * functions is [n_k, c_k]: both per atom, as in the finite-element fields.
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
        typedef Eigen::Array<bool,1,nF> ArrayFB;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;

        const ClusterDynamicsParameters<dim>& cdp;
        // The constants below are computed in ClusterDynamicsParameters; short names for the equations
        const ArrayM Dbar; // orientation-averaged diffusion coefficient of each mobile species
        const ArrayM polarityM; // +1 interstitial, -1 vacancy
        const ArrayF polarityF;
        const ArrayFM polarityFM; // p_k*p_m: +1 when family k and mobile species m have the same sign
        const ArrayF Iv; // vacancy selector of the families (report I_v)
        const ArrayF Ii; // interstitial selector of the families (report I_i)
        const int vacancyIndex; // index of the single vacancy in the mobile species, or -1
        const int interstitialIndex; // index of the single interstitial in the mobile species, or -1
        const double epsilonV; // fraction of the cascade vacancies born in immobile clusters
        const double epsilonI; // fraction of the cascade interstitials born in immobile clusters
        const double G; // defect production rate per atom
        const std::vector<std::tuple<int,int,double,int>> nucleationChannels; // (species a, species b, rate coefficient, 0 vacancy/1 interstitial)

        /**********************************************************************/
        ImmobileRateEquations(const ClusterDynamicsParameters<dim>& cdp_in) :
        /* init */ cdp(cdp_in)
        /* init */,Dbar(cdp.mobileMeanDiffusivity)
        /* init */,polarityM(cdp.polarityM)
        /* init */,polarityF(cdp.polarityF)
        /* init */,polarityFM(cdp.polarityFM)
        /* init */,Iv(cdp.Iv)
        /* init */,Ii(cdp.Ii)
        /* init */,vacancyIndex(cdp.vacancyIndex)
        /* init */,interstitialIndex(cdp.interstitialIndex)
        /* init */,epsilonV(cdp.vacancyLoopCascadeFraction)
        /* init */,epsilonI(cdp.interstitialLoopCascadeFraction)
        /* init */,G(cdp.G0) // nominal dpa rate (Eq. 3 and 42)
        /* init */,nucleationChannels(cdp.nucleationChannels)
        {

        }

        /**********************************************************************/
        ArrayF shrinkGate(const ArrayF& R) const
        {/*!\returns the factor in [0,1] applied to the absorption of defects of
          * opposite sign, which makes a cluster shrink. It closes smoothly as the
          * radius approaches r_min, so that a family is not driven to negative
          * content. One when r_min=0. The same factor multiplies the loss of the
          * mobile species, so that the exchange conserves the defects.
          */
            const ArrayF x((R/cdp.r_min-1.0)/0.3);
            const ArrayF smoothStep((x<=0.0).select(ArrayF::Zero(),(x>=1.0).select(ArrayF::Ones(),x*x*(3.0-2.0*x))));
            return (cdp.r_min>0.0).select(smoothStep,1.0);
        }

        /**********************************************************************/
        ArrayFM absorptionCoefficients(const ArrayF& n,const ArrayF& c,const ArrayF& R) const
        {/*!@param[in] R the radii of the families, cdp.effectiveClusterRadius(n,c)
          *\returns the coefficients A_km such that Gamma_km=A_km*cM_m is the
          * number of defects of mobile species m absorbed by family k, per atom
          * and per unit time (Eq. 44):
          *
          *   A_km = gate_km * D_m * 2*pi*Z_km*R_k*n_k/Omega
          *
          * gate_km is shrinkGate for species of sign opposite to the family, and 1 otherwise.
          */
            const ArrayFM gate((polarityFM<0.0).select(shrinkGate(R).transpose().replicate(1,mSize),1.0));
            return gate*Dbar.replicate(nF,1)*2.0*M_PI*cdp.loopBias*R.transpose().replicate(1,mSize)*n.transpose().replicate(1,mSize)/cdp.omega;
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
//              const double m(std::max(cdp.nNuc(k),1.0));
                const double m(std::max(cdp.nucleatedClusterSize(k),1.0));
                const double m0(0.5*(cdp.shapeTransitionMinSize(k)+cdp.shapeTransitionMaxSize(k)));
                const double w(cdp.shapeTransitionWidth*(cdp.shapeTransitionMaxSize(k)-cdp.shapeTransitionMinSize(k)));
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
//      Eigen::Array<double,1,2> clusteringContent(const ArrayM& cM,const bool& events=false) const
        Eigen::Array<double,1,2> clusteringContent(const ArrayM& cM) const
        {/*!@param[in] cM the mobile concentrations, in defects per atom
          *\returns the content transferred from the mobile species to the immobile
          * vacancy (0) and interstitial (1) families by the reactions whose product
          * is not mobile, in defects per atom and per unit time.
          */
            Eigen::Array<double,1,2> temp(Eigen::Array<double,1,2>::Zero());
            for(const auto& channel : nucleationChannels)
            {
                const int a(std::get<0>(channel));
                const int b(std::get<1>(channel));
                const double Ca(std::max(cM(a),0.0)); // defects per atom (content space)
                const double Cb(std::max(cM(b),0.0));
                const double loss(std::get<2>(channel)*Ca*Cb); // clusters of each reactant lost per atom and unit time
//              temp(std::get<3>(channel))+= events ? loss
//              /*                               */ : a==b ? std::fabs(cdp.msVector(a))*loss
//              /*                               */ : (std::fabs(cdp.msVector(a))+std::fabs(cdp.msVector(b)))*loss;
                temp(std::get<3>(channel))+= a==b ? std::fabs(cdp.msVector(a))*loss
                /*                               */ : (std::fabs(cdp.msVector(a))+std::fabs(cdp.msVector(b)))*loss;
            }
            return temp;
        }

        /**********************************************************************/
        /* Terms of the cluster density equation (36) and of the content       */
        /* equation (41). Each returns the vector over the immobile families,  */
        /* per atom and per unit time. n,c>=0 (clusters and defects per atom). */

        /**********************************************************************/
        ArrayF PI(const ArrayChi& chi) const
        {/*!\returns P_I, the content of the clusters born in cascades (Eq. 42):
          * G*(epsilon_iL*chi_i + epsilon_vL*chi_v)
          */
            return G*(epsilonV*chi.row(0)*Iv+epsilonI*chi.row(1)*Ii);
        }

        /**********************************************************************/
        ArrayF snuc(const ArrayM& cM,const ArrayChi& chi) const
        {/*!\returns s_nuc, the content of the clusters born from the reactions
          * between mobile species whose product is not mobile (Eq. 47):
          * R_2i*chi_i + R_2v*chi_v
          */
            const Eigen::Array<double,1,2> R2(clusteringContent(cM.max(0.0))); // R_2v, R_2i (Eqs. 49, 50)
            return R2(0)*chi.row(0)*Iv+R2(1)*chi.row(1)*Ii;
        }

        /**********************************************************************/
        ArrayF Nnuc(const ArrayM& cM,const ArrayChi& chi) const
        {/*!\returns the nucleation rate (Eq. 38): (snuc + PI)/n_nuc.
          */
//             if(cdp.nucleationPerReaction)
//             { // REMOVE THIS IF STATEMENT 
//                 const Eigen::Array<double,1,2> events(clusteringContent(cM.max(0.0),true));
// //              return chi.row(0)*Iv*(G*epsilonV/cdp.nNuc.max(1.0)+events(0))+chi.row(1)*Ii*(G*epsilonI/cdp.nNuc.max(1.0)+events(1));
//                 return chi.row(0)*Iv*(G*epsilonV/cdp.nucleatedClusterSize.max(1.0)+events(0))+chi.row(1)*Ii*(G*epsilonI/cdp.nucleatedClusterSize.max(1.0)+events(1));
//             }
//          return (snuc(cM,chi)+PI(chi))/cdp.nNuc.max(1.0);
            return (snuc(cM,chi)+PI(chi))/cdp.nucleatedClusterSize.max(1.0);
        }

        /**********************************************************************/
        ArrayF an(const ArrayF& n,const ArrayF& c) const
        {/*!\returns a_n, the clusters lost by dissolution, n*nu_a (Eq. 40),
          * with nu_a=clusterDissolutionRate
          */
            return cdp.populated(n,c).select(n*cdp.clusterDissolutionRate,0.0);
        }

        /**********************************************************************/
        ArrayF Ac(const ArrayF& n,const ArrayF& c) const
        {/*!\returns A_c, the content lost by dissolution. The clusters that
          * dissolve are embryos, which hold nucleatedClusterSize defects, and the
          * loops that have grown keep their content: min(nucleatedClusterSize*n,c)*nu_a,
          * as in the 0D code that fitted the parameters (Eq. 51 of the report removes
          * the mean content c*nu_a instead). The content returns to the single mobile
          * vacancy or interstitial (Eq. 55).
          */
//          const ArrayF dissolved(cdp.dissolveEmbryosOnly? (cdp.nNuc.max(1.0)*n).min(c).eval() : c);
//          const ArrayF dissolved(cdp.dissolveEmbryosOnly? (cdp.nucleatedClusterSize.max(1.0)*n).min(c).eval() : c);
            const ArrayF dissolved((cdp.nucleatedClusterSize.max(1.0)*n).min(c));
            return cdp.populated(n,c).select(dissolved*cdp.clusterDissolutionRate,0.0);
        }

        /**********************************************************************/
        ArrayF Aemit(const ArrayF& n,const ArrayF& c,const ArrayF& R) const
        {/*!@param[in] R the radii of the families, cdp.effectiveClusterRadius(n,c)
          *\returns A_emit, the content lost by the vacancy clusters through thermal
          * emission of single vacancies (Eq. 54). Detailed balance with the
          * absorption of vacancies gives, for a vacancy family,
          *
          *   alpha_k*c_k = D_v * 2*pi*Z_kv*R_k*n_k/Omega * exp(-Eb_k/kB/T)
          *
          * with Eb_k the binding energy of a vacancy to the clusters of family k
          * (clusterBindingEnergy_eV). Zero for the interstitial families.
          */
            if(vacancyIndex<0)
            {
                return ArrayF::Zero();
            }
            const ArrayF emission(Dbar(vacancyIndex)*2.0*M_PI*cdp.loopBias.col(vacancyIndex).transpose()*R*n/cdp.omega*(-cdp.clusterBindingEnergy/cdp.kB/cdp.T).exp());
            return (cdp.populated(n,c) && (polarityF<0.0)).select(emission,0.0);
        }

        /**********************************************************************/
        ArrayF Gc(const ArrayM& cM,const ArrayF& n,const ArrayF& c,const ArrayFM& A) const
        {/*!@param[in] A the absorption coefficients, absorptionCoefficients(n,c,R)
          *\returns Gdot_c, the net content gained by absorption (Eqs. 44-46):
          * sum_m p_k p_m Gamma_km. Like defects add to the content, opposite
          * defects remove from it.
          */
            const ArrayFM Gamma(A*cM.max(0.0).replicate(nF,1)); // Eq. 44
            return cdp.populated(n,c).select((polarityFM*Gamma).rowwise().sum().transpose(),0.0);
        }

        /**********************************************************************/
        Eigen::Array<double,2,nF> coalescenceRates(const ArrayM& cM,const ArrayF& n,const ArrayF& c,const ArrayFM& A) const
        {/*!\returns nu_LL*phi_LL (row 0) and nu_LN*phi_LN (row 1) of each family
          * (Eqs. 58, 59 of the 0D model). The rate scale is the climb speed of
          * the loop from the absorbed like defects.
          */
            const ArrayFM Gamma(A*cM.max(0.0).replicate(nF,1)); // Eq. 44
            const ArrayF gain((polarityFM>0.0).select(Gamma,0.0).rowwise().sum().transpose()); // like defects absorbed
            const ArrayFB active(cdp.populated(n,c) && ((cdp.cLL>0.0) || (cdp.cLN>0.0)));
            const ArrayF m(active.select(c/n,1.0)); // defects per cluster (1 where not used)
            const ArrayF r(cdp.loopRadius(m));
            const ArrayF N(n/cdp.omega); // loops per unit volume
//          const ArrayF vAbs(0.5*r*gain/c*(cdp.coalescenceClimbBurgers? ArrayF::Ones() : cdp.immobileSpeciesBurgersMagnitude));
            const ArrayF vAbs(0.5*r*gain/c); // dr/dt of a loop of radius r (pi r^2 b_k = m Omega) growing by the like defects it absorbs
            const ArrayF rLL(r.min(N.pow(-1.0/3.0)));
            const ArrayF nuLL(cdp.cLL*vAbs*N.pow(1.0/3.0));
            const ArrayF phiLL(1.0-(-cdp.kappaLL*4.0/3.0*M_PI*rLL*rLL*rLL*N).exp());
            ArrayF nuLNphiLN(ArrayF::Zero());
//          if(cdp.rhoNetwork>0.0)
            if(cdp.dislocationSinkDensity>0.0)
            {
//              const ArrayF rLN(r.min(1.0/std::sqrt(cdp.rhoNetwork)));
//              nuLNphiLN=cdp.cLN*vAbs*std::sqrt(cdp.rhoNetwork)*(1.0-(-cdp.kappaLN*M_PI*rLN*rLN*cdp.rhoNetwork).exp());
                const ArrayF rLN(r.min(1.0/std::sqrt(cdp.dislocationSinkDensity)));
                nuLNphiLN=cdp.cLN*vAbs*std::sqrt(cdp.dislocationSinkDensity)*(1.0-(-cdp.kappaLN*M_PI*rLN*rLN*cdp.dislocationSinkDensity).exp());
            }
            Eigen::Array<double,2,nF> temp;
            temp.row(0)=active.select(nuLL*phiLL,0.0);
            temp.row(1)=active.select(nuLNphiLN,0.0);
            return temp;
        }

        /**********************************************************************/
        ArrayF DeltaN(const ArrayM& cM,const ArrayF& n,const ArrayF& c,const ArrayFM& A) const
        {/*!@param[in] A the absorption coefficients, absorptionCoefficients(n,c,R)
          *\returns DeltaN, the clusters lost by coalescence (Eq. 37):
          * (nu_LL*phi_LL + nu_LN*phi_LN)*n
          */
            const Eigen::Array<double,2,nF> nuPhi(coalescenceRates(cM,n,c,A));
            return (nuPhi.row(0)+nuPhi.row(1))*n;
        }

        /**********************************************************************/
        ArrayF DeltaC(const ArrayM& cM,const ArrayF& n,const ArrayF& c,const ArrayFM& A) const
        {/*!@param[in] A the absorption coefficients, absorptionCoefficients(n,c,R)
          *\returns DeltaC, the content lost by loop-network coalescence (Eq. 43):
          * nu_LN*phi_LN*c. Like-loop coalescence conserves the content.
          */
            return coalescenceRates(cM,n,c,A).row(1)*c;
        }

        /**********************************************************************/
        void rates(const ArrayM& cM,const ArrayChi& chi,const ArrayF& nIn,const ArrayF& cIn,ArrayF& nDot,ArrayF& cDot) const
        {/*!@param[in] cM the mobile concentrations, in defects per atom
          * @param[in] chi the nucleation fractions of the families
          * @param[in] nIn the clusters per atom of each family
          * @param[in] cIn the defects per atom stored in each family
          * @param[out] nDot,cDot their rates, Eqs. (36) and (41)
          */
            const ArrayF n(nIn.max(0.0));
            const ArrayF c(cIn.max(0.0));
            const ArrayF R(cdp.effectiveClusterRadius(n,c));
            const ArrayFM A(absorptionCoefficients(n,c,R)); // computed once, shared by Gc, DeltaN and DeltaC

            nDot=Nnuc(cM,chi)-an(n,c)-DeltaN(cM,n,c,A); // Eq. (36)
            cDot=PI(chi)+Gc(cM,n,c,A)+snuc(cM,chi)-Ac(n,c)-Aemit(n,c,R)-DeltaC(cM,n,c,A); // Eq. (41)
        }


        /**********************************************************************/
        /* Previous versions of the code above, commented out. Kept for        */
        /* reference; none of it is compiled.                                  */
        /**********************************************************************/

        // ---- class members ----
//      const double floorValue; // below this number per atom a family is treated as absent

        // ---- in ImmobileRateEquations() ----
//      /* init */,floorValue(1.0e-30)
//      /* init */,Dbar(getDbar())
//      /* init */,polarityM(cdp.msVector.sign())
//      /* init */,polarityF(cdp.immobileSpeciesVector.sign())
//      /* init */,polarityFM((polarityF.transpose().matrix()*polarityM.matrix()).array())
//      /* init */,Iv((polarityF<0.0).template cast<double>())
//      /* init */,Ii((polarityF>0.0).template cast<double>())
//      /* init */,vacancyIndex(getVacancyIndex())
//      /* init */,interstitialIndex(getInterstitialIndex())
//      /* init */,epsilonV(std::max(1.0-(cdp.msCascadeFractions.array()*(polarityM<0.0).template cast<double>()).sum(),0.0))
//      /* init */,epsilonI(std::max(1.0-(cdp.msCascadeFractions.array()*(polarityM>0.0).template cast<double>()).sum(),0.0))

        // ---- in ImmobileRateEquations() ----
//      /* init */,G(cdp.G0*cdp.msSurvivingEfficiency)

        // ---- in ImmobileRateEquations() ----
//      /* init */,nucleationChannels(cdp.clusteringNucleationInODE? getNucleationChannels() : std::vector<std::tuple<int,int,double,int>>())

        // ---- previous getDbar, getVacancyIndex, getInterstitialIndex, getNucleationChannels, populated() ----
//      /**********************************************************************/
//      ArrayM getDbar() const
//      {
//          ArrayM temp(ArrayM::Zero());
//          if(!cdp.detD.empty())
//          {
//              for(int m=0;m<mSize;++m)
//              {
//                  temp(m)=std::pow(std::max(cdp.detD.begin()->second(m),0.0),1.0/3.0);
//              }
//          }
//          return temp;
//      }

//      /**********************************************************************/
//      int getVacancyIndex() const
//      {
//          for(int m=0;m<mSize;++m)
//          {
//              if(std::fabs(cdp.msVector(m)+1.0)<FLT_EPSILON)
//              {
//                  return m;
//              }
//          }
//          return -1;
//      }

//      /**********************************************************************/
//      int getInterstitialIndex() const
//      {
//          for(int m=0;m<mSize;++m)
//          {
//              if(std::fabs(cdp.msVector(m)-1.0)<FLT_EPSILON)
//              {
//                  return m;
//              }
//          }
//          return -1;
//      }

//      /**********************************************************************/
//      std::vector<std::tuple<int,int,double,int>> getNucleationChannels() const
//      {/*!\returns the reactions between mobile species a and b whose product is
//        * not mobile: m_a+m_b larger than the largest mobile cluster (the product
//        * joins the interstitial families) or smaller than the smallest one (the
//        * vacancy families). The rate coefficient is that of the mobile equations,
//        * K_ab=p_ab*4*pi*(r_a+r_b)*(D_a+D_b)/Omega, divided here by |m_a||m_b| so that
//        * it multiplies the contents c_a c_b (content space).
//        */
//          std::vector<std::tuple<int,int,double,int>> temp;
//          if(cdp.detD.empty() || nF==0)
//          {
//              return temp;
//          }
//          const double mMax(cdp.msVector.maxCoeff());
//          const double mMin(cdp.msVector.minCoeff());
//          ArrayM rn(ArrayM::Zero());
//          for(int k=0;k<mSize;k++)
//          {// capture radii, as in ClusterDynamicsParameters::getR2
//              rn(k)= std::fabs(std::fabs(cdp.msVector(k))-1.0)<FLT_EPSILON ? std::pow(3.0*cdp.omega/4.0/M_PI,1.0/3.0)
//              /*                                                        */ : std::sqrt(std::fabs(cdp.msVector(k))*cdp.omega/cdp.b/M_PI);
//          }
//          for(const auto& pair : cdp.reactionMap)
//          {
//              const int a(pair.first.first);
//              const int b(pair.first.second);
//              const double product(cdp.msVector(a)+cdp.msVector(b));
//              if(pair.second>0.0 && (product>mMax+FLT_EPSILON || product<mMin-FLT_EPSILON))
//              {
//                  temp.emplace_back(a,b,pair.second*4.0*M_PI*(rn(a)+rn(b))*(Dbar(a)+Dbar(b))/cdp.omega,product>0.0? 1 : 0);
//                  temp.emplace_back(a,b,pair.second*4.0*M_PI*(rn(a)+rn(b))*(Dbar(a)+Dbar(b))/cdp.omega/std::fabs(cdp.msVector(a))/std::fabs(cdp.msVector(b)),product>0.0? 1 : 0); // content space
//              }
//          }
//          return temp;
//      }

//      /**********************************************************************/
//      bool populated(const double& n,const double& c) const
//      {
//          return n>floorValue && c>floorValue;
//      }

//      /**********************************************************************/
//      ArrayFB populated(const ArrayF& n,const ArrayF& c) const
//      {/*!\returns true for the families that hold clusters. The terms of the
//        * equations that involve the clusters are zero for the other families.
//        */
//          return (n>floorValue) && (c>floorValue);
//      }

        // ---- previous effectiveRadius, loopRadius() ----
//      ArrayF effectiveRadius(const ArrayF& n,const ArrayF& c) const
//      {/*!\returns the radius R_k of the clusters of each family: that of a
//        * bi-pyramid for small clusters and of a loop for large ones, with the
//        * sigmoidal interpolation of ClusterDynamicsParameters. Zero for an absent family.
//        */
//          ArrayF temp(ArrayF::Zero());
//          for(int k=0;k<nF;++k)
//          {
//              if(populated(n(k),c(k)))
//              {
//                  const double m(c(k)/n(k)); // defects per cluster
//                  const double m0(0.5*(cdp.nmin(k)+cdp.nmax(k))-cdp.n_s);
//                  const double w(cdp.w0*(cdp.nmax(k)-cdp.nmin(k)));
//                  const double S(w>0.0? 1.0/(1.0+std::exp(-(m-m0)/w)) : (m>=m0? 1.0 : 0.0));
//                  temp(k)=(1.0-S)*std::pow(m*cdp.omega/std::sqrt(8.0),1.0/3.0)+S*loopRadius(k,m);
//              }
//          }
//          return temp;
//      }

//      ArrayF effectiveRadius(const ArrayF& n,const ArrayF& c) const
//      {/*!\returns the radius R_k of the clusters of each family: that of a
//        * bi-pyramid for small clusters and of a loop for large ones, with the
//        * sigmoidal interpolation of ClusterDynamicsParameters. Zero for an absent family.
//        */
//          const ArrayFB pop(populated(n,c));
//          const ArrayF m(pop.select(c/n,1.0)); // defects per cluster (1 for an absent family, not used)
//          const ArrayF m0(0.5*(cdp.nmin+cdp.nmax)-cdp.n_s);
//          const ArrayF m0(0.5*(cdp.shapeTransitionMinSize+cdp.shapeTransitionMaxSize));
//          const ArrayF w(cdp.shapeTransitionWidth*(cdp.shapeTransitionMaxSize-cdp.shapeTransitionMinSize));
//          const ArrayF S((w>0.0).select(1.0/(1.0+(-(m-m0)/w).exp()),(m>=m0).select(ArrayF::Ones(),ArrayF::Zero())));
//          return pop.select((1.0-S)*(m*cdp.omega/std::sqrt(8.0)).pow(1.0/3.0)+S*loopRadius(m),0.0);
//      }

//      ArrayF effectiveRadius(const ArrayF& n,const ArrayF& c) const
//      {/*!\returns the effective radius R_k of the clusters of each family
//        * (ClusterDynamicsParameters::effectiveClusterRadius, Eq. 30). Zero for an absent family.
//        */
//          const ArrayFB pop(populated(n,c));
//          const ArrayF m(pop.select(c/n,1.0)); // defects per cluster (1 for an absent family, not used)
//          return pop.select(cdp.effectiveClusterRadius(m),0.0);
//      }

        // /**********************************************************************/
        // double loopRadius(const int& k,const double& m) const
        // {/*!\returns the radius of a planar loop of family k that stores m defects
        //   */
        //     return std::sqrt(m*cdp.omega/(M_PI*cdp.b*cdp.immobileSpeciesBurgersMagnitude(k)));
        // }

//      /**********************************************************************/
//      ArrayF loopRadius(const ArrayF& m) const
//      {/*!\returns the radius of a planar loop of each family that stores m defects
//        */
//          return (m*cdp.omega/(M_PI*cdp.b*cdp.immobileSpeciesBurgersMagnitude)).sqrt();
//      }

        // ---- previous shrinkGate() ----
//      ArrayF shrinkGate(const ArrayF& R) const
//      {/*!\returns the factor in [0,1] applied to the absorption of defects of
//        * opposite sign, which makes a cluster shrink. It closes smoothly as the
//        * radius approaches r_min, so that a family is not driven to negative
//        * content. One when r_min=0. The same factor multiplies the loss of the
//        * mobile species, so that the exchange conserves the defects.
//        */
//          ArrayF temp(ArrayF::Ones());
//          for(int k=0;k<nF;++k)
//          {
//              if(cdp.r_min(k)>0.0)
//              {
//                  const double x((R(k)/cdp.r_min(k)-1.0)/0.3);
//                  temp(k)= x<=0.0? 0.0 : (x>=1.0? 1.0 : x*x*(3.0-2.0*x));
//              }
//          }
//          return temp;
//      }

        // ---- previous absorptionCoefficients() ----
//      ArrayFM absorptionCoefficients(const ArrayF& n,const ArrayF& c) const
//      {/*!\returns the coefficients A_km such that Gamma_km=A_km*cM_m is the
//        * number of defects of mobile species m absorbed by family k, per atom
//        * and per unit time:
//        *
//        *   A_km = gate_km * D_m * 2*pi*Z_km*R_k*n_k/Omega
//        *
//        * gate_km is shrinkGate for species of sign opposite to the family, and 1 otherwise.
//        */
//          const ArrayF R(effectiveRadius(n,c));
//          const ArrayF gate(shrinkGate(R));
//          ArrayFM temp(ArrayFM::Zero());
//          for(int k=0;k<nF;++k)
//          {
//              if(R(k)>0.0)
//              {
//                  for(int m=0;m<mSize;++m)
//                  {
//                      const double g(polarityF(k)*polarityM(m)<0.0? gate(k) : 1.0);
//                      temp(k,m)=g*Dbar(m)*2.0*M_PI*cdp.loopBias(k,m)*R(k)*n(k)/cdp.omega;
//                  }
//              }
//          }
//          return temp;
//      }

        // ---- in absorptionCoefficients() ----
//          return absorptionCoefficients(n,c,effectiveRadius(n,c));

        // ---- previous emissionCoefficient, dissolvedContent, vacancySource, interstitialSource() ----
//      Replaced by the terms of the equations, below clusteringContent: emissionCoefficient by Aemit,
//      dissolvedContent by Ac (and an), vacancySource and interstitialSource by (Ac+Aemit).Iv and .Ii
//      in ImmobileClusterSource (Eq. 55).
//      /**********************************************************************/
//      double emissionCoefficient(const int& k,const double& n,const double& c) const
//      {/*!\returns the rate alpha_k*c_k at which the vacancy family k loses content
//        * by thermal emission of single vacancies, per atom and per unit time.
//        * Detailed balance with the absorption of vacancies gives
//        *
//        *   alpha_k*c_k = D_v * 2*pi*Z_kv*R_k*n_k/Omega * exp(-Eb/kB/T)
//        *
//        * with Eb the binding energy of a vacancy to the cluster
//        * (first entry of mobileSpeciesBindingEnergy_eV for the vacancy).
//        */
//          if(vacancyIndex<0 || polarityF(k)>0.0 || !populated(n,c))
//          {
//              return 0.0;
//          }
//          ArrayF nk(ArrayF::Zero());
//          ArrayF ck(ArrayF::Zero());
//          nk(k)=n;
//          ck(k)=c;
//          const double R(effectiveRadius(nk,ck)(k));
//          return Dbar(vacancyIndex)*2.0*M_PI*cdp.loopBias(k,vacancyIndex)*R*n/cdp.omega*std::exp(-cdp.Eb(vacancyIndex)/cdp.kB/cdp.T);
//      }

//      /**********************************************************************/
//      double dissolvedContent(const int& k,const double& n,const double& c) const
//      {/*!\returns the content that family k loses by the dissolution of its
//        * clusters, per atom and per unit time. Clusters dissolve at the rate
//        * nu*n, with nu=clusterDissolutionRate(k). By default each takes the mean content
//        * of the family, c/n, so that the content is lost at the rate nu*c. With
//        * dissolveEmbryosOnly the clusters that dissolve are embryos, which hold nNuc
//        * defects, and the loops that have grown keep their content: the rate is nu*nNuc*n.
//        * The content goes to the single vacancy (vacancySource) or the single
//        * interstitial (interstitialSource).
//        */
//          if(cdp.tauVac<=0.0 || polarityF(k)>0.0)
//          if(cdp.clusterDissolutionRate(k)<=0.0)
//          {
//              return 0.0;
//          }
//          return cdp.dissolveEmbryosOnly ? std::min(std::max(cdp.nNuc(k),1.0)*n,c)/cdp.tauVac : c/cdp.tauVac;
//          return cdp.dissolveEmbryosOnly ? std::min(std::max(cdp.nNuc(k),1.0)*n,c)*cdp.clusterDissolutionRate(k) : c*cdp.clusterDissolutionRate(k);
//      }

//      /**********************************************************************/
//      double vacancySource(const ArrayF& n,const ArrayF& c) const
//      {/*!\returns the vacancies returned to the mobile species by the vacancy
//        * families, per atom and per unit time: dissolution nu_k*c_k plus thermal
//        * emission. It is, term by term, what the content equation removes.
//        */
//          double temp(0.0);
//          for(int k=0;k<nF;++k)
//          {
//              if(polarityF(k)<0.0 && populated(n(k),c(k)))
//              {
//                  temp+=dissolvedContent(k,n(k),c(k));
//                  temp+=emissionCoefficient(k,n(k),c(k));
//              }
//          }
//          return temp;
//      }

//      /**********************************************************************/
//      double interstitialSource(const ArrayF& n,const ArrayF& c) const
//      {/*!\returns the interstitials returned to the mobile species by the
//        * interstitial families, per atom and per unit time: dissolution nu_k*c_k.
//        * It is, term by term, what the content equation removes.
//        */
//          double temp(0.0);
//          for(int k=0;k<nF;++k)
//          {
//              if(polarityF(k)>0.0 && populated(n(k),c(k)))
//              {
//                  temp+=dissolvedContent(k,n(k),c(k));
//              }
//          }
//          return temp;
//      }

        // ---- in nucleationWeights() ----
//              const double m0(0.5*(cdp.nmin(k)+cdp.nmax(k))-cdp.n_s);

        // ---- in clusteringContent() ----
//              const double Ca(std::max(cM(a),0.0)/std::fabs(cdp.msVector(a))); // clusters per atom

        // ---- in clusteringContent() ----
//              const double Cb(std::max(cM(b),0.0)/std::fabs(cdp.msVector(b)));

        // ---- in an() ----
//          return populated(n,c).select(n*cdp.clusterDissolutionRate,0.0);

        // ---- in Ac() ----
//          return populated(n,c).select(dissolved*cdp.clusterDissolutionRate,0.0);

        // ---- in Aemit() ----
//          return Aemit(n,c,effectiveRadius(n,c));

        // ---- in Aemit() ----
//          const ArrayF emission(Dbar(vacancyIndex)*2.0*M_PI*cdp.loopBias.col(vacancyIndex).transpose()*R*n/cdp.omega*std::exp(-cdp.Eb(vacancyIndex)/cdp.kB/cdp.T));

        // ---- in Aemit() ----
//          return (populated(n,c) && (polarityF<0.0)).select(emission,0.0);

        // ---- in Gc() ----
//          return populated(n,c).select((polarityFM*Gamma).rowwise().sum().transpose(),0.0);

        // ---- in coalescenceRates() ----
//          const ArrayFB active(populated(n,c) && ((cdp.cLL>0.0) || (cdp.cLN>0.0)));

        // ---- in coalescenceRates() ----
//          const ArrayF r(loopRadius(m));

        // ---- previous rates() ----
//      void rates(const ArrayM& cM,const ArrayChi& chi,const ArrayF& nIn,const ArrayF& cIn,ArrayF& nDot,ArrayF& cDot) const
//      {/*!@param[in] cM the mobile concentrations, in defects per atom
//        * @param[in] chi the nucleation fractions of the families
//        * @param[in] nIn the clusters per atom of each family
//        * @param[in] cIn the defects per atom stored in each family
//        * @param[out] nDot,cDot their rates
//        */
//          const ArrayF n(nIn.max(0.0));
//          const ArrayF c(cIn.max(0.0));
//          const ArrayM cMp(cM.max(0.0));

//          const ArrayFM A(absorptionCoefficients(n,c));
//          const Eigen::Array<double,1,2> sNuc(clusteringContent(cMp));
//          const Eigen::Array<double,1,2> eventsNuc(cdp.nucleationPerReaction? clusteringContent(cMp,true) : Eigen::Array<double,1,2>::Zero());

//          for(int k=0;k<nF;++k)
//          {
//              const int row(polarityF(k)<0.0? 0 : 1);

                // Nucleation: cascades and reactions between mobile species
//              const double source(chi(row,k)*(G*(row==0? epsilonV : epsilonI)+sNuc(row)));
//              nDot(k)=cdp.nucleationPerReaction? chi(row,k)*(G*(row==0? epsilonV : epsilonI)/std::max(cdp.nNuc(k),1.0)+eventsNuc(row))
//              /*                             */ : source/std::max(cdp.nNuc(k),1.0);
//              cDot(k)=source;

//              if(populated(n(k),c(k)))
//              {
                    // Absorption: like defects add to the content, opposite defects remove from it
//                  double gain(0.0);
//                  for(int m=0;m<mSize;++m)
//                  {
//                      const double Gamma(A(k,m)*cMp(m));
//                      cDot(k)+=polarityF(k)*polarityM(m)*Gamma;
//                      if(polarityF(k)*polarityM(m)>0.0)
//                      {
//                          gain+=Gamma;
//                      }
//                  }

//                  // Dissolution and thermal emission of the vacancy clusters
//                  if(row==0)
//                  {
//                      if(cdp.tauVac>0.0)
//                      {
//                          nDot(k)-=n(k)/cdp.tauVac;
//                      }
//                      cDot(k)-=dissolvedContent(k,n(k),c(k));
//                      cDot(k)-=emissionCoefficient(k,n(k),c(k));
//                  }
                    // Dissolution of the clusters (Eqs. 39, 40, 51). Zero for the families with no dissolution rate
//                  nDot(k)-=cdp.clusterDissolutionRate(k)*n(k);
//                  cDot(k)-=dissolvedContent(k,n(k),c(k));
                    // Thermal emission of the vacancy clusters
//                  if(row==0)
//                  {
//                      cDot(k)-=emissionCoefficient(k,n(k),c(k));
//                  }

                    // Coalescence. The rate scale is the climb speed of the loop from the absorbed like defects
//                  if(cdp.cLL(k)>0.0 || cdp.cLN(k)>0.0)
//                  {
//                      const double m(c(k)/n(k));
//                      const double r(loopRadius(k,m));
//                      const double N(n(k)/cdp.omega); // loops per unit volume
//                      const double vAbs(0.5*r*gain/c(k)*(cdp.coalescenceClimbBurgers? 1.0 : cdp.immobileSpeciesBurgersMagnitude(k)));
//                      const double rLL(std::min(r,std::pow(N,-1.0/3.0)));
//                      const double nuLL(cdp.cLL(k)*vAbs*std::pow(N,1.0/3.0));
//                      const double phiLL(1.0-std::exp(-cdp.kappaLL*4.0/3.0*M_PI*rLL*rLL*rLL*N));
//                      double nuLNphiLN(0.0);
//                      if(cdp.rhoNetwork>0.0)
//                      {
//                          const double rLN(std::min(r,1.0/std::sqrt(cdp.rhoNetwork)));
//                          nuLNphiLN=cdp.cLN(k)*vAbs*std::sqrt(cdp.rhoNetwork)*(1.0-std::exp(-cdp.kappaLN*M_PI*rLN*rLN*cdp.rhoNetwork));
//                      }
//                      nDot(k)-=(nuLL*phiLL+nuLNphiLN)*n(k);
//                      cDot(k)-=nuLNphiLN*c(k); // like-loop coalescence conserves the content
//                  }
//              }
//          }
//      }

        // ---- in rates() ----
//          const ArrayF R(effectiveRadius(n,c));

    };

}
#endif
