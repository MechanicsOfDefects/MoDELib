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

        /**********************************************************************/
        ImmobileRateEquations(const ClusterDynamicsParameters<dim>& cdp_in) :
        /* init */ cdp(cdp_in)
        {

        }

        /**********************************************************************/
//      ArrayF shrinkGate(const ArrayF& R) const
//      {/*!\returns the factor in [0,1] applied to the absorption of defects of
//        * opposite sign, which makes a cluster shrink. It closes smoothly as the
//        * radius approaches r_min, so that a family is not driven to negative
//        * content. One when r_min=0. The same factor multiplies the loss of the
//        * mobile species, so that the exchange conserves the defects.
//        */
//          const ArrayF x((R/cdp.r_min-1.0)/0.3);
//          const ArrayF smoothStep((x<=0.0).select(ArrayF::Zero(),(x>=1.0).select(ArrayF::Ones(),x*x*(3.0-2.0*x))));
//          return (cdp.r_min>0.0).select(smoothStep,1.0);
//      }

        ArrayF shrinkGate(const ArrayF& m) const
        {/*!@param[in] m the mean number of defects per cluster of each family
          *\returns the factor in [0,1] applied to the absorption of defects of
          * opposite sign, which makes a cluster shrink. It closes smoothly as m
          * approaches the smallest immobile cluster of the family (cdp.minimumImmobileClusterSize),
          * and is fully open at twice that size, so that a family is not driven to negative
          * content. The same factor multiplies the loss of the mobile species, so that
          * the exchange conserves the defects.
          */
            const ArrayF x(m/cdp.minimumImmobileClusterSize-1.0);
            return (x<=0.0).select(ArrayF::Zero(),(x>=1.0).select(ArrayF::Ones(),x*x*(3.0-2.0*x)));
        }

        /**********************************************************************/
        ArrayFM absorptionCoefficients(const ArrayF& n,const ArrayF& c,const ArrayF& L) const
        {/*!@param[in] L the capture lengths of the families, cdp.captureLength(n,c)
          *\returns the coefficients A_km such that Gamma_km=A_km*cM_m is the
          * number of defects of mobile species m absorbed by family k, per atom
          * and per unit time (Eq. 44, with 2*pi*R_k replaced by the capture length L_k):
          *
          *   A_km = gate_km * D_m * Z_km*L_k*n_k/Omega
          *
          * gate_km is shrinkGate for species of sign opposite to the family, and 1 otherwise.
          */
            const ArrayFM gate((cdp.polarityFM<0.0).select(shrinkGate(cdp.populated(n,c).select(c/n,1.0)).transpose().replicate(1,mSize),1.0));
//          return gate*cdp.mobileMeanDiffusivity.replicate(nF,1)*2.0*M_PI*cdp.loopBias*R.transpose().replicate(1,mSize)*n.transpose().replicate(1,mSize)/cdp.omega; // loop sink 2*pi*R with the blended radius
            return gate*cdp.mobileMeanDiffusivity.replicate(nF,1)*cdp.loopBias*L.transpose().replicate(1,mSize)*n.transpose().replicate(1,mSize)/cdp.omega;
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
            const ArrayF S(cdp.shapeSigmoid(cdp.nucleatedClusterSize.max(1.0))); // shape of the nucleus of each family (Eqs. 32, 33)
            for(int k=0;k<nF;++k)
            {
                const double m(std::max(cdp.nucleatedClusterSize(k),1.0));
                const double b2(burgers.col(k).squaredNorm());
                const MatrixDim loopTensor(b2>0.0? (burgers.col(k)*burgers.col(k).transpose()/b2).eval() : MatrixDim::Zero());
                const MatrixDim formationVolume((S(k)*loopTensor+(1.0-S(k))*MatrixDim::Identity()/3.0)*m*cdp.omega);
                const int row(cdp.polarityF(k)<0.0? 0 : 1);
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
        Eigen::Array<double,1,2> clusteringContent(const ArrayM& cM) const
        {/*!@param[in] cM the mobile concentrations, in defects per atom
          *\returns the content transferred from the mobile species to the immobile
          * vacancy (0) and interstitial (1) families by the reactions whose product
          * is not mobile, in defects per atom and per unit time.
          */
            Eigen::Array<double,1,2> temp(Eigen::Array<double,1,2>::Zero());
            for(const auto& channel : cdp.nucleationChannels)
            {
                const int a(std::get<0>(channel));
                const int b(std::get<1>(channel));
                const double Ca(std::max(cM(a),0.0)); // defects per atom (content space)
                const double Cb(std::max(cM(b),0.0));
                const double loss(std::get<2>(channel)*Ca*Cb); // clusters of each reactant lost per atom and unit time
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
            return cdp.G0*(cdp.vacancyLoopCascadeFraction*chi.row(0)*cdp.Iv+cdp.interstitialLoopCascadeFraction*chi.row(1)*cdp.Ii);
        }

        /**********************************************************************/
        ArrayF snuc(const ArrayM& cM,const ArrayChi& chi) const
        {/*!\returns s_nuc, the content of the clusters born from the reactions
          * between mobile species whose product is not mobile (Eq. 47):
          * R_2i*chi_i + R_2v*chi_v
          */
            const Eigen::Array<double,1,2> R2(clusteringContent(cM.max(0.0))); // R_2v, R_2i (Eqs. 49, 50)
            return R2(0)*chi.row(0)*cdp.Iv+R2(1)*chi.row(1)*cdp.Ii;
        }

        /**********************************************************************/
        ArrayF Nnuc(const ArrayM& cM,const ArrayChi& chi) const
        {/*!\returns the nucleation rate (Eq. 38): (snuc + PI)/n_nuc.
          */
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
            const ArrayF dissolved((cdp.nucleatedClusterSize.max(1.0)*n).min(c));
            return cdp.populated(n,c).select(dissolved*cdp.clusterDissolutionRate,0.0);
        }

        /**********************************************************************/
        ArrayF Aemit(const ArrayF& n,const ArrayF& c,const ArrayF& L) const
        {/*!@param[in] L the capture lengths of the families, cdp.captureLength(n,c)
          *\returns A_emit, the content lost by the vacancy clusters through thermal
          * emission of single vacancies (Eq. 54). Detailed balance with the
          * absorption of vacancies gives, for a vacancy family,
          *
          *   alpha_k*c_k = D_v * Z_kv*L_k*n_k/Omega * exp(-Eb_k/kB/T)
          *
          * with Eb_k the binding energy of a vacancy to the clusters of family k
          * (clusterBindingEnergy_eV). Zero for the interstitial families.
          */
            if(cdp.vacancyIndex<0)
            {
                return ArrayF::Zero();
            }
//          const ArrayF emission(cdp.mobileMeanDiffusivity(cdp.vacancyIndex)*2.0*M_PI*cdp.loopBias.col(cdp.vacancyIndex).transpose()*R*n/cdp.omega*(-cdp.clusterBindingEnergy/cdp.kB/cdp.T).exp()); // loop sink 2*pi*R with the blended radius
            const ArrayF emission(cdp.mobileMeanDiffusivity(cdp.vacancyIndex)*cdp.loopBias.col(cdp.vacancyIndex).transpose()*L*n/cdp.omega*(-cdp.clusterBindingEnergy/cdp.kB/cdp.T).exp());
            return (cdp.populated(n,c) && (cdp.polarityF<0.0)).select(emission,0.0);
        }

        /**********************************************************************/
        ArrayF Gc(const ArrayM& cM,const ArrayF& n,const ArrayF& c,const ArrayFM& A) const
        {/*!@param[in] A the absorption coefficients, absorptionCoefficients(n,c,R)
          *\returns Gdot_c, the net content gained by absorption (Eqs. 44-46):
          * sum_m p_k p_m Gamma_km. Like defects add to the content, opposite
          * defects remove from it.
          */
            const ArrayFM Gamma(A*cM.max(0.0).replicate(nF,1)); // Eq. 44
            return cdp.populated(n,c).select((cdp.polarityFM*Gamma).rowwise().sum().transpose(),0.0);
        }

        /**********************************************************************/
        Eigen::Array<double,2,nF> coalescenceRates(const ArrayM& cM,const ArrayF& n,const ArrayF& c,const ArrayFM& A) const
        {/*!\returns nu_LL*phi_LL (row 0) and nu_LN*phi_LN (row 1) of each family
          * (Eqs. 58, 59 of the 0D model). The rate scale is the climb speed of
          * the loop from its net growth: like defects absorbed minus opposite defects absorbed minus thermal emission.
          */
            const ArrayFM Gamma(A*cM.max(0.0).replicate(nF,1)); // Eq. 44
//          const ArrayF gain((cdp.polarityFM>0.0).select(Gamma,0.0).rowwise().sum().transpose()); // like defects absorbed
            // Coalescence happens when the loops grow into each other, so its climb speed is that of the NET growth:
            // like defects absorbed - opposite defects absorbed - thermal emission (vacancy families), in content.
            // A family that shrinks (net<0) does not coalesce by this channel. The report and the 0D model use the
            // gain of like defects only, which keeps coalescence active at saturation but is not a growth speed:
            // a loop absorbing 101 v and 100 i would move 50 times faster than one absorbing 2 v and 1 i.
//          const ArrayF net(((cdp.polarityFM*Gamma).rowwise().sum().transpose()-Aemit(n,c,cdp.effectiveClusterRadius(n,c))).max(0.0)); // net content gained by growth, per atom and unit time
            const ArrayF net(((cdp.polarityFM*Gamma).rowwise().sum().transpose()-Aemit(n,c,cdp.captureLength(n,c))).max(0.0)); // net content gained by growth, per atom and unit time
            const ArrayFB active(cdp.populated(n,c) && ((cdp.cLL>0.0) || (cdp.cLN>0.0)));
            const ArrayF m(active.select(c/n,1.0)); // defects per cluster (1 where not used)
            const ArrayF r(cdp.shapeInterpolation(m,cdp.pyramidRadius(m),cdp.loopRadius(m))); // effective radius R_k (Eq. 30), as in the absorption
            const ArrayF N(n/cdp.omega); // loops per unit volume
//          const ArrayF vAbs(0.5*r*gain/c); // loop-form growth speed dr/dt = r*gain/(2c), with r the effective radius (as in MoDELib3); the exact blended derivative would add a dS/dm term
//          const ArrayF vAbs(0.5*r*net/c); // loop-form growth speed dr/dt = r*net/(2c), with r the effective radius; the exact blended derivative would add a dS/dm term
            const ArrayF vAbs(cdp.effectiveRadiusDerivative(m)*m*net/c); // dR/dt = dR/dm * dm/dt, with dm/dt = m*net/c the net growth per cluster and dR/dm that of the effective radius, change of shape included
            const ArrayF rLL(r.min(N.pow(-1.0/3.0)));
            const ArrayF nuLL(cdp.cLL*vAbs*N.pow(1.0/3.0));
            const ArrayF phiLL(1.0-(-cdp.kappaLL*4.0/3.0*M_PI*rLL*rLL*rLL*N).exp());
            ArrayF nuLNphiLN(ArrayF::Zero());
            if(cdp.dislocationSinkDensity>0.0)
            {
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
//          const ArrayF R(cdp.effectiveClusterRadius(n,c));
//          const ArrayFM A(absorptionCoefficients(n,c,R)); // computed once, shared by Gc, DeltaN and DeltaC
            const ArrayF L(cdp.captureLength(n,c)); // capture lengths, shared by the absorption and the emission
            const ArrayFM A(absorptionCoefficients(n,c,L)); // computed once, shared by Gc, DeltaN and DeltaC

            nDot=Nnuc(cM,chi)-an(n,c)-DeltaN(cM,n,c,A); // Eq. (36)
//          cDot=PI(chi)+Gc(cM,n,c,A)+snuc(cM,chi)-Ac(n,c)-Aemit(n,c,R)-DeltaC(cM,n,c,A); // Eq. (41)
            cDot=PI(chi)+Gc(cM,n,c,A)+snuc(cM,chi)-Ac(n,c)-Aemit(n,c,L)-DeltaC(cM,n,c,A); // Eq. (41)
        }


    };

}
#endif
