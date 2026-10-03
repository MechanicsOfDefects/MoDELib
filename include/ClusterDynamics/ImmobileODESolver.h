/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ImmobileODESolver_H_
#define model_ImmobileODESolver_H_

#include <vector>
#include <Eigen/Dense>
#include <ImmobileRateEquations.h>

namespace model
{

    /*!\brief Integrates the rate equations of the immobile clusters over a time
     * step, at every node of the finite-element mesh.
     *
     * With the mobile concentrations held fixed during the step, the equations
     * of a node involve the state of that node only (ImmobileRateEquations).
     * The mesh therefore gives as many independent systems of iSize equations
     * as it has nodes. They are stiff: the lifetime of the vacancy clusters and
     * the coalescence times are much shorter than a dose step. Each system is
     * integrated implicitly with CVODE (variable-order BDF, dense Newton),
     * and the nodes are distributed over the OpenMP threads.
     *
     * The solver works on the vectors of nodal values held by the trial
     * functions of ClusterDynamicsFEM: it reads the mobile values and
     * overwrites the immobile ones. Nothing goes through a file.
     *
     * The class is compiled without CVODE when SUNDIALS is not found
     * (see available()); integrate() then throws.
     */
    template <int dim>
    class ImmobileODESolver
    {

    public:

        typedef ImmobileRateEquations<dim> RateEquationsType;
        typedef typename RateEquationsType::ArrayChi ArrayChi;
        static constexpr int mSize=RateEquationsType::mSize;
        static constexpr int iSize=RateEquationsType::iSize;
        static constexpr int nF=RateEquationsType::nF;

        const RateEquationsType& rateEquations;
        const double relativeTolerance;
        const double absoluteTolerance; // on the numbers and contents per atom

        ImmobileODESolver(const RateEquationsType& rateEquations_in,const double& relativeTolerance_in,const double& absoluteTolerance_in);

        /*!\returns true if MoDELib was built with SUNDIALS
         */
        static bool available();

        /*!@param[in] mobileDofs nodal values of the mobile species, mSize per node, in defects per atom
         * @param[in,out] immobileDofs nodal values of the immobile families, iSize per node:
         *   number densities per unit volume, then defect contents per atom
         * @param[in] chi nucleation fractions of the families: one entry per node, or a single entry for all nodes
         * @param[in] dt the time step
         */
        void integrate(const Eigen::VectorXd& mobileDofs,Eigen::VectorXd& immobileDofs,const std::vector<ArrayChi>& chi,const double& dt);

        size_t lastInternalSteps() const; // sum over the nodes
        size_t lastMaxInternalSteps() const; // largest number at one node

    private:

        size_t internalSteps;
        size_t maxInternalSteps;

    };

}
#endif
