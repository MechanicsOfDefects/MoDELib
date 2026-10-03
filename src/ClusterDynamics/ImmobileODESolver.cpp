/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ImmobileODESolver_cpp_
#define model_ImmobileODESolver_cpp_

#ifdef _OPENMP
#include <omp.h>
#endif

#include <stdexcept>
#include <string>
#include <ImmobileODESolver.h>

#ifdef _MODEL_SUNDIALS_
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#endif

namespace model
{

    template <int dim>
    ImmobileODESolver<dim>::ImmobileODESolver(const RateEquationsType& rateEquations_in,const double& relativeTolerance_in,const double& absoluteTolerance_in) :
    /* init */ rateEquations(rateEquations_in)
    /* init */,relativeTolerance(relativeTolerance_in)
    /* init */,absoluteTolerance(absoluteTolerance_in)
    /* init */,internalSteps(0)
    /* init */,maxInternalSteps(0)
    {

    }

    template <int dim>
    size_t ImmobileODESolver<dim>::lastInternalSteps() const
    {
        return internalSteps;
    }

    template <int dim>
    size_t ImmobileODESolver<dim>::lastMaxInternalSteps() const
    {
        return maxInternalSteps;
    }

#ifdef _MODEL_SUNDIALS_

    /*!\brief The data of one node passed to the right-hand side called by CVODE
     */
    template <int dim>
    struct ImmobileODENodeData
    {
        const ImmobileRateEquations<dim>* rateEquations;
        typename ImmobileRateEquations<dim>::ArrayM cM;
        typename ImmobileRateEquations<dim>::ArrayChi chi;
        double dt;
    };

    /*!\brief Right-hand side in the form required by CVODE.
     * The unknowns are [n_k, c_k], both per atom, and the time is scaled by the
     * step, s=t/dt in [0,1], so that all the numbers handled by the integrator
     * are of order one or less whatever the units.
     */
    template <int dim>
    int immobileODErhs(sunrealtype, N_Vector y, N_Vector ydot, void* userData)
    {
        typedef ImmobileRateEquations<dim> RateEquationsType;
        constexpr int nF(RateEquationsType::nF);
        const ImmobileODENodeData<dim>* const data(static_cast<const ImmobileODENodeData<dim>*>(userData));
        const sunrealtype* const yData(N_VGetArrayPointer(y));
        sunrealtype* const ydotData(N_VGetArrayPointer(ydot));
        typename RateEquationsType::ArrayF n,c,nDot,cDot;
        for(int k=0;k<nF;++k)
        {
            n(k)=yData[k];
            c(k)=yData[nF+k];
        }
        data->rateEquations->rates(data->cM,data->chi,n,c,nDot,cDot);
        for(int k=0;k<nF;++k)
        {
            ydotData[k]=nDot(k)*data->dt;
            ydotData[nF+k]=cDot(k)*data->dt;
        }
        return 0;
    }

    template <int dim>
    bool ImmobileODESolver<dim>::available()
    {
        return true;
    }

    template <int dim>
    void ImmobileODESolver<dim>::integrate(const Eigen::VectorXd& mobileDofs,Eigen::VectorXd& immobileDofs,const std::vector<ArrayChi>& chi,const double& dt)
    {
        internalSteps=0;
        maxInternalSteps=0;
        if(iSize==0 || dt<=0.0)
        {
            return;
        }
        const long int nNodes(immobileDofs.size()/iSize);
        if(mobileDofs.size()!=nNodes*mSize)
        {
            throw std::runtime_error("ImmobileODESolver: mobile and immobile fields have different numbers of nodes.");
        }
        if(chi.size()!=1 && long(chi.size())!=nNodes)
        {
            throw std::runtime_error("ImmobileODESolver: one set of nucleation fractions is needed for each node, or one for all nodes.");
        }
        const double omega(rateEquations.cdp.omega);

        long int failedNodes(0);
        long int firstFailedNode(-1);
        int firstFailedFlag(0);
        size_t stepsSum(0);
        size_t stepsMax(0);

#ifdef _OPENMP
#pragma omp parallel
#endif
        {// one integrator per thread, reinitialized at each node
            SUNContext context;
            SUNContext_Create(SUN_COMM_NULL,&context);
            N_Vector y(N_VNew_Serial(iSize,context));
            N_Vector absTol(N_VNew_Serial(iSize,context));
            N_VConst(absoluteTolerance,absTol);
            N_VConst(0.0,y);
            ImmobileODENodeData<dim> data;
            data.rateEquations=&rateEquations;
            data.cM.setZero();
            data.chi=chi[0];
            data.dt=dt;
            void* cvodeMemory(CVodeCreate(CV_BDF,context));
            CVodeInit(cvodeMemory,immobileODErhs<dim>,0.0,y);
            CVodeSVtolerances(cvodeMemory,relativeTolerance,absTol);
            SUNMatrix jacobian(SUNDenseMatrix(iSize,iSize,context));
            SUNLinearSolver linearSolver(SUNLinSol_Dense(y,jacobian,context));
            CVodeSetLinearSolver(cvodeMemory,linearSolver,jacobian); // Jacobian by finite differences
            CVodeSetUserData(cvodeMemory,&data);
            CVodeSetMaxNumSteps(cvodeMemory,200000);
            SUNContext_ClearErrHandlers(context); // failures are reported once, below

#ifdef _OPENMP
#pragma omp for schedule(dynamic,64) reduction(+:failedNodes,stepsSum) reduction(max:stepsMax)
#endif
            for(long int node=0;node<nNodes;++node)
            {
                sunrealtype* const yData(N_VGetArrayPointer(y));
                for(int m=0;m<mSize;++m)
                {
                    data.cM(m)=mobileDofs(node*mSize+m);
                }
                if(chi.size()>1)
                {
                    data.chi=chi[node];
                }
                for(int k=0;k<nF;++k)
                {
                    yData[k]=std::max(immobileDofs(node*iSize+k),0.0)*omega; // clusters per atom
                    yData[nF+k]=std::max(immobileDofs(node*iSize+nF+k),0.0);
                }
                CVodeReInit(cvodeMemory,0.0,y);
                sunrealtype tReached(0.0);
                const int flag(CVode(cvodeMemory,1.0,y,&tReached,CV_NORMAL));
                if(flag<0)
                {
                    failedNodes+=1;
#ifdef _OPENMP
#pragma omp critical (ImmobileODESolverFailure)
#endif
                    {
                        if(firstFailedNode<0 || node<firstFailedNode)
                        {
                            firstFailedNode=node;
                            firstFailedFlag=flag;
                        }
                    }
                }
                else
                {
                    for(int k=0;k<nF;++k)
                    {
                        immobileDofs(node*iSize+k)=std::max(yData[k],0.0)/omega;
                        immobileDofs(node*iSize+nF+k)=std::max(yData[nF+k],0.0);
                    }
                }
                long int steps(0);
                CVodeGetNumSteps(cvodeMemory,&steps);
                stepsSum+=steps;
                stepsMax=std::max(stepsMax,size_t(steps));
            }

            CVodeFree(&cvodeMemory);
            SUNLinSolFree(linearSolver);
            SUNMatDestroy(jacobian);
            N_VDestroy(y);
            N_VDestroy(absTol);
            SUNContext_Free(&context);
        }

        internalSteps=stepsSum;
        maxInternalSteps=stepsMax;
        if(failedNodes>0)
        {
            throw std::runtime_error("ImmobileODESolver: CVODE failed at "+std::to_string(failedNodes)+" nodes. First node "+std::to_string(firstFailedNode)+", CVODE flag "+std::to_string(firstFailedFlag)+".");
        }
    }

#else

    template <int dim>
    bool ImmobileODESolver<dim>::available()
    {
        return false;
    }

    template <int dim>
    void ImmobileODESolver<dim>::integrate(const Eigen::VectorXd&,Eigen::VectorXd&,const std::vector<ArrayChi>&,const double&)
    {
        throw std::runtime_error("ImmobileODESolver: MoDELib was built without SUNDIALS. Install SUNDIALS (CVODE) and configure again, or set immobileIntegrator=euler in ClusterDynamics.txt.");
    }

#endif

    template class ImmobileODESolver<3>;

}
#endif
