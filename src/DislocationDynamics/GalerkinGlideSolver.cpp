/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GalerkinGlideSolver_cpp_
#define model_GalerkinGlideSolver_cpp_

#include <GalerkinGlideSolver.h>
#include <PeriodicVelocityConstraints.h>
#include <TerminalColors.h>

#include <DislocationNetwork.h>

namespace model
{
    
    template <int dim>
    GalerkinGlideSolver<dim>::GalerkinGlideSolver(const DislocationNetwork<dim>& DN_in) :
    /* init */ DislocationGlideSolverBase<dim>(DN_in)
    {
        std::cout<<greenBoldColor<<"Creating GalerkinGlideSolver"<<defaultColor<<std::endl;
    }
        
    template <int dim>
    Eigen::VectorXd GalerkinGlideSolver<dim>::getNodeVelocities() const
    {
        return lumpedSolve();
    }
    
    
    template <int dim>
    size_t GalerkinGlideSolver<dim>::assembleNCtriplets(TripletContainerType& kqqT, Eigen::VectorXd& Fq) const
    {
        const size_t Ndof(this->DN.networkNodes().size()*NdofXnode); // the total number of dof in the subnetwork
        // Rebuild after every topology change; never retain indices on nodes.
        std::map<size_t,size_t> nodeIndices;
        for(const auto& entry : this->DN.networkNodes())
        {
            nodeIndices.emplace(entry.second.lock()->sID,nodeIndices.size());
        }
        kqqT.clear();
        Fq.setZero(Ndof);
        for (const auto& networkLink : this->DN.networkLinks())
        {
            if (!networkLink.second.lock()->isBoundarySegment())
            {
                const auto link=networkLink.second.lock();
                link->addToGlobalAssembly(kqqT,Fq,nodeIndices.at(link->source->sID),nodeIndices.at(link->sink->sID)); // loop over each segment and add segment contributions to kqqT and Fq
            }
        }
        return Ndof;
    }
    
// Shared periodic crossings enforce every participating loop interpolation.
    template <int dim>
    size_t GalerkinGlideSolver<dim>::assemblePeriodicConstraints(TripletContainerType &zT) const
    {
        return periodicVelocityConstraints<dim>(this->DN,zT);
    }

    template <int dim>
    Eigen::VectorXd GalerkinGlideSolver<dim>::lumpedSolve() const
    {
        std::cout<<", glideSolver "<<std::flush;
        TripletContainerType kqqT;
        Eigen::VectorXd Fq;
        const size_t Ndof=assembleNCtriplets(kqqT,Fq);

        // Physical boundary nodes keep independent velocities. Only actual
        // same-grain periodic subdivisions interpolate their source endpoints.
        // With no subdivisions Z is the identity, including nonperiodic domains.
        TripletContainerType zT;
        const size_t nUnconstrained=assemblePeriodicConstraints(zT);
        SparseMatrixType K(Ndof,Ndof);
        K.setFromTriplets(kqqT.begin(),kqqT.end());
        SparseMatrixType Z(Ndof,dim*nUnconstrained);
        Z.setFromTriplets(zT.begin(),zT.end());
        const SparseMatrixType reduced(Z.transpose()*K*Z);
        const Eigen::VectorXd force(Z.transpose()*Fq);

        // Symmetric row-sum lumping after eliminating subdivision velocities.
        Eigen::VectorXd diagonal=Eigen::VectorXd::Zero(dim*nUnconstrained);
        for(int k=0;k<reduced.outerSize();++k)
        {
            for(typename SparseMatrixType::InnerIterator it(reduced,k);it;++it)
            {
                if(it.row()==it.col())
                {
                    diagonal(it.row())+=it.value();
                }
                else
                {
                    diagonal(it.row())+=0.5*it.value();
                    diagonal(it.col())+=0.5*it.value();
                }
            }
        }
        Eigen::VectorXd velocity=Eigen::VectorXd::Zero(dim*nUnconstrained);
        for(Eigen::Index k=0;k<diagonal.size();++k)
        {
            if(std::fabs(diagonal(k))>FLT_EPSILON)
            {
                velocity(k)=force(k)/diagonal(k);
            }
            else if(std::fabs(force(k))>FLT_EPSILON)
            {
                throw std::runtime_error("Nonzero glide force on a zero-mass degree of freedom.");
            }
        }
        return Z*velocity;
    }

    template class GalerkinGlideSolver<3>;
    
}
#endif
