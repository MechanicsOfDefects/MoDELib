#include <DislocationNetwork.h>
#include <DislocationQuadraturePoint.h>
#include <DefectiveCrystal/DefectiveCrystalConfigIO.h>
#include <DefectiveCrystal/DefectiveCrystalAuxIO.h>

namespace model
{
template<short unsigned int dim>
    DislocationNodeIO<dim>::DislocationNodeIO(const DislocationNode<dim>& dn) :
        /* init */ sID(dn.sID),
        /* init */ P(dn.get_P()),
        /* init */ V(dn.get_V()),
        /* init */ climbVelocityScalar(dn.climbVelocityScalar),
        /* init */ velocityReduction(dn.velocityReduction()),
//        /* init */ snID(dn.pSN()->sID),
        /* init */ meshLocation(dn.meshLocation())
//        /* init */ masterID(dn.masterNode? dn.masterNode->sID : sID)
        {
         
            
        }

template<short unsigned int dim>
    DislocationLoopIO<dim>::DislocationLoopIO(const DislocationLoop<dim>& dL) :
        /* init */ sID(dL.sID)
        /* init */,B(dL.flow().cartesian())
        /* init */,N(dL.glidePlane? (dL.slippedArea()>FLT_EPSILON? dL.rightHandedUnitNormal() : dL.glidePlane->unitNormal) : VectorDim::Zero())
        /* init */,P(dL.glidePlane? dL.glidePlane->P : (*dL.loopLinks().begin())->source->get_P() )
        /* init */,grainID(dL.grain.grainID)
        /* init */,loopType(dL.loopType)
        /* init */,loopLength(dL.network().outputLoopLength? dL.loopLength() : std::make_tuple(0.0,0.0,0.0,0.0))
        /* init */,slippedArea(dL.slippedArea())
        {
            
        }

template<short unsigned int dim>
    DislocationLoopNodeIO<dim>::DislocationLoopNodeIO(const DislocationLoopNode<dim>& dn) :
        /* init */ loopID(dn.loop()->sID)
        /* init */,sID(dn.sID)
        /* init */,P(dn.get_P())
        /* init */,networkNodeID(dn.networkNode->sID)
        /* init */,periodicShift(dn.periodicPlanePatch()? dn.periodicPlanePatch()->shift : VectorDim::Zero())
        /* init */,edgeIDs(dn.periodicPlaneEdge.first? (dn.periodicPlaneEdge.second ? std::make_pair(dn.periodicPlaneEdge.first->edgeID,dn.periodicPlaneEdge.second->edgeID) : std::make_pair(dn.periodicPlaneEdge.first->edgeID,(short int) -1))  : std::make_pair((short int) -1,(short int) -1))
        {
         
            
        }

template<short unsigned int dim>
    DislocationLoopLinkIO<dim>::DislocationLoopLinkIO(const DislocationLoopLink<dim>& ll) :
        /* init */ loopID(ll.loop->sID),
        /* init */ sourceID(ll.source->sID),
        /* init */ sinkID(ll.sink->sID),
        /* init */ hasNetworkLink(ll.hasNetworkLink()),
        /* init */ meshLocation(ll.networkLink()? ll.networkLink()->meshLocation() : TypeTraits<DislocationLoopLink<dim>>::onMeshBoundary)
        {
         
            
        }

template<int dim>
    DislocationQuadraturePointIO<dim>::DislocationQuadraturePointIO(const DislocationQuadraturePoint<dim>& qPoint) :
        /* init */ sourceID(qPoint.sourceID)
        /* init */,sinkID(qPoint.sinkID)
        /* init */,qID(qPoint.qID)
        /* init */,r(qPoint.r)
        /* init */,j(qPoint.j)
        /* init */,rl(qPoint.rl)
        /* init */,stress(qPoint.stress)
        /* init */,pkForce(qPoint.pkForce)
        /* init */,stackingFaultForce(qPoint.stackingFaultForce)
        /* init */,lineTensionForce(qPoint.lineTensionForce)
        /* init */,velocity(qPoint.velocity)
        /* init */,elasticEnergyPerLength(qPoint.elasticEnergyPerLength)
        /* init */,coreEnergyPerLength(qPoint.coreEnergyPerLength)
        /* init */,cCD(qPoint.cCD)
        /* init */,cDD(qPoint.cDD)
        {
            
        }
    template DislocationNodeIO<3>::DislocationNodeIO(const DislocationNode<3>&);
    template DislocationLoopIO<3>::DislocationLoopIO(const DislocationLoop<3>&);
    template DislocationLoopNodeIO<3>::DislocationLoopNodeIO(const DislocationLoopNode<3>&);
    template DislocationLoopLinkIO<3>::DislocationLoopLinkIO(const DislocationLoopLink<3>&);
    template DislocationQuadraturePointIO<3>::DislocationQuadraturePointIO(const DislocationQuadraturePoint<3>&);
}
