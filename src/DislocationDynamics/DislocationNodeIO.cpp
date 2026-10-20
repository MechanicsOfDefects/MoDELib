#include <DislocationDynamics/DislocationNodeIO.h>
#include <cfloat>

namespace model
{
    

    template<short unsigned int dim>
    DislocationNodeIO<dim>::DislocationNodeIO(const size_t& sID_in,          // sID
                          const VectorDim& P_in,          // position
                          const VectorDim& V_in,          // velocity
                          const ClimbVelocityScalarType& climbVelocityScalar_in,
                          const double& velocityReduction_in,             // velocity reduction factor
//                          const size_t& snID_in,          // component ID
                          const int&  meshLocation_in) :
        /* init */ sID(sID_in),
        /* init */ P(P_in),
        /* init */ V(V_in),
        /* init */ climbVelocityScalar(climbVelocityScalar_in),
        /* init */ velocityReduction(velocityReduction_in),
//        /* init */ snID(snID_in),
        /* init */ meshLocation(meshLocation_in)
//        /* init */ masterID(sID)
        {// Constructor for MicrostructureGenerator
        }

    template<short unsigned int dim>
    DislocationNodeIO<dim>::DislocationNodeIO() :
        /* init */ sID(0),
        /* init */ P(VectorDim::Zero()),
        /* init */ V(VectorDim::Zero()),
        /* init */ climbVelocityScalar(ClimbVelocityScalarType::Zero()),
        /* init */ velocityReduction(1.0),
//        /* init */ snID(0),
        /* init */ meshLocation(0)
//        /* init */ masterID(0)
        {
        }

    template<short unsigned int dim>
    DislocationNodeIO<dim>::DislocationNodeIO(std::stringstream& ss) :
        /* init */ sID(0),
        /* init */ P(VectorDim::Zero()),
        /* init */ V(VectorDim::Zero()),
        /* init */ climbVelocityScalar(ClimbVelocityScalarType::Zero()),
        /* init */ velocityReduction(1.0),
//        /* init */ snID(0),
        /* init */ meshLocation(0)
//        /* init */ masterID(0)
        {
            ss>>sID;
            for(int d=0;d<dim;++d)
            {
                ss>>P(d);
            }
            for(int d=0;d<dim;++d)
            {
                ss>>V(d);
            }
            for(int k=0;k<ClusterDynamicsParameters<dim>::mSize;++k)
            {
                ss>>climbVelocityScalar(k);
            }
            ss>>velocityReduction;
//            ss>>snID;
            ss>>meshLocation;
//            ss>>masterID;
        }

    template<short unsigned int dim>
    std::ostream& DislocationNodeIO<dim>::write(std::ostream& os) const
    {
        const auto& ds=*this;
            os  << ds.sID<<"\t"
            /**/<< std::setprecision(15)<<std::scientific<<ds.P.transpose()<<"\t"
            /**/<< ds.V.transpose()<<"\t"
            /**/<< ds.climbVelocityScalar<<"\t"
            /**/<< ds.velocityReduction<<"\t"
//            /**/<< ds.snID<<"\t"
            /**/<< ds.meshLocation;
//            /**/<< ds.masterID;
            return os;
        
    }

    template struct DislocationNodeIO<3>;
}
