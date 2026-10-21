#include <DislocationDynamics/DislocationSegmentIO.h>

namespace model
{
    template<short unsigned int dim>
    DislocationSegmentIO<dim>::DislocationSegmentIO(const size_t& sourceID_in,          // position
                          const size_t& sinkID_in) :
        /* init */ sourceID(sourceID_in),
        /* init */ sinkID(sinkID_in),
        /* init */ b(VectorDim::Zero()),
        /* init */ n(VectorDim::Zero()),
        /* init */ meshLocation(-1)
//        /* init */ grainID(-2) // -2 means never assigned
//        loopCounter(0),
//        isGlissile(true)
//        /* init */ meshLocation(meshLocation_in)
        {
            
            //            assert(0 && "FINISH HERE, THIS MUST BE COMPATIBLE WITH ID READER");
            
        }

    template<short unsigned int dim>
    std::ostream& DislocationSegmentIO<dim>::write(std::ostream& os) const
    {
        const auto& ds=*this;
            os  << ds.sourceID<<" "
            /**/<< ds.sinkID<<" "
            /**/<< std::setprecision(15)<<std::scientific<<ds.b.transpose()<<" "
            /**/<< std::setprecision(15)<<std::scientific<<ds.n.transpose()<<" "
            /**/<< ds.meshLocation;
            return os;
        
    }

    template struct DislocationSegmentIO<3>;
}
