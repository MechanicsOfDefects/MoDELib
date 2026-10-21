#include <DislocationDynamics/DislocationLoopLinkIO.h>
#include <cfloat>

namespace model
{
    

    template<short unsigned int dim>
    DislocationLoopLinkIO<dim>::DislocationLoopLinkIO(const size_t& loopID_in,          // sID
                          const size_t& sourceID_in,          // position
                          const size_t& sinkID_in,          // velocity
                              const bool& hasNetworkLink_in,
                          const int& meshLocation_in) :
        /* init */ loopID(loopID_in),
        /* init */ sourceID(sourceID_in),
        /* init */ sinkID(sinkID_in),
        /* init */ hasNetworkLink(hasNetworkLink_in),
        /* init */ meshLocation(meshLocation_in)
        {
            
            //            assert(0 && "FINISH HERE, THIS MUST BE COMPATIBLE WITH ID READER");
            
        }

    template<short unsigned int dim>
    DislocationLoopLinkIO<dim>::DislocationLoopLinkIO() :
        /* init */ loopID(0),
        /* init */ sourceID(0),
        /* init */ sinkID(0),
        /* init */ hasNetworkLink(true),
        /* init */ meshLocation(0)
        {
            
            
        }

    template<short unsigned int dim>
    DislocationLoopLinkIO<dim>::DislocationLoopLinkIO(std::stringstream& ss) :
        /* init */ loopID(0),
        /* init */ sourceID(0),
        /* init */ sinkID(0),
        /* init */ hasNetworkLink(true),
        /* init */ meshLocation(0)
        {
            ss>>loopID;
            ss>>sourceID;
            ss>>sinkID;
            ss>>hasNetworkLink;
            ss>>meshLocation;
        }

    template<short unsigned int dim>
    std::ostream& DislocationLoopLinkIO<dim>::write(std::ostream& os) const
    {
        const auto& ds=*this;
            os  << ds.loopID<<"\t"
            /**/<< ds.sourceID<<"\t"
            /**/<< ds.sinkID<<"\t"
            /**/<< ds.hasNetworkLink<<"\t"
            /**/<< ds.meshLocation;
            return os;
        
    }

    template struct DislocationLoopLinkIO<3>;
}
