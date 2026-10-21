#include <DislocationDynamics/DislocationLoopNodeIO.h>
#include <cfloat>

namespace model
{
    

    template<short unsigned int dim>
    DislocationLoopNodeIO<dim>::DislocationLoopNodeIO(const size_t& sID_in,          // sID
                                const size_t& loopID_in,
                          const VectorDim& P_in,          // position
                          const size_t& nnID,
                            const VectorDim& shift,
                              const  std::pair<short int,short int>& edgeID_in) :
        /* init */ loopID(loopID_in),
        /* init */ sID(sID_in),
        /* init */ P(P_in),
        /* init */ networkNodeID(nnID),
        /* init */ periodicShift(shift),
        /* init */ edgeIDs(edgeID_in)
        {// Constructor for MicrostructureGenerator
        }

    template<short unsigned int dim>
    DislocationLoopNodeIO<dim>::DislocationLoopNodeIO() :
        /* init */ loopID(0),
        /* init */ sID(0),
        /* init */ P(VectorDim::Zero()),
        /* init */ networkNodeID(0),
        /* init */ periodicShift(VectorDim::Zero()),
        /* init */ edgeIDs(std::make_pair((short int) -1,(short int) -1))
        {
        }

    template<short unsigned int dim>
    DislocationLoopNodeIO<dim>::DislocationLoopNodeIO(std::stringstream& ss) :
        /* init */ loopID(0),
        /* init */ sID(0),
        /* init */ P(VectorDim::Zero()),
        /* init */ networkNodeID(0),
        /* init */ periodicShift(VectorDim::Zero()),
        /* init */ edgeIDs(std::make_pair((short int)-1,(short int)-1))
        {
            ss>>loopID;
            ss>>sID;
            for(int d=0;d<dim;++d)
            {
                ss>>P(d);
            }
            ss>>networkNodeID;
            for(int d=0;d<dim;++d)
            {
                ss>>periodicShift(d);
            }
            ss>>edgeIDs.first;
            ss>>edgeIDs.second;
        }

    template<short unsigned int dim>
    std::ostream& DislocationLoopNodeIO<dim>::write(std::ostream& os) const
    {
        const auto& ds=*this;
            os  << ds.loopID<<"\t"
            /**/<< ds.sID<<"\t"
            /**/<< std::setprecision(15)<<std::scientific<<ds.P.transpose()<<"\t"
            /**/<< ds.networkNodeID<<"\t"
            /**/<< std::setprecision(15)<<std::scientific<<ds.periodicShift.transpose()<<"\t"
            /**/<< ds.edgeIDs.first<<"\t"
            /**/<< ds.edgeIDs.second<<"\t";
            return os;
        
    }

    template struct DislocationLoopNodeIO<3>;
}
