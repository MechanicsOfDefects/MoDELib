#include <DislocationDynamics/DislocationLoopIO.h>
#include <cfloat>

namespace model
{
    

    template<short unsigned int dim>
    DislocationLoopIO<dim>::DislocationLoopIO(const size_t& sID_in,         // sID
                          const VectorDim& B_in,          // position
                          const VectorDim& N_in,          // velocity
                          const VectorDim& P_in,          // velocity
                          const size_t& grainID_in,
                          const int& loopType_in
                          ) :
        /* init */ sID(sID_in)
        /* init */,B(B_in)
        /* init */,N(N_in)
        /* init */,P(P_in)
        /* init */,grainID(grainID_in)
        /* init */,loopType(loopType_in)
        /* init */,loopLength(std::make_tuple(0.0,0.0,0.0,0.0))
        /* init */,slippedArea(0.0)
        {
            
        }

    template<short unsigned int dim>
    DislocationLoopIO<dim>::DislocationLoopIO() :
        /* init */ sID(0)
        /* init */,B(VectorDim::Zero())
        /* init */,N(VectorDim::Zero())
        /* init */,P(VectorDim::Zero())
        /* init */,grainID(0)
        /* init */,loopType(0)
        /* init */,loopLength(std::make_tuple(0.0,0.0,0.0,0.0))
        /* init */,slippedArea(0.0)
        {
            
        }

    template<short unsigned int dim>
    DislocationLoopIO<dim>::DislocationLoopIO(std::stringstream& ss) :
        /* init */ sID(0)
        /* init */,B(VectorDim::Zero())
        /* init */,N(VectorDim::Zero())
        /* init */,P(VectorDim::Zero())
        /* init */,grainID(0)
        /* init */,loopType(0)
        /* init */,loopLength(std::make_tuple(0.0,0.0,0.0,0.0))
        /* init */,slippedArea(0.0)
        {
            ss>>sID;
            for(int d=0;d<dim;++d)
            {
                ss>>B(d);
            }
            for(int d=0;d<dim;++d)
            {
                ss>>N(d);
            }
            for(int d=0;d<dim;++d)
            {
                ss>>P(d);
            }
            ss>>grainID;
            ss>>loopType;
            double l1,l2,l3,A;
            ss>>l1>>l2>>l3>>A;
            loopLength=std::make_tuple(l1,l2,l3,A);
            ss>>slippedArea;
        }

    template<short unsigned int dim>
    std::ostream& DislocationLoopIO<dim>::write(std::ostream& os) const
    {
        const auto& ds=*this;
            os  << ds.sID<<"\t"
            /**/<< std::setprecision(15)<<std::scientific<<ds.B.transpose()<<"\t"
            /**/<< std::setprecision(15)<<std::scientific<<ds.N.transpose()<<"\t"
            /**/<< std::setprecision(15)<<std::scientific<<ds.P.transpose()<<"\t"
            /**/<< ds.grainID<<"\t"
            /**/<< ds.loopType<<"\t"
            /**/<< std::get<0>(ds.loopLength)<<"\t"<< std::get<1>(ds.loopLength)<<"\t"<< std::get<2>(ds.loopLength)<<"\t"<< std::get<3>(ds.loopLength)<<"\t"
            /**/<< ds.slippedArea;
            return os;
        
    }

    template struct DislocationLoopIO<3>;
}
