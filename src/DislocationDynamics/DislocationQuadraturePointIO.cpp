#include <DislocationDynamics/DislocationQuadraturePointIO.h>
#include <cfloat>

namespace model
{
    template<int dim>
    DislocationQuadraturePointIO<dim>::DislocationQuadraturePointIO() :
        /* init */ sourceID(0)
        /* init */,sinkID(0)
        /* init */,qID(0)
        /* init */,r(VectorDim::Zero())
        /* init */,j(0.0)
        /* init */,rl(VectorDim::Zero())
        /* init */,stress(MatrixDim::Zero())
        /* init */,pkForce(VectorDim::Zero())
        /* init */,stackingFaultForce(VectorDim::Zero())
        /* init */,lineTensionForce(VectorDim::Zero())
        /* init */,velocity(VectorDim::Zero())
        /* init */,elasticEnergyPerLength(0.0)
        /* init */,coreEnergyPerLength(0.0)
        /* init */,cCD(Eigen::Matrix<double,1,ClusterDynamicsParameters<dim>::mSize>::Zero())
        /* init */,cDD(Eigen::Matrix<double,1,ClusterDynamicsParameters<dim>::mSize>::Zero())
        {
        }

    

    template<int dim>
    DislocationQuadraturePointIO<dim>::DislocationQuadraturePointIO(std::stringstream& ss) :
        /* init */ sourceID(0)
        /* init */,sinkID(0)
        /* init */,qID(0)
        /* init */,r(VectorDim::Zero())
        /* init */,j(0.0)
        /* init */,rl(VectorDim::Zero())
        /* init */,stress(MatrixDim::Zero())
        /* init */,pkForce(VectorDim::Zero())
        /* init */,stackingFaultForce(VectorDim::Zero())
        /* init */,lineTensionForce(VectorDim::Zero())
        /* init */,velocity(VectorDim::Zero())
        /* init */,elasticEnergyPerLength(0.0)
        /* init */,coreEnergyPerLength(0.0)
        /* init */,cCD(Eigen::Matrix<double,1,ClusterDynamicsParameters<dim>::mSize>::Zero())
        /* init */,cDD(Eigen::Matrix<double,1,ClusterDynamicsParameters<dim>::mSize>::Zero())
        {
            
            ss>>sourceID;
            ss>>sinkID;
            ss>>qID;
            for(int d=0;d<dim;++d)
            {
                ss>>r(d);
            }
            ss>>j;
            for(int d=0;d<dim;++d)
            {
                ss>>rl(d);
            }
            for(int i=0;i<dim;++i)
            {
                for(int j=0;j<dim;++j)
                {
                    ss>>stress(i,j);
                }
            }
            for(int d=0;d<dim;++d)
            {
                ss>>pkForce(d);
            }
            for(int d=0;d<dim;++d)
            {
                ss>>stackingFaultForce(d);
            }
            for(int d=0;d<dim;++d)
            {
                ss>>lineTensionForce(d);
            }
            for(int d=0;d<dim;++d)
            {
                ss>>velocity(d);
            }
            ss>>elasticEnergyPerLength;
            ss>>coreEnergyPerLength;
            for(int d=0;d<ClusterDynamicsParameters<dim>::mSize;++d)
            {
                ss>>cCD(d);
            }
            for(int d=0;d<ClusterDynamicsParameters<dim>::mSize;++d)
            {
                ss>>cDD(d);
            }
        }

    template<int dim>
    std::ostream& DislocationQuadraturePointIO<dim>::write(std::ostream& os) const
    {
        const auto& ds=*this;
            os  << ds.sourceID<<"\t"
            /**/<< ds.sinkID<<"\t"
            /**/<< ds.qID<<"\t"
            /**/<<std::scientific
            /**/<<std::setprecision(15)
            /**/<< ds.r.transpose()<<"\t"
            /**/<< ds.j<<"\t"
            /**/<< ds.rl.transpose()<<"\t";
            for(int d=0;d<dim;++d)
            {
                os  << ds.stress.row(d)<<" ";
            }
            os  << ds.pkForce.transpose()<<"\t"
            /**/<< ds.stackingFaultForce.transpose()<<"\t"
            /**/<< ds.lineTensionForce.transpose()<<"\t"
            /**/<< ds.velocity.transpose()<<"\t"
            /**/<< ds.elasticEnergyPerLength<<"\t"
            /**/<< ds.coreEnergyPerLength<<"\t"
            /**/<< ds.cCD<<"\t"
            /**/<< ds.cDD;
            return os;
        
    }

    template struct DislocationQuadraturePointIO<3>;
}
