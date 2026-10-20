#include <DislocationDynamics/DislocationLinkingNumber.h>
#include <DislocationNetwork.h>

namespace model
{
template<int dim>
typename LinkingNumber<dim>::VectorDim LinkingNumber<dim>::unitCross(const VectorDim& v1,
                                   const VectorDim& v2)
        {
            VectorDim temp(v1.cross(v2));
            const double tempNorm(temp.norm());
            return (tempNorm<tol)? VectorDim::Zero() : (temp/tempNorm).eval();
        }

template<int dim>
double LinkingNumber<dim>::myAsin(const double& x)
        {
            if(x>=1.0)
            {
                return 0.5*std::numbers::pi;
            }
            else if(x<=-1.0)
            {
                return -0.5*std::numbers::pi;
            }
            else
            {
                return asin(x);
            }
        }

template<int dim>
double LinkingNumber<dim>::segmentPairLN(const VectorDim& R1,
                                    const VectorDim& R2,
                                    const VectorDim& R3,
                                    const VectorDim& R4)
        {
            //            const VectorDim& R1(link1.source->get_P());
            //            const VectorDim& R2(link1.  sink->get_P());
            //            const VectorDim& R3(link2.source->get_P());
            //            const VectorDim& R4(link2.  sink->get_P());
            
            const VectorDim R13(R3-R1);
            const VectorDim R12(R2-R1);
            const VectorDim R34(R4-R3);
            
            const double R13norm(R13.norm());
            const double R12norm(R12.norm());
            const double R34norm(R34.norm());
            
            if(R13norm>tol && R12norm>tol && R34norm>tol)
            {
                const double prod(R34.cross(R12).dot(R13)/R13norm/R12norm/R34norm);
                if(fabs(prod)>tol)
                {
                    const VectorDim R14(R4-R1);
                    const VectorDim R24(R4-R2);
                    const VectorDim R23(R3-R2);
                    
                    const VectorDim n1(unitCross(R13,R14));
                    const VectorDim n2(unitCross(R14,R24));
                    const VectorDim n3(unitCross(R24,R23));
                    const VectorDim n4(unitCross(R23,R13));
                    
                    
                    const double Os= myAsin(n1.dot(n2))
                    /*            */+myAsin(n2.dot(n3))
                    /*            */+myAsin(n3.dot(n4))
                    /*            */+myAsin(n4.dot(n1));
                    
                    return Os/(4.0*std::numbers::pi)*prod/fabs(prod);
                    
                }
                else
                {
                    return 0.0;
                }
            }
            else
            {
                return 0.0;
            }
        }

template<int dim>
double LinkingNumber<dim>::loopPairLN(const std::vector<VectorDim>& loop1,
                                 const std::vector<VectorDim>& loop2)
        {
            
            double temp(0.0);
            
            for(size_t i=0;i<loop1.size();++i)
            {
                
                const size_t i1(i==loop1.size()-1? 0 : i+1);
                
                for(size_t j=0;j<loop2.size();++j)
                {
                    
                    const size_t j1(j==loop2.size()-1? 0 : j+1);
                    
                    temp+=segmentPairLN(loop1[i],
                                        loop1[i1],
                                        loop2[j],
                                        loop2[j1]);
                }
            }
            return temp;
        }

template<int dim>
double LinkingNumber<dim>::loopPairHelicity(const std::vector<VectorDim>& loop1,
                                       const VectorDim& b1,
                                       const std::vector<VectorDim>& loop2,
                                       const VectorDim& b2)
        {
            
            
            return loopPairLN(loop1,loop2)*b1.dot(b2);
            
        }

template<int dim>
double DislocationLinkingNumber<dim>::segmentPairLN(const LoopLink& link1,
                                    const LoopLink& link2)
        {
            
            const VectorDim& R1(link1.source->get_P());
            const VectorDim& R2(link1.  sink->get_P());
            const VectorDim& R3(link2.source->get_P());
            const VectorDim& R4(link2.  sink->get_P());
            
            if(   (link1.source->sID==link2.source->sID && link1.sink->sID==link2.sink->sID)
               || (link1.source->sID==link2.sink->sID && link1.sink->sID==link2.source->sID))
            {
                return 0.0;
            }
            else
            {
                return LinkingNumber<dim>::segmentPairLN(R1,R2,R3,R4);
            }
        }

template<int dim>
std::pair<double,bool> DislocationLinkingNumber<dim>::linkingNumber(const LoopType& loop1,
                                                    const LoopType& loop2)
        {
            double Ln=0.0;
            
            bool loopsTouch=false;
            
            for(const auto& link1 : loop1.loopLinks())
            {
                for(const auto& link2 : loop2.loopLinks())
                {
                    if(SegmentSegmentDistance<dim>(link1->source->get_P(),
                                                   link1->sink->get_P(),
                                                   link2->source->get_P(),
                                                   link2->sink->get_P()).dMin<FLT_EPSILON)
                    {
                        loopsTouch=true;
                    }
                    
                    Ln +=  segmentPairLN(*link1,*link2);
                    
                }
            }
            
            return std::make_pair(Ln,loopsTouch);
        }

template<int dim>
DislocationLinkingNumber<dim>::DislocationLinkingNumber(const DislocationNetwork<dim>& DN_in) :
        /* init */ DN(DN_in)
        {
            
            
        }

template<int dim>
std::ostream& DislocationLinkingNumber<dim>::write(std::ostream& os) const
{
    const auto& ln=*this;
            
            for(auto loopIter1=ln.DN.loops().begin();loopIter1!=ln.DN.loops().end();++loopIter1)
            {
                for(auto loopIter2=loopIter1;loopIter2!=ln.DN.loops().end();++loopIter2)
                {
                    if(loopIter1!=loopIter2)
                    {
                        const std::pair<double,bool> Ln(linkingNumber(*loopIter1->second.lock(),*loopIter2->second.lock()));
                        
                        os  <<loopIter1->first<<" "<<loopIter2->first;
                        os  <<" "<<Ln.second<<" "<<Ln.first;
                        os  <<" "<<loopIter1->second.lock()->burgers().dot(loopIter2->second.lock()->burgers());
                        //                        os  <<" "<<loopIter1->second->burgers().transpose();
                        //                        os  <<" "<<loopIter2->second->burgers().transpose();
                        os  <<"\n";
                    }
                }
            }
            
            return os;
        
}

 template struct LinkingNumber<3>;
 template class DislocationLinkingNumber<3>;
}
