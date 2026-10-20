#include <DislocationDynamics/DislocationNodeIO.h>
#include <DislocationDynamics/DislocationLoopIO.h>
#include <DislocationDynamics/DislocationLoopNodeIO.h>
#include <DislocationDynamics/DislocationLoopLinkIO.h>
#include <DislocationDynamics/DislocationQuadraturePointIO.h>
#include <DislocationDynamics/DislocationSegmentIO.h>
#include <DislocationNetwork.h>
#include <DefectiveCrystal/DefectiveCrystalAuxIO.h>
#include <DislocationDynamics/DislocationLinkingNumber.h>
#include <type_traits>
#include <sstream>
#include <stdexcept>

namespace
{
    struct Unrelated {};
    using namespace model;
    static_assert(std::is_constructible_v<DislocationNodeIO<3>,const DislocationNode<3>&>);
    static_assert(std::is_constructible_v<DislocationLoopIO<3>,const DislocationLoop<3>&>);
    static_assert(std::is_constructible_v<DislocationLoopNodeIO<3>,const DislocationLoopNode<3>&>);
    static_assert(std::is_constructible_v<DislocationLoopLinkIO<3>,const DislocationLoopLink<3>&>);
    static_assert(std::is_constructible_v<DislocationQuadraturePointIO<3>,const DislocationQuadraturePoint<3>&>);
    static_assert(!std::is_constructible_v<DislocationNodeIO<3>,const Unrelated&>);
    static_assert(!std::is_constructible_v<DislocationLoopIO<3>,const Unrelated&>);
    static_assert(!std::is_constructible_v<DislocationLoopNodeIO<3>,const Unrelated&>);
    static_assert(!std::is_constructible_v<DislocationLoopLinkIO<3>,const Unrelated&>);
    static_assert(!std::is_constructible_v<DislocationQuadraturePointIO<3>,const Unrelated&>);
    static_assert(std::is_same_v<decltype(std::declval<DislocationNetwork<3>>().networkRemesher),DislocationNetworkRemesh<3>>);
    static_assert(std::is_same_v<decltype(std::declval<DislocationNetwork<3>>().junctionsMaker),DislocationJunctionFormation<3>>);
    static_assert(std::is_same_v<decltype(std::declval<DislocationNetwork<3>>().nodeContractor),DislocationNodeContraction<3>>);
    static_assert(std::is_same_v<decltype(std::declval<DislocationNetwork<3>>().crossSlipMaker),DislocationCrossSlip<3>>);
    template<class Record>
    void roundTrip(const Record& record)
    {
        std::ostringstream written;
        written<<record;
        std::stringstream input(written.str());
        Record restored(input);
        std::ostringstream repeated;
        repeated<<restored;
        if(written.str()!=repeated.str())
        {
            throw std::runtime_error("DD IO text round trip changed serialized fields");
        }
    }
}
int main()
{
    model::DefectiveCrystalAuxIO<3> aux;
    if(aux.getBinFilename(7)!="config/ddAux_7.bin")
    {
        throw std::runtime_error("Default auxiliary folder is wrong");
    }
    model::DislocationNodeIO<3> node;
    node.sID=71;
    node.P=Eigen::Vector3d(1.25,-2.5,3.75);
    node.V=Eigen::Vector3d(-4,5,-6);
    node.meshLocation=2;
    roundTrip(node);
    model::DislocationLoopIO<3> loop;
    loop.sID=9;
    loop.B=Eigen::Vector3d(1,2,3);
    roundTrip(loop);
    roundTrip(model::DislocationLoopNodeIO<3>());
    roundTrip(model::DislocationLoopLinkIO<3>(9,71,72,true,2));
    model::DislocationQuadraturePointIO<3> point;
    point.sourceID=71;
    point.sinkID=72;
    point.r=node.P;
    point.pkForce=node.V;
    roundTrip(point);
    const Eigen::Vector3d a(0,0,0),b(1,0,0),c(0,1,0),d(1,1,0);
    if(model::LinkingNumber<3>::segmentPairLN(a,b,c,d)!=0.0)
    {
        throw std::runtime_error("Separated coplanar segments have nonzero linking number");
    }
}
