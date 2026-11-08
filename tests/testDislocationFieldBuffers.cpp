#include <DislocationDynamics/DislocationFieldBuffers.h>
#include <iostream>
#include <stdexcept>
#include <type_traits>

void require(bool condition, const char* message)
{
    if(!condition)
    {
        throw std::runtime_error(message);
    }
}
int main()
{
    using Buffers=model::DislocationFieldBuffers<3>;
    using Geometry=Buffers::Geometry;
    using Purpose=model::DislocationTargetPurpose;
    static_assert(std::is_trivially_copyable<Geometry>::value,"Geometry must be pointer-free numerical data");
    const Geometry point{{1,2,3},{1,0,0},0.5,2};
    Buffers buffers;
    buffers.addSegment({1,2},{0,3,0},{true,false,false},{point},{point,point},Purpose::glide);
    require(buffers.sources().size()==1 && buffers.targets().size()==2,"Independent quadrature counts");
    require(buffers.sources().weightedBurgersTensors.values[3][0]==6,"Weighted Burgers tensor");
    // Sessile, external boundary, and GB all remain sources during glide.
    buffers.addSegment({2,3},{0,3,0},{false,false,false},{point},{point},Purpose::glide);
    buffers.addSegment({3,4},{0,3,0},{true,true,false},{point},{point},Purpose::glide);
    buffers.addSegment({4,5},{0,3,0},{true,false,true},{point},{point},Purpose::glide);
    require(buffers.sources().size()==4 && buffers.targets().size()==2,"Stationary stress sources");
    // An endpoint on a face does not classify the segment as a face segment.
    buffers.addSegment({5,6},{0,3,0},{true,false,false},{point},{point},Purpose::glide);
    require(buffers.targets().size()==3,"Interior segment target retained");
    buffers.addSegment({6,7},{0,0,0},{true,false,false},{},{},Purpose::glide);
    require(buffers.sources().size()==5,"Cancelled junction must not emit stress");
    const auto sources=buffers.sources().size();
    try
    {
        buffers.addSegment({2,1},{0,-3,0},{true,false,false},{point},{point},Purpose::glide);
        throw std::runtime_error("Duplicate was accepted");
    }
    catch(const std::invalid_argument&)
    {
        require(buffers.sources().size()==sources,"Duplicate changed snapshot");
    }
    auto invalid=point;
    invalid.lineWeight=-1;
    try
    {
        buffers.addSegment({8,9},{0,3,0},{true,false,false},{invalid},{point},Purpose::glide);
        throw std::runtime_error("Invalid quadrature was accepted");
    }
    catch(const std::invalid_argument&)
    {
        require(buffers.sources().size()==sources && !buffers.segments().count({8,9}),"Invalid input changed snapshot");
    }
    buffers.addSegment({8,9},{0,3,0},{false,false,false},{point},{point},Purpose::climb);
    require(buffers.segments().at({8,9}).targets.size==1,"Climb eligibility");
    buffers.addSegment({9,10},{0,3,0},{false,false,true},{point},{point},Purpose::diagnostic);
    require(buffers.targetBindings().back().point.purpose==Purpose::diagnostic,"Diagnostics are not mobility");
    const auto external=buffers.addExternalTarget({4,5,6});
    require(buffers.targetBindings()[external].external,"External target binding");
    model::DislocationStressArrays<3> result;
    result.reset(buffers.targets().size());
    require(result.stresses.values[8].size()==buffers.targets().size(),"Output extent");
    // Subdivision conserves the integrated tensor; orientation reversal is invariant.
    Buffers reversed;
    auto reverse=point;
    reverse.tangent={-1,0,0};
    reverse.lineWeight=1;
    reversed.addSegment({2,1},{0,-3,0},{true,false,false},{reverse,reverse},{},Purpose::glide);
    require(reversed.sources().weightedBurgersTensors.values[3][0]+
            reversed.sources().weightedBurgersTensors.values[3][1]==6,"Orientation/subdivision conservation");
    const auto capacity=buffers.sources().positions.values[0].capacity();
    buffers.clear();
    require(buffers.sources().size()==0 && buffers.targets().size()==0 && buffers.segments().empty(),"Clear metadata");
    require(buffers.sources().positions.values[0].capacity()==capacity,"Reuse capacity");
    std::cout<<"Dislocation source/target buffer tests passed\n";
}
