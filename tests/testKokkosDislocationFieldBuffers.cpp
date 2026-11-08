#include <DislocationDynamics/KokkosDislocationFieldBuffers.h>
#include <iostream>
#include <stdexcept>
void require(bool ok,const char* message)
{
    if(!ok)
    {
        throw std::runtime_error(message);
    }
}
template<class ExecutionSpace>
void check()
{
    using namespace model;
    KokkosDislocationFieldBuffers<3,ExecutionSpace> device;
    DislocationFieldBuffers<3> cpu;
    const DislocationQuadratureGeometry<3> point{{1,2,3},{1,0,0},0.5,2};
    cpu.addSegment({1,2},{1,2,3},{true,false,false},{point},{point},DislocationTargetPurpose::glide);
    cpu.addExternalTarget({4,5,6});
    device.uploadSources(cpu.sources());
    device.uploadTargets(cpu.targets());
    // Manufactured tensor operation exercises every input/output component;
    // it is deliberately not an elastic field kernel.
    const auto x=device.targetPositions();
    const auto y=device.sourcePositions();
    const auto q=device.sourceTensors();
    const auto output=device.stresses();
    const auto count=device.sourceCount();
    Kokkos::parallel_for("field buffer transfer test",
        Kokkos::RangePolicy<ExecutionSpace>(device.executionSpace(),0,device.targetCount()),
        KOKKOS_LAMBDA(const int t)
        {
            for(int c=0;c<9;++c)
            {
                double value=0;
                for(std::size_t s=0;s<count;++s)
                {
                    value+=(x(t,c%3)-y(s,c%3))*q(s,c);
                }
                output(t,c)=value;
            }
        });
    DislocationStressArrays<3> result;
    device.downloadStresses(result);
    for(std::size_t t=0;t<cpu.targets().size();++t)
    {
        for(int c=0;c<9;++c)
        {
            const double expected=(cpu.targets().positions.values[c%3][t]-
                cpu.sources().positions.values[c%3][0])*cpu.sources().weightedBurgersTensors.values[c][0];
            require(result.stresses.values[c][t]==expected,"CPU/backend result mismatch");
        }
    }
    const auto sourceAddress=device.sourcePositions().data();
    const auto targetAddress=device.targetPositions().data();
    DislocationTargetArrays<3> smaller;
    smaller.positions.append({7,8,9});
    device.uploadTargets(smaller);
    require(device.sourcePositions().data()==sourceAddress,"Target batch reallocated sources");
    require(device.targetPositions().data()==targetAddress,"Smaller target batch reallocated");
    device.downloadStresses(result);
    require(result.stresses.values[0].size()==1 && result.stresses.values[0][0]==0,"Stale results after upload");
    auto malformed=smaller;
    malformed.positions.values[1].clear();
    try
    {
        device.uploadTargets(malformed);
        throw std::runtime_error("Malformed target accepted");
    }
    catch(const std::invalid_argument&)
    {
        require(device.targetCount()==1,"Invalid upload changed active count");
    }
    cpu.clear();
    device.uploadSources(cpu.sources());
    device.uploadTargets(cpu.targets());
    device.downloadStresses(result);
    require(device.sourceCount()==0 && device.targetCount()==0 && result.stresses.values[0].empty(),"Empty snapshot");
    // Grow past capacity, then verify every uploaded component on the host.
    for(int i=0;i<17;++i)
    {
        cpu.addSegment({std::size_t(2*i),std::size_t(2*i+1)},{1,2,3},
            {false,false,true},{point},{},DislocationTargetPurpose::glide);
        cpu.addExternalTarget({double(i),2,3});
    }
    device.uploadSources(cpu.sources());
    device.uploadTargets(cpu.targets());
    const auto copy=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{},device.sourceTensors());
    for(int i=0;i<17;++i)
    {
        require(copy(i,6)==6.0,"Source growth lost tensor data");
    }
    require(device.sourceCount()==17 && device.targetCount()==17,"Growth counts");
    std::cout<<ExecutionSpace::name()<<" buffer tests passed\n";
}
int main(int argc,char** argv)
{
    try
    {
        model::KokkosDislocationFieldBuffers<3> uninitialized;
        throw std::runtime_error("Missing runtime initialization accepted");
    }
    catch(const std::logic_error&)
    {
    }
    Kokkos::ScopeGuard runtime(argc,argv);
    check<Kokkos::DefaultExecutionSpace>();
#ifdef KOKKOS_ENABLE_SERIAL
    check<Kokkos::Serial>();
#endif
}
