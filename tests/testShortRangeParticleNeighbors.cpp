#include <DislocationDynamics/LazarStressParticleKernel.h>
#include <DislocationDynamics/KokkosDislocationFieldBuffers.h>
#include <FMM/KokkosShortRangeParticleNeighbors.h>
#include <random>
#include <iostream>
void require(bool b,const char* message) { if(!b) { throw std::runtime_error(message); } }
int main()
{
    try
    {
        Kokkos::ScopeGuard guard;
        std::mt19937 gen(143); std::uniform_real_distribution<double> u(.01,.99);
        for(int rank=0;rank<=3;++rank) { for(double rc:{.13,.4,1.4})
        {
            model::DislocationSourceArrays<3> source;
            model::DislocationTargetArrays<3> target;
            // Explicit opposite corners plus random points, including sources without targets.
            for(int i=0;i<40;++i)
            {
                source.positions.append(i==0?std::array<double,3>{.99,.99,.99}:std::array<double,3>{u(gen),u(gen),u(gen)});
                std::array<double,9> q; for(double& v:q) { v=u(gen)-.5; } source.weightedBurgersTensors.append(q);
            }
            for(int i=0;i<17;++i) { target.positions.append(i==0?std::array<double,3>{.01,.01,.01}:std::array<double,3>{u(gen),u(gen),u(gen)}); }
            model::ShortRangeParticleNeighbors::Basis basis(3,rank);
            if(rank>0) { basis.col(0)<<1,0,0; } if(rank>1) { basis.col(1)<<.2,1,0; } if(rank>2) { basis.col(2)<<.1,.15,1; }
            model::ShortRangeParticleNeighbors list(source.positions.values,target.positions.values,rc,basis);
            model::KokkosDislocationFieldBuffers<3> buffers;
            buffers.uploadSources(source); buffers.uploadTargets(target);
            model::KokkosShortRangeParticleNeighbors<> device(buffers.executionSpace(),list);
            model::LazarStressParticleKernel<model::LazarPart::shortRange> kernel(1,.3,.02,.08);
            device.evaluate(buffers.executionSpace(),kernel,buffers.sourcePositions(),buffers.sourceTensors(),buffers.targetPositions(),buffers.stresses(),source.size(),target.size());
            model::DislocationStressArrays<3> actual; buffers.downloadStresses(actual);
            std::size_t expectedCount=0;
            for(std::size_t t=0;t<target.size();++t)
            {
                double stress[6]={}; std::size_t perTarget=0;
                // Independent fixed replica cube, larger than required for these fixtures.
                for(int a=rank>0?-4:0;a<=(rank>0?4:0);++a) { for(int b=rank>1?-4:0;b<=(rank>1?4:0);++b) { for(int c=rank>2?-4:0;c<=(rank>2?4:0);++c)
                {
                    Eigen::Vector3d shift=Eigen::Vector3d::Zero(); int n[3]={a,b,c};
                    for(int j=0;j<rank;++j) { shift+=n[j]*basis.col(j); }
                    for(std::size_t s=0;s<source.size();++s)
                    {
                        double r[3],q[9],d2=0;
                        for(int j=0;j<3;++j) { r[j]=(target.positions.values[j][t]-source.positions.values[j][s])-shift[j]; d2+=r[j]*r[j]; }
                        if(d2>rc*rc) { continue; }
                        ++expectedCount; ++perTarget;
                        std::size_t matches=0;
                        for(std::size_t p=list.offsets[t];p<list.offsets[t+1];++p)
                        { if(list.pairs[p].source==s&&(list.shifts[list.pairs[p].image]-shift).norm()<1e-12) { ++matches; } }
                        require(matches==1,"Missing or duplicate source image");
                        for(int j=0;j<9;++j) { q[j]=source.weightedBurgersTensors.values[j][s]; }
                        kernel.addSymmetric(r,q,stress);
                    }
                } } }
                require(perTarget==list.offsets[t+1]-list.offsets[t],"Unexpected neighbor");
                int k=0;
                for(int i=0;i<3;++i) { for(int j=i;j<3;++j)
                { require(std::abs(actual.stresses.values[3*i+j][t]-stress[k])<1e-10*std::max(1.0,std::abs(stress[k])),"Neighbor stress mismatch"); ++k; } }
            }
            require(expectedCount==list.pairs.size(),"Neighbor count mismatch");
            if(rank==0)
            {
                auto moved=target.positions.values; moved[0][0]+=10;
                model::ShortRangeParticleNeighbors rebuilt(source.positions.values,moved,rc);
                require(rebuilt.offsets[1]==0,"Moved target rebuild");
            }
        } }
        using List=model::ShortRangeParticleNeighbors;
        List::Positions one; for(auto& c:one) { c={0}; }
        auto boundary=one; boundary[0]={1};
        require(List(one,boundary,1).pairs.size()==1,"Inclusive exact cutoff");
        boundary[0][0]=std::nextafter(1.,2.);
        require(List(one,boundary,1).pairs.empty(),"Outside exact cutoff");
        require(List(List::Positions{},one,1).offsets.size()==2,"Empty sources");
        require(List(one,List::Positions{},1).pairs.empty(),"Empty targets");
        bool rejected=false; try { List bad(one,one,0); } catch(const std::invalid_argument&) { rejected=true; }
        require(rejected,"Reject zero cutoff");
        List::Basis invalid=List::Basis::Zero(3,1); rejected=false;
        try { List bad(one,one,1,invalid); } catch(const std::invalid_argument&) { rejected=true; }
        require(rejected,"Reject dependent periodic generators");
        std::cout<<"Neighbor identities and Kokkos stresses match independent brute force for 0/1/2/3 periodic directions, skew cells, corners and multiple replicas."<<std::endl;
    }
    catch(const std::exception& e) { std::cerr<<e.what()<<std::endl; return 1; }
}
