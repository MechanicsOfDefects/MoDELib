#ifdef MODELIB_TEST_KOKKOS
#include <FMM/KokkosBoxParticleIndex.h>
#endif
#include <FMM/BoxParticleIndex.h>
#include <iostream>
#include <limits>
void require(bool ok,const char* message)
{
    if(!ok) { throw std::runtime_error(message); }
}
template<int dim> void check()
{
    using Tree=model::BoxTree<dim>;
    Tree tree(Tree::Vector::Constant(-2),Tree::Vector::Constant(6),2);
    std::array<std::vector<double>,dim> points;
    for(const auto& leaf:tree.boxes(tree.depth()))
    {
        for(int a=0;a<dim;++a) { points[a].push_back(leaf->center()[a]); }
    }
    for(double v:{-2.,0.,2.,4.,6.,2.})
    {
        for(int a=0;a<dim;++a) { points[a].push_back(v); }
    }
    model::BoxParticleIndex<dim> index;
    std::array<std::vector<double>,dim> empty;
    index.rebuild(tree,points,empty);
    require(index.matches(tree),"Revision mismatch");
    const auto& data=index.sources();
    require(data.offsets.back()==points[0].size(),"Missing sources");
    std::vector<bool> seen(points[0].size(),false);
    for(std::size_t leaf=0;leaf<tree.boxes(tree.depth()).size();++leaf)
    {
        std::size_t previous=0;
        for(auto j=data.offsets[leaf];j<data.offsets[leaf+1];++j)
        {
            const auto p=data.permutation[j];
            require(!seen[p],"Duplicate point"); seen[p]=true;
            require(j==data.offsets[leaf] || p>previous,"Unstable ordering"); previous=p;
            typename Tree::Vector x;
            for(int a=0;a<dim;++a) { x[a]=points[a][p]; }
            require(tree.locate(x,tree.depth())==tree.boxes(tree.depth())[leaf],"Wrong leaf");
            require(data.leafOfPoint[p]==leaf,"Reverse lookup mismatch");
        }
    }
    const auto old=index.sources().permutation;
    for(double bad:{7.,std::numeric_limits<double>::quiet_NaN()})
    {
        auto invalid=points; invalid[0][0]=bad;
        bool failed=false;
        try { index.rebuild(tree,empty,invalid); }
        catch(const std::out_of_range&) { failed=true; }
        require(failed && index.sources().permutation==old,"Failed update mutated state");
    }
    index.rebuild(tree,points,points);
#ifdef MODELIB_TEST_KOKKOS
    model::KokkosBoxParticleIndex<dim> device;
    device.upload(index,tree);
    auto copy=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),device.targetPermutation());
    for(std::size_t i=0;i<old.size();++i) { require(copy(i)==old[i],"Device permutation differs"); }
#endif
#ifdef MODELIB_TEST_KOKKOS
    auto offsets=device.sourceOffsets();
    auto permutation=device.sourcePermutation();
    Kokkos::View<std::size_t*> visits("Point visits",points[0].size());
    Kokkos::deep_copy(visits,std::size_t(0));
    Kokkos::parallel_for("Visit leaf sources",tree.boxes(tree.depth()).size(),KOKKOS_LAMBDA(const int leaf)
    {
        for(auto j=offsets(leaf);j<offsets(leaf+1);++j)
        {
            Kokkos::atomic_increment(&visits(permutation(j)));
        }
    });
    auto hostVisits=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),visits);
    for(std::size_t p=0;p<points[0].size();++p) { require(hostVisits(p)==1,"Device traversal lost or duplicated source"); }
#endif
    tree.setDepth(0);
    require(!index.matches(tree),"Stale tree accepted");
#ifdef MODELIB_TEST_KOKKOS
    bool rejected=false;
    try { device.upload(index,tree); } catch(const std::logic_error&) { rejected=true; }
    require(rejected,"Stale device upload accepted");
#endif
    index.rebuild(tree,points,empty);
    require(index.sources().offsets.size()==2,"Depth zero");
    index.rebuild(tree,empty,empty);
    require(index.sources().permutation.empty(),"Empty rebuild");
}
int main(int argc,char** argv)
{
#ifdef MODELIB_TEST_KOKKOS
    Kokkos::ScopeGuard guard(argc,argv);
#endif
    check<1>(); check<2>(); check<3>();
    std::cout<<"Particle leaf assignment passed\n";
}
