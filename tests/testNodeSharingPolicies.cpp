#include <FiniteElement.h>
#include <ClusterBoundarySinks.h>
#include <filesystem>
#include <iostream>
#include <set>
using namespace model;
void require(bool value,const char* message)
{
    if(!value) { throw std::runtime_error(message); }
}
template<int order> void check(const std::filesystem::path& root,bool periodic)
{
    SimplicialMesh<3> mesh((root/"Library/Meshes/bicrystal_48.msh").string(),Eigen::Matrix3d::Identity(),Eigen::Vector3d::Zero(),periodic?std::set<int>{-1}:std::set<int>{});
    using FE=FiniteElement<LagrangeElement<3,order>>;
    FE continuous(mesh),grain(mesh,NodeSharingPolicy::grainContinuous),discontinuous(mesh,NodeSharingPolicy::discontinuous);
    require(grain.nodeSize()>continuous.nodeSize(),"GB nodes were not split");
    require(discontinuous.nodeSize()==discontinuous.elementSize()*FE::nodesPerElement,"Discontinuous node count");
    require(discontinuous.periodicNodePairs().empty(),"Discontinuous policy added strong periodic constraints");
    std::map<size_t,std::map<size_t,size_t>> copies;
    for(const auto& node:grain.nodes())
    {
        std::set<size_t> regions;
        for(const auto* element:node) { regions.insert(element->simplex.region->regionID); }
        require(regions.size()==1,"Grain node is shared across GB");
        const auto& geometric=continuous.correspondingNode(node);
        copies[geometric.gID].emplace(*regions.begin(),node.gID);
        require(grain.nodeFinder(*regions.begin()).at(node.P0)==&node,"Grain node lookup");
    }
    for(const auto& node:discontinuous.nodes()) { require(node.size()==1,"Discontinuous node shared by elements"); }
    for(const auto& pair:grain.periodicNodePairs())
    {
        require((*grain.node(pair.source).begin())->simplex.region->regionID==(*grain.node(pair.sink).begin())->simplex.region->regionID,"Periodic tie crosses grain boundary");
    }
    auto scalar=grain.template trial<1>();
    auto smooth=continuous.template trial<1>();
    for(const auto& node:continuous.nodes()) { smooth->dofVector()[node.gID]=node.P0.sum(); }
    for(const auto& node:grain.nodes()) { scalar->dofVector()[node.gID]=10*(*node.begin())->simplex.region->regionID+node.P0.sum(); }
    for(const auto& item:grain.elements())
    {
        const auto& element=item.second;
        const Eigen::Vector4d bary=Eigen::Vector4d::Constant(.25);
        const double position=element.position(bary).sum();
        require(std::abs((*eval(smooth))(element,bary)(0)-position)<1e-12,"Foreign element read wrong continuous DOFs");
        const auto& other=continuous.elements().at(item.first);
        require(std::abs((*eval(scalar))(other,bary)(0)-(position+10*element.simplex.region->regionID))<1e-12,"Foreign element lost grain side");
    }
    const auto sinks=clusterBoundarySinks(grain);
    bool shared=false;
    for(const auto& item:copies)
    {
        if(item.second.size()>1)
        {
            shared=true;
            for(const auto& side:item.second) { require(sinks.count(side.second),"GB side missing Dirichlet sink"); }
            bool rejected=false;
            try { (*eval(scalar))(continuous.node(item.first)); } catch(const std::invalid_argument&) { rejected=true; }
            require(rejected,"Ambiguous GB node silently averaged correction");
        }
    }
    require(shared,"No shared GB nodes tested");
    // Standard mass assembly must not connect separate grain DOFs.
    const auto domain=grain.template domain<EntireDomain,4,GaussLegendre>();
    const auto mass=(test(scalar),scalar)*domain;
    const auto triplets=mass->globalTriplets();
    for(const auto& entry:triplets)
    {
        if(entry.value()!=0)
        {
            require((*grain.node(entry.row()).begin())->simplex.region->regionID==(*grain.node(entry.col()).begin())->simplex.region->regionID,"Mass assembly couples different grains");
        }
    }
}
int main(int argc,char** argv)
{
    try
    {
        require(argc==2,"Expected repository root");
        for(bool periodic:{false,true}) { check<1>(argv[1],periodic);check<2>(argv[1],periodic); }
        std::cout<<"PASS node-sharing policies, grain traces, sinks and assembly"<<std::endl;
    }
    catch(const std::exception& error) { std::cerr<<error.what()<<std::endl;return 1; }
}
