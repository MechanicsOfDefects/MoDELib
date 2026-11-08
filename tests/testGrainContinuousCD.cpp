#include "ClusterPeriodicFixture.h"
// A known grain-local contribution exercises the actual boundary-correction path.
struct GrainField : MicrostructureBase<3>
{
    GrainField(MicrostructureContainer<3>& c):MicrostructureBase<3>("ClusterDynamics",c) {}
    void initializeConfiguration(const DefectiveCrystalConfigIO<3>&,const std::ofstream&,const std::ofstream&) override {}
    void solve() override {}
    double getDt() const override { return 1.; }
    void output(DefectiveCrystalConfigIO<3>&,DefectiveCrystalAuxIO<3>&,std::ofstream&,std::ofstream&) const override {}
    void updateConfiguration() override {}
    MatrixDim averagePlasticDistortion() const override { return MatrixDim::Zero(); }
    MatrixDim averagePlasticDistortionRate() const override { return MatrixDim::Zero(); }
    MatrixDim averageStress() const override { return MatrixDim::Zero(); }
    VectorDim displacement(const VectorDim&,const NodeType* const,const ElementType* const,const SimplexDim* const) const override { return VectorDim::Zero(); }
    MatrixDim stress(const VectorDim&,const NodeType* const,const ElementType* const,const SimplexDim* const) const override { return MatrixDim::Zero(); }
    VectorDim inelasticDisplacementRate(const VectorDim&,const NodeType* const,const ElementType* const,const SimplexDim* const) const override { return VectorDim::Zero(); }
    VectorMSize mobileConcentration(const VectorDim&,const NodeType* const n,const ElementType* const,const SimplexDim* const) const override
    {
        require(n && n->size(),"Boundary evaluation must provide grain-side node");
        const auto grain=(*n->begin())->simplex.region->regionID;
        for(const auto* e:*n) { require(e->simplex.region->regionID==grain,"Mixed grain boundary target"); }
        return VectorMSize::Constant((grain+1)*1.e-9);
    }
    MobileClusterFluxType mobileClusterFlux(const VectorDim&,const NodeType* const,const ElementType* const,const SimplexDim* const) const override { return MobileClusterFluxType::Zero(); }
    VectorISize immobileClusters(const VectorDim&,const NodeType* const,const ElementType* const,const SimplexDim* const) const override { return VectorISize::Zero(); }
};
int main(int argc,char** argv)
{
    fs::path folder;
    try
    {
        require(argc==2,"Expected repository root");
        folder=fixture(argv[1],0,0,false,"bicrystal_48.msh");
        DefectiveCrystalBase<3> base(folder.string());
        MicrostructureContainer<3> container(base);
        auto cd=std::make_shared<ClusterDynamics<3>>(container);
        auto source=std::make_shared<GrainField>(container);
        container.microstructures().push_back(cd);
        container.microstructures().push_back(source);
        auto& fem=*cd->clusterDynamicsFEM;
        require(fem.fe->nodeSharing==NodeSharingPolicy::grainContinuous,"CD policy");
        require(base.fe->nodeSharing==NodeSharingPolicy::continuous,"Elasticity policy changed");
        fem.initializeSolver();
        cd->applyBoundaryConditions();
        fem.solveMobileClusters(true);
        constexpr int m=ClusterDynamicsParameters<3>::mSize;
        constexpr int i=ClusterDynamicsParameters<3>::iSize;
        const GrainField::VectorMSize equilibrium=fem.cdp.boundaryMobileConcentration(0.,0.);
        for(const auto& entry:fem.sinkNormals)
        {
            const auto& node=fem.fe->node(entry.first);
            const GrainField::VectorMSize total=fem.mobileClusters->dofVector().template segment<m>(m*node.gID)+source->mobileConcentration(node.P0,&node,nullptr,nullptr);
            require((total-equilibrium).norm()<1.e-16,"Total concentration differs from prescribed GB value");
        }
        DefectiveCrystalConfigIO<3> config((folder/"config").string());
        config.cdMatrix().resize(fem.fe->nodeSize(),m+i);
        config.cdMatrix().leftCols(m)=fem.mobileClusters->dofVector().reshaped(m,fem.fe->nodeSize()).transpose();
        config.cdMatrix().rightCols(i).setConstant(1.e-12);
        const auto saved=fem.mobileClusters->dofVector().eval();
        std::ofstream unused;
        fem.initializeConfiguration(config,unused,unused);
        require((saved-fem.mobileClusters->dofVector()).norm()==0.,"Grain-side restart changed mobile values");
        config.cdMatrix().resize(base.fe->nodeSize(),m+i);
        for(const auto& node:base.fe->nodes()) { config.cdMatrix().row(node.gID).setConstant((node.gID+1)*1.e-12); }
        fem.initializeConfiguration(config,unused,unused);
        for(const auto& node:fem.fe->nodes())
        {
            const auto oldID=base.fe->correspondingNode(node).gID;
            require(fem.mobileClusters->dofVector()(m*node.gID)==config.cdMatrix()(oldID,0),"Legacy restart expansion");
        }
        fs::remove_all(folder);
        std::cout<<"PASS grain-side CD boundary solve and restart migration"<<std::endl;
    }
    catch(const std::exception& e)
    {
        std::cerr<<e.what()<<std::endl;
        if(!folder.empty()) { fs::remove_all(folder); }
        return 1;
    }
}
