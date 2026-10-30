#include <ClusterDynamics.h>
#include <FixedDirichletSolver.h>
#include <filesystem>
#include <fstream>
#include <regex>
#include <chrono>
#include <iostream>
#include <iomanip>
using namespace model;
namespace fs=std::filesystem;
void require(bool ok,const std::string& message)
{
    if(!ok)
    {
        throw std::runtime_error(message);
    }
}
std::string read(const fs::path& path)
{
    std::ifstream stream(path);
    return {std::istreambuf_iterator<char>(stream),std::istreambuf_iterator<char>()};
}
void parameter(std::string& text,const std::string& key,const std::string& value)
{
    const std::regex pattern("(^|\\n)\\s*"+key+"\\s*=[^;]*;");
    require(std::regex_search(text,pattern),"Missing parameter "+key);
    text=std::regex_replace(text,pattern,"$1"+key+"="+value+";");
}
fs::path fixture(const fs::path& root,int rank,int reactions,bool source,const std::string& name="unitCube24.msh")
{
    const auto folder=fs::temp_directory_path()/("cluster-periodic-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    fs::create_directories(folder/"inputFiles");
    fs::create_directories(folder/"config");
    fs::create_directories(folder/"F");
    fs::copy_file(root/"Library/Meshes"/name,folder/"inputFiles/mesh.msh");
    auto material=read(root/"Library/Materials/Zr4_Fitted.txt");
    parameter(material,"computeReactions",reactions? "1" : "0");
    parameter(material,"doseRate_dpaPerSec",source? "1e-7" : "0");
    parameter(material,"otherSinks_SI",reactions==1? "1e14 1e14 1e14 1e14" : "0 0 0 0");
    std::string map;
    const int species[4]={-1,1,2,3};
    for(int i=0;i<4;++i)
    {
        for(int j=i;j<4;++j)
        {
            map+=std::to_string(species[i])+" "+std::to_string(species[j])+" "+((reactions==2 && i==0 && j==1)? "1\n" : "0\n");
        }
    }
    parameter(material,"reactionPrefactorMap",map);
    if(reactions==2)
    {
        parameter(material,"mobileSpeciesCascadeFractions","1 1 0 0");
    }
    std::ofstream(folder/"inputFiles/material.txt")<<material;
    std::ofstream(folder/"inputFiles/ClusterDynamics.txt")<<"useClusterDynamicsFEM=1;\n";
    const Eigen::Matrix3d F=(Eigen::Matrix3d()<<100,50,0,0,50*std::sqrt(3.0),0,0,0,100*std::sqrt(8.0/3.0)).finished();
    StaticID<Polyhedron<3,3-1>>::force_count(0);
    std::set<int> seeds;
    {
        SimplicialMesh<3> mesh((root/"Library/Meshes"/name).string(),F,Eigen::Vector3d::Zero(),{-1});
        for(const auto& region : mesh.regions())
        {
            for(const auto& entry : region.second->faces())
            {
                if(entry.second->periodicFacePair.second)
                {
                    for(int d=0;d<rank;++d)
                    {
                        if(entry.second->periodicFacePair.first.normalized().cross(F.col(d).normalized()).norm()<1.e-8)
                        {
                            seeds.insert(entry.first);
                        }
                    }
                }
            }
        }
    }
    StaticID<Polyhedron<3,3-1>>::force_count(0);
    std::ofstream poly(folder/"inputFiles/polycrystal.txt");
    poly<<std::setprecision(17)<<"materialFile=material.txt;\nmeshFile="<<"mesh.msh"<<";\nabsoluteTemperature=600;\nF=";
    for(int i=0;i<3;++i)
    {
        poly<<F.row(i)<<"\n";
    }
    poly<<";\nX0=0 0 0;\nperiodicFaceIDs=";
    for(int seed : seeds)
    {
        poly<<seed<<" ";
    }
    poly<<";\n";
    const int grains=name=="bicrystal_48.msh"? 2 : (name.find("polycrystal")!=std::string::npos? 10 : 1);
    for(int g=1;g<=grains;++g)
    {
        poly<<"C2G"<<g<<"=1 0 0 0 1 0 0 0 1;\n";
    }
    std::ofstream dc(folder/"inputFiles/DefectiveCrystal.txt");
    dc<<"physics=ClusterDynamics;\nfmmTreeDepth=0;\nuseFEM=1;\nNsteps=1;\nmaxResolveSteps=0;\ndtMax=1;\noutputFrequency=1;\noutputBinary=0;\nstartAtTimeStep=0;\nperiodicImageSize=";
    for(int d=0;d<rank;++d)
    {
        dc<<"0 ";
    }
    dc<<";\nEwaldLengthFactor=0;\n";
    return folder;
}
template<class Trial>
void periodic(const Trial& field,double relative=1.e-10)
{
    const auto& fe=field.fe();
    for(const auto& pair : fe.periodicNodePairs())
    {
        constexpr int count=Trial::dofPerNode;
        Eigen::VectorXd jump=field.dofVector().segment(count*pair.sink,count)-field.dofVector().segment(count*pair.source,count);
        if(field.periodicNodalOffset().size())
        {
            jump-=field.periodicNodalOffset().segment(count*pair.sink,count)-field.periodicNodalOffset().segment(count*pair.source,count);
        }
        require(jump.norm()<=relative*std::max(1.e-40,field.dofVector().norm()),"Periodic field mismatch");
    }
}
void diffusion(const fs::path& root,int rank)
{
    const auto folder=fixture(root,rank,0,false);
    {
        DefectiveCrystalBase<3> base(folder.string());
        ClusterDynamicsParameters<3> parameters(base);
        ClusterDynamicsFEM<3> fem(base,parameters);
        fem.initializeSolver();
        require(base.periodicLatticeBasis.cols()==rank,"Wrong active rank");
        auto u=fem.mobileClusters;
        Eigen::VectorXd target(u->gSize()),offset=Eigen::VectorXd::Zero(u->gSize());
        for(const auto& node : base.fe->nodes())
        {
            for(int c=0;c<4;++c)
            {
                offset(4*node.gID+c)=rank>0? 1.e-9*(c+1)*node.P0.x() : 0.0;
                target(4*node.gID+c)=1.e-6*(c+1)+1.e-8*(c+1)*node.P0.z()+offset(4*node.gID+c);
            }
        }
        u->setPeriodicNodalOffset(offset);
        for(auto& bc : u->dirichletConditions())
        {
            bc.second=target(bc.first);
        }
        Eigen::VectorXd legacy;
        if(rank==0)
        {
            FixedDirichletSolver<Eigen::SparseLU<Eigen::SparseMatrix<double>>,Eigen::BiCGSTAB<Eigen::SparseMatrix<double>>> solver(true,1.e-10);
            solver.compute(fem.mBWF);
            legacy=solver.solve(fem.cascadeGlobalProduction);
        }
        fem.solveMobileClusters(true);
        require((u->dofVector()-target).norm()<1.e-8*target.norm(),"Manufactured diffusion solution changed");
        if(rank==0)
        {
            require((u->dofVector()-legacy).norm()<1.e-9*target.norm(),"Nonperiodic result differs from legacy solver");
        }
        periodic(*u);
        const auto count=fem.mSolver->factorizationCount();
        for(auto& bc : u->dirichletConditions())
        {
            bc.second+=1.e-6;
        }
        fem.solveMobileClusters(true);
        require(fem.mSolver->factorizationCount()==count,"Dirichlet RHS change refactorized diffusion");
        if(rank>0)
        {
            // Two constrained copies must not silently prescribe incompatible values.
            for(const auto& pair : base.fe->periodicNodePairs())
            {
                if(u->dirichletConditions().count(4*pair.source) && u->dirichletConditions().count(4*pair.sink))
                {
                    const auto before=u->dofVector().eval();
                    u->dirichletConditions()[4*pair.sink]+=1.e-4;
                    bool rejected=false;
                    try
                    {
                        fem.solveMobileClusters(true);
                    }
                    catch(const std::runtime_error&)
                    {
                        rejected=true;
                    }
                    require(rejected && (u->dofVector()-before).norm()==0.0,"Failed boundary update modified concentration");
                    break;
                }
            }
        }
    }
    fs::remove_all(folder);
}
void closed(const fs::path& root,int reaction,bool source)
{
    const auto folder=fixture(root,3,reaction,source);
    {
        DefectiveCrystalBase<3> base(folder.string());
        ClusterDynamicsParameters<3> parameters(base);
        ClusterDynamicsFEM<3> fem(base,parameters);
        fem.initializeSolver();
        require(fem.sinkNormals.empty(),"Same-grain periodic seams became sinks");
        auto u=fem.mobileClusters;
        u->setConstant(Eigen::Vector4d::Constant(1.e-7));
        const auto before=u->dofVector().eval();
        if(reaction==0 && source)
        {
            bool rejected=false;
            try
            {
                fem.solve(true);
            }
            catch(const std::runtime_error& error)
            {
                rejected=std::string(error.what()).find("no steady state")!=std::string::npos;
            }
            require(rejected && (u->dofVector()-before).norm()==0.0,"Closed source-only solve was not rejected transactionally");
        }
        else
        {
            fem.solve(true);
            periodic(*u);
            Eigen::Vector4d expected=Eigen::Vector4d::Constant(1.e-7);
            if(reaction==1)
            {
                expected=-parameters.R1cd.fullPivLu().solve(parameters.G);
            }
            if(reaction==2)
            {
                expected(0)=expected(1)=std::sqrt(-parameters.G(0)/parameters.R2[0](0,1));
            }
            for(const auto& node : base.fe->nodes())
            {
                require((u->dofVector().segment<4>(4*node.gID)-expected).norm()<2.e-6*expected.norm(),"Closed reaction balance or conserved mean is wrong");
            }
            const auto immobile=fem.immobileClusters->dofVector().eval();
            fem.solveImmobileClusters();
            require((fem.immobileClusters->dofVector()-immobile).norm()==0.0,"Rate solve overwrote immobile concentrations");
            periodic(*fem.immobileRateTrial);
            const auto count=fem.iSolver->factorizationCount();
            fem.solveImmobileClusters();
            require(fem.iSolver->factorizationCount()==count,"Immobile mass factorization not reused");
            fem.updateImmobileClusters(0.01);
            periodic(*fem.immobileClusters);
            DefectiveCrystalConfigIO<3> config((folder/"config").string());
            config.cdMatrix().resize(base.fe->nodeSize(),12);
            config.cdMatrix().leftCols(4)=u->dofVector().reshaped(4,base.fe->nodeSize()).transpose();
            config.cdMatrix().rightCols(8)=fem.immobileClusters->dofVector().reshaped(8,base.fe->nodeSize()).transpose();
            std::ofstream unused;
            fem.initializeConfiguration(config,unused,unused);
            const auto state=fem.immobileClusters->dofVector().eval();
            config.cdMatrix()(base.fe->periodicNodePairs().front().sink,4)+=1.0;
            bool rejected=false;
            try
            {
                fem.initializeConfiguration(config,unused,unused);
            }
            catch(const std::runtime_error&)
            {
                rejected=true;
            }
            require(rejected && (fem.immobileClusters->dofVector()-state).norm()==0.0,"Inconsistent restart was not rejected transactionally");
        }
    }
    fs::remove_all(folder);
}
void sinks(const fs::path& root,const std::string& name,int rank)
{
    const auto folder=fixture(root,rank,0,false,name);
    {
        DefectiveCrystalBase<3> base(folder.string());
        MicrostructureContainer<3> container(base);
        auto physics=std::make_shared<ClusterDynamics<3>>(container);
        container.microstructures().push_back(physics);
        auto& fem=*physics->clusterDynamicsFEM;
        require(!fem.sinkNormals.empty(),"GB and free-surface sinks missing");
        bool internal=false,periodicGB=false;
        for(const auto& node : base.fe->nodes())
        {
            bool expected=false;
            for(const auto* element : node)
            {
                const auto bary=element->simplex.pos2bary(node.P0);
                for(int face=0;face<4;++face)
                {
                    if(std::abs(bary(face))>1.e-8)
                    {
                        continue;
                    }
                    const auto& simplex=element->simplex.child(face);
                    if(simplex.isRegionBoundarySimplex())
                    {
                        expected=internal=true;
                    }
                    if(simplex.isBoundarySimplex())
                    {
                        const auto found=base.mesh.periodicSimplexPairs().find(&simplex);
                        if(found==base.mesh.periodicSimplexPairs().end())
                        {
                            expected=true;
                        }
                        else if(found->second.partner->parents().begin()->second->region->regionID!=element->simplex.region->regionID)
                        {
                            expected=periodicGB=true;
                        }
                    }
                }
            }
            require(bool(fem.sinkNormals.count(node.gID))==expected,"Boundary sink classification mismatch");
        }
        require(internal && periodicGB,"Fixture must exercise internal and periodic grain boundaries");
        physics->solve();
        periodic(*fem.mobileClusters);
        require(!fem.mobileClusters->dirichletConditions().empty(),"Periodic polycrystal lost GB concentrations");
    }
    fs::remove_all(folder);
}
int main(int argc,char** argv)
{
    try
    {
        require(argc==2,"Expected repository path");
        for(int rank=0;rank<3;++rank)
        {
            diffusion(argv[1],rank);
        }
        closed(argv[1],0,false);
        closed(argv[1],0,true);
        closed(argv[1],1,true);
        closed(argv[1],2,true);
        sinks(argv[1],"bicrystal_48.msh",3);
        sinks(argv[1],"polycrystalCube_10_periodicXY.msh",2);
        std::cout<<"PASS: ClusterDynamics periodic diffusion, reactions, sinks, rates and restart\n";
    }
    catch(const std::exception& error)
    {
        std::cerr<<error.what()<<'\n';
        return 1;
    }
}
