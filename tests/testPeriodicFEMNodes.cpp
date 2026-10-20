#include <FiniteElement.h>
#include <PeriodicMeshFacePairing.h>
#include <Eigen/Geometry>
#include <filesystem>
#include <iostream>
#include <numeric>
#include <set>

void require(bool value,const std::string& message)
{
    if(!value)
    {
        throw std::runtime_error(message);
    }
}

void checkSimplices(const model::SimplicialMesh<3>& mesh)
{
    size_t expected=0;
    std::set<const model::PlanarMeshFace<3>*> seen;
    const double tol=1.e-7*(mesh.xMax()-mesh.xMin()).norm();
    for(const auto& region : mesh.regions())
    {
        for(const auto& entry : region.second->faces())
        {
            const auto* face=entry.second.get();
            if(!seen.insert(face).second)
            {
                continue;
            }
            for(const auto* triangle : face->simplices())
            {
                const auto found=mesh.periodicSimplexPairs().find(triangle);
                if(!face->periodicFacePair.second)
                {
                    require(found==mesh.periodicSimplexPairs().end(),"Inactive or internal simplex was paired");
                    continue;
                }
                ++expected;
                require(found!=mesh.periodicSimplexPairs().end(),"Active face triangle missing");
                const auto& pair=found->second;
                require(triangle->isBoundarySimplex() && pair.partner->isBoundarySimplex(),"Internal simplex paired");
                require(face->periodicFacePair.second->simplices().count(pair.partner),"Partner on wrong face");
                const auto& reverse=mesh.periodicSimplexPairs().at(pair.partner);
                require(reverse.partner==triangle && (reverse.shift+pair.shift).norm()<tol,"Bad reciprocal simplex pair");
                require((triangle->outNormal()+pair.partner->outNormal()).norm()<1.e-7,"Normals not antiparallel");
                for(size_t i=0;i<3;++i)
                {
                    const size_t j=pair.vertexPermutation[i];
                    require(j<3 && reverse.vertexPermutation[j]==i,"Bad inverse vertex permutation");
                    require((triangle->vertices()[i]->P0+pair.shift-pair.partner->vertices()[j]->P0).norm()<tol,"Bad vertex translation");
                }
            }
        }
    }
    require(mesh.periodicSimplexPairs().size()==expected,"Unexpected simplex pairs");
}

template<int order>
void checkNodes(model::FiniteElement<model::LagrangeElement<3,order>>& fe,int rank)
{
    const auto& mesh=fe.mesh;
    const double tol=1.e-7*(mesh.xMax()-mesh.xMin()).norm();
    using Key=std::pair<size_t,size_t>;
    std::set<Key> expected;
    // Independent geometric oracle: use plane distances and translated nodal
    // positions, not barycentric integer permutations used by production code.
    for(const auto& entry : mesh.periodicSimplexPairs())
    {
        const auto* triangle=entry.first;
        const auto& pair=entry.second;
        const auto& a=fe.elements().at(triangle->parents().begin()->first);
        const auto& b=fe.elements().at(pair.partner->parents().begin()->first);
        size_t count=0;
        for(size_t i=0;i<a.nodesPerElement;++i)
        {
            const auto& source=a.node(i);
            if(std::abs(triangle->outNormal().dot(source.P0-triangle->center()))>tol)
            {
                continue;
            }
            ++count;
            size_t matches=0;
            for(size_t j=0;j<b.nodesPerElement;++j)
            {
                const auto& sink=b.node(j);
                if((source.P0+pair.shift-sink.P0).norm()<tol)
                {
                    expected.emplace(std::minmax(source.gID,sink.gID));
                    ++matches;
                }
            }
            require(matches==1,"Geometric FEM partner is missing or ambiguous");
        }
        require(count==size_t((order+1)*(order+2)/2),"Wrong number of face nodes for polynomial order");
    }
    require(fe.periodicNodePairs().size()==expected.size(),"FEM pairing coverage or deduplication failure");
    std::vector<size_t> parents(fe.nodeSize());
    std::iota(parents.begin(),parents.end(),0);
    const auto root=[&](size_t i)
    {
        while(parents[i]!=i)
        {
            i=parents[i];
        }
        return i;
    };
    for(const auto& pair : expected)
    {
        const size_t a=root(pair.first),b=root(pair.second);
        parents[std::max(a,b)]=std::min(a,b);
    }
    std::set<size_t> vertexNodes;
    for(const auto& entry : fe.mesh2femIDmap())
    {
        vertexNodes.insert(entry.second->gID);
    }
    size_t extraNodes=0;
    for(const auto& pair : fe.periodicNodePairs())
    {
        require(expected.count({pair.source,pair.sink}),"Spurious periodic FEM pair");
        require((fe.node(pair.source).P0+pair.shift-fe.node(pair.sink).P0).norm()<tol,"Wrong FEM pair offset");
        require(model::periodicLatticeContains<3>(mesh.periodicBasis(),pair.shift,tol),"Inactive periodic direction used");
        require((fe.periodicShift(pair.source)+pair.shift-fe.periodicShift(pair.sink)).norm()<tol,"Edge/corner path inconsistency");
        extraNodes+=vertexNodes.count(pair.source)==0;
    }
    std::map<size_t,size_t> orbitSizes;
    for(size_t n=0;n<fe.nodeSize();++n)
    {
        require(fe.periodicRepresentative(n)==root(n),"Wrong periodic equivalence class");
        require((fe.node(root(n)).P0+fe.periodicShift(n)-fe.node(n).P0).norm()<tol,"Wrong representative offset");
        ++orbitSizes[root(n)];
    }
    if(rank==0)
    {
        require(expected.empty(),"Nonperiodic domain has paired FEM nodes");
    }
    else
    {
        require(!expected.empty(),"Periodic domain has no FEM pairs");
        require(order==1 || extraNodes>0,"Higher-order nodes were not paired");
    }
    // Only the box-shaped bicrystal guarantees four/eight coincident copies.
    // Non-convex cells are checked against the independent pairing graph above.
    if(rank>=2 && mesh.regions().size()==2)
    {
        size_t largest=0;
        for(const auto& orbit : orbitSizes)
        {
            largest=std::max(largest,orbit.second);
        }
        require(largest>=4,"No edge/corner equivalence class exercised");
        if(rank==3)
        {
            require(largest>=8,"Bicrystal corner equivalence class is incomplete");
        }
    }
}

template<int order>
void exerciseOrder(model::SimplicialMesh<3>& mesh,const std::set<int>& partial)
{
    model::FiniteElement<model::LagrangeElement<3,order>> fe(mesh);
    checkNodes(fe,mesh.periodicBasis().cols());
    // Exercise a two-direction subset of a three-dimensional lattice as well.
    if(mesh.periodicBasis().cols()==3)
    {
        auto twoDirections=partial;
        for(const auto& region : mesh.regions())
        {
            for(const auto& entry : region.second->faces())
            {
                const auto& pair=entry.second->periodicFacePair;
                if(pair.second && pair.first.cross(mesh.periodicBasis().col(1)).norm()<1.e-7*pair.first.norm()*mesh.periodicBasis().col(1).norm())
                {
                    twoDirections.insert(int(entry.first));
                    break;
                }
            }
            if(twoDirections.size()==2)
            {
                break;
            }
        }
        require(twoDirections.size()==2,"No second independent seed face");
        mesh.identifyParallelFaces(twoDirections);
        require(mesh.periodicBasis().cols()==2,"Two-face selection did not select two vectors");
        checkSimplices(mesh);
        fe.updatePeriodicNodePairs();
        checkNodes(fe,2);
    }
    // Explicit refresh must clear obsolete pairs and preserve FE node IDs.
    const auto nodeCount=fe.nodeSize();
    mesh.identifyParallelFaces(partial);
    checkSimplices(mesh);
    require(mesh.periodicBasis().cols()==1,"Partial selection did not select one vector");
    fe.updatePeriodicNodePairs();
    checkNodes(fe,1);
    mesh.identifyParallelFaces({});
    checkSimplices(mesh);
    fe.updatePeriodicNodePairs();
    checkNodes(fe,0);
    mesh.identifyParallelFaces({-1});
    checkSimplices(mesh);
    fe.updatePeriodicNodePairs();
    checkNodes(fe,mesh.periodicBasis().cols());
    require(fe.nodeSize()==nodeCount,"Periodicity changed FEM node numbering");
}

int main(int argc,char** argv)
{
    try
    {
        require(argc==2,"Expected repository root");
        const std::filesystem::path meshes=std::filesystem::path(argv[1])/"Library/Meshes";
        for(const auto* name : {"bicrystal_48.msh","poly30_periodicXYZ.msh","polycrystalCube_10_periodicX.msh","polycrystalCube_10_periodicXY.msh"})
        {
            // Rotated, translated coordinates ensure pairing does not assume
            // global axes or matching local tetrahedron vertex order.
            const Eigen::Matrix3d rotation=Eigen::AngleAxisd(.37,Eigen::Vector3d(1,2,3).normalized()).toRotationMatrix();
            model::SimplicialMesh<3> mesh((meshes/name).string(),rotation,Eigen::Vector3d(2,-3,1),{-1});
            checkSimplices(mesh);
            const int expectedRank=std::string(name).find("periodicXY.msh")!=std::string::npos ? 2 :
                (std::string(name).find("periodicX.msh")!=std::string::npos ? 1 : 3);
            require(mesh.periodicBasis().cols()==expectedRank,"Unexpected tutorial lattice rank");
            std::set<int> partial;
            for(const auto& region : mesh.regions())
            {
                for(const auto& entry : region.second->faces())
                {
                    const auto& pair=entry.second->periodicFacePair;
                    if(pair.second && pair.first.cross(mesh.periodicBasis().col(0)).norm()<1.e-7*pair.first.norm()*mesh.periodicBasis().col(0).norm())
                    {
                        partial={int(entry.first)};
                        break;
                    }
                }
                if(!partial.empty())
                {
                    break;
                }
            }
            require(!partial.empty(),"No single-vector seed face");
            exerciseOrder<1>(mesh,partial);
            exerciseOrder<2>(mesh,partial);
            exerciseOrder<3>(mesh,partial);
            std::cout<<"PASS periodic FEM orders 1,2,3: "<<name<<std::endl;
        }
    }
    catch(const std::exception& error)
    {
        std::cerr<<error.what()<<std::endl;
        return 1;
    }
    return 0;
}
