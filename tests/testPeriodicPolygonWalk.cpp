#include <DislocationDynamics/DislocationDynamicsBase.h>
#include <DefectiveCrystal/DefectiveCrystalBase.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <iostream>

void require(bool condition, const char* message)
{
    if(!condition)
    {
        throw std::runtime_error(message);
    }
}

struct OriginalNode {};

void exercise(model::DefectiveCrystalBase<3>& base,
              const std::vector<Eigen::Vector2d>& vertices, double expectedArea, bool reverse=false)
{
    const double length=base.mesh.xMax()[0]-base.mesh.xMin()[0];
    const auto& grain=base.poly.grains.begin()->second;
    const auto normal=grain->reciprocalLatticeDirection(Eigen::Vector3d::UnitZ());
    auto plane=base.dislocationDynamics().periodicGlidePlaneFactory.get(model::GlidePlaneKey<3>(Eigen::Vector3d::Zero(),normal));
    std::vector<OriginalNode> nodes(vertices.size());
    std::vector<std::pair<Eigen::Vector3d,const OriginalNode* const>> polygon;
    for(size_t k=0; k<vertices.size(); ++k)
    {
        const size_t i=reverse ? vertices.size()-1-k : k;
        polygon.emplace_back(length*Eigen::Vector3d(vertices[i][0],vertices[i][1],0),&nodes[i]);
    }
    const auto records=plane->polygonPatchIntersection(polygon,true);
    require(records.size()>=3,"Walk returned no contour");
    double area=0.0;
    for(size_t i=0; i<records.size(); ++i)
    {
        const Eigen::Vector3d p=plane->referencePlane->globalPosition(std::get<0>(records[i]));
        const Eigen::Vector3d q=plane->referencePlane->globalPosition(std::get<0>(records[(i+1)%records.size()]));
        area+=p.cross(q).z()/2.0;
        const auto& provenance=std::get<6>(records[i]);
        using Kind=model::PolygonPatchNodeKind;
        using Arc=model::PolygonPatchArcKind;
        const auto verifySegment=[&](const model::PolygonOriginalSegment<OriginalNode>& segment)
        {
            require(segment.sourceIndex<polygon.size() && segment.sinkIndex==(segment.sourceIndex+1)%polygon.size(),
                    "Invalid original directed segment indices");
            require(segment.source==polygon[segment.sourceIndex].second && segment.sink==polygon[segment.sinkIndex].second,
                    "Segment provenance lost its original node pointers");
        };
        if(provenance.polygonLocation)
        {
            const auto& location=*provenance.polygonLocation;
            verifySegment(location.segment);
            require(location.parameter>=0 && location.parameter<=1,"Invalid original segment parameter");
            const Eigen::Vector3d expected=(1-location.parameter)*polygon[location.segment.sourceIndex].first+
                                         location.parameter*polygon[location.segment.sinkIndex].first;
            require((expected-p).norm()<1.e-6,"Original segment provenance does not reconstruct the intersection");
        }
        const auto metadataPatch=plane->getPatch(std::get<1>(records[i]));
        bool physical=false;
        size_t incidentEdges=0;
        for(const auto& edge : metadataPatch->edges())
        {
            const Eigen::Vector2d d=*edge->sink-*edge->source;
            const Eigen::Vector2d delta=std::get<0>(records[i])-*edge->source;
            const double t=delta.dot(d)/d.squaredNorm();
            if(t>=-1.e-8 && t<=1+1.e-8 && (delta-t*d).norm()<1.e-6)
            {
                ++incidentEdges;
                for(const auto* face : edge->meshIntersection->faces)
                {
                    physical=physical || !face->periodicFacePair.second ||
                        face->regionIDs!=face->periodicFacePair.second->regionIDs;
                }
            }
        }
        require(provenance.onPhysicalBoundary==physical && provenance.atPatchCorner==(incidentEdges>1),
                "Boundary metadata disagrees with actual mesh faces");
        bool seam=false;
        for(size_t neighbor : {(i+records.size()-1)%records.size(),(i+1)%records.size()})
        {
            seam=seam || ((std::get<0>(records[i])-std::get<0>(records[neighbor])).norm()<1.e-6 &&
                         (std::get<1>(records[i])-std::get<1>(records[neighbor])).norm()>1.e-6);
        }
        require(provenance.crossesPeriodicSeam==seam,"Seam classification disagrees with contour patch transitions");
        if(provenance.originalVertex)
        {
            require(provenance.kind==Kind::OriginalVertex &&
                    std::get<5>(records[i])==polygon[*provenance.originalVertex].second,
                    "Original vertex identity changed");
        }
        else if(provenance.onPhysicalBoundary)
        {
            require(provenance.kind==(provenance.polygonLocation ? Kind::PhysicalBoundaryIntersection : Kind::PhysicalBoundaryCorner),
                    "Physical boundary node misclassified");
        }
        else
        {
            require(provenance.kind==Kind::PeriodicCrossing && provenance.crossesPeriodicSeam,
                    "Periodic crossing misclassified");
        }
        if(provenance.incomingArc==Arc::OriginalSegment)
        {
            require(provenance.incomingOriginalSegment.has_value(),"Polygon arc has no original segment");
            verifySegment(*provenance.incomingOriginalSegment);
            const auto& origin=*provenance.incomingOriginalSegment;
            const Eigen::Vector3d a=polygon[origin.sourceIndex].first;
            const Eigen::Vector3d d=polygon[origin.sinkIndex].first-a;
            const Eigen::Vector3d previous=plane->referencePlane->globalPosition(std::get<0>(records[(i+records.size()-1)%records.size()]));
            require((p-a).cross(d).norm()<1.e-6 && (previous-a).cross(d).norm()<1.e-6,
                    "Incoming polygon arc attributed to wrong segment");
        }
        else
        {
            require(!provenance.incomingOriginalSegment,"Boundary or seam arc falsely attributed to original segment");
            if(provenance.incomingArc==Arc::PeriodicSeam)
            {
                require(provenance.crossesPeriodicSeam &&
                        (std::get<0>(records[i])-std::get<0>(records[(i+records.size()-1)%records.size()])).norm()<1.e-6,
                        "Seam arc is not a zero-length periodic transition");
            }
        }

        const Eigen::Vector3d mapped=p+std::get<1>(records[i]);
        require((mapped.array()>=base.mesh.xMin().array()-1.e-6).all() &&
                (mapped.array()<=base.mesh.xMax().array()+1.e-6).all(),"Mapped node outside RVE");
        if(!std::get<5>(records[i]))
        {
            require(std::get<2>(records[i]).first>=0,"Inserted boundary node lacks edge ID");
            const auto patch=plane->getPatch(std::get<1>(records[i]));
            for(short id : {std::get<2>(records[i]).first,std::get<2>(records[i]).second})
            {
                if(id>=0)
                {
                    const auto edge=patch->edges().at(id);
                    const Eigen::Vector2d direction=*edge->sink-*edge->source;
                    const Eigen::Vector2d delta=std::get<0>(records[i])-*edge->source;
                    const double t=delta.dot(direction)/direction.squaredNorm();
                    require(t>=-1.e-8 && t<=1+1.e-8 && (delta-t*direction).norm()<1.e-6,
                            "Intersection record references the wrong patch edge");
                }
            }
        }
    }
    require(std::abs(area/(length*length)-(reverse ? -expectedArea : expectedArea))<1.e-8,
            "Wrong retained area or orientation");
    const bool retainsOriginal=std::any_of(records.begin(),records.end(),[](const auto& record)
    {
        return std::get<5>(record)!=nullptr;
    });
    require(!retainsOriginal || std::get<5>(records.back())!=nullptr,"Contour must finish at an original node");
    // Geometry must also be checked in the runtime mode; this dummy deliberately
    // has no periodicPlanePatch() accessor, so null metadata is safe too.
    const auto runtime=plane->polygonPatchIntersection(polygon,false);
    require(runtime.size()==records.size(),"Runtime walk differs from generator walk");
    std::vector<std::pair<Eigen::Vector2d,const OriginalNode* const>> localPolygon;
    for(const auto& vertex : polygon)
    {
        localPolygon.emplace_back(plane->referencePlane->localPosition(vertex.first),vertex.second);
    }
    const auto localRecords=plane->polygonPatchIntersection(localPolygon,true);
    for(size_t i=0; i<records.size(); ++i)
    {
        const auto& a=std::get<6>(records[i]);
        const auto& b=std::get<6>(localRecords.at(i));
        require(a.kind==b.kind && a.originalVertex==b.originalVertex && a.incomingArc==b.incomingArc,
                "2D and 3D provenance disagree");
    }
}

// No polygon segment remains inside a fully enclosed cross-section.
void exerciseEnclosure(model::DefectiveCrystalBase<3>& base,bool reverse)
{
    using V=Eigen::Vector3d;
    const double length=base.mesh.xMax().x()-base.mesh.xMin().x();
    const auto grain=base.poly.grains.begin()->second;
    const auto plane=base.dislocationDynamics().periodicGlidePlaneFactory.get(model::GlidePlaneKey<3>(V::Zero(),grain->reciprocalLatticeDirection(V::UnitZ())));
    for(int winding : {1,2})
    {
        std::vector<V> vertices;
        for(int layer=0; layer<winding; ++layer)
        {
            for(const V& v : {V(-.8,-.8,0),V(.8,-.8,0),V(.8,.8,0),V(-.8,.8,0)})
            {
                vertices.push_back(length*v);
            }
        }
        if(reverse)
        {
            std::reverse(vertices.begin(),vertices.end());
        }
        std::vector<std::pair<V,const OriginalNode* const>> polygon;
        for(const auto& v : vertices)
        {
            polygon.emplace_back(v,nullptr);
        }
        const auto records=plane->polygonPatchIntersection(polygon,false);
        require(records.size()>=4,"Enclosed grain cross-section was discarded");
        double area=0;
        for(size_t i=0; i<records.size(); ++i)
        {
            const auto& record=records[i];
            const V a=plane->referencePlane->globalPosition(std::get<0>(record));
            const V b=plane->referencePlane->globalPosition(std::get<0>(records[(i+1)%records.size()]));
            area+=0.5*a.cross(b).z();
            const auto& info=std::get<6>(record);
            require(info.onPhysicalBoundary && info.atPatchCorner && !info.originalVertex && !info.polygonLocation,
                    "Enclosed boundary has incorrect provenance");
            require(info.incomingArc==model::PolygonPatchArcKind::PhysicalBoundary,
                    "Enclosed boundary falsely attributed to a polygon segment");
        }
        require(std::abs(area/(length*length)-(reverse ? -winding : winding))<1.e-8,
                "Enclosed boundary lost its orientation or winding multiplicity");
        for(auto& v : polygon)
        {
            v.first.x()+=3*length;
        }
        require(plane->polygonPatchIntersection(polygon,false).empty(),"Disjoint polygon retained a slipped area");
    }
}

void exercisePeriodicEnclosure(model::DefectiveCrystalBase<3>& base,bool reverse)
{
    using V=Eigen::Vector3d;
    const double length=base.mesh.xMax().x()-base.mesh.xMin().x();
    const auto grain=base.poly.grains.begin()->second;
    const auto plane=base.dislocationDynamics().periodicGlidePlaneFactory.get(model::GlidePlaneKey<3>(V::Zero(),grain->reciprocalLatticeDirection(V(1,0,3))));
    std::vector<V> vertices={V(-2,-.8,2.0/3),V(2,-.8,-2.0/3),V(2,.8,-2.0/3),V(-2,.8,2.0/3)};
    if(reverse)
    {
        std::reverse(vertices.begin(),vertices.end());
    }
    std::vector<std::pair<V,const OriginalNode* const>> polygon;
    for(const auto& v : vertices)
    {
        polygon.emplace_back(length*v,nullptr);
    }
    const auto records=plane->polygonPatchIntersection(polygon,false);
    require(records.size()>=4,"Enclosed periodic component was discarded");
    V area=V::Zero();
    size_t seams=0;
    for(size_t i=0; i<records.size(); ++i)
    {
        const auto& record=records[i];
        const V a=plane->referencePlane->globalPosition(std::get<0>(record));
        const V b=plane->referencePlane->globalPosition(std::get<0>(records[(i+1)%records.size()]));
        area+=0.5*a.cross(b);
        const auto& info=std::get<6>(record);
        require(info.onPhysicalBoundary && !info.originalVertex && !info.polygonLocation,
                "Enclosed periodic boundary has incorrect provenance");
        seams+=info.crossesPeriodicSeam;
        const V mapped=a+std::get<1>(record);
        require((mapped.array()>=base.mesh.xMin().array()-1.e-6).all() &&
                (mapped.array()<=base.mesh.xMax().array()+1.e-6).all(),"Enclosed boundary mapped outside the RVE");
    }
    require(seams>0,"Enclosed boundary did not cross same-grain periodic twins");
    const V expected=(reverse ? -1.0 : 1.0)*length*length*V(1,0,3);
    require((area-expected).norm()<1.e-8*expected.norm(),"Wrong enclosed periodic area or orientation");
}

// A convex hexagonal prism whose X-periodic copies have a zigzag surface.
// Fan triangulation about each face center ensures matching opposite faces.
void writeHexagonalPrism(const std::filesystem::path& file)
{
    const std::vector<Eigen::Vector2d> hex={{0,.5},{-.5,.25},{-.5,-.25},{0,-.5},{.5,-.25},{.5,.25}};
    std::vector<Eigen::Vector3d> nodes;
    for(double z : {-.5,.5})
    {
        for(const auto& v : hex)
        {
            nodes.emplace_back(v.x(),v.y(),z);
        }
    }
    std::vector<std::vector<int>> faces={{0,1,2,3,4,5},{6,7,8,9,10,11}};
    for(int i=0; i<6; ++i)
    {
        faces.push_back({i,(i+1)%6,(i+1)%6+6,i+6});
    }
    const int center=nodes.size();
    nodes.emplace_back(0,0,0);
    std::vector<std::array<int,4>> tetrahedra;
    for(const auto& face : faces)
    {
        Eigen::Vector3d mean=Eigen::Vector3d::Zero();
        for(int i : face)
        {
            mean+=nodes[i];
        }
        mean/=face.size();
        const int faceCenter=nodes.size();
        nodes.push_back(mean);
        for(size_t i=0; i<face.size(); ++i)
        {
            tetrahedra.push_back({center,faceCenter,face[i],face[(i+1)%face.size()]});
        }
    }
    std::ofstream out(file);
    out<<"$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n"<<nodes.size()<<'\n';
    for(size_t i=0; i<nodes.size(); ++i)
    {
        out<<i+1<<' '<<nodes[i].transpose()<<'\n';
    }
    out<<"$EndNodes\n$Elements\n"<<tetrahedra.size()<<'\n';
    for(size_t i=0; i<tetrahedra.size(); ++i)
    {
        out<<i+1<<" 4 2 1 1";
        for(int j : tetrahedra[i])
        {
            out<<' '<<j+1;
        }
        out<<'\n';
    }
    out<<"$EndElements\n";
}

int main(int argc,char** argv)
{
    try
    {
        require(argc==2,"Expected repository root");
        const std::filesystem::path root(argv[1]);
        const auto folder=std::filesystem::temp_directory_path()/
            ("periodic-walk-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        const auto input=folder/"inputFiles";
        std::filesystem::create_directories(input);
        std::filesystem::create_directories(folder/"F");
        std::filesystem::copy_file(root/"Library/Materials/W.txt",input/"W.txt");
        std::filesystem::copy_file(root/"Library/Meshes/unitCube24.msh",input/"mesh.msh");
        int seed=-1;
        {
            model::SimplicialMesh<3> mesh((input/"mesh.msh").string(),Eigen::Matrix3d::Identity(),Eigen::Vector3d::Zero(),{-1});
            for(const auto& region:mesh.regions())
            {
                for(const auto& face:region.second->faces())
                {
                    if(face.second->periodicFacePair.first.x()>0.9)
                    {
                        seed=face.first;
                    }
                }
            }
        }
        require(seed>=0,"No X face");
        for(int rank : {0,1,3})
        {
            const bool periodic=rank>0;
            {
                std::ofstream f(input/"DefectiveCrystal.txt");
                f<<"physics=DislocationDynamics;\nfmmTreeDepth=0;\nuseFEM=0;\nNsteps=0;\nmaxResolveSteps=0;\ndtMax=1;\noutputFrequency=1;\noutputBinary=0;\nstartAtTimeStep=0;\nEwaldLengthFactor=0;\n";
                f<<"periodicImageSize="<<(rank==3 ? "0 0 0" : "0")<<";\n";
                std::ofstream p(input/"polycrystal.txt");
                const double length=200/std::sqrt(3.0);
                p<<std::setprecision(17)<<"materialFile=W.txt;\nmeshFile=mesh.msh;\nabsoluteTemperature=300;\nC2G1=1 0 0 0 1 0 0 0 1;\nF="
                 <<length<<" 0 0 0 "<<length<<" 0 0 0 "<<length<<";\nX0=0 0 0;\nperiodicFaceIDs=";
                if(periodic)
                {
                    p<<(rank==3 ? -1 : seed);
                }
                p<<";\n";
            }
            model::DefectiveCrystalBase<3> base(folder.string());
            for(bool reverse:{false,true})
            {
                exercise(base,{{-.2,-.2},{.2,-.2},{.2,.2},{-.2,.2}},.16,reverse);
                std::vector<Eigen::Vector2d> balanced{{-.2,-.2},{.2,.2},{-.2,.2},{.2,-.2}};
                for(size_t start=0;start<balanced.size();++start)
                {
                    exercise(base,balanced,0.0,reverse);
                    auto clippedBalanced=balanced;
                    for(auto& vertex : clippedBalanced)
                    {
                        vertex*=4.0;
                    }
                    exercise(base,clippedBalanced,0.0,reverse);
                    std::rotate(balanced.begin(),balanced.begin()+1,balanced.end());
                }
                // Tangency at a corner must not start a spurious boundary walk.
                exercise(base,{{.1,.3},{.5,.5},{.1,.4}},.02,reverse);
                if(rank==3)
                {
                    exercise(base,{{-.2,.2},{1.2,.2},{1.2,.8},{-.2,.8}},.84,reverse);
                    exercise(base,{{.1,.1},{.9,.9},{-.1,.3}},.16,reverse);
                    continue;
                }
                // Two overlapping continuations meet the physical face at
                // the same point. Keep the next branch in contour order,
                // rather than jumping to the later overlapping return branch.
                for(double epsilon : {1.e-3,1.e-4,1.e-5})
                {
                    std::vector<Eigen::Vector2d> overlap={{0,.5},{.1,.4},{0,.4},
                        {-.1*epsilon,.5+.1*epsilon},{.3,.2},{-.3,.2},{-.3,.5}};
                    // Exact area of the contour clipped by y<=0.5.
                    const double expected=-.14-.005*epsilon/(1+epsilon);
                    for(size_t start=0; start<overlap.size(); ++start)
                    {
                        exercise(base,overlap,expected,reverse);
                        std::rotate(overlap.begin(),overlap.begin()+1,overlap.end());
                    }
                }
                // A later boundary-following return passes an earlier inward
                // branch. Follow the return itself instead of cycling forever.
                std::vector<Eigen::Vector2d> surfaceReturn={{0,.4},{.1,.3},{.1,.5001},
                    {.08,.5},{.05,.4},{.3,.3},{.3,.5001},{.28,.5},{0,.5}};
                for(size_t start=0; start<surfaceReturn.size(); ++start)
                {
                    exercise(base,surfaceReturn,.051,reverse);
                    std::rotate(surfaceReturn.begin(),surfaceReturn.begin()+1,surfaceReturn.end());
                }
                // Physical face exit and re-entry, including an outside first vertex.
                exercise(base,{{-.2,.8},{-.2,.2},{.2,.2},{.2,.8}},.12,reverse);
                exercise(base,{{-.2,.2},{.2,.2},{.2,.5},{-.2,.5}},.12,reverse);
                if(periodic)
                {
                    exercisePeriodicEnclosure(base,reverse);
                    // A self-intersection gives the top lobe the opposite
                    // winding to the overall contour. Its physical-boundary
                    // connection must follow predecessors on the open chain.
                    // Translate it across a seam to exercise twin crossings
                    // in that reverse walk, and vary the starting vertex.
                    for(double offset : {0.0,0.4,0.5,0.6,1.5})
                    {
                        std::vector<Eigen::Vector2d> folded={{-.4,0},{.2,.6},{-.2,.6},{.4,0}};
                        for(auto& point : folded)
                        {
                            point.x()+=offset;
                        }
                        for(size_t start=0; start<folded.size(); ++start)
                        {
                            exercise(base,folded,-.15,reverse);
                            std::rotate(folded.begin(),folded.begin()+1,folded.end());
                        }
                    }
                    // Original polygon crosses same-grain seams without clipping.
                    exercise(base,{{-.2,-.2},{1.2,-.2},{1.2,.2},{-.2,.2}},.56,reverse);
                    // Boundary walk crosses a seam at a physical-face corner.
                    exercise(base,{{-.2,.2},{1.2,.2},{1.2,.8},{-.2,.8}},.42,reverse);
                    // Two surviving components: retain only the starting left one.
                    exercise(base,{{-.4,.1},{-.2,.1},{-.2,.7},{.2,.7},{.2,.1},{.4,.1},{.4,.9},{-.4,.9}},.08,reverse);
                }
                else
                {
                    exerciseEnclosure(base,reverse);
                    // A negative local lobe at a closed grain boundary must
                    // not acquire an extra winding around the entire grain.
                    std::vector<Eigen::Vector2d> folded={{-.4,0},{.2,.6},{-.2,.6},{.4,0}};
                    for(size_t start=0;start<folded.size();++start)
                    {
                        exercise(base,folded,-.15,reverse);
                        std::rotate(folded.begin(),folded.begin()+1,folded.end());
                    }
                    // Same-face ports can connect around the whole boundary,
                    // rather than directly backwards along their shared edge.
                    exercise(base,{{0,.2},{.2,.8},{1,1},{1,-1},{-1,-1},{-1,1},{-.2,.8}},-.97,reverse);
                    // Follow two physical faces through their shared corner.
                    exercise(base,{{0,0},{.8,0},{.8,.8},{0,.8}},.25,reverse);
                }
            }
        }
        writeHexagonalPrism(input/"mesh.msh");
        {
            // Select the vertical X faces directly; the test isolates patch
            // traversal from lattice-basis discovery and its disambiguation.
            std::ofstream p(input/"polycrystal.txt");
            const double length=200/std::sqrt(3.0);
            p<<std::setprecision(17)<<"materialFile=W.txt;\nmeshFile=mesh.msh;\nabsoluteTemperature=300;\nC2G1=1 0 0 0 1 0 0 0 1;\nF="
             <<length<<" 0 0 0 "<<length<<" 0 0 0 "<<length<<";\nX0=0 0 0;\nperiodicFaceIDs=;\n";
        }
        {
            model::DefectiveCrystalBase<3> base(folder.string());
            model::Polyhedron<3,3-1>* left=nullptr;
            model::Polyhedron<3,3-1>* right=nullptr;
            for(const auto& region : base.mesh.regions())
            {
                for(const auto& face : region.second->faces())
                {
                    if(face.second->outNormal().x() < -.99)
                    {
                        left=face.second.get();
                    }
                    if(face.second->outNormal().x() > .99)
                    {
                        right=face.second.get();
                    }
                }
            }
            require(left && right,"Missing prism seam");
            const Eigen::Vector3d shift((base.mesh.xMax()[0]-base.mesh.xMin()[0]),0,0);
            left->periodicFacePair={shift,right};
            right->periodicFacePair={-shift,left};
            for(bool reverse : {false,true})
            {
                exercise(base,{{-.2,.1},{1.2,.1},{1.2,.7},{-.2,.7}},.415,reverse);
                // All vertices lie in a distant copy; straight-line location
                // from the reference center would leave the zigzag domain.
                exercise(base,{{1.95,.4},{2.05,.4},{2.05,.45},{1.95,.45}},.005,reverse);
            }
        }
        // Cross-grain periodic faces are walls for a loop belonging to grain 1.
        std::filesystem::copy_file(root/"Library/Meshes/bicrystal_48.msh",input/"mesh.msh",
                                   std::filesystem::copy_options::overwrite_existing);
        {
            std::ofstream p(input/"polycrystal.txt");
            const double length=200/std::sqrt(3.0);
            p<<std::setprecision(17)<<"materialFile=W.txt;\nmeshFile=mesh.msh;\nabsoluteTemperature=300;\nC2G1=1 0 0 0 1 0 0 0 1;\nC2G2=1 0 0 0 1 0 0 0 1;\nF="
             <<length<<" 0 0 0 "<<length<<" 0 0 0 "<<length<<";\nX0=0 0 0;\nperiodicFaceIDs=-1;\n";
        }
        {
            std::ofstream f(input/"DefectiveCrystal.txt");
            f<<"physics=DislocationDynamics;\nfmmTreeDepth=0;\nuseFEM=0;\nNsteps=0;\nmaxResolveSteps=0;\ndtMax=1;\noutputFrequency=1;\noutputBinary=0;\nstartAtTimeStep=0;\nperiodicImageSize=0 0 0;\nEwaldLengthFactor=0;\n";
        }
        {
            model::DefectiveCrystalBase<3> base(folder.string());
            for(bool reverse : {false,true})
            {
                exercise(base,{{-.4,-.1},{-.2,-.1},{-.2,.1},{-.4,.1}},.04,reverse);
                exercise(base,{{-.4,-.1},{-.6,-.1},{-.6,.1},{-.4,.1}},-.02,reverse);
            }
        }
        std::filesystem::remove_all(folder);
        std::cout<<"PASS: periodic polygon boundary walking\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr<<error.what()<<'\n';
        return 1;
    }
}
