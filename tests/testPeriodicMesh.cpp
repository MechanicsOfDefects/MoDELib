#include <SimplicialMesh.h>
#include <PeriodicMeshFacePairing.h>
#include <Eigen/Geometry>
#include <iostream>
#include <stdexcept>

void require(bool condition, const std::string& message)
{
    if(!condition)
    {
        throw std::runtime_error(message);
    }
}

void check(const std::string& path, int expected, const Eigen::Matrix3d& rotation)
{
    model::SimplicialMesh<3> mesh(path,rotation,Eigen::Vector3d::Zero(),{-1});
    const auto basis = mesh.periodicBasis();
    require(basis.cols() == expected,"Wrong detected periodic rank: " + path);
    std::vector<model::PlanarMeshFace<3>*> faces;
    for(const auto& region : mesh.regions())
    {
        for(const auto& entry : region.second->faces())
        {
            if(entry.second->isExternal())
            {
                faces.push_back(entry.second.get());
            }
        }
    }
    const double tol = 1.e-7*(mesh.xMax()-mesh.xMin()).norm();
    std::set<int> seeds;
    size_t pairedFaces=0;
    size_t diagonalFaces=0;
    for(const auto* face : faces)
    {
        const auto& pair=face->periodicFacePair;
        if(!pair.second)
        {
            continue;
        }
        ++pairedFaces;
        require(pair.second->periodicFacePair.second==face,"Nonreciprocal pairing");
        require((pair.first+pair.second->periodicFacePair.first).norm()<tol,"Nonreciprocal translation");
        Eigen::Vector3d shift;
        require(model::isPeriodicFacePair(*face,*pair.second,shift,tol),"Invalid triangle pairing");
        require(model::periodicLatticeContains<3>(basis,shift,tol),"Pair is outside the integer lattice");
        const Eigen::Vector3d local=rotation.transpose()*shift;
        require(std::abs(local.z())<tol,"Unexpected Z periodicity");
        if(std::abs(local.x())>tol && std::abs(local.y())>tol)
        {
            ++diagonalFaces;
        }
    }
    require(pairedFaces==size_t(expected==1?24:50),"Wrong number of paired faces");
    require(diagonalFaces==size_t(expected==1?0:12),"Missing X+Y or X-Y diagonal pairs");
    for(int axis=0; axis<expected; ++axis)
    {
        require(model::periodicLatticeContains<3>(basis,Eigen::Vector3d(rotation.col(axis)),tol),
                "Detected basis does not generate an expected unit translation");
        const auto seed=std::find_if(faces.begin(),faces.end(),[&](const auto* face)
        {
            return face->periodicFacePair.second
                && face->periodicFacePair.first.normalized().cross(rotation.col(axis)).norm()<1.e-6;
        });
        require(seed!=faces.end(),"Missing axis seed face");
        seeds.insert((*seed)->sID);
        const int oppositeSeed=(*seed)->periodicFacePair.second->sID;
        if(expected==2)
        {
            require(std::abs((*seed)->outNormal().dot(rotation.col(0)))>1.e-6
                    && std::abs((*seed)->outNormal().dot(rotation.col(1)))>1.e-6,
                    "XY selection fixture must exercise inclined face normals");
        }
        mesh.identifyParallelFaces({static_cast<int>((*seed)->sID)});
        require(mesh.periodicBasis().cols()==1,"Single face enabled extra directions");
        size_t activeFaces=0;
        for(const auto* face : faces)
        {
            if(face->periodicFacePair.second)
            {
                ++activeFaces;
                require(face->periodicFacePair.first.normalized().cross(rotation.col(axis)).norm()<1.e-6,
                        "Single-axis selection enabled another direction or diagonal");
            }
        }
        require(activeFaces==size_t(expected==1?24:(axis==0?20:18)),"Incomplete single-axis activation");
        const auto rediscovered=model::discoverPeriodicMeshLattice(faces,tol);
        require(rediscovered.basis.cols()==expected && rediscovered.pairs.size()==pairedFaces,
                "Active selection changed discovery of the full available lattice");
        mesh.identifyParallelFaces({oppositeSeed});
        require(mesh.periodicBasis().cols()==1
                && mesh.periodicBasis().col(0).normalized().cross(rotation.col(axis)).norm()<1.e-6,
                "Selecting the opposite face changed the active direction");
        mesh.identifyParallelFaces({-1});

    }
    // Explicit X/Y seeds must activate the same complete pairing, including diagonals.
    std::map<size_t,size_t> original;
    for(const auto* face : faces)
    {
        if(face->periodicFacePair.second)
        {
            original.emplace(face->sID,face->periodicFacePair.second->sID);
        }
    }
    mesh.identifyParallelFaces(seeds);
    require(mesh.periodicBasis().cols()==expected,"Wrong explicitly selected rank");
    for(const auto* face : faces)
    {
        require(bool(face->periodicFacePair.second)==bool(original.count(face->sID)),"Selected face coverage changed");
        if(face->periodicFacePair.second)
        {
            require(original.at(face->sID)==face->periodicFacePair.second->sID,"Selected partner changed");
        }
    }
    if(expected==2)
    {
        std::set<int> diagonalSeeds;
        for(const int sign : {-1,1})
        {
            const Eigen::Vector3d diagonal=rotation*Eigen::Vector3d(1,sign,0).normalized();
            const auto face=std::find_if(faces.begin(),faces.end(),[&](const auto* item)
            {
                return item->periodicFacePair.second
                    && item->periodicFacePair.first.normalized().cross(diagonal).norm()<1.e-6;
            });
            require(face!=faces.end(),"Missing diagonal face pair");
            diagonalSeeds.insert((*face)->sID);
            mesh.identifyParallelFaces({static_cast<int>((*face)->sID)});
            require(mesh.periodicBasis().cols()==1,"One diagonal face enabled both axes");
            size_t activeCount=0;
            for(const auto* item : faces)
            {
                if(item->periodicFacePair.second)
                {
                    ++activeCount;
                    require(item->periodicFacePair.first.normalized().cross(diagonal).norm()<1.e-6,
                            "One diagonal face activated a different translation");
                }
            }
            require(activeCount==size_t(sign==1?8:4),"Wrong single-diagonal pair count");
            mesh.identifyParallelFaces(seeds);
        }
        // Selecting diagonal generators activates their index-two sublattice,
        // without implicitly enabling the X or Y translations.
        mesh.identifyParallelFaces(diagonalSeeds);
        require(mesh.periodicBasis().cols()==2,"Wrong diagonal sublattice rank");
        size_t activeDiagonals=0;
        for(const auto* face : faces)
        {
            if(face->periodicFacePair.second)
            {
                const Eigen::Vector3d local=rotation.transpose()*face->periodicFacePair.first;
                require(std::abs(local.x())>tol && std::abs(local.y())>tol,"Diagonal selection enabled X or Y");
                ++activeDiagonals;
            }
        }
        require(activeDiagonals==12,"Missing diagonal sublattice pairs");
        diagonalSeeds.insert(*seeds.begin());
        mesh.identifyParallelFaces(diagonalSeeds);
        for(const auto* face : faces)
        {
            if(original.count(face->sID))
            {
                require(face->periodicFacePair.second
                        && original.at(face->sID)==face->periodicFacePair.second->sID,
                        "Alternative generators lost a periodic pair");
            }
        }
        mesh.identifyParallelFaces(seeds);
    }
    // A locally matching pair is insufficient when another required face has no partner.
    auto incomplete=faces;
    const auto removed=std::find_if(incomplete.begin(),incomplete.end(),[](const auto* face)
    {
        return face->periodicFacePair.second!=nullptr;
    });
    require(removed!=incomplete.end(),"Missing test face");
    incomplete.erase(removed);
    bool incompleteRejected=false;
    try
    {
        model::discoverPeriodicMeshLattice(incomplete,tol);
    }
    catch(const std::exception&)
    {
        incompleteRejected=true;
    }
    require(incompleteRejected,"Accepted lattice despite a missing required face partner");
    // The hull is unchanged when a contained triangle is removed; pairing must still fail.
    for(auto* face : faces)
    {
        if(face->periodicFacePair.second)
        {
            const auto* triangle = *face->simplices().begin();
            face->simplices().erase(triangle);
            Eigen::Vector3d shift;
            const bool matched = model::isPeriodicFacePair(*face,*face->periodicFacePair.second,shift,tol);
            face->simplices().insert(triangle);
            require(!matched,"Accepted differing triangulations with matching hulls");
            break;
        }
    }
    // Invalid selections must not change the current pairing.
    const auto oldBasis=mesh.periodicBasis();
    bool rejected=false;
    try
    {
        mesh.identifyParallelFaces({2147483647});
    }
    catch(const std::exception&)
    {
        rejected=true;
    }
    require(rejected,"Unknown face ID accepted");
    require((oldBasis-mesh.periodicBasis()).norm()<tol,"Failed selection changed mesh");
    for(const auto* face : faces)
    {
        if(!face->periodicFacePair.second)
        {
            rejected=false;
            try
            {
                mesh.identifyParallelFaces({static_cast<int>(face->sID)});
            }
            catch(const std::exception&)
            {
                rejected=true;
            }
            require(rejected,"Nonperiodic face accepted");
            break;
        }
    }
    mesh.identifyParallelFaces({});
    require(mesh.periodicBasis().cols()==0,"Empty selection did not clear periodicity");
}

void checkBicrystal(const std::string& path,const Eigen::Matrix3d& rotation)
{
    model::SimplicialMesh<3> mesh(path,rotation,Eigen::Vector3d::Zero(),{-1});
    const auto basis=mesh.periodicBasis();
    require(basis.cols()==3,"Bicrystal must have three periodic directions");
    std::set<int> axes;
    for(int k=0; k<3; ++k)
    {
        const Eigen::Vector3d local=rotation.transpose()*basis.col(k);
        Eigen::Index axis;
        local.cwiseAbs().maxCoeff(&axis);
        Eigen::Vector3d expected=Eigen::Vector3d::Zero();
        expected[axis]=local[axis]>0.0 ? (axis==0?2.0:1.0) : -(axis==0?2.0:1.0);
        require((local-expected).norm()<1.e-7,"Bicrystal basis is not box aligned");
        axes.insert(static_cast<int>(axis));
    }
    require(axes.size()==3,"Repeated bicrystal direction");
    size_t pairedFaces=0;
    for(const auto& region : mesh.regions())
    {
        for(const auto& entry : region.second->faces())
        {
            const auto& face=*entry.second;
            if(face.periodicFacePair.second)
            {
                ++pairedFaces;
                require(face.isExternal() && face.periodicFacePair.second->isExternal(),"Internal face paired");
                require(face.periodicFacePair.second->periodicFacePair.second==&face,"Nonreciprocal bicrystal pairing");
                const Eigen::Vector3d shift=rotation.transpose()*face.periodicFacePair.first;
                require((shift.cwiseAbs().array()>1.e-7).count()==1,"Diagonal bicrystal pairing survived");
                require(face.regionIDs==face.periodicFacePair.second->regionIDs || std::abs(shift.x())>1.e-7,
                        "Y/Z pairing crosses grains");
            }
        }
    }
    require(pairedFaces==10,"Expected all five bicrystal face pairs");
}

int main(int argc,char** argv)
{
    try
    {
        require(argc==2,"Repository path required");
        const std::string root=argv[1];
        const std::string meshes=root+"/Library/Meshes/";
        const Eigen::Matrix3d rotation=Eigen::AngleAxisd(0.71,Eigen::Vector3d(1,2,3).normalized()).toRotationMatrix();
        check(meshes+"polycrystalCube_10_periodicX.msh",1,Eigen::Matrix3d::Identity());
        check(meshes+"polycrystalCube_10_periodicX.msh",1,rotation);
        check(meshes+"polycrystalCube_10_periodicXY.msh",2,Eigen::Matrix3d::Identity());
        check(meshes+"polycrystalCube_10_periodicXY.msh",2,rotation);
        checkBicrystal(meshes+"bicrystal_48.msh",Eigen::Matrix3d::Identity());
        checkBicrystal(meshes+"bicrystal_48.msh",rotation);
        // Membership must reject fractional coefficients and out-of-plane shifts.
        Eigen::Matrix<double,3,Eigen::Dynamic> xy(3,2);
        xy.col(0)=Eigen::Vector3d::UnitX();
        xy.col(1)=Eigen::Vector3d::UnitY();
        require(model::periodicLatticeContains<3>(xy,Eigen::Vector3d(1,-1,0),1.e-8),"Diagonal membership failed");
        require(!model::periodicLatticeContains<3>(xy,Eigen::Vector3d(0.5,0,0),1.e-8),"Fractional translation accepted");
        require(!model::periodicLatticeContains<3>(xy,Eigen::Vector3d(1,0,1),1.e-8),"Out-of-span translation accepted");
        std::cout<<"Periodic mesh tests passed\n";
    }
    catch(const std::exception& error)
    {
        std::cerr<<error.what()<<std::endl;
        return 1;
    }
    return 0;
}
