#include <SimplicialMesh.h>
#include <iostream>
#include <stdexcept>

void require(bool condition,const char* message)
{
    if(!condition)
    {
        throw std::runtime_error(message);
    }
}

int main(int argc,char** argv)
{
    try
    {
        require(argc==2,"Expected repository root");
        model::SimplicialMesh<3> mesh(std::string(argv[1])+"/Library/Meshes/poly30_periodicXYZ.msh",
                                    Eigen::Matrix3d::Identity(),Eigen::Vector3d::Zero(),{});
        const auto* guess=&mesh.simplices().begin()->second;
        // Node 398 in the reported tutorial frame. The original global walk
        // stopped at an external face, despite this point belonging to grain 26.
        const Eigen::Vector3d point=Eigen::Vector3d(-311.5153825721246,816.6921641196241,2568.115086056312)/3666.1742093541243;
        require(mesh.searchRegion(26,point).first,"Regression point is not in grain 26");
        require(mesh.search(point).first,"Global search missed node 398");
        require(mesh.searchWithGuess(point,guess).first,"Guessed search missed node 398");
        require(mesh.searchWithGuess(point,nullptr).first,"Null-guess search missed node 398");
        require(!mesh.searchRegion(1,point).first,"Region-restricted search escaped its grain");
        for(const auto& entry : mesh.simplices())
        {
            const Eigen::Vector3d center=entry.second.bary2pos(Eigen::Vector4d::Constant(0.25));
            const auto found=mesh.searchWithGuess(center,guess);
            require(found.first,"Global search missed a tetrahedron center");
            require(found.second->pos2bary(center).minCoeff()>-1.e-7,"Search returned a non-containing tetrahedron");
        }
        const Eigen::Vector3d outside=mesh.xMax()+Eigen::Vector3d::Ones()*(mesh.xMax()-mesh.xMin()).norm();
        const auto missed=mesh.search(outside);
        require(!missed.first && missed.second && missed.second->isBoundarySimplex(),
                "Outside search must retain a boundary simplex and report failure");
        std::cout<<"PASS: non-convex polycrystal mesh lookup\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr<<error.what()<<'\n';
        return 1;
    }
}
