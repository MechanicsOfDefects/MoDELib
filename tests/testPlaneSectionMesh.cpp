#include <PlaneSectionMesh.h>
#include <MeshPlane.h>
#include <iostream>

void require(bool value,const char* message)
{
    if(!value)
    {
        throw std::runtime_error(message);
    }
}

// Independent area oracle: intersect every tetrahedron, then sum its convex cut.
double tetrahedralArea(const model::SimplicialMesh<3>& mesh,const model::Plane<3>& plane)
{
    double area=0.0;
    for(const auto& entry : mesh.simplices())
    {
        const auto vertices=entry.second->vertexPositionMatrix();
        std::vector<Eigen::Vector2d> points;
        for(int i=0;i<4;++i)
        {
            for(int j=i+1;j<4;++j)
            {
                const double a=(vertices.col(i)-plane.P).dot(plane.unitNormal);
                const double b=(vertices.col(j)-plane.P).dot(plane.unitNormal);
                if(a*b<0.0)
                {
                    points.push_back(plane.localPosition(vertices.col(i)+a/(a-b)*(vertices.col(j)-vertices.col(i))));
                }
            }
        }
        if(points.size()<3)
        {
            continue;
        }
        Eigen::Vector2d center=Eigen::Vector2d::Zero();
        for(const auto& point : points)
        {
            center+=point/double(points.size());
        }
        std::sort(points.begin(),points.end(),[&](const auto& a,const auto& b)
        {
            return std::atan2(a.y()-center.y(),a.x()-center.x())<std::atan2(b.y()-center.y(),b.x()-center.x());
        });
        for(size_t k=0;k<points.size();++k)
        {
            const Eigen::Vector2d a=points[k]-center;
            const Eigen::Vector2d b=points[(k+1)%points.size()]-center;
            area+=0.5*(a.x()*b.y()-a.y()*b.x());
        }
    }
    return area;
}

double check(const model::SimplicialMesh<3>& mesh,const model::Plane<3>& plane,double maxArea)
{
    const auto section=model::planeSectionMesh(mesh,plane,maxArea);
    double area=0.0;
    for(const auto& triangle : section.triangles())
    {
        const Eigen::Vector2d a=section.vertices().at(triangle[0]);
        const Eigen::Vector2d b=section.vertices().at(triangle[1]);
        const Eigen::Vector2d c=section.vertices().at(triangle[2]);
        const double triArea=0.5*((b.x()-a.x())*(c.y()-a.y())-(b.y()-a.y())*(c.x()-a.x()));
        require(triArea>0.0 && triArea<=maxArea*1.00001,"Invalid triangle or mesh area not respected");
        area+=triArea;
        require(mesh.search(plane.globalPosition((a+b+c)/3.0)).first,"Triangle bridges empty space outside mesh");
    }
    return area;
}

int main(int argc,char** argv)
{
    try
    {
        require(argc==2,"Expected repository path");
        const double scale=14693.56435087598;
        model::SimplicialMesh<3> mesh(std::string(argv[1])+"/Library/Meshes/polycrystalCube_10_periodicXY.msh",
                                    scale*Eigen::Matrix3d::Identity(),Eigen::Vector3d::Zero(),{});
        // Reconstruct the reported cut from three of its printed edge endpoints.
        const Eigen::Vector3d a(-4724.304124390713,15395.35594008324,0);
        const Eigen::Vector3d b(6756.931196017786,15784.62596385352,0);
        const Eigen::Vector3d c(-4183.466787360541,11902.59340322532,3411.985017726859);
        const model::Plane<3> plane(a,(b-a).cross(c-a));
        const double expected=tetrahedralArea(mesh,plane);
        for(const double maxArea : {1.e6,2.e5})
        {
            require(std::abs(check(mesh,plane,maxArea)-expected)<1.e-8*expected,"Section area differs from tetrahedron cuts");
        }
        require(check(mesh,model::Plane<3>(mesh.xMax()+Eigen::Vector3d::Ones()*scale,Eigen::Vector3d::UnitZ()),1.e6)==0.0,
                "Outside plane should produce an empty section");
        model::SimplicialMesh<3> cube(std::string(argv[1])+"/Library/Meshes/bicrystal_48.msh",
                                    Eigen::Matrix3d::Identity(),Eigen::Vector3d::Zero(),{});
        const Eigen::Vector3d center=0.5*(cube.xMin()+cube.xMax());
        // Coincident internal GB must be displayed once, not twice.
        const double gbArea=check(cube,model::Plane<3>(center,Eigen::Vector3d::UnitX()),0.1);
        const double boxArea=(cube.xMax().y()-cube.xMin().y())*(cube.xMax().z()-cube.xMin().z());
        require(std::abs(gbArea-boxArea)<1.e-9,"Coincident internal face duplicated or missing");
        require(std::abs(check(cube,model::Plane<3>(cube.xMin(),Eigen::Vector3d::UnitZ()),0.1)
                -(cube.xMax().x()-cube.xMin().x())*(cube.xMax().y()-cube.xMin().y()))<1.e-9,"External coincident face missing");
        require(check(cube,model::Plane<3>(cube.xMin(),Eigen::Vector3d::Ones()),0.1)==0.0,"Vertex tangent should be empty");
        bool rejected=false;
        try
        {
            model::planeSectionMesh(cube,plane,0.0);
        }
        catch(const std::invalid_argument&)
        {
            rejected=true;
        }
        require(rejected,"Zero mesh area accepted");
        std::cout<<"PASS: non-convex field plane, area, refinement, outside and coincident cuts\n";
    }
    catch(const std::exception& error)
    {
        std::cerr<<error.what()<<'\n';
        return 1;
    }
}
