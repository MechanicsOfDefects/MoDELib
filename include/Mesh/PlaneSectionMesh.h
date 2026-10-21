#ifndef model_PlaneSectionMesh_H_
#define model_PlaneSectionMesh_H_

#include <SimplicialMesh.h>
#include <Plane.h>
#include <TriangularMesh.h>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace model
{
    // Convex grains may form a non-convex or disconnected overall section.
    // Triangulate their sections independently rather than ordering the whole
    // external perimeter by polar angle. Coordinates are local to plane.
    inline TriangularMesh planeSectionMesh(const SimplicialMesh<3>& mesh,
                                          const Plane<3>& plane,
                                          const double maxArea)
    {
        if(!std::isfinite(maxArea) || maxArea<=0.0)
        {
            throw std::invalid_argument("Plane mesh area must be finite and positive");
        }
        const double scale=(mesh.xMax()-mesh.xMin()).norm();
        const double tolerance=1.e-10*scale;
        TriangularMesh result;
        std::vector<std::deque<Eigen::Vector2d>> sections;
        for(const auto& region : mesh.regions())
        {
            std::vector<Eigen::Vector2d> points;
            const auto addPoint=[&](const Eigen::Vector3d& point)
            {
                const Eigen::Vector2d local=plane.localPosition(point);
                for(const auto& existing : points)
                {
                    if((existing-local).norm()<=tolerance)
                    {
                        return;
                    }
                }
                points.push_back(local);
            };
            for(const auto& entry : region.second->faces())
            {
                const auto& hull=entry.second->convexHull();
                for(size_t k=0;k<hull.size();++k)
                {
                    const auto& a=hull[k]->P0;
                    const auto& b=hull[(k+1)%hull.size()]->P0;
                    const double da=(a-plane.P).dot(plane.unitNormal);
                    const double db=(b-plane.P).dot(plane.unitNormal);
                    if(std::abs(da)<=tolerance)
                    {
                        addPoint(a);
                    }
                    if((da>tolerance && db< -tolerance) || (da< -tolerance && db>tolerance))
                    {
                        addPoint(a+(da/(da-db))*(b-a));
                    }
                }
            }
            if(points.size()<3)
            {
                continue;
            }
            std::sort(points.begin(),points.end(),[](const auto& a,const auto& b)
            {
                return a.x()<b.x() || (a.x()==b.x() && a.y()<b.y());
            });
            const auto cross=[](const Eigen::Vector2d& a,const Eigen::Vector2d& b,const Eigen::Vector2d& c)
            {
                return (b.x()-a.x())*(c.y()-a.y())-(b.y()-a.y())*(c.x()-a.x());
            };
            std::deque<Eigen::Vector2d> hull;
            for(const auto& point : points)
            {
                while(hull.size()>1 && cross(hull[hull.size()-2],hull.back(),point)<=tolerance*tolerance)
                {
                    hull.pop_back();
                }
                hull.push_back(point);
            }
            const auto lowerSize=hull.size();
            for(auto it=points.rbegin()+1;it!=points.rend();++it)
            {
                while(hull.size()>lowerSize && cross(hull[hull.size()-2],hull.back(),*it)<=tolerance*tolerance)
                {
                    hull.pop_back();
                }
                hull.push_back(*it);
            }
            hull.pop_back();
            if(hull.size()<3)
            {
                continue;
            }
            // A cut coincident with an internal face belongs to both grains.
            // Keep just one copy to avoid overlapping rendered triangles.
            bool duplicate=false;
            for(const auto& previous : sections)
            {
                if(previous.size()!=hull.size())
                {
                    continue;
                }
                bool same=true;
                for(const auto& point : hull)
                {
                    same=same && std::any_of(previous.begin(),previous.end(),[&](const auto& other)
                    {
                        return (point-other).norm()<=tolerance;
                    });
                }
                duplicate=duplicate || same;
            }
            if(duplicate)
            {
                continue;
            }
            sections.push_back(hull);
            TriangularMesh part;
            part.reMesh(hull,std::deque<Eigen::Vector2d>(),maxArea,"pazqQ");
            const int offset=result.vertices().size();
            result.vertices().insert(result.vertices().end(),part.vertices().begin(),part.vertices().end());
            for(const auto& triangle : part.triangles())
            {
                result.triangles().push_back(triangle+Eigen::Vector3i::Constant(offset));
            }
        }
        return result;
    }
}
#endif
