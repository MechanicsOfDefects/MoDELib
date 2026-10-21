#ifndef model_LoopAreaTriangulation_h_
#define model_LoopAreaTriangulation_h_

#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

namespace model
{
    // Rendering tessellation of a directed planar contour, using nonzero winding.
    // Split at every vertex and crossing height. Within each open horizontal slab
    // edges cannot exchange order, so filled intervals are trapezoids. Retraced
    // bridges cancel, and self-intersections do not require a simple polygon.
    inline std::vector<std::array<Eigen::Vector2d,3>> triangulateLoopArea(
        const std::vector<Eigen::Vector2d>& contour)
    {
        using V=Eigen::Vector2d;
        std::vector<std::array<V,3>> triangles;
        if(contour.size()<3)
        {
            return triangles;
        }
        const V origin=contour.front();
        double scale=0.0;
        for(const auto& point : contour)
        {
            scale=std::max(scale,(point-origin).norm());
        }
        if(scale==0.0)
        {
            return triangles;
        }
        const double eps=128.0*std::numeric_limits<double>::epsilon();
        const auto cross=[](const V& a,const V& b)
        {
            return a.x()*b.y()-a.y()*b.x();
        };
        struct Edge
        {
            V a,b;
            double x(double y) const
            {
                return a.x()+(y-a.y())*(b.x()-a.x())/(b.y()-a.y());
            }
        };
        std::vector<Edge> edges;
        std::vector<double> levels;
        for(size_t i=0;i<contour.size();++i)
        {
            const V a=(contour[i]-origin)/scale;
            const V b=(contour[(i+1)%contour.size()]-origin)/scale;
            levels.push_back(a.y());
            if(std::abs(b.y()-a.y())>eps)
            {
                edges.push_back({a,b});
            }
        }
        for(size_t i=0;i<edges.size();++i)
        {
            for(size_t j=i+1;j<edges.size();++j)
            {
                const V u=edges[i].b-edges[i].a;
                const V v=edges[j].b-edges[j].a;
                const V d=edges[j].a-edges[i].a;
                const double denominator=cross(u,v);
                if(std::abs(denominator)>eps*u.norm()*v.norm())
                {
                    const double t=cross(d,v)/denominator;
                    const double s=cross(d,u)/denominator;
                    if(t>0.0 && t<1.0 && s>0.0 && s<1.0)
                    {
                        levels.push_back(edges[i].a.y()+t*u.y());
                    }
                }
            }
        }
        std::sort(levels.begin(),levels.end());
        const auto append=[&](const V& a,const V& b,const V& c)
        {
            if(cross(b-a,c-a)>eps)
            {
                triangles.push_back({origin+scale*a,origin+scale*b,origin+scale*c});
            }
        };
        for(size_t k=1;k<levels.size();++k)
        {
            const double low=levels[k-1];
            const double high=levels[k];
            if(high-low<=eps)
            {
                continue;
            }
            const double mid=(low+high)*0.5;
            std::vector<const Edge*> active;
            for(const auto& edge : edges)
            {
                if(mid>std::min(edge.a.y(),edge.b.y()) && mid<std::max(edge.a.y(),edge.b.y()))
                {
                    active.push_back(&edge);
                }
            }
            std::stable_sort(active.begin(),active.end(),[&](const Edge* a,const Edge* b)
            {
                return a->x(mid)<b->x(mid);
            });
            int winding=0;
            for(size_t i=0;i+1<active.size();++i)
            {
                winding+=active[i]->b.y()>active[i]->a.y()?1:-1;
                if(winding!=0 && active[i+1]->x(mid)-active[i]->x(mid)>eps)
                {
                    const V a(active[i]->x(low),low);
                    const V b(active[i+1]->x(low),low);
                    const V c(active[i+1]->x(high),high);
                    const V d(active[i]->x(high),high);
                    append(a,b,c);
                    append(a,c,d);
                }
            }
        }
        return triangles;
    }
}
#endif
