#include <SutherlandHodgman.h>
#include <LoopAreaTriangulation.h>
#include <Polygon2D.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

using V=Eigen::Vector2d;
void require(bool condition,const char* message)
{
    if(!condition)
    {
        throw std::runtime_error(message);
    }
}
void check(const std::vector<V>& subject,const std::vector<V>& boundary)
{
    const auto clipped=model::SutherlandHodgman::clip(subject,boundary);
    V low=subject.front(),high=low;
    for(const auto& p : subject)
    {
        low=low.cwiseMin(p);
        high=high.cwiseMax(p);
    }
    const double tolerance=1.e-8*std::max(1.0,(high-low).norm());
    const auto bounded=[&](const V& p)
    {
        require(p.allFinite(),"Clipping produced a nonfinite point");
        require((low-p).maxCoeff()<=tolerance && (p-high).maxCoeff()<=tolerance,
                "Clipping extrapolated outside the subject bounds");
        for(size_t k=0;k<boundary.size();++k)
        {
            const V edge=boundary[(k+1)%boundary.size()]-boundary[k];
            const V relative=p-boundary[k];
            require(edge.x()*relative.y()-edge.y()*relative.x()>=-tolerance*edge.norm(),
                    "Clipped point lies outside its grain patch");
        }
    };
    for(const auto& p : clipped)
    {
        bounded(p);
    }
    const auto triangles=model::triangulateLoopArea(clipped);
    for(const auto& t : triangles)
    {
        for(const auto& p : t)
        {
            bounded(p);
        }
    }
    // Check rendered coverage against the original polygon and grain patch,
    // independently of the clipped contour supplied to the triangulator.
    for(int i=0;i<29;++i)
    {
        for(int j=0;j<29;++j)
        {
            const V p=low+(high-low).cwiseProduct(V((i+.371)/29,(j+.613)/29));
            const bool expected=model::Polygon2D::windingNumber(p,subject)!=0 &&
                                model::Polygon2D::windingNumber(p,boundary)!=0;
            size_t hits=0;
            for(const auto& t : triangles)
            {
                hits+=model::Polygon2D::windingNumber(p,std::vector<V>{t[0],t[1],t[2]})!=0;
            }
            require(hits==size_t(expected),"Shaded coverage disagrees with polygon/patch intersection");
        }
    }
}
int main(int argc,char** argv)
{
    try
    {
        require(argc==2,"Expected repository root");
        const std::vector<V> square{{0,0},{1,0},{1,1},{0,1}};
        require(model::SutherlandHodgman::clip(std::vector<V>{{-3,0},{-2,0},{-2,1},{-3,1}},square).empty(),
                "Disjoint clipping should remain empty through subsequent edges");
        require(model::SutherlandHodgman::clip(std::vector<V>{},square).empty(),"Empty subject was not retained");
        check(square,square);
        check({{-.2,.2},{1.2,.2},{1.2,.8},{-.2,.8}},square);
        check({{0,0},{1,1},{0,1},{1,0}},square);
        std::ifstream in(std::filesystem::path(argv[1])/"tests/data/polycrystalClipSpikes.txt");
        size_t cases=0;
        in>>cases;
        require(cases==9,"Missing saved spike contours");
        for(size_t k=0;k<cases;++k)
        {
            size_t n,m;
            in>>n>>m;
            std::vector<V> subject(n),boundary(m);
            for(auto& p : subject)
            {
                in>>p.x()>>p.y();
            }
            for(auto& p : boundary)
            {
                in>>p.x()>>p.y();
            }
            require(bool(in),"Incomplete saved spike contour");
            for(int direction=0;direction<2;++direction)
            {
                for(size_t start=0;start<subject.size();++start)
                {
                    check(subject,boundary);
                    std::rotate(subject.begin(),subject.begin()+1,subject.end());
                }
                std::reverse(subject.begin(),subject.end());
            }
        }
        std::cout<<"PASS: bounded clipping and spike-free loop-area coverage\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr<<error.what()<<'\n';
        return 1;
    }
}
