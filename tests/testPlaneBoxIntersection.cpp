#include <FMM/PlaneBoxIntersection.h>
#include <iostream>
#include <random>
#include <stdexcept>
#include <set>
#include <limits>

using Cut = model::PlaneBoxIntersection<3>;
using Tree = model::BoxTree<3>;
void require(bool condition, const char* message)
{
    if(!condition)
    {
        throw std::runtime_error(message);
    }
}

double check(const Tree& tree, const Cut& cut)
{
    double area = 0;
    std::map<const Cut::Segment*, std::vector<int>> uses;
    std::map<Cut::SegmentKey, const Cut::Segment*> identities;
    for(const auto& [key, polygon] : cut.polygons())
    {
        require(polygon.vertices.size() >= 3 && polygon.vertices.size() == polygon.boundary.size(), "Invalid polygon size");
        require(polygon.area > 0, "Nonpositive area");
        Eigen::Vector3d orientedArea = Eigen::Vector3d::Zero();
        const auto origin = polygon.vertices[0]->position;
        for(size_t i = 0; i < polygon.vertices.size(); ++i)
        {
            const auto a = polygon.vertices[i];
            const auto b = polygon.vertices[(i + 1) % polygon.vertices.size()];
            const auto side = polygon.boundary[i];
            require((side.orientation == 1 ? side.segment->first : side.segment->second) == a &&
                    (side.orientation == 1 ? side.segment->second : side.segment->first) == b,
                    "Contour connectivity does not follow shared segments");
            require(a != b && (a->position - b->position).norm() > 0, "Zero-length contour side");
            require(std::abs(cut.normal().dot(a->position - cut.origin())) < 1e-10, "Point is off plane");
            require((a->position.array() >= polygon.box->lower().array() - 1e-10).all() &&
                    (a->position.array() <= polygon.box->upper().array() + 1e-10).all(), "Point outside cell");
            orientedArea += (a->position - origin).cross(b->position - origin) * .5;
            uses[side.segment.get()].push_back(side.orientation);
            const Cut::SegmentKey segmentKey(side.segment->first->support->key(), side.segment->second->support->key());
            auto [it, inserted] = identities.emplace(segmentKey, side.segment.get());
            require(inserted || it->second == side.segment.get(), "Same segment has multiple identities");
        }
        require(std::abs(cut.normal().dot(orientedArea) - polygon.area) < 1e-10, "Wrong polygon orientation/area");
        area += polygon.area;
    }
    for(const auto& [segment, signs] : uses)
    {
        require(signs.size() <= 2, "Segment used more than twice");
        if(signs.size() == 2)
        {
            require(signs[0] + signs[1] == 0, "Shared edges do not cancel");
        }
        else
        {
            bool outer = false;
            for(int axis = 0; axis < 3; ++axis)
            {
                for(double boundary : {tree.lower()[axis], tree.upper()[axis]})
                {
                    outer = outer || (std::abs(segment->first->position[axis] - boundary) < 1e-10 &&
                                      std::abs(segment->second->position[axis] - boundary) < 1e-10);
                }
            }
            require(outer, "Unpaired internal cut edge");
        }
    }
    return area;
}

void checkLevels(const Tree& tree, const Eigen::Vector3d& p, const Eigen::Vector3d& n)
{
    const Cut root(tree, 0, p, n);
    const double expected = check(tree, root);
    for(unsigned level = 1; level <= tree.depth(); ++level)
    {
        const Cut cut(tree, level, p, n);
        const double area = check(tree, cut);
        require(std::abs(area - expected) < 1e-10 * std::max(1.0, expected), "Area not conserved across levels");
        // Every child polygon belongs to the cut of its geometric parent.
        const Cut parent(tree, level - 1, p, n);
        std::map<Cut::Key, double> childAreas;
        for(const auto& [key, polygon] : cut.polygons())
        {
            childAreas[polygon.box->refinementParent()->key()] += polygon.area;
        }
        for(const auto& [key, polygon] : parent.polygons())
        {
            require(std::abs(childAreas[key] - polygon.area) < 1e-10, "Child areas do not recover parent area");
        }
    }
}

int main()
{
    try
    {
        Tree tree(Eigen::Vector3d::Zero(), Eigen::Vector3d::Ones(), 3);
        for(int axis = 0; axis < 3; ++axis)
        {
            for(double position : {0.0, 0.125, 0.5, 0.743, 1.0, -0.1, 1.1})
            {
                Eigen::Vector3d p = Eigen::Vector3d::Zero();
                p[axis] = position;
                for(double sign : {-1.0, 1.0})
                {
                    checkLevels(tree, p, sign * Eigen::Vector3d::Unit(axis));
                }
            }
        }
        for(const Eigen::Vector3d& n : {Eigen::Vector3d(1,1,1), Eigen::Vector3d(1,1,0), Eigen::Vector3d(1,-1,0), Eigen::Vector3d(1,2,3)})
        {
            for(const Eigen::Vector3d& p : {Eigen::Vector3d::Zero().eval(), Eigen::Vector3d::Constant(.5).eval(), Eigen::Vector3d::Ones().eval()})
            {
                checkLevels(tree, p, n);
                checkLevels(tree, p, -n);
            }
        }
        Cut coincident(tree, 3, Eigen::Vector3d(0,0,.5), Eigen::Vector3d::UnitZ());
        require(coincident.polygons().size() == 64, "Coincident shared face counted twice");
        require(std::abs(check(tree, coincident) - 1) < 1e-12, "Wrong analytic square area");
        Cut diagonal(tree, 3, Eigen::Vector3d::Zero(), Eigen::Vector3d(1,-1,0));
        require(std::abs(check(tree, diagonal) - std::sqrt(2.0)) < 1e-12, "Wrong analytic diagonal area");
        Cut triangle(tree, 3, Eigen::Vector3d(1,0,0), Eigen::Vector3d::Ones());
        require(std::abs(check(tree, triangle) - std::sqrt(3.0)/2) < 1e-12, "Wrong analytic triangle area");
        std::mt19937 generator(42);
        std::uniform_real_distribution<double> random(-1, 1);
        for(int trial = 0; trial < 50; ++trial)
        {
            const Eigen::Vector3d normal(random(generator), random(generator), random(generator));
            const Eigen::Vector3d point(.5 + random(generator), .5 + random(generator), .5 + random(generator));
            checkLevels(tree, point, normal);
        }
        // Near-coincidence, anisotropic boxes and translated coordinates.
        checkLevels(tree, Eigen::Vector3d(0,0,.5 + 1e-10), Eigen::Vector3d::UnitZ());
        Tree stretched(Eigen::Vector3d(-3,2,5), Eigen::Vector3d(4,3,9), 2);
        checkLevels(stretched, Eigen::Vector3d(0,2.5,6), Eigen::Vector3d(2,-3,1));
        bool rejected = false;
        try
        {
            Cut invalid(tree, 0, Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero());
        }
        catch(const std::invalid_argument&)
        {
            rejected = true;
        }
        require(rejected, "Invalid normal accepted");
        const auto oldArea = check(tree, coincident);
        tree.setDepth(1);
        require(coincident.revision() != tree.revision() && std::abs(check(tree, coincident) - oldArea) < 1e-12,
                "Retained intersection snapshot changed after tree rebuild");
        std::cout << "Plane/box topology, ownership, orientation and area conservation passed\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
