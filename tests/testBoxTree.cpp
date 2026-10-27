#include <FMM/BoxTree.h>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>

void require(bool condition, const char* message)
{
    if(!condition)
    {
        throw std::runtime_error(message);
    }
}

template<int dim>
void checkTopology()
{
    using Tree = model::BoxTree<dim>;
    using Vector = typename Tree::Vector;
    using Key = typename Tree::Key;
    Vector lo;
    Vector hi;
    for(int axis = 0; axis < dim; ++axis)
    {
        lo[axis] = -3.25 + axis;
        hi[axis] = lo[axis] + 2.0 * (axis + 1);
    }
    Tree tree(lo, hi, 2);
    for(unsigned level = 0; level <= tree.depth(); ++level)
    {
        const unsigned n = 1u << level;
        require(tree.boxes(level).size() == (1u << (dim * level)), "Wrong number of boxes");
        double volume = 0.0;
        size_t external = 0;
        for(const auto& box : tree.boxes(level))
        {
            require(box->boundary().size() == 2 * dim, "Missing box faces");
            require(tree.locate(box->center(), level) == box, "Point lookup disagrees with geometry");
            volume += box->measure();
            for(const auto& face : box->boundary())
            {
                const auto neighbor = tree.neighbor(box, face.axis, face.side);
                if(neighbor)
                {
                    const auto& other = neighbor->face(face.axis, 1 - face.side);
                    require(face.entity == other.entity, "Adjacent boxes do not share their face object");
                    require(face.orientation == -other.orientation, "Shared face orientation does not cancel");
                    require(tree.neighbor(neighbor, face.axis, 1 - face.side) == box, "Adjacency is not reciprocal");
                }
                else
                {
                    ++external;
                    const int axis = face.axis;
                    require(face.entity->key()[axis] == 0 || face.entity->key()[axis] == Tree::latticeSize,
                            "An interior face has no neighbor");
                }
            }
            if(level < tree.depth())
            {
                double sum = 0.0;
                require(box->subdivisions().size() == (1u << dim), "Wrong refinement child count");
                for(const auto& child : box->subdivisions())
                {
                    require(child->refinementParent() == box, "Missing refinement parent");
                    sum += child->measure();
                }
                require(std::abs(sum - box->measure()) < 1.e-12, "Refinement does not conserve volume");
            }
        }
        require(std::abs(volume - (hi - lo).prod()) < 1.e-11, "Level does not partition root volume");
        require(external == static_cast<size_t>(2 * dim * std::pow(n, dim - 1)), "Wrong exterior face count");
        require(tree.edges(level).size() == static_cast<size_t>(dim * n * std::pow(n + 1, dim - 1)),
                "Unique edge count is incorrect");
    }
    for(const auto& entry : tree.entities())
    {
        const auto& entity = entry.second;
        require(entity->boundary().size() == static_cast<size_t>(2 * entity->order()), "Invalid cell boundary");
        std::map<Key, int> secondBoundary;
        for(const auto& face : entity->boundary())
        {
            require(face.entity->order() + 1 == entity->order(), "Boundary has wrong order");
            require(face.entity->parents().at(entity->key()).lock() == entity, "Missing incidence back-link");
            for(const auto& edge : face.entity->boundary())
            {
                secondBoundary[edge.entity->key()] += face.orientation * edge.orientation;
            }
        }
        for(const auto& term : secondBoundary)
        {
            require(term.second == 0, "Boundary of a boundary is nonzero");
        }
        // Refined boundaries must equal the oriented union of fine-cell boundaries.
        if(!entity->subdivisions().empty())
        {
            std::map<Key, int> fineBoundary;
            for(const auto& child : entity->subdivisions())
            {
                for(const auto& face : child->boundary())
                {
                    fineBoundary[face.entity->key()] += face.orientation;
                }
            }
            for(const auto& face : entity->boundary())
            {
                if(face.entity->order() == 0)
                {
                    fineBoundary[face.entity->key()] -= face.orientation;
                }
                else
                {
                    require(!face.entity->subdivisions().empty(), "Coarse face has no subdivision");
                    for(const auto& subface : face.entity->subdivisions())
                    {
                        fineBoundary[subface->key()] -= face.orientation;
                    }
                }
            }
            for(const auto& term : fineBoundary)
            {
                require(term.second == 0, "Internal refinement faces do not cancel");
            }
        }
    }
    require(tree.locate(lo, 2) && tree.locate(hi, 2), "Root corners are not included");
    auto mid = tree.locate((lo + hi) * 0.5, 2);
    for(int axis = 0; axis < dim; ++axis)
    {
        require(mid->key()[axis] == Tree::latticeSize / 2, "Half-open split ownership is inconsistent");
    }
    Vector outside = lo;
    outside[0] -= 1.0;
    require(!tree.locate(outside, 2), "Outside point assigned to a box");
    outside[0] = std::numeric_limits<double>::quiet_NaN();
    require(!tree.locate(outside, 2), "NaN point assigned to a box");
    const auto root = tree.root();
    const auto revision = tree.revision();
    tree.setDepth(2);
    require(tree.root() == root && tree.revision() == revision, "No-op rebuild invalidates identities");
    bool rejected = false;
    try
    {
        tree.setDepth(Tree::maxDepth + 1);
    }
    catch(const std::invalid_argument&)
    {
        rejected = true;
    }
    require(rejected && tree.root() == root && tree.revision() == revision, "Failed rebuild changed the tree");
    tree.setDepth(1);
    require(tree.depth() == 1 && tree.revision() == revision + 1, "Coarsening did not replace hierarchy");
    rejected = false;
    try
    {
        tree.neighbor(root, 0, 0);
    }
    catch(const std::invalid_argument&)
    {
        rejected = true;
    }
    require(rejected, "Stale entity accepted by new topology");
    std::weak_ptr<typename Tree::Entity> weak;
    {
        Tree temporary(lo, hi, 2);
        weak = temporary.root()->face(0, 0).entity;
    }
    require(weak.expired(), "Shared topology contains an ownership cycle");
    rejected = false;
    try
    {
        Tree invalid(hi, lo);
    }
    catch(const std::invalid_argument&)
    {
        rejected = true;
    }
    require(rejected, "Reversed bounds accepted");
}

int main()
{
    try
    {
        checkTopology<1>();
        checkTopology<2>();
        checkTopology<3>();
        std::cout << "Box topology, orientation, refinement, lookup and lifetime checks passed in 1D/2D/3D\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
