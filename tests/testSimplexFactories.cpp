#include <SimplicialMesh.h>
#include <iostream>
#include <stdexcept>

namespace
{
    void require(bool condition, const char *message)
    {
        if (!condition)
        {
            throw std::runtime_error(message);
        }
    }

    template <short order> void check(model::SimplicialMesh<3> &mesh)
    {
        for (const auto &entry : mesh.simplexFactory<order>().liveSimplices())
        {
            const auto &simplex = entry.second;
            require(mesh.simplexFactory<order>().get(simplex->xID) == simplex, "Factory identity mismatch");
            if constexpr (order > 0)
            {
                for (const auto &child : simplex->children())
                {
                    require(mesh.simplexFactory<order - 1>().get(child->xID) == child, "Child is not shared");
                    require(child->parents().at(simplex->xID) == simplex, "Missing reciprocal parent");
                }
            }
        }
        if constexpr (order > 0)
        {
            check<order - 1>(mesh);
        }
    }

    void exercise(const std::string &path)
    {
        std::weak_ptr<model::Simplex<3, 0>> weakVertex;
        std::shared_ptr<model::Simplex<3, 2>> retainedFace;
        {
            model::SimplicialMesh<3> mesh(path, Eigen::Matrix3d::Identity(), Eigen::Vector3d::Zero(), {-1});
            check<3>(mesh);
            auto root = mesh.simplices().begin()->second;
            const auto key = root->xID;
            auto reversed = key;
            std::reverse(reversed.begin(), reversed.end());
            require(mesh.simplexFactory<3>().getOrCreate(reversed, root->region->regionID) == root,
                    "Canonical factory lookup duplicated a simplex");
            const auto before = mesh.simplexFactory<3>().liveSimplices().size();
            bool rejected = false;
            try
            {
                mesh.simplexFactory<3>().getOrCreate(key, root->region->regionID + 100000);
            }
            catch (const std::runtime_error &)
            {
                rejected = true;
            }
            require(rejected, "Conflicting region was accepted");
            require(mesh.simplexFactory<3>().liveSimplices().size() == before, "Failed creation changed the hierarchy");
            check<3>(mesh);
            // Keep children alive while replacing their parent under the same canonical key.
            const auto regionID = root->region->regionID;
            {
                const auto children = root->children();
                std::weak_ptr<model::Simplex<3, 3>> previous = root;
                mesh.simplices().erase(key);
                root.reset();
                require(previous.expired(), "Children retain their deleted parent");
                root = mesh.simplexFactory<3>().getOrCreate(key, regionID);
                mesh.simplices().emplace(key, root);
                for (const auto &child : children)
                {
                    require(child->parents().at(key) == root, "Recreation did not replace expired parent link");
                }
            }
            check<3>(mesh);
            retainedFace = root->children().front();
            weakVertex = retainedFace->children().front()->children().front();
            std::weak_ptr<model::Simplex<3, 3>> weakRoot = root;
            root.reset();
            mesh.simplices().clear();
            require(weakRoot.expired(), "Parent/child ownership cycle retains volume simplex");
            require(retainedFace->parents().empty(), "Expired parents remain in topology queries");
            require(mesh.simplexFactory<3>().liveSimplices().empty(), "Factory owns volume simplices");
            mesh.simplexFactory<3>().pruneExpired();
            require(mesh.simplexFactory<3>().empty(), "Expired factory entries were not pruned");
            require(!weakVertex.expired(), "Retained face lost its children");
        }
        require(retainedFace->parents().empty(), "Retained face has dangling parent links");
        for (const auto* vertex : retainedFace->vertices())
        {
            require(vertex != nullptr && vertex->P0.allFinite(), "Retained simplex lost its vertex geometry");
        }
        retainedFace.reset();
        require(weakVertex.expired(), "Simplex hierarchy leaked after final owner release");
    }
} // namespace

int main(int argc, char **argv)
{
    try
    {
        require(argc == 2, "Expected repository root");
        for (const auto *name : {"bicrystal_48", "polycrystalCube_10_periodicX", "polycrystalCube_10_periodicXY"})
        {
            exercise(std::string(argv[1]) + "/Library/Meshes/" + name + ".msh");
        }
        std::cout << "PASS: shared simplex identity, weak parent links, failed creation, and lifetime cleanup\n";
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
