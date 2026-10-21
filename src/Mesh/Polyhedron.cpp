#include <MeshModule.h>
#include <PolyhedralMesh.h>
#include <ConvexHull.h>
namespace model
{
    template <int dim, int order>
    void
    Polyhedron<dim, order>::buildSingleFace(const Simplex<dim, dim - 1> *newStart,
                                            std::shared_ptr<Polyhedron<dim, dim - 1>> newFace,
                                            std::map<std::vector<size_t>, const Simplex<dim, dim - 1> *> &allSimplices)
        requires(order == dim)
    {
        allSimplices.erase(std::vector<size_t>(newStart->xID.begin(), newStart->xID.end()));
        const auto boundaryNeighbors(newStart->boundaryNeighbors());
        for (const auto &neighbor : boundaryNeighbors)
        {
            const auto neighborRgnIDs(neighbor->regionIDs());
            if (neighborRgnIDs.find(this->regionID) != neighborRgnIDs.end() // neighbor in region
                && allSimplices.find(std::vector<size_t>(neighbor->xID.begin(), neighbor->xID.end())) !=
                       allSimplices.end()) // neighbor found in allSimplices
            {
                if ((newStart->outNormal() - neighbor->outNormal()).norm() < FLT_EPSILON)
                { // same plane
                    allSimplices.erase(std::vector<size_t>(neighbor->xID.begin(), neighbor->xID.end()));
                    newFace->insert(neighbor);
                    buildSingleFace(neighbor, newFace, allSimplices);
                }
            }
        }
        faces().emplace(newFace->sID, newFace);
    }

    /**********************************************************************/
    template <int dim, int order>
    void Polyhedron<dim, order>::buildFaces()
        requires(order == dim)
    { /*!Constructs and stores the external PlanarMeshFace(s) of this MeshRegion
       * This function supports non-convex regions.
       */
        const auto t0 = std::chrono::system_clock::now();
        std::cout << "MeshRegion " << this->regionID << " buiding PlanarMeshFace's " << std::flush;

        std::map<std::vector<size_t>, const Simplex<dim, dim - 1> *> allSimplices;
        for (const auto &simplex : simplices())
        { // collect each possible boundary/region_boundary simplices
            for (const auto &child : simplex->children())
            {
                if (child->isBoundarySimplex())
                {
                    allSimplices.emplace(std::vector<size_t>(child->xID.begin(), child->xID.end()), child.get());
                }
            }
        }

        while (allSimplices.size())
        {
            auto newFace = cellFactory.owner->createFace(allSimplices.begin()->second);
            buildSingleFace(allSimplices.begin()->second, newFace, allSimplices);
        }

        for (auto &face : faces())
        {
            face.second->finalize();
        }
        std::cout << magentaColor << " ["
                  << (std::chrono::duration<double>(std::chrono::system_clock::now() - t0)).count() << " sec]"
                  << defaultColor << std::endl;
    }

    template <int dim, int order>
    void Polyhedron<dim, order>::update()
        requires(order == dim)
    {
        buildFaces();
    }
    template <int dim, int order>
    typename Polyhedron<dim, order>::VectorDim Polyhedron<dim, order>::outNormal(const std::set<size_t> &ids) const
        requires(order == dim)
    {
        VectorDim result = VectorDim::Zero();
        for (const auto id : ids)
        {
            result += faces().at(id)->outNormal();
        }
        const double norm = result.norm();
        return norm > FLT_EPSILON ? (result / norm).eval() : VectorDim::Zero();
    }

    template <int dim, int order>
    void Polyhedron<dim, order>::initialize(const SimplexType *simplex)
        requires(order == dim - 1)
    {
        this->regionIDs = {*simplex->regionIDs().begin(), *simplex->regionIDs().rbegin()};
        this->n = simplex->outNormal(this->regionIDs.first);
        this->c = simplex->center();
        this->insert(simplex);
    }
    template <int dim, int order>
    const std::vector<const Simplex<dim, 0> *> &Polyhedron<dim, order>::convexHull() const
        requires(order == dim - 1)
    {
        return this->hull;
    }

    /**********************************************************************/
    template <int dim, int order>
    void Polyhedron<dim, order>::finalize()
        requires(order == dim - 1)
    { /*!\todo This function uses ConvexHull and therefore will yield a wrong
       * face boundary for non-convex faces. In order to support non-convex faces
       * we need to have a structure with face edges.
       */

        //            buildEdges();

        // Recompute this->c and this->n to minimize floating point errors
        this->n.setZero();
        this->c.setZero();
        for (const auto &simplex : *this)
        {
            this->n += simplex->outNormal(this->regionIDs.first) / this->simplices().size();
            this->c += simplex->center() / this->simplices().size();
        }
        this->n.normalize();

        for (const auto &simplex : *this)
        { // Make sure all simplices belong to face
            if (fabs((this->c - simplex->center()).dot(this->n)) > FLT_EPSILON)
            {
                throw std::runtime_error("simplex does not belong to PlanarMeshFace.");
            }
        }

        // Compute ConvexHull of vertices
        std::set<const Simplex<dim, 0> *> vertices;
        for (const auto &simplex : *this)
        { // Collect vertices on this face
            for (const auto &v : simplex->vertices())
            {
                vertices.insert(v);
            }
        }

        //            PlanarMeshFaceHull<dim>::makeHull(vertices,outNormal(),center,this->hull);

        if (vertices.size())
        {

            //            assert(dim==3 && "ALGORITHM ONLY VALID IN dim=3");

            const Plane<dim> plane(this->c, this->n);
            ConvexHull<2, Simplex<dim, 0>> hull;
            for (const auto &v : vertices)
            { // compute local position in the plane as insert in ConvexHull
                const auto x(plane.localPosition(v->P0));
                hull.emplace(std::array<double, 2>{x[0], x[1]}, v);
            }

            const auto hullPts = hull.getPoints();
            this->hull.clear();
            for (const auto &hp : hullPts)
            {
                this->hull.push_back(hp.t);
            }
            if (this->hull.size() < 3)
            {
                std::cout << "vertices.size()=" << vertices.size() << std::endl;
                throw std::runtime_error("PlanarMeshFace hull has size " + std::to_string(this->hull.size()));
            }

            // Recompute center using Hull
            this->c.setZero();
            for (const auto &v : convexHull())
            {
                this->c += v->P0 / convexHull().size();
            }
        }
    }

    /**********************************************************************/
    template <int dim, int order>
    Eigen::Matrix<double, dim, 1> Polyhedron<dim, order>::outNormal() const
        requires(order == dim - 1)
    {
        return isExternal() ? this->n : Eigen::Matrix<double, dim, 1>::Zero();
    }

    /**********************************************************************/
    template <int dim, int order>
    Eigen::Matrix<double, dim, 1> Polyhedron<dim, order>::outNormal(const int &k) const
        requires(order == dim - 1)
    {
        return k == this->regionIDs.first
                   ? this->n
                   : (k == this->regionIDs.second ? (-this->n).eval() : Eigen::Matrix<double, dim, 1>::Zero());
    }

    /**********************************************************************/
    template <int dim, int order>
    const Eigen::Matrix<double, dim, 1> &Polyhedron<dim, order>::center() const
        requires(order == dim - 1)
    {
        return this->c;
    }

    /**********************************************************************/
    template <int dim, int order>
    Plane<dim> Polyhedron<dim, order>::asPlane() const
        requires(order == dim - 1)
    {
        return Plane<dim>(center(), this->n);
    }

    template struct Polyhedron<3, 3>;
    template struct Polyhedron<3, 2>;
} // namespace model
