#include <PolyhedralMesh.h>
#include <MeshModule.h>
namespace model
{
    template <int dim>
    Eigen::Matrix<double, dim, 1> PolyhedralFaceIncidence<dim>::outNormal(size_t regionID) const
        requires(dim == 3)
    {
        const auto owner = face.lock();
        if (!owner || !owner->parents.count(regionID))
        {
            throw std::runtime_error("Face incidence does not belong to this region");
        }
        return owner->outNormal(regionID);
    }
    template <int dim>
    PolyhedralFaceIncidence<dim> PolyhedralFaceIncidence<dim>::otherFace(size_t regionID) const
        requires(dim == 3)
    {
        const auto owner = face.lock();
        if (!owner || !owner->parents.count(regionID))
        {
            throw std::runtime_error("Face incidence does not belong to this region");
        }
        std::shared_ptr<Polyhedron<dim, 2>> other;
        for (const auto &entry : edge->parents)
        {
            const auto candidate = entry.second.lock();
            if (candidate && candidate != owner && candidate->parents.count(regionID))
            {
                if (other)
                {
                    throw std::runtime_error("Non-manifold grain edge");
                }
                other = candidate;
            }
        }
        if (!other)
        {
            throw std::runtime_error("Open grain boundary edge");
        }
        return other->perimeter.at(edge->faceIncidences.at(other->key));
    }
    template <> void PolyhedralMesh<3>::orientMeshFaces()
    {
        for (const auto &entry : factory<2>())
        {
            const auto &face = entry.second;
            Eigen::Vector3d area = Eigen::Vector3d::Zero();
            const auto origin = face->perimeter.front().source()->position;
            for (const auto &incidence : face->perimeter)
            {
                area += (incidence.source()->position - origin).cross(incidence.sink()->position - origin);
            }
            if (area.dot(face->asPlane().unitNormal) < 0)
            {
                std::reverse(face->perimeter.begin(), face->perimeter.end());
                for (auto &incidence : face->perimeter)
                {
                    incidence.reversed = !incidence.reversed;
                }
            }
            for (size_t k = 0; k < face->perimeter.size(); ++k)
            {
                face->perimeter[k].index = k;
                face->perimeter[k].edge->faceIncidences.at(face->key) = k;
            }
        }
    }
    template <> void PolyhedralMesh<3>::validate() const
    {
        validateLevel<3>();
        for (const auto &entry : factory<2>())
        {
            const auto &face = *entry.second;
            if (face.simplices().empty())
            {
                continue; // Topology-only input has no underlying mesh geometry.
            }
            for (const auto &parent : face.parents)
            {
                const auto region = parent.second.lock();
                for (const auto &incidence : face.perimeter)
                {
                    const auto other = incidence.otherFace(region->key);
                    const auto otherOwner = other.face.lock();
                    if (incidence.sink() != incidence.next().source() || incidence.next().prev().edge != incidence.edge)
                    {
                        throw std::runtime_error("Disconnected polyhedral perimeter");
                    }
                    const bool flipped = face.regionIDs.first != int(region->key);
                    const bool otherFlipped = otherOwner->regionIDs.first != int(region->key);
                    if ((incidence.reversed != flipped) == (other.reversed != otherFlipped))
                    {
                        throw std::runtime_error("Inconsistent oriented grain incidences");
                    }
                }
            }
        }
    }
    template struct PolyhedralFaceIncidence<3>;
    template class PolyhedralMesh<3>;
} // namespace model
