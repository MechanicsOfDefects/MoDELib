// Compile the real SimplicialMesh adapters in dimensions without library instantiations.
#include <MeshModule.h>
#include <PolyhedralMesh.h>

template <int dim> void compilePolyhedralAdapter(const model::SimplicialMesh<dim> &mesh)
{
    const auto& topology = mesh.polyhedralMesh();
    topology.validate();
}
template void compilePolyhedralAdapter<1>(const model::SimplicialMesh<1> &);
template void compilePolyhedralAdapter<2>(const model::SimplicialMesh<2> &);
template void compilePolyhedralAdapter<4>(const model::SimplicialMesh<4> &);
