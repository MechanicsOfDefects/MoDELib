#include "MeshBindingsImplementation.h"
void bindMesh2(pybind11::module_& module)
{
    mesh_python::bindMesh<2>(module);
}
