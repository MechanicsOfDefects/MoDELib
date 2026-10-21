#include "MeshBindingsImplementation.h"
void bindMesh1(pybind11::module_& module)
{
    mesh_python::bindMesh<1>(module);
}
