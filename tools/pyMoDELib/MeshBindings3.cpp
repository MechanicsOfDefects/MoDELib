#include "MeshBindingsImplementation.h"
void bindMesh3(pybind11::module_ &module)
{
    mesh_python::bindMesh<3>(module);
}
