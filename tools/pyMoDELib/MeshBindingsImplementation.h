/* This file is part of MoDELib, distributed under GNU GPL v2. */
#pragma once
#include "MeshBindings.h"
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <MeshModule.h>
#include <PolyhedralMesh.h>
namespace mesh_python
{
    namespace py = pybind11;
    using namespace model;

    template <class Key> py::object keyObject(const Key &key)
    {
        if constexpr (std::is_integral_v<Key>)
        {
            return py::int_(key);
        }
        else
        {
            py::tuple result(key.size());
            for (size_t k = 0; k < key.size(); ++k)
            {
                result[k] = py::int_(key[k]);
            }
            return result;
        }
    }
    // Keep the complete owning mesh alive, rather than making Python parent/child
    // references retain each other in cycles. Dynamic attributes participate in GC.
    template <class T> py::object retained(const T &value, const py::object &owner)
    {
        auto result = py::cast(value, py::return_value_policy::reference);
        if (!result.is_none() && !result.is(owner))
        {
            result.attr("_mesh_owner") = py::hasattr(owner, "_mesh_owner") ? owner.attr("_mesh_owner") : owner;
        }
        return result;
    }
    template <class Map> py::dict retainedMap(const Map &map, const py::object &owner)
    {
        py::dict result;
        for (const auto &entry : map)
        {
            result[keyObject(entry.first)] = retained(entry.second, owner);
        }
        return result;
    }
    template <int dim, int order> void bindSimplex(py::module_ &module)
    {
        using S = Simplex<dim, order>;
        const auto name = "Simplex" + std::to_string(dim) + std::to_string(order);
        auto b = py::class_<S, std::shared_ptr<S>>(module, name.c_str(), py::dynamic_attr())
                     .def_property_readonly_static("dim",
                                                   [](py::object)
                                                   {
                                                       return dim;
                                                   })
                     .def_property_readonly_static("order",
                                                   [](py::object)
                                                   {
                                                       return order;
                                                   })
                     .def_property_readonly("xID",
                                            [](const S &s)
                                            {
                                                return keyObject(s.xID);
                                            })
                     .def_property_readonly("sID",
                                            [](const S &s)
                                            {
                                                return s.sID;
                                            })
                     .def("center",
                          [](const S &s)
                          {
                              return s.center().eval();
                          })
                     .def("isBoundarySimplex", &S::isBoundarySimplex)
                     .def("isRegionBoundarySimplex", &S::isRegionBoundarySimplex)
                     .def("regionIDs", &S::regionIDs)
                     .def("vertices",
                          [](py::object self)
                          {
                              py::list result;
                              for (const auto *v : self.cast<const S &>().vertices())
                              {
                                  result.append(retained(v, self));
                              }
                              return result;
                          });
        if constexpr (order > 0)
        {
            b.def_property_readonly("vol0",
                                    [](const S &s)
                                    {
                                        return s.vol0;
                                    })
                .def("vertexPositionMatrix", &S::vertexPositionMatrix)
                .def("bary2pos", &S::bary2pos)
                .def("children",
                     [](py::object self)
                     {
                         py::list result;
                         for (const auto &child : self.cast<const S &>().children())
                         {
                             result.append(retained(child, self));
                         }
                         return result;
                     });
        }
        else
        {
            b.def_readonly("P0", &S::P0);
        }
        if constexpr (order < dim)
        {
            b.def("parents",
                  [](py::object self)
                  {
                      return retainedMap(self.cast<const S &>().parents(), self);
                  })
                .def("outNormal",
                     [](const S &s)
                     {
                         return s.outNormal();
                     })
                .def("outNormal",
                     [](const S &s, size_t id)
                     {
                         return s.outNormal(id);
                     });
        }
        else
        {
            b.def_property_readonly("region",
                                    [](py::object self)
                                    {
                                        return retained(self.cast<const S &>().region, self);
                                    })
                .def("pos2bary", &S::pos2bary);
        }
    }
    template <int dim, int order> void bindPolyhedron(py::module_ &module)
    {
        using C = Polyhedron<dim, order>;
        const auto name = "Polyhedron" + std::to_string(dim) + std::to_string(order);
        auto b = py::class_<C, std::shared_ptr<C>>(module, name.c_str(), py::dynamic_attr())
                     .def_property_readonly_static("dim",
                                                   [](py::object)
                                                   {
                                                       return dim;
                                                   })
                     .def_property_readonly_static("order",
                                                   [](py::object)
                                                   {
                                                       return order;
                                                   })
                     .def_property_readonly("key",
                                            [](const C &c)
                                            {
                                                return keyObject(c.key);
                                            })
                     .def("simplexIDs",
                          [](const C &c)
                          {
                              py::list result;
                              for (const auto *s : c.simplices())
                              {
                                  result.append(keyObject(s->xID));
                              }
                              return result;
                          })
                     .def("simplices",
                          [](py::object self)
                          {
                              py::list result;
                              for (const auto *s : self.cast<const C &>().simplices())
                              {
                                  result.append(retained(s, self));
                              }
                              return result;
                          });
        if constexpr (order > 0)
        {
            b.def_property_readonly("children",
                                    [](py::object self)
                                    {
                                        return retainedMap(self.cast<const C &>().children, self);
                                    });
        }
        if constexpr (order < dim)
        {
            b.def_property_readonly("parents",
                                    [](py::object self)
                                    {
                                        py::dict result;
                                        for (const auto &entry : self.cast<const C &>().parents)
                                        {
                                            if (const auto p = entry.second.lock())
                                            {
                                                result[keyObject(entry.first)] = retained(p, self);
                                            }
                                        }
                                        return result;
                                    });
        }
        if constexpr (order == dim)
        {
            b.def_readonly("regionID", &C::regionID)
                .def("faces",
                     [](py::object self)
                     {
                         return retainedMap(self.cast<const C &>().faces(), self);
                     });
        }
        if constexpr (order == dim - 1)
        {
            b.def_readonly("sID", &C::sID)
                .def_readonly("regionIDs", &C::regionIDs)
                .def("isExternal", &C::isExternal)
                .def("center", &C::center)
                .def("outNormal",
                     [](const C &c)
                     {
                         return c.outNormal();
                     })
                .def("outNormal",
                     [](const C &c, int id)
                     {
                         return c.outNormal(id);
                     })
                .def_property_readonly("periodicShift",
                                       [](const C &c)
                                       {
                                           return c.periodicFacePair.first;
                                       })
                .def_property_readonly("periodicPartner",
                                       [](py::object self)
                                       {
                                           return retained(self.cast<const C &>().periodicFacePair.second, self);
                                       });
        }
        if constexpr (order == 0)
        {
            b.def_readonly("position", &C::position);
        }
    }
    template <int dim, int order> void bindOrders(py::module_ &module)
    {
        bindSimplex<dim, order>(module);
        bindPolyhedron<dim, order>(module);
        if constexpr (order > 0)
        {
            bindOrders<dim, order - 1>(module);
        }
    }
    template <int dim, int order = dim>
    py::dict simplicesAt(const SimplicialMesh<dim> &mesh, int requested, const py::object &owner)
    {
        if (requested == order)
        {
            return retainedMap(mesh.template simplexFactory<order>().liveSimplices(), owner);
        }
        if constexpr (order > 0)
        {
            return simplicesAt<dim, order - 1>(mesh, requested, owner);
        }
        throw py::value_error("Simplex order must be between zero and dim");
    }
    template <int dim, int order = dim>
    py::dict cellsAt(const PolyhedralMesh<dim> &mesh, int requested, const py::object &owner)
    {
        if (requested == order)
        {
            return retainedMap(mesh.template factory<order>(), owner);
        }
        if constexpr (order > 0)
        {
            return cellsAt<dim, order - 1>(mesh, requested, owner);
        }
        throw py::value_error("Polyhedron order must be between zero and dim");
    }
    template <int dim> void bindMesh(py::module_ &module)
    {
        using M = SimplicialMesh<dim>;
        using P = PolyhedralMesh<dim>;
        using R = SimplexReader<dim>;
        using V = Eigen::Matrix<double, dim, 1>;
        using A = Eigen::Matrix<double, dim, dim>;
        bindOrders<dim, dim>(module);
        const auto suffix = std::to_string(dim);
        py::class_<R>(module, ("SimplexReader" + suffix).c_str())
            .def(py::init<>())
            .def(
                "read",
                [](R &reader, const std::string &file, const A &transform, const V &origin)
                {
                    reader.read(file, transform, origin);
                },
                py::arg("filename"), py::arg("A") = A::Identity().eval(), py::arg("x0") = V::Zero().eval())
            .def("nodes",
                 [](const R &reader)
                 {
                     return reader.nodes();
                 })
            .def("elements",
                 [](const R &reader)
                 {
                     py::dict result;
                     for (const auto &entry : reader.elements())
                     {
                         result[py::int_(entry.first)] =
                             py::make_tuple(keyObject(entry.second.first), entry.second.second);
                     }
                     return result;
                 });
        using Boundary = MeshRegionBoundary<dim>;
        py::class_<Boundary>(module, ("MeshRegionBoundary" + suffix).c_str(), py::dynamic_attr())
            .def_readonly("regionBndID", &Boundary::regionBndID)
            .def("faces",
                 [](py::object self)
                 {
                     return retainedMap(self.cast<const Boundary &>().faces(), self);
                 })
            .def("simplices",
                 [](py::object self)
                 {
                     py::list result;
                     for (const auto *simplex : self.cast<const Boundary &>().simplices())
                     {
                         result.append(retained(simplex, self));
                     }
                     return result;
                 });
        py::class_<P, std::shared_ptr<P>>(module, ("PolyhedralMesh" + suffix).c_str(), py::dynamic_attr())
            .def(
                "cells",
                [](py::object self, int order)
                {
                    return cellsAt(self.cast<const P &>(), order, self);
                },
                py::arg("order") = dim)
            .def("region",
                 [](py::object self, size_t id)
                 {
                     return retained(&self.cast<const P &>().region(id), self);
                 })
            .def("face",
                 [](py::object self, size_t id)
                 {
                     return retained(&self.cast<const P &>().face(id), self);
                 })
            .def("validate", &P::validate);
        py::class_<M, std::shared_ptr<M>>(module, ("SimplicialMesh" + suffix).c_str(), py::dynamic_attr())
            .def(py::init<>())
            .def(py::init<const std::string &, const A &, const V &, const std::set<int> &>(), py::arg("filename"),
                 py::arg("A") = A::Identity().eval(), py::arg("x0") = V::Zero().eval(),
                 py::arg("periodicFaceIDs") = std::set<int>{})
            .def_property_readonly_static("dim",
                                          [](py::object)
                                          {
                                              return dim;
                                          })
            .def("xMin",
                 [](const M &mesh)
                 {
                     return mesh.xMin().eval();
                 })
            .def("xMax",
                 [](const M &mesh)
                 {
                     return mesh.xMax().eval();
                 })
            .def("xCenter", &M::xCenter)
            .def("volume", &M::volume)
            .def("periodicBasis",
                 [](const M &mesh)
                 {
                     return Eigen::MatrixXd(mesh.periodicBasis());
                 })
            .def("identifyParallelFaces", &M::identifyParallelFaces, py::arg("periodicFaceIDs"))
            .def("region",
                 [](py::object self, size_t id)
                 {
                     return retained(self.cast<const M &>().regions().at(id), self);
                 })
            .def("regions",
                 [](py::object self)
                 {
                     return retainedMap(self.cast<const M &>().regions(), self);
                 })
            .def(
                "simplices",
                [](py::object self, int order)
                {
                    return simplicesAt(self.cast<const M &>(), order, self);
                },
                py::arg("order") = dim)
            .def("regionBoundaries",
                 [](py::object self)
                 {
                     py::dict result;
                     for (const auto &entry : self.cast<const M &>().regionBoundaries())
                     {
                         result[py::make_tuple(entry.first.first, entry.first.second)] = retained(&entry.second, self);
                     }
                     return result;
                 })
            .def("polyhedralMesh",
                 [](py::object self)
                 {
                     return retained(&self.cast<const M &>().polyhedralMesh(), self);
                 })
            .def("search",
                 [](py::object self, const V &point)
                 {
                     if (self.cast<const M &>().simplices().empty())
                     {
                         return py::object(py::none());
                     }
                     const auto found = self.cast<const M &>().search(point);
                     return found.first ? retained(found.second, self) : py::none();
                 })
            .def("periodicSimplexPairs",
                 [](py::object self)
                 {
                     py::dict result;
                     for (const auto &entry : self.cast<const M &>().periodicSimplexPairs())
                     {
                         result[keyObject(entry.first->xID)] = py::make_tuple(
                             retained(entry.second.partner, self), entry.second.shift, entry.second.vertexPermutation);
                     }
                     return result;
                 });
    }
} // namespace mesh_python
