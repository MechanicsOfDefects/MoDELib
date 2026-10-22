"""Dimension-specific mesh bindings, connectivity, periodicity, and ownership."""
import gc
import itertools
from pathlib import Path
import sys
import tempfile
import weakref
import numpy as np

sys.path.insert(0, sys.argv[1])
import pyMoDELib as model


def box_mesh(path, dim, split=False):
    vertices = list(itertools.product((0, 1), repeat=dim))
    indices = {point: k + 1 for k, point in enumerate(vertices)}
    cells = []
    for axes in itertools.permutations(range(dim)):
        point = [0] * dim
        cell = [indices[tuple(point)]]
        for axis in axes:
            point[axis] = 1
            cell.append(indices[tuple(point)])
        cells.append(cell)
    if split:
        # Two translated copies share a facet, with distinct region IDs.
        original_count = len(cells)
        for cell in list(cells):
            translated = []
            for vertex_id in cell:
                point = list(vertices[vertex_id - 1])
                point[0] += 1
                point = tuple(point)
                if point not in indices:
                    indices[point] = len(vertices) + 1
                    vertices.append(point)
                translated.append(indices[point])
            cells.append(translated)
    lines = ["$MeshFormat", "2.2 0 8", "$EndMeshFormat", "$Nodes", str(len(vertices))]
    for point, index in indices.items():
        xyz = list(point) + [0] * (3 - dim)
        lines.append(f"{index} " + " ".join(map(str, xyz)))
    lines += ["$EndNodes", "$Elements", str(len(cells))]
    for index, cell in enumerate(cells, 1):
        region = 2 if split and index > original_count else 1
        lines.append(f"{index} { {1: 1, 2: 2, 3: 4}[dim]} 2 1 {region} " + " ".join(map(str, cell)))
    path.write_text("\n".join(lines + ["$EndElements", ""]))


def check(dim, path):
    Mesh = getattr(model, f"SimplicialMesh{dim}")
    mesh = Mesh(str(path), periodicFaceIDs={-1})
    assert Mesh.dim == dim
    assert Mesh().search(np.zeros(dim)) is None
    np.testing.assert_allclose(mesh.volume(), 1.)
    assert mesh.periodicBasis().shape == (dim, dim)
    topology = mesh.polyhedralMesh()
    assert isinstance(topology, getattr(model, f"PolyhedralMesh{dim}"))
    topology.validate()
    for order in range(dim, -1, -1):
        simplices = mesh.simplices(order)
        assert simplices
        for key, simplex in simplices.items():
            assert isinstance(simplex, getattr(model, f"Simplex{dim}{order}"))
            assert simplex.xID == key
            assert simplex.order == order and simplex.dim == dim
            assert len(simplex.vertices()) == order + 1
            if order:
                assert len(simplex.children()) == order + 1
                for child in simplex.children():
                    assert child.parents()[key] is simplex
                bary = np.full(order + 1, 1 / (order + 1))
                np.testing.assert_allclose(simplex.bary2pos(bary), simplex.center())
            else:
                np.testing.assert_allclose(simplex.P0, simplex.center())
        for cell in topology.cells(order).values():
            assert isinstance(cell, getattr(model, f"Polyhedron{dim}{order}"))
            if order:
                for child in cell.children.values():
                    assert child.parents[cell.key] is cell
    for source, (partner, shift, permutation) in mesh.periodicSimplexPairs().items():
        original = mesh.simplices(dim - 1)[source]
        np.testing.assert_allclose(original.center() + shift, partner.center(), atol=1e-12)
        for i, vertex in enumerate(original.vertices()):
            np.testing.assert_allclose(vertex.P0 + shift, partner.vertices()[permutation[i]].P0, atol=1e-12)
    assert mesh.periodicSimplexPairs()
    root_simplex = mesh.search(np.full(dim, 0.3))
    assert root_simplex is not None
    assert mesh.search(np.full(dim, 2.)) is None
    np.testing.assert_allclose(root_simplex.bary2pos(root_simplex.pos2bary(np.full(dim, 0.3))), np.full(dim, 0.3))
    assert root_simplex.region is mesh.region(1)
    physical = Mesh(str(path))
    assert physical.periodicBasis().shape == (dim, 0)
    assert not physical.periodicSimplexPairs()
    faces = physical.polyhedralMesh().cells(dim - 1)
    selected = next(iter(faces))
    partial = physical
    partial.identifyParallelFaces({selected})
    assert partial.periodicBasis().shape == (dim, 1)
    reader = getattr(model, f"SimplexReader{dim}")()
    reader.read(str(path))
    assert len(reader.elements()) == len(mesh.simplices())
    scaled = Mesh(str(path), 2 * np.eye(dim), np.full(dim, 0.25))
    np.testing.assert_allclose(scaled.xMin(), np.full(dim, -0.5))
    np.testing.assert_allclose(scaled.xMax(), np.full(dim, 1.5))
    np.testing.assert_allclose(scaled.volume(), 2.**dim)
    for method in (mesh.simplices, topology.cells):
        try:
            method(dim + 1)
        except ValueError:
            pass
        else:
            raise AssertionError("Invalid order accepted")


def lifetime(dim, path):
    mesh = getattr(model, f"SimplicialMesh{dim}")(str(path))
    ref = weakref.ref(mesh)
    cell = mesh.region(1)
    del mesh
    gc.collect()
    assert ref() is not None
    assert cell.simplexIDs()
    retained = cell.simplices()[0].vertices()[0]
    del cell
    gc.collect()
    assert ref() is not None
    assert retained.parents()
    del retained
    gc.collect()
    assert ref() is None, "Mesh ownership cycle"


with tempfile.TemporaryDirectory(prefix="modelib-python-mesh-") as directory:
    for dim in (3, 2, 1):
        path = Path(directory) / f"box{dim}.msh"
        box_mesh(path, dim)
        check(dim, path)
        lifetime(dim, path)
        box_mesh(path, dim, split=True)
        mesh = getattr(model, f"SimplicialMesh{dim}")(str(path))
        assert len(mesh.regions()) == 2
        boundary = mesh.regionBoundaries()[(1, 2)]
        assert isinstance(boundary, getattr(model, f"MeshRegionBoundary{dim}"))
        assert len(boundary.faces()) == 1
        for face in boundary.faces().values():
            assert not face.isExternal()
            assert set(face.parents) == {1, 2}
            assert mesh.region(1).faces()[face.key] is mesh.region(2).faces()[face.key]
        assert boundary.simplices()
assert model.SimplicialMesh is model.SimplicialMesh3
assert model.Polyhedron2 is model.Polyhedron32
print("PASS: Python mesh dimensions 1–3, all orders, periodicity, connectivity, ownership")
