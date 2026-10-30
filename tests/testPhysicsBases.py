"""Python ownership and shared-context API regression. Takes module directory and repository root."""
import gc
import math
from pathlib import Path
import shutil
import sys
import tempfile
import numpy as np

sys.path.insert(0, sys.argv[1])
import pyMoDELib as model
root = Path(sys.argv[2])

with tempfile.TemporaryDirectory(prefix="modelib-python-bases-") as directory:
    folder = Path(directory)
    inputs = folder / "inputFiles"
    inputs.mkdir()
    (folder / "F").mkdir()
    shutil.copy(root / "Library/Materials/W.txt", inputs / "W.txt")
    shutil.copy(root / "Library/Meshes/unitCube24.msh", inputs / "mesh.msh")

    def configure(physics, length):
        (inputs / "DefectiveCrystal.txt").write_text(
            f"physics={physics};\nfmmTreeDepth=2;\nuseFEM=0;\nNsteps=0;\nmaxResolveSteps=0;\ndtMax=1;\n"
            "outputFrequency=1;\noutputBinary=0;\nstartAtTimeStep=0;\n"
            "periodicImageSize=0 0 0;\nEwaldLengthFactor=0;\n")
        (inputs / "polycrystal.txt").write_text(
            "materialFile=W.txt;\nmeshFile=mesh.msh;\nabsoluteTemperature=300;\n"
            "C2G1=1 0 0 0 1 0 0 0 1;\n"
            f"F={length:.17g} 0 0 0 {length:.17g} 0 0 0 {length:.17g};\n"
            "X0=0 0 0;\nperiodicFaceIDs=-1;\n")

    configure("ElasticDeformation", 123.456)
    traits = model.DefectiveCrystalTraitsIO(str(folder))
    assert traits.configFolder == str(folder / "config")
    assert traits.auxFolder == traits.configFolder
    default_config = model.DefectiveCrystalConfigIO()
    custom_config = model.DefectiveCrystalConfigIO(str(folder / "custom"))
    base = model.DefectiveCrystalBase(str(folder))
    assert base.fmmContext.fmmTreeDepth == 2
    assert base.dislocationDynamicsBase is None
    assert base.clusterDynamicsBase is None
    assert base.elasticDeformationBase is base.elasticDeformation()
    assert base.hasPhysics("ElasticDeformation")
    assert base.periodicLatticeBasis.shape == (3, 3)
    assert base.periodicCoordinates(np.zeros(3)).shape == (3,)
    generator = model.MicrostructureGenerator(base)
    assert isinstance(generator.configIO, model.DefectiveCrystalConfigIO)
    del base
    gc.collect()
    try:
        generator.addPrismaticLoopDensity(model.PrismaticLoopDensitySpecification())
    except RuntimeError as error:
        assert "not enabled" in str(error)
    else:
        raise AssertionError("DD generation was allowed without DD")
    del generator

    configure("DislocationDynamics", 200 / math.sqrt(3))
    base = model.DefectiveCrystalBase(str(folder))
    dd = base.dislocationDynamicsBase
    assert dd is base.dislocationDynamics()
    context = base.fmmContext
    assert context.fmmTreeDepth == 2 and context.treeDepth == 2
    assert context.boxCount(2) == 64
    assert not hasattr(dd, "fmmTreeDepth")
    plane = dd.glidePlane(np.zeros(3), np.array([0., 0., 1.]), 1)
    periodic = dd.periodicGlidePlane(np.zeros(3), np.array([0., 0., 1.]), 1)
    assert dd.glidePlaneCount() >= 1
    assert dd.periodicGlidePlaneCount() == 1
    del base
    gc.collect()
    assert context.boxCount(2) == 64
    assert np.isfinite(dd.EwaldLength)
    assert plane.planeBaseID() == periodic.referencePlane.planeBaseID()
    del dd
    gc.collect()
    assert plane.planeBaseID() == periodic.referencePlane.planeBaseID()
    del plane, periodic

    # CD context methods and ownership, without constructing DD factories.
    shutil.copy(root / "tests/spatialCDtest/inputFiles/Zr4_Fitted.txt", inputs / "Zr4_Fitted.txt")
    configure("ClusterDynamics", 123.456)
    poly = inputs / "polycrystal.txt"
    poly.write_text(poly.read_text().replace("W.txt", "Zr4_Fitted.txt").replace("=300", "=553"))
    base = model.DefectiveCrystalBase(str(folder))
    cd = base.clusterDynamicsBase
    assert cd is base.clusterDynamics()
    assert base.fmmContext.fmmTreeDepth == 2
    assert base.dislocationDynamicsBase is None
    del base
    gc.collect()
    assert np.isfinite(cd.equilibriumMobileConcentration(0.)).all()
    assert np.isfinite(cd.getDlocal()).all()
print("PASS: Python base API, activation, and lifetimes")

# Region/face APIs expose the same factory-owned polyhedra at every level.
mesh = model.SimplicialMesh(str(root / "Library/Meshes/bicrystal_48.msh"),
                           np.eye(3), np.zeros(3), set())
regions = mesh.regions()
assert len(regions) == 2
shared_faces = {}
for region_id, region in regions.items():
    assert mesh.region(region_id) is region
    assert region.regionID == region_id
    assert region.simplexIDs()
    for face_id, face in region.faces().items():
        assert region.children[face_id] is face
        assert face.key == face_id == face.sID
        assert face.parents[region_id] is region
        assert face.simplexIDs()
        if face_id in shared_faces:
            assert shared_faces[face_id] is face
            assert not face.isExternal()
        shared_faces[face_id] = face
        for edge_key, edge in face.children.items():
            assert tuple(edge.key) == edge_key
            assert edge.parents[face_id] is face
            assert len(edge.simplexIDs()) == 1
            for vertex_key, vertex in edge.children.items():
                assert tuple(vertex.key) == vertex_key
                assert vertex.parents[edge_key] is edge
                assert len(vertex.simplexIDs()) == 1
assert any(not face.isExternal() for face in shared_faces.values())
print("PASS: Python unified polyhedral hierarchy")
