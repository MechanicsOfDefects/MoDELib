# Box topology and DDqt Boxes tab

This is the geometry foundation for loop-FMM. It contains no moments, plane
intersections, source assignment, periodic summation, force or energy kernels.

## Geometry API

`BoxTree<dim>` is implemented for dimensions 1, 2 and 3, with declarations in
`include/FMM/BoxTree.h` and explicit instantiations in `src/FMM/BoxTree.cpp`.
It builds a uniform dyadic hierarchy within finite, strictly increasing lower
and upper bounds. Boxes may be rectangular. Depth zero contains just the root;
each additional level divides every active coordinate in two. The initial
implementation limits depth to five to bound uniform-tree memory use.

Each `BoxEntity<dim>` has a runtime `order()` from zero (vertex) to `dim` (box).
The tree's keyed shared-pointer factory gives adjacent entities the same face,
edge and vertex objects. Keys are exact dyadic integer bounds, independent of
floating-point comparisons. There are two distinct relationships:

- `boundary()` stores lower-order incidences with axis, side, and cubical-chain
  orientation. `parents()` contains weak reverse links. Opposite box incidences
  on a shared face have opposite orientations, and the boundary of a boundary
  is zero.
- `subdivisions()` contains same-order refinement children;
  `refinementParent()` is a weak link. Coarse-face subdivisions are the same
  objects as the fine-box faces they cover. Vertices do not subdivide.

The factory owns all entities. Downward links are shared, upward links weak;
there is no observer registration or raw-pointer ownership. A retained old
entity remains readable after a rebuild, but is no longer part of the live tree.

Use `boxes(level)`, `edges(level)`, `root()`, and `entities()` for inspection.
`neighbor(box, axis, side)` follows incidence and returns the same-level neighbor
or null at the root boundary. `locate(point, level)` uses half-open ownership at
interior splits and includes the root's upper corner in the last cell. Points
outside the root or with nonfinite coordinates return null.

`setDepth()` builds a replacement before publishing it, and increments
`revision()` only after success. Repeating the current depth is a no-op. Cache
users must invalidate entity-dependent data when the revision changes. This is
uniform refinement, not an adaptive tree: mixed-level leaf neighbors and
hanging-face traversal are deliberately not exposed yet.

## Ownership and viewer

`DislocationDynamicsBase<3>::boxTree` owns a shared tree initialized from the
simulation mesh's axis-aligned bounds. Its initial depth is read from `inputFiles/DD.txt` using `fmmTreeDepth`
(required integer 0..5 whenever DD is enabled). Missing settings or a missing
DD file cause an input error. Explicitly specify zero for one root box. The configured value is available as
`DislocationDynamicsBase::fmmTreeDepth`, including in Python.
The Boxes tab displays this tree, rather than creating a second
viewer-only hierarchy. If DD is disabled, the tab explains why it is unavailable.

- Enable **Show boxes** to display the wireframe (off initially).
- **Tree depth** is read-only. Change `fmmTreeDepth` in `inputFiles/DD.txt` and
  reload the simulation to change it.
- **Display level** selects an existing level without modifying topology.
- **Opacity** and **Line width** control the appearance.
- Counts refer to the displayed level; shared edges are rendered only once.

Visibility toggles retain settings, and frame changes retain box identities and
reuse the VTK geometry. A revision change refreshes the displayed geometry.
The standard periodic-image viewer can replicate this actor for inspection;
that is only a rendering operation, not periodic identification in the tree.

The root encloses the full mesh, including empty space around non-convex grains.
Boxes are not clipped to grains and their axis alignment does not select the
physical periodic directions. Intersections with GlidePlanes will be a later
layer sharing the entities introduced here.

## Tests

`boxTree` checks dimensions 1/2/3, shared topology, incidence orientation,
boundary-of-boundary cancellation, coarse/fine boundary consistency, partition
measure, lookup on splits and corners, failed rebuilds, stale handles, and
ownership lifetime. `boxTreeActor` checks VTK wireframe counts and bounds,
controls, frame reuse, revisions, disabled-DD behavior and renderer cleanup.
Both are registered with CTest under the existing numerical/integration labels.

## Plane cuts

Enable **Glide-plane cuts** in the Boxes tab and select an existing glide plane.
The display level selects its box-cut polygons; orange segments are emitted once
per shared edge. Filled polygons show the infinite plane inside the root box,
not the slipped area or the grain boundary. A configuration without glide planes
has no cuts to display. Plane selection is retained across frames while that plane
exists. Removing a plane clears its geometry.

`include/FMM/PlaneBoxIntersection.h` and `src/FMM/PlaneBoxIntersection.cpp`
implement immutable cuts for dimension 3. `GlidePlane::boxIntersection(tree, level)`
returns a shared cached cut. Cache identity includes tree ownership, tree revision,
and level. Calls on the same glide plane are serialized by a mutex; rebuilding a
tree must not race its readers. Old returned cuts remain valid after rebuilds.

Each unique box edge is intersected once per cached level. A cut node is identified
by its supporting edge, or by its vertex for endpoint intersections. Each face cut
is computed once; segments share their endpoint identities. Cell polygons follow
the resulting degree-two adjacency graph, without angular sorting. Their boundary
incidences are oriented about the supplied plane normal. Neighboring polygons
reference identical segments with opposite incidences.

Coincident internal faces use half-open cell ownership: the box on the positive
coordinate side owns the area, independently of normal sign. Outer upper faces
are retained. Point/edge tangencies have no area. Distances within 64 machine
epsilons times the coordinate scale are treated as coincident; this avoids
independent rounding decisions at shared vertices.

`planeBoxIntersection` tests analytic sections, random planes, reversed normals,
vertex/edge/face coincidences, ownership, shared identities and orientations,
and area conservation at every refinement level. `physicsBases` tests cache
reuse/invalidation and failed requests. `boxTreeActor` additionally tests actual
glide-plane cut rendering, frame reuse, level selection and disappearing planes.

Loop clipping, moments, force kernels and energy kernels remain unimplemented.
