# SurfaceMeshActor

The SurfaceMesh tab, implemented by `SurfaceMeshActor`, refines the mesh's exterior boundary triangles and displays
`X + amplification * u(X)` using the same displacement field as Fields, summed
across the selected microstructures. Enable **Show deformed surface** after loading
a frame. The surface is opaque gray with lighting, recomputed normals, and no
scalar coloring. Hide other actors when you want an unobstructed surface view.

**Displacement physics** lists the same physics contributions as the Fields tab.
All are selected initially. Changing a checkbox recomputes the visible surface
using only the selected contributions, including both sides of slip walls.
Selecting none displays the undeformed surface. Choices persist across frame
updates; changes made while hidden take effect when the surface is enabled.
These selections are independent of the Fields tab.

**Mesh size** is the maximum reference triangle area in Burgers-vector units
squared, matching the area convention in Fields. Press Enter in the mesh-size field to apply it. A common subdivision count keeps shared boundary edges
conforming, including nonconvex boundaries and holes represented by the input
mesh. Internal grain faces are excluded. Refinement is limited to 500000 triangles.

**Displacement amplification** changes the displayed displacement without
recomputing the field. Zero displays the reference surface. Reference coordinates
are retained, so neither changing frames nor changing amplification accumulates
deformation. Each successful frame load, including playback, refreshes the field
before image capture. Hidden surfaces defer computation until enabled.

The exterior surface is cut along the live glide planes in GlidePlaneFactory,
using their mesh-intersection faces to restrict cuts to the appropriate grain
boundary portions. Separate vertices on the two sides retain the same reference
position but sample displacement just inside their respective surface regions.
Wall triangles join the displaced trace edges, giving sharp slip steps, including
where traces cross. Zero-height walls collapse in the reference configuration.
The wall follows the full displacement jump, including tangential components;
it is normal to the surface only when the jump is normal to that surface.

Cuts are reconstructed from the current factory on each frame update. Only
external faces are displayed; grain interfaces are excluded. Coplanar traces
are merged. The 500000-triangle limit also applies after cutting and adding walls.
One-sided samples move toward the cut polygon's centroid, with offsets scaled
to the domain to avoid the solid-angle evaluator's on-plane tolerance and capped
to stay inside each polygon. Two samples on the same side extrapolate to the
trace, removing offset errors for affine fields. The limit remains a numerical
approximation, particularly for very closely spaced traces. The
mesh size still controls resolution of the smooth field between traces.

This is geometric surface topography, not an electron-scattering or SEM detector
model. Cuts follow the factory's glide planes; discontinuities from other sources
(such as nonplanar loop facets) are not explicitly cut. Field capabilities and
boundary-condition limitations are those of the loaded simulation.

**Clip surface** trims the displayed, displaced surface using a plane. Enter its
center as global `x y z` coordinates in b and its normal as `nx ny nz`, then press
Enter. The normal is normalized automatically and must be nonzero. The retained
half-space satisfies `normal · (displayedPosition - center) >= 0`; reverse the
normal to retain the other side. The default plane passes through the box center
with normal `(1,0,0)`, and clipping starts disabled. Settings persist across frames
and apply to slip walls too. Clipping updates the display without rebuilding the
mesh or recomputing displacement. The cut is open, with no artificial cap.

The surface uses Phong interpolation with matte diffuse shading and ambient
fill to reveal changes in slope gently. An oblique camera-relative light is added
while the surface is visible and removed when hidden. It also illuminates other
actors sharing the renderer while active. Normal splitting at 30 degrees keeps
sharp edges and steep slip walls crisp while smoothing shading within connected,
gently curved regions. Lighting highlights actual relief; a truly planar region
still has uniform diffuse shading. Increase displacement amplification when the
physical relief is too small to distinguish at the current view scale.

The clip-position slider translates the plane along its applied unit normal,
without changing the entered center. Its midpoint is zero offset; its limits
are one reference-box diagonal in either direction. The signed offset is shown
in b. Dragging updates clipping immediately without recomputing displacement.
Frame changes and disabling/re-enabling clipping preserve the offset. Applying
a different center or normal with Enter resets the slider to zero. Unapplied
text edits do not affect dragging.

Displacement sampling uses OpenMP when available, with one result per surface
point and both one-sided samples evaluated by the same worker. Surfaces with
fewer than 256 sample points use serial sampling to avoid thread overhead.
Physics selections and FEM element pointers are collected before launching
workers. The configuration remains fixed during evaluation; Qt and VTK updates,
slip-trace splitting, and rendering stay on the main thread. Worker errors are
reported after the parallel loop completes, without displaying partial results.
The existing OpenMP thread settings (including `OMP_NUM_THREADS`) apply. Builds
without OpenMP retain serial execution. This does not add GPU field evaluation
or make the UI asynchronous during computation.

Surface shading uses a matte gray material with ambient fill and no specular
highlight, so small changes in triangle orientation do not dominate the slip
relief. Exactly coincident display vertices are merged before normal calculation
to smooth adjacent cut polygons. This changes only rendering connectivity, not
field sampling or displacement. No distance tolerance is used, so distinct
positions on opposite sides of even a small step remain separate. The 30-degree
normal-splitting threshold preserves sharp walls and box edges.

The Shading group adjusts Ambient, Diffuse, and Specular coefficients independently
from 0 to 1. Defaults are 0.25, 0.75, and 0 respectively. Changes render immediately
without rebuilding the mesh or evaluating displacement. Values are retained across
frame updates and surface visibility changes for the current viewer session.
