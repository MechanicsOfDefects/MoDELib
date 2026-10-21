# Thermal cross-slip

Model 2 (HEX basal–prismatic cross-slip) first collects continuous eligible
screw branches. Physical faces, grain boundaries, shared junction segments,
and junction nodes delimit branches; same-grain periodic seam connectors do
not. Eligibility does not use the remeshing minimum length.

For each target slip system, integrate the existing stress-dependent rate
over the physical branch. The probability of an event over the accepted
motion interval is `-expm1(-dt * sum(rates))`. Conditional on an event, select
the target in proportion to its integrated rate, retaining equal-rate targets.
At most one event is proposed per branch per step. Small timesteps are needed
when multiple events on one branch would be physically significant.

`csLb_SI` (basal source) and `csLp_SI` (prismatic source) are positive physical
nucleation lengths in metres. They retain their role in the rate prefactor
and now also set the seed length. A shorter branch cannot nucleate. The
position is selected using the chosen target's segment-integrated rates,
with uniform position within each segment and the center clamped so the
whole nucleus fits. Spatial rate resolution therefore follows the field
quadrature/segment resolution; uniform-rate results are invariant under
subdivision, and nonuniform results converge with field resolution.

The complete screw branch is aligned transactionally to the selected plane
before inserting nucleus endpoints and a midpoint. Alignment uses a
length-weighted branch center. Rejected alignment leaves connectivity and
positions unchanged. The seed occupies only the prescribed physical length;
existing interior nodes and periodic crossings are retained. Shared nodes
remain shared, and the reversed parent path cancels Burgers vectors on the
seed's return side. New free nodes allow a single long segment to nucleate.

Thermal sampling runs immediately before motion, after the current timestep
has been selected. Velocities are recomputed after successful seed insertion.
Athermal model 1 keeps its existing timing and selection rule.

Slip families are identified by their reciprocal-lattice plane normals,
independently of mobility model names or crystal rotation. Material
slip-family enablement comes from `enabledSlipSystems`; the remaining
thermal attempt-frequency, barrier, and activation-volume parameters retain
their existing names. Thermal parameters must be supplied for model 2; the
stock material files do not supply calibrated thermal cross-slip parameters.

Tests:
- `thermalCrossSlipSampling`: identical draws across refinements, competing
  target frequencies, equal rates, nonuniform piecewise rates, short branches,
  and zero rate/time.
- `thermalCrossSlipBranches`: actual HEX rate integration and branch discovery
  from one to 200 segments, including segments below the remeshing threshold.
- `thermalCrossSlipTopology`: fixed 300b nuclei with different segmentations,
  physical/partial/full periodic domains, preserved junction membership,
  Burgers cancellation, failed alignment, and seed expansion.
- `crossSlipTopology`: existing athermal topology regression.
