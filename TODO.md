# Development TODO

## Ewald FMM accuracy and performance (deferred)

Use `multipoleMode=FMM` with `coreType=Cai` for current work. Keep
`ewaldFMM` available for controlled experiments; its optimization is deferred.

- [ ] Separate long-field interpolation error into self-box, neighbor-box,
  and well-separated contributions before choosing an optimization.
- [ ] Test one extra level of internal long-field refinement at fixed Ewald
  width and interpolation order 4. Preserve the original short-field 27-box
  stencil. Integrate fine source coefficients directly from segments and
  evaluate fine target locals without projecting back onto the coarse basis.
- [ ] Measure translation singular values; evaluate tolerance-controlled SVD
  compression of cached operators. Do not assume self/neighbor blocks have
  low rank. Compression must not be confused with interpolation refinement.
- [ ] Benchmark at matched total accuracy, counting initialization, geometry
  updates, translation-cache memory, field evaluation, and GMRES applications.
- [ ] Choose settings that control both short-field truncation and long-field
  interpolation errors before reconsidering Ewald FMM as the default.

### Current evidence (2026-11-08)

Shared full/near/far interfaces and cached Ewald self/neighbor translations
are implemented and tested. Full kernels use their direct expressions.
The remaining accuracy issue is independent of that refactor.

Frozen annealing benchmark: 53 loops, 1060 segments, 5300 targets; Helmholtz
core width 2; depth 3, interpolation order 4, Ewald width factor 0.75.
Relative L2 errors, normalized by the direct full field:

| Field | Short truncation | Long interpolation | Combined |
|---|---:|---:|---:|
| DD stress | 1.150% | 0.964% | 1.510% |
| DD concentration | 0.893% | 0.145% | 0.900% |

Increasing order does not reduce the truncation floor. Decreasing Ewald
width suppresses truncation but increases interpolation error. Field-error
norms do not add; concentration source strengths were held fixed.

Diagnostic scripts, frozen-input hashes and full results are in the sibling
workspace `../MoDELib2/benchmarks/ewald_error_decomposition_20261108/`.
Timing comparisons are in `../MoDELib2/benchmarks/common_near_far_20261108/`.

## Kokkos grouped M2L products (deferred)

Keep Eigen grouped matrix products on CPU. Device backends retain the
ungrouped Kokkos translation path until a faster grouped implementation is validated.

- [ ] Evaluate KokkosBlas/batched or vendor GEMM for grouped M2L on GPU.
  Preserve cached grouping and translations, cache-budget fallbacks, and
  target-owned accumulation in deterministic pair order.
- [ ] Benchmark on actual GPU backends and check stress, concentration, and
  climb velocities against the Eigen reference before enabling the new path.

The handwritten tiled Kokkos prototype was removed after the 2026-11-09
CPU benchmark: M2L was 3.44x slower than Eigen and the complete timestep was
6.6% slower. Output differences were at roundoff level; GPU performance was
not tested. Benchmark inputs and results remain in the sibling workspace
`../MoDELib2/benchmarks/portable_m2l_cpu_20261109/`.

## True periodic FMM: exterior operator and physical conventions (deferred)

Step 1 is implemented as an internal explicit-image neighborhood: interaction
entries reference a central source box and a physical lattice shift. Shifted
P2P/segment evaluation and M2L reuse the original source geometry and moments.
`localImageShifts` is an internal assembly option for this component, not a
user-selected approximation to a periodic solution. It is not populated from
`periodicImageSize`. Complete periodic FMM requests remain rejected until the
following work is finished. `periodicImageSize` remains a direct (`none`) control.

- [ ] **Step 2: exterior periodic root operator.** Construct and cache the map
  from root moments to root locals for all images outside the explicit
  neighborhood. Subtract that exact neighborhood from the periodic kernel to
  avoid double counting. Combine the implemented physical-distance root
  neighborhood with the exterior correction before enabling periodic solves;
  allow a larger explicit neighborhood if required by the requested accuracy.
- [ ] **Step 3: accurate operator construction.** Implement a convergent
  periodic Green-function/lattice-summation evaluator (initialization only).
  Define the summation convention and requested accuracy. Cache by lattice,
  active directions, explicit neighborhood, interpolation order, kernel/core/
  material parameters, tolerance and boundary convention. Verify independence
  of the explicit neighborhood size once the exterior correction is included.
- [ ] **Step 4: stress compatibility.** Check net Burgers content and admissible
  periodic dislocation sources; separate fluctuating stress from the
  macroscopic stress/strain controlled by the load controller. Define zero-mode
  and nonperiodic-direction conditions explicitly.
- [ ] **Step 4: concentration compatibility.** For full periodicity, account
  explicitly for zero net source in the steady diffusion equation, the mean
  concentration, and the zero-mean Green function solving delta minus 1/V.
  Integrate source balance into CD/climb rather than silently subtracting the
  source mean. Specify the far-boundary conditions for partial periodicity.
- [ ] Validate full periodicity across opposite faces, conservation and
  compatibility, invariance under lattice translations, and convergence with
  operator accuracy. The existing explicit-neighborhood comparisons only
  validate step 1 and are not evidence of infinite-sum convergence.

## Anisotropic elasticity with grain-continuous corrections (deferred)

For grain-specific anisotropic DD reference fields, use the existing
`NodeSharingPolicy::grainContinuous` for the FEM displacement correction.
Separate correction DOFs at GBs must remain coupled: the total displacement
and traction must satisfy the perfectly bonded interface conditions.

- [ ] Implement grain-oriented anisotropic elastic tensors and consistent
  grain-specific DD displacement/stress reference fields.
- [ ] Assemble the correction problem with separate grain-side GB DOFs and
  impose `[u_correction] = -[u_DD]`, so `[u_total] = 0`.
- [ ] Enforce total traction balance using outward normals:
  `(sigma_DD^+ + sigma_correction^+) n^+ +
   (sigma_DD^- + sigma_correction^-) n^- = 0`.
  Include the known DD traction mismatch in the interface loading.
- [ ] Select and validate the interface coupling formulation and its rigid-mode
  constraints. Do not turn this into independent per-grain elasticity solves;
  stress transmission across GBs must be retained.
- [ ] Test displacement-jump cancellation, traction balance, patch tests,
  mesh convergence, and recovery of the common-isotropic continuous-correction
  result when the reference field is shared across grains.

This discontinuity concerns the correction field, not a physical opening or
sliding of the GB. Anisotropy alone does not require discontinuous total
displacement; a common-reference formulation remains an alternative.
