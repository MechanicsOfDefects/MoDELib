# Coupled solved quantities

MicrostructureContainer owns SolvedQuantityRegistry, including its ordered SolvedQuantityGroups. DefectiveCrystal calls finalizeSolvedQuantities after adding all physics. Registration is virtual on MicrostructureBase; custom physics should override it. The default adapter publishes stress and calls the existing solve once.

A contribution is one update operation and may publish multiple quantity names. A dependency selects every other operation publishing a quantity, or one named physics provider. Required missing providers are rejected at finalization. Conditional dependencies are evaluated at the start of each solve (DD climb vs glide); strongly connected components are rebuilt then, before any update. Physics registration order does not set execution order.

ED stress is solved before DD velocities and CD concentrations. DD stress evaluation uses current geometry directly. During climb, DD velocity/concentration and CD concentration form a cyclic group. ED displacement rates and CD immobile rates follow convergence. During glide, no concentration feedback edge exists. Geometry, time, and immobile populations do not advance inside these groups.

The group's maximum componentwise normalized unrelaxed update must be at most one. Each denominator is absoluteTolerance + relativeTolerance * max(abs(old),abs(new)). CD measures its FEM correction relative to each species' zero-stress equilibrium concentration, preventing the equilibrium background from masking convergence. Absolute floors are species-specific. DD monitors per-species scalar rates and vector nodal velocities. Relaxation applies after measuring the unrelaxed update; it cannot manufacture convergence. Cyclic state callbacks restore iterated fields/rates on failure; no time step is accepted. Prepared matrices and boundary caches may remain allocated and are refreshed on the next solve.

Optional DefectiveCrystal.txt settings (defaults shown):

```
solvedQuantityMaxIterations=100;
solvedQuantityRelativeTolerance=1.e-3;
solvedQuantityRelaxation=1.0;
solvedQuantityVelocityAbsoluteTolerance=1.e-24;
solvedQuantityConcentrationAbsoluteFraction=1.e-10;
```

Velocity tolerance uses MoDELib velocity units. Concentration fraction multiplies each species' equilibrium concentration. maxResolveSteps is retained only as a deprecated C++ compatibility field (zero) and is no longer read or used. Old inputs remain readable.

The DD climb interaction matrix is constructed on the first group iteration and reused subsequently. Quadrature stress and dislocation equilibrium concentration remain fixed; only CD concentration is refreshed. CD sink equilibrium concentrations are cached for the duration of the group, then invalidated at its next first iteration. Existing weak-problem factorization reuse remains intact. The DD-to-CD concentration influence matrices H are assembled once per coupled solve, including retries, using the same segment kernels, periodic shifts, grain filtering, and shared network-node endpoints as direct evaluation. Each later iteration computes H times the current scalar climb velocities. Only sink/periodic FEM nodes are stored unless reactions require all nodes. Other concentration providers remain evaluated directly; standalone CD solves also use direct evaluation. The dense matrices retain no ownership of network nodes and are invalidated at the next first iteration.

Extending this system to concentration-dependent stress requires registering that stress dependency explicitly, adding ED to the coupled component, and invalidating stress caches accordingly. The current adapters describe the current stress-independent CD implementation.

Regression coverage includes a non-periodic DD/CD/ED calculation with actual dislocation stress, plus partial/full periodic calculations with a controlled zero reference stress and line-tension-driven climb. The latter isolate coupling from the finite-image DD stress approximation: its non-periodic residual can otherwise prescribe incompatible sink concentrations at equivalent periodic nodes. This is an existing field-evaluation limitation, not relaxed by the scheduler. Tests also force non-convergence, check restoration of concentration and climb state, then retry successfully.
