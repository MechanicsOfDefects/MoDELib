from generate_mesh import MeshSpec, generate

spec = MeshSpec(
    geometry="hexagonal",
    size_nm=500.0,
    height_nm=800.0,
    target_elements=0,
    lc_nm=60.0,
    element_order=1,
    boundary_layer=True,
    boundary_thickness_nm=100.0,
    boundary_lc_nm=20.0,
    boundary_layers_across=None,
    optimize=True,
)

result = generate(spec, force=True)

print(result.path)
print(result.n_nodes)
print(result.n_tets)
print(result.F_b)