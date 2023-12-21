using LoopFieldCalc

# Define grid on which we compute magnetic field strength.
# This is a 2D grid in the y-z plane
xs = 0.0:0.0
ys = -0.4:0.005:0.4
zs = -0.3:0.005:0.3

# Define number of wire turns per coil
nturns = 10

# Define geometry of inner coil
inner_coil_radius = 0.1
inner_coil_width = 0.15
inner_coil_offset = -0.1
inner_coil_current = 1.0

# Generate loops for inner coil
inner_coil = [
    LoopFieldCalc.CurrentLoop(
        LoopFieldCalc.CylindricalPoint(inner_coil_offset + z, 0.0, 0.0);
        radius = inner_coil_radius,
        current = inner_coil_current / nturns,
    )
    for z in LinRange(0.0, inner_coil_width, nturns)
]

# Define geometry for outer coil
outer_coil_radius = 0.25
outer_coil_width = 0.15
outer_coil_offset = -0.15
outer_coil_current = 0.5 * inner_coil_current

# Generate loops for outer coil
outer_coil = [
    LoopFieldCalc.CurrentLoop(
        LoopFieldCalc.CylindricalPoint(outer_coil_offset + z, 0.0, 0.0);
        radius = outer_coil_radius,
        current = outer_coil_current / nturns,
    )
    for z in LinRange(0.0, outer_coil_width, nturns)
]

# Combine inner and outer coil loops into a single vector and compute field
loops = [outer_coil; inner_coil]
Bx, By, Bz = LoopFieldCalc.field_on_grid(loops, xs, ys, zs)

# Write output. By reversing z and y and transposing the magnetic field matrices, we
# can rotate the output so that the loop axis is aligned with the x axis
LoopFieldCalc.write_field("examples/example_2.dat", zs, ys, Bz[1, :, :]', By[1, :, :]')
