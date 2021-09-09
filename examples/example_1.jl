using LoopFieldCalc

# Define current loop with radius of 50 cm, current of 1 A, centered at the origin
loop = LoopFieldCalc.CurrentLoop(
  radius = 0.5,         # radius in meters
  current = 1.0,        # current in Amperes
  center = LoopFieldCalc.CartesianPoint(0.0, 0.0, 0.0)    # cartesian coordinates of the loop center
)

# Define grid points along x, y, and z axes. We have a 2D grid of points in the y-z plane.
xs = 0.0:0.0
ys = -1.0:0.01:1.0
zs = -2.0:0.01:2.0

# Define current loops. Here, we only have the single loop, defined above
loops = [loop]

# Compute the magnetic field components
Bx, By, Bz = LoopFieldCalc.field_on_grid(loops, xs, ys, zs)

LoopFieldCalc.write_field("examples/example_1.dat", ys, zs, By[1, :, :], Bz[1, :, :])