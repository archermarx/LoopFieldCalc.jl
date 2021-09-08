using LoopFieldCalc
using Test

function test_many_loops()
	xs = -0.4:0.005:0.4
	ys = -0.4:0.005:0.4
	zs = -0.3:0.005:0.3

	angles = 0:π/3:5π/3
	r1 = 0.2

	field_strength = 1.0
	loop_radius = 0.04
	loop_current = 2 * loop_radius * field_strength / LoopFieldCalc.μ₀

	#outer_loops =
	#	[LoopFieldCalc.CurrentLoop(loop_radius, loop_current, LoopFieldCalc.CylindricalPoint(-0.01, r2, θ)) for θ in angles]

	inner_loops =
		[LoopFieldCalc.CurrentLoop(loop_radius, -loop_current, LoopFieldCalc.CylindricalPoint(-0.01, r1, θ)) for θ in angles]

	#loops = [outer_loops; inner_loops]

	Bx, By, Bz = LoopFieldCalc.field_on_box(inner_loops, xs, ys, zs)
	LoopFieldCalc.write_field("many_loops.dat", xs, ys, zs, Bx, By, Bz)
end

test_many_loops()
