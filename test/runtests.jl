using LoopFieldCalc
using Test


function test_many_loops()
	xs = 0.0:0.0
	ys = -0.4:0.005:0.4
	zs = -0.3:0.005:0.3

	nturns = 10

	inner_bobbin_radius = 0.1
	inner_bobbin_width = 0.15
	inner_bobbin_offset = -0.1
	inner_bobbin_current = 1.0

	outer_bobbin_radius = 0.25
	outer_bobbin_width = 0.15
	outer_bobbin_offset = -0.15
	outer_bobbin_current = 0.3 * inner_bobbin_current

	inner_bobbin = [
		LoopFieldCalc.CurrentLoop(
			radius = inner_bobbin_radius,
			current = inner_bobbin_current / nturns,
			center = LoopFieldCalc.CylindricalPoint(inner_bobbin_offset + z, 0.0, 0.0)
		)
		for z in LinRange(0.0, inner_bobbin_width, nturns)
	]

	outer_bobbin = [
		LoopFieldCalc.CurrentLoop(
			radius = outer_bobbin_radius,
			current = outer_bobbin_current / nturns,
			center = LoopFieldCalc.CylindricalPoint(outer_bobbin_offset + z, 0.0, 0.0)
		)
		for z in LinRange(0.0, outer_bobbin_width, nturns)
	]

	loops = [outer_bobbin; inner_bobbin]
	Bx, By, Bz = LoopFieldCalc.field_on_box(loops, xs, ys, zs)
	LoopFieldCalc.write_field("many_loops.dat", zs, ys, Bz[1, :, :]', By[1, :, :]')
end
test_many_loops()

