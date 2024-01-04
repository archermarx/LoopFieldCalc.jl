using LoopFieldCalc
using Test

using LoopFieldCalc: CurrentLoop, CartesianPoint, CylindricalPoint, field_on_grid, field_at_point, μ0

@testset "Many loops" begin
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
		CurrentLoop(
			CylindricalPoint(inner_bobbin_offset + z, 0.0, 0.0);
			radius = inner_bobbin_radius,
			current = inner_bobbin_current / nturns,
		)
		for z in LinRange(0.0, inner_bobbin_width, nturns)
	]

	outer_bobbin = [
		CurrentLoop(
			CylindricalPoint(outer_bobbin_offset + z, 0.0, 0.0);
			radius = outer_bobbin_radius,
			current = outer_bobbin_current / nturns,
		)
		for z in LinRange(0.0, outer_bobbin_width, nturns)
	]

	loops = [outer_bobbin; inner_bobbin]
	Bx, By, Bz = field_on_grid(loops, xs, ys, zs)

	# check that this equals the sum of the fields of the individual loops
	Bx_indiv = zeros(size(Bx))
	By_indiv = zeros(size(By))
	Bz_indiv = zeros(size(Bz))

	for loop in loops
		Bx_i, By_i, Bz_i, _, _ = field_on_grid([loop], xs, ys, zs)
		@. Bx_indiv += Bx_i
		@. By_indiv += By_i
		@. Bz_indiv += Bz_i
	end

	@test all(Bx_indiv .≈ Bx)
	@test all(By_indiv .≈ By)
	@test all(Bz_indiv .≈ Bz)

	LoopFieldCalc.write_field("many_loops.dat", zs, ys, Bz[1, :, :]', By[1, :, :]')
	rm(joinpath("test", "many_loops.dat"), force = true)
	rm("many_loops.dat", force = true)
end

@testset "Centerline field" begin
	# Check that on-axis field gives the correct result
	R = 1.0
	I = 1.0
	center = CylindricalPoint(z = 0.0, r = 0.0, θ = 0.0)

	loop = CurrentLoop(center; radius = R, current = I)

	for z in 0.0:0.1:10.0
		Bx, By, Bz, _, _ = LoopFieldCalc.field_at_point(loop, CylindricalPoint(z, center.r, center.θ))
		@test Bx ≈ 0.0
		@test By ≈ 0.0
		@test Bz ≈ μ0 * I * R^2 / (2 * (z^2 + R^2)^1.5)
	end
end

include("runtests_elliptic.jl")
