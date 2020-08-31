using LoopFieldCalc
using Test, Plots

@testset "Field calculation tests" begin
	rL = 0.075
	rT = 0.5 * rL
	zmin, zmax = (-3.0, 7) .* rL
	rmin, rmax = (0.0, 6) .* rL

  	npts_z = 50
	npts_r = 40
	IB = 20.0

	#loop_z = [-0.3, 0.3] * rL
	loop_z = [-0.3, 0.3] * rL
	#loop_z = [0.0] * rL
	nloops = length(loop_z)

	zT = 0.63 * rL
	zmin = zT

	zrange = range(zmin, length = npts_z, stop = zmax)
	rrange = range(0, length = npts_r, stop = rmax)
	zgrid = [z for r in rrange, z in zrange]
	rgrid = [r for r in rrange, z in zrange]

	Bz, Br, λ, χ = zero(zgrid), zero(zgrid), zero(zgrid), zero(zgrid)

	for i in 1:nloops
		myloop = CurrentLoop(rL, IB/nloops, loop_z[i])
		Bzi, Bri, λi, χi = calculatefieldonmesh(zgrid, rgrid, myloop)
		Bz += Bzi
		Br += Bri
		λ += λi
		χ += χi
	end

	myloop = CurrentLoop(rL, IB)
	throatline = LoopFieldCalc.extractfieldline(zT, rT, zgrid, rgrid, λ)
	polyline = throatline ./rL

	#LoopFieldCalc.replace_in_polygon(zgrid, rgrid, (Bz, Br, λ, χ), polyline, NaN)

	#pt1 = (20, 10)
	#pt2 = (10, 20)

	#@show Bz[20, 10]
	#@show Bz[10, 20]

	#polygon = Polygon(polyline)
	#@show inpolygon(polygon, Point2D(pt1...))
	#@show inpolygon(polygon, Point2D(pt2...))
	B = hypot.(Bz, Br)
	p1 = plot(
		yaxis = ("r", (rmin, rmax)./rL),
		xaxis = ("z", (zmin, zmax)./rL),
		aspect_ratio = 1,
		minorticks = 4,
		minorgrid = true,
		tick_direction = :out
	)
	heatmap!(zrange./ rL, rrange ./ rL, Br; cb = :none)
	plot!(polyline[:, 1], polyline[:, 2]; label = "", lw = 1, ls = :dash, color = :red)
	#plot!(line2[:, 1]./rL, line2[:, 2]./rL, label = "", lw = 1, ls = :dash, color = :red)
	#plot!(line3[:, 1]./rL, line3[:, 2]./rL, label = "", lw = 1, ls = :dash, color = :red)
	hline!([rT / rL], linestyle = :dash, color = :yellow, label = "")
	#scatter!([pt1, pt2], label = "")
	display(p1)
	casename = "Little2019"

	Bscale = 21 * IB / 10000 / maximum(B)
	Bz = Bscale * Bz
	Br = Bscale * Br
	LoopFieldCalc.save_field(zgrid .- zT, rgrid, Bz, Br, abspath("BField_" * casename * ".dat"))

	boundarypts = [
		0.0 0.0;
		0.0 rT;
		0.0 rmax;
		zmax.-zT rmax;
		zmax.-zT 0.0
	]
	#rmax_T = maximum(throatline[:, 2]) / rT
	#@show rmax, rmax_T
	#if rmax_T < rmax
	#	rmax = rmax_T
	#end

	#@show rmax, rmax_T

	#if throatline[2, end] > throatline[1, end]
	#	boundarypts = [
	#		throatline;
	#		[zmax rmax; zmax 0; 0 0] .* rT;
	#	]
	#else
	#	boundarypts = [
	#		[0 0; zmax 0; zmax rmax] .* rT;
	#		throatline;
	#	]
	#end


	LoopFieldCalc.generate_boundary_file(boundarypts, casename)
	#p1 = plot(boundarypts[:, 1], boundarypts[:, 2], aspect_ratio = 1, label = "")

	#display(p1)
end
