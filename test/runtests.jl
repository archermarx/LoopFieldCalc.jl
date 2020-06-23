using LoopFieldCalc
using Test, Plots, GeometricalPredicates

@testset "Field calculation tests" begin
	rT = 0.0275
	zmax = 80.
	rmax = 40.
  	h = 1

	zgrid = [z for r in 0:h:rmax, z in 0:h:zmax] * rT
	rgrid = [r for r in 0:h:rmax, z in 0:h:zmax] * rT

	B0 = 0.01

	myloop = CurrentLoop(3.5 * rT; Bfield = B0)
	Bz, Br, λ, χ = calculatefieldonmesh(zgrid, rgrid, myloop)

	throatline = LoopFieldCalc.extractfieldline(0, rT, zgrid, rgrid, λ, myloop)
	polyline = [throatline; -1 rmax+0.1; -1 1]
	#LoopFieldCalc.replace_in_polygon(zgrid, rgrid, (Bz, Br, λ, χ), polyline, NaN)

	#pt1 = (20, 10)
	#pt2 = (10, 20)

	#@show Bz[20, 10]
	#@show Bz[10, 20]

	#polygon = Polygon(polyline)
	#@show inpolygon(polygon, Point2D(pt1...))
	#@show inpolygon(polygon, Point2D(pt2...))

	p1 = plot(yaxis = ("r", (0, rmax)), xaxis = ("z", (0, zmax)), aspect_ratio = 1)
	heatmap!(0:h:zmax, 0:h:rmax, Bz; cb = :none)
	plot!(polyline[:, 1], polyline[:, 2]; label = "", lw = 2, color = :black)

	#scatter!([pt1, pt2], label = "")
	display(p1)

	LoopFieldCalc.save_field(zgrid, rgrid, Bz, Br, abspath("test.dat"))

	boundarypts = [
		throatline;
		[zmax rmax; zmax 0; 0 0] .* rT;
	]
	casename = "LoopNozzle"
	generate_boundary_file(boundarypts, casename)
	p1 = plot(boundarypts[:, 1], boundarypts[:, 2], aspect_ratio = 1, label = "")

	display(p1)
end

function generate_boundary_file(nodes, casename)
	filename = "Boundary_" * casename * ".dat"
	npts = size(nodes, 1)
    open(filename, "w") do io
        write(io, string(npts))
        for i = 1:npts
            write(io, "\n", string(nodes[i, 1]), " ", string(nodes[i, 2]))
        end
    end
end
