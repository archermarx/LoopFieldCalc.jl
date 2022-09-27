module LoopFieldCalc

using Elliptic2, Printf, Contour

const μ₀ = 4π * 1e-7

include("geometry.jl")

Base.@kwdef struct CurrentLoop{R, I, X, Y, Z}
	radius::R
	current::I
	center::CartesianPoint{X, Y, Z}
end

CurrentLoop(;radius, current, x, y, z) =
	CurrentLoop(radius, current, CartesianPoint(x, y, z))

CurrentLoop(;radius, current, z, r, θ) =
	CurrentLoop(radius, current, cylindrical2cartesian(CylindricalPoint(z, r, θ)))

CurrentLoop(radius, current, cyl::CylindricalPoint) =
	CurrentLoop(radius, current, cylindrical2cartesian(cyl))

B0(loop) = μ₀ * loop.current / 2 / loop.radius


field_at_point(loop, point::CylindricalPoint) = field_at_point(loop, cylindrical2cartesian(point))

"""
	field_at_point(loop, point::CartesianPoint)
Compute the field strength due to a current loop at a point in the cartesian plane. Returns
Bx, By, and Bz (three components of magnetic field vector at this point) as well as the
magnetic stream function λ and magnetic scalar potential χ.
"""
function field_at_point(loop, point::CartesianPoint)

	# Compute (z, r, θ), the coordinates of the point relative to the loop center in cylindrical coordinates
	x, y, z	= point.x, point.y, point.z
	x′, y′, z′ = loop.center.x, loop.center.y, loop.center.z
	z, r, θ = cartesian2cylindrical(x - x′, y - y′, z - z′)

	# Compute the field strength at z, r (using the fact that the loop is centered at z = 0, r= 0)
	I = loop.current
	rL = loop.radius

	Rmax² = z^2 + (r+rL)^2
	Rmin² = z^2 + (r-rL)^2
	Rmax = sqrt(Rmax²)
	Rmin = sqrt(Rmin²)

	k² = 4*r*rL / Rmax²
	ξ = atan(z, abs(r-rL))
	kp² = 1-k²

	K, E = Elliptic.ellipke(k²)
	Ei = Elliptic.E(ξ, kp²)
	Fi = Elliptic.F(ξ, kp²)
	Λ = 2/pi * (E*Fi + K*(Ei - Fi))

	Ω = 0.0
	# solid angle subtended by current loop at (z,r)
	if r < rL
		Ω = 2*π - 2*z/Rmax*K - π*Λ
	elseif r == rL
		Ω = π - 2*z/Rmax*K
	else
		Ω = - 2*z/Rmax*K + π*Λ
	end

	B₀ = μ₀ * I / 2/ rL;

	Bz = μ₀ * I / π ./ (2 * Rmin² * Rmax) * ((rL^2 - r^2 - z^2)*E + Rmin² * K)
	Br = μ₀ * I * z /π ./ (2 * Rmin² * Rmax * r) * ((rL^2 + r^2 + z^2)*E - Rmin² * K)
	λ = B₀*rL/2/pi * Rmax * ((2-k²)*K - 2*E)
	χ = I/4/pi * Ω

	if isnan(Br) || isinf(Br)
		Br = 0
	end

	if isnan(Bz) || isinf(Bz)
		Bz = 0
	end

	# Convert Bz and Br to Bx, By, Bz
	Bx, By, Bz = cylindrical2cartesian(Bz, Br, θ)

	return Bx, By, Bz, λ, χ
end

function field_on_grid(loops, xs, ys, zs)
	Nx, Ny, Nz = length(xs), length(ys), length(zs)

	Bx = zeros(Nx, Ny, Nz)
	By = zeros(Nx, Ny, Nz)
	Bz = zeros(Nx, Ny, Nz)

	for (i, x) in enumerate(xs), (j, y) in enumerate(ys), (k, z) in enumerate(zs)
		for loop in loops
			bx, by, bz, _, _ = field_at_point(loop, CartesianPoint(x, y, z))
			Bx[i, j, k] += bx
			By[i, j, k] += by
			Bz[i, j, k] += bz
		end
	end

	return Bx, By, Bz
end

"""
	write_field(filename, xs, ys, Bx, By)
Write magnetic field to Tecplot ASCII datafile for two-dimensional data
"""
function write_field(filename, xs, ys, Bx, By)
	split_fname = splitpath(filename)
	if length(split_fname) > 1 # If file is in folder, make sure that folder exists
		mkpath(joinpath(split_fname[1:end-1]...))
	end
	Nx, Ny= length(xs), length(ys)
	open(filename, "w") do f
		println(f, "VARIABLES = \"X\" \"Y\" \"Bx\" \"By\"")
		println(f, "ZONE I = $Nx, J = $Ny")
		for (j, y) in enumerate(ys), (i, x) in enumerate(xs)
			row = "$(x)\t$(y)\t$(Bx[i, j])\t$(By[i, j])"
			println(f, row)
		end
	end
end
"""

	write_field(filename, xs, ys, zs, Bx, By, Bz)
Write magnetic field to Tecplot ASCII datafile for three-dimensional data
"""
function write_field(filename, xs, ys, zs, Bx, By, Bz)
	split_fname = splitpath(filename)
	if length(split_fname) > 1 # If file is in folder, make sure that folder exists
		mkpath(joinpath(split_fname[1:end-1]...))
	end
	Nx, Ny, Nz = length(xs), length(ys), length(zs)
	open(filename, "w") do f
		println(f, "VARIABLES = \"X\" \"Y\" \"Z\" \"Bx\" \"By\" \"Bz\"")
		println(f, "ZONE I = $Nx, J = $Ny, K = $Nz")
		for (k, z) in enumerate(zs), (j, y) in enumerate(ys), (i, x) in enumerate(xs)
			row = "$(x)\t$(y)\t$(z)\t$(Bx[i, j, k])\t$(By[i, j, k])\t$(Bz[i, j, k])"
			println(f, row)
		end
	end
end

#=
function Λ₀(ξ::Float64, k::Float64)
	kprime = 1 - k
	K, E = Elliptic.ellipke(k)
	return 2/π * (E * Elliptic.F(ξ, kprime) +
				  K * (Elliptic.E(ξ, kprime) - Elliptic.F(ξ, kprime)))
end

function fieldatpoint(z, r, rL, I)
	Bz = 0.0
	Br = 0.0
	λ = 0.0
	χ = 0.0

	Rmax² = z^2 + (r+rL)^2
	Rmin² = z^2 + (r-rL)^2
	Rmax = sqrt(Rmax²)
	Rmin = sqrt(Rmin²)

	k² = 4*r*rL / Rmax²
	ξ = atan(z, abs(r-rL))
	kp² = 1-k²

	K, E = Elliptic.ellipke(k²)
	Ei = Elliptic.E(ξ, kp²)
	Fi = Elliptic.F(ξ, kp²)
	Λ = 2/pi * (E*Fi + K*(Ei - Fi))

	Ω = 0.0
	# solid angle subtended by current loop at (z,r)
	if r<rL
		Ω = 2*π - 2*z/Rmax*K - π*Λ
	elseif r==rL
		Ω = π - 2*z/Rmax*K
	else
		Ω = - 2*z/Rmax*K + π*Λ
	end

	μ₀ = VACUUM_PERMEABILITY

	B₀ = μ₀ * I / 2/ rL;

	Bz = μ₀ * I / π ./ (2 * Rmin² * Rmax) * ((rL^2 - r^2 - z^2)*E + Rmin² * K)
	Br = μ₀ * I * z /π ./ (2 * Rmin² * Rmax * r) * ((rL^2 + r^2 + z^2)*E - Rmin² * K)
	λ = B₀*rL/2/pi * Rmax * ((2-k²)*K - 2*E)
	χ = I/4/pi * Ω

	if isnan(Br) || isinf(Br)
		Br = 0
	end

	if isnan(Bz) || isinf(Bz)
		Bz = 0
	end

	return Bz, Br, λ, χ
end

#=function loopstreamfunc(z, r, rL, I)
	Rmax² = z^2 + (r+rL)^2
	Rmax = sqrt(Rmax²)
    k² = 4*r*rL / Rmax²
	B₀ = VACUUM_PERMEABILITY * I / 2/ rL

    E = Elliptic.E(k²)
	K = Elliptic.K(k²)

	λ = B₀*rL/2/pi * Rmax * ((2-k²)*K - 2*E)
    return λ
end

function loopscalarpotential(z, r, rL, I)
	Rmax² = z^2 + (r+rL)^2
	k² = 4*r*rL / Rmax²
	ξ = atan(z, abs(r-rL))

	Rmax = sqrt(Rmax²)

	Ω = 0.0

	if r<rL
		Ω = 2*π - 2*z/Rmax*Elliptic.K(k²) - π*Λ₀(ξ, k²)
	elseif r==rL
		Ω = π - 2*z/Rmax*Elliptic.K(k²)
	else
		Ω = - 2*z/Rmax*Elliptic.K(k²) + π*Λ₀(ξ, k²)
	end

	return I/4/π * Ω
end

function loopscalarpotential(points, loop::CurrentLoop)
	@views χ = loopscalarpotential.(
		points[:, 1] .- loop.z_coord,
		points[:, 2],
		loop.radius,
		loop.current,
	)
	return χ
end

function loopstreamfunc(points, loop::CurrentLoop)
	@views λ = loopstreamfunc.(
		points[:, 1] .- loop.z_coord,
		points[:, 2],
		loop.radius,
		loop.current,
	)
	return λ
end=#




function calculatefieldonmesh(zgrid, rgrid, loop::CurrentLoop)
	Bz = zero(zgrid)
	Br = zero(zgrid)
	λ = zero(zgrid)
	χ = zero(zgrid)

	for i in eachindex(zgrid)
		Bz[i], Br[i], λ[i], χ[i] = fieldatpoint.(
			zgrid[i] .- loop.z_coord,
			rgrid[i],
			loop.radius,
			loop.current,
		)
	end

	return Bz, Br, λ, χ
end

function extractfieldline(z, r, zgrid, rgrid, λgrid)
	zvals = zgrid[1, :]
	rvals = rgrid[:, 1]
	zind = findfirst(zval -> zval > z, zvals)
	rind = findfirst(rval -> rval > r, rvals)

	if zind == 1 && rind == 1
		λT = λ[1,1]
	elseif zind == 1
		λ_bot = λgrid[rind - 1, 1]
		λ_top = λgrid[rind, 1]
		rwt = (r - rgrid[rind - 1, 1]) / (rgrid[rind, 1] - rgrid[rind - 1, 1])
		λT = (λ_top - λ_bot) * rwt + λ_bot
	elseif rind == 1
		λ_left = λgrid[1, zind - 1]
		λ_right = λgrid[1, zind]
		zwt = (z - zgrid[1, zind - 1]) / (zgrid[1, zind] - zgrid[1, zind - 1])
		λT = (λ_right - λ_left) * zwt + λ_left
	else
		λ_br = λgrid[rind - 1, zind]
		λ_tr = λgrid[rind, zind]
		λ_bl = λgrid[rind - 1, zind - 1]
		λ_tl = λgrid[rind, zind - 1]
		zwt = (z - zgrid[1, zind - 1]) / (zgrid[1, zind] - zgrid[1, zind - 1])
		rwt = (r - rgrid[rind - 1, 1]) / (rgrid[rind, 1] - rgrid[rind - 1, 1])

		λ_top = zwt * (λ_tr - λ_tl) + λ_tl
		λ_bot = zwt * (λ_br - λ_bl) + λ_bl

		λT = (λ_top - λ_bot) * rwt + λ_bot
	end

	throatline = Contour.contour(zgrid, rgrid, λgrid, λT)
	zs, rs = coordinates(lines(throatline)[1])
	return [zs rs]
end

function save_field(z, r, Bz, Br, path)
	ni, nj = size(Bz)
	#for i in eachindex(Br)
	#	if Br[i] == 0.0
	#		Br[i] = 1e-16
	#	end
	#end
	open(path, "w") do io
		write(io, "VARIABLES = \"z\" \"r\" \"Bz\" \"Br\"\n")
		zonestr = @sprintf "ZONE I=%5d, J=%5d, F = POINT\n" ni nj
		write(io, zonestr)
		for j in 1:nj, i in 1:ni
			line = @sprintf "%e\t%e\t%e\t%e\n" z[i, j] r[i, j] Bz[i, j] Br[i, j]
			write(io, line)
		end
	end
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
end=#

end
