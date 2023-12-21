module LoopFieldCalc

using Elliptic2, Printf, Contour

const μ₀ = π * 4e-7
const μ0 = μ₀

include("geometry.jl")

struct CurrentLoop{R, I, X, Y, Z}
	radius::R
	current::I
	center::CartesianPoint{X, Y, Z}
end

CurrentLoop(center::CartesianPoint; radius, current) = CurrentLoop(radius, current, center)
CurrentLoop(center::CylindricalPoint; radius, current) = CurrentLoop(radius, current, center)

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

	K, E = Elliptic2.ellipke(k²)
	Ei = Elliptic2.E(ξ, kp²)
	Fi = Elliptic2.F(ξ, kp²)
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
	streamfunc = zeros(Nx, Ny, Nz)
	potential = zeros(Nx, Ny, Nz)

	for (i, x) in enumerate(xs), (j, y) in enumerate(ys), (k, z) in enumerate(zs)
		for loop in loops
			bx, by, bz, s, p = field_at_point(loop, CartesianPoint(x, y, z))
			Bx[i, j, k] += bx
			By[i, j, k] += by
			Bz[i, j, k] += bz
			streamfunc[i, j, k] += s
			potential[i, j, k] += p
		end
	end

	return Bx, By, Bz, streamfunc, potential
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

end
