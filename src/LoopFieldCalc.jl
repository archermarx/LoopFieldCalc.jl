module LoopFieldCalc

using Elliptic, Printf
include("meshtools.jl")
export 	CurrentLoop,
		calculatefieldonmesh,
		loopstreamfunc,
		loopscalarpotential,
		extractfieldline,
		save_field,
		B0

const VACUUM_PERMITTIVITY = 4 * pi * 1e-7

struct CurrentLoop{T <: Real}
	radius::T
	current::T
	z_coord::T
end

function CurrentLoop(radius::T, current::T; z_coord::T = zero(T)) where {T<:Real}
	return CurrentLoop(radius, current, z_coord)
end

function CurrentLoop(radius::T; Bfield::T, z_coord::T = zero(T)) where {T<:Real}
	current = Bfield / VACUUM_PERMITTIVITY * 2 * radius
	return CurrentLoop(radius, current, z_coord)
end

B0(loop) = VACUUM_PERMITTIVITY * loop.current / 2 / loop.radius;

function loopstreamfunc(z, r, rL, I)
	Rmax² = z^2 + (r+rL)^2
	Rmax = sqrt(Rmax²)
    k² = 4*r*rL / Rmax²
	B₀ = VACUUM_PERMITTIVITY * I / 2/ rL

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
end

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

	μ₀ = VACUUM_PERMITTIVITY

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

function extractfieldline(z, r, zgrid, rgrid, λgrid, loop::CurrentLoop)
	λT = loopstreamfunc([z r], loop)
	throatline = Contour.contour(zgrid, rgrid, λgrid, λT[1])
	zs, rs = coordinates(lines(throatline)[1])
	return [zs rs]
end

function save_field(z, r, Bz, Br, path)
	ni, nj = size(Bz)
	open(path, "w") do io
		write(io, "VARIABLES = \"z\" \"r\" \"Bz\" \"Br\"\n")
		zonestr = @sprintf "ZONE I=%5d, J=%5d, F = POINT\n" ni nj
		write(io, zonestr)
		for  i in 1:ni, j in 1:nj
			line = @sprintf "%f\t%f\t%f\t%f\n" z[i, j] r[i, j] Bz[i, j] Br[i, j]
			write(io, line)
		end
	end
end

end
