module LoopFieldCalc

using Elliptic, Printf, Contour
export 	CurrentLoop,
		calculatefieldonmesh,
		loopstreamfunc,
		loopscalarpotential,
		extractfieldline,
		save_field,
		B0

VACUUM_PERMEABILITY = 4 * pi * 1e-7

struct CurrentLoop{T <: Real}
	radius::T
	current::T
	z_coord::T
end

function CurrentLoop(radius::T, current::T; z_coord::T = zero(T)) where {T<:Real}
	return CurrentLoop(radius, current, z_coord)
end

function CurrentLoop(radius::T; Bfield::T, z_coord::T = zero(T)) where {T<:Real}
	current = Bfield / VACUUM_PERMEABILITY * 2 * radius
	return CurrentLoop(radius, current, z_coord)
end

B0(loop) = VACUUM_PERMEABILITY * loop.current / 2 / loop.radius;

function loopstreamfunc(z, r, rL, I)
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
end

end
