abstract type Point end

Base.@kwdef struct CartesianPoint{X, Y, Z} <: Point
	x::X
	y::Y
	z::Z
end

Base.@kwdef struct CylindricalPoint{Z, R, T} <: Point
	z::Z
	r::R
	θ::T
end

"""
	cartesian2cylindrical(x, y, z)
Transform a point or vector from cartesian (x, y, z)  to cylindrical (z, r, θ) coordinates,
assuming the z-axis is the axis of the cylinder
"""
function cartesian2cylindrical(x, y, z)
	r = √(x^2 + y^2)
	θ = atan(y, x)
	return (z, r, θ)
end

"""
	cylindrical2cartesian(z, r, θ)
Transform a point or vector (z, r, θ) from cylindrical  to cartesian (x, y, z) coordinates,
assuming the z-axis is the axis of the cylinder
"""
function cylindrical2cartesian(z, r, θ)
	x = r * cos(θ)
	y = r * sin(θ)
	return (x, y, z)
end

function cartesian2cylindrical(cart::CartesianPoint)
	z, r, θ = cartesian2cyclindrical(cart.x, cart.y, cart.z)
	return CylindricalPoint(z, r, θ)
end

function cylindrical2cartesian(cyl::CylindricalPoint)
	x, y, z = cylindrical2cartesian(cyl.z, cyl.r, cyl.θ)
	return CartesianPoint(x, y, z)
end
