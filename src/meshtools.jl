using Contour, GeometricalPredicates
import GeometricalPredicates.Polygon

function Polygon(pts::Array{T, 2}) where T<:Real
	npts = size(pts, 1)
	pointarr = Array{Point2D, 1}(undef, npts)
	for p in 1:npts
		pointarr[p] = Point(pts[p, 1], pts[p, 2])
	end
	return Polygon(pointarr...)
end

function replace_in_polygon(xgrid, ygrid, fields, points, replace_val)
	polygon = Polygon(points)
	for i in eachindex(xgrid)
		point = Point(xgrid[i], ygrid[i])
		if inpolygon(polygon, point)
			for f in fields
				f[i] = NaN
			end
		end
	end
end
