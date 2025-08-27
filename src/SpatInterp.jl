module SpatInterp

using GeometryBasics
using Parameters
using Distances: Haversine
using NearestNeighbors
using Base.Threads
using ProgressMeter
using DocStringExtensions

using SpatRasters: SpatRaster, rast,
  st_coords, st_dims, meshgrid,
  bbox, st_bbox, make_rast

export SpatRaster, rast, bbox, make_rast
export Point, st_points

st_points(X::AbstractMatrix{T}) where {T} = map(p -> Point{2,T}(p[1], p[2]), eachrow(X))

# @inline Base.getproperty(p::Point{T}, ::Symbol{Val{:x}}) where {T} = p[1]
# @inline function Base.getproperty(p::Point{T}, s::Symbol) where {T}
#   s === :x && return p[1]
#   s === :y && return p[2]
#   s === :z && return p[3]
#   return getfield(p, s)  # fallback for other fields
# end

include("angle.jl")
include("distance.jl")
include("find_neighbors.jl")
include("weights.jl")

include("bilinear/bilinear_helper.jl")
include("bilinear/bilinear.jl")
include("bilinear/bilinear_irregular.jl")
include("bilinear/frac_dist.jl")
include("bilinear/find_quad.jl")

include("interp.jl")
include("interp_tps.jl")

# include("utilize.jl")
# include("get_fractional.jl")
# include("solve_quadratic.jl")
# include("get_fractional.jl")

end # module SpatInterp
