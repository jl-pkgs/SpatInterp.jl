module Bilinear

using GeometryBasics
using Parameters
using SpatRasters: SpatRaster, st_dims
using Distances: Haversine
using NearestNeighbors


export Point

# @inline Base.getproperty(p::Point{T}, ::Symbol{Val{:x}}) where {T} = p[1]
# @inline function Base.getproperty(p::Point{T}, s::Symbol) where {T}
#   s === :x && return p[1]
#   s === :y && return p[2]
#   s === :z && return p[3]
#   return getfield(p, s)  # fallback for other fields
# end


include("bilinear_irregular.jl")
include("angle.jl")
include("find_neighbors.jl")

# include("utilize.jl")
# include("get_fractional.jl")
# include("solve_quadratic.jl")
# include("get_fractional.jl")

end # module Bilinear
