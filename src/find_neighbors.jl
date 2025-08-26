export Neighbor, find_neighbor


@with_kw struct Neighbor{FT,N}
  nmax::Int = 24
  dims::Tuple = ()

  source::Union{Nothing,Vector{Point{2,FT}}} = nothing # donor points
  lon::AbstractVector{FT} = 1:dims[1]                # target raster
  lat::AbstractVector{FT} = 1:dims[2]

  count::AbstractArray{Integer} = zeros(Int, dims...)
  index::AbstractArray{Integer,N} = zeros(Int, dims..., nmax)

  "distance (km)"
  dist::AbstractArray{FT,N} = zeros(FT, dims..., nmax)

  "azimuth angle (in radian, degree, 0-360), N: 0°, E: 90°"
  angle::AbstractArray{FT,N} = zeros(FT, dims..., nmax)

  weight::AbstractArray{FT,N} = zeros(FT, dims..., nmax)
end

function Neighbor(nmax::Int, dims::Tuple=(); FT=Float64, source=nothing, kw...)
  N = length(dims) + 1
  Neighbor{FT,N}(; source, nmax, dims, kw...)
end


function find_neighbor(point, tree; nmax::Int=20, radius::Real=100)
  inds, dists = knn(tree, point, nmax)
  I = dists .<= radius
  inds = @view inds[I]
  dists = @view dists[I]
  (; index=inds, distance=dists)
end


"""
    find_neighbor(ra::SpatRaster, X; nmax::Int=20, radius::Real=100, do_angle=true) 

已经按照距离进行了排序。
"""
function find_neighbor(ra::SpatRaster{FT}, X;
  fast=false, nmax::Int=20, radius::Real=100, do_angle=true) where {FT<:Real}

  Xt = collect(X')
  if fast
    tree = KDTree(Xt; reorder=true, leafsize=128*2)
    printstyled("Note: radius should be in degree in fast mode!\n", color=:blue, bold=true)
  else
    tree = BallTree(Xt, Haversine(6371.0); leafsize=64 * 1, reorder = true)
  end

  lon, lat = st_dims(ra)
  lon, lat = FT.(lon), FT.(lat)

  nlon, nlat = length(lon), length(lat)
  neighbor = Neighbor(nmax, (nlon, nlat); source=st_points(X), lon, lat, FT)
  (; count, index, dist, angle) = neighbor

  p = Progress(nlat)
  @inbounds @threads for j in 1:nlat
    next!(p)
    for i in 1:nlon

      p0 = [lon[i], lat[j]]
      _inds, _dists = knn(tree, p0, nmax) ## 这里可以更快
      _I = sortperm(_dists) # 已经排序
      _I = _I[_dists[_I].<=radius]

      _inds = @view _inds[_I]
      _dists = @view _dists[_I]

      n_control = length(_dists)
      n_control == 0 && continue

      count[i, j] = n_control
      index[i, j, 1:n_control] .= _inds
      dist[i, j, 1:n_control] .= _dists

      if do_angle
        for c in 1:n_control
          _i = _inds[c]
          p1 = @view X[_i, :]
          angle[i, j, c] = angle_azimuth_sphere(p0, p1; to_degree=true)
        end
      end

    end
  end
  neighbor
end
