export interp

"""
    interp(ra::SpatRaster, target::SpatRaster;
        nmax::Int=20, radius::Real=100, do_angle=false, method="adw", kw...)
    
    interp(x::AbstractMatrix, y::AbstractArray{FT}, target::SpatRaster;
        nmax::Int=20, radius::Real=200, do_angle=false, method="adw", kw...)
    
    interp(neighbor::Neighbor{FT,N}, y::AbstractArray{FT}, target::SpatRaster; ignored...)

## Arguments
- `method`: one of "adw", "idw", "nearest", "bilinear"

- `kw`: other parameters to `wfun`
  + `weight_idw!`: `m`
  + `weight_adw!`: `cdd`, `m`
"""
function interp(ra::SpatRaster, target::SpatRaster;
  nmax::Int=20, radius::Real=100, do_angle=false, method="adw", kw...)
  dims = size(ra)
  nlon, nlat = dims[1:2]
  ntime = length(dims) >= 3 ? dims[3] : 1

  if method == "bilinear"
    x, y = st_dims(ra)
    xx, yy = st_dims(target)
    Z = ra.A
    R = bilinear(x, y, Z, xx, yy)
    rast(R, st_bbox(target))
  else
    X = st_coords(ra) # [lon, lat]
    Y = reshape(ra.A, nlon * nlat, ntime)
    interp(X, Y, target; nmax, radius, do_angle, method, kw...)
  end
end


function interp(x::AbstractMatrix, y::AbstractArray{FT}, target::SpatRaster;
  nmax::Int=20, radius::Real=100, do_angle=false,
  method="adw", kw...) where {FT}

  neighbor = find_neighbor(target, x; nmax, radius, do_angle)

  if method == "nearest"
    interp_nearest(neighbor, y)
  else
    wfun = WFUNS[method]
    wfun(neighbor; kw...)
    interp_weight(neighbor, y)
  end
end

function interp_nearest(neighbor::Neighbor{FT,N}, y::AbstractArray{FT};
  progress=true, ignored...) where {FT,N}
  ntime = size(y, 2)
  (; lon, lat) = neighbor
  nlon, nlat = length(lon), length(lat)
  b = st_bbox(lon, lat)
  R = zeros(FT, nlon, nlat, ntime)

  (; count, index) = neighbor

  p = Progress(nlat)
  @inbounds for j in 1:nlat
    progress && next!(p)
    for i in 1:nlon

      n_control = count[i, j]
      inds = @view index[i, j, 1:n_control]

      @threads for k in 1:ntime
        for l in 1:n_control # control points
          _i = inds[l]
          yᵢ = y[_i, k]
          notnan = yᵢ == yᵢ
          if notnan # 最邻近原则，只要不是Nan，就可以跳出
            R[i, j, k] = yᵢ
            break # 跳出for循环
          end
        end # end n_control
      end # end time
    end # end lat
  end # end lon
  rast(R, b)
end



# - `Y`: [n_source, ntime]
function interp_weight(neighbor::Neighbor{T,N}, Y::AbstractArray{T};
  progress=true, ignored...) where {T<:AbstractFloat,N}
  
  (; lon, lat) = neighbor
  nlon, nlat = length(lon), length(lat)
  ntime = size(Y, 2)
  b = st_bbox(lon, lat)
  
  # 综合考虑 T1 和 T2，选择合适的计算类型
  # T = promote_type(T1, T2)
  R = zeros(T, nlon, nlat, ntime)

  (; count, index, weight) = neighbor
  p = Progress(nlat)

  @inbounds @threads for j in 1:nlat
    progress && next!(p)
    for i in 1:nlon

      n_control = count[i, j]
      inds = @view index[i, j, 1:n_control]
      ws = @view weight[i, j, 1:n_control]
      _Y = @view Y[inds, :] # target
      _Z = @view R[i, j, :]
      _weighted_mean!(_Z, ws, _Y)
    end
  end
  rast(R, b)
end


"""
- `Z` : [ntime]
- `Y` : [n_control, ntime]
- `ws`: [n_control]
"""
function _weighted_mean!(Z::AbstractVector{T}, ws::AbstractVector{T}, Y::AbstractMatrix{T}) where {T<:AbstractFloat}
  n_control, ntime = size(Y)
  # Φ::T = T(0)
  @inbounds @simd for k in 1:ntime
    ∑::T = T(0)
    ∑w::T = T(0)

    for l in 1:n_control # control points
      yᵢ::T = Y[l, k]
      w::T = ws[l]
      
      mask = yᵢ == yᵢ  # 非NaN检查
      ∑w += ifelse(mask, w, T(0))
      ∑ = ifelse(mask, muladd(yᵢ, w, ∑), ∑)
    end
    Z[k] = ifelse(∑w > T(0), ∑ / ∑w, T(0))
  end
end

export interp_weight
