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
