function interp_weight(neighbor::Neighbor{T1,N}, y::AbstractArray{T2}, target::SpatRaster;
  progress=true, ignored...) where {T1, T2<:Real,N}

  ntime = size(y, 2)
  lon, lat = st_dims(target)
  nlon, nlat = length(lon), length(lat)
  
  # 综合考虑 T1 和 T2，选择合适的计算类型
  T = promote_type(T1, T2)
  R = zeros(T, nlon, nlat, ntime)

  (; count, index, weight) = neighbor
  ∅ = zero(T)
  p = Progress(nlat)
  @inbounds @threads for j in 1:nlat
    progress && next!(p)
    for i in 1:nlon

      n_control = count[i, j]
      inds = @view index[i, j, 1:n_control]
      ws = @view weight[i, j, 1:n_control]

      for k in 1:ntime
        ∑ = zero(T)
        ∑w = zero(T)

        for l in 1:n_control # control points
          _i = inds[l]
          yᵢ = y[_i, k]
          notnan = yᵢ == yᵢ
          w_val = ws[l]
          ∑ += ifelse(notnan, yᵢ * w_val, ∅)
          ∑w += ifelse(notnan, w_val, ∅)
        end
        R[i, j, k] = ∑ / ∑w
      end
    end
  end
  rast(R, target)
end

export interp_weight
