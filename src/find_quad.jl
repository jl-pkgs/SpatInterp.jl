export angle2quad, find_quad

# --- 角度→象限（度制，0°含北，90°含东，半开区间） ---
# 映射： [0,90)→q1(NE), [90,180)→q4(SE), [180,270)→q3(SW), [270,360)→q2(NW)
# 这样 0° 归 NE，90° 归 SE，180° 归 SW，270° 归 NW（轴点各归顺时针的那个象限）
@inline function angle2quad(θ::Real)
  s = fld(mod(float(θ), 360.0), 90.0)        # 0,1,2,3
  return (2, 4, 3, 1)[Int(s)+1]              # 1→1, 2→4, 3→3, 4→2
end

# nearest_per_quadrant
function find_quad(
  source::AbstractVector{Point{2,T}}, out_x::T, out_y::T,
  index::AbstractVector, angle::V, dist::V, weight::V; m::Int=2) where {
  T<:Real,V<:AbstractVector{T}}

  idxq = zeros(Int, 4)         # 没有点 → 0
  distq = fill(T(Inf), 4)     # 没有点 → Inf
  angq = fill(T(NaN), 4)      # 没有点 → NaN
  wq = zeros(T, 4)            # 没有点 → 0.0

  @inbounds for j in eachindex(index)
    q = angle2quad(angle[j])
    dj = dist[j]

    if dj < distq[q]
      distq[q] = dj
      idxq[q] = index[j]
      angq[q] = angle[j]
      wq[q] = weight[j]
    end
  end
  
  n = sum(idxq .!== 0)
  t = s = T(0)
  if n == 4
    # 采用bilinear, 计算权重，点的顺序要对
    _points = @view source[idxq]
    t, s = frac_dist(_points..., out_x, out_y) ## 之后，转换为权重
    wq[1] = (1 - s) * (1 - t)
    wq[2] = s * (1 - t)
    wq[3] = (1 - s) * t
    wq[4] = s * t
    # if sum(isnan.(wq)) >=2 
    #   println("$out_x, $out_y")
    #   jldsave("debug.jld2"; _points, out_x, out_y, t, s, wq)
    #   error("stop here")
    # end
  end

  if n < 4 || isnan(t + s) #
    wq = @.(1 / distq^m) # 不足4点，或bilinear求解失败时; 退化为idw插值
    wq .= wq ./ sum(wq)
  end

  # 最后需要进行排序
  if n < 4
    I = sortperm(idxq, rev=true) # 查找到的在前
    idxq = idxq[I]
    distq = distq[I]
    angq = angq[I]
    wq = wq[I]
  end
  return idxq, distq, angq, wq  # 长度为4
end


function find_quad(neighbor::Neighbor{FT,3}; progress::Bool=true) where {FT}
  (; source, lon, lat) = neighbor
  neighbor4 = Neighbor(4, neighbor.dims; FT, source, lon, lat) # 4象限邻居
  (; count, index, dist, angle, weight) = neighbor4
  nlon, nlat = size(count)

  p = Progress(nlat)

  @views @inbounds for j in 1:nlat
    progress && next!(p)
    for i in 1:nlon

      n = neighbor.count[i, j]
      n == 0 && continue

      out_x, out_y = lon[i], lat[j]

      _index = neighbor.index[i, j, 1:n]
      _angle = neighbor.angle[i, j, 1:n]
      _dist = neighbor.dist[i, j, 1:n]
      _weight = neighbor.weight[i, j, 1:n]

      idxq, distq, angq, wq = find_quad(source, out_x, out_y, _index, _angle, _dist, _weight)

      count[i, j] = sum(idxq .!== 0)
      index[i, j, :] .= idxq
      angle[i, j, :] .= angq
      dist[i, j, :] .= distq
      weight[i, j, :] .= wq
    end
  end
  neighbor4
end
