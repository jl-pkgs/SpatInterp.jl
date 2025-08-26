

# points: 2×N，每列一个点；q: 长度2；k: 近邻数（默认32）
function quad32(points::AbstractMatrix, q::AbstractVector; k::Int=32)
  tree = KDTree(points)
  kk = min(k, size(points, 2))
  idxs, dists = knn(tree, q, kk, true)
  if !isempty(dists) && dists[1] == 0.0           # q恰是样本点，去掉自身
    idxs, dists = idxs[2:end], dists[2:end]
  end
  P = points[:, idxs]
  dx, dy = P[1, :] .- q[1], P[2, :] .- q[2]
  quads = Dict(
    1 => idxs[(dx.>=0).&(dy.>=0)],        # 第一象限（含轴）
    2 => idxs[(dx.<0).&(dy.>=0)],
    3 => idxs[(dx.<0).&(dy.<0)],
    4 => idxs[(dx.>=0).&(dy.<0)],
  )
  return (; quads, idxs, dists)
end
