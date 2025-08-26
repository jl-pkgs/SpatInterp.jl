export weight_idw!


function weight_idw!(neighbor::Neighbor{FT,N}; m::Int=2) where {FT,N}
  (; count, dist, weight, dims) = neighbor
  weight .= FT(0.0)
  nlon, nlat = dims[1:2]
  for i in 1:nlon, j in 1:nlat
    n = count[i, j]
    # _w = @view weight[i, j, 1:n]
    _dist = @view dist[i, j, 1:n]
    wk = @.(1 / _dist^m)
    wk .= wk ./ sum(wk)
    weight[i, j, 1:n] .= wk
  end
end
