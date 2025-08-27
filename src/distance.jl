function distance_norm(p1::AbstractArray{FT}, p2::AbstractArray{FT}; kw...) where {FT<:AbstractFloat}
  # @assert length(p1) == length(p2)
  dist2 = FT(0.0)
  @inbounds for i in eachindex(p1)
    dist2 += (p1[i] - p2[i])^2
  end
  sqrt(dist2)
end
