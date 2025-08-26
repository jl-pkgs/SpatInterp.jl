using SpatRasters: bbox, make_rast, st_dims
using Statistics
using RTableTools

@testset "find_neighbor" begin
  indir = "$(@__DIR__)/.." |> abspath
  d = fread("$indir/data/prcp_st174_shiyan.csv")

  X = [d.lon d.lat] # d.alt
  Y = d.prcp
  # Y = repeat(y, outer=(1, 24 * 30))
  b = bbox(109.5, 31.5, 112 - 0.5, 33.5)
  target = make_rast(; b, cellsize=0.01)

  neighbor = find_neighbor(target, X; nmax=24, radius=20, do_angle=true)
  @test mean(neighbor.count) >= 5

  neighbor = find_neighbor(target, X; nmax=24, radius=50, do_angle=true)
  @test mean(neighbor.count) >= 19
end

# # mean(neighbor.angle)
# using GLMakie, MakieLayers
# neighbor = find_neighbor(target, X; nmax=20, radius=20, do_angle=true)

# begin
#   lon, lat = st_dims(target)
#   imagesc(lon, lat, neighbor.count)
# end
