using SpatInterp, Test
import SpatRasters
using SpatRasters: bbox, make_rast, st_dims
using RTableTools, Statistics


begin
  indir = "$(@__DIR__)/.." |> abspath
  d = fread("$indir/data/prcp_st174_shiyan.csv")

  X = [d.lon d.lat] # d.alt
  Y = d.prcp
  # Y = repeat(y, outer=(1, 24 * 30))
  b = bbox(109.5, 31.5, 112 - 0.5, 33.5)
  target = make_rast(; b, cellsize=0.01)
end
