using NetCDFTools: nc_read, st_dims
using GLMakie, MakieLayers
using SpatInterp
using SpatRasters: SpatRaster, bbox, make_rast, bbox2dims, rast

import MakieLayers: imagesc, imagesc!
function imagesc!(handle, ra::SpatRaster, args...; kw...)
  lon, lat = st_dims(ra)
  imagesc!(handle, lon, lat, ra.A, args...; kw...)
end

f = "data/FY4B_20250714014500_20250714015959.nc"
# lon, lat = st_dims(f)

I, J = nc_read(f, "i"), nc_read(f, "j")
lon = nc_read(f, "lon")
lat = nc_read(f, "lat")
Prcp = nc_read(f, "P")

# 采用bilinear_irregular，对数据进行插值
inds = findall(.!isnan.(lon[:]))
X = [lon[inds] lat[inds]] # 76.6%, lon, lat
Z = Prcp[inds]

function make_rast2(; b::bbox, cellsize, FT=Float64)
  lon, lat = bbox2dims(b; cellsize)
  nlon, nlat = length(lon), length(lat)
  rast(zeros(FT, nlon, nlat), b)
end

b = bbox(70., 15., 140., 55.)
target = make_rast2(; b, cellsize=1/120*4, FT=Float32)

# radius = 4km * 4
@profview neighbor = find_neighbor(target, X; nmax=24, radius=4*4, do_angle=true)
neighbor4 = find_quad(neighbor)

@profview res = interp_weight(neighbor4, Z, target)


fig = Figure(; size=(1500, 600))
imagesc!(fig, res)
fig

begin
  fig = Figure(; size=(1500, 600))
  imagesc!(fig[1, 1], I, J, lon; title="lon")
  imagesc!(fig[1, 2], I, J, lat; title="lat")
  imagesc!(fig[1, 3], I, J, P; title="P")
  fig
end
