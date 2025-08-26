using NetCDFTools
using GLMakie, MakieLayers
# using SpatRasters: 

f = "data/FY4B_20250714014500_20250714015959.nc"
# lon, lat = st_dims(f)

I, J = nc_read(f, "i"), nc_read(f, "j")
lon = nc_read(f, "lon")
lat = nc_read(f, "lat")
P = nc_read(f, "P")



begin
  fig = Figure(;size=(1500, 600))
  imagesc!(fig[1, 1], I, J, lon; title = "lon")
  imagesc!(fig[1, 2], I, J, lat; title="lat")
  imagesc!(fig[1, 3], I, J, P; title="P")
  fig
end
