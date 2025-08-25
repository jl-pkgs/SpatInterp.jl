using NearestNeighbors
using LinearAlgebra

# Base class for bilinear interpolation
mutable struct BilinearBase
  bilinear_t::Union{Nothing,Vector{Float64}}
  bilinear_s::Union{Nothing,Vector{Float64}}
  slices_x::Union{Nothing,Matrix{Int32}}
  slices_y::Union{Nothing,Matrix{Int32}}
  mask_slices::Union{Nothing,BitMatrix}
  out_coords_x::Union{Nothing,Vector{Float64}}
  out_coords_y::Union{Nothing,Vector{Float64}}

  _valid_input_index::Union{Nothing,BitVector}
  _index_array::Union{Nothing,Matrix{Int32}}
  _distance_array::Union{Nothing,Matrix{Float64}}
  _neighbours::Int
  _epsilon::Float64
  _reduce_data::Bool
  _radius_of_influence::Float64
  _resample_kdtree::Union{Nothing,KDTree}
  _target_lons::Union{Nothing,Array{Float64}}
  _target_lats::Union{Nothing,Array{Float64}}

  function BilinearBase(source_geo_def, target_geo_def, radius_of_influence;
    neighbours=32, epsilon=0.0, reduce_data=true)
    new(nothing, nothing, nothing, nothing, nothing, nothing, nothing,
      nothing, nothing, nothing, neighbours, epsilon, reduce_data,
      radius_of_influence, nothing,
      nothing, nothing)
  end
end

function resample!(resampler::BilinearBase, data; fill_value=nothing, nprocs=1)
  throw(ErrorException("Not implemented - should be overridden by subclass"))
end

function get_bil_info!(resampler::BilinearBase, source_geo_def, target_geo_def; kdtree_class=KDTree, nprocs=1)
  """Calculate bilinear neighbour info."""
  if size(source_geo_def)[1] < resampler._neighbours
    @warn "Searching for $(resampler._neighbours) neighbours in $(size(source_geo_def)[1]) data points"
  end

  _get_valid_input_index_and_kdtree!(resampler, source_geo_def, target_geo_def; kdtree_class=kdtree_class, nprocs=nprocs)

  if resampler._resample_kdtree === nothing
    return
  end

  resampler._target_lons, resampler._target_lats = get_lonlats(target_geo_def)
  _get_index_array!(resampler)

  # Calculate vertical and horizontal fractional distances t and s
  _get_fractional_distances!(resampler, source_geo_def, target_geo_def)
  _get_target_proj_vectors!(resampler, target_geo_def)
  _get_slices!(resampler, source_geo_def)
end

function _get_valid_input_index_and_kdtree!(resampler::BilinearBase, source_geo_def, target_geo_def; kdtree_class=KDTree, nprocs=1)
  valid_input_index, resample_kdtree = _create_resample_kdtree(resampler, source_geo_def, target_geo_def; kdtree_class=kdtree_class, nprocs=nprocs)

  if resample_kdtree !== nothing
    resampler._valid_input_index = valid_input_index
    resampler._resample_kdtree = resample_kdtree
  else
    # Handle if all input data is reduced away
    _create_empty_bil_info!(resampler, source_geo_def, target_geo_def)
  end
end

function _create_empty_bil_info!(resampler::BilinearBase, source_geo_def, target_geo_def)
  source_size = size(source_geo_def)[1]
  target_size = size(target_geo_def)[1]

  resampler._valid_input_index = trues(source_size)
  resampler._index_array = ones(Int32, target_size, 4)
  resampler.bilinear_s = fill(NaN, target_size)
  resampler.bilinear_t = fill(NaN, target_size)
  resampler.slices_x = zeros(Int32, target_size, 4)
  resampler.slices_y = zeros(Int32, target_size, 4)
  resampler.out_coords_x, resampler.out_coords_y = get_proj_vectors(target_geo_def)
  resampler.mask_slices = resampler._index_array .>= source_size
end

function _get_valid_output_indices(resampler::BilinearBase)
  return vec((resampler._target_lons .>= -180) .& (resampler._target_lons .<= 180) .&
             (resampler._target_lats .<= 90) .& (resampler._target_lats .>= -90))
end

function _get_index_array!(resampler::BilinearBase)
  valid_output_indices = _get_valid_output_indices(resampler)
  index_array = _query_no_distance(resampler._target_lons, resampler._target_lats,
    valid_output_indices, resampler._resample_kdtree,
    resampler._neighbours, resampler._epsilon,
    resampler._radius_of_influence)
  resampler._index_array = _reduce_index_array(resampler, index_array)
end

function _reduce_index_array(resampler::BilinearBase, index_array)
  input_size = sum(resampler._valid_input_index)
  index_mask = index_array .== input_size
  return ifelse.(index_mask, 0, index_array)
end

function _get_fractional_distances!(resampler::BilinearBase, source_geo_def, target_geo_def)
  out_x, out_y = _get_output_xy(target_geo_def, resampler)
  # Get the four closest corner points around each output location
  corner_points, resampler._index_array = _get_four_closest_corners(
    _get_input_xy(source_geo_def, target_geo_def, resampler._valid_input_index, resampler._index_array)..., out_x, out_y,
    resampler._neighbours, resampler._index_array)
  resampler.bilinear_t, resampler.bilinear_s = _get_fractional_distances(
    corner_points, out_x, out_y)
end

function _get_output_xy(target_geo_def, resampler::BilinearBase)
  out_x, out_y = _get_output_xy(target_geo_def)
  valid_output_indices = _get_valid_output_indices(resampler)
  return out_x[valid_output_indices], out_y[valid_output_indices]
end

function _get_input_xy(source_geo_def, target_geo_def, valid_input_index, index_array)
  return _get_input_xy(source_geo_def,
    get_proj(target_geo_def),
    valid_input_index,
    index_array)
end

function _get_target_proj_vectors!(resampler::BilinearBase, target_geo_def)
  try
    resampler.out_coords_x, resampler.out_coords_y = get_proj_vectors(target_geo_def)
  catch
    # AttributeError equivalent - do nothing
  end
end

function _get_slices!(resampler::BilinearBase, source_geo_def)
  shp = size(source_geo_def)
  try
    cols = repeat(1:shp[2], shp[1])
    lines = repeat(1:shp[1], inner=shp[2])
    data = (vec(lines), vec(cols))
  catch BoundsError
    data = (zeros(UInt32, shp[1]), collect(1:shp[1]))
  end

  valid_lines_and_columns = array_slice_for_multiple_arrays(
    resampler._valid_input_index, data)

  resampler.slices_y, resampler.slices_x = array_slice_for_multiple_arrays(
    resampler._index_array, valid_lines_and_columns)
  resampler.mask_slices = resampler._index_array .>= size(source_geo_def)[1]
end

function _create_resample_kdtree(resampler::BilinearBase, source_geo_def, target_geo_def; kdtree_class=KDTree, nprocs=1)
  """Set up kd tree on input."""
  valid_input_index, input_coords = _get_valid_input_index_and_input_coords(resampler, source_geo_def, target_geo_def)
  kdtree = nothing
  if !isempty(input_coords)
    kdtree = KDTree(input_coords')  # Julia KDTree expects features as rows
  end
  return valid_input_index, kdtree
end

function _get_valid_input_index_and_input_coords(resampler::BilinearBase, source_geo_def, target_geo_def)
  valid_input_index, source_lons, source_lats = _get_valid_input_index(
    source_geo_def, target_geo_def,
    resampler._reduce_data, resampler._radius_of_influence)

  input_coords = lonlat2xyz(source_lons, source_lats)
  valid_input_index = vec(valid_input_index)
  input_coords = input_coords[valid_input_index, :]

  return valid_input_index, input_coords
end

function get_sample_from_bil_info(resampler::BilinearBase, data; fill_value=nothing, output_shape=nothing)
  """Resample using pre-computed resampling LUTs."""
  fill_value = _check_fill_value(fill_value, eltype(data))

  res = _resample(_slice_data(resampler, data, fill_value),
    (resampler.bilinear_s, resampler.bilinear_t))

  return _finalize_output_data(resampler, data, res, fill_value)
end

# Helper functions

function _check_fill_value(fill_value, dtype)
  """Check that fill value is usable for the data."""
  if fill_value === nothing
    if dtype <: Integer
      fill_value = 0
    else
      fill_value = NaN
    end
  elseif dtype <: Integer
    if isnan(fill_value)
      fill_value = 0
    elseif isa(fill_value, AbstractFloat)
      fill_value = Int(fill_value)
    end
  end
  return fill_value
end

function _get_output_xy(target_geo_def)
  out_x, out_y = get_proj_coords(target_geo_def)
  return vec(out_x), vec(out_y)
end

function _get_input_xy(source_geo_def, proj, valid_input_index, index_array)
  """Get x/y coordinates for the input area and reduce the data."""
  input_xy_coordinates = mask_coordinates(_get_raveled_lonlats(source_geo_def)...)
  valid_xy_coordinates = array_slice_for_multiple_arrays(valid_input_index, input_xy_coordinates)
  # Expand input coordinates for each output location
  expanded_coordinates = array_slice_for_multiple_arrays(index_array, valid_xy_coordinates)

  return proj_transform(proj, expanded_coordinates...)
end

function array_slice_for_multiple_arrays(idxs, data)
  """Slices multiple arrays using the same indices."""
  return [d[idxs] for d in data]
end

function mask_coordinates(lons, lats)
  """Mask invalid coordinate values."""
  idxs = (find_indices_outside_min_and_max(lons, -180.0, 180.0) .|
          find_indices_outside_min_and_max(lats, -90.0, 90.0))
  return _np_where_for_multiple_arrays(idxs, (NaN, NaN), (lons, lats))
end

function find_indices_outside_min_and_max(data, min_val, max_val)
  """Return array indices outside the given minimum and maximum values."""
  return (data .< min_val) .| (data .> max_val)
end

function _get_raveled_lonlats(geo_def)
  lons, lats = get_lonlats(geo_def)
  if isempty(lons) || isempty(lats)
    throw(ArgumentError("Cannot resample empty data set"))
  elseif length(lons) != length(lats) || size(lons) != size(lats)
    throw(ArgumentError("Mismatch between lons and lats"))
  end
  return vec(lons), vec(lats)
end

function _get_four_closest_corners(in_x, in_y, out_x, out_y, neighbours, index_array)
  """Get bounding corners.

  Get four closest locations from (in_x, in_y) so that they form a
  bounding rectangle around the requested location given by (out_x,
  out_y).
  """
  # Find four closest pixels around the target location
  stride, valid_corners = _get_stride_and_valid_corner_indices(
    out_x, out_y, in_x, in_y, neighbours)
  res = []
  indices = []
  for valid in valid_corners
    x__, y__, idx = _get_corner(stride, valid, in_x, in_y, index_array)
    push!(res, hcat(x__, y__))
    push!(indices, idx)
  end

  return res, hcat(indices...)
end



function _invalid_s_and_t_to_nan(t__, s__)
  return _np_where_for_multiple_arrays(
    (find_indices_outside_min_and_max(t__, 0, 1) .|
     find_indices_outside_min_and_max(s__, 0, 1)),
    (NaN, NaN), (t__, s__))
end


function _ensure_array(val)
  """Ensure that val is an array-like object."""
  if isa(val, Number)
    val = [val]
  end
  return val
end


function _update_fractional_distances(func, t__, s__, corner_points, out_x, out_y)
  idxs = vec(isnan.(t__) .| isnan.(s__))
  if any(idxs)
    new_values = func(corner_points, out_y, out_x)
    updated_t_and_s = _np_where_for_multiple_arrays(idxs, new_values, (t__, s__))
    t__, s__ = _invalid_s_and_t_to_nan(updated_t_and_s...)
  end
  return t__, s__
end

function _get_stride_and_valid_corner_indices(out_x, out_y, in_x, in_y, neighbours)
  out_x_tile = _tile_output_coordinate_vector(out_x, neighbours)
  out_y_tile = _tile_output_coordinate_vector(out_y, neighbours)

  # Get differences in both directions
  x_diff = out_x_tile .- in_x
  y_diff = out_y_tile .- in_y

  return (1:size(x_diff, 1), (
    (x_diff .> 0) .& (y_diff .< 0),  # Upper left corners
    (x_diff .< 0) .& (y_diff .< 0),  # Upper right corners
    (x_diff .> 0) .& (y_diff .> 0),  # Lower left corners
    (x_diff .< 0) .& (y_diff .> 0))  # Lower right corners
  )
end

function _tile_output_coordinate_vector(vector, neighbours)
  return repeat(vector', neighbours, 1)
end

function _get_corner(stride, valid, in_x, in_y, index_array)
  """Get closest set of coordinates from the valid locations."""
  closest_indices = [argmax(valid[i, :]) for i in 1:size(valid, 1)]
  x__ = [in_x[stride[i], closest_indices[i]] for i in eachindex(stride)]
  y__ = [in_y[stride[i], closest_indices[i]] for i in eachindex(stride)]
  idx = [index_array[stride[i], closest_indices[i]] for i in eachindex(stride)]

  # Replace invalid points with NaN
  invalid_mask = .!vec(maximum(valid, dims=2))
  x__[invalid_mask] .= NaN
  y__[invalid_mask] .= NaN

  return x__, y__, idx
end

function _np_where_for_multiple_arrays(idxs, values_for_idxs, otherwise_arrays)
  return [ifelse.(idxs, values_for_idxs[i], arr) for (i, arr) in enumerate(otherwise_arrays)]
end

function _get_valid_input_index(source_geo_def, target_geo_def, reduce_data, radius_of_influence)
  """Find indices of reduced input data."""
  source_lons, source_lats = _get_raveled_lonlats(source_geo_def)

  valid_input_index = .!(find_indices_outside_min_and_max(source_lons, -180.0, 180.0) .|
                         find_indices_outside_min_and_max(source_lats, -90.0, 90.0))

  if reduce_data && is_swath_to_grid_or_grid_to_grid(source_geo_def, target_geo_def)
    valid_input_index .&= get_valid_indices_from_lonlat_boundaries(
      target_geo_def, source_lons, source_lats, radius_of_influence)
  end

  return valid_input_index, source_lons, source_lats
end

function is_swath_to_grid_or_grid_to_grid(source_geo_def, target_geo_def)
  """Check whether the resampling is from swath or grid to grid."""
  # This would need to be implemented based on your geometry types
  # For now, return true as a placeholder
  return true
end

function get_valid_indices_from_lonlat_boundaries(target_geo_def, source_lons, source_lats, radius_of_influence)
  """Get valid indices from lonlat boundaries."""
  # Resampling from swath to grid or from grid to grid
  lonlat_boundary = get_boundary_lonlats(target_geo_def)

  # This would need to be implemented based on your data_reduce module
  # For now, return all true as a placeholder
  return trues(length(source_lons))
end

function get_slicer(data)
  """Get slicing function for 2D or 3D arrays, depending on the data."""
  if ndims(data) == 2
    return _slice2d
  elseif ndims(data) == 3
    return _slice3d
  else
    throw(ArgumentError("Only 2D and 3D arrays are supported"))
  end
end

function _slice2d(values, sl_x, sl_y, mask, fill_value)
  # Slice 2D data
  arr = values[sl_y, sl_x]
  arr[mask] .= fill_value
  return arr[:, 1], arr[:, 2], arr[:, 3], arr[:, 4]
end

function _slice3d(values, sl_x, sl_y, mask, fill_value)
  # Slice 3D data
  arr = values[:, sl_y, sl_x]
  arr[:, mask] .= fill_value
  return arr[:, :, 1], arr[:, :, 2], arr[:, :, 3], arr[:, :, 4]
end

function _resample(corner_points, fractional_distances)
  p_1, p_2, p_3, p_4 = corner_points
  s__, t__ = fractional_distances
  return (p_1 .* (1 .- s__) .* (1 .- t__) .+
          p_2 .* s__ .* (1 .- t__) .+
          p_3 .* (1 .- s__) .* t__ .+
          p_4 .* s__ .* t__)
end


function _query_no_distance(target_lons, target_lats, valid_output_index, kdtree, neighbours, epsilon, radius)
  """Query the kdtree. No distances are returned."""
  target_lons_valid = vec(target_lons)[valid_output_index]
  target_lats_valid = vec(target_lats)[valid_output_index]

  target_coords = lonlat2xyz(target_lons_valid, target_lats_valid)

  # Use NearestNeighbors.jl knn function
  indices, distances = knn(kdtree, target_coords', neighbours)

  # Convert to matrix format similar to Python version
  index_array = hcat(indices...)' # Transpose to match Python kdtree output format

  # Filter by distance threshold
  distance_matrix = hcat(distances...)'
  index_array[distance_matrix.>radius] .= size(kdtree.data, 2) + 1

  return Int32.(index_array)
end

