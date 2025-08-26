function _get_fractional_distances(corner_points, out_x, out_y)
  """Calculate vertical and horizontal fractional distances t and s."""
  # General case, ie. where the corners form an irregular rectangle
  t__, s__ = _invalid_s_and_t_to_nan(
    _get_fractional_distances_irregular(corner_points, out_y, out_x)...)

  # Update where verticals are parallel
  t__, s__ = _update_fractional_distances(
    _get_fractional_distances_uprights_parallel,
    t__, s__, corner_points, out_x, out_y)

  # Update where both verticals and horizontals are parallel
  corner_points_reduced = (corner_points[1], corner_points[2], corner_points[3])
  t__, s__ = _update_fractional_distances(
    _get_fractional_distances_parallellogram,
    t__, s__, corner_points_reduced, out_x, out_y)

  return t__, s__
end


function _solve_another_fractional_distance(f__, y_corners, out_y)
  """Solve parameter t__ from s__, or vice versa."""
  y_1, y_2, y_3, y_4 = y_corners

  y_21 = y_2 .- y_1
  y_43 = y_4 .- y_3

  g__ = ((out_y .- y_1 .- y_21 .* f__) ./
         (y_3 .+ y_43 .* f__ .- y_1 .- y_21 .* f__))

  # Limit values to interval [0, 1]
  g__ = ifelse.(find_indices_outside_min_and_max(g__, 0, 1), NaN, g__)

  return g__
end



function _get_fractional_distances_uprights_parallel(corner_points, out_y, out_x)
  """Get parameters for the case where uprights are parallel."""
  pt_1, pt_2, pt_3, pt_4 = corner_points
  # Get the valid roots from interval [0, 1]. Note the different order needed here.
  corner_points_reordered = (pt_1, pt_3, pt_2, pt_4)
  s__ = _solve_quadratic(_calc_abc(corner_points_reordered, out_y, out_x)...,
    min_val=0.0, max_val=1.0)

  # Calculate parameter t
  y_corners = (pt_1[:, 2], pt_2[:, 2], pt_3[:, 2], pt_4[:, 2])
  t__ = _solve_another_fractional_distance(s__, y_corners, out_y)

  return t__, s__
end

function _get_fractional_distances_parallellogram(points, out_y, out_x)
  """Get parameters for the case where all sides are parallel (parallelogram)."""
  pt_1, pt_2, pt_3 = points
  # Pairwise longitudinal separations between reference points
  x_21 = pt_2[:, 1] .- pt_1[:, 1]
  x_31 = pt_3[:, 1] .- pt_1[:, 1]

  # Pairwise latitudinal separations between reference points
  y_21 = pt_2[:, 2] .- pt_1[:, 2]
  y_31 = pt_3[:, 2] .- pt_1[:, 2]

  t__ = (x_21 .* (out_y .- pt_1[:, 2]) .- y_21 .* (out_x .- pt_1[:, 1])) ./
        (x_21 .* y_31 .- y_21 .* x_31)
  t__ = ifelse.(find_indices_outside_min_and_max(t__, 0.0, 1.0), NaN, t__)

  s__ = (out_x .- pt_1[:, 1] .+ x_31 .* t__) ./ x_21
  s__ = ifelse.(find_indices_outside_min_and_max(s__, 0.0, 1.0), NaN, s__)

  return t__, s__
end


function _get_fractional_distances_irregular(corner_points, out_y, out_x)
  """Get parameters for the case where none of the sides are parallel."""
  # Get the valid roots from interval [0, 1]
  t__ = _solve_quadratic(_calc_abc(corner_points, out_y, out_x)...,
    min_val=0.0, max_val=1.0)

  # Calculate parameter s
  pt_1, pt_2, pt_3, pt_4 = corner_points
  y_corners = (pt_1[:, 2], pt_3[:, 2], pt_2[:, 2], pt_4[:, 2])
  s__ = _solve_another_fractional_distance(t__, y_corners, out_y)

  return t__, s__
end
