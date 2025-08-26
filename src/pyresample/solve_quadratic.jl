function _solve_quadratic(a__, b__, c__; min_val=0.0, max_val=1.0)
  """Solve quadratic equation.

  Solve quadratic equation and return the valid roots from interval
  [min_val, max_val].
  """
  a__ = _ensure_array(a__)
  b__ = _ensure_array(b__)
  c__ = _ensure_array(c__)
  discriminant = b__ .* b__ .- 4 .* a__ .* c__

  # Solve the quadratic polynomial
  x_1 = (-b__ .+ sqrt.(discriminant)) ./ (2 .* a__)
  x_2 = (-b__ .- sqrt.(discriminant)) ./ (2 .* a__)
  # Linear case
  x_3 = -c__ ./ b__

  # Find valid solutions, ie. 0 <= t <= 1
  x__ = ifelse.(find_indices_outside_min_and_max(x_1, min_val, max_val) .| isnan.(x_1),
    x_2, x_1)

  x__ = ifelse.(find_indices_outside_min_and_max(x__, min_val, max_val) .| isnan.(x__),
    x_3, x__)

  x__ = ifelse.(find_indices_outside_min_and_max(x__, min_val, max_val),
    NaN, x__)

  return x__
end


function _calc_abc(corner_points, out_y, out_x)
  """Calculate coefficients for quadratic equation."""
  pt_1, pt_2, pt_3, pt_4 = corner_points
  # Pairwise longitudinal separations between reference points
  x_21 = pt_2[:, 1] .- pt_1[:, 1]
  x_31 = pt_3[:, 1] .- pt_1[:, 1]
  x_42 = pt_4[:, 1] .- pt_2[:, 1]

  # Pairwise latitudinal separations between reference points
  y_21 = pt_2[:, 2] .- pt_1[:, 2]
  y_31 = pt_3[:, 2] .- pt_1[:, 2]
  y_42 = pt_4[:, 2] .- pt_2[:, 2]

  a__ = x_31 .* y_42 .- y_31 .* x_42
  b__ = out_y .* (x_42 .- x_31) .- out_x .* (y_42 .- y_31) .+
        x_31 .* pt_2[:, 2] .- y_31 .* pt_2[:, 1] .+
        y_42 .* pt_1[:, 1] .- x_42 .* pt_1[:, 2]
  c__ = out_y .* x_21 .- out_x .* y_21 .+
        pt_1[:, 1] .* pt_2[:, 2] .- pt_2[:, 1] .* pt_1[:, 2]

  return a__, b__, c__
end
