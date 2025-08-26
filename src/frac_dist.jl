function solve_another_frac_dist(t::T, y_1::T, y_2::T, y_3::T, y_4::T, out_y::T) where T<:Real
  y_21 = y_2 - y_1
  y_43 = y_4 - y_3

  denominator = y_43 * t - y_21 * t + y_3 - y_1
  numerator = out_y - y_1 - y_21 * t
  s = numerator / denominator
  (abs(denominator) < eps(T) || not_valid_root(s)) && (s = T(NaN))
  return s
end

"""
获取不规则情况的分数距离
"""
function frac_dist_irregular(p1::P, p2::P, p3::P, p4::P, out_x::T, out_y::T) where {T<:Real,N,P<:Point{N,T}}
  a, b, c = calc_abc(p1, p2, p3, p4, out_y, out_x)
  t = solve_quadratic(a, b, c)
  s = solve_another_frac_dist(t, p1.y, p3.y, p2.y, p4.y, out_y) # 注意p2和p3交换
  return t, s
end


# vertical parallel
function frac_dist_uprights_parallel(p1::P, p2::P, p3::P, p4::P, out_x::T, out_y::T) where {T<:Real,N,P<:Point{N,T}}
  a, b, c = calc_abc(p1, p3, p2, p4, out_y, out_x)
  s = solve_quadratic(a, b, c) # note diff from irregular
  t = solve_another_frac_dist(s, p1.y, p2.y, p3.y, p4.y, out_y) # note diff
  return t, s
end

# vertical and horizontal are parallel
function frac_dist_parallellogram(p1::P, p2::P, p3::P, out_x::T, out_y::T) where {T<:Real,N,P<:Point{N,T}}
  x_21 = p2.x - p1.x
  x_31 = p3.x - p1.x
  y_21 = p2.y - p1.y
  y_31 = p3.y - p1.y

  t = (x_21 * (out_y - p1.y) - y_21 * (out_x - p1.x)) / (x_21 * y_31 - y_21 * x_31)
  not_valid_root(t) && (t = T(NaN))

  s = (out_x - p1.x + x_31 * t) / x_21
  not_valid_root(s) && (s = T(NaN))
  return t, s
end


function frac_dist(p1::P, p2::P, p3::P, p4::P, out_x::T, out_y::T) where {T<:Real,N,P<:Point{N,T}}
  t, s = frac_dist_irregular(p1, p2, p3, p4, out_x, out_y)

  # vertical parallel
  if isnan(t + s)
    t, s = frac_dist_uprights_parallel(p1, p2, p3, p4, out_x, out_y)
  end
  # vertical and horizonal
  if isnan(t + s)
    t, s = frac_dist_parallellogram(p1, p2, p3, out_x, out_y)
  end
  return t, s
end


export frac_dist_irregular, frac_dist_uprights_parallel, frac_dist_parallellogram
export frac_dist
