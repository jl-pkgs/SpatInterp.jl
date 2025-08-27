# using NearestNeighbors
# using LinearAlgebra

"""
计算二次方程系数
"""
function calc_abc(p1::P, p2::P, p3::P, p4::P, out_x::T, out_y::T) where {T<:Real,P<:Point{T}}
  # 参考点之间的经向分离
  x_21 = p2.x - p1.x
  x_31 = p3.x - p1.x
  x_42 = p4.x - p2.x

  y_21 = p2.y - p1.y
  y_31 = p3.y - p1.y
  y_42 = p4.y - p2.y

  a = x_31 * y_42 - y_31 * x_42
  b = out_y * (x_42 - x_31) - out_x * (y_42 - y_31) +
      x_31 * p2.y - y_31 * p2.x +
      y_42 * p1.x - x_42 * p1.y
  c = out_y * x_21 - out_x * y_21 + p1.x * p2.y - p2.x * p1.y
  return a, b, c
end


function solve_quadratic(a::T, b::T, c::T; min::T=T(0), max::T=T(1)) where {T<:Real}
  # a x^2 + b x + c = 0
  Δ = b * b - 4 * a * c

  sqrt_Δ = sqrt(Base.max(Δ, T(0)))
  # 求解二次多项式
  x_1 = (-b + sqrt_Δ) / (2 * a)
  x_2 = (-b - sqrt_Δ) / (2 * a)
  x_3 = -c / b # (线性情况, a = 0)

  # 找到有效解，即 0 <= t <= 1
  # @show x_1, x_2, x_3
  t = not_valid_root(x_1) ? x_2 : x_1
  not_valid_root(t) && (t = x_3)
  not_valid_root(t) && (t = T(NaN))
  return t
end


function not_valid_root(x::T; min::T=T(0), max::T=T(1)) where {T}
  isnan(x) && return true
  (x < min || x > max) && return true
end
