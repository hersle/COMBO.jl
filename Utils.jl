using Printf
using Dates

# Kronecker δ
δ(i::Integer, j::Integer) = i == j ? 1 : 0

function quadroots(a::Real, b::Real, c::Real)
    d = b^2 - 4*a*c
    if d >= 0
        x1 = (-b + √(d)) / (2*a)
        x2 = (-b - √(d)) / (2*a)
    else
        x1 = NaN
        x2 = NaN
    end
    return min(x1, x2), max(x1, x2) # order roots in ascending order
end

# Add nextra points between each point in x
function extend(x::AbstractVector, nextra::Integer)
    dx = diff(x)
    return sort(vcat(x, (x[1:end-1] .+ i/(nextra+1)*dx for i in 1:nextra)...))
end

integration_points(sol::ODESolution; nextra=0) = extend(sol.t, nextra)

spline(x, y) = scale(interpolate(y, BSpline(Cubic(Line(OnGrid())))), x) # TODO: OOB BC

#=
function multirange(posts, lengths)
    return vcat((range(posts[i], posts[i+1], length=lengths[i])[1:end-1] for i in range(1, length(lengths)))..., posts[end])
end
=#

function format_time_variations(bg::Background, x::Real)
    return @sprintf("x = %+4.2f, a = %6.4f, z = %7.2f, η = %4.1f Gyr, t = %8.5f Gyr",
                     x,          a(x),      z(x),      η(bg,x) / Gyr, t(bg,x) / Gyr)
end

function fixed_point_iterate(xfunc::Function, x0; tol::Float64=1e-10, maxiters::Integer=100)
    converged = false
    iters = 0
    x = x0
    while !converged && iters < maxiters
        x_new = xfunc(x)
        converged = maximum(abs.(x_new .- x)) < tol # TODO: or vector 2-norm?
        x = x_new
        iters += 1
    end
    @assert iters < maxiters
    return x
end
