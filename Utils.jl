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

# internal function that computes y(x) on a cubic spline,
# given dy/dx, from x=x1 with y(x1)=y1, to x=x2
# TODO: generalize to take (x0, y0), integrate in both directions, and combine splines?
function _spline_integral(dy_dx::Function, x1::Float64, x2::Float64, y1::Float64, abstol::Float64)
    sol = solve(ODEProblem((y, p, x) -> dy_dx(x, y), y1, (x1, x2)), Tsit5(); reltol=1e-15, abstol=abstol)
    return sol # use dense output to interpolate
end

function multirange(posts, lengths)
    return vcat((range(posts[i], posts[i+1], length=lengths[i])[1:end-1] for i in range(1, length(lengths)))..., posts[end])
end
