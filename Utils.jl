using Printf

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
# TODO: write about choice of solver
# CLASS uses "rkck4" (non-stiff) or "ndf15" (stiff) evolver (https://github.com/lesgourg/class_public/blob/aa92943e4ab86b56970953589b4897adf2bd0f99/include/common.h)
# Hans says "RK2" is sufficient (https://cmb.wintherscoming.no/theory_numerical.php)
# DifferentialEquations recommends (for high accuracy (low tolerance)) auto-switching stiffness solver AutoVern7(Rodas5()) solver (https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/#Unknown-Stiffness-Problems)
# TODO: autodiff?
function _spline_integral(dy_dx::Function, x1::Float64, x2::Float64, y1; solver=AutoVern7(Rodas5(autodiff=false, diff_type = Val{:central})), reltol::Float64=1e-10, abstol::Float64=1e-10)
    sol = solve(ODEProblem((y, p, x) -> dy_dx(x, y), y1, (x1, x2)), solver, reltol=reltol, abstol=abstol)
    return sol # use dense output to interpolate
end

function multirange(posts, lengths)
    return vcat((range(posts[i], posts[i+1], length=lengths[i])[1:end-1] for i in range(1, length(lengths)))..., posts[end])
end

function format_time_variations(co::ΛCDM, x::Real)
    return @sprintf("x = %+4.2f, a = %6.4f, z = %7.2f, η = %4.1f Gyr, t = %8.5f Gyr",
                     x,          a(x),      z(x),      η(co,x) / Gyr, t(co,x) / Gyr)
end
