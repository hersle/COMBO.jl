using Printf
using Dates

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

# internal function that computes y(x) on a cubic spline,
# given dy/dx, from x=x1 with y(x1)=y1, to x=x2
# TODO: generalize to take (x0, y0), integrate in both directions, and combine splines?
# TODO: write about choice of solver
# CLASS uses "rkck4" (non-stiff) or "ndf15" (stiff) evolver (https://github.com/lesgourg/class_public/blob/aa92943e4ab86b56970953589b4897adf2bd0f99/include/common.h)
# Hans says "RK2" is sufficient (https://cmb.wintherscoming.no/theory_numerical.php)
# DifferentialEquations recommends (for high accuracy (low tolerance)) auto-switching stiffness solver AutoVern7(Rodas5()) solver (https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/#Unknown-Stiffness-Problems)
# TODO: autodiff?
#

spline(x, y) = Spline1D(x, y; bc="error")
#spline(x, y) = scale(interpolate(y, BSpline(Cubic(Line(OnGrid())))), x) # TODO: OOB BC

splinex(spline::ODESolution) = spline.t # TODO: does this make sense to use?

splinex(spline::Spline1D) = Dierckx.get_knots(spline) # TODO: does this make sense to use?
spliney(spline::Spline1D) = spline(splinex(spline))

splinex(spline::Spline2D) = Dierckx.get_knots(spline)[1]
spliney(spline::Spline2D) = Dierckx.get_knots(spline)[2]
splinez(spline::Spline2D) = spline(splinex(spline), spliney(spline))

function splinejoin(x1, x2, spl1, spl2)
    x = unique(vcat(x1, x2)) # don't duplicate midpoint
    y = vcat(spl1.(x1), spl2.(x2)[2:end])
    return x, spline(x, y)
end

function splinejoin(x1, x2, spl1s::Vector{Spline1D}, spl2s::Vector{Spline1D})
    x, spls = [], []
    for (spl1, spl2) in zip(spl1s, spl2s)
        x, spl = splinejoin(x1, x2, spl1, spl2)
        push!(spls, spl)
    end
    return x, spls
end

# use spline points for plotting,
# but add nextra points between each of them
# TODO: make accessible to plotter?
# take an array of x values (e.g. spline points),
# then add nextra points between each of them
function extendx(x::Vector{Float64}, nextra::Integer)
    dx = diff(x)
    return sort(vcat(x, (x[1:end-1] .+ i/(nextra+1)*dx for i in 1:nextra)...))
end

function multirange(posts, lengths)
    return vcat((range(posts[i], posts[i+1], length=lengths[i])[1:end-1] for i in range(1, length(lengths)))..., posts[end])
end

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
