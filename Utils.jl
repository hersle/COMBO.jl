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

function _spline_integral_generic(f::Function, x1::Float64, x2::Float64, y1; solver=Tsit5(), name="unnamed quantity", abstol=1e-8, reltol=1e-8, maxiters=1e7, xskip=1, benchmark=false, verbose=false, kwargs...)
    if benchmark
        sol = solve(ODEProblem(f, y1, (x1, x2)), solver; maxiters=maxiters, kwargs..., abstol=abstol, reltol=reltol) # pre-compile before measuring
    end
    t1 = now()
    sol = solve(ODEProblem(f, y1, (x1, x2)), solver; maxiters=maxiters, kwargs..., abstol=abstol, reltol=reltol)
    t2 = now()
    dt = t2 - t1

    success = sol.retcode == SciMLBase.ReturnCode.Success
    @assert success "failed integrating $name"

    # print some statistics
    if verbose
        print("Integrated $name ")
        print("on [$(sol.t[1]), $(sol.t[end])] ")
        print("with $(length(y1)) variables and $(length(sol.t)) points ")
        print("using $(typeof(sol.alg)), abstol $abstol and reltol $reltol ")
        print("in $dt\n")
    end

    # spline wants points with ascending x values,
    # while the integrator can output them in a different order
    sortinds = sortperm(sol.t)
    x = sol.t[sortinds]
    y = sol.u[sortinds]

    # let caller skip points before splining (to save memory)
    filter = 1:xskip:length(x)
    if filter[end] != length(x)
        filter = vcat(filter, length(x)) # include endpoint regardless of "xskip divisibility"
    end
    x = x[filter]
    y = y[filter]

    return x, y
end

# integrate systems of equations with in-place RHS
function _spline_integral(dy_dx!::Function, x1::Float64, x2::Float64, y1::Vector{Float64}; kwargs...)
    f!(dy, y, p, x) = dy_dx!(x, y, dy)
    x, y = _spline_integral_generic(f!, x1, x2, y1; kwargs...)
    y = [y[ix][iy] for iy in 1:length(y1), ix in 1:length(x)] # convert vector of vectors to 2D matrix
    return x, [spline(x, y[i,:]) for i in 1:length(y1)]
end

# integrate scalar equations with out-of-place (scalar) RHS
function _spline_integral(dy_dx::Function, x1::Float64, x2::Float64, y1::Float64; kwargs...)
    f(y, p, x) = dy_dx(x, y)
    x, y = _spline_integral_generic(f, x1, x2, y1; kwargs...)
    return x, spline(x, y)
end

spline(x, y) = Spline1D(x, y; bc="error")
#spline(x, y) = scale(interpolate(y, BSpline(Cubic(Line(OnGrid())))), x) # TODO: OOB BC

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
