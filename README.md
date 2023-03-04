AST5220 project: Einstein-Boltzmann solver
==========================================

Try to be a hipster and use Julia...

Instructions
------------

[Download and install Julia](https://julialang.org/downloads/), add `julia` to your `$PATH`, then
```
git clone https://github.com/hersle/AST5220-project
cd AST5220-project/
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.resolve(); Pkg.precompile()' # install, resolve and precompile dependencies
julia --project=. Milestone1.jl # produce output for first milestone
```

Julia is a precompiled language,
so dependency installation, precompilation and the first run can take long!

For the plots to look as intended, PGFPlots must be installed.
Otherwise, a different plotting backend is used,
and the plots may not look as nice.
