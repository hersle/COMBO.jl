AST5220 project: Einstein-Boltzmann solver
==========================================

Try to be a hipster and use Julia...

Instructions
------------

[Download and install Julia](https://julialang.org/downloads/), add `julia` to your `$PATH`, then
```
git clone https://github.com/hersle/AST5220-project
cd AST5220-project/
julia --project=. -e 'import Pkg; Pkg.instantiate()' # install dependencies
julia --project=. Milestone1.jl # produce output for first milestone
```

Since Julia is a precompiled language,
so it can take a long time to install all dependencies and run the program for the first time!
Run `Milestone1.jl` a second time to better assess the actual runtime.
