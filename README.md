# BootstrapVirasoro.jl

[![Build Status](https://github.com/Paul-Roux-Phys/BootstrapVirasoro.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Paul-Roux-Phys/BootstrapVirasoro.jl/actions/workflows/CI.yml?query=branch%3Amain)

This is a Julia package for performing numerical bootstrap computation in 2-dimensional conformal field theories with Virasoro symmetry, written by Paul Roux. Part of the code is based on a previous [python code](https://gitlab.com/s.g.ribault/Bootstrap_Virasoro/).

The best way to learn how to use the package is to have a look at the [examples](./examples)

If you use this program in your scientific work, you should cite it as

```text
BibTex: 
@misc{BVJL,
	author = "Roux, Paul",
	title = {BootstrapVirasoro.jl},
	type = {code},
	url = {https://github.com/Paul-Roux-Phys/BootstrapVirasoro.jl},
	year = {2025},
}
```

## Installation

To install this package, run

```julia
using Pkg; Pkg.add("BootstrapVirasoro")
```

or
```julia
using Pkg; Pkg.add("https://github.com/Paul-Roux-Phys/BootstrapVirasoro.jl.git")
```

in a Julia script, or

```julia-repl
julia> ]add BootstrapVirasoro
# or
julia> ]add https://github.com/Paul-Roux-Phys/BootstrapVirasoro.jl.git
```

in a Julia REPL.
Installing directly from git will ensure you have the most up-to-date version, while the version on the Julia registry is kept more stable.
