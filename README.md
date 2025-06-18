# BarnesDoubleGamma

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://Paul-Roux-Phys.github.io/BootstrapVirasoro.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://Paul-Roux-Phys.github.io/BootstrapVirasoro.jl/dev/)
[![Build Status](https://github.com/Paul-Roux-Phys/BootstrapVirasoro.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Paul-Roux-Phys/BootstrapVirasoro.jl/actions/workflows/CI.yml?query=branch%3Amain)

This is a Julia package for performing numerical bootstrap computation in 2-dimensional conformal field theories with Virasoro symmetry, written by Paul Roux. Part of the code is based on a previous [python code](https://gitlab.com/s.g.ribault/Bootstrap_Virasoro/).

For more information, please refer to [the documentation](https://paul-roux-phys.github.io/BootstrapVirasoro.jl/stable/).

BibTex: 
@misc{roux251,
	author = "Roux, Paul",
	title = {BarnesDoubleGamma.jl},
	type = {code},
	url = {https://juliapackages.com/p/barnesdoublegamma},
	year = {2025},
}

## Installation

To install this package, run

```julia
using Pkg; Pkg.add("BootstrapVirasoro")
# or
# using Pkg; Pkg.add("https://github.com/Paul-Roux-Phys/BootstrapVirasoro.jl.git")
```

in a Julia script, or

```julia-repl
julia> ]add BarnesDoubleGamma
julia> ]add https://github.com/Paul-Roux-Phys/BootstrapVirasoro.jl.git
```

in a Julia REPL.
