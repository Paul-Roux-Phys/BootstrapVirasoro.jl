# `BootstrapVirasoro.jl` documentation

```@meta
CurrentModule = BootstrapVirasoro
DocTestSetup = quote
    using BootstrapVirasoro
	setprecision(BigFloat, 256)
end
```

This is the documentation page for the `BootstrapVirasoro` package. `BootstrapVirasoro` is a package for doing bootstrap computations in theories with Virasoro symmetry.

## Contents

* [Installation](installation.md)

For any new user, we advise taking a look at the two following files, which contain simple self-contained code examples:

* [compute blocks](conformalblocks.md): learn how to compute various types of virasoro conformal blocks,
* [solve bootstrap equations](bootstrapeqs.md): learn how to use the package's utilities to solve bootstrap equations and pretty-print their solutions.

For more detailed documentation of the public API of the package, refer to 

* [public API](public_api.md)

