# Computing conformal blocks

The library uses the [Zamolodchikov recursion](https://en.wikipedia.org/wiki/Virasoro_conformal_block#Zamolodchikov's_recursive_representation)
to compute four-point conformal blocks on the sphere, and one-point conformal blocks on the torus.
Because it is common in bootstrap computations to have to evaluate a given block at hundreds of positions,
the library first creates a `Block` object which stores the coefficients of the block's series expansion, which can then be evaluated at any number of positions.

## Examples

The two blocks of code below can be copy-pasted and ran as is. More detailed documentation of the package's functionality can be found [below](@ref Documentation).

### Chiral conformal blocks

```julia
using BootstrapVirasoro
setprecision(BigFloat, 20, base=10) # work with 20 digits of precision
# define a central charge. big"..." is the syntax for abitrary precision numbers.
# b is the Liouville parameter
c = CC(b = big"0.523"+big"0.1"*im)
# if only standard (16-digit) precision is needed, just do c = CC(b=0.523+0.1*im)
# standard precision is *much* faster.
d1 = ConformalDimension(c, r=1, s=0) # Conformal dimension with Kac indices (r, s) = (1, 0).
d2 = CD(c, P=big"0.54") # CD is synonymous to ConformalDimension, P is the momentum.
d3 = CD(c, Δ=big"4.32") # Δ is the conformal dimension. Type it as \Delta, then tab or enter.
d4 = d3
d = CD(c, P=big"0.32"+big"0.1"*im)

# four-point block:
co1 = Correlation(d1, d2, d3, d4, 10) # Create a correlation object that stores the R_{m, n} up to Nmax=10.
b1 = ChiralBlock(co1, :s, d) # Create a chiral s-channel four-point block with dimension d in the channel.
b2 = CBlock(co1, :t, d) # t-channel block. CBlock is synonymous of ChiralBlock.
x = big"0.3" + big"0.2" * im
println(b1(x)) # print the value of the block at position `x`.

# one-point block:
co2 = Correlation(d1, 10)
b3 = ChiralBlock(co2, d) # chiral one-point block
τ = big"1.1" + big"0.3" * im
println(b3(τ))
```

### Non-chiral conformal blocks

```julia
using BootstrapVirasoro
setprecision(BigFloat, 20, base=10) # work with 20 digits of precision
# define a central charge. big"..." is the syntax for abitrary precision numbers.
# b is the Liouville parameter
c = CC(b = big"0.523"+big"0.1"*im)
V2 = Field(c, P=big"0.54") # Field with momentum P=0.54.
V1 = Field(c, r=1, s=1) # non-diagonal field of dimensions (Δ_(1, 1), Δ_(1, -1))
V1 = Field(c, r=1, s=2, diagonal=true) # degenerate field with Kac indices (r, s) = (1, 2).
V3 = Field(c, Δ=big"4.32") # Δ is the conformal dimension. Type it as \Delta, then tab or enter.
V4 = V3
V = Field(c, P=big"0.32"+big"0.1"*im)
V_degenerate = Field(c, r=2, s=3, diagonal=true)

# four-point block:
co1 = Correlation(V1, V2, V3, V4, 10) # Create a correlation object that stores the R_{m, n} up to Nmax=10.
b1 = NonChiralBlock(co1, :s, V) # Create a non-chiral s-channel four-point block with field V in the channel.
b2 = NCBlock(co1, :t, V_degenerate) # regularised t-channel four-point block with degenerate field V_degenerate in the channel.
                                    # NCBlock is synonymous of NonChiralBlock
x = big"0.3" + big"0.2" * im
println(b1(x)) # print the value of the block at position `x`.

# one-point block:
co2 = Correlation(d1, 10)
b2 = NCBlock(co2, d) # non chiral one-point block
τ = big"1.1" + big"0.3" * im
println(b2(τ))
```

## Documentation

The library provides a hierarchy of conformal block types:

```@docs
Block
ChiralBlock
NonChiralBlock
LinearCombinationBlock
```

Defining and computing blocks requires defining a central charge, and conformal dimensions of the external and channel fields.

```@docs
CentralCharge
ConformalDimension
Field
```

The parameters of the external fields are stored in a `Correlation` object, that also stores the
[residues ``R_{m, n}``](https://en.wikipedia.org/wiki/Virasoro_conformal_block#Zamolodchikov's_recursive_representation) for each channel,
and their sums 

$$C^N_{m,n} = R_{m,n}\left(\delta_{N-mn,0} +
\sum_{m'n'\leq N-mn}
\frac{C^{N-mn}_{m',n'}}{\delta_{(m,-n)}-\delta_{(m',n')}} \right)$$

used in the Zamolodchikov recursion.

They are all created via a common interface 

```@docs
Correlation
```
