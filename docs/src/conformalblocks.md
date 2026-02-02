# Computing conformal blocks

The library uses the [Zamolodchikov recursion](https://en.wikipedia.org/wiki/Virasoro_conformal_block#Zamolodchikov's_recursive_representation)
to compute four-point conformal blocks on the sphere, and one-point conformal blocks on the torus.

Since the library is optimized for bootstrap computations in arbitrary precision, it caches any intermediary computations that can be reused.
To this end, we use several types, each storing its own cache: `CentralCharge`, `ConformalDimension`, `Field`, `Correlation`, `Block`.

```julia
using BootstrapVirasoro

# set the precision of high-precision floats
setprecision(BigFloat, 30, base=10)

# Define a high precision number (BigFloat type) from float literals:
β = big"0.3" + big"0.6"*im

# Define a corresponding central charge
c = CC(β = β)

# Define ConformalDimension objects.
D1 = CD(c, P = big"1.2")
D2 = CD(c, r = 2, s = 1//2) # a // b is the Julia syntax for rational numbers.
# we can use regular floats as well, but beware that the will not be correctly rounded
# when they are converted to a BigFloat type.
D3 = CD(c, P=1.2) # D3 not strictly equal to D2.

# Define Field objects.
V1 = Field(c, P = big"1" + big"0.8"*im)           # diagonal field of momentum P
VΔ = Field(c, Δ = big"0.4" + big"0.3"*im)         # diagonal field of conformal dimension Δ
V2 = Field(c, r = 2 , s= 1//2)                    # non-diagonal field defined from Kac indices
V2diag = Field(c, r = 2 , s= 1//2, diagonal=true) # we can force the field to be diagonal
V3 = Field(c, r = 2, s = 1, diagonal = true)      # degenerate field
V4 = Field(D1, D2)                                # define a field from a pair of conformal dimensions
V2.dims                                           # inspect left and right conformal dimensions
V2.dims.left
swap_lr(V2)                                       # swap the left and right conformal dimensions (space parity)
swap_lr(V4)

#=
define a correlation with Correlation(fields, Nmax)
This precomputes the coefficients R_mn that are the residues of conformal blocks
for this correlation, up to m * n = Nmax,
and also the coefficients C^N_mn obtained by solving the
Zamolodchikov recursion relation. =#

# non-chiral, four-point
Cor = Correlation(V1, V1, V1, V1, 20)
# chiral, four-point
bndrycor = Correlation(D1, D1, D1, D1, 20)
# non-chiral, one-point (torus)
toruscor = Correlation(V1, 20)

bndrycor.Rmn.s[2, 1] # the correlation stores the value of Rmn, which can be inspected

# define a non-chiral s-channel block. This stores the series coefficients
b4 = NCBlock(Cor, :s, V4)
b2 = NCBlock(Cor, :s, V2)
# evaluate the block at a value of the cross-ratio z
b4(0.95 + 2.04im)

# define a chiral s-channel block
bndry = CBlock(bndrycor, :s, D2)
x = big"0.84"+big"0.3"*im
# evaluate it at a position
@time bndry(x)

# we can linearly combine block and easily evaluate the linear combination
# the first argument is a list of blocks, the second argument a list of coefficients
# here we compute 2 b2 + 3 b4.
Lc = LinearCombinationBlock([b2, b4], [2, 3])
# evaluate
Lc(big"1"+big"1"*im)

# when evaluating many blocks at the same position, we can compute once and for all
# all the data about the positions that is common to all blocks, i.e. the nome q, powers
# of q, and the prefactors x^E1 (1-x)^E3 theta_3(q)^E_4. The evaluation of the blocks is
# several orders of magnitude faster if we input the cache instead of the raw position.
# the cache is created with BootstrapVirasoro.PosCache.
cache = BootstrapVirasoro.PosCache(x, bndrycor.fields, :s, 20)
@time bndry(cache)
```
