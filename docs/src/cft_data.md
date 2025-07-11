# Basic types

The package has a few types to represent objects in a Virasoro CFT:

* `CentralCharge`s
* `ConformalDimension`s
* `Field`s
* `Correlation`s
* `Block`s

## Central charges

We parametrise the central charge of our theories in terms of
variables `c`, `B`, `b` or `β` related by

$$c = 13 + 6B + 6 B^{-1} \quad , \quad B = b^2 = -\beta^2,
\quad B = \frac{c-13 \pm \sqrt{(c-1)(c-25)}}{12}$$

By convention we keep $\beta = ib$. In $O(n)$ and $U(n)$ models
the central charge is related to $n$ via

$$n = - 2 \cos(\pi \beta^2)$$

```@docs
CentralCharge
```

## Conformal Dimensions

We parametrise the conformal dimensions
in terms of variables $\Delta, P, p, \delta$, related by

$$\Delta = \frac{c-1}{24} + \delta  \quad , \quad \delta = P^2 = -p^2$$

The variable $P$ is called the momentum, and $\Delta$ is the eigenvalue
of the Virasoro generator $L_0$. By convention, we always keep
$P=ip$. Moreover, we introduce the following parametrisation of
dimensions in terms of Kac indices $r, s$:

$$P_{(r,s)}=\frac{1}{2}(\beta r - \beta^{-1}s)$$

where $r,s$ are arbitrary numbers. This convention is different from the
one in [Sylvain\'s
code](https://gitlab.com/s.g.ribault/Bootstrap_Virasoro.git), but
similar to our more recent conventions, such as in [Sylvain's review on
solvable CFTs](https://github.com/ribault/CFT-Review).

```@docs
ConformalDimension
P_rs
```

`ConformalDimension`s can be shifted by units of $\beta^{-1}/2$:

```@docs
shift(::ConformalDimension, i)
```

## Fields

We see a field as the data of left and right conformal dimensions, plus labels
indicating whether the field is diagonal and degenerate.
We denote $V_{(r, s)}$ a field with left and right dimensions
$(\Delta_{(r, s)}, \Delta_{(r, -s)})$.

```@docs
Field
isdiagonal
isdegenerate
shift(V::Field, i)
```

## Correlation functions

The residues of conformal blocks only depend on the external insertion point.
We store this data in a `Correlation` object, which can represent any concrete
(one-point or four-point) correlation.

```@docs
Correlation
```

## Conformal blocks

We compute conformal blocks via the
[Zamolodchikov Spectra](https://en.wikipedia.org/wiki/Virasoro_conformal_block).
Blocks can be created via a unified interface, which does a best-effort dispatch
to the right type of conformal block depending on its inputs. A block object
contains the values of the coefficients of the block as a series in the nome `q`.

```@docs
Block
```

Blocks can be evaluated at a position. For four-point blocks, this is the cross-ratio `x`
for one-point blocks this is the modular parameter $\tau$.

```jl
c = CentralCharge(c=0.54)
V1 = Field(c, r=2, s=3//2)
co = Correlation(V1, Δmax=15)
b = Block(co, :s, V)
τ = 0.4 + 2.1im
b(z)
```

## Spectra

The program exposes two types to deal with CFT spectra.

The type `Spectrum` simply contains a list of conformal dimensions or fields, with conformal dimensions bounded by some ``\Delta_{\mathrm{max}}``.

The type `ChannelSpectrum` constructs all the blocks corresponding to a correlation and a list of channel fields in a given channel.

```@docs
Spectrum
ChannelSpectrum
```

## BootstrapMatrix
