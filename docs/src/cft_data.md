# Basic types

The package has a few types to represent objects in a Virasoro CFT:

* `CentralCharge`s
* `ConformalDimension`s
* `Field`s

## Central charges

We parametrise the central charge of our theories in terms of
variables `c`, `B`, `b` or `β` related by

$$c = 13 + 6B + 6 B^{-1} \quad , \quad B = b^2 = -\beta^2,
\quad B = \frac{c-13 \pm \sqrt{(c-1)(c-25)}}{12}$$

By convention we keep $\beta = ib$. In $O(n)$ and $U(n)$ models
the central charge is related to $n$ via

$$n = - 2 \cos(\pi \beta^2)$$

The program allows to conveniently create central charges from any
of these four parameters, and to retrieve the value of any parameter:

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

The program lets us define these objects and access the various
parametrisations:

```@docs
ConformalDimension
P_rs
```

`ConformalDimension`s can be shifted by units of $\beta^{-1}/2$:

```@docs
shift
```

## Fields

For our purposes, a field is the data of left and right conformal dimensions.
We denote $V_{(r, s)}$ a field with left and right dimensions
$(\Delta_{(r, s)}, \Delta_{(r, -s)})$.

The program exposes a `Field` struct and convenient constructors:

```@docs
Field
isdiagonal
isdegenerate
shift
```

## Spectra

The program exposes two types to deal with CFT spectra.

The type `Spectrum` simply contains a list of conformal dimensions or fields, with conformal dimensions bounded by some ``\Delta_{\mathrm{max}}``.

The type `ChannelSpectrum` contains all the data needed to evaluate blocks corresponding to a correlation and a list of channel fields in a given channel.

```@docs
Spectrum
ChannelSpectrum
```
