# Correlation functions and conformal blocks

The program is capable of computing chiral and non-chiral Virasoro
blocks for four-point functions on the sphere and for one-point
functions on the torus, by Zamalodchikov's recursive formula.

## Chiral blocks and Zamolodchikov's recursion formula

Blocks can be computed efficiently thanks to
[Zamolodchikov's recursion](https://en.wikipedia.org/wiki/Virasoro_conformal_block).

They are expressed as series in a variable $q$ related to $x$ through

$$\begin{cases}
x = \frac{\theta_2(q)^4}{\theta_3(q)^4}, \quad q = e^{-\pi\frac{K(1-x)}{ K(x)}} &
    \text{ for four-point blocks on the sphere} \\
q = e^{2i\pi \tau} &\text{ for one-point blocks on the torus}
\end{cases}$$

where

$$\theta_3(q) = \sum_{n\in\mathbb{Z}} q^{n^2} \quad ,
\quad \theta_2(q) = 2q^\frac14\sum_{n=0}^\infty q^{n(n+1)}$$

are Jacobi special $\theta$-functions, $K(x)$ is the elliptic $K$
function, and $\tau$ is the modulus of the torus.

```@docs
xfromq
```

```@docs
qfromx
```

### Four-point blocks on the sphere

In terms of these variables, the chiral $s$-channel sphere four point conformal block is

$$\begin{align}
\mathcal{F}^{(s)}_{\delta}(c | \Delta_{1}, \dots, \Delta_{4} | x) =
x^{E_0} (1-x)^{E_1} \theta_3(q)^{-4E_2}
(16q)^{\delta} H_{\delta}(c | \Delta_{1},\dots, \Delta_{4} | q)
\end{align}$$

where we use the exponents

$$E_0 = -\delta_1-\delta_2-\frac{c-1}{24} \quad , \quad E_1 = 
-\delta_1-\delta_4-\frac{c-1}{24} \quad ,
\quad E_2 = \delta_1+\delta_2+\delta_3+\delta_4+\frac{c-1}{24}$$

The non-trivial coefficient is the series

$$H_{\delta}(q) = 1 + \sum_{N=1}^{N_{max}} \sum_{mn\leq N} C_{m,n}^N
 \frac{(16q)^N}{\delta-\delta_{(m,n)}}$$

Where the coefficient $C_{m,n}^N$ is defined by the recursive formula

$$C^N_{m,n} = R_{m,n}\left(\delta_{N-mn,0} +
\sum_{m'n'\leq N-mn}
\frac{C^{N-mn}_{m',n'}}{\delta_{(m,-n)}-\delta_{(m',n')}} \right)$$

And the coefficients $R_{m,n}$ can be written

$$\begin{align}
 R_{m,n} = \frac{1}{2}\frac{1}{D_{mn}}
\prod_{r\overset{2}{=} 1-m}^{m-1}
\prod_{s\overset{2}{=}1-n}^{n-1}
&\sqrt{(\delta_2-\delta_1)^2 -2\delta_{(r,s)}(\delta_1+\delta_2) + \delta_{(r,s)}^2}\nonumber\\
&\sqrt{(\delta_3-\delta_4)^2 -2\delta_{(r,s)}(\delta_3+\delta_4) + \delta_{(r,s)}^2}
\end{align}$$

We do not actually take square roots, because each factor appears twice,
except the $(r,s)=(0,0)$ factor which is however a perfect square. The
normalization factor is

$$\begin{equation}
D_{m,n} = mn \prod_{r=1}^{m-1} r^2B \left(r^2B - \frac{n^2}{B}\right)
\prod_{s=1}^{n-1} \frac{s^2}{B}\left(\frac{s^2}{B} - m^2B\right)
\prod_{r=1}^{m-1} \prod_{s=1}^{n-1} \left(r^2B -\frac{s^2}{B} \right)^2.
\end{equation}$$

### One-point blocks on the torus

A similar formula holds for torus one-point blocks:

$$\begin{align}
    \mathcal F_{\Delta}(\tau, c, \Delta_{1} | x) = \frac{q^{\delta}}{\eta(q)}
    H^{\text{torus}}_{\Delta}(\tau, c, \Delta_{1} | q),
\end{align}$$

The recursion formula for
$H^{\text{torus}}_{\Delta}(\tau, c, \Delta_{1} | q)$ is

$$\begin{align}
  H_{\Delta}^{\text{torus}} (\tau, c, \Delta_{1} | q) = 1 + \sum_{N=1}^{N_{\text{max}}}
  \sum C^{N, \text{torus}}_{m,n} \frac{q^N}{\delta - \delta_{(m,n)}}
\end{align}$$

The coefficients $C_{m,n}^{N,\text{torus}}$ have the recursive
representation

$$\begin{equation}
C^{N,\text{torus}}_{m,n} = R^{\text{torus}}_{m,n}\left(\delta_{N-mn,0} +
\sum_{m'n'\leq N-mn} \frac{C^{N-mn}_{m',n'}}{\delta_{(m,-n)}-\delta_{(m',n')}} \right)
\end{equation}$$

An expression for the $R_{m,n}^{\text{torus}}$ can be found on 
[this wikipedia article](https://en.wikipedia.org/wiki/Virasoro_conformal_block).
It can be rewritten

$$R_{m,n}^{\text{torus}} = \frac{1}{2 D_{m,n}} \prod_{r\overset2=1-2m}^{2m-1}
\prod_{s\overset2=1-2n}^{2n-1} \sqrt{\delta_{(r,s)} - \delta_1}$$

where we do not actually take square roots, because each factor appears
twice. The normalization factor is the same $D_{m,n}$ as in the
four-point case.

## The `Correlation` type

To prevent recomputation of the blocks' residues, which are
independent on the channel field, the user must create a struct
`Correlation`. If passed one (resp. four) `ConformalDimension`
object(s), and an integer `Nmax`, the `Correlation` object will contain all residues of
the corresponding blocks, up to order `Nmax`.
The `Correlation` struct also serves as a convenient way to bundle
together external parameters of blocks.

```@docs
Correlation
Correlation()
```

### Regularisation of blocks

The blocks have poles at $P=P_{(r, s)}$ when $r, s \in \mathbb N^*$. We define regularised block via:

$$\mathcal F_{P_{(r, s)}} \underset{P \to P_{(r, s)}}= \frac{R_{r, s}}{P - P_{(r, s)}} \mathcal F_{P_{(r, -s)}} + \mathcal F^{\text{reg}}_{P_{(r, s)}}.$$

## Logarithmic blocks

### Logarithmic blocks on the sphere

The expression of logarithmic four-point blocks on the sphere can be
found by assuming the holomorphicity of the 4-point function

$$\begin{align}
 Z(P) = \sum_{k\in\mathbb{Z}} D_{P+k\beta^{-1}} \left|\mathcal{F}_{P+k\beta^{-1}}\right|^2 +\sum_{r=1}^\infty \sum_{s\in\frac{1}{r}\mathbb{Z}} D_{(r,s)}(P) \mathcal{G}_{(r,s)}\ .
\end{align}.$$

The coefficient $D_P$ has a double pole at $P_{(r,-s)}$. The blocks
$\mathcal F_{P}$ have a simple pole at $P_{(r,s)}$, and we write

$$\begin{align}
  \mathcal{F}_{P} = \frac{R_{r,s}}{P-P_{(r,s)}} \mathcal{F}_{P_{(r,-s)}}
   + \mathcal{F}^{\text{reg}}_{P_{(r,s)}}.
\end{align}$$

Explicitly, using Zamolodchikov's recursion, $\mathcal F^{\text{reg}}$
is written as

$$\begin{align}
  \mathcal{F}^{\text{reg}}_{P_{(r,s)}} = (\text{prefactor}) H^{\text{reg}}_{P_{(r,s)}},
\end{align}$$

where the prefactor is the prefactor in Zamolodchikov's recursion, and

$$\begin{align}
  H^{\text{reg}}_{P_{(r,s)}} = 1 + \sum_{m,n} \left( \frac{1}{P^{2}_{(r,s)} - P^{2}_{(m,n)}} \right)^{\text{reg}} (16q)^{mn} R_{m,n} H_{P_{(m,-n)}}
\end{align}$$

and

$$\begin{align}
\left(  \frac{(16q)^{P^{2}}}{P^{2}_{(r,s)} - P^{2}_{(m,n)}} \right)^{\text{reg}} =
(16q)^{P^{2}} \times
\begin{cases}
\log 16q - \frac{1}{4P_{(r,s)}^{2}} &\text{  if  } (m,n)=(r,s) \\
\frac{1}{P^{2}_{(r,s)} - P^{2}_{(m,n)}} &\text{  otherwise}
\end{cases}.
\end{align}$$

Analysing the poles of this expression (there are double poles and
simple ones), one arrives at the following expression for the
logarithmic blocks: for $(r, s) \in \mathbb{N}^{*}$,

$$\begin{align}
\mathcal{G}_{(r,s)} = (\mathcal{F}_{P_{(r,s)}}^{\text{reg}} - R_{r,s}& \mathcal{F}^{'}_{P_{(r,-s)}}) \bar{\mathcal{F}}_{P_{(r,-s)}} + \frac{R_{r,s}}{\bar R_{r,s}} \mathcal{F}_{P_{(r,-s)}} (\bar{\mathcal{F}}_{P_{(r,s)}}^{\text{reg}} - \bar{R}_{r,s} \bar{\mathcal{F}}^{'}_{P_{(r,-s)}})\nonumber \\
& +R_{r,s} \underbrace{\left( \frac{D^{'}_{P_{(r,s)}}}{D_{P_{(r,s)}}} - \lim_{P \to P_{(r,-s)}} \left[ \frac{2}{P-P_{(r,-s)}} + \frac{D_{P}^{'}}{D_{P}} \right] \right)}_{-\ell^{(1)-}_{(r,s)}}\left|\mathcal{F}_{P_{(r,-s)}}\right|^{2},
\end{align}$$

in which the primes denote derivatives with respect to the momentum $P$ (see [arXiv:1503.02067](https://arxiv.org/abs/1503.02067)).

The derivative of the block is

$$\begin{align}
  \mathcal{F}_{P_{(r,-s)}}^{'} = (\text{prefactor}) H^{\text{der}}_{P_{(r,-s)}}, \quad \text{where} \quad H^{\text{der}}_{P} = 2P\log(16q) H_{P} + H_{P}^{'}.
\end{align}$$

The term $\ell^{(1)-}_{(r,s)}$ can be computed as the order 1 term in
the Taylor expansion of

$$\begin{align}
  \log \left( \epsilon^{2} \frac{D_{P_{(r,-s)}+\epsilon}}{D_{P_{(r,s)+\epsilon}}} \right) = \sum_{n\geq 0} \ell^{(n)-}_{(r,s)} \epsilon^{n}.
\end{align}$$

Explicitly,

$$\begin{align}
 \beta\ell^{(1)-}_{(r,s)} = 4\sum_{j=1-s}^s &\Big\{ \psi(-2\beta^{-1}P_{(r,j)}) +\psi(2\beta^{-1}P_{( r,-j)}) \Big\}
 -4\pi \cot(\pi s \beta^{-2}) \nonumber
 \\
 &-\sum_{j\overset{2}{=}1-s}^{s-1}\sum_{\pm,\pm}\Big\{
 \psi\left(\tfrac12-\beta^{-1}(P_{( r,j)}\pm P_1\pm P_2)\right)
 + \psi\left(\tfrac12+\beta^{-1}(P_{( r,j)}\pm \bar P_1\pm \bar P_2)\right)
 \Big\} \nonumber
 \\
 &-\sum_{j\overset{2}{=}1-s}^{s-1}\sum_{\pm,\pm}\Big\{
 \psi\left(\tfrac12-\beta^{-1}(P_{( r,j)}\pm P_3\pm P_4)\right)
 + \psi\left(\tfrac12+\beta^{-1}(P_{( r,j)}\pm \bar P_3\pm \bar P_4)\right)
 \Big\}
\end{align}$$

For $(r, s) \in \mathbb{N}^{*}$, $\mathcal G_{(r,s)}$ can actually be
non-logarithmic, due to residues $R_{(r,s)}$ and $\bar R_{(r,s)}$
vanishing.

### Logarithmic blocks on the torus

The argument we used for computing logarithmic blocks on the sphere can
be transferred verbatim to the case of one point blocks on the torus. In
particular, the expression of the logarithmic block as a residue
is also valid for the torus one-point block, if we replace $D_P$ by the
corresponding structure constant on the torus, namely

$$\begin{align}
  D_{P} \to \frac{C^{\text{ref}}_{P,P, P_1}}{B_{P}}
\end{align}$$

where $P_1$ is the momentum of the external field.
