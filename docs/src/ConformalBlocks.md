# Virasoro conformal blocks

The file [ConformalBlocks.jl](../../src/ConformalBlocks.jl) implements Zamolodchikov's recursion formula for computing four-point conformal blocks on the sphere and one-point conformal blocks on the torus.

It uses structures from [CFTdata.jl](CFTData.md).

## Four-point conformal blocks on the sphere

The module `FourPointBlocksSphere` exports

* a struct `FourPointBlockSphere` that encapsulates the data needed to compute a 4pt conformal block, namely a channel, four external fields and the field propagating in the channel
* a function `F_four_point_sphere(block, charge, x)` which computes the value of the non-chiral block $\mathcal F_{\Delta}^{(s)}(\Delta_i | x)$ as defined in the paragraph [Background Material](#four-point-blocks-on-the-sphere) below.

## One-point conformal blocks on the torus

The module `OnePointBlocksTorus` exports

* a struct `OnePointBlockTorus` that encapsulates the data needed to compute a 4pt conformal block, namely an external field.
* a function `F_one_point_torus(block, q)` which computes the value of the non-chiral block $\mathcal F_{\Delta}^{\text{torus}}(\Delta | q)$ as defined in the paragraph [Background Material](#one-point-blocks-on-the-torus) below.

## Background material

### Four-point blocks on the sphere

Computation of any four-point correlation function reduces to the computation of

$$ \mathcal F(x) = \langle V_{1}(x) V_{2}(0) V_{3}(\infty) V_{4}(1) \rangle $$

Conformal blocks encode the universal part of correlation functions. More precisely, by performing the OPE of $V_{1}$ with $V_{2}$ and $V_{3}$ with $V_{4}$ (s-channel), or $1\leftrightarrow 4, 2\leftrightarrow 3$ (t-channel) or $1\leftrightarrow 3, 2\leftrightarrow 4$ (u-channel), $\mathcal F(x)$ can be written as

$$ \mathcal F(x) = \sum_{s \in \mathcal S^{(s)}} \frac{C_{12s}C_{s34}}{B_s} \mathcal F_\Delta^{(s)}(\Delta_i | x)$$

where $\sum_{\mathcal S^{(s)}}$ is a sum over a basis of the $s$-channel spectrum $\mathcal S^{(s)}$, $B_s$ are two-point structure constants and $C_{ijk}$ are three-point structure constants. Analogous expressions can be written for the $t$- and $u$- channels. The functions $\mathcal F_{\Delta}^{(s)}(\Delta_i | x)$ are called $s$-channel conformal blocks (resp. t,u).

Conformal blocks factorise into holomorphic and anti-holomorphic parts, and they are characterized by the normalisation conditions

$$
\begin{align}
 \mathcal{F}^{(s)}_\Delta(x) & \underset{x\to 0}{=} \left| x^{\Delta-\Delta_1-\Delta_2}\right|^2 \left(1+O(x)\right) \nonumber
 \\
 \mathcal{F}^{(t)}_\Delta(x) & \underset{x\to 1}{=} \left|(1-x)^{\Delta-\Delta_1-\Delta_4}\right|^2 \left(1+O(1-x)\right) \nonumber
 \\
 \mathcal{F}^{(u)}_\Delta(x) & \underset{x\to \infty}{=} \left|\left(\frac{1}{x}\right)^{\Delta+\Delta_1-\Delta_3} \right|^2\left(1+O\left(\frac{1}{x}\right)\right)
\end{align}
$$
(we omit the $\Delta_i$ dependence in the notation $\mathcal{F}^{(u)}_\Delta(x)$).

Together with the invariance of $\left\langle \prod_{i=1}^4 V_{\Delta_i}(z_i) \right\rangle$ under permutations, this leads to the relations
$$
\begin{align*}
\mathcal{F}^{(t)}_{\Delta}(\Delta_1,\Delta_2,\Delta_3,\Delta_4|x)
&= (-1)^{S_1+S_2+S_3+S_4}
\mathcal{F}^{(s)}_{\Delta}(\Delta_1,\Delta_4,\Delta_3,\Delta_2|1-x)
\\
\mathcal{F}^{(u)}_\Delta(\Delta_1,\Delta_2,\Delta_3,\Delta_4|x)
&= (-1)^{S_1+S_2+S_3+S_4}
\left|x^{-2\Delta_1}\right|^2 \mathcal{F}^{(s)}_\Delta(\Delta_1,\Delta_3,\Delta_2,\Delta_4|\tfrac{1}{x})
\end{align*}
$$
where $S=\Delta-\bar\Delta$ is the conformal spin, which we assume to be integer.
<!-- If some dimensions and spectra are equal, we may assume that some spectra are purely even/odd spin, and impose relations between four-point structure constants:
$$
 \begin{array}{|l|l|l|}
 \hline
 \text{Equal dimensions} & \text{Even spin sector} & \text{Odd spin sector}
 \\
 \hline 
 1,2 \text{ or }  3, 4& D^{(t)} = (-)^{S+S_2+S_3} D^{(u)} & D^{(t)} = -(-)^{S+S_2+S_3} D^{(u)}
  \\
  \hline 
  2,4 \text{ or } 1, 3 & D^{(s)} = (-)^{S_\text{tot}} D^{(t)} & D^{(s)} = -(-)^{S_\text{tot}} D^{(t)}
  \\
  \hline 
  2,3 & D^{(s)} = (-)^{S_\text{tot}} D^{(u)} & D^{(s)} = -(-)^{S_\text{tot}} D^{(u)}
  \\
  \hline 
  1,4 & D^{(s)} = D^{(u)} & D^{(s)} = - D^{(u)}
  \\
  \hline 
  1,2,3 \text{ or } 2,3,4 & D^{(t)}=D^{(u)}= (-)^{S_\text{tot}} D^{(s)} & D^{(t)}=D^{(u)}= -(-)^{S_\text{tot}}D^{(s)}
  \\
  \hline 
  1,3,4 \text{ or } 1,2,4 & D^{(s)} = D^{(u)} = (-)^{S_\text{tot}} D^{(t)} & D^{(s)} = -D^{(u)} = - (-)^{S_\text{tot}} D^{(t)}
  \\
  \hline 
  1,2,3,4 & D^{(s)} = D^{(t)} = D^{(u)} & D^{(s)}=-D^{(t)} = -D^{(u)} 
  \\
  \hline 
 \end{array}
 $$ -->

#### Notations

In addition to $x$, we use the elliptic nome $q$ such that
$$
x = \frac{\theta_2(q)^4}{\theta_3(q)^4}, \quad q = e^{-\pi\frac{K(1-x)}{ K(x)}}
$$
where
$$
\theta_3(q) = \sum_{n\in\mathbb{Z}} q^{n^2} \quad , \quad \theta_2(q) = 2q^\frac14\sum_{n=0}^\infty q^{n(n+1)}
$$
are Jacobi special $\theta$-functions, and $K(x)$ is the elliptic $K$ function.

#### Expression

Our chiral $s$-channel conformal block is
$$
 \mathcal{F}^{(s)}_{\delta}(x) =  x^{E_0} (1-x)^{E_1} \theta_3(q)^{-4E_2}
(16q)^{\delta} H(q,\delta)
$$
where we use the exponents
$$
E_0 = -\delta_1-\delta_2-\frac{c-1}{24} \quad , \quad E_1 = -\delta_1-\delta_4-\frac{c-1}{24} \quad ,
\quad E_2 = \delta_1+\delta_2+\delta_3+\delta_4+\frac{c-1}{24}
$$
The non-trivial coefficient is the series
$$
H(q,\delta) = 1 + \sum_{N=1}^{N_{max}} \sum_{mn\leq N} C_{m,n}^N \frac{(16q)^N}{\delta-\delta_{(m,n)}}
$$

<!-- If $R_{m,n}=0$, we compute a finite regularization of $R_{m,n}$. This is never used for computing chiral conformal blocks, but only for computing non-chiral logarithmic blocks. The regularization of vanishing factors is
$$
\left(\delta_2-\delta_1\right)_\text{reg} = 2p_2
$$
$$
\left((\delta_2-\delta_1)^2 -2\delta_{(r,s)}(\delta_1+\delta_2) + \delta_{(r,s)}^2\right)_\text{reg} = 8p_1p_2p_{(r,s)}
$$ -->

#### Degenerate blocks

In the case $\delta_1 = \delta_{(r_1,s_1)}$ and $\delta_2 = \delta_{(r_2,s_2)}$ with $r_i,s_i\in\mathbb{N}^*$, we have

$$
\left\{\begin{array}{l} m\in |r_1-r_2|+1+2\mathbb{N}
\\ n \in |s_1-s_2| + 1+2\mathbb{N} \end{array} \right. \quad
\Rightarrow \quad R_{m,n} = 0
$$
and similarly if the fields with numbers $3$ and $4$ are degenerate. Thanks to $R_{m,n}=0$ thus $C_{m,n}^N=0$, the block $\mathcal{F}^{(s)}_{\delta_{(m,n)}}(x)$ is finite and can be computed exactly.

This can be generalized to fractional indices $r_i,s_i$. In this case, we have to add the following restriction, which was redundant for positive integer indices:
$$
\left\{\begin{array}{l} m\in |r_1+r_2|+1+2\mathbb{N}
\\ n \in |s_1+s_2| + 1+2\mathbb{N} \end{array} \right. \quad
\Rightarrow \quad R_{m,n} = 0
$$
and similarly if the fields with numbers $3$ and $4$ are degenerate. In particular, for $\delta_i = \delta_{(0,\frac12)}$, we have $m\in 2\mathbb{N}+1\Rightarrow R_{m,n}=0$.

#### Derivative and regularization

For the purpose of computing conformal blocks for logarithmic channel representations, we need to compute derivatives of conformal blocks with respect to the channel dimension, and regularized values of blocks at their poles. Taking the derivative amounts to
$$
H(q, \delta) \to \log(16q) H(q,\delta) +\frac{\partial}{\partial\delta} H(q, \delta)
$$
And the regularization we are interested in amounts to
$$
\left.\frac{1}{\delta-\delta_{(m,n)}}\right|_{\delta=\delta_{(m,n)}} = \log(16q)-\frac{1}{4\delta_{(m,n)}}
$$
The code can formally compute a regularization of the block's derivative, but this regularization is a priori not meaningful.

### One point blocks on the torus

A chiral one-point conformal block on the torus can be written as

$$
\mathcal{F}^\text{torus}_{\Delta}(\Delta_1|q) = q^{\Delta-\frac{c-1}{24}}\eta(q)^{-1}H^\text{torus}_{\Delta}(\Delta_1|q)
$$

where $\eta$ is the Dedekind eta function, $q$ is related to the torus' modular parameter $\tau$ via $q=e^{2i\pi \tau}$, $c$ is the central charge and the non-trivial factor $H^{\text{torus}}_{\Delta}(\Delta_1|q)$ can be computed by a recursion formula:

$$
H^\text{torus}_{\Delta}(\Delta_1|q) = 1 + \sum_{N=1}^{N_{max}} \sum_{mn\leq N} C_{m,n}^{N,\text{torus}} \frac{q^N}{\delta-\delta_{(m,n)}}
$$
