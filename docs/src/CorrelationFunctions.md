# Virasoro conformal blocks

The file [CorrelationFunctions.jl](../../src/CorrelationFunctions.jl) provides structs and methods for representing and computing correlation functions.

It uses structures from [CFTdata.jl](CFTData.md).

## Four-point correlation functions on the sphere

The module `FourPointCorrelationFunctions` defines

* a struct `FourPointCorrelation` that represents a four point function
    $$
    < V_1(0) V_2(1) V_3(\infty) V_4(x)>
    $$
* a method `computeCNmn` that computes the coefficients $C^N_{m,n}$ which serve to compute the conformal blocks that enter the expansion of the 4-pt function.

The coefficient $C_{m,n}^N$ has the recursive representation
$$
C^N_{m,n} = R_{m,n}\left(\delta_{N-mn,0} + \sum_{m'n'\leq N-mn} \frac{C^{N-mn}_{m',n'}}{\delta_{(m,-n)}-\delta_{(m',n')}} \right)
$$
An expression for the $R_{m,n}$ can be found on [this wikipedia article](https://en.wikipedia.org/wiki/Virasoro_conformal_block). It can be rewritten
$$
 R_{m,n} = \frac{1}{2}\frac{1}{D_{mn}}
 \prod_{r\overset{2}{=} 1-m}^{m-1}
 \prod_{s\overset{2}{=}1-n}^{n-1}
 \sqrt{(\delta_2-\delta_1)^2 -2\delta_{(r,s)}(\delta_1+\delta_2) + \delta_{(r,s)}^2}
\sqrt{(\delta_3-\delta_4)^2 -2\delta_{(r,s)}(\delta_3+\delta_4) + \delta_{(r,s)}^2}
$$
where we do not actually take square roots, because each factor appears twice, except the $(r,s)=(0,0)$ factor which is however a perfect square. The normalization factor is
$$
D_{m,n} = mn \prod_{r=1}^{m-1} r^2B \left(r^2B - \frac{n^2}{B}\right)
\prod_{s=1}^{n-1} \frac{s^2}{B}\left(\frac{s^2}{B} - m^2B\right)
\prod_{r=1}^{m-1} \prod_{s=1}^{n-1} \left(r^2B -\frac{s^2}{B} \right)^2.
$$

### TODO: logarithmic blocks

## One-point correlation functions on the torus

The module `OnePointCorrelationFunctions` defines

* a struct `OnePointCorrelation` that represents a one point function
    $$
    < V >_\tau,
    $$
    where $\tau$ is the modulus of the torus,
* a method `computeCNmn` that computes the coefficients $C^{N,\text{torus}}_{m,n}$ which serve to compute the conformal blocks that enter the expansion of the 4-pt function.

The coefficient $C_{m,n}^{N,\text{torus}}$ has the recursive representation
$$
C^{N,\text{torus}}_{m,n} = R^{\text{torus}}_{m,n}\left(\delta_{N-mn,0} + \sum_{m'n'\leq N-mn} \frac{C^{N-mn}_{m',n'}}{\delta_{(m,-n)}-\delta_{(m',n')}} \right)
$$

where the coefficients $R_{mn}^{\text{torus}}$ can be computed from the formula in this [wikipedia article](https://en.wikipedia.org/wiki/Virasoro_conformal_block):

An expression for the $R_{m,n}$ can be found on [this wikipedia article](https://en.wikipedia.org/wiki/Virasoro_conformal_block). It can be rewritten
$$
R_{m,n}^{\text{torus}} = \frac{1}{2 D_{m,n}} \prod_{r\overset2=1-2m}^{2m-1} \prod_{s\overset2=1-2n}^{2n-1} \sqrt{\delta_{(r,s)} - \delta_1}
$$
where we do not actually take square roots, because each factor appears twice. The normalization factor is the same $D_{m,n}$ as in the [four-point](#four-point-correlation-functions-on-the-sphere) case.
