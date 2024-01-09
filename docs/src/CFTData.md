# Conformal Field Theory data

The file CFTdata.jl contains types representing objects in CFTs:
* the central charge
* fields

## Central charge

The central charge $c$ can be parametrized by the variables $B$, $b$ or $\beta$ such that
$$
c = 13 + 6B + 6 B^{-1} \quad , \quad B = b^2 = -\beta^2 
$$
Conversely, we have 
$$
B = \frac{c-13 \pm \sqrt{(c-1)(c-25)}}{12}
$$

We represent a central charge as a dictionary with 4 keys $b$,$\beta$,$B$,$c$, and the associated values of the parameters.

## Fields

Fields in 2D CFTs are elements of representations of the Virasoro algebra.

Non-logarithmic fields are elements of representations where the generators $L_0$ and $\bar L_0$ of the Virasoro algebra act diagonally, with eigenvalues called conformal weights or dimensions and respectively denoted $(\Delta,\bar \Delta)$. 
Logarithmic fields are elements of representations where $L_0$ and $\bar L_0$ act as triangular matrices; see [arXiv:2007.04190](https://arxiv.org/abs/2007.04190).

We parametrise the dimensions in terms of $P,p,\Delta,\delta$, related by

$$
\Delta = \frac{c-1}{24} + \delta  \quad , \quad \delta = -P^2 = p^2
$$

The Kac parametrisation for conformal weights is
$$ p_{(r,s)}=\frac{1}{2}(b r - b^{-1}s)$$

where $r,s$ are arbitrary numbers. We say the field is degenerate if $r,s\in \mathbb Z$ and rs>0. In terms of $r,s$, the dimension $\Delta$ is written


$$\Delta_{(r,s)} = \frac14 B (1-r^2) + \frac12 (1-rs) + \frac14\frac{1-s^2}{B}
$$

This convention is consistent with the one in [this code](https://gitlab.com/s.g.ribault/Bootstrap_Virasoro.git)  but differs from the one in [this paper](https://arxiv.org/abs/2208.14298).

In our models, non-diagonal fields are written $V_{(r,s)}$ and are parametrised by Kac indices $r$,$s$, with left and right conformal dimension $(p_{(r,s)},p_{(r,-s)})$.

We represent a field as the data of 
* a dictionary giving the values of its left and right conformal dimensions ($\Delta$, $\delta$, $p$, $P$) and ($\bar\Delta$, $\bar\delta$, $\bar p$, $\bar P$)
* a boolean saying wether the field is given in the Kac parametrisation
* the values of r and s, which are by convention set to zero if the field is not given in the Kac parametrisation
* a boolean saying whether the field is degenerate
* a boolean saying wether the field is logarithmic
* a boolean saying wether the field is diagonal