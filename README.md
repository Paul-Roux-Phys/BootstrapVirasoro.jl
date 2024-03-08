#+title: JuliVirBootstrap Documentation
#+author: Paul ROUX
#+setupfile: https://fniessen.github.io/org-html-themes/org/theme-readtheorg.setup
#+options: toc:3 num:1
#+startup: latexpreview
#+property: header_args:julia :session JuliVirBootstrap :eval no-export

* Table of contents :toc:noexport:
- [[#main-module][Main module]]
- [[#cft-data-central-charges-and-fields][CFT data: central charges and fields]]
  - [[#parametrisations][Parametrisations]]
  - [[#the-cftdata-module][The ~CFTData~ module]]
- [[#correlation-functions][Correlation Functions]]
  - [[#four-point-correlation-functions-on-the-sphere][Four-point correlation functions on the sphere]]
  - [[#one-point-correlation-functions-on-the-torus][One-point correlation functions on the torus]]
- [[#virasoro-conformal-blocks][Virasoro conformal blocks]]
  - [[#four-point-conformal-blocks-on-the-sphere][Four-point conformal blocks on the sphere]]
  - [[#one-point-conformal-blocks-on-the-torus][One-point conformal blocks on the torus]]
- [[#special-functions][Special functions]]
- [[#unit-testing][Unit testing]]
- [[#development-tests][Development tests]]
  - [[#relation-between-four-point-blocks-on-the-sphere-and-one-point-blocks-on-the-torus][Relation between four-point blocks on the sphere and one-point blocks on the torus]]

* Main module
:PROPERTIES:
:header-args:julia: :tangle ./src/JuliVirBootstrap.jl
:END:

The module ~JuliVirBootstrap~ is the main module of this package, and it includes the sub-modules.

- ~CFTData~ provides types for central charges and fields.
- ~CorrelationFunctions~ provides types for one-point and four-point correlation functions, as well as methods for computing coefficients appearing in their conformal blocks.
- ~VirasoroConformalBlocks~ provides types for representing four-point conformal blocks on the sphere and one-point conformal blocks on the torus, as well as methods for computing them.

#+begin_src julia
#===========================================================================================

Written by Paul Roux, adapting a Python code written by Sylvain Ribault & Rongvoram
Nivesvivat

===========================================================================================#

module JuliVirBootstrap

#===========================================================================================
Central charges and fields
===========================================================================================#
include("CFTData.jl")
using .CFTData
export CentralCharge, Field

#===========================================================================================
Correlation functions
===========================================================================================#
include("CorrelationFunctions.jl")
using .FourPointCorrelationFunctions
export FourPointCorrelation

using .OnePointCorrelationFunctions
export OnePointCorrelation

#===========================================================================================
Conformal blocks
===========================================================================================#
include("ConformalBlocks.jl")
using .FourPointBlocksSphere
export FourPointBlockSphere, F_four_point_sphere

using .OnePointBlocksTorus
export OnePointBlockTorus, F_one_point_torus


#===========================================================================================
Special functions
===========================================================================================#
include("SpecialFunctions.jl")
export log_double_gamma, double_gamma


end
#+end_src

* CFT data: central charges and fields
:PROPERTIES:
:header-args:julia: :tangle ./src/CFTData.jl
:END:

The file [[file:src/CFTData.jl][CFTData.jl]] defines

- a struct ~CentralCharge~ that represents a central charge $c$ and contains the value of the four corresponding parameters $b, B, \beta, c$
- a struct ~Field~ that represents a field $V$. The field can be defined from its Kac indices $r, s$, be diagonal, logarithmic, or degenerate. The struct contains booleans for these three characteristics, as well as rationals for $r$ and $s$, and the pairs of (left, right) values $(\Delta, \bar \Delta)$, $(p, \bar p)$, $(\delta, \bar \delta)$, $(P, \bar P)$.

** Parametrisations

*** Central charge

The central charge \(c\) can be parametrized by variables \(B\), \(b\) or \(\beta\) such that

\[
c = 13 + 6B + 6 B^{-1} \quad , \quad B = b^2 = -\beta^2
\]

Conversely, we have

\[
B = \frac{c-13 \pm \sqrt{(c-1)(c-25)}}{12}
\]

The central charge can be given directly or as a function of b,β,B. No matter what parameter is passed we compute B from it, then b,β,c. We ensure that the relation between b and \beta never changes (square root branch).

*** Fields

Fields in 2D CFTs are elements of representations of the Virasoro algebra.
Non-logarithmic fields are elements of representations where the generators \(L_0\) and \(\bar L_0\) of the Virasoro algebra act diagonally, with eigenvalues called conformal weights or dimensions and respectively denoted \((\Delta,\bar \Delta)\). Logarithmic fields are elements of representations where \(L_0\) and \(\bar L_0\) act as triangular matrices; see [[https://arxiv.org/abs/2007.04190][this paper]].
We parametrise the dimensions in terms of \(P,p,\Delta,\delta\), related by

\[
\Delta = \frac{c-1}{24} + \delta  \quad , \quad \delta = -P^2 = p^2
\]

The Kac parametrisation for conformal weights is

\[P_{(r,s)}=\frac{1}{2}(b r + b^{-1}s)\]

Or equivalently

\[p_{(r,s)} = -\frac{1}{2} (\beta r - \beta^{-1}s)\]

where \(r,s\) are arbitrary numbers. We say the field is degenerate if \(r,s\in \mathbb Z\) and $rs > 0$. In terms of \(r,s\), the dimension \(\Delta\) is written

\[\Delta_{(r,s)} = \frac14 B (1-r^2) + \frac12 (1-rs) + \frac14\frac{1-s^2}{B}\]

This convention is consistent with the one in [[https://gitlab.com/s.g.ribault/Bootstrap_Virasoro.git][this code]] but differs from the one in [[https://arxiv.org/abs/2208.14298][this paper]].
In our models, non-diagonal fields are written \(V_{(r,s)}\) and are parametrised by Kac indices \(r\),\(s\), with left and right conformal dimension \((P_{(r,s)},P_{(r,-s)})\).

** The ~CFTData~ module

*** Header

#+begin_src julia

#===========================================================================================

CFTData.jl contains a module CFTData that provides types representing
central charges and fields in 2D CFTs with Virasoro symmetry.

Written by Paul Roux, adapting a Python code written by Sylvain Ribault & Rongvoram
Nivesvivat

============================================================================================#

"""
Provides types representing central charges and fields in CFT.
"""
module CFTData

using Match;

export CentralCharge, Field

"""print complex numbers in latex format"""
function Base.show(io::IO,::MIME"text/latex",z::Complex)
    print("$(real(z)) + $(imag(z))i")
end

#+end_src

*** Central charge

#+begin_src julia
"""Get B from given parameter"""
function Bfrom(parameter, value)
    @match parameter begin
        "c" => (value-13+sqrt(complex((value-1)*(value-25))))/12
        "b" => value^2
        "β" => -value^2
        "B" => value
    end
end

"""Get asked parameter from B"""
function Bto(parameter, value)
    @match parameter begin
        "c" => 13+6*value+6/value
        "b" => sqrt(complex(value))
        "β" => im*sqrt(complex(value))
        "B" => value
    end
end

"""
    CentralCharge{T}
Object representing the central charge.
Contains the values of the 4 parameters representing it.
"""
struct CentralCharge{T}

    #= T is the type of the parameters; either Complex{Float64} or Complex{BigFloat}
    for arbitrary precision. =#
    values::Dict{String, T}

end

"""
    CentralCharge(parameter, value)

Constructor function for the CentralCharge type.

Given one of the four parameters `"c"`, `"b"`, `"β"`, `"B"` and its value,
creates an object CentralCharge{T} where T is the type of `value`.

# Example
```julia-repl
julia> setprecision(BigFloat, 20, base=10)
julia> charge = CentralCharge("β", sqrt(big(2)))
Central charge :
B = -2.0 + 0.0im
c = -2.0 + 0.0im
b = 0.0 + 1.414213562373095048804im
β = -1.414213562373095048804 + 0.0im
```
"""
function CentralCharge(parameter = "c", value = 1)
    # Constructor
    T=typeof(AbstractFloat(real(value)))
    B=Bfrom(parameter, value)
    dict=Dict(key => Bto(key, B) for key in ("c", "b", "β", "B"))
    CentralCharge{complex(T)}(dict)
end
#+end_src

**** Pretty printing

#+begin_src julia
"""Display an object of type CentralCharge"""
function Base.show(io::IO, charge::CentralCharge)
    println("Central charge:")
    for (key, value) in charge.values
        println(io, "$key = $value")
    end
end

"""Display the value of the central charge in LaTeX format"""
function Base.show(io::IO, ::MIME"text/latex", charge::CentralCharge, parameter)
    if parameter=="β"
        print("\\beta = ")
    else
        print(parameter," = ")
    end
    show(io, MIME("text/latex"), charge[parameter])
end

"""Overload of [] to access values in charge"""
Base.getindex(charge::CentralCharge, key) = charge.values[key];
#+end_src

*** Fields

Fields can be given from any of the four parameters $\Delta, \delta, P, p$. Optional keyword arguments lets us choose whether the field is diagonal, degenerate, logarithmic. The field can also be defined from its r and s indices using the keyword argument Kac = true.

#+begin_src julia
"""Get p from any given parameter"""
function p_from(parameter, value, charge::CentralCharge)
    @match parameter begin
        "Δ" => sqrt(complex(value - (charge["c"]-1)/24))
        "δ" => sqrt(complex(value))
        "P" => -im*value
        "p" => value
    end
end

"""Get all parameters from p"""
function p_to(parameter, value, charge::CentralCharge)
    @match parameter begin
        "Δ" => value^2 + (charge["c"]-1)/24
        "δ" => value^2
        "P" => im*value
        "p" => value
    end
end

"""
    Field{T}
Object representing a conformal field.
Contains the values of the 4 parameters `"Δ"`,`"δ"`,`"P"`,`"p"` for its conformal dimension,
and flags saying whether the field is in the Kac table, degenerate, logarithmic or diagonal.
"""
struct Field{T}

    values::Dict{String, Vector{T}}
    isKac::Bool
    r::Rational
    s::Rational
    isdegenerate::Bool
    islogarithmic::Bool
    isdiagonal::Bool

end

"""
    Field(charge, parameter, leftvalue, rightvalue; kwargs...)

Constructor function for the Field type.

Given a charge `charge`, one of the four parameters `"Δ"`, `"δ"`, `"P"`, `"p"` and two values,
create an object Field{T} (where T is the type of the values in `charge`) that represents a
field of left and right dimensions given by leftvalue and rightvalue in the chosen
parametrisation.

# keyword arguments:

- `Kac::Bool`: if set to true, the field can be constructed from the values of its r and s
indices,
- `r::Rational`,`s::Rational`: used in conjunction to `Kac=true`, must be given rational
values,
- `logarithmic::Bool`: set to True if the field is logarithmic,
- `degenerate::Bool`: set to True if the field is degenerate,
- `diagonal::Bool`: set to True to get a diagonal field ; only the leftvalue needs to be
given.

# Examples
```julia-repl
julia> charge = CentralCharge("b", big(0.5));
julia> field = Field(charge, Kac=true, r=0, s=1)
Non-diagonal field with Kac indices r = 0//1, s = 1//1 and (left,right) dimensions:
Δ = ( 2.5625 + 0.0im, 2.5625 + 0.0im )
P = ( -0.0 - 1.0im, 0.0 + 1.0im )
δ = ( 1.0 - 0.0im, 1.0 + 0.0im )
p = ( -1.0 + 0.0im, 1.0 + 0.0im )
```
```julia-repl
julia> charge = CentralCharge("β", 1.5+im);
julia> Field(charge, "δ", 2, 3)
Non-diagonal field with (left, right) dimensions:
Δ = ( 2.1579142011834325 - 0.6789940828402367im, 3.1579142011834316 - 0.6789940828402367im )
P = ( 0.0 + 1.4142135623730951im, 0.0 + 1.7320508075688772im )
δ = ( 2.0000000000000004 + 0.0im, 2.9999999999999996 + 0.0im )
p = ( 1.4142135623730951 + 0.0im, 1.7320508075688772 + 0.0im )
```
```julia-repl
julia> charge = CentralCharge();
julia> Field(charge, "δ", 1, diagonal=true)
Diagonal field of dimension:
Δ = 1.0 + 0.0im
P = 0.0 + 1.0im
δ = 1.0 + 0.0im
p = 1.0 + 0.0im
```
"""
function Field(
    charge::CentralCharge = CentralCharge("c", 1),
    parameter = "Δ",
    leftvalue = 0, rightvalue = 0;
    Kac = false, r = 0, s = 0,
    logarithmic = false, degenerate = false, diagonal = false
    )

    T=typeof(charge.values["c"])   #dimensions have the same type as central charges
    if degenerate
        Kac = true
    end
    if Kac
        pleft = -1/2*(charge["β"]*r - 1/charge["β"]*s)
        pright = -1/2*(charge["β"]*r - 1/charge["β"]*s)
    else
        pleft, pright = p_from.(parameter, [leftvalue, rightvalue], Ref(charge))
    end
    if diagonal
        pright = pleft
    end
    values = Dict(key => p_to.(key, [pleft, pright], Ref(charge))
                  for key in ("Δ", "δ", "P", "p")
                      )
    Field{complex(T)}(values, Kac, r, s, degenerate, logarithmic, diagonal)
end

"""Compute the spin Δleft - Δright of a field."""
function spin(field::Field)
    #Computes the spin Δ-Δbar
    if field.isdiagonal
        return 0
    else
        return field["Δ"][1] - field["Δ"][2]
    end
end
#+end_src

**** Pretty printing

#+begin_src julia
"""Display field"""
function Base.show(io::IO,field::Field)
    #Print fields
    if field.isdiagonal
        println("Diagonal field of dimension:")
        for (key, value) in field.values
            println(io, "  $key = $(value[1])")
        end
    else
        print("Non-diagonal field ")
        if field.isKac
            print("with Kac indices r = $(field.r), s = $(field.s) and ")
        else
            print("with ")
        end
        println("(left, right) dimensions:")
        for (key, value) in field.values
            println(io, "  $key = ($(value[1]), $(value[2]))")
        end
    end
end

"""Display dimension of field in latex format"""
function Base.show(io::IO,::MIME"text/latex", field::Field,parameter)
    if field.isdiagonal
        if parameter == "Δ"
            print("\\Delta = ")
        elseif parameter == "δ"
            print("\\delta = ")
        else
            print(parameter," = ")
        end
        show(io, MIME("text/latex"), field[parameter][1])
    else
        if parameter=="Δ"
            print("(\\Delta, \\bar\\Delta) = ")
        elseif parameter=="δ"
            print("(\\delta, \\bar\\delta) = ")
        else
            print("($parameter, \\bar$parameter) = ")
        end
        print("("); show(io, MIME("text/latex"), field[parameter][1]); print(", ");
        show(io, MIME("text/latex"), field[parameter][2]); print(")")
    end
end

# function Base.show(io::IO, arr::Vector{Field{T}}) where {T}
#     println(io, "Vector{Field{$T}} with $(length(arr)) elements:")
#     for (index, field) in enumerate(arr)
#         print(io, "$(index): ")
#         show(io, field)
#         println()
#     end
# end

"""Overload []"""
Base.getindex(field::Field,key) = field.values[key];
#+end_src

*** End of module

#+begin_src julia
end # end module
#+end_src

* Correlation Functions
:PROPERTIES:
:header-args:julia: :tangle ./src/CorrelationFunctions.jl
:END:

The file [[./src/CorrelationFunctions.jl][CorrelationFunctions.jl]] provides structs and methods for representing and computing correlation functions.

It uses the types defined in [[The ~CFTData~ module]].

** Four-point correlation functions on the sphere

*** Four point functions in CFT

Because of conformal invariance, computation of any four-point correlation function reduces to the computation of

\[ \mathcal F(x) = \langle V_{1}(x) V_{2}(0) V_{3}(\infty) V_{4}(1) \rangle \]

*** Coefficients $C^{N}_{m,n}$

We define a coefficient \(C_{m,n}^N\) by the recursive formula

\[
C^N_{m,n} = R_{m,n}\left(\delta_{N-mn,0} + \sum_{m'n'\leq N-mn} \frac{C^{N-mn}_{m',n'}}{\delta_{(m,-n)}-\delta_{(m',n')}} \right)
\]
 
An expression for the \(R_{m,n}\) can be found in [[https://en.wikipedia.org/wiki/Virasoro_conformal_block][this wikipedia article]]. It can be rewritten

\[
R_{m,n} = \frac{1}{2}\frac{1}{D_{mn}}
\prod_{r\overset{2}{=} 1-m}^{m-1}
\prod_{s\overset{2}{=}1-n}^{n-1}
\sqrt{(\delta_2-\delta_1)^2 -2\delta_{(r,s)}(\delta_1+\delta_2) + \delta_{(r,s)}^2}
\sqrt{(\delta_3-\delta_4)^2 -2\delta_{(r,s)}(\delta_3+\delta_4) + \delta_{(r,s)}^2}
\]

where we do not actually take square roots, because each factor appears twice, except the \((r,s)=(0,0)\) factor which is however a perfect square. The normalization factor is

#+name: Dmn
\[
D_{m,n} = mn \prod_{r=1}^{m-1} r^2B \left(r^2B - \frac{n^2}{B}\right)
\prod_{s=1}^{n-1} \frac{s^2}{B}\left(\frac{s^2}{B} - m^2B\right)
\prod_{r=1}^{m-1} \prod_{s=1}^{n-1} \left(r^2B -\frac{s^2}{B} \right)^2.
\]

The $C^{N}_{m,n}$ only depend on the external fields and the central charge, hence we compute them in the ~FourPointCorrelationFunctions~ module.

*** The =FourPointCorrelationFunctions= module

The module =FourPointCorrelationFunctions= defines

- a struct =FourPointCorrelation= that represents a four point function
  
  \[
  < V_1(0) V_2(1) V_3(\infty) V_4(x)>
  \]

- a method =computeCNmn= that computes the coefficients \(C^N_{m,n}\) which serve to compute the conformal blocks that enter the expansion of the 4-pt function.

**** Header

#+begin_src julia
#===========================================================================================

Written by Paul Roux, adapting a Python code written by Sylvain Ribault & Rongvoram
Nivesvivat

===========================================================================================#

module FourPointCorrelationFunctions

export FourPointCorrelation, computeCNmn

using ..CFTData
using Match
import Memoization: @memoize
#+end_src

**** Four-point function type

We create a struct ~FourPointCorrelation~ for representing a four-point function on the sphere, that is, a central charge and four external fields.

#+begin_src julia
"""Struct representing a four-point function. Contains
- a central charge
- 4 external fields
"""
struct FourPointCorrelation{T}
    charge::CentralCharge{T}
    fields::Vector{Field{T}}
end

function FourPointCorrelation(charge::CentralCharge{T}, V1, V2, V3, V4) where {T}
    return FourPointCorrelation{T}(charge, [V1, V2, V3, V4])
end

"""Display a four-point function"""
function Base.show(io::IO, corr::FourPointCorrelation)
    println("Four-point correlation function: < V_1 V_2 V_3 V_4 > where ")
    print("V_1 = "); show(corr.fields[1])
    print("V_2 = "); show(corr.fields[2])
    print("V_3 = "); show(corr.fields[3])
    print("V_4 = "); show(corr.fields[4])
end

# explicit names for the indices of left and right dimensions
const left = 1
const right = 2
#+end_src

**** Compute $C^N_{m,n}$

The function ~permute_ext_fields~ permutes the external fields such that the first two and last two are fused together in the channel.

The ~@memoize~ macro stores the result of the function such that subsequent calls with the same arguments only require a memory access.

The function ~Rmn_zero_order~ computes the order of a zero of R, to avoid computing 0/0 in $\frac{R_{m,n}}{\delta - \delta_{r,s}}$. At generic central charge (non-rational) $R_{m,n}$ is zero iff one of the two pairs of fused fields have Kac indices such that $r_1 \pm r_2 \in \{1-m, 3-m, \dots, m-1\}$ or $s_1 \pm s_2 \in \{1-n, 3-n, \dots, n-1\}$.

#+begin_src julia
double_prod_in_Dmn(m, n, B) = prod(prod((r^2*B - s^2/B)^2 for s in 1:n-1) for r in 1:m-1)

δrs(r, s, B) = -1/4 * (B*r^2 + 2*r*s + s^2/B)

function Dmn(m, n, B)
    if m == 1 && n == 1 # treat cases m = 1, n=1 separately
        return 1
    elseif m == 1
        return n * prod(s^2/B * (s^2/B - m^2*B) for s in 1:n-1)
    elseif n == 1
        return m * prod(r^2*B * (r^2*B - n^2/B) for r in 1:m-1)
    else
        f1 = prod(r^2*B * (r^2*B - n^2/B) for r in 1:m-1)
        f2 = prod(s^2/B * (s^2/B - m^2*B) for s in 1:n-1)
        f3 = double_prod_in_Dmn(m, n, B)
        return m*n*f1*f2*f3
    end
end

"""Permute the external fields to get t- or u-channels from s-channel"""
function permute_ext_fields(corr::FourPointCorrelation, channel)
    Vs=corr.fields
    Vs = @match channel begin
        "s" => [Vs[1], Vs[2], Vs[3], Vs[4]]
        "t" => [Vs[1], Vs[4], Vs[3], Vs[2]]
        "u" => [Vs[1], Vs[3], Vs[2], Vs[4]]
    end
    return FourPointCorrelation(corr.charge, Vs)
end

"""Order of a pole of Rmn, assuming the central charge is generic"""
function Rmn_zero_order(m, n, B, corr::FourPointCorrelation, channel)

    order=0
    V=permute_ext_fields(corr, channel).fields

    if !((V[1].isKac && V[2].isKac) || (V[3].isKac && V[4].isKac))
        return 0
    end

    r=[V[i].r for i in 1:4]
    s=[V[i].s for i in 1:4]

    #= Rmn is zero if r1 \pm r2 or r3 \pm r4 is an integer in 1-m:2:m-1, and
    s1 \pm s2 or s3 \pm s4 is an integer in 1-n:2:n-1.
    equivalently, if (|r1 \pm r2| <= m-1 and r1-r2 - (m-1) % 2 == 0)
    and (|s1 \pm s2| <= n-1 and s1-s2 - (n-1) % 2 == 0)
    =#
    for pm in (-1,1)
        for (i,j) in ((1,2),(3,4))
            if V[i].isdegenerate && V[j].isdegenerate
                if (abs(r[i]+pm*r[j]) <= m-1 && (r[i]+pm*r[j]-(m-1))%2 == 0) &&
                    (abs(s[i]+pm*s[j]) <= n-1 && (s[i]+pm*s[j]-(n-1))%2 == 0)
                    order += 1
                end
            end
        end
    end

    return order
end

function helper_Rmn(δ1, δ2, δ3, δ4, r, s, B)
    if r == 0 && s == 0
        return (δ2-δ1)*(δ3-δ4)
    else
        return (((δ2-δ1)^2 - 2*δrs(r, s, B)*(δ1+δ2) + δrs(r, s, B)^2)
                ,*((δ3-δ4)^2 - 2*δrs(r, s, B)*(δ3+δ4) + δrs(r, s, B)^2))
    end
end

"""
Compute `Rmn`.
lr indicates the left or right moving parts of the fields
Cache the result.
TODO: value of regularisation
"""
@memoize function Rmn(m, n, corr::FourPointCorrelation, channel, lr)

    B = corr.charge["B"]
    Vs = permute_ext_fields(corr, channel).fields
    δ1 = Vs[1]["δ"][lr]
    δ2 = Vs[2]["δ"][lr]
    δ3 = Vs[3]["δ"][lr]
    δ4 = Vs[4]["δ"][lr]

    if Rmn_zero_order(m, n, B, corr, channel) > 0
        return 0
    else
        if m == 1
            res = prod(helper_Rmn(δ1, δ2, δ3, δ4, 0, s, B) for s in 1-n:2:0)
        else # m > 1
            res = prod(prod(helper_Rmn(δ1, δ2, δ3, δ4, r, s, B)
                            for s in 1-n:2:n-1) for r in 1-m:2:-1)
            if m%2 == 1 # m odd -> treat r=0 term separately
                res *= prod(helper_Rmn(δ1, δ2, δ3, δ4, 0, s, B) for s in 1-n:2:0)
            end
        end
    end

    return res/(2*Dmn(m, n, B))
end

@memoize function computeCNmn(N, m, n, corr::FourPointCorrelation, channel, lr)
    B = corr.charge["B"]
    if Rmn_zero_order(m, n, B, corr, channel) > 0
        return 0
    elseif m*n > N
        return 0
    elseif m*n == N
        return Rmn(m, n, corr, channel, lr)
    else
        res = sum(sum(computeCNmn(N-m*n, mp, np, corr, channel, lr)/(δrs(m, -n, B) - δrs(mp, np, B))
                      for mp in 1:N-m*n if mp*np <= N-m*n)
                  for np in 1:N-m*n)
        return Rmn(m, n, corr, channel, lr) * res
    end
end
#+end_src

**** End module

#+begin_src julia
end # end module
#+end_src



** One-point correlation functions on the torus

*** One point functions in CFT

Because of translation invariance, one-point functions on the torus do not depend on the field's position, but only on the modulus of the torus.

*** Coefficients $C^{N,\text{torus}}_{m,n}$

The coefficients \(C_{m,n}^{N,\text{torus}}\) have the recursive representation

#+name: CNmn-torus 
\begin{equation}
C^{N,\text{torus}}_{m,n} = R^{\text{torus}}_{m,n}\left(\delta_{N-mn,0} + \sum_{m'n'\leq N-mn} \frac{C^{N-mn}_{m',n'}}{\delta_{(m,-n)}-\delta_{(m',n')}} \right)
\end{equation}


An expression for the \(R_{m,n}^{\text{torus}}\) can be found on [[https://en.wikipedia.org/wiki/Virasoro_conformal_block][this wikipedia article]]. It can be rewritten

\[
R_{m,n}^{\text{torus}} = \frac{1}{2 D_{m,n}} \prod_{r\overset2=1-2m}^{2m-1} \prod_{s\overset2=1-2n}^{2n-1} \sqrt{\delta_{(r,s)} - \delta_1}
\]

where we do not actually take square roots, because each factor appears twice. The normalization factor is the same \(D_{m,n}\) as in the [[Dmn][four point]] case.

*** The =OnePointCorrelationFunctions= module

The module =OnePointCorrelationFunctions= defines

- a struct =OnePointCorrelation= that represents a one point function \[
  < V >,
  \]
- a method =computeCNmn= that computes the coefficients \(C^{N,\text{torus}}_{m,n}\) which serve to compute the conformal blocks that enter the expansion of the 1-pt function.

**** Header

#+begin_src julia
module OnePointCorrelationFunctions

export OnePointCorrelation, computeCNmn

using ..CFTData
import ..FourPointCorrelationFunctions: Dmn, δrs # re-use the Dmn from four-point functions
#+end_src

**** One-point function type

#+begin_src julia
struct OnePointCorrelation{T}
    charge::CentralCharge{T}
    field::Field{T}
end

"""Display a one-point function"""
function Base.show(io::IO, corr::OnePointCorrelation)
    println("One-point correlation function: < V > where ")
    print("V = "); show(corr.field)
end
#+end_src

**** Compute $C^{N,\text{torus}}_{m,n}$

The computation of the $C^{N,\text{torus}}_{m,n}$ is very similar to that of the [[*Compute $C^N_{m,n}$][coefficients $C^{N}_{m,n}$]]. We re-use much of the code.

#+begin_src julia
"""Order of a pole of Rmn^torus, assuming the central charge is generic"""
function Rmn_zero_order(m, n, corr::OnePointCorrelation)
    B = corr.charge["B"]
    V = corr.field
    if V.isKac && V.r%2==1 && V.s%2==1 && abs(V.r) <= 2*m-1 && abs(V.s) <= 2*n-1
        return 1
    end
    return 0
end

"""
Compute `Rmn^torus`.
lr indicates the left or right moving parts of the fields
TODO: value of regularisation
"""
function Rmn(m, n, corr::OnePointCorrelation, lr)
    B = corr.charge["B"]
    V = corr.field
    δ1 = V["δ"][lr]
    if Rmn_zero_order(m, n, corr) > 0
        return 0
    else
        res = prod(prod(δrs(r, s, B) - δ1 for r in 1:2:2*m-1) for s in 1-2n:2:2n-1)
        return res/(2*Dmn(m, n, B))
    end
end

function computeCNmn(N, m, n, corr::OnePointCorrelation, lr)
    B = corr.charge["B"]
    if Rmn_zero_order(m, n, corr) > 0
        return 0
    elseif m*n > N
        return 0
    elseif m*n == N
        return Rmn(m, n, corr, lr)
    else
        res = sum(sum(computeCNmn(N-m*n, mp, np, corr, lr)/(δrs(m, -n, B)-δrs(mp, np, B))
                      for mp in 1:N-m*n if mp*np <= N-m*n)
                  for np in 1:N-m*n)
        return Rmn(m, n, corr, lr) * ((N-m*n==0)+res)
    end
end
#+end_src

#+RESULTS:

**** End module

#+begin_src julia
end # end module
#+end_src

* Virasoro conformal blocks
:PROPERTIES:
:header-args:julia: :tangle ./src/ConformalBlocks.jl
:END:

The file [[./src/ConformalBlocks.jl][ConformalBlocks.jl]] implements Zamolodchikov's recursion formula for computing four-point conformal blocks on the sphere and one-point conformal blocks on the torus.

It uses types from [[*The ~CFTData~ module]] and [[*Correlation Functions]].

** Four-point conformal blocks on the sphere

*** Four-point conformal blocks

Conformal blocks encode the universal part of correlation functions. More precisely, by performing the OPE of \(V_{1}\) with \(V_{2}\) and \(V_{3}\) with \(V_{4}\) (s-channel), or \(1\leftrightarrow 4, 2\leftrightarrow 3\) (t-channel) or \(1\leftrightarrow 3, 2\leftrightarrow 4\) (u-channel), \(\mathcal F(x)\) can be written as

\[ \mathcal F(x) = \sum_{s \in \mathcal S^{(s)}} \frac{C_{12s}C_{s34}}{B_s} \mathcal F_\Delta^{(s)}(\Delta_i | x)\]

where \(\sum_{\mathcal S^{(s)}}\) is a sum over a basis of the \(s\)-channel spectrum \(\mathcal S^{(s)}\), \(B_s\) are two-point structure constants and \(C_{ijk}\) are three-point structure constants. Analogous expressions can be written for the \(t\)- and \(u\)- channels. The functions \(\mathcal F_{\Delta}^{(s)}(\Delta_i | x)\) are called \(s\)-channel conformal blocks (resp. t,u).

**** Normalisation of the blocks

Conformal blocks factorise into holomorphic and anti-holomorphic parts, and they are characterized by the normalisation conditions

#+name: block-normalisation
\begin{align}
 \mathcal{F}^{(s)}_\Delta(x) & \underset{x\to 0}{=} \left| x^{\Delta-\Delta_1-\Delta_2}\right|^2 \left(1+O(x)\right) \nonumber
 \\
 \mathcal{F}^{(t)}_\Delta(x) & \underset{x\to 1}{=} \left|(1-x)^{\Delta-\Delta_1-\Delta_4}\right|^2 \left(1+O(1-x)\right) \nonumber
 \\
 \mathcal{F}^{(u)}_\Delta(x) & \underset{x\to \infty}{=} \left|\left(\frac{1}{x}\right)^{\Delta+\Delta_1-\Delta_3} \right|^2\left(1+O\left(\frac{1}{x}\right)\right)
\end{align}

(we omit the \(\Delta_i\) dependence in the notation \(\mathcal{F}^{(u)}_\Delta(x)\)).

**** Notations

\[
x = \frac{\theta_2(q)^4}{\theta_3(q)^4}, \quad q = e^{-\pi\frac{K(1-x)}{ K(x)}}
\]

where

\[
\theta_3(q) = \sum_{n\in\mathbb{Z}} q^{n^2} \quad , \quad \theta_2(q) = 2q^\frac14\sum_{n=0}^\infty q^{n(n+1)}
\] 

are Jacobi special \(\theta\)-functions, and \(K(x)\) is the elliptic \(K\) function.

**** Expression

Our chiral \(s\)-channel conformal block is 

\[
\mathcal{F}^{(s)}_{\delta}(x) =  x^{E_0} (1-x)^{E_1} \theta_3(q)^{-4E_2}
(16q)^{\delta} H(q,\delta)
\] 

where we use the exponents
\[
E_0 = -\delta_1-\delta_2-\frac{c-1}{24} \quad , \quad E_1 = -\delta_1-\delta_4-\frac{c-1}{24} \quad ,
\quad E_2 = \delta_1+\delta_2+\delta_3+\delta_4+\frac{c-1}{24}
\] 
The non-trivial coefficient is the series

\[
H(q,\delta) = 1 + \sum_{N=1}^{N_{max}} \sum_{mn\leq N} C_{m,n}^N \frac{(16q)^N}{\delta-\delta_{(m,n)}}
\]

# If $R_{m,n}=0$, we compute a finite regularization of $R_{m,n}$. This is never used for computing chiral conformal blocks, but only for computing non-chiral logarithmic blocks. The regularization of vanishing factors is
# $$
# \left(\delta_2-\delta_1\right)_\text{reg} = 2p_2
# $$
# $$
# \left((\delta_2-\delta_1)^2 -2\delta_{(r,s)}(\delta_1+\delta_2) + \delta_{(r,s)}^2\right)_\text{reg} = 8p_1p_2p_{(r,s)}
# $$

**** Changing channels

In practice, we compute $s$-channel blocks, and then use these to compute $t$- and $u$- channel blocks, from the relations

#+name: tu-from-s
\begin{align*}
\mathcal{F}^{(t)}_{\Delta}(\Delta_1,\Delta_2,\Delta_3,\Delta_4|x)
&= (-1)^{S_1+S_2+S_3+S_4}
\mathcal{F}^{(s)}_{\Delta}(\Delta_1,\Delta_4,\Delta_3,\Delta_2|1-x)
\\
\mathcal{F}^{(u)}_\Delta(\Delta_1,\Delta_2,\Delta_3,\Delta_4|x)
&= (-1)^{S_1+S_2+S_3+S_4}
\left|x^{-2\Delta_1}\right|^2 \mathcal{F}^{(s)}_\Delta(\Delta_1,\Delta_3,\Delta_2,\Delta_4|\frac{1}{x})
\end{align*}

where \(S=\Delta-\bar\Delta\) is the conformal spin, which we assume to be integer. These relations are obtained from [[block-normalisation][the block normalisation]].

**** TODO Degenerate blocks

In the case \(\delta_1 = \delta_{(r_1,s_1)}\) and \(\delta_2 = \delta_{(r_2,s_2)}\) with \(r_i,s_i\in\mathbb{N}^*\), we have

\[
\left\{\begin{array}{l} m\in |r_1-r_2|+1+2\mathbb{N}
\\ n \in |s_1-s_2| + 1+2\mathbb{N} \end{array} \right. \quad
\Rightarrow \quad R_{m,n} = 0
\] 

and similarly if the fields with numbers \(3\) and \(4\) are degenerate. Thanks to \(R_{m,n}=0\) thus \(C_{m,n}^N=0\), the block \(\mathcal{F}^{(s)}_{\delta_{(m,n)}}(x)\) is finite and can be computed exactly.

This can be generalized to fractional indices \(r_i,s_i\). In this case, we have to add the following restriction, which was redundant for positive integer indices:
 
\[
\left\{\begin{array}{l} m\in |r_1+r_2|+1+2\mathbb{N}
\\ n \in |s_1+s_2| + 1+2\mathbb{N} \end{array} \right. \quad
\Rightarrow \quad R_{m,n} = 0
\] 

and similarly if the fields with numbers \(3\) and \(4\) are degenerate. In particular, for \(\delta_i = \delta_{(0,\frac12)}\), we have \(m\in 2\mathbb{N}+1\Rightarrow R_{m,n}=0\).

**** TODO Derivative and regularization

For the purpose of computing conformal blocks for logarithmic channel representations, we need to compute derivatives of conformal blocks with respect to the channel dimension, and regularized values of blocks at their poles. Taking the derivative amounts to \[
H(q, \delta) \to \log(16q) H(q,\delta) +\frac{\partial}{\partial\delta} H(q, \delta)
\] And the regularization we are interested in amounts to \[
\left.\frac{1}{\delta-\delta_{(m,n)}}\right|_{\delta=\delta_{(m,n)}} = \log(16q)-\frac{1}{4\delta_{(m,n)}}
\] The code can formally compute a regularization of the block's derivative, but this regularization is a priori not meaningful.

*** The ~FourPointBlocksSphere~ module

The module ~FourPointBlocksSphere~ exports

- a struct ~FourPointBlockSphere~ that encapsulates the data needed to compute a 4pt conformal block, namely a channel, four external fields and the field propagating in the channel
- a function ~F_four_point_sphere(block, charge, x)~ which computes the value of the non-chiral block \(\mathcal F_{\Delta}^{(s)}(\Delta_i | x)\) as defined in this paragraph.

**** Header

#+begin_src julia
#===========================================================================================

ConformalBlocks.jl contains modules that compute series expansions for
Virasoro four-point conformal blocks on the sphere and Virasoro one-point conformal blocks
on the torus.

Written by Paul Roux, adapting a Python code written by Sylvain Ribault & Rongvoram
Nivesvivat

===========================================================================================#


"""
Series expansion of four-point blocks on the sphere.
"""
module FourPointBlocksSphere

export FourPointBlockSphere, F_four_point_sphere

using ..CFTData, ..FourPointCorrelationFunctions
import ..FourPointCorrelationFunctions: permute_ext_fields
using Match, EllipticFunctions, Memoization

#+end_src

**** Four-point block sphere type

#+begin_src julia
#===========================================================================================
Struct FourPointBlockSphere
===========================================================================================#
"""
    FourPointBlockSphere{T}

Composite type that represents the list of arguments of a four-point conformal block:
a channel and a field propagating in the channel. The external fields and central charge are
provided in a `FourPointCorrelation` object.

# Example

```julia-repl
julia> c = CentralCharge("c",0.5); V = Field(c, "δ", 0.6, diagonal = true);
julia> FourPointBlockSphere("s", V)
Four-point block
Channel:        s
Channel Field:
Diagonal field of dimension:
  Δ = 0.5791666666666667 + 0.0im
  P = 0.0 + 0.7745966692414834im
  δ = 0.6000000000000001 + 0.0im
  p = 0.7745966692414834 + 0.0im
```
"""
struct FourPointBlockSphere{T}

    channel::String
    channelField::Field{T}

end

"""Display blocks"""
function Base.show(io::IO, block::FourPointBlockSphere)
    println("Four-point block")
    println("Channel:\t$(block.channel)")
    println("Channel Field:")
    show(block.channelField)
    # println("External Fields:")
    # print("1. "); show(block.extFields[1])
    # print("2. "); show(block.extFields[2])
    # print("3. "); show(block.extFields[3])
    # print("4. "); show(block.extFields[4])
end

# explicit names for the indices of left and right dimensions
const left = 1
const right = 2
#+end_src

**** Change of channel

The $t$- and $u$-channel blocks are computed from the $s$-channel one, using [[tu-from-s][the relation]] described above.

#+begin_src julia
#===========================================================================================
Get t- and u- channel blocks from s-channel block
===========================================================================================#
"""Prefactor to get t- or u-channel blocks from the s-channel block"""
function channelprefactor(block::FourPointBlockSphere, corr::FourPointCorrelation, x)
    @match block.channel begin
        "s" => 1
        "t" => (-1)^(sum(spin(corr.fields)))
        "u" => (-1)^(sum(spin.(corr.fields)))*abs2(x)^(-2*corr.fields[1]["Δ"])
    end
end

"""Cross-ratio at which to evaluate the s-channel block to get t- or u-channel block"""
function crossratio(channel, x)
    @match channel begin
        "s" => x
        "t" => 1-x
        "u" => 1/x
    end
end
#+end_src

**** Prefactors, elliptic nome

The nome $q$ is related to $x$ via

\[
q(x) = \exp(-\pi \frac{K(1-x)}{K(x)})
\]

where $K$ is the elliptic $K$ function. The inverse of this relation is

\[
x(q) = \left(\frac{\theta_{4}(q)}{\theta_{3}(q)}\right)^{2}
\]

#+begin_src julia
#===========================================================================================
Set prefactors, relate the cross-ratio x and the elliptic nome q
===========================================================================================#
"""Nome `q` from the cross-ratio `x`"""
@memoize qfromx(x) = exp(-π*ellipticK(1-x) / ellipticK(x))

""""Cross ratio `x` from the nome `q`"""
xfromq(q) = jtheta2(0,q)^4 / jtheta3(0,q)^4

"""Prefactor for getting the block F from H. The argument `lr` indicates if we are working
with a left or right moving block"""
function blockprefactor(block::FourPointBlockSphere, corr::FourPointCorrelation, x, lr)

    c = corr.charge["c"]
    e0 = - corr.fields[1]["δ"][lr] - corr.fields[2]["δ"][lr] - (c-1)/24
    e1 = - corr.fields[1]["δ"][lr] - corr.fields[4]["δ"][lr] - (c-1)/24
    e2 = sum(corr.fields[i]["δ"][lr] for i in 1:4) + (c-1)/24
    q=qfromx(x)

    return x^e0 * (1-x)^e1 * jtheta3(0,q)^(-4*e2) * (16*q)^block.channelField["δ"][1]
end

"""Degenerate dimensions"""
δrs(r, s, B) = -1/4 * (B*r^2 + 2*r*s + s^2/B)
#+end_src

**** Computation of the block

#+begin_src julia
#===========================================================================================
Compute the conformal block
===========================================================================================#
"""Compute the function ``H(q,δ)``."""
function H(q, Nmax, block::FourPointBlockSphere, corr::FourPointCorrelation, lr)
    δ = block.channelField["δ"][lr]
    B = corr.charge["B"]
    sq = 16*q
    res=1
    pow = 1
    for N in 1:Nmax
        sum_mn = sum(sum(computeCNmn(N, m, n, corr, block.channel, lr)/(δ-δrs(m, n, B))
                         for n in 1:N if m*n <= N) for m in 1:N)
        pow *= sq
        res += pow * sum_mn
    end
    return res
end

"""
    Fs_chiral(block::FourPointBlockSphere, corr::FourPointCorrelation, x, lr)

Compute the chiral conformal block

``\\mathcal F^{(s)}_{\\delta}(x)``

"""
function Fs_chiral(x, Nmax, block::FourPointBlockSphere, corr::FourPointCorrelation, lr)
    blockprefactor(block, corr, x, lr) * H(qfromx(x), Nmax, block, corr, lr)
end

"""Compute the chiral conformal block

``\\mathcal F^{(\\text{chan})}_{\\delta}(x)``

where `chan` is `s`, `t`, or `u`."""
function F_chiral(x, Nmax, block::FourPointBlockSphere, corr::FourPointCorrelation, lr)
    chan = block.channel
    Fs_chiral(crossratio(chan, x), Nmax, block, permute_ext_fields(corr, chan), lr)
end

"""
Compute the non-chiral conformal block

``\\mathcal F^{(\\text{chan})}_{\\delta}(x) \\overline{\\mathcal F}^{(\\text{chan})}_{\\delta}( \bar x )``

where `chan` is `s`,`t` or `u`.

TODO: logarithmic blocks
"""
function F_four_point_sphere(x, Nmax, block::FourPointBlockSphere, corr::FourPointCorrelation)
    channelprefactor(block, corr, x) * \
        F_chiral(x, Nmax, block::FourPointBlockSphere, corr::FourPointCorrelation, left) * \
        conj(F_chiral(conj(x), Nmax, block::FourPointBlockSphere, corr::FourPointCorrelation, right))
end
#+end_src

**** End of module

#+begin_src julia
end # end module
#+end_src

** One-point conformal blocks on the torus

*** One point blocks on the torus

A chiral one-point conformal block on the torus can be written as

#+name: F-chiral-torus
\[
\mathcal{F}^\text{torus}_{\Delta}(\Delta_1|q) = q^{\delta}\eta(q)^{-1}H^\text{torus}_{\Delta}(\Delta_1|q)
\]

where \(\eta\) is the Dedekind eta function, \(q\) is related to the torus' modular parameter \(\tau\) via \(q=e^{2i\pi \tau}\), and the non-trivial factor \(H^{\text{torus}}_{\Delta}(\Delta_1|q)\) can be computed by the following recursion formula:

\[
H^\text{torus}_{\Delta}(\Delta_1|q) = 1 + \sum_{N=1}^{N_{max}} \sum_{mn\leq N} C_{m,n}^{N,\text{torus}} \frac{q^N}{\delta-\delta_{(m,n)}}
\]

where the coefficient $C^{N,\text{torus}}_{m,n}$ is given [[CNmn-torus][here]].

*** The ~OnePointBlocksTorus~ module

The module ~OnePointBlocksTorus~ exports

- a struct ~OnePointBlockTorus~ that encapsulates the data needed to compute a 4pt conformal block, namely an external field.
- a function ~F_one_point_torus(block, charge, x)~ which computes the value of the non-chiral block \(\mathcal F_{\Delta}^{\text{torus}}(\Delta | q(x))\) as defined in [[F-chiral-torus][this equation]].

**** Header

#+begin_src julia
"""
Series expansion of one-point blocks on the torus
"""
module OnePointBlocksTorus

using ..CFTData, ..OnePointCorrelationFunctions
import EllipticFunctions: etaDedekind as η

export OnePointBlockTorus, F_one_point_torus

#===========================================================================================
Struct containing the data required to compute a block: an external field
===========================================================================================#
struct OnePointBlockTorus{T}
    channelField::Field{T}
end

# explicit names for the indices of left and right dimensions
const left = 1
const right = 2
#+end_src

**** Computation of the block

#+begin_src julia

qfromtau(τ) = exp(2im*π*τ)
δrs(r, s, B) = -1/4 * (B*r^2 + 2*r*s + s^2/B)

#===========================================================================================
Compute the conformal block
===========================================================================================#
"""
    H(q, Nmax, block, corr, leftright)
Compute the function  ``H^{\\text{torus}}(q,δ)``."""
function H(q, Nmax, block::OnePointBlockTorus, corr::OnePointCorrelation, lr)
    δ = block.channelField["δ"][lr]
    B = corr.charge["B"]
    res = 1
    pow = 1
    for N in 1:Nmax
        sum_mn = sum(sum(computeCNmn(N, m, n, corr, lr)/(δ-δrs(m, n, B))
                         for n in 1:N if m*n <= N) for m in 1:N)
        pow *= q
        res += pow * sum_mn
    end
    return res
end

"""
    Fs_chiral(block::FourPointBlockSphere, corr::FourPointCorrelation, x, lr)

Compute the chiral conformal block

``\\mathcal F^{\text{torus}}_{\\delta}(x)``

"""
function F_chiral(τ, Nmax, block::OnePointBlockTorus, corr::OnePointCorrelation, lr)
    δ = block.channelField["δ"][lr]
    return q^δ/η(τ) * H(qfromtau(τ), Nmax, block, corr, lr)
end

"""
Compute the non-chiral conformal block

`` \\mathcal F_{\\Delta}^{(\\text{chan})}(\\Delta_i| x)``

where ``\\text{chan}`` is `s`,`t` or `u`.

TODO: logarithmic blocks
"""
function F_one_point_torus(τ, Nmax, block::OnePointBlockTorus, corr::OnePointCorrelation)
    F_chiral(τ, Nmax, block, corr, left) * conj(F_chiral(conj(τ), Nmax, block, corr, right))
end
#+end_src

**** End of module

#+begin_src julia
end # end module
#+end_src

* Special functions
:PROPERTIES:
:header-args:julia: :tangle ./src/SpecialFunctions.jl
:END:

#+begin_src julia
#==================

SpecialFunctions.jl computes the special functions relevant for our applications in 2D CFT.

==================#

using EllipticFunctions

function log_double_gamma(beta, w)
end

function double_gamma(beta, w)
    return exp(log_double_gamma(beta, w))
end
#+end_src

* Unit testing
:PROPERTIES:
:header-args:julia: :tangle ./test/runtests.jl
:END:

#+begin_src julia
using JuliVirBootstrap
using Test

@testset "CFTData.jl" begin

    #ensure the relation between b and β does not change
    c1 = CentralCharge("c", -1.1+.2im)
    b = c1["b"]
    c2 = CentralCharge("b", b)
    @test c1["c"] == c2["c"]
    @test c1["β"] == c2["β"]

    #ensure the relation between p and P does not change
    left = 1
    right = 2
    V1 = Field(c1, "P", 0.5, diagonal=true)
    p = V1["p"][left]
    V2 = Field(c1, "p", p, diagonal=true)
    @test V1["P"] == V2["P"]

    #ensure the keyword diagonal also works for fields given from Kac indices
    V1 = Field(c1, Kac=true, r=3, s=4, diagonal=true)
    @test V1["Δ"][left] == V1["Δ"][right]


    #ensure degenerate and diagonal work well together
    V1 = Field(c1, Kac=true, degenerate=true, r=2, s=5, diagonal=true)
    @test V1["Δ"][left] == V1["Δ"][right]

end

@testset "FourPointCorrelationFunctions" begin

    left=1
    right=2

    c = CentralCharge("β", 1.2+.1*1im)
    V1 = Field(c, "Δ", 0.23+.11im, diagonal=true)
    V2 = Field(c, "Δ", 3.43, diagonal=true)
    V3 = Field(c, "Δ", 0.13, diagonal=true)
    V4 = Field(c, "Δ", 1.3, diagonal=true)
    corr = FourPointCorrelation(c, V1, V2, V3, V4)

    @test isapprox(JuliVirBootstrap.FourPointCorrelationFunctions.Rmn(2, 1, corr, "s", left),
                   0.31097697185245077-0.70523695127635733im, # value taken from Sylvain's code
                   atol=1e-8)

    @test isapprox(JuliVirBootstrap.FourPointCorrelationFunctions.Rmn(3, 3, corr, "t", left),
                   4.3964194233662846e-5-1.1534661157146291e-5im, # value taken from Sylvain's code
                   atol=1e-8)

    @test isapprox(JuliVirBootstrap.FourPointCorrelationFunctions.computeCNmn(7, 2, 3, corr, "s", left),
                   0.0019498393368877166+0.0026353877950837049im, # value taken from Sylvain's code
                   atol=1e-8)

end
#+end_src

* Development tests
:PROPERTIES:
:header-args:julia: :tangle ./test/devtests.jl :session test
:END:

#+begin_src julia :results silent
import Pkg; Pkg.activate(".")
using JuliVirBootstrap
#+end_src

#+begin_src julia
using JuliVirBootstrap, BenchmarkTools, EllipticFunctions

left=1;
right=2;

c = CentralCharge("β", big(1.2+.1*1im));
V1 = Field(c, "Δ", 0.23+.11im, diagonal=true);
V2 = Field(c, "Δ", 3.43, diagonal=true);
V3 = Field(c, "Δ", 0.13, diagonal=true);
V4 = Field(c, "Δ", 1.3, diagonal=true);
V = Field(c, "Δ", 0.1, diagonal = true);

x = BigFloat("0.05", RoundUp);
function test()
    corr = FourPointCorrelation(c, V1, V2, V3, V4)
    block = FourPointBlockSphere("s", V)
    calc = JuliVirBootstrap.FourPointBlocksSphere.Fs_chiral(x, 5, block, corr, left);
end;
#+end_src

#+RESULTS:

#+begin_src julia :results output
@btime test()
#+end_src

#+RESULTS:
:   1.046 ms (14047 allocations: 779.80 KiB)
: 2337.403601860713925958391009719678809055546365124340617423994215164154638728732 + 4771.392251761078219452023733430973991671048305871845481326523306698785531079742im


** Relation between four-point blocks on the sphere and one-point blocks on the torus

Four point blocks on the sphere are related to one-point blocks on the torus through the relation

\[
\mathcal H^{\text{torus}}_{c, P}(P_{1} | q^{2}) = \mathcal H_{c', \sqrt{2}P'}\left(\left. P'_{(0,\frac12)}, \left(\frac{P_{1}}{\sqrt{2}}\right)', P'_{(0,\frac12)}, P'_{(0,\frac12)} \right| q \right)
\]

where
+ $c'$ is related to $c$ via $\beta'=\frac\beta{\sqrt 2}$.
+ $P'$ denotes the Virasoro module with primary field of dimension $\Delta'(P') = \frac{c'-1}{24} - P'^{2}$

#+begin_src julia :results silent
import Pkg; Pkg.activate(".")
using JuliVirBootstrap, BenchmarkTools, EllipticFunctions

left=1;
right=2;

import JuliVirBootstrap.FourPointBlocksSphere.qfromx
c_torus = CentralCharge("b", 1.2+.1*1im);
c_sphere = CentralCharge("b", (1.2+.1*1im)/sqrt(2))

q = JuliVirBootstrap.FourPointBlocksSphere.qfromx(0.05)

P = 0.23+.11im
P1 = 0.41+1.03im
V_torus_chan = Field(c_torus, "P", P, diagonal=true)
δ_torus = V_torus_chan["δ"][left]
δ11_torus = Field(c_torus, Kac=true, r=1, s=1, diagonal=true)["δ"][left]
V_torus_ext = Field(c_torus, "P", P1, diagonal=true)
corr_torus = OnePointCorrelation

V_sphere_chan = Field(c_sphere, "P", sqrt(2)*P, diagonal=true)
δ_sphere = V_sphere_chan["δ"][left]
δ21_sphere = Field(c_sphere, Kac=true, r=2, s=1, diagonal=true)["δ"][left]
δ12_sphere = Field(c_sphere, Kac=true, r=1, s=2, diagonal=true)["δ"][left]
V_sphere_ext = Field(c_sphere, "P", P1/sqrt(2), diagonal=true)
VKac_sphere = Field(c_sphere, Kac=true, r=0, s=1//2, diagonal=true)

corr_torus = OnePointCorrelation(c_torus, V_torus_ext)
block_torus = OnePointBlockTorus(V_torus_chan)

corr_sphere = FourPointCorrelation(c_sphere, [VKac_sphere, V_sphere_ext, VKac_sphere,VKac_sphere])
block_sphere = FourPointBlockSphere("s", V_sphere_chan)


C111_torus = JuliVirBootstrap.OnePointCorrelationFunctions.computeCNmn(1, 1, 1, corr_torus, left)

C212_sphere = JuliVirBootstrap.FourPointCorrelationFunctions.computeCNmn(2, 1, 2, corr_sphere, "s", left)
C221_sphere = JuliVirBootstrap.FourPointCorrelationFunctions.computeCNmn(2, 2, 1, corr_sphere, "s", left)

16^2 * (C212_sphere/(δ_sphere-δ12_sphere) + C221_sphere/(δ_sphere-δ21_sphere))
C111_torus/(δ_torus-δ11_torus)

h1 = JuliVirBootstrap.OnePointBlocksTorus.H(q^2, 5, block_torus, corr_torus, left)
h2 = JuliVirBootstrap.FourPointBlocksSphere.H(q, 5, block_sphere, corr_sphere, left)
#+end_src

#+begin_src julia :results output
println("torus block = $h1 \nsphere block = $h2")
#+end_src

#+RESULTS:
: torus block = 1.0000059915273005 - 1.1912765043504052e-5im
: sphere block = 1.000005991527301 - 1.1912765042311957e-5im
