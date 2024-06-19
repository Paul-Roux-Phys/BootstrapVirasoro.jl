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

using Match, Latexify;

export CentralCharge
export ConformalDimension
export Field, spin

"""
    CentralCharge{T}

Type representing a central charge.
T is expected to be a real or complex number, of standard or arbitrary precision
"""
struct CentralCharge{T <: Union{AbstractFloat, Complex{Float64}, Complex{BigFloat}}}

    β::T

end

"""Get B from given parameter"""
function Bfrom(s::Symbol, x)
    a = (x-1)*(x-25)
    @match s begin
        :β => -x^2
        :c => if isreal(a) && a > 0
                 (x-13+sqrt((x-1)*(x-25)))/12
              else # a is complex
                 (x-13+sqrt(complex((x-1)*(x-25))))/12
              end
        :b => x^2
        :B => x
    end
end

"""Get asked parameter from B"""
function Bto(s::Symbol, x)
    rx = sqrt(complex(x))
    res = @match s begin
        :β => -im*rx
        :c => 13+6*x+6/x
        :b => -rx
        :B => x
    end
    return isreal(res) ? real(res) : res
end


function Base.getproperty(c::CentralCharge, s::Symbol)
    β = Bto(:β, Bfrom(:β, getfield(c, :β)))
    β = isreal(β) ? real(β) : β
    if s === :β
        β
    elseif s === :c
        13 - 6*β^2 - 6/β^2
    elseif s === :B
        -β^2
    elseif s === :b
        -im*β
    elseif s === :n
        -2*cos(oftype(β, π)*β^2)
    else
        error("$s is not a supported parametrisation of the central charge")
    end
end

"""
    CentralCharge(parameter, value)

Constructor function for the CentralCharge type.

Given one of the four parameters `c`, `b`, `β`, `B` and its value,
creates an object CentralCharge{T} where T is real if `β` is real.

# Example
```julia-repl
julia> setprecision(BigFloat, 20, base=10)
julia> CentralCharge(big"1.2")
c = 0.1933333333333333332741, β = 1.200000000000000000003

```
"""
function CentralCharge(s::Symbol, x)
    β = Bto(:β, Bfrom(s, x))
    CentralCharge(β)
end

"""Display an object of type CentralCharge"""
function Base.show(io::IO, c::CentralCharge)
    println("c = $(c.c), β = $(c.β)")
end

"""Get P from any given parameter"""
function Pfrom(s::Symbol, x, c::CentralCharge)
    res = @match s begin
        :Δ => sqrt(complex(x - (c.c-1)/24))
        :δ => sqrt(complex(x))
        :P => x
        :p => im*x
    end
    return isreal(res) ? real(res) : res
end

"""Get all parameters from P"""
function Pto(s::Symbol, x, c::CentralCharge)
    @match s begin
        :Δ => x^2 + (c.c-1)/24
        :δ => x^2
        :P => x
        :p => -im*x
        :w => -2*cos(oftype(c.β, π)*c.β*x)
    end
end

"""
    ConformalDimension{T}
Type for encoding a conformal dimension, and conveniently access its values in all parametrisations
"""
struct ConformalDimension{T <: Union{AbstractFloat, Complex{Float64}, Complex{BigFloat}}}

    c::CentralCharge{T}
    P::T
    isKac::Bool
    r::Rational
    s::Rational

end

function ConformalDimension(c::CentralCharge{T}, sym::Symbol=:P, P=0; Kac=false, r=0, s=0) where {T}
    if Kac
        P = (r*c.β-s/c.β)/2
    else
        P = Pto(:P, Pfrom(sym, P, c), c)
    end
    ConformalDimension{T}(c, P, Kac, r, s)
end

function Base.getproperty(d::ConformalDimension, s::Symbol)
    c = getfield(d, :c)
    P = Pto(:P, Pfrom(:P, getfield(d, :P), c), c)
    P = isreal(P) ? real(P) : P
    if s in (:P, :p, :Δ, :δ, :w)
        return Pto(s, P, c)
    else
        return getfield(d, s)
    end
    res = isreal(res) ? real(res) : res
end

const left, right = 1, 2


"""
    Field{T}
Object representing a conformal field.
Contains the conformal dimensions, and flags saying whether the field has (rational) Kac indices, is degenerate, or diagonal.
"""
struct Field{T <: Union{AbstractFloat, Complex{Float64}, Complex{BigFloat}}}

    dim::Tuple{ConformalDimension{T}, ConformalDimension{T}}
    isdiagonal::Bool
    isdegenerate::Bool

end

"""
   TODO: update the examples
    Field(charge, parameter, leftvalue, rightvalue; kwargs...)

Constructor function for the Field type.

Given a charge `charge`, one of the four parameters `Δ`, `δ`, `P`, `p` and two values,
create an object `Field{T}` (where T is the type of the values in `charge`) that represents a
field of left and right dimensions given by leftvalue and rightvalue in the chosen
parametrisation.
If given only one value for the parameters `Δ`, `δ`, `P` or `p`, the field is diagonal by default

# keyword arguments:

- `Kac::Bool`: if set to true, the field can be constructed from the values of its r and s
indices. By convention V_(r,s) has left and right momenta (P_(r,s), P_(r,-s))
- `r::Rational`,`s::Rational`: used in conjunction to `Kac=true`, must be given rational
values,
- `degenerate::Bool`: set to True if the field is degenerate,
- `diagonal::Bool`: set to True to get a diagonal field ; only the leftvalue needs to be
given.

# Examples
```julia-repl
julia> charge = CentralCharge(:b, big(0.5));
julia> field = Field(charge, Kac=true, r=0, s=1)
Non-diagonal field with Kac indices r = 0//1, s = 1//1 and (left,right) dimensions:
Δ = ( 2.5625 + 0.0im, 2.5625 + 0.0im )
  P = ( -0.0 - 1.0im, 0.0 + 1.0im )
δ = ( 1.0 - 0.0im, 1.0 + 0.0im )
p = ( -1.0 + 0.0im, 1.0 + 0.0im )
```
```julia-repl
julia> charge = CentralCharge(:β, 1.5+im);
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
function Field(c::CentralCharge{T}, sym::Symbol=:P, dim=0;
               Kac=false, r=0, s=0, degenerate=false, diagonal=false) where {T}

    if !Kac
        # diagonal = true # a field not given from Kac indices is diagonal
    end
    if degenerate # degenerate fields are diagonal and must be given from Kac indices
        Kac = true
        diagonal = true
    end
    dim_left = ConformalDimension(c, sym, dim, Kac=Kac, r=r, s=s)
    if diagonal
        dim_right = dim_left
    else
        @assert Kac==true "A non-diagonal field must be given from Kac indices"
        dim_right = ConformalDimension(c, sym, dim_left, Kac=Kac, r=r, s=-s)
    end

    Field{T}((dim_left, dim_right), diagonal, degenerate)
end

function Base.getproperty(V::Field, s::Symbol)
    ds = getfield(V, :dim)
    if s === :P
        return ds[left].P, ds[right].P
    elseif s === :Δ
        return ds[left].Δ, ds[right].Δ
    elseif s === :p
        return ds[left].p, ds[right].p
    elseif s === :δ
        return ds[left].δ, ds[right].δ
    elseif s in (:r, :s)
        return getfield(ds[left], s) # by convention V_(r,s) denotes the field with left right dimension P_(r, s), P_(r, -s)
    elseif s === :isKac
        return (V.dim[left].isKac && V.dim[right].isKac && V.dim[left].r == V.dim[left].r && V.dim[left].s == -V.dim[right].s)
    else
        return getfield(V, s)
    end
end

# Overload the == operator
function Base.:(==)(V1::Field, V2::Field)
    return V1.Δ == V2.Δ
end

"""Compute the spin Δleft - Δright of a field."""
function spin(V::Field)::Rational
    if V.isdiagonal
        return 0
    elseif V.isKac
        return V.r*V.s
    else # this should never happen
        return V.Δ[1] - V.Δ[2]
    end
end

function Base.show(io::IO, d::ConformalDimension)
    if d.isKac
        print(io, "Kac indices r = $(d.r), s=$(d.s)")
    else
        print(io, "Δ = $(d.Δ), P = $(d.P)")
    end
end

function Base.show(io::IO, V::Field)
    if V.isdiagonal
        print(io, "Diagonal $(typeof(V)) with ")
        show(V.dim[left])
    else
        println(io, "Non-diagonal $(typeof(V))")
        print(io, "left: ")
        show(V.dim[left])
        print(io, "\nright: ")
        show(V.dim[right])
    end
end

end # end module
