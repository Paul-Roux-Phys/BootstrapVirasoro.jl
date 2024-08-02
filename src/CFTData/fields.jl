using ..JuliVirBootstrap: left, right

"""
    Field{T}
Object representing a conformal field.
Contains the conformal dimensions, and flags saying whether the field has (rational) Kac indices, is degenerate, or diagonal.
"""
struct Field{T <: Union{AbstractFloat, Complex{Float64}, Complex{BigFloat}}}

    dims::Tuple{ConformalDimension{T}, ConformalDimension{T}}
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
julia> Field(charge, :δ, 2, 3)
Non-diagonal field with (left, right) dimensions:
Δ = ( 2.1579142011834325 - 0.6789940828402367im, 3.1579142011834316 - 0.6789940828402367im )
P = ( 0.0 + 1.4142135623730951im, 0.0 + 1.7320508075688772im )
δ = ( 2.0000000000000004 + 0.0im, 2.9999999999999996 + 0.0im )
p = ( 1.4142135623730951 + 0.0im, 1.7320508075688772 + 0.0im )
```
```julia-repl
julia> charge = CentralCharge();
julia> Field(charge, :δ, 1, diagonal=true)
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
        @assert Kac == true """
          A non-diagonal field must be given from Kac indices.
          If you mean to define a diagonal field, use `diagonal=true`.
        """
        @assert r*s % 1 == 0 "The product r*s of Kac indices must be an integer"
        dim_right = ConformalDimension(c, sym, dim_left, Kac=Kac, r=r, s=-s)
    end

    Field{T}((dim_left, dim_right), diagonal, degenerate)
end

function Base.getproperty(V::Field, s::Symbol)
    ds = getfield(V, :dims)
    s === :P && return ds[left].P, ds[right].P
    s === :Δ && return ds[left].Δ, ds[right].Δ
    s === :p && return ds[left].p, ds[right].p
    s === :δ && return ds[left].δ, ds[right].δ
    s in (:r, :s) && return getfield(ds[left], s) # by convention V_(r,s) denotes the field
                                            # with left right dimension P_(r, s), P_(r, -s)
    s === :isKac && return (
        V.dims[left].isKac && V.dims[right].isKac && 
        V.dims[left].r == V.dims[left].r && V.dims[left].s == -V.dims[right].s
    )
    return getfield(V, s)
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
        @warn "You are computing the spin of a field not defined from Kac indices"
        return V.Δ[1] - V.Δ[2]
    end
end

function swap_lr(V::Field{T}) where {T}
    return Field{T}((V.dims[right], V.dims[left]), V.isdiagonal, V.isdegenerate)
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
        show(V.dims[left])
    else
        println(io, "Non-diagonal $(typeof(V))")
        print(io, "left: ")
        show(V.dims[left])
        print(io, "\nright: ")
        show(V.dims[right])
    end
end
