"""
        LeftRight{T}
Left and right pairs of objects. Can be accessed with
`obj[:left]` and `obj[:right]`.
"""
const LeftRight{T} = Tuple{T,T} # left and right pairs of things

function Base.getindex(x::LeftRight, s::Symbol)
    s === :left && return x[1]
    s === :right && return x[2]
    error("tried to access pair element other than 1, 2, :left or :right")
end

include("CFTData/central_charges.jl")
include("CFTData/conformal_dimensions.jl")
include("CFTData/fields.jl")

# convenience aliases
const CC = CentralCharge
const CD = ConformalDimension
