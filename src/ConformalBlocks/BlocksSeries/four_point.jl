struct FourPointBlock{T} <: Block{T}

    Nmax::Int
    channel::Symbol
    channel_field::Field{T}
    _CNmn::CNmnTable{T}
    _coefficients::Tuple{Vector{T}, Vector{T}}

end

function FourPointBlock(
    Nmax::Int,
    c::CentralCharge,
    corr::FourPointCorrelation,
    channel::Symbol,
    channel_field::Field{T};
    der=false,
    reg=false
) where {T}
    CNmn = corr._CNmn[channel]
    tmp = FourPointBlock{T}(
        Nmax, channel, channel_field, CNmn, (Vector{T}(), Vector{T}())
    )
    coeffs = Tuple(series_H(c, tmp, lr, der, reg) for lr in (left, right))

    FourPointBlock{T}(
        Nmax, channel, channel_field, CNmn, coeffs
    )
end