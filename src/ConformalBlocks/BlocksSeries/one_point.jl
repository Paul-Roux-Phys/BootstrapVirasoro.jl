struct OnePointBlock{T} <: Block{T}

    Nmax::Int
    channel_field::Field{T}
    _CNmn::CNmnTable{T}
    _coefficients::Tuple{Vector{T}, Vector{T}}

end

function OnePointBlock(
    Nmax::Int,
    c::CentralCharge,
    corr::OnePointCorrelation,
    channel_field::Field{T};
    der=false,
    reg=false
) where {T}
    CNmn = corr._CNmn
    tmp = OnePointBlock{T}(
        Nmax, channel_field, CNmn, (Vector{T}(), Vector{T}())
    )

    coeffs = Tuple(series_H(c, tmp, lr, der, reg) for lr in (left, right))

    OnePointBlock{T}(
        Nmax, channel, channel_field, CNmn, coeffs
    )
end