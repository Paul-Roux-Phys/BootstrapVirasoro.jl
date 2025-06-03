struct BlockChiral{T} <: Block{T}

    corr::CorrelationChiral{T}
    channel::Symbol
    channel_dimension::ConformalDimension{T}
    Nmax::Int
    _coefficients::Vector{T}
    _coefficients_der::Vector{T}
    _missing_terms::Vector{T}

end

function BlockChiral(
    corr::CorrelationChiral,
    chan::Symbol,
    d::ConformalDimension{T},
    Nmax::Int;
    der=false
) where {T}
    CNmn = corr._CNmn[chan]
    coeffs = series_H(d, Nmax, CNmn)
    coeffs_der = Vector{T}()
    if der
        coeffs_der = series_H_der(d, Nmax, CNmn)
    end
    if !isdegenerate(d)
        missing_terms = []
    else
        r, s = d.indices
        missing_terms = [
            (N, r, s) in CNmn.keys ? CNmn[N, r, s] : zero(T)
            for N in 0:Nmax
        ]
    end

    BlockChiral{T}(corr, chan, d, Nmax, coeffs, coeffs_der, missing_terms)
end

BlockChiral(corr::CorrelationChiral, chan, d::ConformalDimension; der=false) = 
    BlockChiral(corr, chan, d, corr.Nmax, der=der)

BlockChiral(corr::CorrelationChiral, chan, V::Field, lr::Symbol, Nmax::Int; der=false) = 
    BlockChiral(corr, chan, V.dims[lr], Nmax, der=der)

BlockChiral(corr::CorrelationNonChiral, chan, V::Field, lr, Nmax; der=false) = 
    BlockChiral(corr[lr], chan, V.dims[lr], Nmax, der=der)

BlockChiral(corr::CorrelationNonChiral, chan, d::ConformalDimension, lr, Nmax; der=false) = 
    BlockChiral(corr[lr], chan, d, Nmax, der=der)

BlockChiral(corr, chan, V_or_d, lr::Symbol; der=false) = 
    BlockChiral(corr, chan, V_or_d, lr, corr.Nmax, der=der)

function Base.getproperty(b::BlockChiral, s::Symbol)
    s in (:_Rmn, :_CNmn) && return getfield(b.corr, s)[b.channel]
    s in (:c, :dims) && return getproperty(b.corr, s)
    return getfield(b, s)
end

function series_H(d::ConformalDimension{T}, Nmax, CNmn) where {T}
    B = d.c.B
    P = d.P
    P2 = P^2
    isKac = d.isKac
    r, s = d.r, d.s
    δKac  = isKac ? δrs(r, s, B) : zero(T)

    H = zeros(T, Nmax+1)
    H[1] = one(T)
    for (N, m, n) in CNmn.keys
        N > Nmax && continue
        δ = δrs(m, n, B)
        if isKac && m == r &&  n == s
            coeff = -inv(4δKac)
        else
            coeff = inv(P2-δ)
        end
        H[N+1] += CNmn[N, m, n] * coeff
    end

    return H
end

function series_H_der(d::ConformalDimension{T}, Nmax, CNmn) where {T}
    B = d.c.B
    P = d.P
    P2 = P^2

    H = zeros(T, Nmax+1)
    for (N, m, n) in CNmn.keys
        N > Nmax && continue
        δ = δrs(m, n, B)
        coeff = -2P / (P2 - δ)^2
        H[N+1] += CNmn[N, m, n] * coeff
    end

    return H
end

function Base.show(io::IO, ::MIME"text/plain", b::BlockChiral)
    println(io, "Chiral block for the correlation")
    println(io, b.corr)
    println(io, "Channel: $(b.channel), $(b.channel_dimension)")
end

function Base.show(io::IO, b::BlockChiral)
    println(io, "F^($(b.channel))($(b.channel_dimension))")
end
