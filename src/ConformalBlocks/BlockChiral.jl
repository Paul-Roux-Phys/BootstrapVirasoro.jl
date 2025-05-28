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
    coeffs = series_H(d, Nmax, CNmn, false)
    coeffs_der = Vector{T}()
    if der
        coeffs_der = series_H(d, Nmax, CNmn, true)
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

function Base.show(io::IO, ::MIME"text/plain", b::BlockChiral)
    println(io, "Chiral block for the correlation")
    println(io, b.corr)
    println(io, "Channel: $(b.channel), $(b.channel_dimension)")
end

function Base.show(io::IO, b::BlockChiral)
    println(io, "F^($(b.channel))($(b.channel_dimension))")
end

function series_H_N(N, d::ConformalDimension, CNmn, der=false)
    B = d.c.B
    (N == 0 && !der) && return 1
    (N == 0 &&  der) && return 0
    P = d.P
    ind = [(m, n) for m in 1:N for n in 1:N if (N, m, n) in CNmn.keys]
    isempty(ind) && return 0
    function coeff(m, n)
        # for the derivative series
        der && return -2P / (P^2 - δrs(m, n, B))^2
        # for the regularised series
        d.isKac && (m, n) == (d.r, d.s) && return -inv(4*δrs(d.r, d.s, B))
        # vanilla series
        return inv(P^2 - δrs(m, n, B))
    end
    return sum(CNmn[N, m, n] * coeff(m, n) for (m, n) in ind)
end

function series_H(d::ConformalDimension, Nmax, CNmn, der=false)
    [series_H_N(N, d, CNmn, der) for N in 0:Nmax]
end
