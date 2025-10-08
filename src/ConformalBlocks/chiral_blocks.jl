struct ChiralBlock{T} <: Block{T}
    corr::CCo{T}
    chan_dim::CD{T}
    fields::Tuple{Vararg{CD{T}}}
    coeffs::Vector{T}
    coeffs_der::Vector{T}
    missing_terms::Vector{T}
    chan::Symbol
    Δmax::Int
end
const CBlock = ChiralBlock

function ChiralBlock(co::CCo{T}, chan::Symbol, d::CD, Δmax::Int, der) where {T}
    CNmn = getCNmn(co, chan)
    coeffs = series_H(d, Δmax, CNmn)
    coeffs_der = Vector{T}()
    if der
        coeffs_der = series_H_der(d, Δmax, CNmn)
    end
    if !d.degenerate
        missing_terms = []
    else
        r, s = indices(d)
        missing_terms =
            [(N, r, s) in CNmn.keys ? CNmn[N, r, s] : zero(T) for N = 0:Δmax]
    end

    CBlock{T}(co, d, co.fields, coeffs, coeffs_der, missing_terms, chan, Δmax)
end

getRmn(b::CBlock) = getRmn(b.corr, b.chan)
getRmnreg(b::CBlock) = getRmnreg(b.corr, b.chan)
getCNmn(b::CBlock) = getCNmn(b.corr, b.chan)
getc(b::CBlock) = b.corr.c

function series_H(d::CD{T}, Δmax, CNmn) where {T}
    B = d.c.B
    P = d.P
    P2 = P^2
    isKac = d.isKac
    r, s = d.r, d.s
    δKac = isKac ? δrs(r, s, B) : zero(T)

    H = zeros(T, Δmax + 1)
    H[1] = one(T)
    for (N, m, n) in CNmn.keys
        N > Δmax && continue
        δ = δrs(m, n, B)
        if isKac && m == r && n == s
            coeff = -inv(4δKac)
        else
            coeff = inv(P2 - δ)
        end
        H[N+1] += CNmn[N, m, n] * coeff
    end

    return H
end

function series_H_der(d::CD{T}, Δmax, CNmn) where {T}
    B = d.c.B
    P = d.P
    P2 = P^2

    H = zeros(T, Δmax + 1)
    for (N, m, n) in CNmn.keys
        N > Δmax && continue
        δ = δrs(m, n, B)
        coeff = -2P / (P2 - δ)^2
        H[N+1] += CNmn[N, m, n] * coeff
    end

    return H
end

function Base.show(io::IO, ::MIME"text/plain", b::CBlock)
    println(io, "Chiral block for the correlation")
    println(io, b.corr)
    println(io, "Channel: $(b.chan), $(b.chan_dim)")
end

function Base.show(io::IO, b::CBlock)
    println(io, "F^($(b.chan))($(b.chan_dim))")
end
