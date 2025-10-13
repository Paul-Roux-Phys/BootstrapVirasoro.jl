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
        missing_terms = Vector{T}()
    else
        r, s = indices(d)
        missing_terms =
            [(N > 0 && (r, s) in CNmn.keys[N]) ? CNmn[N, r, s] : zero(T) for N = 0:Δmax]
    end

    CBlock{T}(co, d, co.fields, coeffs, coeffs_der, missing_terms, chan, Δmax)
end

getRmn(b::CBlock) = getRmn(b.corr, b.chan)
getRmnreg(b::CBlock) = getRmnreg(b.corr, b.chan)
getCNmn(b::CBlock) = getCNmn(b.corr, b.chan)
getc(b::CBlock) = b.corr.c

# function series_H(d::CD{T}, Δmax, CNmn) where {T}
#     B = d.c.B
#     P = d.P
#     P2 = P^2
#     isKac = d.isKac
#     r, s = d.r, d.s
#     δKac = isKac ? δrs(r, s, B) : zero(T)
#     p = precision(real(T))
#     ϵ = big"2."^(-2p)
#     δs = [δrs(m, n, B) for m in 0:Δmax+1, n in 0:Δmax+1]
#     coeffs = zeros(T, Δmax+1, Δmax+1)
#     mns = Set((k[1], k[2]) for k in CNmn.keys)
#     for m in 0:Δmax+1, n in 0:Δmax+1
#         if (m, n) in mns
#             coeffs[m,n] = 1/(P2-δrs(m, n, B))
#         end
#     end
#     isKac ? coeffs[r, s] = -1/(4δrs(r, s, B)) : nothing

#     H = zeros(T, Δmax + 1)
#     H[1] = one(T)
#     for N in 1:Δmax+1
#         for (m, n) in mns
#             if (N, m, n) in CNmn.keys
#                 # if two coefficients in a row are below numerical precision,
#                 # the rest will be too. stop.
#                 (N > Δmax || (N > 2 && abs2(H[N]) < ϵ && abs2(H[N-1]) < ϵ)) ||
#                     !((N, m, n) in CNmn.keys) && continue
#                 δ = δs[m, n]
#                 if isKac && m == r && n == s
#                     coeff = -inv(4δKac)
#                 else
#                     coeff = inv(P2 - δ)
#                 end
#                 H[N+1] += CNmn[N, m, n] * coeffs[m, n]
#             end
#         end
#     end

#     return H
# end

function series_H(d::CD{T}, Δmax, CNmn) where {T}
    P = d.P
    P2 = P^2
    isKac = d.isKac
    r, s = d.r, d.s
    ϵ = big"2.0"^(-2precision(BigFloat))

    coeffs = Matrix{T}(undef, Δmax + 2, Δmax + 2)
    all_mns = union([CNmn.keys[N] for N in 1:Δmax]...)
    for (m, n) in all_mns
        if isKac && m == r && n == s
            coeffs[m, n] = -inv(4CNmn.δs[r, s])
        else
            coeffs[m, n] = inv(P2 - CNmn.δs[m, n])
        end
    end

    H = zeros(T, Δmax + 1)
    H[1] = one(T)
    for N in 1:Δmax
        for (m, n) in CNmn.keys[N]
            N > Δmax || m > Δmax || n > Δmax && continue
            H[N+1] += CNmn[N, m, n] * coeffs[m, n]
        end
    end

    return H
end

function series_H_der(d::CD{T}, Δmax, CNmn) where {T}
    P = d.P
    P2 = P^2

    H = zeros(T, Δmax + 1)
    coeffs = Dict((m, n) => -2P / (P2 - CNmn.δs[m, n])^2
              for (m, n) in union([CNmn.keys[N] for N in 1:Δmax]...))

    for N in 1:Δmax
        for (m, n) in CNmn.keys[N]
            N > Δmax && continue
            H[N+1] += CNmn[N, m, n] * coeffs[(m, n)]
        end
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
