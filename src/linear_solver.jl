function StructureConstants(S::Dict{Symbol,U}) where {T,U<:ChannelSpectrum{T}}
    constants = Dict(chan => Dict{Field{T},T}() for chan in keys(S))
    errors = Dict(chan => Dict{Field{T},T}() for chan in keys(S))
    return StructureConstants{T}(constants, errors)
end

function Base.getproperty(c::StructureConstants, s::Symbol)
    s === :fields && begin
        consts = getfield(c, :constants)
        return vcat([[V for V in keys(consts[chan])] for chan in keys(consts)]...)
    end
    getfield(c, s)
end

function Base.show(io::IO, c::StructureConstants)
    for (chan, vals) in c.constants
        println(io, "Channel $chan:")
        for (field, constant) in vals
            println(io, "$field: $constant")
        end
    end
end

function format_complex(z::Complex{<:Real}; digits=5)
    real_str = @sprintf("%.*f", digits, real(z))
    imag_str = @sprintf("%.*f", digits, abs(imag(z)))
    sign = imag(z) < 0 ? "-" : "+"
    buf = real(z) > 0 ? " " : ""
    return "$buf$real_str $sign $(imag_str)im"
end

function Base.show(io::IO, c::StructureConstants)
    chans = keys(c.constants)
    nondiags = Dict(
        chan => sort(
            [V for V in keys(c[chan]) if !isdiagonal(V)],
            by=V -> (V.r, V.s)
        )
        for chan in chans
    )
    diags = Dict(
        chan => sort(
            [V for V in keys(c[chan]) if isdiagonal(V) && !isdegenerate(V)],
            by=V -> real(total_dimension(V))
        )
        for chan in chans
    )
    degs = Dict(
        chan => sort(
            [V for V in keys(c[chan]) if isdegenerate(V)],
            by=V -> real(total_dimension(V))
        )
        for chan in chans
    )

    for chan in chans
        # Collect all the labels for this channel
        all_Vs = vcat(degs[chan], diags[chan], nondiags[chan])

        # Find max width of label for alignment
        max_label_width = maximum(length(string(V)) for V in all_Vs)

        # Print channel header
        printstyled(io, "Channel $(chan)\n", color=:red, bold=true)
        printstyled(io, repeat('=', 9 + length(string(chan))) * "\n", color=:red, bold=true)

        # Print each row with alignment
        for V in all_Vs
            label = rpad(string(V), max_label_width)
            value = format_complex(c[chan][V])
            println(io, "$label | $value")
        end
    end
end

Base.length(c::StructureConstants) = sum(length(c.constants[chan]) for chan in keys(c.constants))

function solve!(s::BootstrapSystem; precision_factor=1)
    # raw vector format
    sol1, sol2 = setprecision(BigFloat, precision_factor * precision(BigFloat)) do
        (
            s.matrix.LHS[3:end, :] \ s.matrix.RHS[3:end],
            s.matrix.LHS[1:end-2, :] \ s.matrix.RHS[1:end-2]
        )
    end
    errors = @. abs((sol1 - sol2) / sol2)

    # back to dictionary format
    mat = s.matrix
    nb_fields = vcat([0], [length(mat.fields[chan]) for chan in mat.channels])
    for (i, chan) in enumerate(mat.channels)
        for (j, V) in enumerate(mat.fields[chan])
            s.consts[chan][V] = sol1[nb_fields[i] + j]
            s.consts.errors[chan][V] = errors[nb_fields[i] + j]
        end
    end
end




# """
#     nullspace_qr(A; tol=nothing)

# Compute a basis for the (approximate) nullspace of A using a rank-revealing QR with column pivoting.

# Returns Z (size n×k) whose columns span null(A), i.e. A*Z ≈ 0.

# If tol===nothing, it is chosen automatically based on eps and R’s diagonal.
# """
# function nullspace_qr(A::AbstractMatrix{T}; tol=nothing) where {T}
#     m, n = size(A)
#     F = eltype(A)

#     # 1) QR with column pivoting:  A[:,p] = Q*R
#     Q, R, p = qr(A, Val(true))

#     # 2) Extract only the first n diagonal entries of R
#     rdiag = abs.([R[i, i] for i in 1:min(m, n)])

#     # 3) Pick tolerance if none provided
#     tol === nothing && (
#         tol = maximum(rdiag) * eps(real(F)) * max(m, n)
#     )

#     # 4) Numerical rank r and nullity k
#     r = count(x -> x > tol, rdiag)
#     k = n - r
#     if k == 0
#         return k, rdiag, zeros(F, n, 0)
#     end

#     # 5) Partition R into R11 (r×r) and R12 (r×k)
#     R11 = R[1:r, 1:r]
#     R12 = R[1:r, r+1:n]

#     # 6) Build the nullspace in the *pivoted* coordinate:
#     #      [ x1 ]      [ -R11 \ R12 ]
#     #      [ x2 ]  with x2 = I_k
#     Z_piv = vcat(-R11 \ R12, I(k))

#     # 7) Undo the column pivot: rows of Z_piv correspond to permuted variables p
#     Z = similar(Z_piv, F, n, k)
#     for j in 1:n
#         Z[p[j], :] = Z_piv[j, :]
#     end

#     return k, rdiag, Z
# end

# function setup_matrix(B::BootstrapEquations)
#     # M = Matrix(undef, 
# end

# function setup_matrix(S::Dict{Symbol,U}; extrapoints::Int=6) where {T,U<:Spectrum{T}}
#     chans = [chan for (chan, _) in S]
#     fields = Dict(chan => collect(keys(block_values[chan])) for chan in chans)
#     M = Dict(
#         chan => hcat([block_values[chan][V] for V in fields[chan]]...)
#         for chan in chans
#     )
#     zero_block = Dict(
#         chan => zeros(T, size(M[chan])...)
#         for chan in chans
#     )

#     # Form the matrix
#     if length(chans) == 2
#         matrix = [M[chans[1]] (.-M[chans[2]])]
#     elseif length(chans) == 3
#         matrix = [M[chans[1]] (.-M[chans[2]]) zero_block[chans[3]];
#             M[chans[1]] zero_block[chans[2]] (.-M[chans[3]])]
#     else
#         error("there should be either 2 or 3 channels")
#     end

#     dim, vals, sols = nullspace_qr(matrix)
#     sols_vec = [sols[:, i] for i in 1:dim]

#     return BootstrapEquations{T,U}(S, chans, fields, matrix, dim, vals, sols_vec)
# end
