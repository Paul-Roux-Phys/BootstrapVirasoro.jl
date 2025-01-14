#===========================================================================================
Struct FourPointBlockSphere
===========================================================================================#

"""
    FourPointBlockSphere{T}

Composite type that represents a four-point conformal block:
a channel and a field propagating in the channel. The external fields and central charge are
provided in a `FourPointCorrelation` object.

# Example

```julia-repl
julia> c = CentralCharge(:c,0.5); V = Field(c, :δ, 0.6, diagonal = true);
julia> FourPointBlockSphere(:s, V)
Four-point block
Channel:        s
Channel Field:
Diagonal field of dimension:
  Δ = 0.5791666666666667 + 0.0im
  P = 0.0 + 0.7745966692414834im
  δ = 0.6000000000000001 + 0.0im
  p = 0.7745966692414834 + 0.0im
```
"""
struct FourPointBlockSphere{T}

    corr::FourPointCorrelation{T}
    channel::Symbol
    channelfield::Field{T}
    _seriescoeffs_lr::Tuple{Dict{Tuple{Int, Int, Int}, T}, Dict{Tuple{Int, Int, Int}, T}}
    _Nmax::Int

end

"""Permute the external fields to get t- or u-channels from s-channel."""
function permute_ext_fields(corr::FourPointCorrelation, chan::Symbol)::FourPointCorrelation
    Vs=corr.fields
    Vs = @match chan begin
        :s => [Vs[1], Vs[2], Vs[3], Vs[4]]
        :t => [Vs[1], Vs[4], Vs[3], Vs[2]]
        :u => [Vs[1], Vs[3], Vs[2], Vs[4]]
        _ => error("The parameter $chan is not a valid channel")
    end
    return FourPointCorrelation(corr.charge, Vs)
end


function FourPointBlockSphere(corr::FourPointCorrelation{T}, s::Symbol, V::Field{T}; Nmax=10) where {T}
    corr_permuted = permute_ext_fields(corr, s)
    coeff_left = Dict{Tuple{Int,Int,Int}, T}( ((N, m, n) => computeCNmn(N, m, n, corr_permuted, left))
                       for n in 1:Nmax for m in 1:Nmax for N in 1:Nmax
                           if m*n <= N)
    coeff_right = Dict{Tuple{Int,Int,Int}, T}( ((N, m, n) => computeCNmn(N, m, n, corr_permuted, right))
                        for n in 1:Nmax for m in 1:Nmax for N in 1:Nmax
                            if m*n <= N)
    return FourPointBlockSphere(corr_permuted, s, V, (coeff_left, coeff_right), Nmax)
end

function Base.show(io::IO, block::FourPointBlockSphere)
    print("Four-point block, for the ")
    show(block.corr); print("\n")
    println("Channel:\t$(block.channel)")
    println("Channel Field:")
    show(block.channelfield)
end

#===========================================================================================
Get t- and u- channel blocks from s-channel block
===========================================================================================#
"""Prefactor to get t- or u-channel blocks from the s-channel block"""
function channelprefactor_chiral(block::FourPointBlockSphere, x)
    @match block.channel begin
        :s => 1
        :t => 1
        :u => 1/x^(-2*block.corr.fields[1].Δ[left])
    end
end

function channelprefactor_non_chiral(block::FourPointBlockSphere, x)
    return channelprefactor_chiral(block, x)*channelprefactor_chiral(block, conj(x))
end

"""Sign (-1)^{S_1+S_2+S_3+S_4} when changing from s to t or u channels"""
function channel_sign(block::FourPointBlockSphere, x)
    @match block.channel begin
        :s => 1
        :t => 1 # (-1)^(sum(spin.(corr.fields)))
        :u => 1 # (-1)^(sum(spin.(corr.fields)))
    end
end

"""Cross-ratio at which to evaluate the s-channel block to get t- or u-channel block"""
function crossratio(channel, x)
    @match channel begin
        :s => x
        :t => 1-x
        :u => 1/x
    end
end

#===========================================================================================
Set prefactors, relate the cross-ratio x and the elliptic nome q
===========================================================================================#
"""Nome `q` from the cross-ratio `x`"""
qfromx(x) = exp(-oftype(x, π)*ellipticK(1-x)/ellipticK(x))

"""Cross ratio `x` from the nome `q`"""
xfromq(q) = jtheta2(0,q)^4 / jtheta3(0,q)^4

"""Prefactor for getting the block F from H. The argument `lr` indicates if we are working
with a left or right moving block"""
function blockprefactor(block::FourPointBlockSphere, x, lr)

    corr = block.corr
    c = corr.charge.c
    e0 = - corr.fields[1].δ[lr] - corr.fields[2].δ[lr] - (c-1)/24
    e1 = - corr.fields[1].δ[lr] - corr.fields[4].δ[lr] - (c-1)/24
    e2 = sum(corr.fields[i].δ[lr] for i in 1:4) + (c-1)/24
    q=qfromx(x)

    return Complex(x)^e0 * (Complex(1-x))^e1 * jtheta3(0,q)^(-4*e2) * (16*q)^block.channelfield.δ[lr]
end

"""Degenerate dimensions"""
δrs(r, s, B) = -1/4 * (B*r^2 + 2*r*s + s^2/B)
βm1P(B, r, s) = 1/2*(r+s/B) # \beta^{-1}P_{(r,s)}

"""Factor \ell_{(r,s)} that appears in logarithmic blocks"""
function ell(corr, r, s)
    c = corr.charge
    B, β = c.B, c.β
    βm1P_ext = [[corr.fields[i].P[left]/β for i in 1:4], [corr.fields[i].P[right]/β for i in 1:4]]

    term1(j) = digamma_reg(-2*βm1P(B, r, j)) + digamma_reg(2*βm1P(B, r, -j))

    res = -big(4)*oftype(B, π)/tan(oftype(B, π)*s/B)

    term3(j, lr, pm1, pm2, a, b) = digamma_reg(1/2 + (lr == left ? -1 : 1)*βm1P(B, r, j) + pm1*βm1P_ext[lr][a] + pm2*βm1P_ext[lr][b])

    return res + 4*sum(term1(j) for j in 1-s:s) -
        sum(term3(j, lr, pm1, pm2, a, b)
                        for pm1 in (-1,1)
                        for pm2 in (-1,1)
                        for j in 1-s:2:s-1
                        for (a,b) in ((1,2), (3, 4))
                        for lr in (left, right)
        )
end

function H_series_coeffN(block, lr, N;
                         der=false, reg=false)

    V = block.channelfield
    P = V.P[lr]
    β = block.corr.charge.β

    res = 0
    if !reg && !der
        for m in 1:N
            for n in 1:N
                if m*n <= N
                    Pmn = (β*m-n/β)/2
                    res += block._seriescoeffs_lr[lr][(N, m, n)]/(P^2-Pmn^2)
                end
            end
        end
        return res
    elseif !reg && der
        for m in 1:N
            for n in 1:N
                if m*n <= N
                    Pmn = (β*m-n/β)/2
                    res -= block._seriescoeffs_lr[lr][(N, m, n)]/(P^2-Pmn^2)^2
                end
            end
        end
        return 2*P*res
    elseif reg && V.isKac && V.r%1 == V.s%1 == 0 && V.r > 0 && (lr == left && V.s > 0 || lr == right && V.s < 0)
        for m in 1:N
            for n in 1:N
                if m*n <= N
                    Pmn = (β*m-n/β)/2
                    if m != V.r || n != abs(V.s)
                        res += block._seriescoeffs_lr[lr][(N, m, n)]/(P^2-Pmn^2)
                    else
                        res -= block._seriescoeffs_lr[lr][(N, m, n)]/(4*P^2)
                    end
                end
            end
        end
        return res
    else
        error("Trying to compute the derivative of a regularised block")
    end
end

"""
    H_series(block, lr;
      der = false, reg = false)

Compute the coefficients of the series expansion of the function ``H(q,δ)``. If der=true, compute instead the series of the derivative of H with respect to P. If reg=true, compute instead the P dependent part of the coefficients of ``H^{\\text{reg}}``.
"""
@memoize function H_series(block::FourPointBlockSphere, lr;
                  der=false, reg=false)

    @assert !(der && reg) "you should not compute the derivative of a regularised block"

    if !der
        return vcat(1, [H_series_coeffN(block, lr, N, der=der, reg=reg) for N in 1:block._Nmax]) # H = 1 + series
    else
        return [H_series_coeffN(block, lr, N, der=der, reg=reg) for N in 1:block._Nmax]
    end
end

#===========================================================================================
Compute the conformal block
===========================================================================================#
"""
    block_chiral(x, Nmax, block, lr)

Compute the chiral conformal block

``\\mathcal F^{(\\text{chan})}_{\\delta}(x)``

where `chan` is `s`, `t`, or `u`."""
function block_chiral(x, block::FourPointBlockSphere, lr;
                      der=false, reg=false)
    chan = block.channel
    P = block.channelfield.P[lr]
    x_chan = crossratio(chan, x)

    q = qfromx(x_chan)
    h = evalpoly(16*q, H_series(block, lr, der=der, reg=reg))
    if reg
        V = block.channelfield
        # h += log(sq)*sum(block._seriescoeffs_lr[lr](N, V.r, abs(V.s))*(sq)^N for N in V.r*abs(V.s):block.Nmax) # log(16q) term in log(16q) - 1/4P^2
        h += log(16*q)*(16*q)^(V.r*abs(V.s))*evalpoly(16*q, [block._seriescoeffs_lr[lr][(N, V.r, abs(V.s))] for N in V.r*abs(V.s):block._Nmax]) # log(16q) term in log(16q) - 1/4P^2
    elseif der
        h += 2*P*log(16*q)*evalpoly(16*q, H_series(block, lr, der=false)) # H^der = 2Plog(16q)H + H'
    end

    return channelprefactor_chiral(block, x_chan) * blockprefactor(block, x_chan, lr) * h
end

block_chiral(x, block, lr; der=false, reg=false) = block_chiral(x, block, lr; der=der, reg=reg)

"""
    block_non_chiral(x, Nmax, block)

Compute the non-chiral conformal block G_(r,s) in the s channel.

TODO: regularise R_(r,s) / \bar{R}_(r,s)
"""
function block_non_chiral(x, block::FourPointBlockSphere)

    x_chan = crossratio(block.channel, x)
    Vchan = block.channelfield

    if !Vchan.isKac || (Vchan.isKac && (Vchan.r%1 != 0 || Vchan.s%1 != 0 || spin(Vchan) == 0)) # non-logarithmic block

        return block_chiral(x_chan, block, left) * block_chiral(conj(x_chan), block, right)

    elseif 0 == 1 # accidentally non-logarithmic block
        return
    else
        # logarithmic block
        corr = block.corr

        r, s = Vchan.r, Vchan.s

        @assert !(Vchan.r < 0 || Vchan.s < 0) "Trying to compute a logarithmic block with a negative index: r=$(Vchan.r), s=$(Vchan.s) .
                                               This goes against the chosen convention"
        c = corr.charge
        block1 = block
        block2 = FourPointBlockSphere(corr, :s, Field(c, r=r, s=-s), Nmax=block._Nmax) # block with momenta (P_(r,-s), P_(r,s)) in the channel

        F_Prms = block_chiral(x_chan, block2, left) # F_{P_(r,-s)}
        F_Prms_bar = block_chiral(conj(x_chan), block1, right) # \bar F_{P_(r,-s)}
        F_der_Prms = block_chiral(x_chan, block2, left, der=true) # F'_{P_(r,-s)}
        F_der_Prms_bar = block_chiral(conj(x_chan), block1, right, der=true) # \bar F'_{P_(r,-s)}
        F_reg_Prs = block_chiral(x_chan, block1, left, reg=true) # F^reg_{P_(r,s)}
        F_reg_Prs_bar = block_chiral(conj(x_chan), block2, right, reg=true) # \bar F^reg_{P_(r,s)}

        R = Rmn(r, s, corr, left) # Vchan.P[left] = P_(r,s)
        R_bar = Rmn(r, s, corr, right)

        term1 = (F_reg_Prs - R*F_der_Prms)*F_Prms_bar
        term2 = R/R_bar*F_Prms*(F_reg_Prs_bar - R_bar*F_der_Prms_bar)
        term3 = -R*ell(corr, r, s)*F_Prms*F_Prms_bar

        # return F_Prms, F_Prms_bar, F_der_Prms, F_der_Prms_bar, F_reg_Prs, F_reg_Prs_bar, ell(corr, r, s), R, R_bar
        return channel_sign(block, x)*(term1+term2+term3)
    end
end

"""
    block_non_chiral(x, Nmax, block, corr)

Compute the non-chiral conformal block G_(r,s) in the channel indicated in `block`.

TODO: regularise R_(r,s) / \bar{R}_(r,s)
"""
