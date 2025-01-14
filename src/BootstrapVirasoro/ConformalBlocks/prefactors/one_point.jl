qfromτ(τ) = exp(2*im*(π*τ))

function blockprefactor_chiral(d::OneDimension, b::BlockChiral, τ)
    q = qfromτ(τ)
    return q^b.channel_dimension.δ / etaDedekind(τ)
end
