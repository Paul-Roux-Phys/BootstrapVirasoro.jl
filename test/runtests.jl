using JuliVirBootstrap
using Test

@testset "JuliVirBootstrap.jl" begin
    
    #Display a central charge
    charge = CentralCharge("β",sqrt(big(2)))
    field = Field(charge,"Δ",0.5,diagonal=true)

end
