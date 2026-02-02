# Solving bootstrap equations

The package provides functionality to create and solve linear systems of bootstrap equations.

```julia
using BootstrapVirasoro

# set the precision of high-precision floats
setprecision(BigFloat, 40, base=10)

β    = big"0.3" + big"0.6"*im
# define a central charge
c    = CC(β = β)

V12  = Field(c, r=1, s=2, diagonal=true) # the external, degenerate field
P    = big"0.67"
VP   = Field(c, P=P) # the external diagonal field
VPpm = [Field(c, P = P + pm*1/(2β)) for pm in (-1, 1)] # shifted fields
Vd   = [Field(c, r=1, s=1, diagonal=true),
        Field(c, r=1, s=3, diagonal=true),
        Field(c, r=1, s=5, diagonal=true)] # degenerate fields
# list of all fields that can go in channels (vcat concatenates vectors.)
V    = vcat(VPpm, Vd) 

Cor  = Correlation(V12, VP, V12, VP, 70) # define the correlation

# We need to define ChannelSpectrum objects for each channel.
# They are created as ChannelSpectrum(correlation, channel, listoffields, function)
# where the function is a function that takes as input a field and outputs a block object.
fs(V) = NCBlock(Cor, :s, V) # map fields to blocks
Cs    = ChannelSpectrum(Cor, :s, V, fs) # create s-channel spectrum
ft(V) = NCBlock(Cor, :t, V)
Ct    = ChannelSpectrum(Cor, :t, V, ft) # t-channel
fu(V) = NCBlock(Cor, :u, V)
Cu    = ChannelSpectrum(Cor, :u, V, fu) # u-channel

sol = solve_bootstrap(Channels(Cs, Ct, Cu)) # setup and solve the crossing equations
sol.str_cst # return the structure constants

# you can display the structure constants with
# println(sol.str_cst)
# or
# show(stdout, sol.str_cst)
# changing stdout to a file IO will print to a file instead.
```

The output of running this code should be a colored table indicating the values of the structure constants,
and an estimation of the numerical error. You can see this output in the jupyter notebook on the
[github repo](https://github.com/Paul-Roux-Phys/BootstrapVirasoro.jl/blob/main/examples/bootstrap.ipynb).
