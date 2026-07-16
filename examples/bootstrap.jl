using BootstrapVirasoro

# set the precision of high-precision floats
setprecision(BigFloat, 40, base=10)

β    = big"0.8" + big"0.1"*im
# define a central charge
c    = CC(β = β)

V12  = Field(c, r=1, s=2, diagonal=true) # the external, degenerate field
P    = big"0.67"
VP   = Field(c, P=P) # the external diagonal field
VPpm = [Field(c, P = P + pm*1/(2β)) for pm in (-1, 1)] # shifted fields
Vd   = [Field(c, r=1, s=1, diagonal=true),
        Field(c, r=1, s=3, diagonal=true),
        Field(c, r=1, s=5, diagonal=true)] # degenerate fields
# list of all fields that can go in channels
# vcat concatenates arrays.
V    = vcat(VPpm, Vd) 

Cor  = Correlation(V12, VP, V12, VP, 60) # define the correlation

# We need to define ChannelSpectrum objects for each channel.
# They are created as ChannelSpectrum(correlation, channel, listoffields, function)
# where the function is a function that takes as input a field and outputs a block object.
f(V, chan) = NCBlock(Cor, chan, V) # map fields to blocks
Cs    = ChannelSpectrum(Cor, :s, VPpm, V -> f(V, :s)) # create s-channel spectrum
Ct    = ChannelSpectrum(Cor, :t, VPpm, V -> f(V, :t)) # t-channel
Cu    = ChannelSpectrum(Cor, :u, Vd,   V -> f(V, :u)) # u-channel

sol = solve_bootstrap(Channels(Cs, Ct, Cu)) # setup and solve the crossing equations

# you can display the structure constants with
# println(sol)
# or
# show(stdout, sol)
# changing stdout to a file IO will print to a file instead.

# We can impose that the s and t channels structure constants are equal.
# Other options for rels are :su, :tu, :stu
solve_bootstrap(Channels(Cs, Ct, Cu), rels=:st)

return sol # return the structure constants
