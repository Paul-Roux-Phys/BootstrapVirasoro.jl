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

# We can also use only 2 channels
solve_bootstrap(Channels(Dict(:s => Cs, :t => Ct)))

# Rels also works with 2 channels:
solve_bootstrap(Channels(Dict(:s => Cs, :t => Ct)), rels=:st)

# We can fix the value of a certain number of structure constants
# The format is an array of triples (channel, field, value)
fixed = [(:s, VPpm[1], 2), (:t, VPpm[1], 2), (:u, Vd[3], 0)]
solve_bootstrap(Channels(Cs, Ct, Cu), fix=fixed)

# If we need to, we can first create the linear system, and then solve it:
sys = BootstrapSystem(Channels(Cs, Ct, Cu))
# sys is a struct that contains all of the cached data needed to create the linear
# system: values of the blocks, prefactors per field and per channel, cache associated
# to the moduli.
# This is sometimes useful, because we can modify something in the bootstrap data, for instance
# modify the momentum of a diagonal field, without recomputing everything.
solve(sys, rels=:st) # this constructs and solves the linear system, and returns the
                     # corresponding structure constants
# We can remove a certain number of fields:
deleteat!(sys.unknowns.u, 2) # remove the second field in unknowns of u channel, which here is V^d_<1,5>
mysol2 = solve(sys)

return sol # return the structure constants
