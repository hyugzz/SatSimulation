#Include dependencies
using Revise, Debugger
includet("./gravity_utils.jl")
includet("./signal_generator.jl")
using Main.signal_utils, Main.MyGravity, Plots, SatelliteToolbox, FFTW, PlotlyBase

#Variables initiation
const global c = 2.99792458e8
plotly()
tles = read_tle("starlink.txt")
sig = signal_utils.cofdm([true, false, true, true, false, true, false, false], 1, 15000.0, 5.00e9, 1.0)
doppler_ratio = Float64[]

#SGP4 - Generating orbit
orbit = init_orbit_propagator(Val(:sgp4), tles[1])
r, v = propagate!(orbit, collect(0:1:150)*60)

#Initializing cartesian objects
receiver = MyGravity.space_obj([0.0,0.0,0.0], [0.0,0.0,0.0], [0.0,0.0,0.0], 0.0)
sender   = MyGravity.space_obj([1000.0,0.0,0.0], [1000.0,0.0,0.0], [-10000000000.0,0.0,0.0], 0.0)

#Calculate doppler ratio for the two space objects
for i in 1:1000
    append!(doppler_ratio, signal_utils.Doppler(sender, receiver))
    MyGravity.update_pos(sender, 1.0)
end
plot(doppler_ratio)

#Generating signal
a = signal_utils.cofdm([true, false, true, true, false, true, false, true], 4, 100.0, 1000.0, 1.0)
sig = signal_utils.generate_signal(a, 1.0)
plot(sig)

#Calculating DFT
trans = FFTW.r2r(sig[1:10001], FFTW.DHT)
plot(fftshift(-trans))
plot(trans)


x,y,z = MyGravity.DisassembleVector(r)
x1,y1,z1 = MyGravity.DisassembleVector(v)
Main.MyGravity.set_pos(sender, [x[1], y[1], z[1]])
Main.MyGravity.set_spe(sender, [x1[1], y1[1], z1[1]])

a = signal_utils.generate_signal(sig, 1.0)
signal_utils.Doppler(1.0, sender, receiver)

bit_0 = sig[1:10001]
plots = signal_utils.AnalyzeSignal(sig, 8)
plot(plots[1],plots[2],plots[3])

for i in 1:151
    append!(dist, abs(signal_utils.NormeVecteur(sender.p - receiver.p)))
    append!(doppler_ratio, signal_utils.Doppler(1.0, sender, receiver))
    Main.signal_utils.generate_signal(sig, doppler_ratio)
    Main.MyGravity.set_pos(sender, [x[i], y[i], z[i]])
    Main.MyGravity.set_spe(sender, [x1[i], y1[i], z1[i]])
end

doppler_ratio = Float64[]
posx = Float64[]
spex = Float64[]
accx = Float64[]
dist = Float64[]
for i in 1:151
    append!(posx, sender.p[1])
    append!(spex, sender.v[1])
    append!(accx, sender.a[1])
    append!(dist, abs(signal_utils.NormeVecteur(sender.p - receiver.p)))
    append!(doppler_ratio, signal_utils.Doppler(1.0, sender, receiver))
    Main.MyGravity.set_pos(sender, [x[i], y[i], z[i]])
    Main.MyGravity.set_spe(sender, [x1[i], y1[i], z1[i]])
end
Plots.plot(x,y,z)
Plots.plot(doppler_ratio)
popfirst!(doppler_ratio)
