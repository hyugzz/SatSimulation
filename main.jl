#Include dependencies
using Revise, Debugger
includet("./gravity_utils.jl")
includet("./signal_generator.jl")
using Main.signal_utils, Main.MyGravity, Plots, SatelliteToolbox, FFTW, PlotlyBase

#Variables initiation
const global c = 2.99792458e8
plotly()
tles = read_tle("starlink.txt")
sig = signal_utils.cofdm([true, false, true, true, false, true, false, false], 1, 15000.0, 5.00e9, 1.0, 15.0)
doppler_ratio = Float64[]

#SGP4 - Generating orbit
orbit = init_orbit_propagator(Val(:sgp4), tles[1])
r, v = propagate!(orbit, collect(0:1:150)*60)

#Initializing cartesian objects
x,y,z = MyGravity.DisassembleVector(r)
x1,y1,z1 = MyGravity.DisassembleVector(v)
receiver = MyGravity.space_obj([0.0,0.0,0.0], [0.0,0.0,0.0], [0.0,0.0,0.0], 0.0)
sender   = MyGravity.space_obj([0.0,0.0,0.0], [x1[1], y1[1], z1[1]], [x[1], y[1], z[1]], 0.0)

#Calculate doppler ratio for the two space objects
for i in 1:length(x)
    append!(doppler_ratio, signal_utils.Doppler(sender, receiver))
    MyGravity.set_pos(sender, [x[i], y[i], z[i]])
    MyGravity.set_spe(sender, [x1[i], y1[i], z1[i]])
end
plot(doppler_ratio)

#Generating signal
a = signal_utils.cofdm([true, false, true, true, false, true, false, false], 8, 15000.0, 5.00e5, 1.0, 15.0)
sig = signal_utils.generate_signal(a, 1.0)
plot(sig)

plots = Plot[]
fft_data = Vector{Array{Float64}}

#DFT transform of the signal
plots, fft_data, fft_freq = signal_utils.AnalyzeSignal(sig, 8)
for i in 1:length(plots)
    display(plots[i])
end
plot(fft_data[1])

fft_data[1]
plot(fft_data[1])
k = [Int64[],Int64[],Int64[],Int64[],Int64[],Int64[],Int64[],Int64[]]
for i in 1:length(fft_data)
    k[i] = signal_utils.reconstruct_data(a, fft_data[i], fft_freq)
end

# Stats about orbit
length([fft_data[1]])
for i in 1:151
    append!(dist, abs(signal_utils.NormeVecteur(sender.p - receiver.p)))
    append!(doppler_ratio, signal_utils.Doppler(1.0, sender, receiver))
    Main.signal_utils.generate_signal(sig, doppler_ratio)
    Main.MyGravity.set_pos(sender, [x[i], y[i], z[i]])
    Main.MyGravity.set_spe(sender, [x1[i], y1[i], z1[i]])
end

shifted_fft = Float64[]
tmp_fft = fft_data[1:8][trunc(Int, length(fft_data[1:8])/2):length(fft_data[1:8])]
true_shift = signal_utils.ShiftSignal(fft_data, a, fft_freq)
for i in 2:trunc(Int, length(fft_data[1])/2)
    if i < length(tmp_fft)/2
        push!(shifted_fft, tmp_fft[trunc(Int, i*1/(0.5))])
    end
end
plot(shifted_fft)

plot(shifted_fft)
plot(tmp_fft)


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
