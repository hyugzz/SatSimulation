#Include dependencies
using Revise, Debugger
includet("./gravity_utils.jl")
includet("./signal_generator.jl")
using Main.signal_utils, Main.MyGravity, Plots, SatelliteToolbox, FFTW, PlotlyBase

#Variables initiation
const global c = 2.99792458e8
plotly()
                        #8 channel, 5MHz carrier, 15KHz subcarrier, 15 MHz sampling frequency
a = signal_utils.cofdm([[true, false, true, true, false, true, false, false],
                        [true, true, false, true, false, true, false, true],
                        [false, true, false, false, false, true, false, true],
                        [false, false, true, false, true, false, true, true]],
                        8, 15000.0, 5.00e6, 0.001, 1500.0, 1.5e7)
tles = read_tle("starlink.txt")
doppler_ratio = Float64[]

#SGP4 - Generating orbit
orbit = init_orbit_propagator(Val(:sgp4), tles[1])
r, v = propagate!(orbit, collect(0:1:6000)*15)

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
sig = signal_utils.generate_signal(a, 1.1)
plot(sig)

signal_utils.AnalyzeSignal(a, sig, 4)

plots = Plot[]
fft_data = Vector{Array{Float64}}

#DFT transform of the signal
plots, fft_data, fft_freq = signal_utils.AnalyzeSignal(a, sig, 4)
for i in 1:length(plots)
    display(plots[i])
end
plot(fft_freq, fft_data[1])

fft_data[1]
plot(fft_data[1])
k = [Vector{Bool}() for _ in 1:length(a.data)]
for i in 1:4
    k[i] = signal_utils.reconstruct_data(a, fft_data[i], fft_freq)
end


shifted_ffts = [Vector{Bool}()  for _ in 1:length(fft_data)]
bits = Float64[]
for i in 1:length(fft_data)
    shifted_fft = Vector{Float64}()
    for k in 1:length(fft_data[i])    # FFT frequency shift loop
        if trunc(Int,k*1.1) < length(fft_data[i])
            if trunc(Int, k*1.1) == 0
                push!(shifted_fft, fft_data[i][1])
            else
                push!(shifted_fft, fft_data[i][trunc(Int, k*1.1)])
            end
        else
            push!(shifted_fft, 0.0)
        end
    end
    shifted_ffts[i] = signal_utils.reconstruct_data(a, shifted_fft, fft_freq)
end
popfirst!(shifted_fft)
# Is the is shift correct for this bit
bits = signal_utils.reconstruct_data(a, shifted_fft, fft_freq)

plot(fft_freq,shifted_fft)

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
tmp_fft = fft_data[begin:end][trunc(Int, length(fft_data[begin:end])/2):length(fft_data[begin:end])]
true_shift, k = signal_utils.ShiftSignal(fft_data, a, fft_freq)
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
