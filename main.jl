#Include dependencies
using Revise, Debugger
includet("./gravity_utils.jl")
includet("./signal_generator.jl")
using Main.signal_utils, Main.MyGravity, Plots, SatelliteToolbox, FFTW, PlotlyBase

#Variables initiation
const global c = 2.99792458e8
plotly()
sequences = [rand(Bool, 8,4) for _ in 1:10]
                        #8 channel data, 2 channels sync, 5MHz carrier, 15KHz subcarrier, 15 MHz sampling frequency
a = signal_utils.cofdm([[true, false, true, true, false, true, false, false],
                        [true, true, false, true, false, true, false, true],
                        [false, true, false, false, false, true, false, true],
                        [false, false, true, false, true, false, true, true]],
                        [[false,true, true], [true,true, false], [true, false, true], [true, false, false]],
                        8, [1,8, 11], 15000.0, 5.00e6, 0.001, 2, 1.5e7)
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

#Generating signal
sig = signal_utils.generate_signal(a, 0.85)
#DFT transform of the signal
plots, fft_data, fft_freq = signal_utils.AnalyzeSignal(a, sig, 4)

# Is the is shift correct for this bit
for i in 1:4
    bits = signal_utils.reconstruct_data(a, signal_utils.SmoothedZscoreAlgo(fft_data[i], 3, 4, 0.5)[1], fft_freq)
    if 0 == signal_utils.check_difference(a, bits , i)
        print("$i good")
    else
        print("$i - $(signal_utils.check_difference(a, bits , i))")
    end
end
true_shift, k,  = signal_utils.ShiftSignal(fft_data, a, fft_freq)
signal_utils.check_difference(a,k)
