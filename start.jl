using Plots, FFTW, DSP

# Number of points
N = 2^14 - 1
# Sample period
Ts = 1 / (1.1 * N)
# Start time
t0 = 0
tmax = t0 + N * Ts
# time coordinate
t = t0:Ts:tmax

# signal
signal1 = 3*sin.(2π * 10 .* t .+ 1) # sin (2π f t)
signal2 = 5*sin.(2π * 5 .* t)
signal3 = 10*sin.(2π * 54 .* t)


signal = signal1 + signal2 + signal3
# Fourier Transform of it
F = fft(signal) |> fftshift
freqs = fftfreq(length(t), 1.0/Ts) |> fftshift

# plots
sum_td = plot(t[1:5000], signal[1:5000], title="Signal sum")
time_domain = plot(t[1:5000],[signal1[1:5000],signal2[1:5000],signal3[1:5000]], title = "Signal")
freq_domain = plot(freqs, abs.(F), title = "Spectrum", xlim=(0, +100))
plot(time_domain, freq_domain, sum_td, layout=3)

plot_sine_sum(0:0.001:π, [create_sine(3, 10, 0), create_sine(5, 15, 0)])

c

#Block to animate sines
n = 1000
t = range(0, 2π, length = n)
sig = 5cos.(2π*10t)

@gif for i ∈ 1:n
        plot1 = plot(t, [10cos.(i.+2π*5*t), 5cos.(i.+π*t)], xlim=(0,+π))
end every 50



@gif for i ∈ 0:2π
        plot(sin.(2π * 40 .* 0:0.01:2π+0.01 .+ i))
end


function create_sine(amp, freq, phas, inter=0:0.001:π)
        signal = amp*sin.(2π * freq .* inter .+ phas)
        return signal
end

function plot_sine_sum(inter, sigs)
#        for item in sigs
#                signal += item;
#        end
        sum = plot(inter, signal, title="Input sig")
        decor = plot(inter, signal, Title = "Decorelated")
        plot(decor, sum, layout=2)
end
