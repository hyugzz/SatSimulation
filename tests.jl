#Testing norme
using Revise
includet("./signal_generator.jl")
includet("./gravity_utils.jl")

using Main.signal_utils, Main.MyGravity, Test

#Init of variables for test modularity
cofdm = signal_utils.cofdm([], 8, 15000.0, 500000.0, 4e-3, 4500, 2e-6)
byte = [Bool[] for _ in 1:4]
k = [Bool[] for _ in 1:4]
plot = [Plot() for _ in 1:4]

@testset  "Testing large range of Doppler shift" for shift in 0.85:0.1:1.25
    @info "Starting test for Doppler : $shift"
    # Generate random data
    for i in 1:4
      byte[i] = rand(Bool, 8)
    end
    print(byte)
    signal_utils.define_data(cofdm, byte)

    # Generate signal
    sig = signal_utils.generate_signal(cofdm, shift)
    plot , fft_data, fft_freq = signal_utils.AnalyzeSignal(cofdm, sig, 4)
    for i in 1:4
      @test isempty(signal_utils.reconstruct_data(cofdm, fft_data[i], fft_freq))
    end
    true_shift, k = signal_utils.ShiftSignal(fft_data, cofdm, fft_freq)
    @test true_shift != 0.0
    @test k == cofdm.data
    @test isapprox(shift, true_shift ; atol=0.21)
    if(!isapprox(shift, true_shift ; atol=0.21) || k != cofdm.data)
        @warn "The algorithm didn't work for $shift.\nThe returned result is $k\nThe estimated shift is $true_shift"
    end
  end
end;
