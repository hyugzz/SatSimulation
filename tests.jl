#Testing norme
using Revise
include("./signal_generator.jl")
include("./gravity_utils.jl")

using Main.signal_utils, Main.MyGravity, Test

#Init of variables for test modularity
cofdm = signal_utils.cofdm


@testset "Testing large range of Doppler shift" for shift in 0.85:0.00001:1.25 begin
    @info "Starting test for Doppler : $shift"
    # Generate random data
    byte = rand(Bool, 8,4)
    signal_utils.define_data(cofdm, byte)
    # Generate signal
    sig = signal_utils.generate_signal(cofdm, shift)
        plot , fft_data, fft_freq = singal_utils.AnalyzeSignal(cofdm, sig, 4)
        k, true_shift = signal_utils.ShiftSignal(fft_data, cofdm, fft_freq)
        @test k == cofdm.data
        @test isapprox(shift, true_shift ; atol=0.0001)
    if(!isapprox(shift, true_shift ; atol=0.0001) || k != cofdm.data)
        @warn "The algorithm didn't work for $shift.\nThe returned result is $k\nThe estimated shift is $true_shift"
    end
end;
