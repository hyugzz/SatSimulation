include("./gravity_utils.jl")

module signal_utils
    import Main.MyGravity
    using Plots, FFTW
    const global c = 2.99792458e8
    mutable struct cofdm
        data::Vector{Bool}  #Array of boolean data transmitted
        nb_channel::Int64   #Number of channels of the cofdm
        interval::Float64   #Interval between cofdm channels
        frequency::Float64  #Base frequency for transmission
        transmission_rate::Float64  #Time between two signal modulations
        peak_tolerance::Int64 #Accepted shifting in frequency during FFT to reconstruct a bit
    end
    "Defines data of the frame to be transmitted"
    function  define_data(self::cofdm, in::Vector{Bool})
        self.data = in;
    end
    "Generates signal with cofdm, data is modulated only in frequency"
    function generate_signal(self::cofdm, doppler::Float64)
        #Generate sine wave
        signal = Float64[]
        tmp_signal = Float64[]
        tmp_freq = 0.0;
        for i in 1:(2π*5e3*self.transmission_rate)+1
            append!(tmp_signal, 0.0)
        end
        for bit in self.data
            if bit==1           # 1 --> Higher frequency
                for i in 0:self.nb_channel-1
                    tmp_freq = doppler * (self.frequency+self.interval/2+i*self.interval)
                    tmp_signal += create_sine(tmp_freq, 1/self.transmission_rate)
                end
            elseif bit==0       # 0 --> Lower frequency
                for i in 0:self.nb_channel-1
                    tmp_freq = doppler * (self.frequency+i*self.interval)
                    tmp_signal += create_sine(tmp_freq,1/self.transmission_rate)
                end
            end
            append!(signal, tmp_signal)
            empty!(tmp_signal)
            for i in 1:(2π*5e3*self.transmission_rate)+1 #Resizing the vector for the next loop execution
                append!(tmp_signal, 0.0)
            end
        end
        return signal
    end
    "Returns the length of the vector"
    function NormeVecteur(vec::Vector{Float64})
        tmp = 0
        for i in vec
            tmp = tmp + (i*i)
        end
        if abs(sqrt(tmp)) <= c
            return abs(sqrt(tmp))
        else
            return c
        end
    end
    "Returns the projection of sender speed vector onto the distance vector of sender and receiver"
    function ProjectiveRelativeSpeed(sender::MyGravity.space_obj, receiver::MyGravity.space_obj)
        ϕ = (sender.p[1]*receiver.p[1]+sender.p[2]*receiver.p[2]+sender.p[3]*receiver.p[3]) / (NormeVecteur(sender.v) * NormeVecteur(receiver.p - sender.p))
        v_proj = cos(ϕ).* sender.v;
        return v_proj
    end
    "Returns the doppler ratio (f˕receiver/f˕sender )received by an observer in his referential"
    function Doppler(sender::MyGravity.space_obj, receiver::MyGravity.space_obj)
        lorentz = sqrt(1 - ((NormeVecteur(sender.v) ^ 2)/(c ^ 2)))
        doppler = (1 - NormeVecteur(ProjectiveRelativeSpeed(sender, receiver)) / c)
        coeff = lorentz * doppler;
        return doppler
    end
    "Returns a amp*sin(2π*freq*t + phase) with t = [0 ; inter]"
    function create_sine(freq, inter::Float64, amp=1.0, phas=0.0)
            interval = 0:1e-8:0.00005*2π
            signal = amp*sin.(2π * freq .* interval .+ phas)
            return signal
    end
    "Returns the DFTs of the input signal sequence"
    function AnalyzeSignal(sig::Vector{Float64} , nb_bit)
        #number of points
        N = trunc(Int,(length(sig)/nb_bit)) - 1
        #sample period
        Ts = 1e-8
        t0 = 0
        tmax = t0 + N * Ts
        t = t0:Ts:tmax
        plots = Plots.Plot[]
        ffts = [Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],]
        vcat(ffts)
        freqs = fftfreq(N+1, 1.0/Ts) |> fftshift
        for i in 0:nb_bit-1

            bit = sig[trunc(Int,(length(sig)/nb_bit)*i)+1:trunc(Int, (length(sig)/nb_bit)*(i+1))]
            F = FFTW.r2r(bit, FFTW.DHT) |> fftshift
            freqs = fftfreq(trunc(Int,length(bit)), 1/Ts)
            freq_domain = plot(freqs, abs.(F), title = "Bit $i")
            push!(plots, freq_domain)
            append!(ffts[i+1],F)
        end
        return plots, ffts, freqs
    end
    "Reconstruct the cofdm data from the FFT"
    function reconstruct_data(self::cofdm, fft::Vector{Float64}, freqs::Frequencies{Float64})
        sym = fft[trunc(Int, length(fft)/2):end]
        waves = Int64[]
        around = 2
        is_local_peak = false
        is_down = false
        #Sampling bloc
        for i in 1:trunc(Int, length(sym))-1
            if sym[i] > sym[i+1] && !is_down   #Found peak, and value rising
                is_down = true
                for j in i-around:i+around
                    if j > 0 && j < length(sym)
                        if sym[j] > sym[i] #Local peak
                            is_local_peak = true
                            break
                        end
                        is_local_peak = false
                    end
                end
                if !is_local_peak && sym[i] > 50
                    push!(waves,trunc(Int,freqs[i]))
                end
            else
                is_down = false
            end
        end #Peak filled
        #Bits reconstruction
        bits = []
        for channel in 1:length(waves)
            if abs(waves[channel]-(self.frequency+(channel-1)*self.interval)) < (self.interval/2)   #If the peak is in the cofdm channel
                #0 --> False
                push!(bits, false)
            elseif abs(waves[channel]-(self.frequency+(channel-1)*self.interval)) < (self.interval)
                # 1 --> True
                push!(bits, true)

            end
        end
        return bits
    end
    "This function shifts the peaks of the synchronization signal in order to calculate the Doppler effect coefficient"
    function ShiftSignal(ffts::Vector{Vector{Float64}}, self::cofdm, freq::Frequencies{Float64})
        doppler = 0.0;
        is_correct = false
        true_shift = 0.0
        for shift in 0.999973:0.000001:1.35   # Value for maximum shift (8Km/s)
            for bit in 1:length(ffts)
                shifted_fft = Float64[]
                for k in 1:length(ffts[bit])
                    if trunc(Int,k*shift) < length(ffts[bit])
                        if trunc(Int, k*(1/shift)) == 0
                            push!(shifted_fft, ffts[bit][1])
                        else
                            push!(shifted_fft, ffts[bit][trunc(Int, k*(1/shift))])
                        end
                    end
                end
                # Is the is shift correct for this bit
                bits = reconstruct_data(self, shifted_fft, freq)
                if isempty(bits)
                    is_correct = false
                    break
                elseif bits[1] == self.data[bit]  # If the peak is shifted correctly
                    printstyled("Bit $(bit-1) is correctly shifted with $shift\n", color=:green)
                    println("Its value is $(bits[1])")
                    is_correct = true
                    true_shift = shift
                else
                    is_correct = false
                    break
                end
            end
            if is_correct
                break
            end
        end
        return true_shift
    end
end
