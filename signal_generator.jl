include("./gravity_utils.jl")

module signal_utils
    import Main.MyGravity
    using Plots, FFTW, DSP, Statistics
    const global c = 2.99792458e8
    mutable struct cofdm
        "Array of boolean data transmitted"
        data::Vector{Vector{Bool}}
        "Array of boolean data used for synchronization"
        sync_data::Vector{Vector{Bool}}
        "Number of channels of the cofdm"
        nb_channel::Int64
        "channels used for synchronization"
        sync_channels::Vector{Int64}
        "Interval between cofdm channels"
        interval::Float64
        "Base frequency for transmission"
        frequency::Float64
        "Time between two signal modulations"
        transmission_rate::Float64
        "Accepted number of errors in synchronization channels"
        error_tolerance::Int64
        "Sampling rate respecting Nyquist limit"
        sampling_rate::Float64
    end
    "Defines data of the frame to be transmitted"
    function  define_data(self::cofdm, in::Vector{Vector{Bool}}, sync::Vector{Vector{Bool}})
        self.data = in;
        if length(sync) == length(self.sync_channels)
            self.sync_data = sync
        elseif length(sync) >= length(self.sync_channels)
            self.data = sync[1:length(self.sync_channels)]
        else
            printstyled("Sync data is not the right size ($(length(self.sync_channels)))\nIgnore")
        end
    end
    "Generates signal with cofdm, data is modulated only in frequency"
    function generate_signal(self::cofdm, doppler::Float64)
        #Generate sine wave
        signal = Float64[]
        tmp_signal = Float64[]
        tmp_freq = 0.0;
        for i in 0:1/self.sampling_rate:self.transmission_rate
            append!(tmp_signal, 0.0)
        end
        nb_seq = 0
        for seq in self.data
            nb_seq += 1
            sync_count = 0
            for chan_bit in 1:length(seq)+length(self.sync_channels)
                is_done = false
                for i in 1:length(self.sync_channels)
                    if chan_bit == self.sync_channels[i]
                        if self.sync_data[nb_seq][i]==1           # 1 --> Higher frequency
                            tmp_freq = doppler * (self.frequency+self.interval/2+(chan_bit-1)*self.interval+self.interval/4)
                            tmp_signal += create_sine(tmp_freq, self.sampling_rate, self.transmission_rate)
                        elseif self.sync_data[nb_seq][i]==0       # 0 --> Lower frequency
                            tmp_freq = doppler * (self.frequency+(chan_bit-1)*self.interval+self.interval/4)
                            tmp_signal += create_sine(tmp_freq, self.sampling_rate, self.transmission_rate)
                        end
                        sync_count += 1
                        is_done = true
                    end
                end
                if is_done == false
                    if seq[chan_bit-sync_count]==1           # 1 --> Higher frequency
                        tmp_freq = doppler * (self.frequency+self.interval/2+(chan_bit-1)*self.interval+self.interval/4)
                        tmp_signal += create_sine(tmp_freq, self.sampling_rate, self.transmission_rate)
                    elseif seq[chan_bit-sync_count]==0       # 0 --> Lower frequency
                        tmp_freq = doppler * (self.frequency+(chan_bit-1)*self.interval+self.interval/4)
                        tmp_signal += create_sine(tmp_freq, self.sampling_rate, self.transmission_rate)
                    end
                end
            end
            # Channel filled, append the signal
            window = DSP.Windows.rect(length(tmp_signal), padding=trunc(Int, log10(self.frequency)), zerophase=false)
            tmp_signal = tmp_signal .* sin.(π/length(tmp_signal):π/length(tmp_signal):π)
            append!(signal, tmp_signal)
            empty!(tmp_signal)
            for i in 0:1/self.sampling_rate:self.transmission_rate #Resizing the vector for the next loop execution
                append!(tmp_signal, 0.0)
            end
        end
        return signal
    end
    "Returns a amp*sin(2π*freq*t + phase) with t = [0 ; inter]"
    function create_sine(freq, sampling_rate, upper_bound, amp=1.0, phas=0.0)
            interval = 0:1/sampling_rate:upper_bound
            signal = amp*sin.(2π * freq .* interval .+ phas)
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
    "Returns the DFTs of the input signal sequence"
    function AnalyzeSignal(self::cofdm, sig::Vector{Float64} , nb_bit)
        #number of points
        N = trunc(Int,(length(sig)/nb_bit)) - 1
        #sample period
        Ts = 1 / self.sampling_rate
        tmax =  self.transmission_rate
        t = 0.0:Ts:tmax
        plots = Plots.Plot[]
        ffts = [Vector{Float64}() for _ in 1:nb_bit]
        freqs = fftfreq(N+1, Ts) # |> fftshift
        for i in 0:nb_bit-1

            bit = sig[trunc(Int,(length(sig)/nb_bit)*i)+1:trunc(Int, (length(sig)/nb_bit)*(i+1))]
            F = FFTW.r2r(bit, FFTW.DHT) # |> fftshift
            freqs = fftfreq(trunc(Int,length(bit)), 1/Ts)
            freq_domain = plot(freqs, abs.(F), title = "Bit $i")
            push!(plots, freq_domain)
            append!(ffts[i+1],  abs.(F))
        end
        return plots, ffts, freqs
    end
    "Reconstruct the cofdm data from the FFT"
    function reconstruct_data(self::cofdm, fft::Vector{Float64}, freqs::Frequencies{Float64})
        sym = fft[1:trunc(Int, length(fft)/2)]
        sym_freq = freqs[1:trunc(Int, length(freqs)/2)]
        waves = Int64[]
        around = 4
        is_local_peak = false
        is_down = false
        #Peak detection algorithm
#        for i in 1:trunc(Int, length(sym))-1
#            if sym[i] > sym[i+1] && !is_down   #Found peak, and value rising
#                is_down = true
#                for j in i-around:i+around
#                    if j > 0 && j < length(sym)
#                        if sym[j] > sym[i] #Local peak
#                            is_local_peak = true
#                            break
#                        end
#                        is_local_peak = false
#                    end
#                end
#                if !is_local_peak && sym[i] > 800 # Threshold
#                    push!(waves,trunc(Int,sym_freq[i]))
#                end
#            else
#                is_down = false
#            end
#        end #Peak filled
        sig = sym
        #Cleaning peaks
        tmp = Int64[]
        waves = Float64[]
        for i in 2:length(sig)
            if sig[i-1] == 1.0 && sig[i] == 0.0 # If we are on the end of peak
                push!(waves, sym_freq[i])
            end
        end
        println(length(waves))
        #Bits reconstruction
        bits = Bool[]
        for channel in 1:length(waves)
            if abs(waves[channel]-(self.frequency+(channel-1)*self.interval)) < (self.interval) / 2  #If the peak is in the cofdm channel
                #0 --> False
                push!(bits, false)
            elseif abs(waves[channel]-(self.frequency+(channel-1)*self.interval)) < (self.interval)
                # 1 --> True
                push!(bits, true)

            end
        end
        return bits
    end
    "Z-Score algorithm for peak detection"
    function SmoothedZscoreAlgo(y, lag, threshold, influence)
    # Julia implimentation of http://stackoverflow.com/a/22640362/6029703
        n = length(y)
        signals = zeros(n) # init signal results
        filteredY = copy(y) # init filtered series
        avgFilter = zeros(n) # init average filter
        stdFilter = zeros(n) # init std filter
        avgFilter[lag - 1] = Statistics.mean(y[1:lag]) # init first value
        stdFilter[lag - 1] = Statistics.std(y[1:lag]) # init first value

        for i in range(lag, stop=n-1)
            if abs(y[i] - avgFilter[i-1]) > threshold*stdFilter[i-1]
                if y[i] > avgFilter[i-1]
                    signals[i] += 1 # postive signal
                else
                    signals[i] += -1 # negative signal
                end
                # Make influence lower
                filteredY[i] = influence*y[i] + (1-influence)*filteredY[i-1]
            else
            signals[i] = 0
            filteredY[i] = y[i]
            end
            avgFilter[i] = mean(filteredY[i-lag+1:i])
            stdFilter[i] = std(filteredY[i-lag+1:i])
        end
        return (signals = signals, avgFilter = avgFilter, stdFilter = stdFilter)
    end
    "This function shifts the peaks of the synchronization signal in order to calculate the Doppler effect coefficient"
    function ShiftSignal(ffts::Vector{Vector{Float64}}, self::cofdm, freq::Frequencies{Float64})
        doppler = 0.0;
        is_correct = false
        true_shift = 0.0
        true_data = [Vector{Bool}() for _ in 1:length(ffts)]
        zscore = [SmoothedZscoreAlgo(ffts[i], 3, 4, 0.5)[1] for i in 1:length(ffts)]
        correl = [Plots.Plot() for _ in 1:length(ffts)]
        for i in 1:length(ffts)
            correl[i] = Plots.plot(freq, [zscore[i].*1000, ffts[i]], title="Bit $i", labels=["zScore" "FFT"])
        end

        for shift in 0.84:0.00001:0.86   # Value for maximum shift (8Km/s)
            if shift == 0.85
                println('a')
            end
            error_number = 0
            for bit in 1:length(zscore) # For each sequence of cofdm
                empty!(true_data[bit])
                shifted_fft = Vector{Float64}()
                empty!(shifted_fft)
                if bit == 3
                    println("wow")
                end
                for k in 1:length(ffts[bit])    # FFT frequency shift loop
                    if trunc(Int,k*shift) < length(ffts[bit])
                        if trunc(Int, k*shift) == 0
                            push!(shifted_fft, zscore[bit][1])
                        else
                            push!(shifted_fft, zscore[bit][trunc(Int, k*shift)])
                        end
                    else
                        push!(shifted_fft, 0.0)
                    end
                end
                # Is the is shift correct for this bit
                bits = reconstruct_data(self, shifted_fft, freq)
                if isempty(bits) || length(bits) != length(self.data[bit])+length(self.sync_data[bit]) # If is empty or is not the right size
                    is_correct = false
                    for i in 1:bit  # Clear the false data
                        empty!(true_data[i])
                    end
                    printstyled("Wrong size ($(length(bits))) with bit $bit at doppler $shift\n", color=:red)
                    break
                elseif check_difference(self, bits, bit) <= length(self.sync_data[bit]) # If there is some sync data with the shift
                    printstyled("Bit $(bit) is shifted with $shift\n", color=:green)
                    printstyled("The number of wrong bits in sync channels is $(check_difference(self,bits, bit)).\n", color=:red)
                    error_number += check_difference(self,bits, bit)
                    if error_number > self.error_tolerance
                        is_correct = false
                        for i in 1:bit
                            empty(true_data[i])
                        end
                        printstyled("Error overflow for Doppler = $shift at $bit\n$error_number errors.\n")
                        break
                    else
                        is_correct = true
                        true_shift = shift
                        append!(true_data[bit], bits)
                    end
                else
                    is_correct = false
                    for i in 1:bit  # Clear the false data
                        empty!(true_data[i])
                    end
                    printstyled("$shift - $bit\n", color=:red)
                    break
                end
                if length(true_data[1]) > 11
                    println("davai")
                end
                empty!(bits)
            end
            if is_correct
                printstyled("The total number of errors in cyn channels is $error_number.\nOver a max of $(self.error_tolerance)", color=:green)
                break
            end
        end
        return true_shift, true_data, correl
    end
    "This function returns the number of different bit between two packets"
    function check_difference(sync::cofdm, data::Vector{Bool}, seq::Int64)
        error = 0
        for i in 1:length(sync.sync_channels)
            if sync.sync_data[seq][i] != data[sync.sync_channels[i]]
                error +=1
            end
        end
        return error
    end
    "Overload for total difference of the whole transmission"
    function check_difference(self::cofdm, data::Vector{Vector{Bool}})
        error = 0
        sync_count = 0
        done = false
            for i in 1:length(self.data) # For each sequence
                sync_count = 0
                for j in 1:length(self.data[i]) + length(self.sync_channels) # For each channel
                    done = false
                    for k in 1:length(self.sync_channels)
                        if j == self.sync_channels[k] && !done
                            if self.sync_data[i][k] != data[i][j]
                                error +=1
                            end
                            sync_count += 1
                            done = true
                            break
                        end
                    end
                    if !done
                        if self.data[i][j-sync_count] != data[i][j]
                            error += 1
                        end
                    end
                end
            end
        return error
    end
end
