
mutable struct cofdm
    data::Array{Bool, 1}
    interval::Float64
end

function  define_data(self::cofdm, in::Array{Bool,1})
    self.data = in;
end

function generate_signal(self::cofdm)
    #Generate sine wave
    signal = Float64[]
    for bit in data
        if bit==1
            tmp = create_sine(1)
            append!(signal, tmp)
            append!(data2, bit)
        elseif bit==0
            tmp = create_sine(0.5)
            append!(signal, tmp)
            append!(data2, bit)
        end
    end
    return signal
end
 
