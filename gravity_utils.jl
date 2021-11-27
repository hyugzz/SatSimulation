#using Core, Base, Distances, SatelliteToolbox, Plot

mutable struct space_obj
    a::Vector
    v::Vector
    p::Vector
    mass::Float64
end

function set_acc(self::space_obj, forces::Array{Vector})
    for force in forces
        self.a += force
    end
end

function update_spe(self::space_obj, delta::Float64)
    self.v += self.a.*delta;
end

function update_pos(self::space_obj, delta::Float64)
    update_spe(self, delta);
    return self.p += self.v.*delta
end


function calculate_gravity(space_obj1::space_obj, planets::Array{space_obj})
    tmp = [0.0,0.0,0.0]
    for item in planets
        d = evaluate(Euclidean(), space_obj1.p, item.p)
        d = Euclidean()(space_obj1.p, item.p);
        tmp += ((6.67*10^-11) * (item.mass - space_obj1.mass) / d^2) .* (space_obj1.p - item.p)
    end
    return tmp;
end

function calculate_gravity(space_obj1::space_obj, space_obj2::space_obj)
    d = evaluate(Euclidean(), space_obj1.p, space_obj2.p)
    d = Euclidean()(space_obj1.p, space_obj2.p);
    return ((6.67*10^-11) * (space_obj2.mass - space_obj1.mass) / d^2) .* (space_obj1.p - space_obj2.p)
end

function step!(ojb::space_obj)
    delta = 1000.0
    set_acc(a, Vector[calculate_gravity(a,b)]);
    update_pos(a, delta);
end

function DisassembleVector(vectorList)
    x = Float64[]
    y = Float64[]
    z = Float64[]
    for vec in vectorList
        push!(x, vec[1])
        push!(y, vec[2])
        push!(z, vec[3])
    end
    return x,y,z
end

function ShowOrbit(vectorList)
    x,y,z = DisassembleVector(vectorList)
    plot(x,y,z)
end
