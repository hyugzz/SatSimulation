include("gravity_utils.jl")
using Core,Base,Distances, Plots, SatelliteToolbox


tles = read_tle("starlink.txt")
tles[1]

#SGP4
orbit = init_orbit_propagator(Val(:sgp4), tles[1])
#@gif for i=1:500
r, v = propagate!(orbit, collect(0:1:150)*60)
x = Float64[]
y = Float64[]
z = Float64[]

#x,y,z = DisassembleVector(r)

plt = plot3d(r, xlim=(-6000000,6000000),ylim=(-6000000,6000000),zlim=(-6000000,6000000))

fig = figure()
ax = fig.gca(projection="3d")

ax.quiver(r[:][1],r[:][2],r[:][3], v[:][1],v[:][2],v[:][3])
plot(ax)
