# using Pkg; Pkg.activate(".")
using GRASS
using Random
using Revise
using Statistics
using BenchmarkTools

# plotting
using LaTeXStrings
import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
mpl.style.use(GRASS.moddir * "fig.mplstyle")

# set up paramaters for spectrum
N = 75
Nt = 1
disk = DiskParams(N=N, Nt=Nt, inclination=45.0, Nsubgrid=5)

# allocate the memory for keys, velocities, ld, etc.
ϕc = zeros(size(disk.θc))
θc = zeros(size(disk.θc))
μs = zeros(size(disk.θc))
ld = zeros(size(disk.θc))
dA = zeros(size(disk.θc))
xyz = zeros(size(disk.θc)..., 3)
wts = zeros(size(disk.θc))
z_rot = zeros(size(disk.θc))
ax_codes = zeros(Int, size(disk.θc))

# precompute quantities
GRASS.eclipse_compute_quantities!(disk, ϕc, θc, μs, ld, dA, xyz, wts, z_rot, ax_codes)

# disk coordinates
ϕe = disk.ϕe
ϕc = disk.ϕc
θe = disk.θe
θc = disk.θc
R_x = disk.R_x
R_y = disk.R_y
R_z = disk.R_z

# get color scalar mappable
dat = wts ./ maximum(wts)
cmap = plt.cm.afmhot

# dat = z_rot .* 3e8
# cmap = plt.cm.seismic

norm = mpl.colors.Normalize(vmin=minimum(dat), vmax=maximum(dat))
smap = plt.cm.ScalarMappable(cmap=cmap, norm=norm)

# initialize figure
fig, ax = plt.subplots(1,1, figsize=(8,8))

# loop over grid positions
println("\t>>> Plotting!")
for i in 1:length(ϕe)-1
    lat = range(ϕe[i], ϕe[i+1], length=4)
    for j in 1:disk.Nθ[i]
        lon = range(θe[i,j], θe[i,j+1], length=4)

        border = (([lat[1], lat[2], lat[3], lat[4], lat[4], lat[4], lat[4], lat[3], lat[2], lat[1], lat[1], lat[1]]),
                  ([lon[1], lon[1], lon[1], lon[1], lon[2], lon[3], lon[4], lon[4], lon[4], lon[4], lon[3], lon[2]]))


        out = GRASS.sphere_to_cart.(1.0, border...) #update coordiante system here and below rotation
        x = getindex.(out, 1)
        y = getindex.(out, 2)
        z = getindex.(out, 3)

        # rotate it
        for k in eachindex(x)
            x0 = x[k]
            y0 = y[k]
            z0 = z[k]

            x[k] = x0 * R_x[1,1] + y0 * R_x[1,2] + z0 * R_x[1,3]
            y[k] = x0 * R_x[2,1] + y0 * R_x[2,2] + z0 * R_x[2,3]
            z[k] = x0 * R_x[3,1] + y0 * R_x[3,2] + z0 * R_x[3,3]
        end

        idx = z .>= 0
        if any(idx)
            # ax.plot(x[idx], y[idx], color="k", lw=1)
            ax.fill(x[idx], y[idx], c=smap.to_rgba(dat[i,j]))
        end
    end
end

# # get equator coords
# latitude = deg2rad(0.0)
# longitude = deg2rad.(range(0.0, 360.0, length=200))
# x_eq = []
# y_eq = []
# z_eq = []
# for i in eachindex(longitude)
#     out = GRASS.sphere_to_cart.(1.0, latitude, longitude[i])
#     x = getindex(out, 1)
#     y = getindex(out, 2)
#     z = getindex(out, 3)

#     x0 = x
#     y0 = y
#     z0 = z

#     x = x0 * R_x[1,1] + y0 * R_x[1,2] + z0 * R_x[1,3]
#     y = x0 * R_x[2,1] + y0 * R_x[2,2] + z0 * R_x[2,3]
#     z = x0 * R_x[3,1] + y0 * R_x[3,2] + z0 * R_x[3,3]

#     push!(x_eq, x)
#     push!(y_eq, y)
#     push!(z_eq, z)
# end

# # sort the values on increasing x
# idx_eq = sortperm(x_eq)
# x_eq = x_eq[idx_eq]
# y_eq = y_eq[idx_eq]
# z_eq = z_eq[idx_eq]

# idx_eq = z_eq .> 0.0

# # plot the equator
# ax.plot(x_eq[idx_eq], y_eq[idx_eq], color="white", ls="--", zorder=3, alpha=0.75)

# # get meridians
# latitude = deg2rad.(range(-89.0, 89.0, length=200))
# longitude = deg2rad.(range(0.0, 360.0, step=90.0))

# for j in eachindex(longitude)
#     out = GRASS.sphere_to_cart.(1.0, latitude, longitude[j])

#     out = hcat(out...)

#     x = out[1,:]
#     y = out[2,:]
#     z = out[3,:]

#     x0 = x
#     y0 = y
#     z0 = z

#     x = x0 .* R_x[1,1] .+ y0 .* R_x[1,2] .+ z0 .* R_x[1,3]
#     y = x0 .* R_x[2,1] .+ y0 .* R_x[2,2] .+ z0 .* R_x[2,3]
#     z = x0 .* R_x[3,1] .+ y0 .* R_x[3,2] .+ z0 .* R_x[3,3]

#     # plot the meridian
#     idx = z .> 0.0
#     ax.plot(x[idx], y[idx], color="white", ls="--", zorder=3, alpha=0.75)
# end

ax.set_xlabel(L"\Delta {\rm x\ [Stellar\ Radii]}")
ax.set_ylabel(L"\Delta {\rm y\ [Stellar\ Radii]}")
ax.set_aspect("equal")
ax.grid(false)
ax.invert_xaxis()
cb = fig.colorbar(smap, ax=ax, fraction=0.1, shrink=0.8)
cb.set_label(L"{\rm Weighted\ Relative\ Intensity}")
plt.show()
plt.clf(); plt.close();

