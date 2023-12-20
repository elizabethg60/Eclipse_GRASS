# using Pkg; Pkg.activate(".")
using GRASS
using Random
using Revise
using Statistics
using BenchmarkTools
using SPICE

# plotting
#using LaTeXStrings
#import PyPlot; plt = PyPlot; mpl = plt.matplotlib; plt.ioff()
#mpl.style.use(GRASS.moddir * "fig.mplstyle")

# set up paramaters for spectrum
N = 75
Nt = 1
GRASS.get_kernels()
disk = GRASS.DiskParams(N=N, Nt=Nt, Nsubgrid=5)

reiners_timestamps = ["2015-03-20T7:07:57.2", "2015-03-20T7:09:45.8", "2015-03-20T7:11:34.3", "2015-03-20T7:13:22.9", "2015-03-20T7:15:11.5", "2015-03-20T7:17:00.3", "2015-03-20T7:18:49.7", "2015-03-20T7:20:38.5", "2015-03-20T7:22:27.4", "2015-03-20T7:24:16.5", "2015-03-20T7:26:05.5", "2015-03-20T7:27:53.8", "2015-03-20T7:29:42.5", "2015-03-20T7:31:30.7", "2015-03-20T7:33:19.2", "2015-03-20T7:35:09.5", "2015-03-20T7:36:58.0", "2015-03-20T7:38:46.1", "2015-03-20T7:40:34.6", "2015-03-20T7:42:22.8", "2015-03-20T7:44:11.5", "2015-03-20T7:46:00.1", "2015-03-20T7:47:48.8", "2015-03-20T7:49:37.1", "2015-03-20T7:51:25.3", "2015-03-20T7:53:13.9", "2015-03-20T7:55:02.3", "2015-03-20T7:56:50.7", "2015-03-20T7:58:39.0", "2015-03-20T8:00:27.2", "2015-03-20T8:02:15.7", "2015-03-20T8:04:04.0", "2015-03-20T8:05:53.0", "2015-03-20T8:07:41.4", "2015-03-20T8:09:30.0", "2015-03-20T8:11:18.5", "2015-03-20T8:13:07.1", "2015-03-20T8:14:55.3", "2015-03-20T8:16:43.7", "2015-03-20T8:18:32.1", "2015-03-20T8:20:20.4", "2015-03-20T8:22:09.4", "2015-03-20T8:23:57.7", "2015-03-20T8:25:46.1", "2015-03-20T8:27:34.7", "2015-03-20T8:29:23.0", "2015-03-20T8:31:11.3", "2015-03-20T8:32:59.5", "2015-03-20T8:34:47.9", "2015-03-20T8:36:36.3", "2015-03-20T8:38:54.5", "2015-03-20T8:40:43.0", "2015-03-20T8:42:31.4", "2015-03-20T8:44:19.9", "2015-03-20T8:46:08.5", "2015-03-20T8:47:56.8", "2015-03-20T8:49:45.5", "2015-03-20T8:51:34.0", "2015-03-20T8:53:22.9", "2015-03-20T8:55:11.5", "2015-03-20T8:56:59.6", "2015-03-20T8:58:47.9", "2015-03-20T9:00:36.1", "2015-03-20T9:02:24.8", "2015-03-20T9:04:13.1", "2015-03-20T9:06:01.7", "2015-03-20T9:07:50.2", "2015-03-20T9:09:38.7", "2015-03-20T9:11:27.1", "2015-03-20T9:13:15.5", "2015-03-20T9:15:04.0", "2015-03-20T9:16:53.1", "2015-03-20T9:18:42.1", "2015-03-20T9:20:30.9", "2015-03-20T9:22:19.7", "2015-03-20T9:24:08.5", "2015-03-20T9:25:57.3", "2015-03-20T9:27:45.9", "2015-03-20T9:29:34.7", "2015-03-20T9:31:23.1", "2015-03-20T9:33:11.5", "2015-03-20T9:34:59.9", "2015-03-20T9:36:48.4", "2015-03-20T9:38:37.4", "2015-03-20T9:40:26.3", "2015-03-20T9:42:15.0", "2015-03-20T9:44:03.8", "2015-03-20T9:45:52.2", "2015-03-20T9:47:40.7", "2015-03-20T9:49:29.7", "2015-03-20T9:51:18.4", "2015-03-20T9:53:07.5", "2015-03-20T9:54:55.9", "2015-03-20T9:56:44.9", "2015-03-20T9:58:33.3", "2015-03-20T10:00:21.8", "2015-03-20T10:02:10.1", "2015-03-20T10:03:58.6", "2015-03-20T10:05:47.4", "2015-03-20T10:07:36.2", "2015-03-20T10:09:54.5", "2015-03-20T10:11:43.9", "2015-03-20T10:13:33.6", "2015-03-20T10:15:22.6", "2015-03-20T10:17:11.7", "2015-03-20T10:19:00.9", "2015-03-20T10:20:49.9", "2015-03-20T10:22:38.9", "2015-03-20T10:24:27.9", "2015-03-20T10:26:17.0", "2015-03-20T10:28:07.1", "2015-03-20T10:29:56.1", "2015-03-20T10:31:45.1", "2015-03-20T10:33:34.0", "2015-03-20T10:35:22.9", "2015-03-20T10:37:12.0", "2015-03-20T10:39:01.0", "2015-03-20T10:40:49.9", "2015-03-20T10:42:38.7", "2015-03-20T10:44:27.6", "2015-03-20T10:46:16.8", "2015-03-20T10:48:05.8", "2015-03-20T10:49:54.9", "2015-03-20T10:51:43.6", "2015-03-20T10:53:32.7", "2015-03-20T10:55:21.9", "2015-03-20T10:57:10.8", "2015-03-20T10:58:59.7", "2015-03-20T11:00:49.8", "2015-03-20T11:02:38.6", "2015-03-20T11:04:27.7", "2015-03-20T11:06:16.8", "2015-03-20T11:08:05.7", "2015-03-20T11:09:54.6", "2015-03-20T11:11:43.6", "2015-03-20T11:13:33.1", "2015-03-20T11:15:22.0", "2015-03-20T11:17:10.9", "2015-03-20T11:18:59.9", "2015-03-20T11:20:48.7", "2015-03-20T11:22:37.7", "2015-03-20T11:24:26.7", "2015-03-20T11:26:15.7", "2015-03-20T11:28:04.3", "2015-03-20T11:29:53.0", "2015-03-20T11:31:41.7", "2015-03-20T11:33:30.3", "2015-03-20T11:35:19.0", "2015-03-20T11:37:07.7", "2015-03-20T11:38:56.9", "2015-03-20T11:48:37.5", "2015-03-20T11:50:26.6", "2015-03-20T11:52:15.6", "2015-03-20T11:54:04.2", "2015-03-20T11:55:53.1", "2015-03-20T11:57:41.8", "2015-03-20T11:59:30.5", "2015-03-20T12:01:19.1", "2015-03-20T12:03:07.6"]
#convert from utc to et as needed by SPICE
time_stamps = utc2et.(reiners_timestamps)

#Gottingen location
obs_lat = 51.560583 
obs_long = 9.944333
alt = 0.201

for i in 1:length([1])#time_stamps)
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
    GRASS.eclipse_compute_quantities!(disk, time_stamps[i], obs_long, obs_lat, alt, ϕc, θc, μs, ld, dA, xyz, wts, z_rot, ax_codes)
end

# # disk coordinates
# ϕe = disk.ϕe
# ϕc = disk.ϕc
# θe = disk.θe
# θc = disk.θc
# R_x = disk.R_x
# R_y = disk.R_y
# R_z = disk.R_z

# # get color scalar mappable
# dat = wts ./ maximum(wts)
# cmap = plt.cm.afmhot

# # dat = z_rot .* 3e8
# # cmap = plt.cm.seismic

# norm = mpl.colors.Normalize(vmin=minimum(dat), vmax=maximum(dat))
# smap = plt.cm.ScalarMappable(cmap=cmap, norm=norm)

# # initialize figure
# fig, ax = plt.subplots(1,1, figsize=(8,8))

# # loop over grid positions
# println("\t>>> Plotting!")
# for i in 1:length(ϕe)-1
#     lat = range(ϕe[i], ϕe[i+1], length=4)
#     for j in 1:disk.Nθ[i]
#         lon = range(θe[i,j], θe[i,j+1], length=4)

#         border = (([lat[1], lat[2], lat[3], lat[4], lat[4], lat[4], lat[4], lat[3], lat[2], lat[1], lat[1], lat[1]]),
#                   ([lon[1], lon[1], lon[1], lon[1], lon[2], lon[3], lon[4], lon[4], lon[4], lon[4], lon[3], lon[2]]))


#         out = GRASS.sphere_to_cart.(1.0, border...) #update coordiante system here and below rotation
#         x = getindex.(out, 1)
#         y = getindex.(out, 2)
#         z = getindex.(out, 3)

#         # rotate it
#         for k in eachindex(x)
#             x0 = x[k]
#             y0 = y[k]
#             z0 = z[k]

#             x[k] = x0 * R_x[1,1] + y0 * R_x[1,2] + z0 * R_x[1,3]
#             y[k] = x0 * R_x[2,1] + y0 * R_x[2,2] + z0 * R_x[2,3]
#             z[k] = x0 * R_x[3,1] + y0 * R_x[3,2] + z0 * R_x[3,3]
#         end

#         idx = z .>= 0
#         if any(idx)
#             # ax.plot(x[idx], y[idx], color="k", lw=1)
#             ax.fill(x[idx], y[idx], c=smap.to_rgba(dat[i,j]))
#         end
#     end
# end

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

# ax.set_xlabel(L"\Delta {\rm x\ [Stellar\ Radii]}")
# ax.set_ylabel(L"\Delta {\rm y\ [Stellar\ Radii]}")
# ax.set_aspect("equal")
# ax.grid(false)
# ax.invert_xaxis()
# cb = fig.colorbar(smap, ax=ax, fraction=0.1, shrink=0.8)
# cb.set_label(L"{\rm Weighted\ Relative\ Intensity}")
# plt.show()
# plt.clf(); plt.close();

