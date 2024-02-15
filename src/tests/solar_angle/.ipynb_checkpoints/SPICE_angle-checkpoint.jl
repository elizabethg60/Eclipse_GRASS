using SPICE
using LinearAlgebra 

UTC_time = ["2023-10-14T15:00:00", "2023-10-14T15:05:00", "2023-10-14T15:10:00", "2023-10-14T15:15:00", "2023-10-14T15:20:00","2023-10-14T15:25:00", "2023-10-14T15:30:00", "2023-10-14T15:35:00", "2023-10-14T15:40:00", "2023-10-14T15:45:00","2023-10-14T15:50:00","2023-10-14T15:55:00","2023-10-14T16:00:00", "2023-10-14T16:05:00", "2023-10-14T16:10:00", "2023-10-14T16:15:00", "2023-10-14T16:20:00","2023-10-14T16:25:00", "2023-10-14T16:30:00", "2023-10-14T16:35:00", "2023-10-14T16:40:00", "2023-10-14T16:45:00","2023-10-14T16:50:00","2023-10-14T16:55:00", "2023-10-14T17:00:00", "2023-10-14T17:05:00", "2023-10-14T17:10:00", "2023-10-14T17:15:00", "2023-10-14T17:20:00","2023-10-14T17:25:00", "2023-10-14T17:30:00", "2023-10-14T17:35:00", "2023-10-14T17:40:00", "2023-10-14T17:45:00","2023-10-14T17:50:00","2023-10-14T17:55:00", "2023-10-14T18:00:00", "2023-10-14T18:05:00"]
epoch = utc2et.(UTC_time)

sun_radius = bodvrd("SUN","RADII")[1]
sun_axis_sun = [0.0,0.0,sun_radius]
#sun_axis_bary = Vector{Float64}[]

angle_list = Vector{Float64}[]
for i in 1:length(epoch)
    #push!(sun_axis_bary, pxform("IAU_SUN", "J2000", epoch[i])*sun_axis_sun)
    sun_axis_earth = pxform("IAU_SUN", "IAU_Earth", epoch[i])*sun_axis_sun
    bary_axis_earth = pxform("J2000", "IAU_Earth", epoch[i])*sun_axis_sun
    
    angle = acos(dot(sun_axis_earth, bary_axis_earth) / (norm(sun_axis_earth) * norm(bary_axis_earth)))
    push!(angle_list, [rad2deg(angle)])
end

print(angle_list)