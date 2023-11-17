#parallel v2 multi-processing 
function lat_grid_fc_pa2(num_lats, num_lon)
    ϕ = deg2rad.(range(-90.0, 90.0, length=num_lats))
    return repeat(ϕ,num_lon)
end 

function compute_solar_grid(lats::T, lons::T, epoch) where T 
    #determine xyz stellar coordinates for lat/long grid
    return get_xyz_for_surface(sun_radius, num_lats = lats, num_lons = lons)
end

function compute_rv_pa2(lats::T, lons::T, epoch, obs_long, obs_lat, alt, SP_sun, lat_index) where T
#query JPL horizons for E, S, M position (km) and velocities (km/s)
    earth_pv = spkssb(399,epoch,"J2000") 
    sun_pv = spkssb(10,epoch,"J2000")
    moon_pv = spkssb(301,epoch,"J2000")


#determine required position vectors
    #transform xyz stellar coordinates of grid from sun frame to ICRF
    SP_bary = pxform("IAU_SUN", "J2000", epoch)*SP_sun

    #determine vectors from Earth observatory surface to each patch
    #determine xyz earth coordinates for lat/long 
    EO_earth = pgrrec("EARTH", deg2rad(obs_long), deg2rad(obs_lat), alt, earth_radius, (earth_radius - earth_radius_pole) / earth_radius)
    #transform xyz earth coordinates of observatory from earth frame to ICRF
    EO_bary = pxform("IAU_EARTH", "J2000", epoch)*EO_earth
    #get vector from barycenter to observatory on Earth's surface
    BO_bary = earth_pv[1:3] .+ EO_bary

    #get vector from observatory on earth's surface to moon center
    OM_bary = moon_pv[1:3] .- BO_bary
    #get vector from barycenter to each patch on Sun's surface
    BP_bary = sun_pv[1:3] .+ SP_bary

    #get vector from observatory on Earth's surface to Sun's center
    SO_bary = sun_pv[1:3] .- BO_bary
     
    #vectors from observatory on Earth's surface to each patch on Sun's surface
    OP_bary = BP_bary .- BO_bary


#calculate mu for each patch
    mu_grid = calc_mu(SP_bary, OP_bary)


#determine velocity vectors
    #determine velocity scalar for each patch 
    lat_grid = lat_grid_fc_pa2(lats, lons)
    v_scalar_grid = (2*π*sun_radius*cos(lat_grid[lat_index]))/(rotation_period(lat_grid[lat_index]))
    #convert v_scalar to from km/day m/s
    v_scalar_grid = v_scalar_grid/86.4
    #determine velocity vector + projected velocity for each patch 
    #determine pole vector for each patch
    pole_vector_grid = SP_sun - [0.0, 0.0, SP_sun[3]]

    #get velocity vector direction and set magnitude
    velocity_vector_solar = [SP_sun;(cross(pole_vector_grid, [0.0,0.0,sun_radius]) / norm(cross(pole_vector_grid, [0.0,0.0,sun_radius]))).*v_scalar_grid]
    
    #transform into ICRF frame 
    velocity_vector_ICRF = sxform("IAU_SUN", "J2000", epoch) * velocity_vector_solar

    #get projected velocity for each patch
    vel = [velocity_vector_ICRF[4],velocity_vector_ICRF[5],velocity_vector_ICRF[6]]
    angle = dot(OP_bary, vel) / (norm(OP_bary) * norm(vel))
    projected_velocities = -(norm(vel) * angle)

    
#determine patches that are blocked by moon 
    #calculate the distance between tile corner and moon
    distance = calc_proj_dist2(OP_bary, OM_bary)

    #calculate limb darkening weight for each patch 
    LD_all = quad_limb_darkening_optical(mu_grid, 0.4, 0.26)

    #get indices for visible patches
    idx1 = mu_grid > 0.0
    idx3 = (idx1) & (distance > atan(moon_radius/norm(OM_bary))^2)

    #calculating the area element dA for each tile
    ϕ = range(deg2rad(-90.0), deg2rad(90.0), length=lats)
    θ =range(0.0, deg2rad(360.0), length=lons)   
    dA_sub = calc_dA(sun_radius, lat_grid[lat_index], step(ϕ), step(θ))
    #get total projected, visible area of larger tile
    dp_sub = abs(dot(SP_bary, OP_bary))
    dA_total_proj = dA_sub .* dp_sub
    #if no patches are visible, set mu, LD, projected velocity to zero 
    if idx3 == false
        mu_grid = NaN
        LD_all = 0.0
        projected_velocities = NaN
    end

    return LD_all, projected_velocities, idx1, dA_total_proj
end

function parallel_v2()
    println("# Now Julia is using ", nworkers(), " workers.")
    time_stamps = utc2et("2015-03-20T09:42:00") 
    obs_lat = 51.54548 
    obs_long = 9.905548
    alt = 0.15
    #strong scaling
    num_workers_all = nworkers()
    wall_time = zeros(num_workers_all)
    wall_time_mid = zeros(num_workers_all)
    sun_grid_small = compute_solar_grid(50, 100, time_stamps)
    sun_grid_mid = compute_solar_grid(250, 500, time_stamps)
    #week scaling 
    wall_time_weak = zeros(num_workers_all)
    for nw in num_workers_all:-1:1
        #strong scaling
        wall_time[nw] = @elapsed output = @distributed (vcat) for idx ∈ 1:length(sun_grid_small)
            compute_rv_pa2(50, 100, time_stamps, obs_long, obs_lat, alt, sun_grid_small[idx], idx) 
        end 
        wall_time_mid[nw] = @elapsed output = @distributed (vcat) for idx ∈ 1:length(sun_grid_mid)
            compute_rv_pa2(250, 500, time_stamps, obs_long, obs_lat, alt, sun_grid_mid[idx], idx) 
        end 
        #weak scaling
        sun_grid = compute_solar_grid(nw*25, nw*25*2, time_stamps)
        wall_time_weak[nw] = @elapsed output = @distributed (vcat) for idx ∈ 1:length(sun_grid)
            compute_rv_pa2(nw*25, nw*25*2, time_stamps, obs_long, obs_lat, alt, sun_grid[idx], idx) 
        end
        if nw > 1
            rmprocs(last(workers()))            # Remove one worker
        end
        println("# Now Julia is using ", nworkers(), " workers.")
    end

    @save "src/test/scaling_v2.jld2"
    jldopen("src/test/scaling_v2.jld2", "a+") do file
        file["num_workers_all"] = num_workers_all 
        file["wall_time"] = wall_time
        file["wall_time_mid"] = wall_time_mid
        file["wall_time_weak"] = wall_time_weak
    end
end