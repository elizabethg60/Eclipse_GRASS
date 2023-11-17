function compute_rv_pa(lats::T, lons::T, epoch, obs_long, obs_lat, alt, band; moon_r::Float64=moon_radius) where T 
#query JPL horizons for E, S, M position (km) and velocities (km/s)
    earth_pv = spkssb(399,epoch,"J2000") 
    sun_pv = spkssb(10,epoch,"J2000")
    moon_pv = spkssb(301,epoch,"J2000")


#determine required position vectors
    #determine xyz stellar coordinates for lat/long grid    
    SP_sun = get_xyz_for_surface(sun_radius, num_lats = lats, num_lons = lons)                                   
    #transform xyz stellar coordinates of grid from sun frame to ICRF
    SP_bary = Matrix{Vector{Float64}}(undef,size(SP_sun)...)
    start_time = time()
    frame_transfer_pa(pxform("IAU_SUN", "J2000", epoch), SP_sun, SP_bary)

    #determine xyz earth coordinates for lat/long of observatory
    EO_earth = pgrrec("EARTH", deg2rad(obs_long), deg2rad(obs_lat), alt, earth_radius, (earth_radius - earth_radius_pole) / earth_radius)
    #transform xyz earth coordinates of observatory from earth frame to ICRF
    EO_bary = pxform("IAU_EARTH", "J2000", epoch)*EO_earth

    #get vector from barycenter to observatory on Earth's surface
    BO_bary = earth_pv[1:3] .+ EO_bary
    #get vector from observatory on earth's surface to moon center
    OM_bary = moon_pv[1:3] .- BO_bary
    #get vector from barycenter to each patch on Sun's surface
    BP_bary = Matrix{Vector{Float64}}(undef,size(SP_bary)...)
    Threads.@threads for i in eachindex(BP_bary)
        BP_bary[i] = sun_pv[1:3] + SP_bary[i]
    end
    #get vector from observatory on Earth's surface to Sun's center
    SO_bary = sun_pv[1:3] .- BO_bary  
    #vectors from observatory on Earth's surface to each patch on Sun's surface
    OP_bary = Matrix{Vector{Float64}}(undef,size(SP_bary)...)
    earth2patch_vectors_pa(BP_bary, BO_bary, OP_bary)	
        
 
#calculate mu for each patch
    mu_grid = Matrix{Float64}(undef,size(SP_bary)...)
    calc_mu_grid_pa!(SP_bary, OP_bary, mu_grid)
    

#determine velocity vectors
    #determine velocity scalar for each patch 
    lat_grid = lat_grid_fc_pa(size(SP_bary)...)                                                              
    v_scalar_grid = Matrix{Float64}(undef,size(SP_bary)...)
    v_scalar_pa!(lat_grid, v_scalar_grid)
    #convert v_scalar to from km/day m/s
    Threads.@threads for i in eachindex(v_scalar_grid)
        v_scalar_grid[i] = v_scalar_grid[i]/86.4
    end
    
    #determine pole vector for each patch
    pole_vector_grid = Matrix{Vector{Float64}}(undef,size(SP_sun)...)
    pole_vector_grid_pa!(SP_sun, [0.0,0.0,sun_radius], pole_vector_grid)

    #get velocity vector direction and set magnitude
    velocity_vector_solar = Matrix{Vector{Float64}}(undef,size(pole_vector_grid)...)
    v_vector_pa(SP_sun, pole_vector_grid, v_scalar_grid, velocity_vector_solar)
    
    #transform into ICRF frame 
    velocity_vector_ICRF = Matrix{Vector{Float64}}(undef,size(velocity_vector_solar)...)
    frame_transfer_pa(sxform("IAU_SUN", "J2000", epoch), velocity_vector_solar, velocity_vector_ICRF)
    
    #get projected velocity for each patch
    projected_velocities = Matrix{Float64}(undef,size(SP_bary)...)
    projected_pa!(velocity_vector_ICRF, OP_bary, projected_velocities)
    
 
#determine patches that are blocked by moon 
    #calculate the distance between tile corner and moon
    distance = ThreadsX.map(x -> calc_proj_dist2(x, OM_bary), OP_bary)
    
    #calculate limb darkening weight for each patch 
    if band == "NIR"
        LD_all = map(x -> quad_limb_darkening_NIR(x, 0.4, 0.26), mu_grid)
    end
    if band == "optical"
        LD_all = map(x -> quad_limb_darkening_optical(x, 0.4, 0.26), mu_grid)
    end
    
    #get indices for visible patches                                                                    
    idx1 = mu_grid .> 0.0
    idx3 = (idx1) .& (distance .> atan((moon_r)/norm(OM_bary))^2) 

    #calculating the area element dA for each tile
    ϕ = range(deg2rad(-90.0), deg2rad(90.0), length=lats)
    θ =range(0.0, deg2rad(360.0), length=lons)
    dA_sub = ThreadsX.map(x -> calc_dA(sun_radius, x, step(ϕ), step(θ)), lat_grid)
    #get total projected, visible area of larger tile
    dp_sub = ThreadsX.map((x,y) -> abs(dot(x,y)), SP_bary, OP_bary) 
    dA_total_proj = dA_sub .* dp_sub
    
     #if no patches are visible, set mu, LD, projected velocity to zero 
    Threads.@threads for i in 1:length(idx3)
        if idx3[i] == false
            mu_grid[i] = NaN
            LD_all[i] = 0.0
            projected_velocities[i] = NaN
        end
    end
    

#determine mean weighted velocity from sun given blocking from moon 
    mean_weight_v = NaNMath.sum(matrix_multi(LD_all .* dA_total_proj, projected_velocities)) / NaNMath.sum(LD_all .* dA_total_proj)                     
#determine mean intensity 
    mean_intensity = sum(view(LD_all .* dA_total_proj, idx1)) / length(view(LD_all .* dA_total_proj, idx1))
    return mean_weight_v, mean_intensity, time() - start_time
end