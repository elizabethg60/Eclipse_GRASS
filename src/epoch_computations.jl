function compute_rv(lats::T, lons::T, epoch, obs_long, obs_lat, alt, band, index; moon_r::Float64=moon_radius) where T 
    """
    compute rv for a given grid size and timestamp - serial 
    
    lats: grid latitude size
    lons: grid longitude size
    epoch: timestamp
    obs_long: observer longtiude
    obs_lat: observer latitude
    alt: observer altitude
    """
#query JPL horizons for E, S, M position (km) and velocities (km/s)
    earth_pv = spkssb(399,epoch,"J2000")[1:3] 
    sun_pv = spkssb(10,epoch,"J2000")[1:3] 
    moon_pv = spkssb(301,epoch,"J2000")[1:3] 


#determine required position vectors
    #determine xyz stellar coordinates for lat/long grid
    SP_sun = get_xyz_for_surface(sun_radius, num_lats = lats, num_lons = lons)
    #transform xyz stellar coordinates of grid from sun frame to ICRF
    SP_bary = Matrix{Vector{Float64}}(undef,size(SP_sun)...)
    frame_transfer(pxform("IAU_SUN", "J2000", epoch), SP_sun, SP_bary)

    #determine xyz earth coordinates for lat/long of observatory
    EO_earth = pgrrec("EARTH", deg2rad(obs_long), deg2rad(obs_lat), alt, earth_radius, (earth_radius - earth_radius_pole) / earth_radius)
    #transform xyz earth coordinates of observatory from earth frame to ICRF
    EO_bary = pxform("IAU_EARTH", "J2000", epoch)*EO_earth

    #get vector from barycenter to observatory on Earth's surface
    BO_bary = earth_pv .+ EO_bary
    #get vector from observer to Sun's center 
    OS_bary = sun_pv - BO_bary 
    #get vector from observatory on earth's surface to moon center
    OM_bary = BO_bary .- moon_pv
    #get vector from barycenter to each patch on Sun's surface
    BP_bary = Matrix{Vector{Float64}}(undef,size(SP_bary)...)
    for i in eachindex(BP_bary)
        BP_bary[i] = sun_pv + SP_bary[i]
    end
    #vectors from observatory on Earth's surface to each patch on Sun's surface
    OP_bary = Matrix{Vector{Float64}}(undef,size(SP_bary)...)
    earth2patch_vectors(BP_bary, BO_bary, OP_bary)	

#calculate mu for each patch
    mu_grid = Matrix{Float64}(undef,size(SP_bary)...)
    calc_mu_grid!(SP_bary, OP_bary, mu_grid) 


#determine velocity vectors
    #determine velocity scalar for each patch 
    lat_grid = lat_grid_fc(size(SP_bary)...)
    v_scalar_grid = Matrix{Float64}(undef,size(SP_bary)...)
    v_scalar!(lat_grid, v_scalar_grid)
    #convert v_scalar to from km/day m/s
    for i in eachindex(v_scalar_grid)
        v_scalar_grid[i] = v_scalar_grid[i]/86.4  #+ (2*π*earth_radius*cos(deg2rad(obs_lat)))/86.4
    end

    #determine pole vector for each patch
    pole_vector_grid = Matrix{Vector{Float64}}(undef,size(SP_sun)...)
    pole_vector_grid!(SP_sun, pole_vector_grid) 

    #get velocity vector direction and set magnitude
    velocity_vector_solar = Matrix{Vector{Float64}}(undef,size(pole_vector_grid)...)
    v_vector(SP_sun, pole_vector_grid, v_scalar_grid, velocity_vector_solar)

    #transform into ICRF frame 
    velocity_vector_ICRF = Matrix{Vector{Float64}}(undef,size(velocity_vector_solar)...)
    frame_transfer(sxform("IAU_SUN", "J2000", epoch), velocity_vector_solar, velocity_vector_ICRF)

    #get projected velocity for each patch
    convective_velocities = convective_blueshift_interpol.(mu_grid)
    projected_velocities_no_cb = Matrix{Float64}(undef,size(SP_bary)...)
    projected_velocities_cb = Matrix{Float64}(undef,size(SP_bary)...)
    projected!(velocity_vector_ICRF, OP_bary, projected_velocities_no_cb, projected_velocities_cb, convective_velocities, epoch, BP_bary) 

    
#determine patches that are blocked by moon 
    #calculate the distance between tile corner and moon
    distance = map(x -> calc_proj_dist2(x, OM_bary), OP_bary)

    #calculate limb darkening weight for each patch 
    if band == "NIR"
        LD_all = map(x -> quad_limb_darkening_NIR(x), mu_grid)
    end

    zenith_angle_matrix = rad2deg.(map(x -> calc_proj_dist2(x, EO_bary), OP_bary))
    if band == "optical"
        LD_all = map((x,y) -> quad_limb_darkening_optical(x, index, y), mu_grid, zenith_angle_matrix)
    end

    #get indices for visible patches
    idx1 = mu_grid .> 0.0
    idx3 = (idx1) .& (distance .> atan((moon_r)/norm(OM_bary))) 
    
    #calculating the area element dA for each tile
    ϕ = range(deg2rad(-90.0), deg2rad(90.0), length=lats)
    θ =range(0.0, deg2rad(360.0), length=lons)
    dA_sub = map(x -> calc_dA(sun_radius, x, step(ϕ), step(θ)), lat_grid)
    #get total projected, visible area of larger tile
    dp_sub = map((x,y) -> abs(dot(x,y)), OP_bary, SP_bary) / norm(OS_bary)
    dA_total_proj = dA_sub .* dp_sub 

    #if no patches are visible, set mu, LD, projected velocity to zero 
    for i in 1:length(idx3)
        if idx3[i] == false
            LD_all[i] = 0.0
            projected_velocities_no_cb[i] = NaN
        end
    end
#determine mean intensity 
    mean_intensity = sum(view(LD_all .* dA_total_proj, idx3)) / sum(view(dA_total_proj, idx1))  

#determine mean weighted velocity from sun given blocking from moon 
    mean_weight_v_no_cb = sum(view((mean_intensity .* dA_total_proj) .* projected_velocities_no_cb, idx3)) / sum(view(mean_intensity .* dA_total_proj, idx3))
    mean_weight_v_cb =  sum(view(LD_all .* dA_total_proj .* projected_velocities_cb, idx3)) / sum(view(LD_all .* dA_total_proj, idx3))

#get ra and dec of solar grid patches
    OP_ra_dec = SPICE.recrad.(OP_bary)

    return mean_weight_v_no_cb, mean_weight_v_cb, mean_intensity, rad2deg.(getindex.(OP_ra_dec,2)), rad2deg.(getindex.(OP_ra_dec,3)), projected_velocities_no_cb, projected_velocities_cb
end