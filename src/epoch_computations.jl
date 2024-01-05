function compute_rv(lats::T, epoch, obs_long, obs_lat, alt, band, index; moon_r::Float64=moon_radius) where T 
    """
    compute rv for a given grid size and timestamp - serial 
    
    lats: grid latitude size
    lons: grid longitude size
    epoch: timestamp
    obs_long: observer longtiude
    obs_lat: observer latitude
    alt: observer altitude
    """

# disk gridding
    # get latitude grid edges and centers
    disk_ϕe = range(deg2rad(-90.0), deg2rad(90.0), length=lats+1)
    disk_ϕc = get_grid_centers(disk_ϕe)
    # number of longitudes in each latitude slice
    Nθ = get_Nθ.(disk_ϕc, step(disk_ϕe))
    # make longitude grid
    disk_θe = zeros(lats+1, maximum(Nθ)+1)
    disk_θc = zeros(lats, maximum(Nθ))
    for i in eachindex(Nθ)
        edges = range(deg2rad(0.0), deg2rad(360.0), length=Nθ[i]+1)
        disk_θc[i, 1:Nθ[i]] .= get_grid_centers(edges)
        disk_θe[i, 1:Nθ[i]+1] .= collect(edges)
    end


# collect geometry vectors 
    # query JPL horizons for E, S, M position (km) and velocities (km/s)
    BE_bary = spkssb(399,epoch,"J2000")
    BS_bary = spkssb(10,epoch,"J2000")
    BM_bary = spkssb(301,epoch,"J2000")

    #determine xyz earth coordinates for lat/long of observatory
    flat_coeff = (earth_radius - earth_radius_pole) / earth_radius
    EO_earth_pos = georec(deg2rad(obs_long), deg2rad(obs_lat), alt, earth_radius, flat_coeff)
    #set earth velocity vectors
    EO_earth = vcat(EO_earth_pos, [0.0, 0.0, 0.0])
    #transform into ICRF frame
    EO_bary = sxform("ITRF93", "J2000", epoch) * EO_earth

    # get vector from barycenter to observatory on Earth's surface
    BO_bary = BE_bary .+ EO_bary

    # set string for ltt and abberation
    lt_flag = "CN+S"

    # get light travel time corrected OS vector
    OS_bary, OS_lt, OS_dlt = spkltc(10, epoch, "J2000", lt_flag, BO_bary)
    # get light travel time corrected ES vector
    ES_bary, ES_lt, ES_dlt = spkltc(10, epoch, "J2000", lt_flag, BE_bary)

    # get vector from observatory on earth's surface to moon center
    OM_bary, OM_lt, OM_dlt = spkltc(301, epoch, "J2000", lt_flag, BO_bary)

    # get modified epch
    epoch_lt = epoch - OS_lt

    # get vector for sun pole
    sun_lat = deg2rad(90.0)
    sun_lon = deg2rad(0.0)
    sun_pole_sun = pgrrec("SUN", sun_lon, sun_lat, 0.0, sun_radius, 0.0)
    sun_pole_bary = pxform("IAU_SUN", "J2000", epoch_lt) * sun_pole_sun

    # get rotation matrix for sun
    sun_rot_mat = pxform("IAU_SUN", "J2000", epoch_lt)

    # compute inclination
    # v1 = sun_pole_bary
    # v2 = OS_bary[1:3]
    # sun_angle = 90.0 - rad2deg(acos(dot(v1, v2) / (norm(v1) * norm(v2))))
    # @show sun_angle


    # set size of subgrid
    Nsubgrid = 10
# allocate memory
    dA_total_proj_mean = zeros(length(disk_ϕc), maximum(Nθ))
    mean_intensity = zeros(length(disk_ϕc), maximum(Nθ))
    mean_weight_v_no_cb = zeros(length(disk_ϕc), maximum(Nθ))
    mean_weight_v_cb = zeros(length(disk_ϕc), maximum(Nθ))
    mean_weight_v_earth_orb = zeros(length(disk_ϕc), maximum(Nθ))
    ra_mean = zeros(length(disk_ϕc), maximum(Nθ))
    de_mean = zeros(length(disk_ϕc), maximum(Nθ))
    mean_exti = zeros(length(disk_ϕc), maximum(Nθ))
    # vectors
    pole_vector_grid = fill(Vector{Float64}(undef, 3), Nsubgrid, Nsubgrid)
    SP_sun_pos = fill(Vector{Float64}(undef, 3), Nsubgrid, Nsubgrid)
    SP_sun_vel = fill(Vector{Float64}(undef, 3), Nsubgrid, Nsubgrid)
    SP_bary = fill(Vector{Float64}(undef, 6), Nsubgrid, Nsubgrid)
    SP_bary_pos = fill(Vector{Float64}(undef, 3), Nsubgrid, Nsubgrid)
    SP_bary_vel = fill(Vector{Float64}(undef, 3), Nsubgrid, Nsubgrid)
    OP_bary = fill(Vector{Float64}(undef, 6), Nsubgrid, Nsubgrid)
    # scalars
    mu_grid = zeros(Nsubgrid, Nsubgrid)
    convective_velocities = zeros(Nsubgrid, Nsubgrid)
    projected_velocities_no_cb = zeros(Nsubgrid, Nsubgrid)
    projected_velocities_cb = zeros(Nsubgrid, Nsubgrid)
    distance = zeros(Nsubgrid, Nsubgrid)
    v_scalar_grid = zeros(Nsubgrid, Nsubgrid)
    v_earth_orb_proj = zeros(Nsubgrid, Nsubgrid)


# loop over sub-tiles
    for i in eachindex(disk_ϕc)
        for j in 1:Nθ[i]
        #subdivide the tile
            ϕe_sub = range(disk_ϕe[i], disk_ϕe[i+1], length=Nsubgrid+1)
            θe_sub = range(disk_θe[i,j], disk_θe[i,j+1], length=Nsubgrid+1)
            ϕc_sub = get_grid_centers(ϕe_sub)
            θc_sub = get_grid_centers(θe_sub)
            subgrid = Iterators.product(ϕc_sub, θc_sub)

        #determine required position vectors
            #determine xyz stellar coordinates for lat/long grid
            SP_sun_pos .= map(x -> pgrrec("SUN", getindex(x,2), getindex(x,1), 0.0, sun_radius, 0.0), subgrid)

            # get differential rotation velocities
            v_scalar_grid .= map(x -> v_scalar(x...), subgrid)
            #convert v_scalar to from km/day km/s
            v_scalar_grid ./= 86400.0

            #determine pole vector for each patch
            pole_vector_grid!(SP_sun_pos, pole_vector_grid)

            #get velocity vector direction and set magnitude
            v_vector(SP_sun_pos, pole_vector_grid, v_scalar_grid, SP_sun_vel)

            for k in eachindex(SP_sun_pos)
                SP_bary_pos[k] .= (sun_rot_mat * SP_sun_pos[k])
                SP_bary_vel[k] .= (sun_rot_mat * SP_sun_vel[k])
                SP_bary[k] = vcat(SP_bary_pos[k], SP_bary_vel[k])
            end

            #get vector from obs to each patch on Sun's surface
            for k in eachindex(OP_bary)
                OP_bary[k] = OS_bary .+ SP_bary[k]
            end

        #calculate mu for each patch
            calc_mu_grid!(SP_bary, OP_bary, mu_grid)
            # move on if everything is off the grid
            all(mu_grid .< zero(T)) && continue

        #get projected velocity for each patch
            convective_velocities .= convective_blueshift_interpol.(mu_grid)
            projected!(SP_bary, OP_bary, projected_velocities_no_cb, projected_velocities_cb, convective_velocities)
            # convert from km/s to m/s
            projected_velocities_cb .*= 1000.0
            projected_velocities_no_cb .*= 1000.0

        #get relative orbital motion in m/s
            v_delta = OS_bary[4:6]
            for k in eachindex(v_earth_orb_proj)
                B = view(OP_bary[k], 1:3)
                angle = dot(B, v_delta) / (norm(B) * norm(v_delta))
                v_earth_orb_proj[k] = norm(v_delta) * angle
            end
            # convert from km/s to m/s
            v_earth_orb_proj .*= 1000.0

        #determine patches that are blocked by moon
            #calculate the distance between tile corner and moon
            for i in eachindex(OP_bary)
                distance[i] = calc_proj_dist2(OM_bary[1:3], OP_bary[i][1:3])
            end

            #get indices for visible patches
            idx1 = mu_grid .> 0.0
            idx3 = (idx1) .& (distance .> atan((moon_r)/norm(OM_bary))) 

        #calculate limb darkening weight for each patch
            if band == "NIR"
                LD_all = map(x -> quad_limb_darkening_NIR(x), mu_grid)
            end
            
            if band == "optical"
                LD_all = map(x -> quad_limb_darkening_optical(x), mu_grid)
            end

        #calculating the area element dA for each tile
            dϕ = step(ϕe_sub) 
            dθ = step(θe_sub) 
            dA_sub = map(x -> calc_dA(sun_radius, getindex(x,1), dϕ, dθ), subgrid)
            #get total projected, visible area of larger tile
            dA_total_proj = dA_sub .* mu_grid
            dA_total_proj_mean[i,j] = sum(view(dA_total_proj, idx1))

        #get ra and dec
            OP_ra_dec = SPICE.recrad.([x[1:3] for x in OP_bary])
            ra_mean[i,j] = mean(getindex.(OP_ra_dec, 2))
            de_mean[i,j] = mean(getindex.(OP_ra_dec, 3))

        #determine mean intensity
            mean_intensity[i,j] = mean(view(LD_all, idx3)) 
            #extinction 
            zenith_angle_matrix = rad2deg.(map(x -> calc_proj_dist2(x, EO_bary), OP_bary))
            extin = map(x -> 10^(-((1/cosd(x))*ext_coef[index])/2.5), zenith_angle_matrix)
            mean_exti[i,j] = mean(view(extin, idx3)) 

        #determine mean weighted velocity from sun given blocking from moon
            mean_weight_v_no_cb[i,j] = mean(view(projected_velocities_no_cb, idx3))
            mean_weight_v_cb[i,j] = mean(view(projected_velocities_cb, idx3)) 
            mean_weight_v_earth_orb[i,j] = mean(view(v_earth_orb_proj, idx3))
        end
    end

    #index for correct lat / lon disk grid
    idx_grid = mean_intensity .> 0.0

    contrast = (mean_intensity / NaNMath.maximum(mean_intensity)).^0.1
    brightness = mean_intensity .* mean_exti .* dA_total_proj_mean
    cheapflux = sum(view(brightness, idx_grid))

    #determine final mean intensity for disk grid
    final_mean_intensity = cheapflux 

    #determine final mean weighted velocity for disk grid
    final_weight_v_no_cb = sum(view(contrast .* mean_weight_v_no_cb .* brightness, idx_grid)) / cheapflux 
    final_weight_v_no_cb += mean(view(mean_weight_v_earth_orb, idx_grid)) 

    final_weight_v_cb = sum(view(contrast .* mean_weight_v_cb .* brightness, idx_grid)) / cheapflux
    final_weight_v_cb += mean(view(mean_weight_v_earth_orb, idx_grid)) 

    return final_weight_v_no_cb, final_weight_v_cb, final_mean_intensity, rad2deg.(ra_mean), rad2deg.(de_mean), mean_weight_v_no_cb, mean_weight_v_cb
end
