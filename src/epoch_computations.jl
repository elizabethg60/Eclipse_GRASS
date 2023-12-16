using Base.Iterators

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

    #query JPL horizons for E, S, M position (km) and velocities (km/s)
    earth_pv = spkssb(399,epoch,"J2000")[1:3] 
    sun_pv = spkssb(10,epoch,"J2000")[1:3] 
    moon_pv = spkssb(301,epoch,"J2000")[1:3] 

    earth_vel = spkssb(399,epoch,"J2000")[4:6]
    sun_vel = spkssb(10,epoch,"J2000")[4:6]
    moon_vel = spkssb(301,epoch,"J2000")[4:6]

    # set subgridding
    Nsubgrid = 10

    # allocate memory
    dA_total_proj_mean = zeros(length(disk_ϕc), maximum(Nθ))
    mean_intensity = zeros(length(disk_ϕc), maximum(Nθ))
    mean_weight_v_no_cb = zeros(length(disk_ϕc), maximum(Nθ))
    mean_weight_v_cb = zeros(length(disk_ϕc), maximum(Nθ))
    mean_weight_v_earth_rot= zeros(length(disk_ϕc), maximum(Nθ))
    mean_weight_v_earth_orb = zeros(length(disk_ϕc), maximum(Nθ))
    LD_all_mean = zeros(length(disk_ϕc), maximum(Nθ))
    dA_total_proj_out = zeros(length(disk_ϕc), maximum(Nθ))
    ra_mean = zeros(length(disk_ϕc), maximum(Nθ))
    de_mean = zeros(length(disk_ϕc), maximum(Nθ))
    mu_mean = zeros(length(disk_ϕc), maximum(Nθ))

    for i in eachindex(disk_ϕc)
        for j in 1:Nθ[i]
        # subdivide the tile
            ϕe_sub = range(disk_ϕe[i], disk_ϕe[i+1], length=Nsubgrid+1)
            θe_sub = range(disk_θe[i,j], disk_θe[i,j+1], length=Nsubgrid+1)
            ϕc_sub = get_grid_centers(ϕe_sub)
            θc_sub = get_grid_centers(θe_sub)
            subgrid = Iterators.product(ϕc_sub, θc_sub)
            
        #determine required position vectors
            #determine xyz stellar coordinates for lat/long grid
            SP_sun = map(x -> get_xyz.(sun_radius, x...), subgrid) #get_xyz_for_surface(sun_radius, num_lats = lats, num_lons = lons)
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
            OS_bary = BO_bary - sun_pv
            #get vector from observatory on earth's surface to moon center
            OM_bary = moon_pv .- BO_bary
            #get vector from barycenter to each patch on Sun's surface
            BP_bary = Matrix{Vector{Float64}}(undef,size(SP_bary)...)
            for i in eachindex(BP_bary)
                BP_bary[i] = sun_pv + SP_bary[i]
            end
            #vectors from observatory on Earth's surface to each patch on Sun's surface
            OP_bary = Matrix{Vector{Float64}}(undef,size(SP_bary)...)
            earth2patch_vectors(BP_bary, BO_bary, OP_bary)	


            # get ra and dec
            OP_ra_dec = SPICE.recrad.(OP_bary)


        #calculate mu for each patch
            mu_grid = Matrix{Float64}(undef,size(SP_bary)...)
            calc_mu_grid!(SP_bary, OS_bary, mu_grid) 

            all(mu_grid .< zero(T)) && continue

            # plt.pcolormesh(rad2deg.(getindex.(OP_ra_dec, 2)), rad2deg.(getindex.(OP_ra_dec, 3)), mu_grid, cmap="viridis", vmin=0.0, vmax=1.0)

        #determine velocity vectors
            #determine velocity scalar for each patch 
            v_scalar_grid = map(x -> v_scalar(x...), subgrid)

            #convert v_scalar to from km/day m/s
            v_scalar_grid ./= 86.4

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
            projected!(velocity_vector_ICRF, OP_bary, projected_velocities_no_cb, projected_velocities_cb, convective_velocities)

            # plt.pcolormesh(rad2deg.(getindex.(OP_ra_dec, 2)), rad2deg.(getindex.(OP_ra_dec, 3)), projected_velocities_no_cb, cmap="seismic", vmin=-2000.0, vmax=2000.0)

            earth_v_scalar = (2*π*earth_radius*cosd(obs_lat))/86.4
            earth_pole_vector = EO_earth - [0.0, 0.0, EO_earth[3]]
            temp = cross(earth_pole_vector, [0.0,0.0,earth_radius])
            temp ./= norm(temp)
            temp .*= earth_v_scalar
            velocity_vector_earth = [EO_earth;temp]
            velocity_vector_earth_ICRF = sxform("IAU_EARTH", "J2000", epoch) * velocity_vector_earth
            velocity_vector_earth_ICRF = velocity_vector_earth_ICRF[4:6]

            # project earth_rot velocity vector onto patches
            v_earth_rot_proj = zeros(Nsubgrid, Nsubgrid)
            for k in eachindex(v_earth_rot_proj)
                B = OP_bary[k]
                angle = cos(π - acos(dot(B, velocity_vector_earth_ICRF) / (norm(B) * norm(velocity_vector_earth_ICRF))))
                v_earth_rot_proj[k] = norm(velocity_vector_earth_ICRF) * angle
            end

            # plt.pcolormesh(rad2deg.(getindex.(OP_ra_dec, 2)), rad2deg.(getindex.(OP_ra_dec, 3)), v_earth_rot_proj, cmap="viridis", vmin=135, vmax=145)

            # get relative orbital motion in m/s
            v_delta = (sun_vel .- earth_vel) .* 1000
            v_earth_orb_proj = zeros(Nsubgrid, Nsubgrid)
            for k in eachindex(v_earth_orb_proj)
                B = OP_bary[k]
                angle = cos(π - acos(dot(B, v_delta) / (norm(B) * norm(v_delta))))
                v_earth_orb_proj[k] = norm(v_delta) * angle
            end

            # plt.pcolormesh(rad2deg.(getindex.(OP_ra_dec, 2)), rad2deg.(getindex.(OP_ra_dec, 3)), v_earth_orb_proj, cmap="viridis", vmin=35, vmax=65)
            
        #determine patches that are blocked by moon 
            #calculate the distance between tile corner and moon
            distance = map(x -> calc_proj_dist2(x, OM_bary), OP_bary)

            distance_clone = deepcopy(distance)
            distance_clone[distance .<=atan((moon_r)/norm(OM_bary))] .= 0.0

            # plt.pcolormesh(rad2deg.(getindex.(OP_ra_dec, 2)), rad2deg.(getindex.(OP_ra_dec, 3)), distance_clone, cmap="viridis", vmin=0.0, vmax=deg2rad(0.5))

            #calculate limb darkening weight for each patch 
            if band == "NIR"
                LD_all = map(x -> quad_limb_darkening_NIR(x), mu_grid)
            end

            zenith_angle_matrix = rad2deg.(map(x -> calc_proj_dist2(x, EO_bary), OP_bary))
            if band == "optical"
                LD_all = map((x,y) -> quad_limb_darkening_optical(x, index, y), mu_grid, zenith_angle_matrix)
            end

            # plt.pcolormesh(rad2deg.(getindex.(OP_ra_dec, 2)), rad2deg.(getindex.(OP_ra_dec, 3)), LD_all, cmap="afmhot", vmin=0.0, vmax=1.0)

            #get indices for visible patches
            idx1 = mu_grid .> 0.0
            idx3 = (idx1) .& (distance .> atan((moon_r)/norm(OM_bary)))


            ra_mean[i,j] = mean(getindex.(OP_ra_dec, 2))
            de_mean[i,j] = mean(getindex.(OP_ra_dec, 3))
            mu_mean[i,j] = mean(view(mu_grid, idx1))


            # plt.pcolormesh(rad2deg.(getindex.(OP_ra_dec, 2)), rad2deg.(getindex.(OP_ra_dec, 3)), idx3 .* LD_all, cmap="afmhot", vmin=0.0, vmax=1.0)
            
            #calculating the area element dA for each tile
            dϕ = step(ϕe_sub) # mean(diff(ϕc_sub))
            dθ = step(θe_sub) # mean(diff(θc_sub))
            dA_sub = map(x -> calc_dA(sun_radius, getindex(x,1), dϕ, dθ), subgrid)
            #get total projected, visible area of larger tile
            dp_sub = map((x,y) -> abs(dot(x,y)), OP_bary, SP_bary) ./ (norm.(OP_bary) .* norm.(SP_bary))

            dA_total_proj = dA_sub .* dp_sub 

            # plt.pcolormesh(rad2deg.(getindex.(OP_ra_dec, 2)), rad2deg.(getindex.(OP_ra_dec, 3)), dA_total_proj, cmap="viridis")

            LD_all_mean[i,j] = mean(view(LD_all, idx3))
            dA_total_proj_mean[i,j] = sum(view(dA_total_proj, idx1))

        #determine mean intensity
            # mean_intensity[i,j] = sum(view(LD_all .* dA_total_proj, idx3)) / sum(view(dA_total_proj, idx3))
            mean_intensity[i,j] = sum(view(LD_all, idx3)) / sum(idx3)

        #determine mean weighted velocity from sun given blocking from moon 

            mean_weight_v_no_cb[i,j] = sum(view(LD_all .* dA_total_proj .* projected_velocities_no_cb, idx3)) / sum(view(LD_all .* dA_total_proj, idx3))
            mean_weight_v_cb[i,j] =  sum(view(LD_all .* dA_total_proj .* projected_velocities_cb, idx3)) / sum(view(LD_all .* dA_total_proj, idx1))

            mean_weight_v_earth_rot[i,j] = sum(view(v_earth_rot_proj .* LD_all .* dA_total_proj, idx3)) / sum(view(LD_all .* dA_total_proj, idx3))
            mean_weight_v_earth_orb[i,j] = sum(view(v_earth_orb_proj .* LD_all .* dA_total_proj, idx3)) / sum(view(LD_all .* dA_total_proj, idx3))

        end
    end
    idx_grid = mean_intensity .> 0.0

    # @show sum(dA_total_proj_mean) / sun_radius^2.0

    # ax1 = plt.gca()

    # xlims = ax1.get_xlim()
    # ylims = ax1.get_ylim()

    # ax1.set_xlim(xlims[2]+0.01, xlims[1]-0.01)
    # ax1.set_ylim(ylims[1]-0.01, ylims[2]+0.01)
    # plt.colorbar()
    # ax1.set_aspect("equal")
    # plt.show()

    #double weighting occuring? 
    final_mean_intensity = sum(view(mean_intensity .* dA_total_proj_mean, idx_grid)) / sum(view(dA_total_proj_mean, idx_grid)) 


    #### INCLUDE LIMB DARKENING
    final_weight_v_no_cb = sum(view(mean_intensity .* dA_total_proj_mean .* mean_weight_v_no_cb, idx_grid)) / sum(view(mean_intensity .* dA_total_proj_mean, idx_grid))
    final_weight_v_no_cb -= sum(view(mean_weight_v_earth_rot .* mean_intensity .* dA_total_proj_mean, idx_grid)) / sum(view(mean_intensity .* dA_total_proj_mean, idx_grid))
    final_weight_v_no_cb -= sum(view(mean_weight_v_earth_orb .* mean_intensity .* dA_total_proj_mean, idx_grid)) / sum(view(mean_intensity .* dA_total_proj_mean, idx_grid))

    final_weight_v_cb = sum(view(mean_intensity .* dA_total_proj_mean .* mean_weight_v_cb, idx_grid)) / sum(view(mean_intensity .* dA_total_proj_mean, idx_grid))
    final_weight_v_cb -= sum(view(mean_weight_v_earth_rot .* mean_intensity .* dA_total_proj_mean, idx_grid)) / sum(view(mean_intensity .* dA_total_proj_mean, idx_grid))
    final_weight_v_cb -= sum(view(mean_weight_v_earth_orb .* mean_intensity .* dA_total_proj_mean, idx_grid)) / sum(view(mean_intensity .* dA_total_proj_mean, idx_grid))

    #return mean_weight_v_no_cb, mean_weight_v_cb, mean_intensity, rad2deg.(getindex.(OP_ra_dec,2)), rad2deg.(getindex.(OP_ra_dec,3)), projected_velocities_no_cb, projected_velocities_cb
    return final_weight_v_no_cb, final_weight_v_cb, final_mean_intensity

end
