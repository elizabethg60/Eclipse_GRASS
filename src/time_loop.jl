function neid_loop(lats::T, lons::T) where T
    """
    computes RV for each timestamp for the NEID eclipse 

    lats: number of latitude grid cells
    lons: number of longitude grid cells
    """
    #convert from utc to et as needed by SPICE
    time_stamps = utc2et.(neid_timestamps)

    #NEID location
    obs_lat = 31.9583 
    obs_long = -111.5967  
    alt = 2.097938 

    RV_list_no_cb = Vector{Float64}(undef,size(time_stamps)...)
    RV_list_cb = Vector{Float64}(undef,size(time_stamps)...)
    intensity_list = Vector{Float64}(undef,size(time_stamps)...)
    RA_list = Vector{Matrix{Float64}}(undef,size(time_stamps)...)
    dec_list = Vector{Matrix{Float64}}(undef,size(time_stamps)...)
    vel_no_cb = Vector{Matrix{Float64}}(undef,size(time_stamps)...)
    vel_cb = Vector{Matrix{Float64}}(undef,size(time_stamps)...)
    #run compute_rv (serial) for each timestamp
    for i in 1:length(time_stamps)
        RV_no_cb, RV_cb, intensity, ra, dec, projected_v_no_cb, projected_v_cb  = compute_rv(lats, lons, time_stamps[i], obs_long, obs_lat, alt, "optical", i)
        RV_list_no_cb[i] = RV_no_cb
        RV_list_cb[i] = RV_cb
        intensity_list[i] = intensity
        RA_list[i] = ra
        dec_list[i] = dec
        vel_no_cb[i] = projected_v_no_cb
        vel_cb[i] = projected_v_cb
    end

    @save "src/plots/NEID/model_data.jld2"
    jldopen("src/plots/NEID/model_data.jld2", "a+") do file
        file["RV_list_no_cb"] = RV_list_no_cb 
        file["RV_list_cb"] = RV_list_cb 
        file["intensity_list"] = intensity_list
        file["RA_list"] = RA_list
        file["dec_list"] = dec_list
        file["vel_no_cb"] = vel_no_cb
        file["vel_cb"] = vel_cb
    end
end

function gottingen_loop(lats::T, lons::T) where T
    """
    computes RV for each timestamp for the gottingen eclipse 
    """
    #convert from utc to et as needed by SPICE
    time_stamps = utc2et.(reiners_timestamps)

    #Gottingen location
    obs_lat = 51.560583 
    obs_long = 9.944333
    alt = 0.201

    RV_list_no_cb = Vector{Float64}(undef,size(time_stamps)...)
    RV_list_cb = Vector{Float64}(undef,size(time_stamps)...)
    intensity_list = Vector{Float64}(undef,size(time_stamps)...)
    RA_list = Vector{Matrix{Float64}}(undef,size(time_stamps)...)
    dec_list = Vector{Matrix{Float64}}(undef,size(time_stamps)...)
    vel_no_cb = Vector{Matrix{Float64}}(undef,size(time_stamps)...)
    vel_cb = Vector{Matrix{Float64}}(undef,size(time_stamps)...)
    #run compute_rv (serial) for each timestamp
    for i in 1:length(time_stamps)
        print("\r i = " * string(i) * " of " * string(length(time_stamps)))
        # if i != round(Int, length(time_stamps)/2)
        #     continue
        # end
        RV_no_cb, RV_cb, intensity = compute_rv(lats, lons, time_stamps[i], obs_long, obs_lat, alt, "optical", i)
        RV_list_no_cb[i] = RV_no_cb
        RV_list_cb[i] = RV_cb
        intensity_list[i] = intensity
        # RA_list[i] = ra
        # dec_list[i] = dec
        # vel_no_cb[i] = projected_v_no_cb
        # vel_cb[i] = projected_v_cb
        #, ra, dec, projected_v_no_cb, projected_v_cb 

        # rv_bin = Vector{Float64}(undef,12)
        # sample = utc2et.(reiners_finer_sample[i]) 
        # for j in 1:12
        #     rv_bin[j] = (compute_rv(lats, lons, sample[j], obs_long, obs_lat, alt, "optical", i)[1])
        # end
        # RV_list_no_cb[i] = mean(rv_bin)
    end
    println()

    @save "src/plots/Reiners/model_data.jld2"
    jldopen("src/plots/Reiners/model_data.jld2", "a+") do file
        file["RV_list_no_cb"] = RV_list_no_cb 
        file["RV_list_cb"] = RV_list_cb 
        file["intensity_list"] = intensity_list
        # file["RA_list"] = RA_list
        # file["dec_list"] = dec_list
        # file["vel_no_cb"] = vel_no_cb
        # file["vel_cb"] = vel_cb
    end
    return nothing
end

function expres_loop(lats::T, lons::T) where T
    """
    computes RV for each timestamp for the expres eclipse 
    """
    #convert from utc to et as needed by SPICE
    time_stamps = utc2et.(expres_timestamps)

    #EXPRES location
    obs_lat = 34.744444
    obs_long = -111.421944 
    alt = 2.359152

    RV_list_no_cb = Vector{Float64}(undef,size(time_stamps)...)
    RV_list_cb = Vector{Float64}(undef,size(time_stamps)...)
    intensity_list = Vector{Float64}(undef,size(time_stamps)...)
    RA_list = Vector{Matrix{Float64}}(undef,size(time_stamps)...)
    dec_list = Vector{Matrix{Float64}}(undef,size(time_stamps)...)
    vel_no_cb = Vector{Matrix{Float64}}(undef,size(time_stamps)...)
    vel_cb = Vector{Matrix{Float64}}(undef,size(time_stamps)...)
    #run compute_rv (serial) for each timestamp
    for i in 1:length(time_stamps)
        RV_no_cb, RV_cb, intensity, ra, dec, projected_v_no_cb, projected_v_cb  = compute_rv(lats, lons, time_stamps[i], obs_long, obs_lat, alt, "optical", i)
        RV_list_no_cb[i] = RV_no_cb
        RV_list_cb[i] = RV_cb
        intensity_list[i] = intensity
        RA_list[i] = ra
        dec_list[i] = dec
        vel_no_cb[i] = projected_v_no_cb
        vel_cb[i] = projected_v_cb
    end

    @save "src/plots/EXPRES/model_data.jld2"
    jldopen("src/plots/EXPRES/model_data.jld2", "a+") do file
        file["RV_list_no_cb"] = RV_list_no_cb 
        file["RV_list_cb"] = RV_list_cb 
        file["intensity_list"] = intensity_list
        file["RA_list"] = RA_list
        file["dec_list"] = dec_list
        file["vel_no_cb"] = vel_no_cb
        file["vel_cb"] = vel_cb
    end
end

function boulder_loop(lats::T, lons::T) where T
    """
    computes RV for each timestamp for the boulder eclipse 
    """
    #convert from utc to et as needed by SPICE
    time_stamps = utc2et.(boulder_timestamps)

    #Boulder location
    obs_lat = 39.995380
    obs_long = -105.262390
    alt = 1.6523

    RV_list_no_cb = Vector{Float64}(undef,size(time_stamps)...)
    RV_list_cb = Vector{Float64}(undef,size(time_stamps)...)
    intensity_list = Vector{Float64}(undef,size(time_stamps)...)
    RA_list = Vector{Matrix{Float64}}(undef,size(time_stamps)...)
    dec_list = Vector{Matrix{Float64}}(undef,size(time_stamps)...)
    vel_no_cb = Vector{Matrix{Float64}}(undef,size(time_stamps)...)
    vel_cb = Vector{Matrix{Float64}}(undef,size(time_stamps)...)
    #run compute_rv (serial) for each timestamp
    for i in 1:length(time_stamps)
        RV_no_cb, RV_cb, intensity, ra, dec, projected_v_no_cb, projected_v_cb  = compute_rv(lats, lons, time_stamps[i], obs_long, obs_lat, alt, "NIR", i)
        RV_list_no_cb[i] = RV_no_cb
        RV_list_cb[i] = RV_cb
        intensity_list[i] = intensity
        RA_list[i] = ra
        dec_list[i] = dec
        vel_no_cb[i] = projected_v_no_cb
        vel_cb[i] = projected_v_cb
    end

    @save "src/plots/Boulder/model_data.jld2"
    jldopen("src/plots/Boulder/model_data.jld2", "a+") do file
        file["RV_list_no_cb"] = RV_list_no_cb 
        file["RV_list_cb"] = RV_list_cb 
        file["intensity_list"] = intensity_list
        file["RA_list"] = RA_list
        file["dec_list"] = dec_list
        file["vel_no_cb"] = vel_no_cb
        file["vel_cb"] = vel_cb
    end
end
