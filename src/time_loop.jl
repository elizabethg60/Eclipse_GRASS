function gottingen_loop(lats::T) where T
    """
    computes RV for each timestamp for the gottingen eclipse 
    """
    #convert from utc to et as needed by SPICE
    time_stamps = utc2et.(reiners_timestamps)

    #Gottingen location
    obs_lat = 51.560583 
    obs_long = 9.944333
    alt = 0.201

    wavelength = 5798.8

    RV_list_no_cb = Vector{Float64}(undef,size(time_stamps)...)
    RV_list_cb = Vector{Float64}(undef,size(time_stamps)...)
    RV_list_cb_new = Vector{Float64}(undef,size(time_stamps)...)
    intensity_list = Vector{Float64}(undef,size(time_stamps)...)
    #run compute_rv for each timestamp
    ext_coeff = 0.4 .- 0.3 .* (time_stamps .- time_stamps[1]) ./ (time_stamps[length(time_stamps)] - time_stamps[1])

    for i in 1:length(time_stamps)
        RV_no_cb, RV_cb, RV_cb_new, intensity = compute_rv(lats, time_stamps[i], obs_long, obs_lat, alt, wavelength, i, ext_coeff[i], ext = true)
        RV_list_no_cb[i] = RV_no_cb
        RV_list_cb[i] = RV_cb
        RV_list_cb_new[i] = RV_cb_new
        intensity_list[i] = intensity

    end
    
    @save "Confirm_Reiners/data/model_data_new.jld2"
    jldopen("Confirm_Reiners/data/model_data_new.jld2", "a+") do file
        file["RV_list_no_cb"] = RV_list_no_cb 
        file["RV_list_cb"] = RV_list_cb 
        file["RV_list_cb_new"] = RV_list_cb_new
        file["intensity_list"] = intensity_list
    end
    return nothing
end

function figures()
    boulder_timestamps = ["2023-10-14T13:25:37"]
    time_stamps = utc2et.(boulder_timestamps)

    #Boulder location
    obs_lat = 39.995380
    obs_long = -105.262390
    alt = 1.6523
    wavelength = 15652.79

    #NEID location
    # obs_lat = 31.9583 
    # obs_long = -111.5967  
    # alt = 2.097938 
    # wavelength = 5434.5232

    for i in 1:1
        compute_rv(50, time_stamps[i], obs_long, obs_lat, alt, wavelength, i, 1.0, ext = false)
    end
end