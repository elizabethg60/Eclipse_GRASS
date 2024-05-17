function neid_april_loop(lats::T) where T
    #convert from utc to et as needed by SPICE
    time_stamps = utc2et.(neid_april_timestamps)

    #NEID location
    obs_lat = 31.9583 
    obs_long = -111.5967  
    alt = 2.097938 

    wavelength = [5250.2084, 5250.6453, 5379.5734, 5381.0216, 5382.2562, 5383.368, 5432.546, 5432.947, 5434.5232, 5435.8577, 5436.2945, 5436.5875, 5576.0881, 5578.718, 6149.246, 6151.617, 6169.042, 6169.563, 6170.5056, 6173.3344, 6301.5008, 6302.4932]

    RV_list_no_cb_final = Vector{Vector{Float64}}(undef,size(wavelength)...)
    RV_list_cb_final = Vector{Vector{Float64}}(undef,size(wavelength)...)
    intensity_list_final = Vector{Vector{Float64}}(undef,size(wavelength)...)
    for lambda in 1:length(wavelength)
        println(wavelength[lambda]/10.0)
        RV_list_no_cb = Vector{Float64}(undef,size(time_stamps)...)
        RV_list_cb = Vector{Float64}(undef,size(time_stamps)...)
        intensity_list = Vector{Float64}(undef,size(time_stamps)...)
        #run compute_rv for each timestamp
        for i in 1:length(time_stamps)
            RV_no_cb, RV_cb, RV_cb_new, intensity = compute_rv(lats, time_stamps[i], obs_long, obs_lat, alt, "optical", wavelength[lambda]/10.0, i)
            RV_list_no_cb[i] = RV_no_cb
            RV_list_cb[i] = RV_cb
            intensity_list[i] = intensity
        end
        RV_list_no_cb_final[lambda] = RV_list_no_cb
        RV_list_cb_final[lambda] = RV_list_cb 
        intensity_list_final[lambda] = intensity_list
    end

    @save "src/plots/NEID_April/model_data.jld2"
    jldopen("src/plots/NEID_April/model_data.jld2", "a+") do file
        file["RV_list_no_cb"] = RV_list_no_cb_final 
        file["RV_list_cb"] = RV_list_cb_final 
        file["intensity_list"] = intensity_list_final
    end
end

function neid_october_loop(lats::T) where T
    #convert from utc to et as needed by SPICE
    time_stamps = utc2et.(neid_timestamps)

    #NEID location
    obs_lat = 31.9583 
    obs_long = -111.5967  
    alt = 2.097938 

    wavelength = [5250.2084, 5250.6453, 5379.5734, 5381.0216, 5382.2562, 5383.368, 5432.546, 5432.947, 5434.5232, 5435.8577, 5436.2945, 5436.5875, 5576.0881, 5578.718, 6149.246, 6151.617, 6169.042, 6169.563, 6170.5056, 6173.3344, 6301.5008, 6302.4932]

    RV_list_no_cb_final = Vector{Vector{Float64}}(undef,size(wavelength)...)
    RV_list_cb_final = Vector{Vector{Float64}}(undef,size(wavelength)...)
    intensity_list_final = Vector{Vector{Float64}}(undef,size(wavelength)...)
    for lambda in 1:length(wavelength)
        println(wavelength[lambda]/10.0)
        RV_list_no_cb = Vector{Float64}(undef,size(time_stamps)...)
        RV_list_cb = Vector{Float64}(undef,size(time_stamps)...)
        intensity_list = Vector{Float64}(undef,size(time_stamps)...)
        # airmass_list = Vector{Float64}(undef,size(time_stamps)...)
        #run compute_rv for each timestamp
        for i in 1:length(time_stamps)
            # airmass_list[i] = compute_rv(lats, time_stamps[i], obs_long, obs_lat, alt, "optical", wavelength[lambda]/10.0, i)
            RV_no_cb, RV_cb, RV_cb_new, intensity = compute_rv(lats, time_stamps[i], obs_long, obs_lat, alt, "optical", wavelength[lambda]/10.0, i)
            RV_list_no_cb[i] = RV_no_cb
            RV_list_cb[i] = RV_cb
            intensity_list[i] = intensity
        end
        RV_list_no_cb_final[lambda] = RV_list_no_cb
        RV_list_cb_final[lambda] = RV_list_cb 
        intensity_list_final[lambda] = intensity_list
    end

    @save "src/plots/NEID_October/model_data.jld2"
    jldopen("src/plots/NEID_October/model_data.jld2", "a+") do file
        file["RV_list_no_cb"] = RV_list_no_cb_final 
        file["RV_list_cb"] = RV_list_cb_final 
        file["intensity_list"] = intensity_list_final
    end
end

function gottingen_loop(lats::T) where T
    """
    computes RV for each timestamp for the gottingen eclipse 
    """
    #convert from utc to et as needed by SPICE
    #time_stamps = utc2et.(reiners_timestamps)
    time_stamps = utc2et.(reiners_50)

    #Gottingen location
    obs_lat = 51.560583 
    obs_long = 9.944333
    alt = 0.201

    wavelength = 543.4

    RV_list_no_cb = Vector{Float64}(undef,size(time_stamps)...)
    RV_list_cb = Vector{Float64}(undef,size(time_stamps)...)
    RV_list_cb_new = Vector{Float64}(undef,size(time_stamps)...)
    intensity_list = Vector{Float64}(undef,size(time_stamps)...)
    #run compute_rv for each timestamp
    for i in 1:length(time_stamps)
        RV_no_cb, RV_cb, RV_cb_new, intensity = compute_rv(lats, time_stamps[i], obs_long, obs_lat, alt, "optical", wavelength, i)#, ext = true)
        RV_list_no_cb[i] = RV_no_cb
        RV_list_cb[i] = RV_cb
        RV_list_cb_new[i] = RV_cb_new
        intensity_list[i] = intensity

        # rv_bin = Vector{Float64}(undef,12)
        # cb_bin = Vector{Float64}(undef,12)
        # new_cb_bin = Vector{Float64}(undef,12)
        # intensity_bin = Vector{Float64}(undef,12)
        # sample = utc2et.(reiners_finer_sample_50[i]) 
        # for j in 1:12
        #     rv_bin[j], cb_bin[j], new_cb_bin[j], intensity_bin[j]  = (compute_rv(lats, sample[j], obs_long, obs_lat, alt, "optical", wavelength, i, ext = true))
        # end
        # RV_list_no_cb[i] = mean(rv_bin)
        # RV_list_cb[i] = mean(cb_bin)
        # RV_list_cb_new[i] = mean(new_cb_bin)
        # intensity_list[i] = mean(intensity_bin)
    end
    
    @save "src/plots/Reiners/model_data.jld2"
    jldopen("src/plots/Reiners/model_data.jld2", "a+") do file
        file["RV_list_no_cb"] = RV_list_no_cb 
        file["RV_list_cb"] = RV_list_cb 
        file["RV_list_cb_new"] = RV_list_cb_new
        file["intensity_list"] = intensity_list
    end
    return nothing
end

function expres_loop(lats::T) where T
    """
    computes RV for each timestamp for the expres eclipse 
    """
    #convert from utc to et as needed by SPICE
    time_stamps = utc2et.(expres_timestamps)

    #EXPRES location
    obs_lat = 34.744444
    obs_long = -111.421944 
    alt = 2.359152

    wavelength = 543.4

    RV_list_no_cb = Vector{Float64}(undef,size(time_stamps)...)
    RV_list_cb = Vector{Float64}(undef,size(time_stamps)...)
    RV_list_cb_new = Vector{Float64}(undef,size(time_stamps)...)
    intensity_list = Vector{Float64}(undef,size(time_stamps)...)
    #run compute_rv for each timestamp
    for i in 1:length(time_stamps)
        RV_no_cb, RV_cb, RV_cb_new, intensity = compute_rv(lats, time_stamps[i], obs_long, obs_lat, alt, "optical", wavelength, i)
        RV_list_no_cb[i] = RV_no_cb
        RV_list_cb[i] = RV_cb
        RV_list_cb_new[i] = RV_cb_new
        intensity_list[i] = intensity
    end

    @save "src/plots/EXPRES/model_data.jld2"
    jldopen("src/plots/EXPRES/model_data.jld2", "a+") do file
        file["RV_list_no_cb"] = RV_list_no_cb 
        file["RV_list_cb"] = RV_list_cb 
        file["RV_list_cb_new"] = RV_list_cb_new
        file["intensity_list"] = intensity_list
    end
end

function boulder_loop(lats::T) where T
    """
    computes RV for each timestamp for the boulder eclipse 
    """
    #convert from utc to et as needed by SPICE
    time_stamps = utc2et.(boulder_timestamps)

    #Boulder location
    obs_lat = 39.995380
    obs_long = -105.262390
    alt = 1.6523

    wavelength = 543.4

    RV_list_no_cb = Vector{Float64}(undef,size(time_stamps)...)
    RV_list_cb = Vector{Float64}(undef,size(time_stamps)...)
    RV_list_cb_new = Vector{Float64}(undef,size(time_stamps)...)
    intensity_list = Vector{Float64}(undef,size(time_stamps)...)
    #run compute_rv for each timestamp
    for i in 1:length(time_stamps)
        RV_no_cb, RV_cb, RV_cb_new, intensity  = compute_rv(lats, time_stamps[i], obs_long, obs_lat, alt, "NIR", wavelength, i)
        RV_list_no_cb[i] = RV_no_cb
        RV_list_cb[i] = RV_cb
        RV_list_cb_new[i] = RV_cb_new
        intensity_list[i] = intensity
    end

    @save "src/plots/Boulder/model_data.jld2"
    jldopen("src/plots/Boulder/model_data.jld2", "a+") do file
        file["RV_list_no_cb"] = RV_list_no_cb 
        file["RV_list_cb"] = RV_list_cb 
        file["RV_list_cb_new"] = RV_list_cb_new
        file["intensity_list"] = intensity_list
    end
end

# function chi2(lats::T) where T
#     #Gottingen location
#     obs_lat = 51.560583 
#     obs_long = 9.944333
#     alt = 0.201

#     wavelength = 580

#     seconds = range(63,63,step=1)
#     for sec in seconds
#         print(sec)
#         time_stamps_string = Vector{String}(undef,size(reiners_timestamps)...) 
#         for i in 1:length(time_stamps_string)
#             time_stamps_string[i] = Dates.format(Dates.DateTime(reiners_timestamps[i], "yyyy-mm-ddTHH:MM:SS.ss") + Dates.Second(sec), "yyyy-mm-ddTHH:MM:SS.ss")
#         end
#         time_stamps = utc2et.(time_stamps_string)

#         RV_list_no_cb = Vector{Float64}(undef,size(time_stamps)...)
#         RV_list_cb = Vector{Float64}(undef,size(time_stamps)...)
#         RV_list_cb_new = Vector{Float64}(undef,size(time_stamps)...)
#         intensity_arr = Vector{Float64}(undef,size(time_stamps)...)
#         #run compute_rv for each timestamp
#         for i in 1:length(time_stamps)
#             RV_no_cb, RV_cb, RV_cb_new, intensity, = compute_rv(lats, time_stamps[i], obs_long, obs_lat, alt, "optical", wavelength, i, ext = true)
#             RV_list_no_cb[i] = RV_no_cb
#             RV_list_cb[i] = RV_cb
#             RV_list_cb_new[i] = RV_cb_new
#             intensity_arr[i] = intensity
#         end

#         @save "src/plots/Reiners/model_data.jld2"
#         jldopen("src/plots/Reiners/model_data.jld2", "a+") do file
#             file["RV_list_no_cb"] = RV_list_no_cb 
#             file["RV_list_cb"] = RV_list_cb 
#             file["RV_list_cb_new"] = RV_list_cb_new 
#             file["timestamps"] = time_stamps_string
#             file["intensity_list"] = intensity_arr
#         end
#     end
#     return nothing
# end