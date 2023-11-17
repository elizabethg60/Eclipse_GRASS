#parallel v1 multi-threading 
function parallel_loop()
        initial_N = 200
        final_N = 375
        N_steps = range(initial_N, final_N, step = 2)
        N_steps = Int.(N_steps)
        time_vector = Vector{Float64}(undef,size(N_steps)...) 
        for i in 1:length(N_steps)
            time_vector[i] = compute_rv_pa(N_steps[i], N_steps[i]*2, utc2et("2015-03-20T09:42:00"), 9.905548, 51.54548, 0.15, "optical")[3]                       
        end
        
        @save "src/test/grid_parallel.jld2"
        jldopen("src/test/grid_parallel.jld2", "a+") do file
            file["N_steps"] = N_steps 
            file["time_vector"] = time_vector
        end
end
    
function serial_loop()
        initial_N = 200
        final_N = 375
        N_steps = range(initial_N, final_N, step = 2)
        N_steps = Int.(N_steps)
        time_vector = Vector{Float64}(undef,size(N_steps)...)
        for i in 1:length(N_steps)
            time_vector[i] = compute_rv(N_steps[i], N_steps[i]*2, utc2et("2015-03-20T09:42:00"), 0, 9.905548, 51.54548, 0.15, "optical")[3]
        end
        
        @save "src/test/grid_serial.jld2"
        jldopen("src/test/grid_serial.jld2", "a+") do file
            file["N_steps"] = N_steps 
            file["time_vector"] = time_vector
        end
end
     
#-------------------------------------------------------------------------------------
    
function gottingen_loop(lats::T, lons::T) where T
    time_stamps = utc2et.(reiners_timestamps)

    obs_lat = 51.54548 
    obs_long = 9.905548
    alt = 0.15

    RV_list = Vector{Float64}(undef,size(time_stamps)...)
    intensity_list = Vector{Float64}(undef,size(time_stamps)...)
    for i in 1:length(time_stamps)
        rv, intensity = compute_rv(lats, lons, time_stamps[i], i, obs_long, obs_lat, alt, "optical", obs = "Reiners")
        RV_list[i] = rv
        intensity_list[i] = intensity
    end

    @save "src/plots/Reiners/rv_intensity.jld2"
    jldopen("src/plots/Reiners/rv_intensity.jld2", "a+") do file
        file["RV_list"] = RV_list 
        file["intensity_list"] = intensity_list
        file["timestamps"] = et2utc.(time_stamps, "ISOC", 0)
    end
end

function kitt_loop(lats::T, lons::T) where T
    time_stamps = utc2et.(neid_timestamps)

    obs_lat = 31.9583 
    obs_long = 360-111.5967  
    alt = 2.097938

    RV_list = Vector{Float64}(undef,size(time_stamps)...)
    intensity_list = Vector{Float64}(undef,size(time_stamps)...)
    for i in 1:length(time_stamps)
        rv, intensity = compute_rv(lats, lons, time_stamps[i], i, obs_long, obs_lat, alt, "optical", obs = "NEID")
        RV_list[i] = rv
        intensity_list[i] = intensity
    end

    @save "src/plots/NEID/rv_intensity.jld2"
    jldopen("src/plots/NEID/rv_intensity.jld2", "a+") do file
        file["RV_list"] = RV_list 
        file["intensity_list"] = intensity_list
        file["timestamps"] = et2utc.(time_stamps, "ISOC", 0)
    end
end

function low_loop(lats::T, lons::T) where T
    #array of timestamps 
    time_stamps = utc2et.(expres_timestamps)

    obs_lat = 34.744444
    obs_long = 360-111.421944 
    alt = 2.359152

    RV_list = Vector{Float64}(undef,size(time_stamps)...)
    intensity_list = Vector{Float64}(undef,size(time_stamps)...)
    for i in 1:length(time_stamps)
        rv, intensity = compute_rv(lats, lons, time_stamps[i], i, obs_long, obs_lat, alt, "optical", obs = "EXPRES")
        RV_list[i] = rv
        intensity_list[i] = intensity
    end

    @save "src/plots/EXPRES/rv_intensity.jld2"
    jldopen("src/plots/EXPRES/rv_intensity.jld2", "a+") do file
        file["RV_list"] = RV_list 
        file["intensity_list"] = intensity_list
        file["timestamps"] = et2utc.(time_stamps, "ISOC", 0)
    end
end

function boulder_loop(lats::T, lons::T) where T
    #array of timestamps 
    time_stamps = utc2et.(boulder_timestamps)

    obs_lat = 39.995380
    obs_long = 360-105.262390
    alt = 1.6523

    RV_list = Vector{Float64}(undef,size(time_stamps)...)
    intensity_list = Vector{Float64}(undef,size(time_stamps)...)
    for i in 1:length(time_stamps)
        rv, intensity = compute_rv(lats, lons, time_stamps[i], i, obs_long, obs_lat, alt, "NIR", obs = "Boulder")
        RV_list[i] = rv
        intensity_list[i] = intensity
    end

    @save "src/plots/Boulder/rv_intensity.jld2"
    jldopen("src/plots/Boulder/rv_intensity.jld2", "a+") do file
        file["RV_list"] = RV_list 
        file["intensity_list"] = intensity_list
        file["timestamps"] = et2utc.(time_stamps, "ISOC", 0)
    end
end