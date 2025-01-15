path_october = "/storage/group/ebf11/default/pipeline/neid_solar/data/v1.3/L2/2023/10/15/"
files_and_dirs = readdir(path_october)
files = filter(f -> isfile(joinpath(path_october, f)), files_and_dirs)
files = files[4:length(files)-145]

print(files)
println(length(files))
# final_time = Vector{String}(undef, length(files))

# using Dates
# for i in 1:length(final_time)
#     timestamp_str = split(files[i], '_')[2][1:end-5]
#     timestamp = DateTime(timestamp_str, "yyyymmddTHHMMSS")
#     new_timestamp = timestamp + Millisecond(27500)
#     final_time[i] = string(new_timestamp)
# end
# println(files)
# println(final_time)