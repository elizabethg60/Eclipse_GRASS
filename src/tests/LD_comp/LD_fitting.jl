using Revise
using MyProject
using PyPlot
using LsqFit
using JLD2
using DataFrames
using CSV

model(x, p) = 1 .- p[1] .* (1 .-x) .- p[2] .* (1 .- x) .^ 2.0

# NL94 info - nm
lambda_nm = [303.327, 310.843, 320.468, 329.897, 349.947, 365.875, 374.086, 390.915, 401.970, 416.319, 427.930, 443.885, 445.125, 457.345, 477.427, 492.905, 519.930, 541.76, 559.95, 579.88, 610.975, 640.97, 669.4, 700.875, 748.71, 811.76, 869.6, 948.85, 1046.6, 1098.95]
a0 = [0.08011, 0.08160, 0.08833, 0.09188, 0.11012, 0.12828, 0.12579, 0.12995, 0.12323, 0.12814, 0.14249, 0.16220, 0.15248, 0.16604, 0.19571, 0.20924, 0.23695, 0.26073, 0.26892, 0.28392, 0.30854, 0.33644, 0.34685, 0.37885, 0.40627, 0.42977, 0.45396, 0.47855, 0.49870, 0.51149]
a1 = [0.70695, 0.71609, 0.77285, 0.92459, 1.02168, 1.04969, 0.85402, 0.91836, 1.08648, 1.19947, 1.28796, 1.24893, 1.38517, 1.38544, 1.30551, 1.30798, 1.29927, 1.27428, 1.34319, 1.36896, 1.3662, 1.30590, 1.37539, 1.25553, 1.22842, 1.25182, 1.25101, 1.19813, 1.21429, 1.19354]
a2 = [0.4991, 0.69685, 0.65382, 0.19604, -0.10924, 0.17482, 0.54601, -0.07566, -0.43974, -0.84407, -1.19564, -0.92165, -1.49615, -1.52275, -1.25845, -1.20411, -1.28034, -1.30352, -1.58427, -1.75998, -1.83572, -1.79238, -2.04425, -1.70908, -1.67877, -1.85164, -2.02958, -1.86296, -2.06976, -2.00174]
a3 = [-0.31080, -0.87703, -1.04647, -0.39546, -0.00055, -1.13371, -1.15048, 0.19149, 0.45912, 1.07201, 1.68603, 0.89978, 1.99886, 2.00232, 1.50626, 1.21505, 1.37760, 1.47085, 1.91271, 2.22154, 2.33221, 2.4504, 2.70493, 2.19647, 2.05535, 2.31949, 2.7541, 2.36939, 2.80703, 2.66936]
a4 = [-0.02177, 0.47008, 0.72921, 0.23599, -0.08688, 1.23882, 0.88928, -0.28712, -0.32759, -0.79537, -1.36658, -0.50148, -1.48155, -1.45969, -1.05472, -0.67196, -0.85054, -0.96618, -1.3135, -1.56074, -1.63082, -1.89979, -1.9429, -1.59554, -1.39972, -1.59101, -2.02287, -1.64367, -2.05247, -1.94981]
a5 = [0.04642, -0.0876, -0.19775, -0.05303, 0.06487, -0.4599, -0.26462, 0.12298, 0.0985, 0.23982, 0.44572, 0.1122, 0.44119, 0.42864, 0.30570, 0.14381, 0.21706, 0.26384, 0.37295, 0.4463, 0.45959, 0.59943, 0.55999, 0.47378, 0.38845, 0.44155, 0.59338, 0.46056, 0.60221, 0.57715]

# mu array for figures
mu_arr = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
mu_zero_arr = range(0.0, step=0.05, stop=1.0)
# GRASS lines - A
wavelength = [5250.2084, 5250.6453, 5379.5734, 5381.0216, 5382.2562, 5383.368, 5432.546, 5432.947, 5434.5232, 5435.8577, 5436.2945, 5436.5875, 5576.0881, 5578.718, 6149.246, 6151.617, 6169.042, 6169.563, 6170.5056, 6173.3344, 6301.5008, 6302.4932]

best_u1 = Vector{Float64}(undef,size(wavelength)...)
best_u2 = Vector{Float64}(undef,size(wavelength)...)
#for each GRASS wavelength determine NL94 + Kostogryz LD at mu array
for lambda in 1:length(wavelength)
    NL_LD = Vector{Float64}(undef,size(mu_arr)...)
    NL_LD_full = Vector{Float64}(undef,size(mu_zero_arr)...)
    Kos_LD = Vector{Float64}(undef,size(mu_arr)...)
    Kos_LD_full = Vector{Float64}(undef,size(mu_zero_arr)...)

    # wavelength index for NL94
    index = findmin(x->abs(x-wavelength[lambda]/10.0), lambda_nm)[2]

    #iterate through mu + collect corresponding NL94 + Kostogryz LD
    for i in 1:length(mu_arr)
        NL_LD[i] = a0[index] + a1[index]*mu_arr[i] + a2[index]*mu_arr[i]^2 + a3[index]*mu_arr[i]^3 + a4[index]*mu_arr[i]^4 + a5[index]*mu_arr[i]^5
        Kos_LD[i], best_u1[lambda], best_u2[lambda] = MyProject.quad_limb_darkening_optical(mu_arr[i], wavelength[lambda] ./ 10)
    end

    p0 = [1.0, 1.0]
    fit_NL = curve_fit(model, mu_arr, NL_LD, p0)
    p_opt_NL = fit_NL.param

    for i in 1:length(mu_zero_arr)
        Kos_LD_full[i] = model(mu_zero_arr[i], [best_u1[lambda], best_u2[lambda]])
        NL_LD_full[i] = model(mu_zero_arr[i], p_opt_NL)
    end

    #figures
    plt.scatter(mu_arr, NL_LD, label = "NL94", color = "b")
    plt.scatter(mu_arr, Kos_LD, label = "Kostogryz", color = "r")
    plt.plot(mu_zero_arr, NL_LD_full, color = "b")
    plt.plot(mu_zero_arr, Kos_LD_full, color = "r")
    plt.title(wavelength[lambda])
    plt.legend()
    plt.text(mu_arr[4], 0.5, "NL94 - GRASS wavelength $(round(lambda_nm[index] - wavelength[lambda]/10.0; digits = 3)) nm")
    index = findmin(x->abs(x-wavelength[lambda] ./ 10), MyProject.Kostogryz_LD_file[!, "wavelength"])[2]
    plt.text(mu_arr[4], 0.55, "Kostogryz - GRASS wavelength $(round(collect(MyProject.Kostogryz_LD_file[!, "wavelength"])[index] - wavelength[lambda]/10.0; digits = 3)) nm")
    plt.xlabel("mu")
    plt.ylabel("relative intensity")
    plt.savefig("src/tests/LD_comp/LD_comp_SSD/$(wavelength[lambda]).png")
    plt.clf()
end

df = DataFrame()
df[!, "wavelength"] = wavelength 
df[!, "u1"] = best_u1
df[!, "u2"] = best_u2
CSV.write("quad_ld_coeff_SSD.csv", df)