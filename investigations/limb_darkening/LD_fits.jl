using PyPlot
using LsqFit
using JLD2
using DataFrames
using CSV

Kostogryz_300 = DataFrame(CSV.File("Kostogryz_LD_300.csv", header = ["wavelength", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"]))
Kostogryz_SSD = DataFrame(CSV.File("Kostogryz_LD_SSD.csv", header = ["wavelength", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"]))
Kostogryz_HD = DataFrame(CSV.File("Kostogryz_LD_HD.csv", header = ["wavelength", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"]))

quadratic_model(x, p) = 1 .- p[1] .* (1 .-x) .- p[2] .* (1 .- x) .^ 2.0
four_parameter_model(x, p) = 1 .- p[1] .* (1 .-x.^0.5) .- p[2] .* (1 .-x) .- p[3] .* (1 .-x.^1.5) .- p[4] .* (1 .-x.^2.0)

function best_fit(wavelength::T, file) where T
    """
    limb darkening prescription for optical based on mu angle  
    """

    index = findmin(x->abs(x-wavelength), file[!, "wavelength"])[2]
    Kostogryz_LD_array = file[index, ["0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"]]

    p0 = [1.0, 1.0]
    fit_quadratic = curve_fit(quadratic_model, [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], append!(collect(Kostogryz_LD_array)), p0)
    p_opt_quadratic = fit_quadratic.param

    p0 = [1.0, 1.0, 1.0, 1.0]
    fit_four_parameter = curve_fit(four_parameter_model, [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], append!(collect(Kostogryz_LD_array)), p0)
    p_opt_four_parameter = fit_four_parameter.param

    return collect(Kostogryz_LD_array), p_opt_quadratic, p_opt_four_parameter
end

# GRASS lines
wavelength = [5250.2084, 5250.6453, 5379.5734, 5381.0216, 5382.2562, 5383.368, 5432.546, 5432.947, 5434.5232, 5435.8577, 5436.2945, 5436.5875, 5576.0881, 5578.718, 6149.246, 6151.617, 6169.042, 6169.563, 6170.5056, 6173.3344, 6301.5008, 6302.4932, 15652]

quadratic_u1 = Vector{Float64}(undef,size(wavelength)...)
quadratic_u2 = Vector{Float64}(undef,size(wavelength)...)
four_parameter_u1 = Vector{Float64}(undef,size(wavelength)...)
four_parameter_u2 = Vector{Float64}(undef,size(wavelength)...)
four_parameter_u3 = Vector{Float64}(undef,size(wavelength)...)
four_parameter_u4 = Vector{Float64}(undef,size(wavelength)...)
#for each GRASS wavelength determine Kostogryz LD at mu array
for lambda in 1:length(wavelength)
    Kos_data, quadratic_coeff, four_parameter_coeff = best_fit(wavelength[lambda] ./ 10, Kostogryz_SSD)

    quadratic_u1[lambda] = quadratic_coeff[1]
    quadratic_u2[lambda] = quadratic_coeff[2]
    four_parameter_u1[lambda] = four_parameter_coeff[1]
    four_parameter_u2[lambda] = four_parameter_coeff[2]
    four_parameter_u3[lambda] = four_parameter_coeff[3]
    four_parameter_u4[lambda] = four_parameter_coeff[4]
end

df = DataFrame()
df[!, "wavelength"] = wavelength 
df[!, "u1"] = quadratic_u1
df[!, "u2"] = quadratic_u2
df[!, "4u1"] = four_parameter_u1
df[!, "4u2"] = four_parameter_u2
df[!, "4u3"] = four_parameter_u3
df[!, "4u4"] = four_parameter_u4
CSV.write("LD_coeff_SSD.csv", df)