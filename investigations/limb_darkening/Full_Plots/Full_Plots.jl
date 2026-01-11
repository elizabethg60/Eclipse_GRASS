using PyPlot
using LsqFit
using JLD2
using DataFrames
using CSV
using PyCall
py"""
import sys
sys.path.append(".")
"""
mlines = pyimport("matplotlib.lines")

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
    fit_quadratic = curve_fit(quadratic_model, [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 0.0], append!(collect(Kostogryz_LD_array),0.0), p0)
    p_opt_quadratic = fit_quadratic.param

    p0 = [1.0, 1.0, 1.0, 1.0]
    fit_four_parameter = curve_fit(four_parameter_model, [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 0.0], append!(collect(Kostogryz_LD_array),0.0), p0)
    p_opt_four_parameter = fit_four_parameter.param

    return collect(Kostogryz_LD_array), p_opt_quadratic, p_opt_four_parameter
end

# mu array for figures
mu_arr = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 0.0]
mu_zero_arr = range(0.0, step=0.05, stop=1.0)
# GRASS lines
wavelength = [5250.2084, 5250.6453, 5379.5734, 5381.0216, 5382.2562, 5383.368, 5432.546, 5432.947, 5434.5232, 5435.8577, 5436.2945, 5436.5875, 5576.0881, 5578.718, 6149.246, 6151.617, 6169.042, 6169.563, 6170.5056, 6173.3344, 6301.5008, 6302.4932]
claret_coeff = [0.5311, -0.0545, 0.7301, -0.4053] # V band: 500-700 nm

quadratic_u1 = Vector{Float64}(undef,size(wavelength)...)
quadratic_u2 = Vector{Float64}(undef,size(wavelength)...)
four_parameter_u1 = Vector{Float64}(undef,size(wavelength)...)
four_parameter_u2 = Vector{Float64}(undef,size(wavelength)...)
four_parameter_u3 = Vector{Float64}(undef,size(wavelength)...)
four_parameter_u4 = Vector{Float64}(undef,size(wavelength)...)
#for each GRASS wavelength determine Kostogryz LD at mu array
files_arr = [Kostogryz_300, Kostogryz_SSD]#, Kostogryz_HD]
colors_arr = ["b", "r"]#, "r"]
labels_arr = ["300", "SSD"]#, "HD"]
for j in 9:9#1:length(wavelength)
    for i in 1:length(files_arr)
        Kos_LD_quadratic = Vector{Float64}(undef,size(mu_zero_arr)...)
        Kos_LD_four_parameter = Vector{Float64}(undef,size(mu_zero_arr)...)

        Kos_data, quadratic_coeff, four_parameter_coeff = best_fit(wavelength[j] ./ 10, files_arr[i])

        quadratic_u1[j] = quadratic_coeff[1]
        quadratic_u2[j] = quadratic_coeff[2]
        four_parameter_u1[j] = four_parameter_coeff[1]
        four_parameter_u2[j] = four_parameter_coeff[2]
        four_parameter_u3[j] = four_parameter_coeff[3]
        four_parameter_u4[j] = four_parameter_coeff[4]
    
        for i in 1:length(mu_zero_arr)
            Kos_LD_quadratic[i] = quadratic_model(mu_zero_arr[i], [quadratic_u1[j], quadratic_u2[j]])
            Kos_LD_four_parameter[i] = four_parameter_model(mu_zero_arr[i], [four_parameter_u1[j], four_parameter_u2[j], four_parameter_u3[j], four_parameter_u4[j]])
        end

        #figures
        plt.scatter(mu_arr, append!(Kos_data,0.0), color = colors_arr[i], label = labels_arr[i])
        plt.plot(mu_zero_arr, Kos_LD_quadratic, color = colors_arr[i], label = "Quadratic Law", linestyle="--")
        plt.plot(mu_zero_arr, Kos_LD_four_parameter, color = colors_arr[i], label = "Four Parameter Law", linestyle="-.")
    end
    plt.plot(mu_zero_arr, four_parameter_model(mu_zero_arr, claret_coeff), color = "k", label = "Claret 2000")
    plt.title("Fe I 5434 Å", fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    line_legend_b = mlines.Line2D([], [], color="blue", label="300 Gauss", marker="o")
    line_legend_y = mlines.Line2D([], [], color="r", label="SSD", marker="o")
    # line_legend_r = mlines.Line2D([], [], color="red", label="HD", marker="o")
    line_legend_q = mlines.Line2D([], [], color="k", label="Best Fit Quadratic Law", linestyle="--")
    line_legend_4 = mlines.Line2D([], [], color="k", label="Best Fit Four Parameter Law", linestyle="-.")
    line_legend_c = mlines.Line2D([], [], color="k", label="Claret 2000")
    plt.legend(handles=[line_legend_b, line_legend_y, line_legend_q, line_legend_4, line_legend_c], fontsize=12, frameon=false)
    plt.xlabel("μ", fontsize=12)
    plt.ylabel("Relative Intensity", fontsize=12)
    plt.savefig("$(wavelength[j])_comps_zero_constraint.pdf")
    plt.clf()
end