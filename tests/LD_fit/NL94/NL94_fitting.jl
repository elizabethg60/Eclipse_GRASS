using PyPlot
using GRASS
using LsqFit
using DataFrames
using CSV

LD_model(x, p) = 1 .- p[1] .* (1 .-x) .- p[2] .* (1 .- x) .^ 2.0

NL94_coeff = CSV.read("/storage/home/efg5335/work/GRASS/data/NL94_LD.csv", DataFrame)
lambda_nm = NL94_coeff.lambda_nm
a0 = NL94_coeff.a0
a1 = NL94_coeff.a1
a2 = NL94_coeff.a2
a3 = NL94_coeff.a3
a4 = NL94_coeff.a4
a5 = NL94_coeff.a5

# GRASS lines - A
wavelength = [5250.2084, 5250.6453, 5379.5734, 5381.0216, 5382.2562, 5383.368, 5432.546, 5432.947, 5434.5232, 5435.8577, 5436.2945, 5436.5875, 5576.0881, 5578.718, 6149.246, 6151.617, 6169.042, 6169.563, 6170.5056, 6173.3344, 6301.5008, 6302.4932]

mu_limb = range(0.0, step=0.01, stop=1.0)

#for each GRASS wavelength determine NL94 at mu array
for lambda in 1:length(wavelength)
    # wavelength index for NL94
    index = findmin(x->abs(x-wavelength[lambda]/10.0), lambda_nm)[2]

    # mu array for figures
    mu_arr = [0.12, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    for i in 1:length(mu_arr)
        mu_zero_arr = range(mu_arr[i], step=0.01, stop=1.0)
        NL_LD_full = Vector{Float64}(undef,size(mu_zero_arr)...)

        for j in 1:length(mu_zero_arr)    
            NL_LD_full[j] = a0[index] + a1[index]*mu_zero_arr[j] + a2[index]*mu_zero_arr[j]^2 + a3[index]*mu_zero_arr[j]^3 + a4[index]*mu_zero_arr[j]^4 + a5[index]*mu_zero_arr[j]^5
        end

        p0 = [1.0, 1.0]
        fit_NL = curve_fit(LD_model, mu_zero_arr, NL_LD_full, p0)
        p_opt_NL = fit_NL.param

        NL_LD_full_fit = LD_model(mu_limb, p_opt_NL)

        #figures
        plt.scatter(mu_zero_arr, NL_LD_full)
        plt.plot(mu_limb, NL_LD_full_fit)
    end
    plt.title(wavelength[lambda])
    plt.text(mu_arr[4], 0.5, "NL94 - GRASS wavelength $(round(lambda_nm[index] - wavelength[lambda]/10.0; digits = 3)) nm")
    plt.xlabel("mu")
    plt.ylabel("relative intensity")
    plt.savefig("/storage/home/efg5335/work/Eclipse_GRASS/tests/LD_fit/NL94/$(wavelength[lambda]).png")
    plt.clf()
end