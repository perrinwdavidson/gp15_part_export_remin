# ---------------------------------------------------------
# import packages ::
using NCDatasets
using LoopVectorization
using NaNStatistics
using Statistics
using Distances
using Random
using CairoMakie
using LsqFit
using Printf
using JLD2
using LinearAlgebra
using DelimitedFiles

# set plotting background ::
CairoMakie.activate!(type = "png")

# define structs ::
struct Inputs
    coords::Array{Float64, 2}
    data::Array{Float64, 1}
    bin_space::AbstractRange
    bin_time::AbstractRange
end

# ---------------------------------------------------------
# define decorrelation function ::
function calc_decorr(inputs::Inputs)

    # get dimensions ::
    n = size(inputs.data)[1]
    n_space = length(inputs.bin_space)
    n_time = length(inputs.bin_time)

    # preallocate output ::
    autocov, autovar₁, autovar₂, counts = zeros(n_space, n_time), zeros(n_space, n_time), zeros(n_space, n_time), zeros(n_space, n_time)

    # loop through and calculate ::
    for i in 1 : 1 : (n - 1)

        # get current point ::
        p₁ = inputs.coords[i, :]
        d₁ = inputs.data[i]

        for j in (i + 1) : 1 : n
        
                # get next point ::
                p₂ = inputs.coords[j, :]
                d₂ = inputs.data[j]

                # get distance ::
                distₘ = haversine(p₁[1:2], p₂[1:2])  # [m], great circle distance, using volumetric mean radius r = 6,371,000
                # distₘ = euclidean(p₁[1:2], p₂[1:2]) * 111000  # [m], euclidean distance assuming small angle approx
                distₜ = abs(p₁[3] - p₂[3]) * (60 * 60 * 24) # [s]

                # get index ::
                _, idxₘ = histcountindices([distₘ], inputs.bin_space)
                _, idxₜ = histcountindices([distₜ], inputs.bin_time)
                idxₘ, idxₜ = idxₘ[1], idxₜ[1]
                if idxₘ == 0 || idxₜ == 0
                    continue
                end

                # calculate and store correlation ::
                autocov[idxₘ, idxₜ] += d₁ * d₂
                autovar₁[idxₘ, idxₜ] += d₁ ^ 2 
                autovar₂[idxₘ, idxₜ] += d₂ ^ 2 
                counts[idxₘ, idxₜ] += 1

        end

        # print out ::
        if i % 1E2 == 0
            println("Done with $(i)th row of data.")
        end

    end

    # calculate autocorrelation ::
    return autocov ./ sqrt.(autovar₁ .* autovar₂), autocov, autovar₁,  autovar₂, counts

end

# data function to fit ::
function fit_func(x, p)
    a, b, c = p
    return (a .* exp.(-(x ./ b) .^ 2)) .+ c
end

# data fitting function ::
function data_fit(p0, func, xvar, yvar, wt, residual, fit_type_name)
	"""
	purpose :: to fit data (xvar, yvar) to a function func(p0) using the levenburg-marquadt algorithm
	inputs ::
		[1] p0 - the initial parameters to start fitting func() to
		[2] func - function to fit
		[3] xvar - the independent sample data
		[4] yvar - the dependent sample data
		[5] wt - the weights for least squares
		[6] residual - the residual fuction to minimize
		[7] fit_type_name - the name of the function that is being fit for printing
	outputs ::
		[1] pf - the final parameter from the fitting routine
		[2] pferr - the error in the estimation
		[3] chisq - the chisquared of the estimation
		[4] dof - the degrees of freedom from the estimation
	"""
	# initialize fitting object with inputs ::
	fit = curve_fit(func, xvar, yvar, wt, p0);

	# get fit values:
	pf = fit.param;

	# compute the covariance matrix by finding the inverse of the Jacobian times its transpose:
	cov = inv(transpose(fit.jacobian) * fit.jacobian);

	# calculate chi-squared sum:
	chisq = sum(residual(pf, func, xvar, yvar) .^ 2);

	# calculate degrees of freedom of estimate:
	dof = length(xvar) - length(pf);

	# calculate reduced chi-squared:
	red_chisq = chisq / dof;

	# calculate error in parameter estimate by squaring diagonal elements of the covariance matrix:
	pferr = diag(sqrt.(Diagonal(cov)));

	# print out fit options:
	println();
	println(fit_type_name, " Function");
	@printf("Converged with chi-squared %.2f \n", chisq);
	@printf("Number of degrees of freedom, dof = %.2f \n", dof);
	@printf("Reduced chi-squared %2f \n", red_chisq);
	println();
	Columns = ["Parameter #", "Initial guess values:", "Best fit values:", "Uncertainties in the best fit values:"]
	@printf("%11s | %24s | %24s | %24s \n", Columns[1], Columns[2], Columns[3], Columns[4]);
	for num in 1 : 1 : length(pf)
		@printf("%11i | %24.3e | %24.3e | %24.3e \n", num, p0[num], pf[num], pferr[num]);
      	end

	# return outputs:
	return pf, pferr, chisq, dof

end

# residual function ::
function residual(p, func, xvar, yvar)
	"""
	purpose ::
	inputs ::
		[1] p - parameter for function func
		[2] func - function to estimate the residual for
		[3] xvar - the independent data to evaluate the residual at
		[4] yvar - the dependent data to evaluate the residual at
	outputs ::
		[1] residual - the residual from the data
	"""
	return (func(xvar, p) .- yvar) ./ var(yvar);

end

# nan mean function ::
function calc_mean(autocorr)
    time_autocorr_init = Float64[]
    time_autocorr = Float64[]
    for i in 1 : 1 : size(autocorr)[2]
        dat = autocorr[:, i]
        μ_dat = mean(dat[findall(!isnan, dat)])
        time_autocorr_init = append!(time_autocorr_init, μ_dat)
        time_autocorr = append!(time_autocorr, μ_dat)
    end
    space_autocorr = Float64[]
    for i in 1 : 1 : size(autocorr)[1]
        dat = autocorr[i, :]
        μ_dat = mean(dat[findall(!isnan, dat)])
        space_autocorr = append!(space_autocorr, μ_dat)
    end
    return time_autocorr_init, time_autocorr, space_autocorr
end

# ---------------------------------------------------------
# load data ::
npp_data = NCDataset("/Users/perrindavidson/research/whoi/current/gp15/data/data_pro/readData/cbpm/npp.nc")
npp₀ = npp_data["npp"][:, :, :]
lon = npp_data["lon"][:]
lat = npp_data["lat"][:]
time = npp_data["time"][:]

# make data ::
lon_grid, lat_grid, time_grid = (ones(length(lat)) .* lon') .* ones(1, 1, length(time)), (lat .* ones(length(lon))') .* ones(1, 1, length(time)), ones(length(lat), length(lon)) .* reshape(time, 1, 1, length(time))
npp = hcat(lon_grid[:], lat_grid[:], time_grid[:], npp₀[:])
npp = npp[findall(!isnan, npp[:, 4]), :]
μ = mean(npp[:, 4])
npp = hcat(npp, npp[:, 4] .- μ)

# select data ::
idx = sort(randsubseq(collect(1 : 1 : size(npp)[1]), 0.0005))
coords = npp[idx, 1:3]
data = npp[idx, 5]  # 4 is raw, 5 is detrended

# make bins (Sumata (2017, Ocean Sci.)) ::
bin_space = (0 : 10 : 5E2) .* 1000  # [m]
bin_time = (0 : 16 : 4E2) .* (60 * 60 * 24)  # [s]

# make inputs ::
inputs = Inputs(coords,
                data,
                bin_space,
                bin_time)

# run function ::
autocorr, autocov, autovar₁, autovar₂, counts = calc_decorr(inputs)

# remove points where counts less than 1000 ::
autocorr_init = autocorr
autocorr[counts .< 1000] .= NaN  # see pg. 171 from Sumata (2017, Ocean Sci.)

# calculate autocorrelation ::
time_autocorr_init, time_autocorr, space_autocorr = calc_mean(autocorr)

# make grid vectors ::
bs = bin_space ./ 1000
bt = bin_time ./ (60 * 60 * 24)
bt_plot = collect(bt)
bt_plot_init = collect(bt)
bs_plot = collect(bs)

# given previous work, we know the bounds to fit ::
# n.b.: we do not remove points from spatial comp. as good fit 
# given that length scale of pacific ≫ 500 [km] (about 1400 [km], Marshall and Speer) 
# which would cause a secondary peak. the annual cycle is dealt with now.
time_autocorr[bt_plot .> 150] .= NaN  # see pg. 172 from Sumata (2017, Ocean Sci.)

# make plotting vectors without nans ::
bt_plot = bt_plot[findall(x -> !isnan(x), time_autocorr[:])]
ta_plot = time_autocorr[findall(x -> !isnan(x), time_autocorr[:])]
bs_plot = bs_plot[findall(x -> !isnan(x), space_autocorr[:])]
sa_plot = space_autocorr[findall(x -> !isnan(x), space_autocorr[:])]

# fit spatial data ::
p₀ = [0.5, 400, 0.3]
wt_err = diagm(ones(length(sa_plot))) ./ var(sa_plot)
fit_params_space, fit_params_err_space, fit_χ², fit_ν = data_fit(p₀,
                                                                 fit_func,
                                                                 bs_plot,
                                                                 sa_plot,
                                                                 wt_err,
                                                                 residual,
                                                                 "Correlation")
fitted_data_space = fit_func(bs_plot, fit_params_space)

# fit temporal data :: 
p₀ = [0.6, 150, 0.2]
wt_err = diagm(ones(length(ta_plot))) ./ var(ta_plot)
fit_params_time, fit_params_err_time, fit_χ², fit_ν = data_fit(p₀,
                                                               fit_func,
                                                               bt_plot,
                                                               ta_plot,
                                                               wt_err,
                                                               residual,
                                                               "Correlation")
fitted_data_time = fit_func(bt_plot_init, fit_params_time)

# plot contour data ::
f = Figure(; size = (1800, 2400))

ax₁ = Axis(f[1, 1], 
          title = L"$$Autocorrelation",
          xlabel = L"$$Space Lag [km]",
          ylabel = L"$$Time Lag [d]")
ρ = contourf!(bs, bt, autocorr)  
Colorbar(f[1, 2], ρ)

ax₂ = Axis(f[1, 3], 
          title = L"$$Covariance",
          xlabel = L"$$Space Lag [km]",
          ylabel = L"$$Time Lag [d]")
σ²₁₂ = contourf!(bs, bt, autocov)  
Colorbar(f[1, 4], σ²₁₂)

ax₃ = Axis(f[2, 1], 
          title = L"$$Variance 1",
          xlabel = L"$$Space Lag [km]",
          ylabel = L"$$Time Lag [d]")
σ²₁ = contourf!(bs, bt, autovar₁)  
Colorbar(f[2, 2], σ²₁)

ax₄ = Axis(f[2, 3], 
          title = L"$$Variance 2",
          xlabel = L"$$Space Lag [km]",
          ylabel = L"$$Time Lag [d]")
σ²₂ = contourf!(bs, bt, autovar₂)  
Colorbar(f[2, 4], σ²₂)

ax₅ = Axis(f[3, 1], 
          title = L"$$Counts",
          xlabel = L"$$Space Lag [km]",
          ylabel = L"$$Time Lag [d]")
c = contourf!(bs, bt, counts)  
Colorbar(f[3, 2], c)

save("plots/calcNpp/decorrelation/autocorrelation.png", f, px_per_unit = 4)

# calculate and plot line data ::
f = Figure(; size = (800, 400))

ax₆ = Axis(f[1, 1], 
          title = L"$$Time Mean Autocorrelation",
          ylabel = L"$$Autocorrelation",
          xlabel = L"$$Time Lag [d]")
scatterlines!(ax₆, bt_plot_init, time_autocorr_init[:]) 
lines!(ax₆, bt_plot_init, fitted_data_time, color = "red", linestyle = :dash)

ax₇ = Axis(f[1, 2], 
          title = L"$$Space Mean Autocorrelation",
          xlabel = L"$$Space Lag [km]",
          ylabel = L"$$Autocorrelation")
scatterlines!(ax₇, bs_plot, sa_plot, label = L"$$Data") 
lines!(ax₇, bs_plot, fitted_data_space, color = "red", linestyle = :dash, label = L"$$Fit") 
axislegend()

save("plots/calcNpp/decorrelation/autocorrelation_mean.png", f, px_per_unit = 4)

# save data ::
time_bins_array = bt_plot
spatial_bins_array = bs_plot
time_autocorr_mean = ta_plot
spatial_autocorr_mean = sa_plot
save_object("/Users/perrindavidson/research/whoi/current/gp15/data/sim/calcNpp/decorrelation/ecco_autocorr.jld2", autocorr)
save_object("/Users/perrindavidson/research/whoi/current/gp15/data/sim/calcNpp/decorrelation/ecco_autocov.jld2", autocov)
save_object("/Users/perrindavidson/research/whoi/current/gp15/data/sim/calcNpp/decorrelation/ecco_autovar1.jld2", autovar₁)
save_object("/Users/perrindavidson/research/whoi/current/gp15/data/sim/calcNpp/decorrelation/ecco_autovar2.jld2", autovar₂)
save_object("/Users/perrindavidson/research/whoi/current/gp15/data/sim/calcNpp/decorrelation/ecco_counts.jld2", counts)
save_object("/Users/perrindavidson/research/whoi/current/gp15/data/sim/calcNpp/decorrelation/ecco_time_bins_array.jld2", time_bins_array)
save_object("/Users/perrindavidson/research/whoi/current/gp15/data/sim/calcNpp/decorrelation/ecco_spatial_bins_array.jld2", spatial_bins_array)
save_object("/Users/perrindavidson/research/whoi/current/gp15/data/sim/calcNpp/decorrelation/ecco_time_autocorr_mean.jld2", time_autocorr_mean)
save_object("/Users/perrindavidson/research/whoi/current/gp15/data/sim/calcNpp/decorrelation/ecco_spatial_autocorr_mean.jld2", spatial_autocorr_mean)
save_object("/Users/perrindavidson/research/whoi/current/gp15/data/sim/calcNpp/decorrelation/ecco_spatial_autocorr_mean_fit.jld2", fitted_data_space)
save_object("/Users/perrindavidson/research/whoi/current/gp15/data/sim/calcNpp/decorrelation/ecco_spatial_autocorr_mean_fit_params.jld2", fit_params_space)
save_object("/Users/perrindavidson/research/whoi/current/gp15/data/sim/calcNpp/decorrelation/ecco_spatial_autocorr_mean_fit_params_err.jld2", fit_params_err_space)
save_object("/Users/perrindavidson/research/whoi/current/gp15/data/sim/calcNpp/decorrelation/ecco_temporal_autocorr_mean_fit.jld2", fitted_data_time)
save_object("/Users/perrindavidson/research/whoi/current/gp15/data/sim/calcNpp/decorrelation/ecco_temporal_autocorr_mean_fit_params.jld2", fit_params_time)
save_object("/Users/perrindavidson/research/whoi/current/gp15/data/sim/calcNpp/decorrelation/ecco_temporal_autocorr_mean_fit_params_err.jld2", fit_params_err_time)

# write out length scales ::
writedlm("/Users/perrindavidson/research/whoi/current/gp15/data/sim/calcNpp/decorrelation/spatial_decorrelation.csv", 
         hcat(fit_params_space, fit_params_err_space), 
         ',')
writedlm("/Users/perrindavidson/research/whoi/current/gp15/data/sim/calcNpp/decorrelation/temporal_decorrelation.csv", 
         hcat(fit_params_time, fit_params_err_time), 
         ',')

# write out data ::
npp_array = NCDataset("/Users/perrindavidson/research/whoi/current/gp15/data/sim/calcNpp/decorrelation/npp_matrix.nc",
                      "c")
defDim(npp_array, "serial_number", size(npp)[1])
defDim(npp_array, "observation", 1)
npp_lon = defVar(npp_array, "longitude", Float64, ("serial_number", "observation"))
npp_lat = defVar(npp_array, "latitude", Float64, ("serial_number", "observation"))
npp_time = defVar(npp_array, "time", Float64, ("serial_number", "observation"))
npp_data = defVar(npp_array, "npp", Float64, ("serial_number", "observation"))
npp_lon[:] = npp[:, 1]
npp_lat[:] = npp[:, 2]
npp_time[:] = npp[:, 3]
npp_data[:] = npp[:, 4]
close(npp_array)

# end subroutine
