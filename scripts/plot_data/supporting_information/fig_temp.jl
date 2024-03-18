# load packages
using CairoMakie
using XLSX
using NCDatasets
using Interpolations

# set plotting backend ::
CairoMakie.activate!(type = "png")
Makie.to_font("Arial")

# set plotting basepath ::
plot_data_basepath = pwd() * "/"  # "/Users/perrindavidson/research/whoi/current/gp15/" 

# load station data ::
stat_data = XLSX.readxlsx(plot_data_basepath * "data/data_pro/doQc/stations/gp15_stations.xlsx")["Sheet1"][:]
stat_vars = stat_data[1, :]
stat_data = Float64.(stat_data[2:end, 1:end-1])
ctd_data = XLSX.readxlsx(plot_data_basepath * "data/data_pro/doQc/gp15/gp15_ctd.xlsx")["Sheet1"][:]
ctd_vars = ctd_data[1, :]
ctd_data = Float64.(ctd_data[2:end, :])

# loop through all stations and make array ::
function make_data(zq)
    num_stats = size(stat_data)[1]
    ctd_temp = zeros(length(zq), num_stats)
    for i in 1 : 1 : num_stats
        j₁ = findall(x -> x == stat_data[i, 1], ctd_data[:, 1])
        j₂ = findall(x -> x == stat_data[i, 2], ctd_data[:, 2])
        j = intersect(j₁, j₂)
        dat = ctd_data[j, :]
        f_temp = linear_interpolation(dat[:, 9], dat[:, 10]; extrapolation_bc = Linear())
        ctd_temp[:, i] .= f_temp(zq)
    end
    return ctd_temp
end

zq = 0 : 1 : 400
ctd_temp = make_data(zq)

# make scatter plot :: add in sample locations and change opacy
function lat_names(x)
    if Int64(x) == 0
        return "$(Int64(x))°"  # EQ
    elseif Int64(x) > 0
        return "$(Int64(x))°N"
    elseif Int64(x) < 0
        return "$(Int64(-x))°S"
    end
end

regions = [42, 7.5, -11]

δ = 3

fig = Figure(figure_padding = 25; size = (2400/δ, 1000/δ))

ax₁ = Axis(fig[1, 1],
           # xlabel = L"$$Latitude", 
           ylabel = "Depth (m)",
           titlesize = 20, 
           subtitlesize = 20, 
           xticklabelsize = 16,
           xaxisposition = :top,
           xtickformat = values -> [lat_names(value) for value in values],
           yticklabelsize = 16,
           yreversed = true)
cb = contourf!(stat_data[:, 4], zq, ctd_temp', colormap = (:GnBu_9, 0.75)) # , levels = 0:0.01:0.15)
Colorbar(fig₃[1, 2], cb, label = "Temperature (°C)")  # , ticks = 0:0.05:0.15)
vlines!(ax₁, regions, linestyle = :dash, color = :black, linewidth = 1)
xlims!(ax₁, stat_data[end, 4], stat_data[1, 4])
ylims!(ax₁, 400, 0)

save("plots/gp15Model/supporting_information/fig_temp.png", fig, px_per_unit = 4)

# end plotting routine
