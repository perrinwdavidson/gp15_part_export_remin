# load packages
using CairoMakie
using XLSX
using NCDatasets
using ColorSchemes

# set plotting backend ::
CairoMakie.activate!(type = "png")

# set basepaths ::
plot_data_basepath = pwd() * "/"  # "/Users/perrindavidson/research/whoi/current/gp15/" 

# load station data ::
poc_data = XLSX.readxlsx(plot_data_basepath * "data/data_raw/global_poc/global_poc_flux.xlsx")["Sheet1"][:]
poc_vars = poc_data[1, :]
poc_data = Float64.(poc_data[2:end, :])

# load in gp15 100 m data (for comparison) ::
poc_data_gp15 = XLSX.readxlsx(plot_data_basepath * "data/sim/gp15Model/modelOutput/depthData.xlsx")["Sheet1"][:]
poc_gp15_vars = poc_data_gp15[1, :]
poc_data_gp15 = poc_data_gp15[2:end, :]
poc_data_gp15[findall(ismissing, poc_data_gp15)] .= NaN
poc_data_gp15 = Float64.(poc_data_gp15)

# load in stations ::
stat_data = XLSX.readxlsx(plot_data_basepath * "data/data_pro/doQc/stations/gp15_stations.xlsx")["Sheet1"][:]
stat_vars = stat_data[1, :]
stat_data = Float64.(stat_data[2:end, 1:end-1])

# quality control compilation data ::
poc_data = poc_data[findall(x -> x > 0, poc_data[:, 4]), :]
len_poc_dat = size(poc_data)[1]

# make ppz flux data ::
gp15_ppz_flux = poc_data_gp15[:, 21]
gp15_ppz_flux_upwell = poc_data_gp15[:, 23]
gp15_ppz_flux[findall(!isnan, gp15_ppz_flux_upwell)] .= gp15_ppz_flux_upwell[findall(!isnan, gp15_ppz_flux_upwell)]
gp15_ppz_flux = hcat(stat_data[:, 4], stat_data[:, 3], stat_data[:, 5], gp15_ppz_flux)
gp15_ppz_flux = gp15_ppz_flux[findall(!isnan, gp15_ppz_flux[:, 4]), :]

# add in data ::
poc_data = vcat(poc_data, gp15_ppz_flux)

# load bathymetry data ::
Δ = 20
bath_nc = NCDataset(plot_data_basepath * "data/data_raw/gebco/gebco_2023_n65.0_s-25.0_w-180.0_e-60.0.nc")
x_bath = bath_nc["lon"][1:Δ:end]
y_bath = bath_nc["lat"][1:Δ:end]
bath = bath_nc["elevation"][1:Δ:end, 1:Δ:end]

# set text ::
select_stations = [10, 18, 29, 39]
regions = [42, 7.5, -11]

function lat_names(x)
    if Int64(x) == 0
        return "$(Int64(x))°"  # EQ
    elseif Int64(x) > 0
        return "$(Int64(x))°N"
    elseif Int64(x) < 0
        return "$(Int64(-x))°S"
    end
end

# make scatter plot ::
δ = 1.25
fig₈ = Figure(figure_padding = 25; size = (1200/δ, 1000/δ))
ax₁ = Axis(fig₈[1, 1],
           xticklabelsize = 16,
           xtickformat = values -> ["$(Int64(-1 * value))°W" for value in values],
           yticklabelsize = 16,
           ytickformat = values -> [lat_names(value) for value in values])

contourf!(x_bath, y_bath, bath, colormap = (:oslo, 0.25), levels = range(minimum(bath), 0, length = 20), extendhigh = :gray)  # GnBu_9. but also consider: deep, devon

hlines!(ax₁, regions, xmax = [0.75, 1, 1], linestyle = :dash, color = :black, linewidth = 1)

contourf!(x_bath, y_bath, bath, colormap = (:oslo, 0.25), levels = range(0, 0, length = 1), extendhigh = :gray)  # GnBu_9. but also consider: deep, devon

colors = log10.(poc_data[:, 4])
markersizes = (poc_data[:, 3] ./ maximum(poc_data[:, 3])) .* 30

sblabel = rich("log", subscript("10"), rich(" P"), subscript("POC"), "(z) (dpm m", superscript("-2"), " d", superscript("-1"), ")")

scatter!(Point2f.(poc_data[len_poc_dat+1:end, 2], poc_data[len_poc_dat+1:end, 1]), color = colors[len_poc_dat+1:end], markersize = markersizes[len_poc_dat+1:end], colormap = cgrad(:GnBu_9, 12, categorical = true), strokewidth = 1.5, strokecolor = ColorSchemes.RdBu_9[1], colorrange = (-1, 2))

sb = scatter!(Point2f.(poc_data[1:len_poc_dat, 2], poc_data[1:len_poc_dat, 1]), color = colors[1:len_poc_dat], markersize = markersizes[1:len_poc_dat], colormap = cgrad(:GnBu_9, 12, categorical = true), strokewidth = 1, strokecolor = :black, colorrange = (-1, 2))

Colorbar(fig₈[1, 2], sb, label = sblabel, ticks = -1:0.5:2, height = 600/δ, tellheight = false)

group_type = [MarkerElement(marker = :circle, color = :white, strokecolor = sc, strokewidth = 1.5, markersize = 12.48) for sc in [ColorSchemes.RdBu_9[1], :black]]
axislegend(ax₁, group_type, ["This Study", "Previous Studies"], position = :lb)  # tellheight = true)

group_size = [MarkerElement(marker = :circle, color = :black, strokecolor = :transparent, markersize = ms) for ms in (([50, 100, 200] ./ maximum(poc_data[:, 3])) .* 30)]
axislegend(ax₁, group_size, string.([50, 100, 200]), ["Sampling depth, z (m)"], position = :rt)  # tellheight = true)

xlims!(ax₁, -180, -60)
ylims!(ax₁, -25, 65)
save("plots/gp15Model/manuscript/fig8.png", fig₈, px_per_unit = 4)

# end plotting routine
