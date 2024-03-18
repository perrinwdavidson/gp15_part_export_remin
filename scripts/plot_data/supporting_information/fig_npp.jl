# load packages
using CairoMakie
using XLSX
using DelimitedFiles
using NCDatasets
using ColorSchemes

# set plotting backend ::
CairoMakie.activate!(type = "png")
Makie.to_font("Arial")

# set plotting basepath ::
plot_data_basepath = pwd() * "/"  # "/Users/perrindavidson/research/whoi/current/gp15/" 

# load bathymetry data ::
Δ = 20
bath_nc = NCDataset(plot_data_basepath * "data/data_raw/gebco/gebco_2023_n65.0_s-25.0_w-180.0_e-60.0.nc")
x_bath = bath_nc["lon"][1:Δ:end]
y_bath = bath_nc["lat"][1:Δ:end]
bath = bath_nc["elevation"][1:Δ:end, 1:Δ:end]

# load npp data ::
x_npp = readdlm(plot_data_basepath * "data/data_pro/readData/cbpm/npp_lon.csv", ',', Float64)
y_npp = readdlm(plot_data_basepath * "data/data_pro/readData/cbpm/npp_lat.csv", ',', Float64)
npp = readdlm(plot_data_basepath * "data/data_pro/calcNpp/npp_mean.csv", ',', Float64)

# load in stations ::
stat_data = XLSX.readxlsx(plot_data_basepath * "data/data_pro/doQc/stations/gp15_stations.xlsx")["Sheet1"][:]
stat_vars = stat_data[1, :]
stat_data = Float64.(stat_data[2:end, 1:end-1])

# load flux data ::
poc_data_gp15 = XLSX.readxlsx(plot_data_basepath * "data/sim/gp15Model/modelOutput/depthData.xlsx")["Sheet1"][:]
poc_gp15_vars = poc_data_gp15[1, :]
poc_data_gp15 = poc_data_gp15[2:end, :]
poc_data_gp15[findall(ismissing, poc_data_gp15)] .= NaN
poc_data_gp15 = Float64.(poc_data_gp15)
gp15_ppz_flux = poc_data_gp15[:, 21]
gp15_ppz_flux_upwell = poc_data_gp15[:, 23]
gp15_ppz_flux[findall(!isnan, gp15_ppz_flux_upwell)] .= gp15_ppz_flux_upwell[findall(!isnan, gp15_ppz_flux_upwell)]
gp15_ppz_flux = hcat(stat_data[:, 4], stat_data[:, 3], stat_data[:, 5], gp15_ppz_flux)
gp15_ppz_flux = gp15_ppz_flux[findall(!isnan, gp15_ppz_flux[:, 4]), :]
poc_data = gp15_ppz_flux

gp15_ppz_flux = poc_data_gp15[:, 33]
gp15_ppz_flux_upwell = poc_data_gp15[:, 35]
gp15_ppz_flux[findall(!isnan, gp15_ppz_flux_upwell)] .= gp15_ppz_flux_upwell[findall(!isnan, gp15_ppz_flux_upwell)]
gp15_ppz_flux = hcat(stat_data[:, 4], stat_data[:, 3], stat_data[:, 5], gp15_ppz_flux)
gp15_ppz_flux = gp15_ppz_flux[findall(!isnan, gp15_ppz_flux[:, 4]), :]
pn_data = gp15_ppz_flux

# define plotting function ::
function lat_names(x)
    if Int64(x) == 0
        return "$(Int64(x))°"  # "EQ"
    elseif Int64(x) > 0
        return "$(Int64(x))°N"
    elseif Int64(x) < 0
        return "$(Int64(-x))°S"
    end
end

function reg_names(x)
    if x == 52
        return "NPHPZ"  
    elseif x == 24.75
        return "NPG"
    elseif x == -1.75
        return "EP"
    elseif x == -16.5
        return "SPG"
    end
end

regions = [42, 7.5, -11]

# make scatter plot ::
δ = 1.25
fig = Figure(figure_padding = 25; size = (700/δ*1.5, 1000/δ*1.15))
ax₁ = Axis(fig[1, 1:2],
           xticklabelsize = 16,
           xtickformat = values -> ["$(Int64(-1 * value))°W" for value in values],
           yticklabelsize = 16,
           yticks = [52, 24.75, -1.75, -16.5],
           ytickformat = values -> [reg_names(value) for value in values],
           yticklabelrotation = π/2,
           xgridvisible = false,
           ygridvisible = false)

cb = contourf!(x_npp[1, :], y_npp[:, 1], npp', colormap = (:GnBu, 0.75), levels = 0:10:100, extendlow = ColorSchemes.GnBu[1], extendhigh = ColorSchemes.GnBu[end])  # , extendhigh = ColorSchemes.GnBu[end]) # , levels = 0:0.25:3.5)
Colorbar(fig[0, 2:3], cb, label = rich("NPP (mmol C m", superscript("-2"), " d", superscript("-1"), ")"), ticks = 0:20:100, vertical = false)  # ticks = 0:20:100
hlines!(ax₁, regions, linestyle = :dash, color = :black, linewidth = 1)

colors = poc_data[:, 4]
markersizes = (poc_data[:, 3] ./ maximum(poc_data[:, 3])) .* 30
sb = scatter!(Point2f.(poc_data[:, 2], poc_data[:, 1]), color = colors[:], markersize = markersizes[:], colorrange = (0, 2.5), strokewidth = 1, strokecolor = :black, colormap = cgrad(:GnBu_9, 10, categorical = true))
Colorbar(fig[2, 1:2], sb, label = rich("P", subscript("POC"), "(z) (mmol C m", superscript("-2"), " d", superscript("-1"), ")"), ticks = 0:0.5:2.5, vertical = false, flipaxis = false)

contourf!(x_bath, y_bath, bath, colormap = (:oslo, 0.75), levels = 0:0, extendhigh = :gray)  # also consider: deep, devon

xlims!(ax₁, -180, -120)
ylims!(ax₁, -22, 62)

ax₂ = Axis(fig[1, 3:4],
           xticklabelsize = 16,
           yticks = -20 : 20 : 60, 
           xtickformat = values -> ["$(Int64(-1 * value))°W" for value in values],
           yticklabelsize = 16,
           ytickformat = values -> [lat_names(value) for value in values],
           yaxisposition = :right,
           xgridvisible = false,
           ygridvisible = false)

contourf!(x_npp[1, :], y_npp[:, 1], npp', colormap = (:GnBu, 0.75), levels = 0:10:100, extendlow = ColorSchemes.GnBu[1], extendhigh = ColorSchemes.GnBu[end])  # , extendhigh = ColorSchemes.GnBu[end]) # , levels = 0:0.25:3.5)
hlines!(ax₂, regions, linestyle = :dash, color = :black, linewidth = 1)

colors = pn_data[:, 4]
markersizes = (pn_data[:, 3] ./ maximum(pn_data[:, 3])) .* 30
sb = scatter!(Point2f.(pn_data[:, 2], pn_data[:, 1]), color = colors[:], markersize = markersizes[:], colorrange = (0, 0.25), strokewidth = 1, strokecolor = :black, colormap = cgrad(:GnBu_9, 10, categorical = true))
Colorbar(fig[2, 3:4], sb, label = rich("P", subscript("PN"), "(z) (mmol N m", superscript("-2"), " d", superscript("-1"), ")"), ticks = 0:0.05:0.25, vertical = false, flipaxis = false)

contourf!(x_bath, y_bath, bath, colormap = (:oslo, 0.75), levels = 0:0, extendhigh = :gray)  # also consider: deep, devon

xlims!(ax₂, -180, -120)
ylims!(ax₂, -22, 62)

group_size = [MarkerElement(marker = :circle, color = :black, strokecolor = :transparent, markersize = ms) for ms in (([50, 100, 200] ./ maximum(poc_data[:, 3])) .* 30)]
axislegend(ax₂, group_size, string.([50, 100, 200]), ["Sampling depth (m)"], position = :rt)  # tellheight = true)

colgap!(fig.layout, 60)

for (ax, label) in zip([ax₁, ax₂], ["A", "B"])
    text!(
        ax, 0, 1,
        text = label,
        font = :bold,
        align = (:left, :top),
        offset = (4, -2),
        space = :relative,
        fontsize = 24
    )
end

save("plots/gp15Model/supporting_information/fig_npp.png", fig, px_per_unit = 4)

# end plotting routine
