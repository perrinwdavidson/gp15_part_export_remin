# load packages
using CairoMakie
using XLSX
using NCDatasets
using ColorSchemes

# set plotting backend ::
CairoMakie.activate!(type = "pdf")

# set plotting basepath ::
plot_data_basepath = pwd() * "/"  # "/Users/perrindavidson/research/whoi/current/gp15/" 

# load depth data ::
gp15_data = XLSX.readxlsx(plot_data_basepath * "data/sim/gp15Model/modelOutput/gp15_flux.xlsx")["Sheet1"][:]
gp15_data_vars = gp15_data[1, :]
gp15_data_vars = vcat(gp15_data_vars[1:8], gp15_data_vars[10:end])  # remove dates
gp15_data = gp15_data[2:end, :]
gp15_data[findall(ismissing, gp15_data)] .= NaN
gp15_data = Float64.(hcat(gp15_data[:, 1:8], gp15_data[:, 10:end]))

# load npp data ::
stat_data = XLSX.readxlsx(plot_data_basepath * "data/data_pro/doQc/gp15/gp15_npp.xlsx")["Sheet1"][:]
stat_data_vars = stat_data[1, :]
stat_data = stat_data[2:end, :]
stat_data[findall(ismissing, stat_data)] .= NaN
stat_data = Float64.(stat_data)

# make arrays ::
region_data = gp15_data[:, 2]
depth_data = gp15_data[:, 1]
pn_data = gp15_data[:, 28]
poc_data = gp15_data[:, 18]
pn_stat_data = gp15_data[:, 4]
npp_data0 = stat_data[:, 3]
npp_stat_data0 = stat_data[:, 1]

# interpolate npp to all stations ::
npp_data = NaN .* zeros(size(pn_data))
function fill_npp!(npp_data)
    for i in 1 : 1 : size(npp_stat_data0)[1]
        npp_data[findall(j -> pn_stat_data[j] == npp_stat_data0[i], 1 : 1 : size(pn_stat_data)[1])] .= npp_data0[i]
    end
    return npp_data
end
npp_data = fill_npp!(npp_data)

# sort by npp ::
data = hcat(pn_data, depth_data, npp_data, region_data, poc_data)
# data = data[sortperm(data[:, 3]), :]
# data = data[end:-1:1, :]

# make scatter plot ::
δ = 1.75

markersizes = log10.(npp_data) .* 15
# colors = ColorSchemes.GnBu[regions_data]
colors = ColorSchemes.GnBu[data[:, 3]]

fig₉ = Figure(figure_padding = 25; size = (2400/δ, 1000/δ .* 2))

ax₁ = Axis(fig₉[3, 1],
           ylabel = "Depth (m)",
           xlabelsize = 20, 
           ylabelsize = 20, 
           xticklabelsize = 16,
           yticklabelsize = 16,
           xaxisposition = :top, 
           xgridvisible = false,
           ygridvisible = false)
idx_region = 1
idx_plot = findall(i -> data[i, 4] == idx_region, 1 : 1 : size(data)[1])
colors = npp_data
sb = scatter!(Point2f.(log10.(data[idx_plot, 1]), data[idx_plot, 2]), color = colors[idx_plot], markersize = 20, strokewidth = 1, strokecolor = :black, colormap = cgrad(:GnBu_9, 12, categorical = true), colorrange = (0, 120))

ax₂ = Axis(fig₉[3, 2],
           xlabelsize = 20, 
           ylabelsize = 20, 
           xticklabelsize = 16,
           yticklabelsize = 16,
           yticklabelsvisible = :false,
           xaxisposition = :top, 
           xgridvisible = false,
           ygridvisible = false)
idx_region = 2
idx_plot = findall(i -> data[i, 4] == idx_region, 1 : 1 : size(data)[1])
colors = npp_data
sb = scatter!(Point2f.(log10.(data[idx_plot, 1]), data[idx_plot, 2]), color = colors[idx_plot], markersize = 20, strokewidth = 1, strokecolor = :black, colormap = cgrad(:GnBu_9, 12, categorical = true), colorrange = (0, 120))

ax₃ = Axis(fig₉[3, 3],
           xlabelsize = 20, 
           ylabelsize = 20, 
           xticklabelsize = 16,
           yticklabelsize = 16,
           yticklabelsvisible = :false,
           xaxisposition = :top, 
           xgridvisible = false,
           ygridvisible = false)
idx_region = 3
idx_plot = findall(i -> data[i, 4] == idx_region, 1 : 1 : size(data)[1])
colors = npp_data
sb = scatter!(Point2f.(log10.(data[idx_plot, 1]), data[idx_plot, 2]), color = colors[idx_plot], markersize = 20, strokewidth = 1, strokecolor = :black, colormap = cgrad(:GnBu_9, 12, categorical = true), colorrange = (0, 120))

ax₄ = Axis(fig₉[3, 4],
           xlabelsize = 20, 
           ylabelsize = 20, 
           xticklabelsize = 16,
           yticklabelsvisible = :false,
           yticklabelsize = 16,
           xaxisposition = :top, 
           xgridvisible = false,
           ygridvisible = false)
idx_region = 4
idx_plot = findall(i -> data[i, 4] == idx_region, 1 : 1 : size(data)[1])
colors = npp_data
sb = scatter!(Point2f.(log10.(data[idx_plot, 1]), data[idx_plot, 2]), color = colors[idx_plot], markersize = 20, strokewidth = 1, strokecolor = :black, colormap = cgrad(:GnBu_9, 12, categorical = true), colorrange = (0, 120))

Colorbar(fig₉[3, 5], sb, label = rich("NPP (mmol C m", superscript("-2"), " d", superscript("-1"), ")"), ticks = 0 : 20 : 120)  # rich("P", subscript("PN"), "(z) (mmol N m", superscript("-2"), " d", superscript("-1"), ")"), ticks = 0:0.05:0.25, vertical = false, flipaxis = false)

ylims!(ax₁, 420, 0)
ylims!(ax₂, 420, 0)
ylims!(ax₃, 420, 0)
ylims!(ax₄, 420, 0)
linkyaxes!(ax₁, ax₂, ax₃, ax₄)

Label(fig₉[2, 2:3], rich("log", subscript("10"), " LSF PN:", superscript("234"), "Th (μmol N dpm", superscript("-1"), ")"), fontsize = 24)

for (ax, label) in zip([ax₁, ax₂, ax₃, ax₄], ["NPHPZ", "NPG", "EP", "SPG"])
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

ax₅ = Axis(fig₉[1, 1],
           ylabel = "Depth (m)",
           xlabelsize = 20, 
           ylabelsize = 20, 
           xticklabelsize = 16,
           yticklabelsize = 16,
           xaxisposition = :top, 
           xgridvisible = false,
           ygridvisible = false)
idx_region = 1
idx_plot = findall(i -> data[i, 4] == idx_region, 1 : 1 : size(data)[1])
colors = npp_data
sb = scatter!(Point2f.(log10.(data[idx_plot, 5]), data[idx_plot, 2]), color = colors[idx_plot], markersize = 20, strokewidth = 1, strokecolor = :black, colormap = cgrad(:GnBu_9, 12, categorical = true), colorrange = (0, 120))

ax₆ = Axis(fig₉[1, 2],
           xlabelsize = 20, 
           ylabelsize = 20, 
           xticklabelsize = 16,
           yticklabelsize = 16,
           yticklabelsvisible = :false,
           xaxisposition = :top, 
           xgridvisible = false,
           ygridvisible = false)
idx_region = 2
idx_plot = findall(i -> data[i, 4] == idx_region, 1 : 1 : size(data)[1])
colors = npp_data
sb = scatter!(Point2f.(log10.(data[idx_plot, 5]), data[idx_plot, 2]), color = colors[idx_plot], markersize = 20, strokewidth = 1, strokecolor = :black, colormap = cgrad(:GnBu_9, 12, categorical = true), colorrange = (0, 120))

ax₇ = Axis(fig₉[1, 3],
           xlabelsize = 20, 
           ylabelsize = 20, 
           xticklabelsize = 16,
           yticklabelsize = 16,
           yticklabelsvisible = :false,
           xaxisposition = :top, 
           xgridvisible = false,
           ygridvisible = false)
idx_region = 3
idx_plot = findall(i -> data[i, 4] == idx_region, 1 : 1 : size(data)[1])
colors = npp_data
sb = scatter!(Point2f.(log10.(data[idx_plot, 5]), data[idx_plot, 2]), color = colors[idx_plot], markersize = 20, strokewidth = 1, strokecolor = :black, colormap = cgrad(:GnBu_9, 12, categorical = true), colorrange = (0, 120))

ax₈ = Axis(fig₉[1, 4],
           xlabelsize = 20, 
           ylabelsize = 20, 
           xticklabelsize = 16,
           yticklabelsvisible = :false,
           yticklabelsize = 16,
           xaxisposition = :top, 
           xgridvisible = false,
           ygridvisible = false)
idx_region = 4
idx_plot = findall(i -> data[i, 4] == idx_region, 1 : 1 : size(data)[1])
colors = npp_data
sb = scatter!(Point2f.(log10.(data[idx_plot, 5]), data[idx_plot, 2]), color = colors[idx_plot], markersize = 20, strokewidth = 1, strokecolor = :black, colormap = cgrad(:GnBu_9, 12, categorical = true), colorrange = (0, 120))

Colorbar(fig₉[1, 5], sb, label = rich("NPP (mmol C m", superscript("-2"), " d", superscript("-1"), ")"), ticks = 0 : 20 : 120)  # rich("P", subscript("PN"), "(z) (mmol N m", superscript("-2"), " d", superscript("-1"), ")"), ticks = 0:0.05:0.25, vertical = false, flipaxis = false)

ylims!(ax₅, 420, 0)
ylims!(ax₆, 420, 0)
ylims!(ax₇, 420, 0)
ylims!(ax₈, 420, 0)
linkyaxes!(ax₅, ax₆, ax₇, ax₈)

xlims!(ax₁, -6.2, -2.4)
linkxaxes!(ax₁, ax₂, ax₃, ax₄, ax₅, ax₆, ax₇, ax₈)

Label(fig₉[0, 2:3], rich("log", subscript("10"), " LSF POC:", superscript("234"), "Th (μmol C dpm", superscript("-1"), ")"), fontsize = 24)

for (ax, label) in zip([ax₅, ax₆, ax₇, ax₈], ["NPHPZ", "NPG", "EP", "SPG"])
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

# group_size = [MarkerElement(marker = :circle, color = :black, strokecolor = :transparent, markersize = ms) for ms in (log10.([10, 50, 100]) .* 15)]
# group_color = [PolyElement(color = color, strokecolor = :black, strokewidth = 1) for color in allcolors]
# group_color = [MarkerElement(marker = :circle, color = allcolors[1], strokecolor = ColorSchemes.RdBu_9[1], strokewidth = 1.5, markersize = 20),
#                MarkerElement(marker = :circle, color = allcolors[2], strokecolor = ColorSchemes.RdBu_9[1], strokewidth = 1.5, markersize = 20),
#                MarkerElement(marker = :circle, color = allcolors[3], strokecolor = :black, strokewidth = 1, markersize = 20),
#                MarkerElement(marker = :circle, color = allcolors[4], strokecolor = :black, strokewidth = 1, markersize = 20),
#                MarkerElement(marker = :circle, color = allcolors[5], strokecolor = :black, strokewidth = 1, markersize = 20)]
# 
# fig₉[1, 2] = Legend(fig₉, group_color, ["PMT – NPHPZ", "PMT – Other", "EPZT – Shelf", "EPZT – Offshore", "EPZT – Gyre"], ["Region"], tellheight = true, tellwidth = true)
# fig₉[2, 2] = Legend(fig₉, group_size, string.([10, 50, 100]), [rich("NPP (mmol C m", superscript("-2"), " d", superscript("-1"), ")")], tellheight = true, tellwidth = true)

# linkaxes!(ax₁, ax₂, ax₃)

save("plots/gp15Model/supporting_information/fig_pn.pdf", fig₉)

# end plotting routine
