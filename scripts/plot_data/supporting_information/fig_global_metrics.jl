# load packages
using CairoMakie
using XLSX
using NCDatasets
using ColorSchemes

# set plotting backend ::
CairoMakie.activate!(type = "png")

# set plotting basepath ::
plot_data_basepath = pwd() * "/"  # "/Users/perrindavidson/research/whoi/current/gp15/" 

# load depth data ::
gp15_data = XLSX.readxlsx(plot_data_basepath * "data/sim/gp15Model/modelOutput/depthData.xlsx")["Sheet1"][:]
gp15_data_vars = gp15_data[1, :]
gp15_data = gp15_data[2:end, :]
gp15_data[findall(ismissing, gp15_data)] .= NaN
gp15_data = Float64.(gp15_data)
len_gp15 = size(gp15_data)[1]

# load npp data ::
stat_data = XLSX.readxlsx(plot_data_basepath * "data/data_pro/doQc/gp15/gp15_npp.xlsx")["Sheet1"][:]
stat_data_vars = stat_data[1, :]
stat_data = stat_data[2:end, :]
stat_data[findall(ismissing, stat_data)] .= NaN
stat_data = Float64.(stat_data)

# load epzt data ::
epzt_data = XLSX.readxlsx(plot_data_basepath * "data/data_raw/epzt/epzt_metrics.xlsx")["Sheet1"][:]
epzt_data_vars = epzt_data[1, :]
epzt_data = epzt_data[2:end, :]
epzt_data[findall(ismissing, epzt_data)] .= NaN
epzt_data = Float64.(epzt_data)

# load global data ::
global_data = XLSX.readxlsx(plot_data_basepath * "data/data_raw/global_poc/global_database.xlsx")["Global"][:]
global_data_vars = global_data[1, 4:end]
global_data_regions = global_data[:, 1:3]
global_data = global_data[2:end, 4:end]
global_data[findall(ismissing, global_data)] .= NaN
global_data = Float64.(global_data)

# make arrays ::
ppz_data = vcat(gp15_data[:, 4], epzt_data[:, 5], global_data[:, 5])
t100_data = vcat(gp15_data[:, end-3], epzt_data[:, 7], global_data[:, 8])
ez_ratio_data = vcat(gp15_data[:, end], epzt_data[:, 8], global_data[:, 9])
npp_data = vcat(stat_data[:, 3], epzt_data[:, 6], global_data[:, 7])
gp15_data_regions = gp15_data[:, 1]
gp15_data_regions[findall(x -> x > 1, gp15_data_regions)] .= 2
epzt_data_regions = epzt_data[:, 1] .+ 2
regions_data = vcat(Int64.(round.(1.35 .* vcat(gp15_data_regions, epzt_data_regions), digits = 0)), Int64.(9 .* ones(size(global_data)[1])))

# make scatter plot ::
δ = 1.75

markersizes = log10.(npp_data) .* 15
colors = ColorSchemes.GnBu[regions_data]
allcolors = ColorSchemes.GnBu[unique(regions_data)]

fig₉ = Figure(figure_padding = 25; size = (2200/δ, 800/δ))

ax₁ = Axis(fig₉[1:2, 1],
           xlabel = rich("T", subscript("100"), " (%)"),
           ylabel = rich("E", subscript("z"), " Ratio (%)"),
           xlabelsize = 20, 
           ylabelsize = 20, 
           xticklabelsize = 16,
           yticklabelsize = 16,
           xgridvisible = false,
           ygridvisible = false)
scatter!(ax₁, t100_data[len_gp15+1:end], ez_ratio_data[len_gp15+1:end], color = colors[len_gp15+1:end], markersize = markersizes[len_gp15+1:end], strokewidth = 1, strokecolor = :black)
scatter!(ax₁, t100_data[1:len_gp15], ez_ratio_data[1:len_gp15], color = colors[1:len_gp15], markersize = markersizes[1:len_gp15], strokewidth = 1.5, strokecolor = ColorSchemes.RdBu_9[1])

ax₂ = Axis(fig₉[1:2, 2],
           xlabel = rich("E", subscript("z"), " (m)"),
           ylabel = rich("E", subscript("z"), " Ratio (%)"),
           xlabelsize = 20, 
           ylabelsize = 20, 
           xticklabelsize = 16,
           yticklabelsize = 16,
           xgridvisible = false,
           ygridvisible = false)
scatter!(ax₂, ppz_data[len_gp15+1:end], ez_ratio_data[len_gp15+1:end], color = colors[len_gp15+1:end], markersize = markersizes[len_gp15+1:end], strokewidth = 1, strokecolor = :black)
scatter!(ax₂, ppz_data[1:len_gp15], ez_ratio_data[1:len_gp15], color = colors[1:len_gp15], markersize = markersizes[1:len_gp15], strokewidth = 1.5, strokecolor = ColorSchemes.RdBu_9[1])

xlims!(ax₂, 0, 250)

ax₃ = Axis(fig₉[1:2, 3],
           xlabel = rich("E", subscript("z"), " (m)"),
           ylabel = rich("T", subscript("100"), " (%)"),
           xlabelsize = 20, 
           ylabelsize = 20, 
           xticklabelsize = 16,
           yticklabelsize = 16,
           xgridvisible = false,
           ygridvisible = false)
scatter!(ax₃, ppz_data[len_gp15+1:end], t100_data[len_gp15+1:end], color = colors[len_gp15+1:end], markersize = markersizes[len_gp15+1:end], strokewidth = 1, strokecolor = :black)
scatter!(ax₃, ppz_data[1:len_gp15], t100_data[1:len_gp15], color = colors[1:len_gp15], markersize = markersizes[1:len_gp15], strokewidth = 1.5, strokecolor = ColorSchemes.RdBu_9[1])

xlims!(ax₃, 0, 250)

group_size = [MarkerElement(marker = :circle, color = :black, strokecolor = :transparent, markersize = ms) for ms in (log10.([10, 50, 100]) .* 15)]
# group_color = [PolyElement(color = color, strokecolor = :black, strokewidth = 1) for color in allcolors]
group_color = [MarkerElement(marker = :circle, color = allcolors[1], strokecolor = ColorSchemes.RdBu_9[1], strokewidth = 1.5, markersize = 20),
               MarkerElement(marker = :circle, color = allcolors[2], strokecolor = ColorSchemes.RdBu_9[1], strokewidth = 1.5, markersize = 20),
               MarkerElement(marker = :circle, color = allcolors[3], strokecolor = :black, strokewidth = 1, markersize = 20),
               MarkerElement(marker = :circle, color = allcolors[4], strokecolor = :black, strokewidth = 1, markersize = 20),
               MarkerElement(marker = :circle, color = allcolors[5], strokecolor = :black, strokewidth = 1, markersize = 20),
               MarkerElement(marker = :circle, color = allcolors[6], strokecolor = :black, strokewidth = 1, markersize = 20)]

fig₉[1, 4] = Legend(fig₉, group_color, ["PMT – NPHPZ", "PMT – Other", "EPZT – Shelf", "EPZT – Offshore", "EPZT – Gyre", "Global"], ["Region"], tellheight = true, tellwidth = true)
fig₉[2, 4] = Legend(fig₉, group_size, string.([10, 50, 100]), [rich("NPP (mmol C m", superscript("-2"), " d", superscript("-1"), ")")], tellheight = true, tellwidth = true)

# linkaxes!(ax₁, ax₂, ax₃)

for (ax, label) in zip([ax₁, ax₂, ax₃], ["A", "B", "C"])
    text!(
        ax, 0, 1,
        text = label,
        font = :bold,
        align = (:left, :top),
        offset = (250, -2),
        space = :relative,
        fontsize = 24
    )
end

save("plots/gp15Model/supporting_information/fig_global_metrics.png", fig₉, px_per_unit = 4)

# end plotting routine
