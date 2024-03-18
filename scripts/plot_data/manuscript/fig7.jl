# load packages
using CairoMakie
using XLSX
using NCDatasets
using Interpolations

# set plotting backend ::
CairoMakie.activate!(type = "png")

# set plotting basepath ::
plot_data_basepath = pwd() * "/"  # "/Users/perrindavidson/research/whoi/current/gp15/" 

# load station data ::
stat_data = XLSX.readxlsx(plot_data_basepath * "data/sim/gp15Model/modelOutput/depthData.xlsx")["Sheet1"][:]
stat_vars = stat_data[1, :]
stat_data = stat_data[2:end, :]
stat_data[findall(ismissing, stat_data)] .= NaN
stat_data = Float64.(stat_data)

# make scatter plot :: 
function lat_names(x)
    if Int64(x) == 0
        return "$(Int64(x))°" # EQ
    elseif Int64(x) > 0
        return "$(Int64(x))°N"
    elseif Int64(x) < 0
        return "$(Int64(-x))°S"
    end
end

δ = 2

regions = [42, 7.5, -11]

fig₇ = Figure(figure_padding = 25; size = (2100/δ, 1200/δ))

ax₁ = Axis(fig₇[1:2, 1],
           ylabel = rich("P", subscript("Th"), "(100) (dpm m", superscript("-2"), " d", superscript("-1"), ")"),
           xticklabelsize = 16,
           xtickformat = values -> [lat_names(value) for value in values],
           xticklabelsvisible = false,
           yticklabelsize = 16,
           xticks = -20:10:65,
           xgridvisible = false,
           ygridvisible = false)

lat₀ = stat_data[end:-1:1, 3]
lat = collect(LinRange(round(minimum(lat₀), digits = 0), round(maximum(lat₀), digits = 0), 33)) .+ 1

vlines!(ax₁, regions, linestyle = :dash, color = :black, linewidth = 1)

flux_100₀ = stat_data[end:-1:1, 5]
flux_100_err₀ = stat_data[end:-1:1, 6] 
upwell_flux_100₀ = stat_data[end:-1:1, 7]
upwell_flux_100_err₀ = stat_data[end:-1:1, 8] 

barplot!(ax₁, lat, flux_100₀, color = ColorSchemes.GnBu[7], gap = 0.2, strokecolor = :black, strokewidth = 1)
barplot!(ax₁, lat, flux_100₀, fillto = upwell_flux_100₀, color = ColorSchemes.GnBu[9], gap = 0.2, strokecolor = :black, strokewidth = 1)

barplot!(ax₁, lat, flux_100₀ .- flux_100_err₀, fillto = flux_100₀ .+ flux_100_err₀, color = ColorSchemes.GnBu[4], gap = 0.5, strokecolor = :black, strokewidth = 1)
barplot!(ax₁, lat, upwell_flux_100₀ .- upwell_flux_100_err₀, fillto = upwell_flux_100₀ .+ upwell_flux_100_err₀, color = ColorSchemes.GnBu[5], gap = 0.5, strokecolor = :black, strokewidth = 1)

hlines!(ax₁, 0, color = :black, linewidth = 1)

ax₂ = Axis(fig₇[3:4, 1],
           ylabel = rich("P", subscript("Th"), "(E", subscript("z"), ") (dpm m", superscript("-2"), " d", superscript("-1"), ")"),
           xticklabelsize = 16,
           xtickformat = values -> [lat_names(value) for value in values],
           yticklabelsize = 16,
           xticks = -20:10:65,
           xgridvisible = false,
           ygridvisible = false)

flux_ppz₀ = stat_data[end:-1:1, 9]
flux_ppz_err₀ = stat_data[end:-1:1, 10] 
upwell_flux_ppz₀ = stat_data[end:-1:1, 11]
upwell_flux_ppz_err₀ = stat_data[end:-1:1, 12] 

vlines!(ax₂, regions, linestyle = :dash, color = :black, linewidth = 1)

barplot!(ax₂, lat, flux_ppz₀, color = ColorSchemes.GnBu[7], gap = 0.2, strokecolor = :black, strokewidth = 1)
barplot!(ax₂, lat, flux_ppz₀, fillto = upwell_flux_ppz₀, color = ColorSchemes.GnBu[9], gap = 0.2, strokecolor = :black, strokewidth = 1)

barplot!(ax₂, lat, flux_ppz₀ .- flux_ppz_err₀, fillto = flux_ppz₀ .+ flux_ppz_err₀, color = ColorSchemes.GnBu[4], gap = 0.5, strokecolor = :black, strokewidth = 1)
barplot!(ax₂, lat, upwell_flux_ppz₀ .- upwell_flux_ppz_err₀, fillto = upwell_flux_ppz₀ .+ upwell_flux_ppz_err₀, color = ColorSchemes.GnBu[5], gap = 0.5, strokecolor = :black, strokewidth = 1)

scatterlines!(ax₂, lat, stat_data[end:-1:1, end-1], color = :black, linestyle = :dash)

hlines!(ax₂, 0, color = :black, linewidth = 1)

flux = [PolyElement(color = ColorSchemes.GnBu_9[7], strokecolor = :black, strokewidth = 1, points = Point2f[(0, 0.15), (0, 0.85), (1, 0.85), (1, 0.15)]),
        PolyElement(color = ColorSchemes.GnBu_9[4], strokecolor = :black, strokewidth = 1, points = Point2f[(0.25, 0), (0.25, 1), (0.75, 1), (0.75, 0)])]
upwell = [PolyElement(color = ColorSchemes.GnBu_9[9], strokecolor = :black, strokewidth = 1, points = Point2f[(0, 0.15), (0, 0.85), (1, 0.85), (1, 0.15)]),
          PolyElement(color = ColorSchemes.GnBu_9[5], strokecolor = :black, strokewidth = 1, points = Point2f[(0.25, 0), (0.25, 1), (0.75, 1), (0.75, 0)])]
r100 = [LineElement(color = :black, linestyle = :dash), MarkerElement(color = :black, marker = :circle)] 

Legend(fig₇[2:3, 2], [flux, upwell, r100], ["Flux", "Upwell Correction", rich("R", subscript("100"))])  # , tellheight = false, tellwidth = false) # , labelsize = 30, patchsize = (35, 35))

linkaxes!(ax₁, ax₂)

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

save("plots/gp15Model/manuscript/fig7.png", fig₇, px_per_unit = 4)

# end plotting routine
