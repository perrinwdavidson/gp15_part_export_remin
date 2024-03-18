# load packages
using CairoMakie
using XLSX
using NCDatasets
using ColorSchemes

# set plotting backend ::
CairoMakie.activate!(type = "png")

# set plotting basepath ::
plot_data_basepath = pwd() * "/"  # "/Users/perrindavidson/research/whoi/current/gp15/" 

# load station data ::
stat_data = XLSX.readxlsx(plot_data_basepath * "data/sim/interpData/upwelling/w_ecco.xlsx")["Sheet1"][:]
stat_vars = stat_data[1, :]
stat_data = Float64.(stat_data[2:end, :])
ppz_data = XLSX.readxlsx(plot_data_basepath * "data/data_pro/doQc/stations/gp15_stations.xlsx")["Sheet1"][:]
ppz_vars = ppz_data[1, :]
ppz_data = Float64.(ppz_data[2:end, 1:end-1])

# get data ::
stat_100 = stat_data[findall(x -> x == 100, stat_data[:, 4]), :]
stat_ppz = stat_data[findall(x -> x in ppz_data[:, 5], stat_data[:, 4]), :]

# make scatter plot ::
δ = 3

s2d = 60 * 60 * 24

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

fig₆ = Figure(figure_padding = 25; size = (2400/δ, 1000/δ))

ax₁ = Axis(fig₆[1, 1],
           ylabel = rich("Upwelling velocity, w (m d", superscript("-1"), ")"),
           xticklabelsize = 16,
           xtickformat = values -> [lat_names(value) for value in values],
           yticklabelsize = 16,
           xgridvisible = false,
           ygridvisible = false)

hlines!(ax₁, 0, linestyle = :dash, color = :black)
vlines!(ax₁, regions, linestyle = :dash, color = :black, linewidth = 1)

band!(ax₁, 
      Point2.(Float64.(stat_100[:, 3]), Float64.(stat_100[:, 6] .- stat_100[:, 8]) .* s2d),
      Point2.(Float64.(stat_100[:, 3]), Float64.(stat_100[:, 6] .+ stat_100[:, 8]) .* s2d),
      color = (ColorSchemes.GnBu_9[4], 0.25))
scatterlines!(ax₁, stat_100[:, 3], stat_100[:, 6] .* s2d, color = ColorSchemes.GnBu[4])

band!(ax₁, 
      Point2.(Float64.(stat_ppz[:, 3]), Float64.(stat_ppz[:, 6] .- stat_ppz[:, 8]) .* s2d),
      Point2.(Float64.(stat_ppz[:, 3]), Float64.(stat_ppz[:, 6] .+ stat_ppz[:, 8]) .* s2d),
      color = (ColorSchemes.GnBu_9[9], 0.25))
scatterlines!(ax₁, stat_ppz[:, 3], stat_ppz[:, 6] .* s2d, color = ColorSchemes.GnBu[9])

w_ppz = [PolyElement(color = (ColorSchemes.GnBu_9[9], 0.25), strokecolor = :transparent, points = Point2f[(0, 0.25), (0, 0.75), (1, 0.75), (1, 0.25)]),
         LineElement(color = ColorSchemes.GnBu[9])]
w_100 = [PolyElement(color = (ColorSchemes.GnBu_9[4], 0.25), strokecolor = :transparent, points = Point2f[(0, 0.25), (0, 0.75), (1, 0.75), (1, 0.25)]),
         LineElement(color = ColorSchemes.GnBu_9[4])]

axislegend(ax₁, [w_ppz, w_100], [rich("E", subscript("z")), "100 m"], position = :rt)  

save("plots/gp15Model/manuscript/fig6.png", fig₆, px_per_unit = 4)

# end plotting routine
