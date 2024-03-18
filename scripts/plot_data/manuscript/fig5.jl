# load packages
using CairoMakie
using XLSX
using NCDatasets
using Interpolations

# set plotting backend ::
CairoMakie.activate!(type = "png")
Makie.to_font("Arial")

# set basepaths ::
plot_data_basepath = pwd() * "/"  # "/Users/perrindavidson/research/whoi/current/gp15/" 

# load station data ::
stat_data = XLSX.readxlsx(plot_data_basepath * "data/data_pro/delineateRegions/regionalData/delineate_data.xlsx")["Sheet1"][:]
stat_vars = stat_data[1, :]
stat_data = stat_data[2:end, :]

# make scatter plot :: 
function lat_names(x)
    if Int64(x) == 0
        return "$(Int64(x))°"
    elseif Int64(x) > 0
        return "$(Int64(x))°N"
    elseif Int64(x) < 0
        return "$(Int64(-x))°S"
    end
end

δ = 1.75

fig₅ = Figure(figure_padding = 25; size = (1800/δ, 1200/δ))

ax₁ = Axis(fig₅[1, 1],
           ylabel = "Temperature (°C)",
           xticklabelsize = 16,
           xtickformat = values -> [lat_names(value) for value in values],
           xticklabelsvisible = false,
           yticklabelsize = 16)
lat = Float64.(stat_data[findall(!ismissing, stat_data[:, 4]), 2])
dat = Float64.(stat_data[findall(!ismissing, stat_data[:, 4]), 4])

x0 = round(stat_data[end, 2], digits = 0) - 2
xend = round(stat_data[1, 2], digits = 0) + 2
regions = [42, 7.5, -11]

y0 = round(minimum(dat), digits = 0) - 5
yend = round(maximum(dat), digits = 0) + 5

r1 = Float64.([xend y0; xend yend])
r2 = Float64.([regions[1] y0; regions[1] yend])
r3 = Float64.([regions[2] y0; regions[2] yend])
r4 = Float64.([regions[3] y0; regions[3] yend])
r5 = Float64.([x0 y0; x0 yend])

band!(ax₁, Point2.(r1[:, 1], r1[:, 2]), Point2.(r2[:, 1], r2[:, 2]), color = (ColorSchemes.GnBu_9[2], 0.5))
band!(ax₁, Point2.(r2[:, 1], r2[:, 2]), Point2.(r3[:, 1], r3[:, 2]), color = (ColorSchemes.GnBu_9[4], 0.5))
band!(ax₁, Point2.(r3[:, 1], r3[:, 2]), Point2.(r4[:, 1], r4[:, 2]), color = (ColorSchemes.GnBu_9[6], 0.5))
band!(ax₁, Point2.(r4[:, 1], r4[:, 2]), Point2.(r5[:, 1], r5[:, 2]), color = (ColorSchemes.GnBu_9[8], 0.5))

scatter!(ax₁, lat, dat, marker = :rect, color = :black, markersize = 15)
xlims!(ax₁, x0, xend)
ylims!(ax₁, y0, yend)

ax₂ = Axis(fig₅[2, 1],
           ylabel = "Bulk Particles (volts)",
           xticklabelsize = 16,
           xtickformat = values -> [lat_names(value) for value in values],
           xticklabelsvisible = false,
           # yreversed = true,
           yticklabelsize = 16)
lat = Float64.(stat_data[findall(!ismissing, stat_data[:, 5]), 2])
dat = Float64.(stat_data[findall(!ismissing, stat_data[:, 5]), 5])

y0 = round(minimum(dat), digits = 2) - 0.04
yend = round(maximum(dat), digits = 2) + 0.04

r1 = Float64.([xend y0; xend yend])
r2 = Float64.([regions[1] y0; regions[1] yend])
r3 = Float64.([regions[2] y0; regions[2] yend])
r4 = Float64.([regions[3] y0; regions[3] yend])
r5 = Float64.([x0 y0; x0 yend])

band!(ax₂, Point2.(r1[:, 1], r1[:, 2]), Point2.(r2[:, 1], r2[:, 2]), color = (ColorSchemes.GnBu_9[2], 0.5))
band!(ax₂, Point2.(r2[:, 1], r2[:, 2]), Point2.(r3[:, 1], r3[:, 2]), color = (ColorSchemes.GnBu_9[4], 0.5))
band!(ax₂, Point2.(r3[:, 1], r3[:, 2]), Point2.(r4[:, 1], r4[:, 2]), color = (ColorSchemes.GnBu_9[6], 0.5))
band!(ax₂, Point2.(r4[:, 1], r4[:, 2]), Point2.(r5[:, 1], r5[:, 2]), color = (ColorSchemes.GnBu_9[8], 0.5))

scatter!(ax₂, lat, dat, marker = :rect, color = :black, markersize = 15)
xlims!(ax₂, x0, xend)
ylims!(ax₂, y0, yend)

ax₃ = Axis(fig₅[3, 1],
           ylabel = rich("PO", subscript("4"), superscript("3-"), " (μM)"),
           xticklabelsize = 16,
           xtickformat = values -> [lat_names(value) for value in values],
           yticklabelsize = 16)
lat = Float64.(stat_data[findall(!ismissing, stat_data[:, 6]), 2])
dat = Float64.(stat_data[findall(!ismissing, stat_data[:, 6]), 6])

y0 = round(minimum(dat), digits = 2) - 0.2
yend = round(maximum(dat), digits = 2) + 0.2

r1 = Float64.([xend y0; xend yend])
r2 = Float64.([regions[1] y0; regions[1] yend])
r3 = Float64.([regions[2] y0; regions[2] yend])
r4 = Float64.([regions[3] y0; regions[3] yend])
r5 = Float64.([x0 y0; x0 yend])

band!(ax₃, Point2.(r1[:, 1], r1[:, 2]), Point2.(r2[:, 1], r2[:, 2]), color = (ColorSchemes.GnBu_9[2], 0.5))
band!(ax₃, Point2.(r2[:, 1], r2[:, 2]), Point2.(r3[:, 1], r3[:, 2]), color = (ColorSchemes.GnBu_9[4], 0.5))
band!(ax₃, Point2.(r3[:, 1], r3[:, 2]), Point2.(r4[:, 1], r4[:, 2]), color = (ColorSchemes.GnBu_9[6], 0.5))
band!(ax₃, Point2.(r4[:, 1], r4[:, 2]), Point2.(r5[:, 1], r5[:, 2]), color = (ColorSchemes.GnBu_9[8], 0.5))

scatter!(ax₃, lat, dat, marker = :rect, color = :black, markersize = 15)
xlims!(ax₃, x0, xend)
ylims!(ax₃, y0, yend)

ax₄ = Axis(fig₅[1, 2],
           ylabel = rich("NO", subscript("3"), superscript("2-"), " (μM)"),
           xticklabelsize = 16,
           xtickformat = values -> [lat_names(value) for value in values],
           xticklabelsvisible = false,
           yticklabelsize = 16)
lat = Float64.(stat_data[findall(!ismissing, stat_data[:, 7]), 2])
dat = Float64.(stat_data[findall(!ismissing, stat_data[:, 7]), 7])

y0 = round(minimum(dat), digits = 2) - 2
yend = round(maximum(dat), digits = 2) + 2

r1 = Float64.([xend y0; xend yend])
r2 = Float64.([regions[1] y0; regions[1] yend])
r3 = Float64.([regions[2] y0; regions[2] yend])
r4 = Float64.([regions[3] y0; regions[3] yend])
r5 = Float64.([x0 y0; x0 yend])

band!(ax₄, Point2.(r1[:, 1], r1[:, 2]), Point2.(r2[:, 1], r2[:, 2]), color = (ColorSchemes.GnBu_9[2], 0.5))
band!(ax₄, Point2.(r2[:, 1], r2[:, 2]), Point2.(r3[:, 1], r3[:, 2]), color = (ColorSchemes.GnBu_9[4], 0.5))
band!(ax₄, Point2.(r3[:, 1], r3[:, 2]), Point2.(r4[:, 1], r4[:, 2]), color = (ColorSchemes.GnBu_9[6], 0.5))
band!(ax₄, Point2.(r4[:, 1], r4[:, 2]), Point2.(r5[:, 1], r5[:, 2]), color = (ColorSchemes.GnBu_9[8], 0.5))

scatter!(ax₄, lat, dat, marker = :rect, color = :black, markersize = 15)
xlims!(ax₄, x0, xend)
ylims!(ax₄, y0, yend)

labels = ["NPHZ", "NPG", "EP", "SPG"]
elements = [PolyElement(polycolor = (ColorSchemes.GnBu_9[i], 0.5), strokecolor = :black, strokewidth = 1) for i in [2, 4, 6, 8]]

Legend(fig₅[1, 3], elements, labels)

ax₅ = Axis(fig₅[2, 2],
           ylabel = "Particulate N (μM)",
           xticklabelsize = 16,
           xtickformat = values -> [lat_names(value) for value in values],
           xticklabelsvisible = false,
           yticklabelsize = 16)
lat = Float64.(stat_data[findall(!ismissing, stat_data[:, 8]), 2])
dat = Float64.(stat_data[findall(!ismissing, stat_data[:, 8]), 8])

y0 = round(minimum(dat), digits = 0) 
yend = round(maximum(dat), digits = 1) + 0.1

r1 = Float64.([xend y0; xend yend])
r2 = Float64.([regions[1] y0; regions[1] yend])
r3 = Float64.([regions[2] y0; regions[2] yend])
r4 = Float64.([regions[3] y0; regions[3] yend])
r5 = Float64.([x0 y0; x0 yend])

band!(ax₅, Point2.(r1[:, 1], r1[:, 2]), Point2.(r2[:, 1], r2[:, 2]), color = (ColorSchemes.GnBu_9[2], 0.5))
band!(ax₅, Point2.(r2[:, 1], r2[:, 2]), Point2.(r3[:, 1], r3[:, 2]), color = (ColorSchemes.GnBu_9[4], 0.5))
band!(ax₅, Point2.(r3[:, 1], r3[:, 2]), Point2.(r4[:, 1], r4[:, 2]), color = (ColorSchemes.GnBu_9[6], 0.5))
band!(ax₅, Point2.(r4[:, 1], r4[:, 2]), Point2.(r5[:, 1], r5[:, 2]), color = (ColorSchemes.GnBu_9[8], 0.5))

scatter!(ax₅, lat, dat, marker = :rect, color = :black, markersize = 15)
xlims!(ax₅, x0, xend)
ylims!(ax₅, y0, yend)

ax₆ = Axis(fig₅[3, 2],
           ylabel = "Pigment Composition (%)",
           xticklabelsize = 16,
           xtickformat = values -> [lat_names(value) for value in values],
           yticklabelsize = 16)

y0 = 0
yend = 1

r1 = Float64.([xend y0; xend yend])
r2 = Float64.([regions[1] y0; regions[1] yend])
r3 = Float64.([regions[2] y0; regions[2] yend])
r4 = Float64.([regions[3] y0; regions[3] yend])
r5 = Float64.([x0 y0; x0 yend])

band!(ax₆, Point2.(r1[:, 1], r1[:, 2]), Point2.(r2[:, 1], r2[:, 2]), color = (ColorSchemes.GnBu_9[2], 0.5))
band!(ax₆, Point2.(r2[:, 1], r2[:, 2]), Point2.(r3[:, 1], r3[:, 2]), color = (ColorSchemes.GnBu_9[4], 0.5))
band!(ax₆, Point2.(r3[:, 1], r3[:, 2]), Point2.(r4[:, 1], r4[:, 2]), color = (ColorSchemes.GnBu_9[6], 0.5))
band!(ax₆, Point2.(r4[:, 1], r4[:, 2]), Point2.(r5[:, 1], r5[:, 2]), color = (ColorSchemes.GnBu_9[8], 0.5))

lat₀ = Float64.(stat_data[findall(!ismissing, stat_data[:, 9]), 2])
lat = collect(LinRange(round(minimum(lat₀), digits = 0), round(maximum(lat₀), digits = 0), length(lat₀))) 
lat = lat[end:-1:1]

grp = Int64.(vcat(ones(size(lat)), 5 .* ones(size(lat)), 9 .* ones(size(lat))))
lat = vcat(lat, lat, lat)
dat = vcat(Float64.(stat_data[findall(!ismissing, stat_data[:, 9]), 9]), Float64.(stat_data[findall(!ismissing, stat_data[:, 9]), 10]), Float64.(stat_data[findall(!ismissing, stat_data[:, 9]), 11]))

barplot!(ax₆, lat, dat, stack = grp, color = ColorSchemes.RdBu_9[grp], gap = 0.5, strokecolor = :black, strokewidth = 1)

labels = ["Micro", "Pico", "Nano"]
elements = [PolyElement(polycolor = ColorSchemes.RdBu_9[i], strokecolor = :black, strokewidth = 1) for i in [1, 5, 9]]

Legend(fig₅[3, 3], elements, labels)

xlims!(ax₆, x0, xend)
ylims!(ax₆, y0, yend)

linkxaxes!(ax₁, ax₂, ax₃, ax₄, ax₅, ax₆)

for (ax, label) in zip([ax₁, ax₂, ax₃, ax₄, ax₅, ax₆], ["A", "B", "C", "D", "E", "F"])
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

save("plots/gp15Model/manuscript/fig5.png", fig₅, px_per_unit = 4)

# end plotting routine
