# load packages
using CairoMakie
using XLSX
using NCDatasets
using ColorSchemes

# set plotting backend ::
CairoMakie.activate!(type = "png")
Makie.to_font("Arial")

# set basepaths ::
plot_data_basepath = pwd() * "/"  # "/Users/perrindavidson/research/whoi/current/gp15/" 

# load station data ::
gp15_data = XLSX.readxlsx(plot_data_basepath * "data/data_pro/doQc/gp15/gp15_obs.xlsx")["Sheet1"][:]
gp15_vars = gp15_data[1, :]
gp15_data = gp15_data[2:end, :]

# load mld, eqp, and ppz data ::
gp15_stations = XLSX.readxlsx(plot_data_basepath * "data/data_pro/doQc/stations/gp15_stations.xlsx")["Sheet1"][:]
gp15_mld = XLSX.readxlsx(plot_data_basepath * "data/sim/calcMld/gp15_mld.xlsx")["Sheet1"][:]
gp15_eqp = XLSX.readxlsx(plot_data_basepath * "data/sim/calcEqp/gp15_eqp.xlsx")["Sheet1"][:]
gp15_station_vars = gp15_stations[1, :]
gp15_stations = gp15_stations[2:end, :]
gp15_mld_vars = gp15_mld[1, :]
gp15_mld = gp15_mld[2:end, :]
gp15_eqp_vars = gp15_eqp[1, :]
gp15_eqp = gp15_eqp[2:end, :]

# make scatter plot ::
δ = 1.2
num_stat = 33

fig₄ = Figure(figure_padding = 25; size = (1600/δ, 2750/δ))

for i in 1 : 1 : num_stat
    
    j = CartesianIndices((5, 7))
    k = j[i]

    stat_no = Int64(gp15_stations[i, 1])
              
    if k[1] != 1
        see_y_tick_label = false
    else
        see_y_tick_label = true
    end

    if k[2] != 1
        see_x_tick_label = false
    else
        see_x_tick_label = true
    end

    ax = Axis(fig₄[k[2], k[1]],
              # xlabel = x_label,
              # ylabel = y_label,
              # title = L"Station $%$(Int64(gp15_stations[i, 1]))$",
              titlesize = 20, 
              subtitlesize = 20, 
              xticklabelsize = 20,
              xaxisposition = :top,
              xticklabelsvisible = see_x_tick_label,
              yticklabelsize = 20,
              yticklabelsvisible = see_y_tick_label,
              yreversed = true,
              xgridvisible = false,
              ygridvisible = false)
    i_s = findall(x -> x == gp15_stations[i, 1], gp15_data[:, 2])
    band!(ax, 
          Point2.(Float64.(gp15_data[i_s, 11] .- gp15_data[i_s, 12]), Float64.(gp15_data[i_s, 4])), 
          Point2.(Float64.(gp15_data[i_s, 11] .+ gp15_data[i_s, 12]), Float64.(gp15_data[i_s, 4])), 
          color = (ColorSchemes.GnBu_9[7], 0.5))
    scatterlines!(ax, Float64.(gp15_data[i_s, 11]), Float64.(gp15_data[i_s, 4]), color = :black, marker = :rect, markersize = 12.5) 
    band!(ax, 
          Point2.(Float64.(gp15_data[i_s, 9] .- gp15_data[i_s, 10]), Float64.(gp15_data[i_s, 4])), 
          Point2.(Float64.(gp15_data[i_s, 9] .+ gp15_data[i_s, 10]), Float64.(gp15_data[i_s, 4])), 
          color = (ColorSchemes.GnBu_9[5], 0.5))
    scatterlines!(ax, Float64.(gp15_data[i_s, 9]), Float64.(gp15_data[i_s, 4]), color = :black, marker = :rect, markersize = 12.5, linestyle = :dash) 
    hlines!(ax, gp15_stations[i, 5], color = :black, linestyle = :dash)
    hlines!(ax, gp15_mld[i, 2], color = :black, linestyle = :dot)
    xlims!(ax, 1.0, 3.0)
    ylims!(ax, 400, 0)
    if stat_no in [27, 29, 31]
        local uci
        uci = min(gp15_stations[i, 5], gp15_eqp[i, 2])
        # hlines!(ax, gp15_eqp[i, 2], color = :black, linestyle = ".-")
        d1 = Float64.([1.0 gp15_mld[i, 2]; 1.0 uci])
        d2 = Float64.([3.0 gp15_mld[i, 2]; 3.0 uci])
        band!(ax, Point2.(d1[:, 1], d1[:, 2]), Point2.(d2[:, 1], d2[:, 2]), color = (ColorSchemes.GnBu_9[2], 0.5))
    end
    
    text!(ax, 0, 1,
          text = "Station $(stat_no)",
          font = :bold,
          align = (:left, :top),
          offset = (4, -242.5),
          space = :relative,
          fontsize = 24)

end

colgap!(fig₄.layout, 40)
rowgap!(fig₄.layout, 40)

th234 = [PolyElement(color = (ColorSchemes.GnBu_9[7], 0.5), strokecolor = :transparent, points = Point2f[(0, 0.25), (0, 0.75), (1, 0.75), (1, 0.25)]),
         LineElement(color = :black)]
u238 = [PolyElement(color = (ColorSchemes.GnBu_9[5], 0.5), strokecolor = :transparent, points = Point2f[(0, 0.25), (0, 0.75), (1, 0.75), (1, 0.25)]),
        LineElement(color = :black, linestyle = :dash)]
uci_plot = [PolyElement(color = (ColorSchemes.GnBu_9[2], 0.5), strokecolor = :transparent)]

Legend(fig₄[7, 4], [th234, u238, uci_plot], [rich(superscript("234"), "Th"), rich(superscript("238"), "U"), "UCI"], tellheight = false, tellwidth = false, labelsize = 30, patchsize = (35, 35))

Label(fig₄[0, 2:4], rich("Activity (dpm L", superscript("-1"), ")"), fontsize = 30)
Label(fig₄[3:5, 0], "Depth (m)", fontsize = 30, rotation = π/2)

save("plots/gp15Model/manuscript/fig4.png", fig₄, px_per_unit = 4)

# end plotting routine
