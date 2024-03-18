# load packages
using CairoMakie
using XLSX
using NCDatasets
using ColorSchemes

# set plotting backend ::
CairoMakie.activate!(type = "png")
Makie.to_font("Arial")

# set plotting basepath ::
plot_data_basepath = pwd() * "/"  # "/Users/perrindavidson/research/whoi/current/gp15/" 

# load ctd data ::
gp15_data = XLSX.readxlsx(plot_data_basepath * "data/data_pro/doQc/gp15/gp15_ctd.xlsx")["Sheet1"][:]
gp15_data_vars = gp15_data[1, :]
gp15_data = gp15_data[2:end, :]

# load ppz ::
gp15_ppz = XLSX.readxlsx(plot_data_basepath * "data/data_pro/doQc/ppz/gp15_ppz.xlsx")["Sheet1"][:]
gp15_ppz_vars = gp15_ppz[1, :]
gp15_ppz = gp15_ppz[2:end, :]

# load station data ::
gp15_stations = XLSX.readxlsx(plot_data_basepath * "data/data_pro/doQc/stations/gp15_stations.xlsx")["Sheet1"][:]
gp15_station_vars = gp15_stations[1, :]
gp15_stations = gp15_stations[2:end, :]

# make scatter plot ::
δ = 1.2
num_stat = 33

fig = Figure(figure_padding = 25; size = (1600/δ, 2750/δ))

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

    ax = Axis(fig[k[2], k[1]],
              # xlabel = x_label,
              # ylabel = y_label,
              # title = L"Station $%$(Int64(gp15_stations[i, 1]))$",
              titlesize = 20, 
              subtitlesize = 20, 
              xticklabelsize = 20,
              xaxisposition = :top,
              xticklabelsvisible = see_x_tick_label,
              xticks = 0:0.05:0.15,
              yticklabelsize = 20,
              yticklabelsvisible = see_y_tick_label,
              yreversed = true,
              xgridvisible = false,
              ygridvisible = false)
    i_s = findall(x -> x == gp15_stations[i, 1], gp15_data[:, 1])
    i_c = findall(x -> x == gp15_stations[i, 2], gp15_data[:, 2])
    idx = intersect(i_c, i_s)

    scatterlines!(ax, Float64.(gp15_data[idx, end-2]), Float64.(gp15_data[idx, 9]), color = :black, marker = :rect, markersize = 2.5) 
    hlines!(ax, gp15_stations[i, 5], color = ColorSchemes.RdBu[1], linestyle = :dash)
    xlims!(ax, 0.0, 0.15)
    ylims!(ax, 400, 0)
    
    text!(ax, 1, 0,
          text = "Station $(stat_no)",
          font = :bold,
          align = (:right, :bottom),
          offset = (-4, 0),
          space = :relative,
          fontsize = 24)

end

colgap!(fig.layout, 40)
rowgap!(fig.layout, 40)

Label(fig[0, 2:4], "Fluorescence (volts)", fontsize = 30)
Label(fig[3:5, 0], "Depth (m)", fontsize = 30, rotation = π/2)

save("plots/gp15Model/supporting_information/fig_ppz.png", fig, px_per_unit = 4)

# end plotting routine
