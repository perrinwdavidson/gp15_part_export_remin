# load packages
using CairoMakie
using XLSX
using DelimitedFiles
using NCDatasets
using ColorSchemes

# set plotting backend ::
CairoMakie.activate!(type = "png")
Makie.to_font("Arial")

# set basepaths ::
plot_data_basepath = pwd() * "/"  # "/Users/perrindavidson/research/whoi/current/gp15/" 

# load station data ::
station_data = XLSX.readxlsx(plot_data_basepath * "data/data_raw/stations/cruise_stations.xlsx")["Sheet1"][:]
station_data = station_data[2:end, :]

# load bathymetry data ::
Δ = 20
bath_nc = NCDataset(plot_data_basepath * "data/data_raw/gebco/gebco_2023_n65.0_s-25.0_w-180.0_e-60.0.nc")
x_bath = bath_nc["lon"][1:Δ:end]
y_bath = bath_nc["lat"][1:Δ:end]
bath = bath_nc["elevation"][1:Δ:end, 1:Δ:end]

# load npp data ::
Δ = 20
x_npp = readdlm(plot_data_basepath * "data/data_pro/readData/cbpm/npp_lon.csv", ',', Float64)
y_npp = readdlm(plot_data_basepath * "data/data_pro/readData/cbpm/npp_lat.csv", ',', Float64)
npp = readdlm(plot_data_basepath * "data/data_pro/calcNpp/npp_mean.csv", ',', Float64)

# set text ::
select_stations = [10, 18, 29, 39]
regions = [42, 7.5, -11]

function lat_names(x)
    if Int64(x) == 0
        return "$(Int64(x))°"  # "EQ"
    elseif Int64(x) > 0
        return "$(Int64(x))°N"
    elseif Int64(x) < 0
        return "$(Int64(-x))°S"
    end
end

# set characteristics ::
δ = 1.35
station_colors = Dict("port" => :black,  # ColorSchemes.RdBu_9[1], # Makie.wong_colors()[3], 
                      "test" => :black,  # ColorSchemes.RdBu_9[1],  # :black 
                      "intermediate" => ColorSchemes.RdBu_9[5],  # Makie.wong_colors()[7], 
                      "demi" => ColorSchemes.RdBu_9[3],  # :white, 
                      "full" => ColorSchemes.RdBu_9[9],  # Makie.wong_colors()[1], 
                      "super" => ColorSchemes.RdBu_9[1],  # Makie.wong_colors()[6], 
                      "shelf" => ColorSchemes.RdBu_9[7],  # Makie.wong_colors()[2], 
                      "slope" => ColorSchemes.RdBu_9[7])  # Makie.wong_colors()[2])
station_size = Dict("port" => 20/δ, 
                    "test" => 10/δ, 
                    "intermediate" => 15/δ, 
                    "demi" => 10/δ, 
                    "full" => 20/δ, 
                    "super" => 20/δ, 
                    "shelf" => 10/δ, 
                    "slope" => 10/δ)
station_shape = Dict("port" => :rect, 
                     "test" => :circle, 
                     "intermediate" => :dtriangle, 
                     "demi" => :circle, 
                     "full" => :circle, 
                     "super" => :circle, 
                     "shelf" => :circle, 
                     "slope" => :circle)

# make scatter plot ::
δ = 1.25
fig₁ = Figure(figure_padding = 25; size = (700/δ, 1000/δ))
ax₁ = Axis(fig₁[1, 1],
           xticklabelsize = 16,
           xtickformat = values -> ["$(Int64(-1 * value))°W" for value in values],
           yticklabelsize = 16,
           ytickformat = values -> [lat_names(value) for value in values])
cb = contourf!(x_npp[1, :], y_npp[:, 1], npp', colormap = (:GnBu, 0.75), levels = 0:10:100)
Colorbar(fig₁[1, 2], cb, label = rich("NPP (mmol C m", superscript("-2"), " d", superscript("-1"), ")"), ticks = 0:20:100, height = 350/δ, tellheight = false)
hlines!(ax₁, regions, linestyle = :dash, color = :black, linewidth = 1)
contourf!(x_bath, y_bath, bath, colormap = (:oslo, 0.75), levels = 0:0, extendhigh = :gray)  # also consider: deep, devon
for i in size(station_data)[1] : -1 : 1
    scatter!(ax₁, station_data[i, 4], station_data[i, 3], color = station_colors[station_data[i, 2]], markersize = station_size[station_data[i, 2]], marker = station_shape[station_data[i, 2]], strokecolor = :black, strokewidth = 1)
    stat = station_data[i, 1]
    if stat in select_stations
        text!(ax₁, station_data[i, 4] - 4.0, station_data[i, 3] - 1.01, text = "$(stat)", fontsize = 16)
    end
end

port = MarkerElement(color = :black, marker = :rect, markersize = 20/δ, strokecolor = :black, strokewidth = 1)
test = MarkerElement(color = :black, marker = :circle, markersize = 10/δ, strokecolor = :black, strokewidth = 1)
intermediate = MarkerElement(color = ColorSchemes.RdBu_9[5], marker = :dtriangle, markersize = 15/δ, strokecolor = :black, strokewidth = 1)
demi = MarkerElement(color = ColorSchemes.RdBu_9[3], marker = :circle, markersize = 10/δ, strokecolor = :black, strokewidth = 1)
full = MarkerElement(color = ColorSchemes.RdBu_9[9], marker = :circle, markersize = 20/δ, strokecolor = :black, strokewidth = 1)
super = MarkerElement(color = ColorSchemes.RdBu_9[1], marker = :circle, markersize = 20/δ, strokecolor = :black, strokewidth = 1)
shelf = MarkerElement(color = ColorSchemes.RdBu_9[7], marker = :circle, markersize = 10/δ, strokecolor = :black, strokewidth = 1)

axislegend(ax₁, [port, test, intermediate, demi, full, super, shelf], ["Port", "Test", "Intermediate", "Demi", "Full", "Super", "Shelf/Slope"], position = :rb)

xlims!(ax₁, -180, -120)
ylims!(ax₁, -22, 62)
save("plots/gp15Model/manuscript/fig1.png", fig₁, px_per_unit = 4)

# end plotting routine
