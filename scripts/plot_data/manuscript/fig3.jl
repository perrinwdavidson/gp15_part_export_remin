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
stat_data = XLSX.readxlsx(plot_data_basepath * "data/data_pro/doQc/stations/gp15_stations.xlsx")["Sheet1"][:]
stat_vars = stat_data[1, :]
stat_data = Float64.(stat_data[2:end, 1:end-1])
ctd_data = XLSX.readxlsx(plot_data_basepath * "data/data_pro/doQc/gp15/gp15_ctd.xlsx")["Sheet1"][:]
ctd_vars = ctd_data[1, :]
ctd_data = Float64.(ctd_data[2:end, :])
bottle_data = XLSX.readxlsx(plot_data_basepath * "data/data_pro/doQc/gp15/gp15_obs.xlsx")["Sheet1"][:]
bottle_data_vars = bottle_data[1, :]
bottle_data = bottle_data[2:end, :]

# loop through all stations and make array ::
function make_data(zq)
    num_stats = size(stat_data)[1]
    ctd_pdens = zeros(length(zq), num_stats)
    ctd_sal = zeros(length(zq), num_stats)
    ctd_fluoro = zeros(length(zq), num_stats)
    lsf_234th = zeros(length(zq), 0)
    ssf_234th = zeros(length(zq), 0)
    th234 = zeros(length(zq), num_stats)
    u238 = zeros(length(zq), num_stats)
    lsf_lat = []
    ssf_lat = []
    for i in 1 : 1 : num_stats
        
        # ctd ::
        j₁ = findall(x -> x == stat_data[i, 1], ctd_data[:, 1])
        j₂ = findall(x -> x == stat_data[i, 2], ctd_data[:, 2])
        j = intersect(j₁, j₂)
        dat = ctd_data[j, :]
        f_fluoro = linear_interpolation(dat[:, 9], dat[:, 20]; extrapolation_bc = Linear())
        ctd_fluoro[:, i] .= f_fluoro(zq)
        f_pdens = linear_interpolation(dat[:, 9], dat[:, 22]; extrapolation_bc = Linear())
        ctd_pdens[:, i] .= f_pdens(zq)
        f_sal = linear_interpolation(dat[:, 9], dat[:, 12]; extrapolation_bc = Linear())
        ctd_sal[:, i] .= f_sal(zq)

        # bottle ::
        j = findall(x -> x == stat_data[i, 1], bottle_data[:, 2])
        dat = bottle_data[j, :]

        lsf = dat[:, 15]
        ssf = dat[:, 21]
        ilsf = findall(!ismissing, lsf)
        issf = findall(!ismissing, ssf)
        lsf = lsf[ilsf]
        zlsf = dat[ilsf, 4]
        ssf = ssf[issf]
        zssf = dat[issf, 4]
        ilsf = findall(x -> x < 0.2, lsf)
        issf = findall(x -> x < 1.5, ssf)

        if ~isempty(lsf)
            f_lsf = linear_interpolation(zlsf[ilsf], lsf[ilsf]; extrapolation_bc = Linear())  
            lsf_234th = hcat(lsf_234th, f_lsf(zq))
            lsf_lat = append!(lsf_lat, stat_data[i, 4])
        end
        if ~isempty(ssf)
            f_ssf = linear_interpolation(zssf[issf], ssf[issf]; extrapolation_bc = Linear())  
            ssf_234th = hcat(ssf_234th, f_ssf(zq))
            ssf_lat = append!(ssf_lat, stat_data[i, 4])
        end
        f_th = linear_interpolation(dat[:, 4], dat[:, 11]; extrapolation_bc = Linear()) 
        th234[:, i] .= f_th(zq)
        f_u = linear_interpolation(dat[:, 4], dat[:, 9]; extrapolation_bc = Linear())  
        u238[:, i] .= f_u(zq)

    end
    return ctd_pdens, ctd_sal, ctd_fluoro, lsf_234th, ssf_234th, th234, u238, lsf_lat, ssf_lat
end

zq = 0 : 1 : 400
ctd_pdens, ctd_sal, ctd_fluoro, lsf_234th, ssf_234th, th234, u238, lsf_lat, ssf_lat = make_data(zq)

# make scatter plot :: add in sample locations and change opacy
function lat_names(x)
    if Int64(x) == 0
        return "$(Int64(x))°"  # EQ
    elseif Int64(x) > 0
        return "$(Int64(x))°N"
    elseif Int64(x) < 0
        return "$(Int64(-x))°S"
    end
end

regions = [42, 7.5, -11]

δ = 1.5

fig₃ = Figure(figure_padding = 25; size = (2000/δ, 1200/δ))

ax₁ = Axis(fig₃[1, 1],
           # xlabel = L"$$Latitude", 
           # ylabel = L"$$Depth (m)",
           titlesize = 20, 
           subtitlesize = 20, 
           xticklabelsize = 16,
           xaxisposition = :top,
           xtickformat = values -> [lat_names(value) for value in values],
           yticklabelsize = 16,
           yreversed = true)
cb = contourf!(stat_data[:, 4], zq, ctd_fluoro', colormap = (:GnBu_9, 0.75), levels = 0:0.01:0.15)
Colorbar(fig₃[1, 2], cb, label = "Fluorescence (relative volts)", ticks = 0:0.05:0.15)
scatter!(ax₁, stat_data[:, 4], zeros(size(stat_data[:, 4])), marker = :dtriangle, color = :black, markersize = 15)
contour!(ax₁, stat_data[:, 4], zq, ctd_fluoro'; color = :black, levels = [0.03])
vlines!(ax₁, regions, linestyle = :dash, color = :black, linewidth = 1)
xlims!(ax₁, stat_data[end, 4], stat_data[1, 4])
ylims!(ax₁, 400, 0)

ax₂ = Axis(fig₃[2, 1],
           # xlabel = L"$$Latitude", 
           # ylabel = L"$$Depth (m)",
           titlesize = 20, 
           subtitlesize = 20, 
           xticklabelsize = 16,
           xticklabelsvisible = false,
           xaxisposition = :top,
           xtickformat = values -> [lat_names(value) for value in values],
           yticklabelsize = 16,
           yreversed = true)
cb = contourf!(stat_data[:, 4], zq, ctd_pdens', colormap = (:GnBu_9, 0.75), levels = 20:0.5:28)
Colorbar(fig₃[2, 2], cb, label = rich("Potential density (kg m", superscript("-3"), ")"), ticks = 20:2:28)
scatter!(ax₂, stat_data[:, 4], zeros(size(stat_data[:, 4])), marker = :dtriangle, color = :black, markersize = 15)
vlines!(ax₂, regions, linestyle = :dash, color = :black, linewidth = 1)
xlims!(ax₂, stat_data[end, 4], stat_data[1, 4])
ylims!(ax₂, 400, 0)

ax₃ = Axis(fig₃[3, 1],
           # xlabel = L"$$Latitude", 
           # ylabel = L"$$Depth (m)",
           titlesize = 20, 
           subtitlesize = 20, 
           xticklabelsize = 16,
           xticklabelsvisible = false,
           xaxisposition = :top,
           xtickformat = values -> [lat_names(value) for value in values],
           yticklabelsize = 16,
           yreversed = true)
cb = contourf!(stat_data[:, 4], zq, ctd_sal', colormap = (:GnBu_9, 0.75), levels = 32:0.5:37)
Colorbar(fig₃[3, 2], cb, label = "Salinity", ticks = 32:1:37)
scatter!(ax₃, stat_data[:, 4], zeros(size(stat_data[:, 4])), marker = :dtriangle, color = :black, markersize = 15)
vlines!(ax₃, regions, linestyle = :dash, color = :black, linewidth = 1)
xlims!(ax₃, stat_data[end, 4], stat_data[1, 4])
ylims!(ax₃, 400, 0)

ax₄ = Axis(fig₃[1, 3],
           # xlabel = L"$$Latitude", 
           # ylabel = L"$$Depth (m)",
           titlesize = 20, 
           subtitlesize = 20, 
           xticklabelsize = 16,
           xaxisposition = :top,
           xtickformat = values -> [lat_names(value) for value in values],
           yticklabelsize = 16,
           yreversed = true)

lsf_scatter = hcat(bottle_data[:, 5], bottle_data[:, 4], bottle_data[:, 15])
lsf_scatter = Float64.(lsf_scatter[findall(!ismissing, lsf_scatter[:, 3]), :])

cb = contourf!(Float64.(lsf_lat), Float64.(collect(zq)), lsf_234th', colormap = (:GnBu_9, 0.75), levels = 0:0.025:0.2)  # range(0, 0.2, length = 10)
Colorbar(fig₃[1, 4], cb, label = rich("LSF ", superscript("234"), "Th (dpm L", superscript("-1"), ")"), ticks = 0:0.05:0.2)
contour!(ax₄, stat_data[:, 4], zq, ctd_fluoro'; color = :black, levels = [0.03])

ilsf = findall(x -> x < 0.2, lsf_scatter[:, 3])
scatter!(ax₄, lsf_scatter[ilsf, 1], lsf_scatter[ilsf, 2], color = lsf_scatter[ilsf, 3], marker = :rect, markersize = 10, strokecolor = :black, strokewidth = 1, colormap = :GnBu_9, levels = 0:0.025:0.2)

vlines!(ax₄, regions, linestyle = :dash, color = :black, linewidth = 1)

# scatter!(ax₁, stat_data[:, 4], zeros(size(stat_data[:, 4])), marker = :dtriangle, color = :black)
xlims!(ax₄, stat_data[end, 4], stat_data[1, 4])
ylims!(ax₄, 400, 0)

ax₅ = Axis(fig₃[2, 3],
           # xlabel = L"$$Latitude", 
           # ylabel = L"$$Depth (m)",
           titlesize = 20, 
           subtitlesize = 20, 
           xticklabelsize = 16,
           xaxisposition = :top,
           xticklabelsvisible = false,
           xtickformat = values -> [lat_names(value) for value in values],
           yticklabelsize = 16,
           yreversed = true)

ssf_scatter = hcat(bottle_data[:, 5], bottle_data[:, 4], bottle_data[:, 21])
ssf_scatter = Float64.(lsf_scatter[findall(!ismissing, lsf_scatter[:, 3]), :])

cb = contourf!(Float64.(ssf_lat), Float64.(collect(zq)), ssf_234th', colormap = (:GnBu_9, 0.75), levels = 0:0.125:1.5)  # range(0, 1.5, length = 10)
Colorbar(fig₃[2, 4], cb, label = rich("SSF ", superscript("234"), "Th (dpm L", superscript("-1"), ")"), ticks = 0:0.5:1.5)
contour!(ax₅, stat_data[:, 4], zq, ctd_fluoro'; color = :black, levels = [0.03])
issf = findall(x -> x < 1.5, ssf_scatter[:, 3])
scatter!(ax₅, ssf_scatter[issf, 1], ssf_scatter[issf, 2], color = ssf_scatter[issf, 3], marker = :rect, markersize = 10, strokecolor = :black, strokewidth = 1, colormap = :GnBu_9, levels = 0:0.125:1.5)
vlines!(ax₅, regions, linestyle = :dash, color = :black, linewidth = 1)
# scatter!(ax₁, stat_data[:, 4], zeros(size(stat_data[:, 4])), marker = :dtriangle, color = :black)
xlims!(ax₅, stat_data[end, 4], stat_data[1, 4])
ylims!(ax₅, 400, 0)

ax₆ = Axis(fig₃[3, 3],
           # xlabel = L"$$Latitude", 
           # ylabel = L"$$Depth (m)",
           titlesize = 20, 
           subtitlesize = 20, 
           xticklabelsize = 16,
           xaxisposition = :top,
           xtickformat = values -> [lat_names(value) for value in values],
           xticklabelsvisible = false,
           yticklabelsize = 16,
           yreversed = true)
cb = contourf!(stat_data[:, 4], zq, u238' .- th234', colormap = (:GnBu_9, 0.75), levels = -0.5:0.125:1.5)  # = range(-0.5, 1.2, length = 10)
Colorbar(fig₃[3, 4], cb, label = rich(superscript("238"), "U – ", superscript("234"), "Th (dpm L", superscript("-1"), ")"), ticks = -0.5:0.5:1.5)
contour!(ax₆, stat_data[:, 4], zq, ctd_fluoro'; color = :black, levels = [0.03])
scatter!(ax₆, Float64.(bottle_data[:, 5]), Float64.(bottle_data[:, 4]), color =  Float64.(bottle_data[:, 9]) .-  Float64.(bottle_data[:, 11]), marker = :rect, markersize = 10, strokecolor = :black, strokewidth = 1, colormap = :GnBu_9, levels = -0.5:0.5:1.5)
vlines!(ax₆, regions, linestyle = :dash, color = :black, linewidth = 1)
# scatter!(ax₁, stat_data[:, 4], zeros(size(stat_data[:, 4])), marker = :dtriangle, color = :black)
xlims!(ax₆, stat_data[end, 4], stat_data[1, 4])
ylims!(ax₆, 400, 0)

linkaxes!(ax₁, ax₂, ax₃, ax₄, ax₅, ax₆)

# Label(fig₃[0, 2:3], L"LSF POC:$^{234}$Th ($\mu$mol C dpm$^{-1}$)", fontsize = 20)
Label(fig₃[1:3, 0], "Depth (m)", fontsize = 20, rotation = π/2)

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

save("plots/gp15Model/manuscript/fig3.png", fig₃, px_per_unit = 4)

# end plotting routine
