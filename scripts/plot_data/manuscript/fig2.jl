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
region_data = XLSX.readxlsx(plot_data_basepath * "data/sim/regressRegionalRatio/plotting/plot_output_data.xlsx")["Sheet1"][:]
region_fit = XLSX.readxlsx(plot_data_basepath * "data/sim/regressRegionalRatio/plotting/plot_output_fit.xlsx")["Sheet1"][:]
region_data = Float64.(region_data[2:end, :])
region_fit = Float64.(region_fit[2:end, :])

# load in data ::
gp15_data = XLSX.readxlsx(plot_data_basepath * "data/data_pro/doQc/gp15/gp15_obs.xlsx")["Sheet1"][:]
gp15_data_vars = gp15_data[1, :]
gp15_data = gp15_data[2:end, :]
gp15_data[findall(ismissing, gp15_data)] .= NaN

# load mld and ppz data ::
region_metrics = XLSX.readxlsx(plot_data_basepath * "data/sim/gp15Model/tables/table2.xlsx")["Sheet1"][:]
region_metrics = Float64.(region_metrics[2:end, 2:end])

# make scatter plot ::
δ = 1.75

fig₂ = Figure(figure_padding = 25; size = (2400/δ, 1000/δ))

ax₁ = Axis(fig₂[1:2, 1],
           # xlabel = L"LSF POC:$^{234}$Th ($\mu$mol C dpm$^{-1}$)",
           # ylabel = L"$$Depth (m)",
           ylabelsize = 20, 
           subtitlesize = 20, 
           xticklabelsize = 16,
           xaxisposition = :top,
           yticklabelsize = 16,
           yreversed = true,
           xgridvisible = false,
           ygridvisible = false)

i_d = findall(x -> x == 1, region_data[:, 1])
i_f = findall(x -> x == 1, region_fit[:, 1])
i_f_r = i_f[end:-1:1]
scatter!(ax₁, region_data[i_d, 3], region_data[i_d, 2], color = :black, marker = :rect, markersize = 12.5) 

i_gp15 = findall(x -> (x >= 0) & (x <= region_metrics[3, 1]), gp15_data[:, 4])
gp15_dat = gp15_data[i_gp15, :]
gp15_dat_region = gp15_dat[findall(x -> x in 1:1:10, gp15_dat[:, 2]), :]
scatter!(ax₁, Float64.(gp15_dat_region[:, 17]), Float64.(gp15_dat_region[:, 4]), color = :black, marker = :rect, markersize = 12.5) 

hlines!(ax₁, region_metrics[3, 1], color = :black, linestyle = :dot, linewidth = 1)
hlines!(ax₁, region_metrics[4, 1], color = :black, linestyle = :dash, linewidth = 1)
lines!(ax₁, region_fit[i_f, 3], region_fit[i_f, 2], color = ColorSchemes.GnBu_9[7], linewidth = 2) 
band!(ax₁, Point2.(region_fit[i_f, 4], region_fit[i_f, 2]), Point2.(region_fit[i_f, 5], region_fit[i_f, 2]), color = (ColorSchemes.GnBu_9[7], 0.5))
xlims!(ax₁, 0, 3.1)
ylims!(ax₁, 400, 0)

ax₂ = Axis(fig₂[1:2, 2],
           # xlabel = L"LSF POC:$^{234}$Th ($\mu$mol C dpm$^{-1}$)",
           # ylabel = L"$$Depth (m)",
           xlabelsize = 20, 
           xticklabelsize = 16,
           xaxisposition = :top,
           yticklabelsize = 16,
           yticklabelsvisible = false,
           yreversed = true,
           xgridvisible = false,
           ygridvisible = false)

i_d = findall(x -> x == 2, region_data[:, 1])
i_f = findall(x -> x == 2, region_fit[:, 1])
scatter!(ax₂, region_data[i_d, 3], region_data[i_d, 2], color = :black, marker = :rect, markersize = 12.5) 

i_gp15 = findall(x -> (x >= 0) & (x <= region_metrics[3, 3]), gp15_data[:, 4])
gp15_dat = gp15_data[i_gp15, :]
gp15_dat_region = gp15_dat[findall(x -> x in 11:1:23, gp15_dat[:, 2]), :]
scatter!(ax₂, Float64.(gp15_dat_region[:, 17]), Float64.(gp15_dat_region[:, 4]), color = :black, marker = :rect, markersize = 12.5) 

hlines!(ax₂, region_metrics[3, 3], color = :black, linestyle = :dot, linewidth = 1)
hlines!(ax₂, region_metrics[4, 3], color = :black, linestyle = :dash, linewidth = 1)
lines!(ax₂, region_fit[i_f, 3], region_fit[i_f, 2], color = ColorSchemes.GnBu_9[7], linewidth = 2) 
band!(ax₂, Point2.(region_fit[i_f, 4], region_fit[i_f, 2]), Point2.(region_fit[i_f, 5], region_fit[i_f, 2]), color = (ColorSchemes.GnBu_9[7], 0.5))
xlims!(ax₂, 0, 3.1)
ylims!(ax₂, 400, 0)

ax₃ = Axis(fig₂[1:2, 3],
           # xlabel = L"LSF POC:$^{234}$Th ($\mu$mol C dpm$^{-1}$)",
           # ylabel = L"$$Depth (m)",
           titlesize = 20, 
           xticklabelsize = 16,
           xaxisposition = :top,
           yticklabelsize = 16,
           yticklabelsvisible = false,
           yreversed = true,
           xgridvisible = false,
           ygridvisible = false)

i_d = findall(x -> x == 3, region_data[:, 1])
i_f = findall(x -> x == 3, region_fit[:, 1])
scatter!(ax₃, region_data[i_d, 3], region_data[i_d, 2], color = :black, marker = :rect, markersize = 12.5) 

i_gp15 = findall(x -> (x >= 0) & (x <= region_metrics[3, 5]), gp15_data[:, 4])
gp15_dat = gp15_data[i_gp15, :]
gp15_dat_region = gp15_dat[findall(x -> x in 24:1:35, gp15_dat[:, 2]), :]
scatter!(ax₃, Float64.(gp15_dat_region[:, 17]), Float64.(gp15_dat_region[:, 4]), color = :black, marker = :rect, markersize = 12.5) 

hlines!(ax₃, region_metrics[3, 5], color = :black, linestyle = :dot, linewidth = 1)
hlines!(ax₃, region_metrics[4, 5], color = :black, linestyle = :dash, linewidth = 1)
lines!(ax₃, region_fit[i_f, 3], region_fit[i_f, 2], color = ColorSchemes.GnBu_9[7], linewidth = 2) 
band!(ax₃, Point2.(region_fit[i_f, 4], region_fit[i_f, 2]), Point2.(region_fit[i_f, 5], region_fit[i_f, 2]), color = (ColorSchemes.GnBu_9[7], 0.5))
xlims!(ax₃, 0, 3.1)
ylims!(ax₃, 400, 0)

ax₄ = Axis(fig₂[1:2, 4],
           # xlabel = L"LSF POC:$^{234}$Th ($\mu$mol C dpm$^{-1}$)",
           # ylabel = L"$$Depth (m)",
           titlesize = 20, 
           xticklabelsize = 16,
           xaxisposition = :top,
           yticklabelsize = 16,
           yticklabelsvisible = false,
           yreversed = true,
           xgridvisible = false,
           ygridvisible = false)

i_d = findall(x -> x == 4, region_data[:, 1])
i_f = findall(x -> x == 4, region_fit[:, 1])
scatter!(ax₄, region_data[i_d, 3], region_data[i_d, 2], color = :black, marker = :rect, markersize = 12.5) 

i_gp15 = findall(x -> (x >= 0) & (x <= region_metrics[3, 7]), gp15_data[:, 4])
gp15_dat = gp15_data[i_gp15, :]
gp15_dat_region = gp15_dat[findall(x -> x in 36:1:39, gp15_dat[:, 2]), :]
scatter!(ax₄, Float64.(gp15_dat_region[:, 17]), Float64.(gp15_dat_region[:, 4]), color = :black, marker = :rect, markersize = 12.5) 

hlines!(ax₄, region_metrics[3, 7], color = :black, linestyle = :dot, linewidth = 1)
hlines!(ax₄, region_metrics[4, 7], color = :black, linestyle = :dash, linewidth = 1)
lines!(ax₄, region_fit[i_f, 3], region_fit[i_f, 2], color = ColorSchemes.GnBu_9[7], linewidth = 2) 
band!(ax₄, Point2.(region_fit[i_f, 4], region_fit[i_f, 2]), Point2.(region_fit[i_f, 5], region_fit[i_f, 2]), color = (ColorSchemes.GnBu_9[7], 0.5))
xlims!(ax₄, 0, 3.1)
ylims!(ax₄, 400, 0)

Label(fig₂[0, 2:3], rich("LSF POC:", superscript("234"), "Th (μmol C dpm", superscript("-1"), ")"), fontsize = 24)
Label(fig₂[1:2, 0], "Depth (m)", fontsize = 24, rotation = π/2)

linkaxes!(ax₁, ax₂, ax₃, ax₄)

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

save("plots/gp15Model/manuscript/fig2.png", fig₂, px_per_unit = 4)

# end plotting routine
