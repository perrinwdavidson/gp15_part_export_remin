# GP15 Particulate Export & Remineralization — Code Manual

**Project**: Kenyon and Davidson et al. (2026), GP15 Pacific GEOTRACES transect  
**Languages**: MATLAB (pipeline), Julia (figures), Python via `terraformer` CLI (kriging)  
**Last updated**: 2026-06-05 (Session 7: testFitPiecewiseRatio output-order fix; §24.4 call-signature corrected)

---

## Table of Contents

1. [Scientific Overview](#1-scientific-overview)
2. [Repository Structure](#2-repository-structure)
3. [Conventions and Constants](#3-conventions-and-constants)
4. [Pipeline Execution Order](#4-pipeline-execution-order)
5. [Stage 1 — Read Data](#5-stage-1--read-data)
6. [Stage 2 — Particle Flux Proxy Depth (PPZ)](#6-stage-2--particle-flux-proxy-depth-ppz)
7. [Stage 3 — Quality Control](#7-stage-3--quality-control)
8. [Stage 4 — Net Primary Production](#8-stage-4--net-primary-production)
9. [Stage 5 — Mixed Layer Depth](#9-stage-5--mixed-layer-depth)
10. [Stage 6 — Equilibrium Point Depth (EQP)](#10-stage-6--equilibrium-point-depth-eqp)
11. [Stage 7 — Upwelling Velocity](#11-stage-7--upwelling-velocity)
12. [Stage 8 — Upwelling Interpolation](#12-stage-8--upwelling-interpolation)
13. [Stage 9 — Vertical Th-234 Gradient](#13-stage-9--vertical-th-234-gradient)
14. [Stage 10 — Regional Delineation](#14-stage-10--regional-delineation)
15. [Stage 11 — POC:Th and PN:Th Ratio Regression](#15-stage-11--pocth-and-pnth-ratio-regression)
16. [Stage 12 — Data Collation](#16-stage-12--data-collation)
17. [Stage 13 — Th-234 Flux Model](#17-stage-13--th-234-flux-model)
18. [Stage 14 — Error Propagation](#18-stage-14--error-propagation)
19. [Stage 15 — Metrics and Statistics](#19-stage-15--metrics-and-statistics)
20. [Toolboxes](#20-toolboxes)
21. [Numerical Invariants and Fragile Areas](#21-numerical-invariants-and-fragile-areas)
22. [Function Inventory](#22-function-inventory)
23. [Interpolation Dependency Map](#23-interpolation-dependency-map)
24. [Canonical Method Reference](#24-canonical-method-reference)
25. [Tests](#25-tests)

---

## 1. Scientific Overview

This repository computes particulate organic carbon (POC) and particulate nitrogen (PN) export fluxes and remineralization metrics along the GP15 Pacific GEOTRACES transect (Alaska to Tahiti, 2018). The approach uses naturally-occurring thorium-234 (²³⁴Th) as a particle-reactive proxy for gravitational settling flux, calibrated against in-situ POC and PN measurements.

### Physical and chemical basis

²³⁴Th is produced continuously in seawater by decay of dissolved ²³⁸U, and is removed by scavenging onto sinking particles. In steady state, the excess production over in-situ decay equals the flux of ²³⁴Th leaving the water column on particles:

**1D (steady-state) model:**
```
F_Th(z) = ∫₀ᶻ λ · [²³⁸U(z') − ²³⁴Th(z')] dz'   [dpm m⁻² d⁻¹]
```

where λ = ln(2) / t½ is the ²³⁴Th decay constant (t½ = 24.101 d).

**2D model (equatorial stations):** At stations near the equator, horizontal divergence drives significant upwelling, which advects dissolved ²³⁴Th upward and violates the steady-state assumption. The 2D correction adds an advection term:

```
F_Th,2d(z) = ∫₀ᶻ [λ · (²³⁸U − ²³⁴Th) + w · (∂²³⁴Th/∂z)] dz'
```

where w is the vertical velocity (m d⁻¹) and ∂²³⁴Th/∂z is the vertical gradient of ²³⁴Th.

**POC/PN flux** is obtained by multiplying the Th-234 flux by the measured particulate ratio:

```
F_POC(z) = F_Th(z) × R_POC(z, region)      [mmol C m⁻² d⁻¹]
F_PN(z)  = F_Th(z) × R_PN(z, region)       [mmol N m⁻² d⁻¹]
```

where R_POC = POC:²³⁴Th (mmol dpm⁻¹) is determined by depth-binned weighted linear regression within each biogeographic region.

### Key output metrics

- **F_POC(100 m)** — POC export flux at 100 m
- **F_POC(PPZ)** — POC export flux at the particle flux proxy depth
- **T₁₀₀ = F_POC(PPZ+100) / F_POC(PPZ)** — transfer efficiency: fraction of PPZ flux surviving 100 m below
- **R₁₀₀ = F_Th(PPZ+100) − F_Th(PPZ)** — Th-234 remineralization between PPZ and PPZ+100 m
- **Ez Ratio = F_POC(PPZ) / NPP** — export efficiency: fraction of surface production exported

---

## 2. Repository Structure

```
gp15_part_export_remin/
├── scripts/                     ← pipeline entry points and subroutines
│   ├── gp15Model.m              ← MASTER SCRIPT
│   ├── model/                   ← Th-234 flux model (core science)
│   ├── profile_interpolation/   ← spatiotemporal kriging of upwelling
│   ├── quality_control/         ← QC, duplicate averaging, flier removal
│   ├── read_data/               ← raw data ingestion
│   ├── collate_data/            ← assembles final model input table
│   ├── calculate_ppz/           ← particle flux proxy depth
│   ├── calculate_npp/           ← satellite NPP kriging
│   ├── calculate_mld/           ← mixed layer depth (3 methods)
│   ├── calculate_eqp/           ← equilibrium point depth (Th234/U238 = 1)
│   ├── calculate_upwell/        ← vertical velocity from model divergence
│   ├── calculate_gradient/      ← vertical Th-234 gradient
│   ├── regional_delineation/    ← biogeographic province assignment
│   ├── regress_ratio/           ← POC:Th and PN:Th depth regression
│   └── plot_data/               ← Julia plotting scripts (manuscript + SI)
├── src/
│   ├── toolboxes/
│   │   ├── terraformer/         ← ordinary_kriging (Python CLI wrappers)
│   │   ├── interp1dError/       ← piecewise-linear interpolation with uncertainty
│   │   └── regression/          ← ols.m (OLS/WLS regression)
│   └── utils/model/             ← setupModel, printStart, printEnd
├── data/
│   ├── data_raw/                ← immutable raw inputs
│   ├── data_pro/                ← processed intermediates
│   └── sim/                     ← simulation outputs
├── docs/                        ← methodology and code documentation
└── plots/                       ← output figures
```

Each stage script can be run standalone; `setupModel.m` sets all paths via `mfilename('fullpath')` and clears the workspace.

---

## 3. Conventions and Constants

### Unit systems

Three unit systems coexist. Conversions occur **only in `collateStations.m`**:

| Quantity | Storage units | Model units | Conversion factor |
|----------|---------------|-------------|-------------------|
| ²³⁴Th, ²³⁸U | dpm L⁻¹ | dpm m⁻³ | × 1000 (L2M) |
| POC, PN | µmol | mmol | ÷ 1000 (UMOL2MMOL) |
| Upwelling w | m s⁻¹ | m d⁻¹ | × 86400 (SEC2DAY) |
| Gradient ∂Th/∂z | dpm L⁻¹ m⁻¹ | dpm m⁻³ m⁻¹ | × 1000 (L2M) |
| POC:Th ratio | µmol dpm⁻¹ | mmol dpm⁻¹ | ÷ 1000 (UMOL2MMOL) |

Fluxes are reported as dpm m⁻² d⁻¹ (Th-234) and mmol m⁻² d⁻¹ (POC, PN).

### Depth convention

All depths are positive downward. The integration ceiling is `BTMDEPTH = 400 m`, hardcoded in both `configureInterp.m` and `collateStations.m` — these must remain consistent.

### Longitude convention

Station longitudes are in [−180, 180]. Model grids (ECCO, MERCATOR) use [0, 360]. The `+360` offset is applied in `interpUpwell.m` when constructing query coordinates, and removed when saving output (`query_lon − 360`).

### String style

All scripts use single-quoted `char` arrays (`'like_this'`), not double-quoted `string` objects. Per-product loops in `calcError.m`, `doModelCalcs.m`, `calcParticulateFlux.m`, and `collateStations.m` use MATLAB dynamic struct/table field access — `s.(['field_' dataProduct])` — which requires char arrays and does not use `eval()`.

### Critical constants

| Constant | Value | Location | Role |
|----------|-------|----------|------|
| `LAMBDA` | ln(2) / 24.101 d⁻¹ | `setModelCoefficients.m` | ²³⁴Th decay constant |
| `BTMDEPTH` | 400 m | `configureInterp.m`, `collateStations.m` | Integration ceiling |
| `FLIERDATA` | 30 µmol dpm⁻¹ | `qcGp15Poc.m` | POC:Th ratio flier cutoff |
| `DepthBinEdges` | [0, 100, 200, 400, ∞] m | `interpUpwell.m` | Bins for temporal preprocessing |
| `cbpmFrac` | 0.30 | `calcStats.m` | CBPM algorithmic uncertainty (30%) |
| Equatorial band | \|lat\| ≤ 5° | `calcStats.m` | Stations receiving 2D upwelling correction |

---

## 4. Pipeline Execution Order

**Master script**: `scripts/gp15Model.m`

```
setupModel
  │
  ▼
readData          ← raw GP15 observations, ECCO model, CBPM NPP
calcPpz           ← particle flux proxy depth per station
doQc              ← QC, interpolate ratios, remove fliers
calcNpp           ← kriged NPP at station coordinates
calcMld           ← mixed layer depth (3 algorithms)
calcEqp           ← Th-234/U-238 equilibrium point depth
calcUpwell        ← vertical velocity from divergence integral
interpData        ← spatiotemporal kriging of w to station profiles
calcVertGrad      ← vertical Th-234 gradient with uncertainty
delineateRegions  ← biogeographic province assignment
regressRegionalRatio ← depth-binned POC:Th and PN:Th regression
collateData       ← assemble full model input table
  │
  ▼
runModel          ← fluxes, error propagation, metrics, output tables
```

---

## 5. Stage 1 — Read Data

**Script**: `scripts/read_data/`

### GP15 observations (`readGp15.m`)

Ingests raw Excel spreadsheets from `data_raw/gp15/`:

| File | Contents |
|------|----------|
| `gp15_stat_cast.xlsx` | Station metadata (number, lat, lon, date, bottom depth) |
| `gp15_master.xlsx` (sheet: tTh-234_Summary) | Total ²³⁴Th and ²³⁸U activities (dpm L⁻¹) with uncertainties |
| `gp15_master.xlsx` (sheet: pTh234) | Particulate ²³⁴Th by size fraction (SSF, LSF) |
| `gp15_particles.xlsx` | POC and PN concentrations by size fraction (µmol L⁻¹) |
| `gp15_ctd.xlsx` | Potential density, temperature, salinity profiles |
| `gp15_bottles.xlsx` | Nutrient chemistry (NO₃, PO₄) |
| `gp15_hplc.xlsx` | Pigment concentrations (Chl a, accessory pigments) |

**Output**: `gp15_obs` — a table with one row per (station, depth) observation containing all radionuclide and particulate measurements.

### ECCO model velocities (`readEcco.m`)

Reads four quarterly ECCO netCDF files (`ecco_wo_quart1–4.cdf`), concatenates along the time axis, converts ECCO fill values (−10¹⁰) to NaN, and serializes time to MATLAB `datenum` format.

**Output**: `ecco_wo` struct with 4D velocity arrays (longitude × latitude × depth × time).

### CBPM satellite NPP (`readNpp.m`)

Reads daily CBPM (Carbon-based Productivity Model) NPP HDF5 files from `data_raw/cbpm/`. Each file contains a 1080 × 2160 global array at 0.167° × 0.167° resolution in units of mg C m⁻² d⁻¹.

**Operations**:
- Date string parsed from filename to build time axis
- Fill values removed
- Units converted: mg C → mmol C by dividing by the molar mass of carbon (12 g mol⁻¹)
- Coordinate matrices constructed for the native grid (lat centers 90° to −90°, lon centers −180° to 180°)

**Output**: `NPP` (lat × lon × time), `lat`, `lon`, `nppTime` arrays; also saved as netCDF4.

---

## 6. Stage 2 — Particle Flux Proxy Depth (PPZ)

**Script**: `scripts/calculate_ppz/`

The PPZ is the depth below which the ²³⁴Th deficit (excess ²³⁸U production over ²³⁴Th) is largest, indicating active particle export. It serves as the primary integration depth for computing export fluxes.

**Algorithm**: The PPZ is identified from the observed ²³⁴Th profile per station using a smoothing spline or threshold criterion applied to the ²³⁴Th deficit profile. The exact criterion is documented in `calcPpzDepths.m`.

**Output**: `gp15_ppz` table with `stationNo` and `ppzDepth` (m) per station.

---

## 7. Stage 3 — Quality Control

**Script**: `scripts/quality_control/subroutines/subsubroutines/qcGp15Poc.m`

### Interpolation of particulate data to observation depths

Raw POC and particulate ²³⁴Th measurements are collected at different depths than the total ²³⁴Th samples. Piecewise-linear interpolation with uncertainty propagation (`interp1dError.m`) maps all particulate quantities to the common total ²³⁴Th observation depths.

### Ratio calculation and uncertainty propagation

At each matched depth:

```
R = POC / Th_part    [µmol dpm⁻¹]
```

Uncertainty via the delta method:

```
σ_R = |R| · √[(σ_POC / POC)² + (σ_Th / Th_part)²]
```

Performed separately for large size fraction (LSF, >51 µm) and small size fraction (SSF, 1–51 µm).

### Flier removal

Ratios exceeding `FLIERDATA = 30 µmol dpm⁻¹` are flagged as outliers. The flier index is captured **before** NaN-ification so that uncertainty columns are also nulled consistently. The `flierIdxLarge` / `flierIdxSmall` capture pattern must be preserved if this section is ever refactored.

**Output**: `gp15_obs` augmented with `pocTh234RatioLarge`, `uncertPocTh234RatioLarge`, `pocTh234RatioSmall`, `uncertPocTh234RatioSmall` columns (and analogues for PN).

**Diagnostic plots** (`plotQcPoc.m`, `plotQcPn.m`): standalone scripts in `quality_control/subroutines/subsubroutines/`. Each generates one 8-panel PDF per station showing depth profiles of U-238, total Th-234, SSF and LSF particulate Th-234, POC/PN, and POC:Th / PN:Th ratios. Measurement uncertainties are shown as horizontal error bars (symmetric, capSize = 0) centred on each measurement point.

---

## 8. Stage 4 — Net Primary Production

**Scripts**: `scripts/calculate_npp/subroutines/`

### Motivation

NPP is needed to compute the Ez Ratio (export efficiency). Satellite CBPM NPP is available on a global grid, but must be interpolated to exact station coordinates and sampling dates.

### Algorithm

**Step 1 — Temporal preprocessing** (`fit_temporal_preprocessor`, via `terraformer`):

A harmonic regression (annual + semi-annual sinusoid + linear trend) is fit independently per depth bin [0, 100, 200, 400, ∞] m to the full training NPP field. This removes the deterministic seasonal cycle, leaving anomalies `anomV` which are more spatially smooth and easier to krige.

**Step 2 — Spatiotemporal ordinary kriging** (`ordinary_kriging`, via `terraformer`):

Kriging is performed on the anomaly field with a Matérn-3/2 covariance kernel over the joint (latitude, longitude, time) space. Initial hyperparameters:

| Parameter | Value | Justification |
|-----------|-------|---------------|
| Horizontal scale | 500 km (prior) | Pacific basin-scale |
| Temporal scale | 8.0 d | MODIS 8-day composite period |
| Min/max time scale | 1–100 d | Optimizer bounds |
| FitHp | true, NRestarts=3 | Optimizer re-fits all scales |

The optimizer adjusts hyperparameters by maximum likelihood; initial values are physically motivated priors only.

**Step 3 — Temporal post-processing** (`inverse_temporal_preprocess`):

The kriged anomalies are inverted by adding back the seasonal mean at each station's (latitude, longitude, time) coordinate.

### NPP uncertainty

The kriging returns a predictive standard deviation `gp15_npp.stdev`, which quantifies interpolation uncertainty. This is combined in quadrature with a 30% algorithmic uncertainty floor from the CBPM model structure (Behrenfeld et al.; Saba et al. 2011):

```
σ_NPP = √[(σ_kriging)² + (0.30 × NPP)²]
```

This is applied in `calcStats.m` when computing Ez Ratio uncertainty.

### 35-day trailing mean

Because the ²³⁴Th flux integral represents productivity over the prior ~35 d (mean lifetime = 24.101/ln 2 ≈ 34.8 d), a temporally matched NPP is required for a defensible Ez Ratio. `calcGp15Npp.m` therefore expands the kriging query to include all daily time steps in [t_sample − 34, t_sample] for each station — 35 × N_stations query points total in a single kriging call — then averages the predictions per station. The 35-day mean and its RMS kriging SD are stored as `mmolC_35d` / `stdev_35d`. `calcStats.m` uses these fields exclusively for EzRatio; the instantaneous `mmolC` / `stdev` are retained for the diagnostic NPP-at-sampling-date plot.

**Output**: `gp15_npp` table — one row per station, columns `mmolC` / `stdev` (instantaneous), `mmolC_35d` / `stdev_35d` (35-day trailing mean).

---

## 9. Stage 5 — Mixed Layer Depth

**Script**: `scripts/calculate_mld/subroutines/calcMlDepths.m`

Three MLD algorithms are computed from CTD density and temperature profiles:

| Method | Abbreviation | Criterion |
|--------|-------------|-----------|
| Pickart threshold | JAK | Density deviates by more than 1σ from mean in the subjective mixed layer; smoothing spline zero-crossing |
| Temperature difference | dBM | Surface temperature T(10 m) − ΔT = 0; T(10 m) estimated via ordinary kriging |
| Brunt-Väisälä | BW | Density exceeds surface density by fixed Δσ thresholds (multiple variants) |

The JAK MLD (`gp15_mld.('MLD JAK')`) is the primary value used downstream to define the upper boundary of the Th-234 gradient zone and the ratio regression domain.

**Spline parameters**: all three methods use `csaps` with `p = P_SPLINE_DENSITY = 0.99` (dense, precise CTD profiles) and explicit unit weights `ones(size(profile))` (no per-point CTD uncertainty available; unit weights are numerically equivalent to unweighted but make the convention explicit). Both constants are set in `setModelCoefficients.m`.

**Sensitivity analysis** (`sensitivitySplineDensity.m`): runs inside `calcMld.m` after `calcMlDepths`. Refits the JAK density spline for p ∈ {0.495, 0.99, 1.0} (p/2, nominal, 1.0) at representative equatorial and subtropical stations and reports the resulting MLD depth for each p value. Three reasons why the flux is expected to be insensitive: (1) CTD profiles are dense (sub-metre resolution) and precise — the zero-crossing location is nearly identical across this p range; (2) `P_SPLINE_DENSITY` controls a depth boundary (MLD), not a flux term directly, so a small MLD shift has only second-order effect on the integrated Th-234 flux; (3) the three independent MLD algorithms (JAK, dBM, BW) provide an implicit inter-method robustness check whose spread subsumes any p variation. Figures are saved to `plots/calcMld/sensitivity_splineDensity_stn<N>.pdf`.

**Output**: `gp15_mld` table with MLD (m) from each method per station.

---

## 10. Stage 6 — Equilibrium Point Depth (EQP)

**Script**: `scripts/calculate_eqp/subroutine/calcEqpDepths.m`

### Physical meaning

The EQP is the depth at which ²³⁴Th/²³⁸U = 1, i.e., where the ²³⁴Th deficit transitions from positive (scavenging-dominated, active flux) to negative (regeneration-dominated). It provides an upper bound on where the steady-state flux assumption is physically meaningful.

### Algorithm

Per station:

1. Compute the normalized deficit profile: `δ(z) = (²³⁴Th / ²³⁸U) − 1`. NaN and non-positive values are removed before fitting.
2. Fit `csaps` smoothing splines to `δ(z)` and `|δ(z)|` with `p = P_SPLINE_GRADIENT = 0.9` and inverse-variance weights derived from the measurement uncertainties via the delta method for division: `σ_δ = (Th/U) × √[(σ_Th/Th)² + (σ_U/U)²]`.
3. Find all zero crossings of the signed spline (via `fnzeros` on the pp struct) and the minimum of the absolute-value spline (via `fminbnd(@(z) ppval(...), ...)`) within [surface, PPZ + tolerance].
4. Select the candidate point closest to both PPZ and MLD using a weighted distance metric (equal weights W_EZ = W_MLD = 0.5):

```
score(z) = 0.5 · |z − PPZ| + 0.5 · |z − MLD|
```

The EQP is used in `doVertGradCalc.m` to define the lower boundary of the gradient taper zone.

**Diagnostic plots** (`plotEqp.m`): standalone script in `calculate_eqp/subroutine/`; run after `calcEqp.m` has saved `gp15_eqp.mat`. Loads `gp15_eqp`, `gp15_obs`, `gp15_mld`, `gp15_ppz`, and `gp15_stations` from saved files. Plots Th-234 and U-238 activity profiles with horizontal lines marking the EQP (solid), MLD (dashed), and PPZ (dash-dot) per station. One PDF per station saved to `plots/calcEqp/eqpPlots/`.

**Output**: `gp15_eqp` table with `stationNo` and `eqp` (m) per station.

---

## 11. Stage 7 — Upwelling Velocity

**Scripts**: `scripts/calculate_upwell/subroutines/`

### Data products

Six ocean circulation model products are processed (`configureUpwell.m`). Five are independent reanalyses that form the flux model ensemble; GREP is their arithmetic mean and is carried through the pipeline for diagnostic comparison only:

| Abbreviation | Full name | Institution | Role |
|-------------|-----------|-------------|------|
| ECCO | ECCO v4 | JPL/MIT | Reference + ensemble |
| CGLO | C-GLORS05 | CMCC (Italy) | Ensemble |
| FOAM | GloSea5 | Met Office (UK) | Ensemble |
| GLOR | GLORYS2V4 | Mercator Ocean (France) | Ensemble |
| ORAS | ORAS5 | ECMWF | Ensemble |
| GREP | CMEMS Global Ensemble Reanalysis | CMEMS | Diagnostic only — mean of CGLO/FOAM/GLOR/ORAS; excluded from ensemble uncertainty |

### Vertical velocity calculation (`calculateUpwelling.m`)

Vertical velocity w(z) is computed by integrating the horizontal divergence of the model velocity field from the surface down:

```
w(z) = −∫₀ᶻ (∂u/∂x + ∂v/∂y) dz'
```

**Step-by-step** (per grid point):

1. Extract zonal (u) and meridional (v) velocity profiles at the central point and at its four horizontal neighbors
2. Compute horizontal derivatives using finite differences:
   - Central differencing where neighbors exist; forward/backward at boundaries
   - Convert degree spacing to meters via `deg2km(dxDeg, A) × cos(lat)` for zonal spacing (the `cos(lat)` factor converts the zonal arc length to the true east–west distance); meridional spacing uses `deg2km(dyDeg, A)` directly. `A` is Earth's radius in metres.
3. Compute horizontal divergence: `div = ∂u/∂x + ∂v/∂y`
4. Integrate vertically via cumulative trapezoidal rule: `w(z) = cumtrapz(z, −div)`
5. Apply sign convention: positive w is upward

**Output**: 4D upwelling array `w` (longitude × latitude × depth × time) in m s⁻¹.

### Spatial and temporal smoothing (`spatialTimeAverage.m`)

The raw w field is smoothed to reduce small-scale noise before kriging:

1. **Spatial**: Moving mean over longitude axis, then latitude axis, with window size `spatAve = floor(degreeAve / spaceResolution)`, where `degreeAve` is the mean inter-station spacing (Bindoff & Wunsch 1992)
2. **Temporal**: Moving mean over the time axis with window `timeAve = floor(35 d / timeResolution)`. The 35 d window matches the ²³⁴Th mean lifetime (t½/ln 2 = 24.10/0.693 = 34.8 d), ensuring the velocity field is time-averaged over the same period that the ²³⁴Th inventory integrates.

Both operations use MATLAB's `movmean`. The temporal average must be applied to the spatially-averaged field `wSpatAve`, not the raw `w`.

**ECCO processing**: ECCO provides its vertical velocity (`wo`) directly rather than u/v fields; it therefore bypasses `calculateUpwelling.m`. After the main MERCATOR product loop, `calcUpwell.m` runs a dedicated ECCO block that loads `ecco_wo.mat`, populates the `u_mercator` coordinate interface, and passes `w = ecco_wo.wo` through `spatialTimeAverage`, `plotUpwell`, and `saveUpwell` — producing `w_ecco.nc` in the same format as the MERCATOR outputs. ECCO-specific averaging parameters: `spaceResolution = 1.0` deg (native 1° lat–lon grid, confirmed from `ecco_wo_quart*.cdf`); `timeResolution = 10` d (the ECCO `Wave` field is a 10-day average — 240 model time steps at 1 h per step — confirmed from `ecco_wo_quart*.cdf` time diffs of 240 h), giving `timeAve = floor(35/10) = 3` — a 30-day moving average, the nearest feasible approximation to the 35-d Th-234 integration window at 10-day resolution.

**Sensitivity analysis — temporal window** (`sensitivityDeltaDay.m`): runs inside the product loop after `spatialTimeAverage`. Recomputes `wSpatAve` for `deltaDay ∈ {10, 24, 35}` d and plots all three curves as a function of latitude at 100 m and the mean cruise date. Figure title includes the spatial averaging window in physical degrees (`spatAve * spaceResolution` deg) to confirm the spatial filter is held constant. One figure per product is saved to `plots/calcUpwell/sensitivity_deltaDay_<product>.pdf`.

**Sensitivity analysis — spatial window** (`sensitivitySpatialAve.m`): runs inside both the MERCATOR product loop and the ECCO block after `spatialTimeAverage`. Fixes temporal averaging at 35 d and tests four spatial window sizes — none (window = 1), ~half nominal, nominal, and ~double nominal — where nominal is `spatAve = floor(degreeAve / spaceResolution)`. For the no-averaging case (`spatAveTest = 1`), `movmean` is bypassed entirely rather than called with window 1 (mathematically equivalent, but makes the no-averaging intent unambiguous). Legend entries show physical degree extent (`spatAveTest * spaceResolution` deg) rather than grid-index window size, so MERCATOR and ECCO figures are directly comparable. One figure per product saved to `plots/calcUpwell/sensitivity_spatAve_<product>.pdf`.

**Diagnostic plots** (`plotUpwell.m`): runs inside the product loop. Produces two figures per product. Both figure titles include the spatial averaging window in physical degrees (`spatAve * spaceResolution` deg) and the 35 d temporal window. Figure 1: meridional slice of `wSpatAve` at the mean GP15 transect longitude, 100 m depth, and mean cruise date — plotted against all model latitudes; saved to `plots/calcUpwell/w_<product>.pdf`. Figure 2: `wSpatAve` extracted at the nearest grid point to each GP15 station's individual (lon, lat, 100 m, date) coordinates — plotted against station latitude with `markerSize = 10`; saved to `plots/calcUpwell/w_<product>_stations.pdf`. Figure 2 is directly comparable to the kriging output from `interpUpwell.m`.

**Ensemble comparison** (`compareUpwellModels.m`): runs after the product loop and the ECCO block. Loads `w_<product>.nc` for CGLO, FOAM, GLOR, ORAS, and GREP once each, extracting two views in the same loading pass. Grid coordinates (longitude, latitude, depth, time) are loaded directly from `w_grep.nc` rather than from `u_mercator`, which has been overwritten by the ECCO block at this point; `cmems_spaceResolution = 0.25` and `cmems_spatAve = floor(degreeAve / cmems_spaceResolution)` are computed locally for the CMEMS grid. All figure titles display the spatial averaging window in physical degrees (`cmems_spatAve * cmems_spaceResolution` deg) and the 35 d temporal window. Figure 1 (meridional slice): same mean-longitude, mean-date cross-section as `plotUpwell.m` figure 1, overlaid for all four members, their mean, and GREP; saved to `plots/calcUpwell/w_comparison.pdf`. Figure 2 (station positions): `wSpatAve` extracted at the nearest grid point to each GP15 station's individual (lon_i, lat_i, 100 m, date_i) coordinates, overlaid for all members, mean, and GREP; saved to `plots/calcUpwell/w_comparison_stations.pdf`. Figure 2 is directly comparable to the kriging output from `interpUpwell.m`.

**Output**: Smoothed `wSpatAve` in m s⁻¹, saved as netCDF4 (variable name `w`) for all six products including ECCO.

---

## 12. Stage 8 — Upwelling Interpolation

**Script**: `scripts/profile_interpolation/subroutines/interpUpwell.m`

### Purpose

Model w fields are defined on regular grids at fixed time steps. This stage interpolates w to the exact (station, depth, sampling date) coordinates of the GP15 observations using spatiotemporal ordinary kriging, run separately for each of the six model products.

### Algorithm

**Step 1 — Load and subset**

All six products (ECCO and MERCATOR) are loaded from their `calcUpwell.m` netCDF output (`w_<product>.nc`), which contains `wSpatAve` — the spatially and 35-day temporally averaged field. ECCO previously loaded raw `ecco_wo.mat` directly; it now enters via the same averaging pipeline as MERCATOR products, ensuring consistent magnitude and smoothing across all kriging training datasets.

Extract the Pacific domain: longitude [129, 280]° (129°E to 80°W — strictly within the Pacific basin; 280°E sits just west of the South American coast at all latitudes), latitude [−60, 60]° (upper bound ensures all GP15 stations near 59°N are interpolated rather than extrapolated), depth [0, 400] m (matches BTMDEPTH; data deeper than 400 m has near-zero Matérn-3/2 covariance with 0–400 m queries given `GpVerticalScaleM = 100 m`). Remove zero-velocity points (land/unfilled cells).

**Step 2 — Temporal preprocessing**

`fit_temporal_preprocessor` fits a harmonic regression (annual + semi-annual + linear trend) per depth bin [0, 100, 200, 400, ∞] m to the training w field. This removes the seasonal cycle, producing anomalies that are more suitable for kriging. The preprocessor is saved to disk for inversion.

**Step 3 — Ordinary kriging**

A single kriging call covers all (station, depth) query pairs simultaneously:

```
Query array: [lat, lon+360, depth, time]  — one row per (station, depth) pair
```

Kriging parameters (Matérn-3/2 kernel):

| Parameter | Prior value | Bounds | Role |
|-----------|-------------|--------|------|
| Horizontal scale | 500,000 m | — | Pacific spatial autocorrelation |
| Vertical scale | 100 m | — | Depth decorrelation |
| Temporal scale | 30 d | 1–365.25 d | Intraseasonal variability |
| FitHp | true | NRestarts=3 | Optimizer re-fits all scales |

The optimizer adjusts hyperparameters by maximum likelihood; the values above are physically motivated starting points only.

**Step 4 — Temporal post-processing**

`inverse_temporal_preprocess` adds back the seasonal mean at each query (station, time, depth) to recover physical velocities in m s⁻¹.

### Output table structure

```
gp15_w (ECCO) or gp15_wAve (MERCATOR):
  stationNo   — station number
  longitude   — query_lon − 360  [deg, −180 to 180]
  latitude    — [deg]
  depth       — [m]
  time        — [MATLAB datenum]
  w           — kriged vertical velocity  [m s⁻¹]
  wVar        — wStd²  [m² s⁻²]
  wErr        — kriging prediction SD  [m s⁻¹]
  wN          — 1 (single BLUP per query point)
```

**Loading**: `loadCollateData.m` assembles all six products into `gp15_w.<product>` struct fields. ECCO (`gp15_w`) is the reference product used for measurement-propagated 2D error; CGLO, FOAM, GLOR, and ORAS enter the ensemble spread; GREP is loaded for diagnostic comparison only and is excluded from the ensemble uncertainty calculation.

**Diagnostic comparison** (`compareModels.m`): produces two figures, each with two stacked panels (100 m and PPZ).

*Figure 1 — full ensemble* (`all_model_transect.pdf`): CMEMS members (CGLO, FOAM, GLOR, ORAS) in gray with distinct line styles; CMEMS mean (thick black); GREP (thinner red); ECCO (dark blue); all-model mean of five independent products (thick dark green). GREP is excluded from ensemble means because it is the arithmetic mean of CGLO/FOAM/GLOR/ORAS and is not independent.

*Figure 2 — FOAM and GREP excluded* (`all_model_transect_noFOAM.pdf`): same layout and styling, but FOAM and GREP are omitted. CMEMS mean and all-model mean are recomputed over CGLO, GLOR, ORAS only (n=3) and ECCO + CGLO + GLOR + ORAS (n=4) respectively.

**FOAM outlier**: FOAM (GloSea5) is a consistent outlier in equatorial upwelling relative to the other four products. This is not a pipeline artefact — it reflects known properties of GloSea5/NEMO: (1) NEMOVAR 3D-Var assimilation corrects T/S/SSH but not velocities directly; near-equatorial ageostrophic flow is poorly constrained, producing spurious divergence; (2) NEMO has documented Pacific equatorial undercurrent biases (core too deep, too weak) that directly inflate the diagnosed w; (3) z\*-coordinate free-surface corrections can introduce residual errors in the vertically integrated w (Storto et al. 2019, *Ocean Sci.*). GREP is excluded from the sensitivity figure because it incorporates FOAM and is not independent of the outlier. FOAM is the primary driver of σ_model; its inclusion makes the reported uncertainty conservative.

**Ensemble mean columns** (`collateStations.m`): after the per-product join, `w_meanFull` / `wErr_meanFull` (n=5, all independent products) and `w_meanNoFOAM` / `wErr_meanNoFOAM` (n=4, FOAM excluded) are computed with `wErr_mean = √(Σ σ_k²)/n`. Total-error columns `totalXxxFluxError_meanFull` and `_meanNoFOAM` are stored in `calcError.m`. The FOAM-excluded sensitivity is reported by setting `fluxProduct = uncertProduct = 'meanNoFOAM'` in `calcStats.m` alongside `all_model_transect_noFOAM.pdf`.

---

## 13. Stage 9 — Vertical Th-234 Gradient

**Script**: `scripts/calculate_gradient/subroutines/doVertGradCalc.m`

### Purpose

The 2D upwelling correction requires ∂²³⁴Th/∂z, the vertical gradient of total ²³⁴Th activity. This is estimated from the discrete profile measurements using finite differences.

### Algorithm (per station)

**Step 1 — Finite differences**

For each observation depth i:
- Surface (i=1): Forward difference — `grad_i = (Th_{i+1} − Th_i) / (z_{i+1} − z_i)`
- Bottom (i=n): Backward difference — `grad_i = (Th_i − Th_{i-1}) / (z_i − z_{i-1})`
- Interior: Centered difference — `grad_i = (Th_{i+1} − Th_{i-1}) / (z_{i+1} − z_{i-1})`

**Step 2 — Measurement uncertainty propagation**

The standard deviation of each finite-difference gradient is propagated from the ²³⁴Th measurement uncertainties. For all three schemes the result has the same form with B = depth interval:

```
σ_grad,i = √[(σ_Th,a / B)² + (σ_Th,b / B)²]
```

**Step 3 — Hard-zero outside active zone**

Both the gradient and its uncertainty are set to zero outside `[mldVal, zBot]` where `zBot = min(EQP, PPZ)`:

```
vertGrad      = 0   if depth < MLD  or  depth > min(EQP, PPZ)
vertGradError = 0   if depth < MLD  or  depth > min(EQP, PPZ)
```

Zeroing (rather than NaN-ing) allows cumulative sums in `doModelCalcs.m` to work correctly. A smoothing spline with boundary anchors was evaluated and removed: for GP15 Th-234 profiles (8–15 observations per station), the spline reproduced the raw finite differences within errorbars at all stations while the combined residual + measurement uncertainty inflated `vertGradError` near the anchor points with no scientific gain.

**Output**: `gp15_grad` table with columns `vertGrad` and `vertGradError` (dpm L⁻¹ m⁻¹, converted to dpm m⁻³ m⁻¹ in `collateStations.m`).

---

## 14. Stage 10 — Regional Delineation

**Script**: `scripts/regional_delineation/subroutines/calcRegionalData.m`

The GP15 transect is divided into four biogeographic provinces based on hydrographic and biogeochemical properties:

1. North Pacific High Productivity Zone
2. North Pacific Gyre
3. Equatorial Pacific
4. South Pacific Gyre

Province boundaries are defined in `regions.xlsx` and loaded as `regions.first` / `regions.last` station numbers.

**Data used for characterization**: CTD temperature, transmissometer (bulk particle abundance), bottle NO₃ and PO₄ (quality-flagged == 1 or 2), particulate nitrogen, and HPLC pigment fractions (% Micro, Nano, Pico).

**Output**: `gp15Regions` table assigning a region label to each (station, depth) pair; `regions` struct with province boundaries.

---

## 15. Stage 11 — POC:Th and PN:Th Ratio Regression

**Scripts**: `scripts/regress_ratio/subroutines/`

### Motivation

The POC:²³⁴Th and PN:²³⁴Th particulate ratios are measured at discrete depths. A smooth, depth-dependent ratio profile is needed to convert Th-234 fluxes to POC/PN fluxes throughout the water column.

### Fitting function (`src/toolboxes/regression/fitPiecewiseRatio.m`)

Each of the four ratio types (POC LSF, POC SSF, PN LSF, PN SSF) is fit independently per biogeographic region by `fitPiecewiseRatio`, called from the corresponding `calcRegionalRegression*.m` script. The model is:

```
R(z) = α + β · clamp(z, MLD, z*)     β ≤ 0

where clamp(z, lo, hi) = max(min(z, hi), lo)
```

This gives a **constant surface cap** (ratio held at the MLD value above MLD), a **linearly decreasing** active zone from MLD to z\*, and a **constant deep zone** below z\* — all continuous by construction.

**Step 1 — Equal-count depth binning.** `n_bins = max(5, min(10, floor(n_raw / 5)))` bins are formed by rank-based assignment so each bin holds approximately equal numbers of raw observations regardless of depth ties.

**Step 2 — Per-station-per-bin means.** For each (station, bin) combination, an inverse-variance weighted mean and its propagated uncertainty are computed. These `n_stations × n_bins` means are the observations fed to the fit, removing within-station serial correlation.

**Step 3 — Profile likelihood over z\*.** Candidate values are the unique raw measurement depths in [MLD, 400 m]. At each candidate z\*: (a) the design matrix is X = [1, clamp(z, MLD, z\*)]; (b) the unconstrained WLS solution is computed analytically; (c) if β > 0 violates the constraint, the solution is clamped to the boundary β = 0 (flat line, α = weighted mean). The residual sum of squares RSS(z\*) is recorded.

**Step 4 — Optimal parameters and covariance.** The z\* minimising RSS is selected. The parameter covariance is:

```
paramCov = sSquared · (X'WX)⁻¹,    sSquared = RSS / (n_bin_means − 2)
```

`sSquared` is the residual variance — it scales `paramCov` to account for scatter of the bin means around the fitted line, beyond what the per-bin weights alone capture. The flat-line degenerate case (β = 0) uses a 1-parameter form with `sSquared = RSS / (n − 1)`. `paramCov` elements are `alphaVar = paramCov(1,1)`, `betaVar = paramCov(2,2)`, `coVar = paramCov(1,2)`.

**Step 5 — Profile-likelihood CI on z\*.** The 95% CI is the set of z\* where `n · log(RSS(z\*)/RSS_opt) ≤ 3.841` (χ²₀.₉₅, 1 dof). Stored as `zstarCI_lo` / `zstarCI_hi` in the output table. This CI is only on z\* itself; it is not propagated into the flux error.

**Transition depth z\*** is found by the optimizer for all regions.

### Three-zone evaluation (`makeGp15RegressRatio.m`)

For each (station, depth) pair, the predicted ratio and its **1σ delta-method SE** are evaluated:

| Zone | Depth | Predicted R | 1σ uncertainty of R |
|------|-------|-------------|---------------------|
| Surface cap | z ≤ MLD | α + β·MLD | `√(β_var·MLD² + α_var + 2·MLD·cov)` |
| Active | MLD < z ≤ z\* | α + β·z | `√(β_var·z² + α_var + 2z·cov)` |
| Deep | z > z\* | α + β·z\* | `ratioSubsurfaceStd` = `√(β_var·z*² + α_var + 2·z*·cov)` |

All three expressions are the **delta-method SE of the predicted regression mean** at that depth: `σ_R(z) = √(x(z)ᵀ · paramCov · x(z))` where `x(z) = [1, clamp(z, MLD, z*)]`. This is a confidence interval on the fitted mean, not a prediction interval (it does not add a new-observation scatter term).

The result is stored as `uncertPocTh234RatioRegress` (and PN analogue). It is a **1σ standard error**, propagated as such in `calcError.m`.

**Plot label**: the shaded uncertainty band in `calcRegionalRegression*.m` diagnostic figures is labelled `$\pm 1\sigma$`, consistent with the pipeline's 1σ error propagation throughout.

**Output**: `regressPocRatio` and `regressPnRatio` tables — one row per (station, depth) with predicted ratio and 1σ uncertainty `uncertPocTh234RatioRegress`.

---

## 16. Stage 12 — Data Collation

**Scripts**: `scripts/collate_data/subroutines/`

### Loading (`loadCollateData.m`)

Assembles all preprocessed data into memory:
- `gp15_obs` — radionuclide and particulate profiles
- `gp15_w.<product>` — kriged upwelling for all 6 products (struct with one field per product)
- `gp15_grad` — vertical Th-234 gradient
- `regressPocRatio`, `regressPnRatio` — depth-dependent ratio regression
- `gp15Regions` — province assignments
- `gp15_mld`, `gp15_eqp`, `gp15_ppz`, `gp15_stations` — depth diagnostics

### Assembly (`collateStations.m`)

Per station, the following operations are performed and joined on depth as the key:

1. Subset `gp15_obs` to depths ≤ 400 m
2. Extract upwelling w and wErr for ECCO and all five MERCATOR products
3. Extract vertical gradient and gradient error
4. Extract POC:Th and PN:Th regression ratios
5. Extract province assignment

**Unit conversions applied**:

```
Radionuclides:  × 1000   (dpm L⁻¹ → dpm m⁻³)
Particulates:   ÷ 1000   (µmol → mmol)
Upwelling:      × 86400  (m s⁻¹ → m d⁻¹)
Gradient:       × 1000   (dpm L⁻¹ m⁻¹ → dpm m⁻³ m⁻¹)
Ratios:         ÷ 1000   (µmol dpm⁻¹ → mmol dpm⁻¹)
```

**Ensemble mean columns**: after the per-product unit conversions, `collateStations.m` computes `w_meanFull`, `wErr_meanFull` (mean and propagated kriging uncertainty of 5 independent products) and `w_meanNoFOAM`, `wErr_meanNoFOAM` (4 products, FOAM excluded). Kriging uncertainty propagates as `wErr_mean = √(Σ σ_k²)/n`.

**Zero-velocity rule**: Upwelling is set to zero wherever the vertical gradient is zero. A depth-alignment assert (`isequal(statUpwell.depth, statVertGrad.depth)`) fires before this operation to confirm both tables share the same depth grid. The zeroing loop uses named field access over all eight products (six individual + `meanFull` + `meanNoFOAM`). This ensures no spurious 2D correction at depths outside the active gradient zone.

**Duplicate-depth guard**: Before any join executes, six per-station uniqueness assertions (lines 43–49) verify that each input table has no repeated depths. If an assert fires, the pipeline halts immediately with the station number and the offending source table named, allowing the upstream issue to be fixed at the source.

**Output**: `gp15_inputs` — a flat table with one row per (station, depth) pair containing all variables needed by the flux model.

---

## 17. Stage 13 — Th-234 Flux Model

**Scripts**: `scripts/model/`

### Entry point (`runModel.m`)

Calls subroutines in order: `loadModelData` → loop over stations calling `doModelCalcs` → `calcParticulateFlux` → `calcError` → `makeDepthArrays` → `calcStats` → `saveModelData`.

### Coefficients (`setModelCoefficients.m`)

```
LAMBDA = ln(2) / 24.101 d⁻¹ = 0.02876 d⁻¹
```

This is the ²³⁴Th decay constant. The half-life 24.101 d is the literature value and must not be changed.

### Per-station flux calculation (`doModelCalcs.m`)

**Depth cell widths**: A half-open cell scheme is used to ensure the surface integral is correct:

```
dz_i = (z_{i+1} − z_i)/2 + (z_i − z_{i-1})/2 · bc_i

where bc_i = 2 for i=1 (surface), bc_i = 1 otherwise
```

This gives the surface cell width as the full distance to the next depth, and interior cells as the average of the two half-spacings.

**1D Th-234 flux** (accumulated from surface to each depth):

```
F_Th,1d(z_k) = Σᵢ₌₁ᵏ λ · [²³⁸U_i − ²³⁴Th_i] · dz_i     [dpm m⁻² d⁻¹]
```

The deficit [²³⁸U − ²³⁴Th] is in dpm m⁻³; multiplying by dz (m) gives dpm m⁻²; multiplying by λ (d⁻¹) gives dpm m⁻² d⁻¹.

**2D Th-234 flux** (for all 6 model products):

```
F_Th,2d(z_k) = Σᵢ₌₁ᵏ [λ · (²³⁸U_i − ²³⁴Th_i) + w_i · (∂Th/∂z)_i] · dz_i
```

where w (m d⁻¹) × ∂Th/∂z (dpm m⁻³ m⁻¹) × dz (m) = dpm m⁻² d⁻¹.

### Particulate flux calculation (`calcParticulateFlux.m`)

For each model product and for the 1D case:

```
F_POC = F_Th × R_POC(z, region)    [mmol C m⁻² d⁻¹]
F_PN  = F_Th × R_PN(z, region)     [mmol N m⁻² d⁻¹]
```

where R is the depth- and region-specific ratio from the regression (mmol dpm⁻¹).

### Standard depth extraction (`makeDepthArrays.m`)

After model calculation, fluxes are extracted at three standard depths. All matching uses a tolerance of 10⁻⁶ m to guard against floating-point rounding:

| Variable | Depth | Method |
|----------|-------|--------|
| `gp15_flux100` | 100 m | `find(|depth − 100| ≤ 1e−6)` |
| `gp15_fluxPpz` | PPZ | Per-station loop matching `depthPpz` within 1e−6 |
| `gp15_flux100Ppz` | PPZ + 100 m | Per-station loop matching `depthPpz + 100` within 1e−6 |

Each extraction asserts exactly NUMSTAT rows. If an assert fires, the most likely cause is that the PPZ or PPZ+100 depth was not inserted as a query point during station-profile interpolation (`interpData`).

---

## 18. Stage 14 — Error Propagation

**Script**: `scripts/model/subroutines/calcError.m`

Error propagation is performed per station via loops over the six model products using dynamic field access (`s.(['field_' dataProduct])`). All intermediate quantities are **variances** (σ²); square roots are taken only at the final storage step.

### 1D flux variance

```
σ²_1d(z_k) = Σᵢ₌₁ᵏ (λ · dz_i)² · (σ²_U238,i + σ²_Th234,i)
```

Both σ_U238 and σ_Th234 are measured standard deviations from the GP15 cruise data.

### 2D upwelling correction variance (per model product)

At each depth layer i:

```
σ²_upwell,i = (g_i · dz_i)² · σ_w,i² + (w_i · dz_i)² · σ_g,i²
```

where:
- w_i = upwelling velocity (m d⁻¹) for the given product
- σ_w,i = `wErr` from kriging (m d⁻¹)
- g_i = `vertGrad` (dpm m⁻³ m⁻¹)
- σ_g,i = `vertGradError` (dpm m⁻³ m⁻¹)
- dz_i = layer thickness from `calcLayerThickness` (m)

This absolute-error form is always valid and avoids division by zero when w or g is zero. NaN values in `σ²_upwell` propagate through `cumsum` — a missing upwelling correction at one depth renders the cumulative 2D flux variance undefined at all greater depths.

### 2D total flux variance

```
σ²_2d(z_k) = Σᵢ₌₁ᵏ [(λ · dz_i)² · (σ²_U238,i + σ²_Th234,i) + σ²_upwell,i]
```

### Particulate flux variance

For POC (analogous for PN):

```
σ²_POC(z) = F²_POC(z) · [σ²_2d(z) / F²_Th,2d(z) + σ²_R_POC(z) / R²_POC(z)]
```

This is the delta method for a product of two uncertain quantities.

### Upwelling correction uncertainty

The upwelling correction ΔF = F_2D − F_1D has its own uncertainty. Because F_2D = F_1D + ΔF and the error sources are independent — F_1D depends only on σ_U238 and σ_Th234, while ΔF depends only on σ_w and σ_grad — their variances add:

```
Var[F_2D] = Var[F_1D] + Var[ΔF]   →   Var[ΔF] = Var[F_2D] − Var[F_1D]
```

Therefore:

```
σ_correction = sqrt(max(0, sigma_2D² − sigma_1D²))
```

The `max(0, ...)` guard prevents tiny negative values from floating-point rounding. Do **not** use `sigma_2D − sigma_1D` (difference of standard deviations) — this understates the uncertainty and is dimensionally inconsistent with error propagation.

Stored as `uncertUpwellCorrect_<product>`, `uncertUpwellPocCorrect_<product>`, `uncertUpwellPnCorrect_<product>` in `gp15_flux`. These are always non-negative by construction.

### Model ensemble uncertainty

The spread across the five independent 2D flux estimates (ECCO, CGLO, FOAM, GLOR, ORAS) quantifies model-structural uncertainty. GREP is excluded because it is the arithmetic mean of CGLO, FOAM, GLOR, and ORAS and is not an independent realisation. The inter-model standard deviation is used directly — there is no `√n` divisor, because the models are deterministic reanalyses, not i.i.d. random samples:

```
σ_model,full    = std(F_Th,2d over ECCO, CGLO, FOAM, GLOR, ORAS)       n=5
σ_model,noFOAM  = std(F_Th,2d over ECCO, CGLO, GLOR, ORAS)             n=4
```

σ_model is attached **only to the mean-ensemble products** (`meanFull`, `meanNoFOAM`). Individual members carry measurement-propagated uncertainty (σ_2d) only. FOAM is the primary driver of σ_model at equatorial stations; its inclusion makes the full-ensemble uncertainty conservative.

### Total reported uncertainty

The measurement-propagated uncertainty and the model-structural uncertainty are combined in quadrature. Three named total-error columns are stored in `gp15_flux` (and propagated to `gp15_flux100`, `gp15_fluxPpz`, `gp15_flux100Ppz` via `makeDepthArrays.m`):

| Column | σ_2d from | σ_model from | Use |
|--------|-----------|--------------|-----|
| `totalTh234FluxError` (legacy) | ECCO wErr | full n=5 | backward compat only |
| `totalTh234FluxError_ecco` | ECCO wErr | full n=5 | alias for legacy column |
| `totalTh234FluxError_meanFull` | meanFull wErr | full n=5 | **default reported** |
| `totalTh234FluxError_meanNoFOAM` | meanNoFOAM wErr | no-FOAM n=4 | FOAM sensitivity |

Analogous columns exist for POC (`totalPocFluxError_*`) and PN (`totalPnFluxError_*`).

`calcStats.m` reads `fluxProduct` and `uncertProduct` variables (both default `'meanFull'`) and constructs column names as `['th234FluxCumul2d_' fluxProduct]` / `['totalTh234FluxError_' uncertProduct]`. This makes switching the reported product or uncertainty a one-line change.

---

## 19. Stage 15 — Metrics and Statistics

**Script**: `scripts/model/subroutines/calcStats.m`

### Per-station flux table

For each station, fluxes and uncertainties are extracted at 100 m, PPZ, and PPZ+100 m for both the 1D model and (where applicable) the 2D correction.

**Equatorial 2D correction**: The 2D upwelling correction is applied only at stations within the equatorial band (|latitude| ≤ 5°). Outside this band, the 2D upwelling flux is set to NaN and the 1D value is reported. The latitude threshold is data-driven and geophysically defined.

### Derived metrics

**Transfer efficiency** (T₁₀₀):

```text
T₁₀₀ = F_POC(PPZ+100) / F_POC(PPZ)

σ_T100 = √[max(σ²_Th,PPZ+100 − T₁₀₀·(2−T₁₀₀)·σ²_Th,PPZ, 0)] / |F_Th(PPZ)|
````

T₁₀₀ is defined as a POC flux ratio. However, when both depths use the same POC/²³⁴Th ratio `R`,

```text
F_POC(z) = F_Th(z) · R
```

so

```text
T₁₀₀ = F_Th(PPZ+100)·R / [F_Th(PPZ)·R]
     = F_Th(PPZ+100) / F_Th(PPZ).
```

Therefore, uncertainty in `R` cancels from `T₁₀₀`. The uncertainty in `T₁₀₀` should be propagated from the ²³⁴Th flux uncertainties, not from the POC flux uncertainties.

Let

```text
A = F_Th(PPZ+100)
B = F_Th(PPZ)
T₁₀₀ = A / B.
```

The deeper flux is not independent of the PPZ flux. Physically, it is the PPZ flux after additional attenuation over the next 100 m:

```text
A = B − L,
```

where `L` is the flux loss between `PPZ` and `PPZ+100`. Thus the two fluxes share the same dominant PPZ-flux uncertainty. We approximate

```text
Cov(A, B) ≈ Var(B)
```

because the uncertainty associated with the additional attenuation term `L` is treated as a smaller correction than the uncertainty in the PPZ flux itself.

Using the delta-method formula for a ratio,

```text
Var(A/B) ≈ [Var(A) + T₁₀₀² Var(B) − 2T₁₀₀ Cov(A,B)] / B².
```

Substituting `Cov(A,B) ≈ Var(B)` gives

```text
Var(T₁₀₀) ≈ [Var(A) − T₁₀₀·(2−T₁₀₀)·Var(B)] / B².
```

In code notation,

```text
Var(T₁₀₀) =
[σ²_Th,PPZ+100 − T₁₀₀·(2−T₁₀₀)·σ²_Th,PPZ] / F_Th(PPZ)².
```

The `max(..., 0)` guard is included only to prevent tiny negative values from floating-point roundoff. If this guard frequently clamps the value to zero, that indicates a problem with the variance/covariance assumptions or with using the wrong variances. In particular, POC variances should not be used here because they include uncertainty from `R`, which cancels out of `T₁₀₀`.

**Th-234 remineralization** (R₁₀₀):

```
R₁₀₀ = F_Th(PPZ+100) − F_Th(PPZ)

σ_R100 = √[max(σ²_Th,PPZ+100 − σ²_Th,PPZ, 0)]
```

The subtraction rather than addition corrects for the same shared-layer covariance: `Cov[F(PPZ+100), F(PPZ)] = Var[F(PPZ)]`, so `Var[R₁₀₀] = Var[F(PPZ+100)] − Var[F(PPZ)]`. The `max(..., 0)` guard prevents negative values from floating-point rounding.

**Export efficiency** (Ez Ratio):

```
Ez Ratio = F_POC(PPZ) / NPP

σ_Ez = |Ez| · √[(σ_POC,PPZ / F_POC,PPZ)² + (σ_NPP / NPP)²]
```

where NPP uncertainty includes both the kriging predictive SD and the 30% CBPM algorithmic floor:

```
σ_NPP = √[(σ_kriging)² + (0.30 × NPP)²]
```

### Regional aggregation

For each of the four provinces, station-level fluxes and metrics are averaged (nanmean) and the inter-station standard deviation is reported. Regional statistics include MLD, PPZ, NPP, Th-234 fluxes, POC/PN fluxes, T₁₀₀, R₁₀₀, and Ez Ratio.

### Output tables

| Table | Contents |
|-------|----------|
| `depthData` | Per-station fluxes at all standard depths with uncertainties and derived metrics |
| `regionalData` | Province-averaged summary statistics |
| `table1` | Manuscript Table 1 (per-station Th-234 fluxes) |
| `table2` | Manuscript Table 2 (regional summary) |
| `textStats` | Numerical values referenced in results text |

All tables are saved as both `.mat` and `.xlsx`.

---

## 20. Toolboxes

### Ordinary and weighted least-squares (`src/toolboxes/regression/ols.m`)

Fits a simple linear model y = α + β·x with optional inverse-variance weighting.

**Inputs**: `x` (predictor), `y` (response), `y_uncertainty` (optional measurement SD for WLS)

**Algorithm**:

Design matrix `X = [1, x]`. If WLS, weight matrix `W = diag(1/σ²_y)`.

```
θ = (X'WX)⁻¹ X'Wy        (WLS)
θ = (X'X)⁻¹ X'y           (OLS)
Cov(θ) = σ²_resid · (X'WX)⁻¹
```

**Outputs**:
- `beta`, `alpha` — slope and intercept
- `q` — standard error of the fitted mean at each x
- `ci = [varBeta, varAlpha, covAlphaBeta, n]` — for downstream delta-method propagation

### Spatiotemporal kriging (`src/toolboxes/terraformer/`)

Python CLI wrappers called from MATLAB via `system()`. Three entry points:

| Function | Role |
|----------|------|
| `fit_temporal_preprocessor` | Fit depth-binned harmonic regression; return anomalies |
| `ordinary_kriging` | Matérn-3/2 spatiotemporal BLUP with hyperparameter optimization |
| `inverse_temporal_preprocess` | Add back seasonal mean at query coordinates |

Requires a functioning Python environment with `terraformer` installed. Verify with `python -m terraformer._cli.ok --help` before running.

### Piecewise-linear interpolation with uncertainty (`src/toolboxes/interp1dError/interp1dError.m`)

Interpolates a profile y(x) with uncertainties σ_y(x) to new query points using piecewise-linear interpolation. Uncertainties are propagated linearly through the interpolation weights.

---

## 21. Numerical Invariants and Fragile Areas

### Invariants that must remain consistent

| Quantity | Value | Files where it appears |
|----------|-------|----------------------|
| `LAMBDA` | ln(2)/24.101 d⁻¹ | `setModelCoefficients.m` |
| `BTMDEPTH` | 400 m | `configureInterp.m`, `collateStations.m` |
| `FLIERDATA` | 30 µmol dpm⁻¹ | `qcGp15Poc.m` |
| `DepthBinEdges` | [0, 100, 200, 400, ∞] m | `interpUpwell.m`, `calcGp15Npp.m` |
| Longitude offset | +360 on query, −360 on save | `interpUpwell.m` |
| Ensemble products | ECCO, CGLO, FOAM, GLOR, ORAS (n=5); GREP excluded | `calcError.m` lines 114–124 |
| `deltaDay` | 35 d | `configureUpwell.m`; matches ²³⁴Th mean lifetime |
| `P_SPLINE_GRADIENT` | 0.9 | `setModelCoefficients.m`; used in `calcEqpDepths.m` |
| `P_SPLINE_DENSITY` | 0.99 | `setModelCoefficients.m`; used in `calcMlDepths.m` |
| `BTMDEPTH` in `fitPiecewiseRatio` | 400 m | `fitPiecewiseRatio.m` — must equal the pipeline constant |
| `nDaysWindow` in `calcGp15Npp` | 35 d | `calcGp15Npp.m` — matches ²³⁴Th mean lifetime; must equal `deltaDay` |
| `dz` formula | Half-open midpoint rule; `bc=2` at surface | `calcLayerThickness.m` — called from `doModelCalcs.m` and `calcError.m`; must not be reimplemented inline |

### Known fragile areas

| File | Risk | What to check |
|------|------|---------------|
| `makeDepthArrays.m` | Medium | Asserts fire if PPZ or PPZ+100 were not inserted during `stationInterp.m`; run `interpData` first |
| `calcError.m` | Low | Both `yw`/`nw` loading loops and the error-calculation loop must all cover the same 8 products `{ecco, cglo, foam, glor, oras, grep, meanFull, meanNoFOAM}`; correction uncertainty is `sqrt(max(0, sigma_2D^2 - sigma_1D^2))` — do not revert to `sigma_2D - sigma_1D`; `upwellError` is `(g·dz)²·σ_w² + (w·dz)²·σ_g²`; NaN propagates through `cumsum`; ensemble std uses 5 independent products, no `/ sqrt(n)` |
| `collateStations.m` | Low | Table `join()` on depth key; six per-station uniqueness asserts catch duplicates loudly before any join executes; depth-alignment assert (`isequal(statUpwell.depth, statVertGrad.depth)`) before positional upwelling zeroing — do not remove; all unit conversions use named variable access (not column indices); upwelling zeroed where gradient is zero via named-access loop |
| `interpUpwell.m` | Medium | Requires `terraformer` Python environment; hyperparameters are priors — optimizer takes over |
| `calcStats.m` | Low | Equatorial 2D correction uses `abs(statLat) > 5`; all 8 occurrences must remain consistent. `uncertT100` uses Th-234 flux variances — ratio R cancels exactly in the delta method for the deep-zone case; reverting to POC variances produces a negative argument (clamped to zero by `max`). |
| `doVertGradCalc.m` | Low | Raw FD gradient and measurement-propagated SD hard-zeroed outside `[mldVal, zBot]`. `mldVal` looked up by station number (`gp15_mld.('Station No') == sn`) — do not revert to positional `iStat` indexing. |
| `qcGp15Poc.m` | Low | Flier index must be captured before NaN-ification; the `flierIdxLarge` / `flierIdxSmall` capture pattern must be preserved during refactoring |
| `fitPiecewiseRatio.m` | Low | z\* search bounded to [MLD, 400 m]; `flatLine=true` when β = 0 (legitimate degenerate case for depth-independent regions — not an error) |
| `calcGp15Npp.m` | Low | `tqWindowFlat` constructed by station-major flatten of `tqWindowMat.'`; `nppWindowMat = reshape(VhatWindowMean, nDaysWindow, nStations)` columns must align with stations — do not transpose |

### Per-product dynamic field access

Four scripts loop over model products by name using MATLAB dynamic struct/table field access: `calcError.m`, `doModelCalcs.m`, `calcParticulateFlux.m`, and `collateStations.m`. The pattern is:

```matlab
s.(['fieldName_' dataProduct]) = ...
```

When adding variables that require per-product treatment, follow this pattern and add the new product name to the cell array at the top of each loop.

---

*For dev rules, fragile areas, and the verification checklist, see `dev/project_reference.md`.*

---

## 22. Function Inventory

All `.m` files in the project organised by pipeline stage. **Role** column: E = entry point, S = subroutine, U = utility, T = toolbox.

### Top-level entry points

| File | Role | Key outputs |
|------|------|-------------|
| `scripts/gp15Model.m` | E | — (master orchestration) |
| `scripts/model/runModel.m` | E | `gp15_flux*` tables |
| `scripts/profile_interpolation/interpData.m` | E | interpolated `gp15_obs`, `gp15_w.*` |
| `scripts/quality_control/doQc.m` | E | `gp15_obs`, `gp15_obsNoQc` |
| `scripts/quality_control/plotQcPoc.m` | E | `plots/doQc/stationPlots/poc/` — standalone; loads `gp15_obsNoQc`; uncertainties as horizontal error bars |
| `scripts/quality_control/plotQcPn.m` | E | `plots/doQc/stationPlots/pn/` — standalone; loads `gp15_obs`; uncertainties as horizontal error bars |
| `scripts/read_data/readData.m` | E | `gp15_obs`, `ecco_wo` |
| `scripts/collate_data/collateData.m` | E | `gp15_inputs` |
| `scripts/calculate_ppz/calcPpz.m` | E | `gp15_ppz` |
| `scripts/calculate_mld/calcMld.m` | E | `gp15_mld` |
| `scripts/calculate_mld/plotMld.m` | E | `plots/calcMld/mldPlots/` — standalone; loads `gp15_mld` |
| `scripts/calculate_eqp/calcEqp.m` | E | `gp15_eqp` |
| `scripts/calculate_eqp/subroutine/plotEqp.m` | E | `plots/calcEqp/eqpPlots/` — standalone; loads `gp15_eqp` |
| `scripts/calculate_npp/calcNpp.m` | E | `gp15_npp` |
| `scripts/calculate_upwell/calcUpwell.m` | E | `w_<model>.nc` |
| `scripts/calculate_gradient/calcVertGrad.m` | E | `gp15_grad` |
| `scripts/calculate_gradient/plotVertGrad.m` | E | `plots/calcVertGrad/gradPlots/` — standalone; loads `gp15_grad`; horizontal errorbars (raw FD and zeroed FD, no caps); MLD and zBot reference lines; per-station title and legend |
| `scripts/regional_delineation/delineateRegions.m` | E | `regions` |
| `scripts/regress_ratio/regressRegionalRatio.m` | E | `ratioPocLsf`, `ratioPnLsf`, `ratioPocSsf`, `ratioPnSsf` |

### Model subroutines (`scripts/model/subroutines/`)

| File | Key math |
|------|----------|
| `setModelCoefficients.m` | λ = ln2/t½; P_SPLINE_GRADIENT = 0.9; P_SPLINE_DENSITY = 0.99 |
| `loadModelData.m` | — |
| `doModelCalcs.m` | `calcLayerThickness` for dz; F₁d = cumsum(λ·deficit·dz); F₂d adds w·∂Th/∂z·dz; dynamic field access |
| `calcParticulateFlux.m` | F_POC = F_Th × R_POC; dynamic field access |
| `calcError.m` | absolute-error upwelling variance; `yw`/`nw` loading loops + error loop all cover 8 products; correction uncertainty `sqrt(max(0, σ²_2d − σ²_1d))`; cumsum (NaN propagates); inter-model std (n=5, no √n); dynamic field access |
| `makeDepthArrays.m` | tolerance-based depth index; hard asserts |
| `calcStats.m` | T₁₀₀, R₁₀₀ (covariance-corrected), EzRatio (35d NPP); regional means |
| `saveModelData.m` | — |
| `clearEnvironment.m` | — |
| `plotModel.m` | Per-product bar plots: 3 rows × 2 cols (1D flux, upwelling correction, 2D flux) × (Th-234, POC); one figure per product per depth (16 files); white bars, black no-cap error bars. Comparison line plots: meanFull, meanNoFOAM, and CMEMS-vs-GREP ensembles; `compareUpwellModels.m` styling (ECCO blue, CMEMS gray/distinct styles, mean thick black, GREP red); one figure per comparison type per depth (6 files). Zero-reference `yline` has `HandleVisibility='off'` to suppress it from legends. POC flux is NaN at stations without valid POC:Th ratio data; comparison lines use `ok=~isnan` filter so each line is continuous over its valid stations. |

### Interpolation subroutines (`scripts/profile_interpolation/subroutines/`)

| File | Description |
|------|-------------|
| `configureInterp.m` | BTMDEPTH = 400 m; Pacific domain [129,280]°E × [−60,60]° × [0,400] m; product list |
| `stationInterp.m` | `interp1dError` to {MLD, PPZ, PPZ+100, EQP, 100 m} |
| `interpUpwell.m` | spatiotemporal ordinary kriging of w via `terraformer`; all products load from `w_<product>.nc` |
| `compareModels.m` | two figures: (1) full ensemble — CMEMS gray, CMEMS mean black, GREP red, ECCO dark blue, all-model mean dark green; (2) FOAM+GREP excluded sensitivity — same styling over CGLO/GLOR/ORAS only; FOAM is primary outlier and driver of σ_model |

### Quality control (`scripts/quality_control/`)

| File | Description |
|------|-------------|
| `doQc.m` | entry point |
| `subroutines/doGp15Qc.m` | flags, calls sub-subroutines |
| `subroutines/makeStations.m` | builds `gp15_stations` lookup |
| `subroutines/subsubroutines/aveDuplicate.m` | averages duplicates; σ_combined = √(Σσᵢ²)/n |
| `subroutines/subsubroutines/qcGp15Poc.m` | ratio calc, flier removal (> 30 µmol dpm⁻¹) |
| `subroutines/subsubroutines/qcGp15Pn.m` | same for PN |

### Auxiliary calculation subroutines

| File | Method |
|------|--------|
| `calculate_ppz/subroutines/calcPpzDepths.m` | `findppz`; fluorescence threshold |
| `calculate_mld/subroutines/calcMlDepths.m` | JAK/dBM/BW; `csaps(P_SPLINE_DENSITY)` + `fnzeros` |
| `calculate_mld/subroutines/sensitivitySplineDensity.m` | refits JAK density spline for p ∈ {0.495, 0.99, 1.0} at equatorial and subtropical stations |
| `calculate_eqp/subroutine/calcEqpDepths.m` | `csaps(P_SPLINE_GRADIENT, inv-var wts)` + `fnzeros`; weighted distance to PPZ/MLD |
| `calculate_npp/subroutines/calcGp15Npp.m` | spatiotemporal ordinary kriging; 35-day trailing mean |
| `calculate_upwell/subroutines/calculateUpwelling.m` | finite-difference divergence + `cumtrapz`; `cos(lat)` zonal correction |
| `calculate_upwell/subroutines/spatialTimeAverage.m` | `movmean` spatial then temporal (35 d) |
| `calculate_upwell/subroutines/sensitivityDeltaDay.m` | recomputes wSpatAve for deltaDay ∈ {10, 24, 35} d at 100 m, mean cruise date; title shows spatial window in physical degrees |
| `calculate_upwell/subroutines/sensitivitySpatialAve.m` | fixes temporal avg at 35 d; tests 4 spatial windows (none, ~half, nominal, ~double) at 100 m; legend in physical degrees (`spatAveTest * spaceResolution`); no-averaging case bypasses `movmean` |
| `calculate_upwell/subroutines/compareUpwellModels.m` | two figures: meridional slice and per-station nearest-point extraction; members + ensemble mean + GREP; titles in physical degrees; grid loaded from `w_grep.nc` (not `u_mercator`) |
| `calculate_gradient/subroutines/doVertGradCalc.m` | forward/central/backward FD; hard-zero outside `[MLD, min(EQP,PPZ)]` |
| `calculate_gradient/subroutines/sensitivitySplineGrad.m` | refits gradient spline for p ∈ {0.45, 0.9, 1.0} at equatorial and subtropical stations — vestigial; spline no longer in main path |
| `regional_delineation/subroutines/calcRegionalData.m` | per-station means of T, nutrients, HPLC above PPZ |
| `regress_ratio/subroutines/calcRegionalRegression*.m` (×4) | `fitPiecewiseRatio`; per-station-per-bin means; profile-likelihood z\* |
| `regress_ratio/subroutines/makeGp15RegressRatio.m` | three-zone piecewise evaluation at all (station, depth) pairs |
| `collate_data/subroutines/collateStations.m` | table joins on depth key; unit conversions by variable name |

### Toolbox functions (`src/toolboxes/`)

| File | Algorithm |
|------|-----------|
| `regression/ols.m` | OLS/WLS; returns β, α, q (mean SE), ci = [varβ, varα, cov, n] |
| `regression/fitPiecewiseRatio.m` | profile-likelihood piecewise fit; per-station-per-bin means; β ≤ 0 hard constraint |
| `interp1dError/interp1dError.m` | piecewise-linear; propagates σ through basis weights; no extrapolation |
| `terraformer/ordinary_kriging.m` | Matérn-3/2 BLUP via Python CLI; FitHp=true, NRestarts=3 |
| `terraformer/fit_temporal_preprocessor.m` | depth-binned harmonic regression; returns anomalies + preproc struct |
| `terraformer/inverse_temporal_preprocess.m` | adds seasonal mean back at query coordinates |
| `utils/model/calcLayerThickness.m` | half-open midpoint-rule dz; bc=2 at surface; called from `doModelCalcs.m` and `calcError.m` |
| `utils/pipeline/computeFileHash.m` | MD5 hash via Java `MessageDigest`; file read with MATLAB `fread` in 8 kB chunks |
| `utils/pipeline/writeHash.m` | saves `.hash.json` sidecar after `save(...)` |
| `utils/pipeline/checkHash.m` | verifies sidecar hash; hard-errors on mismatch; warns if absent |

### Global constants

| Constant | Value | Set in | Meaning |
|----------|-------|--------|---------|
| `LAMBDA` | ln(2)/24.101 d⁻¹ | `setModelCoefficients.m` | ²³⁴Th decay |
| `P_SPLINE_GRADIENT` | 0.9 | `setModelCoefficients.m` | sparse/noisy profiles |
| `P_SPLINE_DENSITY` | 0.99 | `setModelCoefficients.m` | dense CTD profiles |
| `deltaDay` | 35 d | `configureUpwell.m` | temporal smoothing window |
| `nDaysWindow` | 35 d | `calcGp15Npp.m` | NPP trailing mean window |
| `BTMDEPTH` | 400 m | `configureInterp.m`, `collateStations.m`, `fitPiecewiseRatio.m` | integration ceiling |
| `FLIERDATA` | 30 µmol dpm⁻¹ | `qcGp15Poc.m` | POC:Th ratio flier cutoff |
| `cbpmFrac` | 0.30 | `calcStats.m` | CBPM algorithmic uncertainty floor |

---

## 23. Interpolation Dependency Map

Every interpolation and objective-mapping call through the pipeline, in execution order.

```
Raw CTD / bottle / satellite data
        │
        ▼
[1] PPZ depth        findppz — fluorescence median-filter, 10%-of-max threshold
        │
        ▼
[2] MLD depth        ordinary_kriging (T at 10 m) + csaps(P_SPLINE_DENSITY) + fnzeros
        │
        ▼
[3] EQP depth        csaps(P_SPLINE_GRADIENT, inv-var wts) on Th/U − 1 + fnzeros
        │
        ▼
[4] NPP              ordinary_kriging (lat, lon, time); 35-day trailing mean
        │
        ▼
[5] Profile interp   interp1dError — all obs columns to {MLD, PPZ, PPZ+100, EQP, 100m}
        │
        ▼
[6] Upwelling w      ordinary_kriging (lat, lon, depth, time) on 35-day smoothed wSpatAve
        │
        ▼
[7] Vertical ∂Th/∂z  finite differences + hard-zero outside [MLD, min(EQP,PPZ)]
        │
        ▼
[8] POC:Th ratio     fitPiecewiseRatio — per-station-per-bin means; profile-likelihood z*
        │
        ▼
[9] Flux integration cumsum (not interpolation — listed for completeness)
```

### [4] NPP kriging (`calcGp15Npp.m`)

| Attribute | Detail |
|-----------|--------|
| Input | Satellite CBPM NPP on a regular (lon, lat, time) grid |
| Method | `ordinary_kriging` (Matérn-3/2; FitHp=true; NRestarts=3; RandomState=7) |
| Temporal preprocessing | depth-binned harmonic regression; anomalies kriged |
| Query — instantaneous | one (lat, lon, t_sample) per station → `mmolC`, `stdev` |
| Query — 35-day window | 35 daily query points per station [t−34, t] in the same call → `mmolC_35d`, `stdev_35d` |
| Uncertainty | `√(stdev_35d² + (0.30 × mmolC_35d)²)` in EzRatio calculation |

### [6] Upwelling kriging (`interpUpwell.m`)

| Attribute | Detail |
|-----------|--------|
| Input | `w_<product>.nc` — 35-day temporally averaged, spatially smoothed w (m s⁻¹) |
| Domain | lon [129, 293]°E, lat [−60, 58]°N, depth ≤ 400 m; land/zero cells removed |
| Temporal preprocessing | `fit_temporal_preprocessor` over bins [0, 100, 200, 400, ∞] m |
| Method | `ordinary_kriging`; priors lx = 500 000 m, lz = 100 m, lt = 30 d |
| Post-processing | `inverse_temporal_preprocess` at each (station, depth, time) |
| Output | `gp15_w.<product>`: [stationNo, lon, lat, depth, time, w, wVar, wErr, wN] |

### [8] POC:Th ratio (`fitPiecewiseRatio.m`)

| Attribute | Detail |
|-----------|--------|
| Input | Raw (depth, ratio, σ, stationNo) for all valid obs in [0, 400 m] per region |
| Binning | n_bins = max(5, min(10, ⌊n_raw/5⌋)) equal-count rank-based bins |
| Data for fit | per-station-per-bin inv-variance weighted means (reduces serial correlation) |
| z\* search | unique raw depths in [MLD, 400 m]; profile likelihood; β ≤ 0 hard constraint |
| Covariance | `σ²_resid · (X'WX)⁻¹`; flat-line (β=0) uses 1-parameter form |
| CI on z\* | 95% concentrated profile-likelihood: RSS(z\*) ≤ RSS_opt · exp(3.841/n) |

### Dependency summary

| Consumer | Depends on |
|----------|------------|
| `doModelCalcs.m` | `gp15_inputs` (from `collateStations`) |
| `collateStations.m` | `gp15_obs` (interpolated), `gp15_w.*`, `gp15_grad`, `regressPocRatio`, `regressPnRatio` |
| `stationInterp.m` | `interp1dError`, QC'd `gp15_obs`, `gp15_mld`, `gp15_ppz`, `gp15_eqp` |
| `interpUpwell.m` | `ordinary_kriging`, 35-day-smoothed `w_<product>.nc` |
| `doVertGradCalc.m` | interpolated `gp15_obs`, `gp15_mld`, `gp15_ppz`, `gp15_eqp` |
| `calcGp15Npp.m` | `ordinary_kriging`, satellite NPP grid |
| `calcMlDepths.m` | `ordinary_kriging` (T at 10 m), CTD profiles |
| `calcEqpDepths.m` | QC'd `gp15_obs`, `gp15_mld`, `gp15_ppz` |
| `calcRegionalRegression*.m` | QC'd `gp15_obs`, `gp15_mld`, `gp15_stations`, `regions` |
| `makeGp15RegressRatio.m` | `ratioPocLsf`, `ratioPnLsf`, interpolated `gp15_obsFull` |

---

## 24. Canonical Method Reference

Each section below names the required implementation for a computational process used in multiple places. All call sites must conform. Do not substitute alternatives.

---

### 24.1 1D Profile Interpolation

**Function**: `interp1dError` (`src/toolboxes/interp1dError/interp1dError.m`)

```matlab
[vq, vqStd] = interp1dError(x, v, vStd, xq)
```

Piecewise-linear; propagates uncertainty through Lagrange weights. Returns NaN outside the knot range — never extrapolates. Input `x` must be NaN-free and sorted ascending. Do not use `interp1(..., 'linear')` — it silently extrapolates and discards uncertainty.

| Call site | Variables interpolated |
|---|---|
| `qcGp15Poc.m` | Th234_part LSF/SSF, POC LSF/SSF |
| `qcGp15Pn.m` | PN LSF/SSF |
| `stationInterp.m` | All observation columns |

---

### 24.2 Smoothing Spline

**Function**: `csaps` (MATLAB Curve Fitting Toolbox)

```matlab
splineFit = csaps(x, y, p, [], w);
yHat      = ppval(splineFit, xq);
```

`p` must be set explicitly — never pass `[]`. Both constants are defined in `setModelCoefficients.m` and loaded by every stage that uses splines.

| Constant | Value | Context |
|---|---|---|
| `P_SPLINE_GRADIENT` | 0.9 | Sparse/noisy profiles: Th234/U238 ratio (`calcEqpDepths.m` only) |
| `P_SPLINE_DENSITY` | 0.99 | Dense/precise CTD profiles: potential density, temperature |

Pass `w = 1 ./ sigma^2` when uncertainties are available. For CTD profiles (no per-point uncertainty), pass `w = ones(size(y))`.

| Call site | Variable | `p` |
|---|---|---|
| `calcEqpDepths.m` | Th234/U238 ratio | `P_SPLINE_GRADIENT`, inverse-variance weights (delta method) |
| `calcMlDepths.m` | Potential density, temperature | `P_SPLINE_DENSITY`, unit weights |

---

### 24.3 Ordinary Kriging

**Function**: `ordinary_kriging` (`src/toolboxes/terraformer/ordinary_kriging.m`)

`RandomState` is mandatory at every call site. `'FitHp', true` at all production call sites.

| Call site | Dimensions | H-scale (m) | V-scale (m) | T-scale (d) | RandomState |
|---|---|---|---|---|---|
| `calcGp15Npp.m` | lat, lon, time | 250 000 | — | 8.0 | 7 |
| `calcMlDepths.m` | depth only | — | 100 | — | 7 |
| `interpUpwell.m` | lat, lon, depth, time | 500 000 | 100 | 30.0 | 7 |

`std_pred` is interpolation uncertainty only; algorithmic/model uncertainty must be added separately.

---

### 24.4 Piecewise Ratio Regression

**Function**: `fitPiecewiseRatio` (`src/toolboxes/regression/fitPiecewiseRatio.m`)

Fits R(z) = α + β·clamp(z, MLD, z\*) with β ≤ 0. Profile likelihood over unique depths in [MLD, 400 m]. Per-station-per-bin inverse-variance weighted means reduce within-station serial correlation.

```matlab
[fittedAlpha, fittedBeta, fittedZStar, paramCov, zStarCI, flatLine, nBinMeans] = ...
    fitPiecewiseRatio(z_raw, r_raw, sigma_raw, station_raw, mld)
```

`paramCov = sSquared · (X'WX)⁻¹` where `sSquared = RSS/dof`. The ratio uncertainty at depth z passed downstream is the **1σ delta-method SE of the fitted mean**:

```
σ_R(z) = √( x(z)ᵀ · paramCov · x(z) )    x(z) = [1, clamp(z, MLD, z*)]
```

This is a confidence interval on the fitted mean R(z) — it does not include new-observation scatter. The 95% CI on z\* itself (`zStarCI`) is from the profile likelihood and is not propagated into flux errors.

`flatLine = true` (β = 0) is a legitimate degenerate case, not an error.

| Call site | Ratio | QC bounds |
|---|---|---|
| `calcRegionalRegressionPocLsf.m` | POC LSF : Th | 0.05–3.0 µmol dpm⁻¹ (lower = author-defined noise floor) |
| `calcRegionalRegressionPocSsf.m` | POC SSF : Th | 0–100 µmol dpm⁻¹ |
| `calcRegionalRegressionPnLsf.m` | PN LSF : Th | 0–100 µmol dpm⁻¹ |
| `calcRegionalRegressionPnSsf.m` | PN SSF : Th | 0–100 µmol dpm⁻¹ |

**Plot label**: the shaded band in `calcRegionalRegression*.m` diagnostic figures is labelled `$\pm 1\sigma$`, consistent with the 1σ propagation used throughout the pipeline.

---

### 24.5 Layer Thickness (dz)

**Function**: `calcLayerThickness` (`src/utils/model/calcLayerThickness.m`)

Half-open midpoint rule; `bc = 2` at the surface cell so the surface cell spans the full distance to the first sub-surface observation.

```matlab
dz = calcLayerThickness(z);
```

Do not reimplement inline. Called from `doModelCalcs.m` and `calcError.m`.

---

### 24.6 Finite-Difference Vertical Gradient

```
Forward  (i=1):     dydz = (y(2)   - y(1))   / (z(2)   - z(1))
Central  (1<i<n):   dydz = (y(i+1) - y(i-1)) / (z(i+1) - z(i-1))
Backward (i=n):     dydz = (y(n)   - y(n-1)) / (z(n)   - z(n-1))

σ(dydz) = sqrt( (σ_a/B)² + (σ_b/B)² )   where B = |z2 - z1|
```

Depth spacing `B` must use `abs(z2-z1)`.

| Call site | Quantity | Post-processing |
|---|---|---|
| `doVertGradCalc.m` | Total dissolved ²³⁴Th | hard-zero outside `[MLD, min(EQP,PPZ)]` |
| `calculateUpwelling.m` | u(z), v(z) for divergence | `cumtrapz` → w |

---

### 24.7 Cumulative ²³⁴Th Flux Integration

```
Flux_1D(i) = (U238(i) - Th234(i)) * dz(i) * lambda
Flux_2D(i) = Flux_1D(i) + w(i) * dTh234/dz(i) * dz(i)
F(z_k)     = sum_{i=1}^{k} Flux(i)          [cumsum]
```

Call site: `doModelCalcs.m`. Units after `collateStations` conversions:

| Quantity | Unit |
|---|---|
| U238, Th234 | dpm m⁻³ |
| dz | m |
| lambda | d⁻¹ |
| w | m d⁻¹ |
| dTh234/dz | dpm m⁻⁴ |
| Cumulative flux | dpm m⁻² d⁻¹ |

---

### 24.8 NPP Uncertainty

```matlab
sigma_NPP = sqrt( sigma_krig^2 + (0.30 * NPP)^2 )
```

`f_CBPM = 0.30` (Behrenfeld et al. 2001; Saba et al. 2011). Applied in `calcStats.m` using the 35-day trailing mean fields `mmolC_35d` / `stdev_35d`. The instantaneous fields `mmolC` / `stdev` are retained only for the diagnostic NPP plot.

---

### 24.9 Upwelling Variance Contribution

The **absolute-error** form is required:

```matlab
Var[w * grad * dz] = (grad * dz)^2 * sigma_w^2 + (w * dz)^2 * sigma_grad^2
```

Using relative errors would produce `0/0 = NaN` when w or grad is zero.

---

### 24.10 Flux Uncertainty Propagation

```
Var[F_1D(z_k)] = lambda^2 * sum_{i<=k} dz(i)^2 * (sigma_U(i)^2 + sigma_Th(i)^2)

Var[F_2D(z_k)] = Var[F_1D(z_k)]
               + sum_{i<=k} [ (grad(i)*dz(i))^2 * sigma_w(i)^2
                            + (w(i)*dz(i))^2   * sigma_grad(i)^2 ]

Var[F_POC]     = F_POC^2 * ( Var[F_Th] / F_Th^2 + Var[R_POC] / R_POC^2 )
```

Call site: `calcError.m`. Do not replace NaN variance contributions with zero before `cumsum` — NaN at a depth means the upwelling correction is undefined there.

---

### 24.11 Unit Conversions

Applied once, in `collateStations.m`, by variable name (never by column index):

| Constant | Value | Direction |
|---|---|---|
| `L2M` | 1000 | dpm L⁻¹ → dpm m⁻³ |
| `UMOL2MMOL` | 1000 | µmol → mmol |
| `SEC2DAY` | 86400 | m s⁻¹ → m d⁻¹ |

In `calculateUpwelling.m`, `deg2km(dxDeg, A)` is called with `A = 6378.137 × 1000` m to return metres rather than km. The mandatory `cos(lat)` factor is applied to all zonal distance branches.

---

### 24.12 Ratio Profile Reconstruction

Three-zone piecewise evaluation in `makeGp15RegressRatio.m`:

| Zone | Depth | Ratio | Uncertainty |
|---|---|---|---|
| Surface | z ≤ MLD | β·MLD + α | `sqrt(βVar·MLD² + αVar + 2·MLD·coVar)` |
| Active | MLD < z ≤ z\* | β·z + α | `sqrt(βVar·z² + αVar + 2·z·coVar)` |
| Deep | z > z\* | `ratioSubsurface` | `ratioSubsurfaceStd` (delta-method SE at z\*) |

`ratioSubsurface = alpha + beta * zstar` (continuity condition). `ratioSubsurfaceStd` is the delta-method SE at z\*, not SE of the mean.

---

### 24.13 Model Ensemble Statistics

Two ensemble spreads are computed in `calcError.m` after the per-product loop:

```matlab
% full ensemble — n=5 (ECCO, CGLO, FOAM, GLOR, ORAS)
th234FluxData      = [F_ecco, F_cglo, F_foam, F_glor, F_oras];
modelSpreadFull    = std(th234FluxData, 1, 2);     % no /sqrt(n)

% no-FOAM ensemble — n=4 (ECCO, CGLO, GLOR, ORAS)
th234FluxData_noFOAM   = [F_ecco, F_cglo, F_glor, F_oras];
modelSpreadNoFOAM      = std(th234FluxData_noFOAM, 1, 2);
```

GREP is excluded from both because it is the arithmetic mean of CGLO, FOAM, GLOR, and ORAS and is not independent. Do not divide by `sqrt(n)` — the models are deterministic reanalyses, not i.i.d. samples. σ_model is combined in quadrature with σ_2d only for the mean-ensemble products (`totalXxxFluxError_meanFull`, `_meanNoFOAM`); individual members carry σ_2d only.

---

### 24.14 Equilibrium-Point (EQP) Detection

`csaps(depthV, Th234/U238 − 1, P_SPLINE_GRADIENT, [], wts_eqp)` where `wts_eqp = 1 ./ sigmaRatio.^2` (delta-method ratio uncertainty). Zero crossings via `fnzeros`; minimum via `fminbnd(@(z) ppval(...), ...)`. Selection by weighted distance to MLD and PPZ (equal weights W = 0.5 each).

---

### 24.15 Regional Summary Statistics

```matlab
regionalMean = mean(x, 'omitnan');
regionalStd  = std(x, 0, 'omitnan');
regionalSE   = regionalStd / sqrt(sum(~isnan(x)));
```

Report as mean ± SD in regional tables; SE = SD/√n for the uncertainty of the regional mean.

---

### 24.16 Dependency Hashing

**Save side**: call `writeHash(matPath)` immediately after every `save(matPath, ...)`.

**Load side**: call `checkHash(matPath)` immediately before every `load(matPath, ...)` for pipeline intermediates.

No sidecar → warning (first run). Sidecar present, hash matches → silent pass. Sidecar present, hash mismatch → `error()` (re-run the upstream stage).

Covered files: `gp15_grad.mat`, `gp15_mld.mat`, `gp15_obs.mat`, `gp15_obsNoQc.mat`, `gp15_stations.mat`, `gp15_inputs.mat`, `gp15_flux*.mat`.

---

## 25. Tests

**Location**: `tests/`  
**Runner**: `tests/runTests.m` (calls `setupModel` then runs all four tests in sequence)

| File | What it tests | Data required |
|---|---|---|
| `testCalcError.m` | 1D and 2D flux variance formulae; hand-computed 3-depth synthetic station | None (self-contained) |
| `testFitPiecewiseRatio.m` | `fitPiecewiseRatio` recovers known β, α, z\* from synthetic linear profile | None (self-contained) |
| `testUncertMetrics.m` | T₁₀₀ and R₁₀₀ uncertainty formulae; hand-computed 2-layer synthetic station verifying the Th-234 variance form and ratio-cancellation derivation | None (self-contained) |
| `testGoldenFlux.m` | Full pipeline regression; diffs `gp15_flux100.mat` against a saved golden file; columns checked: `th234FluxCumul1d`, `th234FluxCumul2d_ecco`, `totalPocFluxError`, `pocFluxCumul2d_ecco` | Full pipeline run required; first run creates `tests/golden/gp15_flux100_golden.mat` |

Each test prints `pass` or `FAIL` per check and a summary line.

**Full reproducibility run**: `runGp15.m` in the repository root is the canonical script for a complete run-and-confirm cycle. It calls `gp15Model` (full pipeline), then `runTests` twice. The first `runTests` call generates `tests/golden/gp15_flux100_golden.mat` from the fresh output (the other three tests run and pass immediately). The second `runTests` call validates all four tests — including `testGoldenFlux` — against the newly created golden file. Both runs passing is the confirmation criterion (T4).

**Updating the golden file**: after any intentional model change, delete `tests/golden/gp15_flux100_golden.mat` and run `runGp15.m` from scratch to regenerate and re-confirm it. Do not run `runTests` alone to regenerate — the pipeline must be re-run first so the golden file reflects current output.

**Status (2026-06-05)**: all four tests pass. `gp15_flux100_golden.mat` has been committed as the permanent regression baseline. The pipeline is fully frozen.
