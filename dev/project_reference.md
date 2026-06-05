# GP15 Project Reference

**Project**: Kenyon and Davidson et al. (in prep.), GP15 Pacific GEOTRACES transect  
**Languages**: MATLAB (pipeline), Julia (figures), Python via `terraformer` CLI (kriging)  
**Pipeline status**: FULLY FROZEN. All methodology and submission-blocking issues resolved. All automated tests pass, including golden-file regression (T4 complete, 2026-06-05). No further code changes without explicit author discussion.  
**Last updated**: 2026-06-05 (Session 7 — T4 complete)

---

## 1. Repository Structure

```text
gp15_part_export_remin/
├── scripts/
│   ├── gp15Model.m                    ← MASTER SCRIPT
│   ├── model/                         ← Th-234 flux model (core science)
│   ├── profile_interpolation/         ← spatiotemporal kriging of upwelling
│   ├── quality_control/
│   │   ├── doQc.m                     ← QC entry point
│   │   ├── plotQcPoc.m                ← standalone POC QC diagnostic plots (horizontal error bars)
│   │   └── plotQcPn.m                 ← standalone PN QC diagnostic plots (horizontal error bars)
│   ├── read_data/                     ← raw data ingestion
│   ├── collate_data/                  ← assembles final model input table
│   ├── calculate_ppz/                 ← particle flux proxy depth
│   ├── calculate_npp/                 ← satellite NPP kriging
│   ├── calculate_mld/
│   │   ├── calcMld.m                  ← MLD entry point
│   │   ├── plotMld.m                  ← standalone MLD diagnostic plots
│   │   └── subroutines/sensitivitySplineDensity.m  ← JAK density-spline p sensitivity
│   ├── calculate_eqp/
│   │   ├── calcEqp.m                  ← EQP entry point
│   │   └── subroutine/plotEqp.m       ← standalone EQP diagnostic plots
│   ├── calculate_upwell/
│   │   ├── calcUpwell.m               ← upwelling entry point (MERCATOR loop + ECCO block)
│   │   └── subroutines/
│   │       ├── sensitivityDeltaDay.m  ← temporal window sensitivity; title shows spatial deg
│   │       ├── sensitivitySpatialAve.m ← spatial window sensitivity; legend + title in physical deg; movmean bypassed for window=1
│   │       ├── plotUpwell.m           ← meridional slice + per-station figures; titles in physical deg; markerSize=10
│   │       └── compareUpwellModels.m  ← CMEMS members + mean + GREP; titles in physical deg; grid from w_grep.nc
│   ├── calculate_gradient/
│   │   ├── calcVertGrad.m             ← gradient entry point
│   │   └── plotVertGrad.m             ← standalone gradient diagnostic plots
│   ├── regional_delineation/          ← biogeographic province assignment
│   ├── regress_ratio/                 ← POC:Th and PN:Th depth regression
│   └── plot_data/                     ← Julia plotting (manuscript + SI)
├── src/
│   ├── toolboxes/
│   │   ├── terraformer/               ← ordinary_kriging (Python CLI wrappers)
│   │   ├── interp1dError/             ← piecewise-linear interp with uncertainty
│   │   └── regression/               ← ols.m, fitPiecewiseRatio.m
│   └── utils/
│       ├── model/                     ← setupModel.m, calcLayerThickness.m
│       └── pipeline/                  ← computeFileHash.m, writeHash.m, checkHash.m
├── tests/                             ← assertion-script unit and regression tests
│   ├── runTests.m
│   ├── testCalcError.m
│   ├── testFitPiecewiseRatio.m
│   ├── testGoldenFlux.m
│   └── golden/                        ← gp15_flux100_golden.mat (created on first run)
├── data/
│   ├── data_raw/                      ← immutable raw inputs
│   ├── data_pro/                      ← processed intermediates
│   └── sim/                           ← simulation outputs
├── docs/                              ← manual.md (canonical reference)
└── dev/                               ← this file
```

**Pipeline execution order** (all called by `gp15Model.m`):

```text
readData → calcPpz → doQc → calcNpp → calcMld → calcEqp → calcUpwell
→ interpData → calcVertGrad → delineateRegions → regressRegionalRatio
→ collateData → runModel
```

Each stage can be run standalone. `setupModel.m` sets all paths via `mfilename('fullpath')`, clears the workspace, and calls `addpath(genpath(pwd))` — all subdirectories including `src/utils/pipeline/` are on the path automatically.

---

## 2. Core Methodology

### Th-234 flux model

**1D** (all stations):

$$F_\text{Th}(z) = \sum_i \lambda \cdot \left[{}^{238}\text{U}_i - {}^{234}\text{Th}_i\right] \cdot dz_i \qquad [\text{dpm m}^{-2}\text{ d}^{-1}]$$

**2D** (equatorial stations, $|\varphi| \leq 5°$):

$$F_{\text{Th,2d}}(z) = \sum_i \left[\lambda \cdot \left({}^{238}\text{U}_i - {}^{234}\text{Th}_i\right) + w_i \cdot \left(\frac{\partial\,{}^{234}\text{Th}}{\partial z}\right)_i\right] dz_i$$

**POC/PN flux**:

$$F_\text{POC} = F_\text{Th} \times R_\text{POC}(z,\,\text{region})$$

**Decay constant**: `LAMBDA` $= \ln(2)/24.101\ \text{d}^{-1}$ (`setModelCoefficients.m`). Never change.

**dz scheme** (`calcLayerThickness.m`): half-open cell at surface ($\text{bc}=2$), centered elsewhere. Called from `doModelCalcs.m` and `calcError.m` — do not reimplement inline.

**Total uncertainty** (`calcError.m`):

$$\sigma_\text{total} = \sqrt{\sigma^2_{\text{2d,ECCO}} + \sigma^2_\text{model}}$$

where $\sigma_\text{2d,ECCO}$ is the measurement-propagated 2D flux uncertainty and $\sigma_\text{model} = \text{std}(F_\text{2d}\text{ over 5 independent products})$ (no $\sqrt{n}$ divisor — models are deterministic reanalyses, not i.i.d. samples).

### Upwelling model products

Six products are processed. Five independent products enter the flux model ensemble; GREP is retained for diagnostic plotting only:

| Product | Full name | Institution | Role |
|---------|-----------|-------------|------|
| ECCO | ECCO v4 | JPL/MIT | Reference + ensemble |
| CGLO | C-GLORS05 | CMCC (Italy) | Ensemble |
| FOAM | GloSea5 | Met Office (UK) | Ensemble |
| GLOR | GLORYS2V4 | Mercator Ocean (France) | Ensemble |
| ORAS | ORAS5 | ECMWF | Ensemble |
| GREP | CMEMS Global Ensemble Reanalysis | CMEMS | Diagnostic only — arithmetic mean of CGLO/FOAM/GLOR/ORAS; not independent |

The reported 2D flux and its uncertainty are controlled by two variables at the top of `calcStats.m` (`fluxProduct`, `uncertProduct`; both default to `'meanFull'`). The five independent products (ECCO, CGLO, FOAM, GLOR, ORAS) enter the ensemble spread ($n=5$). GREP is interpolated alongside the others for diagnostic comparison in `compareModels.m` but is excluded from the ensemble uncertainty calculation.

`compareModels.m` produces two figures (100 m and PPZ stacked in each):
- `all_model_transect.pdf` — full ensemble: CMEMS members (gray, distinct line styles), CMEMS mean (thick black), GREP (thinner red), ECCO (dark blue), all-model mean of five independent products (thick dark green).
- `all_model_transect_noFOAM.pdf` — FOAM and GREP excluded: same styling over CGLO, GLOR, ORAS only; ensemble means recomputed without FOAM.

**FOAM outlier**: FOAM (GloSea5) is a strong outlier in equatorial upwelling velocity relative to the other four products. Cause: NEMOVAR 3D-Var does not directly constrain near-equatorial velocities (ageostrophic flow poorly corrected by assimilation increments); NEMO has documented Pacific EUC biases; z\*-coordinate free-surface corrections accumulate in the vertical integral (Storto et al. 2019). FOAM is the primary driver of $\sigma_\text{model}$, making the reported total uncertainty conservative. GREP is excluded from the sensitivity figure because it incorporates FOAM.

**FOAM sensitivity**: `collateStations.m` stores `w_meanNoFOAM` (n=4, FOAM excluded) and `wErr_meanNoFOAM`. `calcError.m` computes `totalXxxFluxError_meanNoFOAM` using the 4-product ensemble spread. Set `fluxProduct = uncertProduct = 'meanNoFOAM'` in `calcStats.m` to run the FOAM-excluded sensitivity; report alongside `all_model_transect_noFOAM.pdf`.

ECCO processing in `calcUpwell.m`: ECCO provides `wo` directly (not u/v), so it bypasses `calculateUpwelling.m`. After the main MERCATOR loop, a dedicated ECCO block loads `ecco_wo.mat`, runs `spatialTimeAverage` with `spaceResolution = 1.0` deg (native 1° grid) and `timeResolution = 10` d (ECCO `Wave` is a 10-day average field, confirmed from `ecco_wo_quart*.cdf`; `timeAve = floor(35/10) = 3`, giving a 30-day moving average), and saves to `w_ecco.nc`. Grid coordinates for `compareUpwellModels.m` are loaded from `w_grep.nc` rather than `u_mercator`, which is overwritten by the ECCO block before `compareUpwellModels` runs. This ensures all six products enter `interpUpwell.m` with identically processed `wSpatAve` training data.

Save naming: all six products → `w_<product>.nc` (netCDF4) via `saveUpwell.m`.

### w-values architecture (implemented)

**Design**: all six per-product kriged w tables are joined as individual columns in `gp15_inputs` by `collateStations.m` (`w_ecco`, `w_cglo`, `w_foam`, `w_glor`, `w_oras`, `w_grep`). After the per-product join, two ensemble mean columns and their propagated kriging uncertainties are computed:

```matlab
% wErr_mean = sqrt(Σ wErr_k²) / n  — uncertainty of mean of n independent quantities
gp15_inputs.w_meanFull      = mean([w_ecco w_cglo w_foam w_glor w_oras], 2, 'omitnan');   % n=5
gp15_inputs.wErr_meanFull   = sqrt(wErr_ecco² + wErr_cglo² + wErr_foam² + wErr_glor² + wErr_oras²) / 5;
gp15_inputs.w_meanNoFOAM    = mean([w_ecco w_cglo w_glor w_oras], 2, 'omitnan');           % n=4
gp15_inputs.wErr_meanNoFOAM = sqrt(wErr_ecco² + wErr_cglo² + wErr_glor² + wErr_oras²) / 4;
```

**Total-error columns** (`calcError.m`): `totalXxxFluxError_ecco` (alias for backward compat), `_meanFull` (σ_2d,meanFull + σ_model,full; n=5), `_meanNoFOAM` (σ_2d,meanNoFOAM + σ_model,noFOAM; n=4). σ_model is attached only to the mean-ensemble products; individual members carry σ_2d only.

**Configurable reporting** (`calcStats.m`): `fluxProduct` and `uncertProduct` variables (both default `'meanFull'`) select which flux and uncertainty columns are written into `depthData` and the manuscript tables. Example: `fluxProduct='ecco', uncertProduct='meanFull'` reports ECCO flux with full-ensemble total uncertainty. Cross-references: `collateData.m` (inline description), `interpData.m` (FOAM outlier note), `compareModels.m` (sensitivity figure).

### Upwelling interpolation (`interpUpwell.m`)

Single spatiotemporal ordinary kriging call via `terraformer`:

1. `fit_temporal_preprocessor` — depth-binned harmonic regression over bins $[0,\,100,\,200,\,400,\,\infty]$ m removes seasonal trend from training $w$ field
2. `ordinary_kriging` — Matérn-3/2 kernel; priors `lxInit = 500 000` m horizontal, `lt = 30` d temporal (`FitHp=true`, `NRestarts=3`, `RandomState=7`)
3. `inverse_temporal_preprocess` — adds back per-depth-bin seasonal mean at query (station, depth, time)

### NPP — 35-day trailing mean

`calcGp15Npp.m` queries the ordinary kriging at 35 daily time steps per station ([t − 34, t]) in a single combined kriging call. The per-station mean and RMS kriging SD are stored as `mmolC_35d` / `stdev_35d`. `calcStats.m` uses these for EzRatio so that numerator (²³⁴Th flux, integrating ~35 d) and denominator (NPP) share the same temporal window.

### Ratio regression (`regressRegionalRatio`)

Profile-likelihood piecewise fit via `fitPiecewiseRatio`. Fits $R(z) = \alpha + \beta \cdot \text{clamp}(z,\,\text{MLD},\,z^*)$ with $\beta \leq 0$. The transition depth $z^*$ is found by sweeping the profile likelihood over unique raw measurement depths in [MLD, 400 m]. Within-station serial correlation is reduced by using per-station inverse-variance weighted means in equal-count depth bins as regression data. Output table includes `zstarCI_lo`, `zstarCI_hi` (95% profile-likelihood CI on $z^*$).

---

## 3. Development Rules

1. **Read before editing.** The pipeline shares workspace variables across scripts — a variable name in one script must match exactly what upstream scripts produce.

2. **Comment-out, then fix.** When correcting a bug, comment the original line with a short explanation, then put the fixed line immediately after. Do not delete original lines.

3. **Unit discipline.** Three unit systems coexist and are only converted in `collateStations.m`:
   - Radionuclides: stored as dpm/L in `gp15_obs`; converted to dpm/m³ (`× L2M = 1000`)
   - Particulates: stored as µmol in `gp15_obs`; converted to mmol (`÷ UMOL2MMOL = 1000`)
   - Upwelling: stored as m/s in `gp15_w`; converted to m/d (`× SEC2DAY`)
   - Fluxes: reported as dpm m⁻² d⁻¹ (Th-234) and mmol m⁻² d⁻¹ (POC, PN)

4. **Depth convention**: all depths are positive downward. `BTMDEPTH = 400 m` is the integration ceiling used consistently across `collateStations.m`, `configureInterp.m`, and the ratio regression.

5. **Longitude convention**: station longitudes are in $[-180,\,180]$. Model grid (ECCO/MERCATOR) uses $[0,\,360]$. The $+360$ offset is applied in `interpUpwell.m` query construction; removed again when saving (`query_lon − 360`).

6. **Dynamic field access**: `calcError.m`, `doModelCalcs.m`, `calcParticulateFlux.m`, and `collateStations.m` loop over the six model products using MATLAB dynamic struct/table field access: `s.(['field_' dataProduct]) = ...`. When adding variables that need per-product treatment, follow this pattern and add the new product name to the cell array at the top of each loop.

7. **No new features without discussion.** The codebase is publication-ready. Changes should be bug fixes or targeted methodology improvements.

8. **String style**: use single-quoted `char` arrays (`'like this'`) throughout. Dynamic field access uses bracket concatenation: `(['prefix_' name])`.

9. **Dependency hashing**: after every `save(path, ...)` for a pipeline intermediate, call `writeHash(path)`. Before every `load(path, ...)` for a pipeline intermediate, call `checkHash(path)`. See `docs/manual.md` §24.16.

---

## 4. Known Fragile Areas

| File | Risk | Detail |
|------|------|--------|
| `model/subroutines/makeDepthArrays.m` | Medium | Tolerance-based logical indexing with hard asserts. Assert fires if PPZ or PPZ+100 depths were not inserted during `stationInterp.m`. If an assert triggers, check `interpData.m` ran successfully before `runModel.m`. |
| `model/subroutines/calcError.m` | Low | Upwelling variance uses absolute-error form `(g·dz)²·σ_w² + (w·dz)²·σ_g²`. Both `yw`/`nw` loading loops and the error-calculation loop must all cover `{ecco, cglo, foam, glor, oras, grep, meanFull, meanNoFOAM}` (8 products). Correction uncertainty: `sqrt(max(0, sigma_2D^2 - sigma_1D^2))` — quadrature difference of variances, not difference of SDs. Ensemble std (5-product and 4-product no-FOAM) computed after the loop. Named total-error columns: `totalXxxFluxError_ecco` (alias), `_meanFull`, `_meanNoFOAM`. No `/ sqrt(n)`. NaN propagates through `cumsum` — do not pre-zero. |
| `profile_interpolation/subroutines/interpUpwell.m` | Medium | If `terraformer` Python environment is not installed, all kriging calls fail. Check `python -m terraformer._cli.ok --help` before running. Hyperparameters are priors — optimizer takes over. All six products (including ECCO) load from `w_<product>.nc`; run `calcUpwell.m` first to generate these files. |
| `model/subroutines/calcStats.m` | Medium | `fluxProduct` and `uncertProduct` at top of file control what is written to `depthData` (default both `'meanFull'`). Equatorial 2D correction is latitude-based: `abs(statLat) > 5`; all 8 occurrences must remain consistent. `uncertT100` uses Th-234 flux variances (not POC) — ratio R cancels exactly in the delta method when both depths share the same deep-zone R; reverting to POC variances produces a negative argument clamped to zero. |
| `regress_ratio/subroutines/calcRegionalRegression*.m` (×4) | Low | Calls `fitPiecewiseRatio`. `flatLine=true` (β=0) is a legitimate degenerate case — not an error. `uncertPocTh234RatioRegress` is the **1σ delta-method SE of the fitted mean** at each depth; used as 1σ in `calcError.m`. The plot legend labels the shaded band `$\pm 1\sigma$`, consistent with the 1σ pipeline propagation. |
| `src/toolboxes/regression/fitPiecewiseRatio.m` | Low | `BTMDEPTH = 400` must equal the pipeline constant. Column j of `nppWindowMat = reshape(VhatWindowMean, nDaysWindow, nStations)` must align with station j — do not transpose. |
| `collate_data/subroutines/collateStations.m` | Low | Six per-station uniqueness asserts catch duplicate depths before any `join()`. Depth-alignment assert before positional upwelling zeroing — do not remove. Zero-gradient zeroing loop covers all 8 products including `meanFull` and `meanNoFOAM`. Unit conversions use named variable access. Mean columns (`w_meanFull`, `wErr_meanFull`, `w_meanNoFOAM`, `wErr_meanNoFOAM`) are computed after the unit conversion loop. |
| `quality_control/subroutines/subsubroutines/qcGp15Poc.m` | Low | Flier index captured before NaN-ification. The `flierIdxLarge` / `flierIdxSmall` capture pattern must be preserved if this section is refactored. |
| `calculate_gradient/subroutines/doVertGradCalc.m` | Low | Raw FD gradient and measurement-propagated SD hard-zeroed outside `[mldVal, zBot]`. `mldVal` looked up by station number (`gp15_mld.('Station No') == sn`) — do not revert to positional `iStat` indexing. |

### Numerical invariants to preserve

| Constant | Value | Files | Impact |
|----------|-------|-------|--------|
| `LAMBDA` | $\ln(2)/24.101\ \text{d}^{-1}$ | `setModelCoefficients.m` | ²³⁴Th decay in flux calculations |
| `BTMDEPTH` | 400 m | `configureInterp.m`, `collateStations.m` | Integration ceiling |
| `FLIERDATA` | 30 µmol dpm⁻¹ | `qcGp15Poc.m` | POC:Th ratio flier cutoff |
| `DepthBinEdges` | $[0,\,100,\,200,\,400,\,\infty]$ m | `interpUpwell.m` | Vertical bins for upwelling preprocessing |
| `cbpmFrac` | 0.30 | `calcStats.m` | CBPM algorithmic uncertainty floor |
| Full ensemble | $n=5$: ECCO, CGLO, FOAM, GLOR, ORAS | `calcError.m` | σ_model,full — inter-model std, no $/ \sqrt{n}$ |
| No-FOAM ensemble | $n=4$: ECCO, CGLO, GLOR, ORAS | `calcError.m` | σ_model,noFOAM — same formula, FOAM excluded |
| `fluxProduct` default | `'meanFull'` | `calcStats.m` (top of file) | Reported 2D flux magnitude — change to swap products |
| `uncertProduct` default | `'meanFull'` | `calcStats.m` (top of file) | Reported total uncertainty — change to swap products |
| `dz` formula | Half-open midpoint rule; bc=2 at surface | `calcLayerThickness.m` | Must not be reimplemented inline |
| Longitude offset | $+360$ / $-360$ | `interpUpwell.m` | Station [−180, 180] ↔ model [0, 360] |
| Equatorial band | $|\varphi| \leq 5°$ | `calcStats.m` (8 occurrences) | 2D upwelling correction reporting filter |

---

## 5. Key Files to Verify After Any Model Change

| File | What to check |
|------|---------------|
| `scripts/model/subroutines/calcError.m` | Both `yw`/`nw` loading loops (lines 15, 26) and the error-calculation loop (line 68) all cover `{ecco, cglo, foam, glor, oras, grep, meanFull, meanNoFOAM}` (8 products — must stay consistent); `uncertUpwellCorrect = sqrt(max(0, sigma_2D^2 - sigma_1D^2))` for Th-234/POC/PN (quadrature variance difference, not SD difference); GREP excluded from ensemble std arrays; full-ensemble std uses $n=5$, no-FOAM std uses $n=4$; no `/ sqrt(n)`; `upwellError` uses absolute-error form with `dz`; NaN propagates through `cumsum`; named total-error columns `_ecco`, `_meanFull`, `_meanNoFOAM` present |
| `scripts/model/subroutines/makeDepthArrays.m` | Uses `find(abs(...) <= tol)` and per-station loop |
| `src/toolboxes/regression/ols.m` | `ci = [Ci.varBeta, Ci.varAlpha, Ci.covAlphaBeta, n]` |
| `scripts/model/subroutines/calcStats.m` | `fluxProduct` / `uncertProduct` at top of file control reported magnitude and uncertainty (default both `'meanFull'`); EzRatio denominator uses `mmolC_35d`; `uncertT100` uses Th-234 variance form (ratio cancels exactly — do not revert to POC variances); `uncertR100` uses `sqrt(max(σ²_Th(PPZ+100) − σ²_Th(PPZ), 0))`; all 8 equatorial guards use `abs(statLat) > 5` |
| `scripts/calculate_npp/subroutines/calcGp15Npp.m` | `GpTimeScale = 8.0` (MODIS period); `nppWindowMat = reshape(VhatWindowMean, nDaysWindow, nStations)` — column j = station j |
| `scripts/regress_ratio/subroutines/calcRegionalRegression*.m` | Calls `fitPiecewiseRatio`; no region-2 override; tables have `zstarCI_lo/hi` columns |
| `scripts/profile_interpolation/subroutines/interpUpwell.m` | `lxInit = 500e3`, `lt = 30` hardcoded; `RandomState = 7`; all products load from `w_<product>.nc` |
| `scripts/profile_interpolation/subroutines/configureInterp.m` | `pacificLon = [129, 280]`; `pacificLat = [-60, 60]`; `pacificDepth = [0, 400]` — see inline comments for rationale |
| `src/utils/model/setupModel.m` | Uses `fileparts(mfilename('fullpath'))`, not a hardcoded path |
| `scripts/calculate_gradient/subroutines/doVertGradCalc.m` | Raw FD hard-zeroed outside `[mldVal, zBot]`; uncertainty is measurement-propagated SD; no spline |
| `src/utils/pipeline/checkHash.m`, `writeHash.m` | Called before/after every pipeline `load`/`save` |
| `tests/runTests.m` | All four tests pass after any model change (`testCalcError`, `testFitPiecewiseRatio`, `testUncertMetrics`, `testGoldenFlux`) |
