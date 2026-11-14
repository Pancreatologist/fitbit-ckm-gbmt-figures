# fitbit-ckm-gbmt-figures

R code that draws the figures for a multi-trajectory analysis of three
device-measured indicators (daily steps, mean heart rate, sleep duration) and
incident cardiometabolic-kidney outcomes.

The repository contains **code only**. No individual-level data, no `.rds`
objects, no imputation objects and no intermediate data files are included.
Every number that appears in a figure is read from a CSV that the script is
given; nothing is typed into the drawing code.

## Contents

| File | Draws |
| --- | --- |
| `R/11a_figure_sample_flow.R` | Figure 1: sample-flow diagram (population boxes, exclusion boxes, two outcome-cohort boxes) |
| `R/11b_figure_ckm_cuminc.R` | Cumulative incidence by trajectory group with a number-at-risk table (Aalen–Johansen estimator, death as a competing event) |
| `R/13_final_figures.R` | Final versions of the group-trajectory figure, the main forest plot, the sensitivity forest plot and the cumulative-incidence figure; reads stored CSVs only |
| `R/_traj_group_contract.R` | Shared definition of the four trajectory groups (codes, reference group, factor levels); sourced by `11b` |
| `R/07d_summary_v33.R` | Included unchanged. Its figure section is the origin of the two forest-plot data tables that `13_final_figures.R` reads. See *Notes* |

## Requirements

* R ≥ 4.2 (tested with R 4.5.3)
* Packages: `readr`, `ggplot2`, `jsonlite`, `survival`, `ragg`
* The **Arial** font must be installed; the PDF device is `cairo_pdf`, so R must
  be built with cairo support
* `ragg::agg_tiff` writes the 600 dpi TIFF versions

```r
install.packages(c("readr", "ggplot2", "jsonlite", "survival", "ragg"))
```

## Paths and configuration

All paths in the scripts are relative to one root directory. The root defaults
to the current working directory and can be changed with three environment
variables, so no absolute path is hard-coded anywhere:

| Variable | Default | Meaning |
| --- | --- | --- |
| `FIG_ROOT` | `.` | Root that contains `output/` and `intermediate/` |
| `FIG_INPUT_RDS` | `<FIG_ROOT>/gate1v3_fitbit_ckm_w8q_min3.rds` | Input `.rds` used by `11b` only |
| `FIG_SCRIPT_DIR` | `R` | Directory that holds `_traj_group_contract.R` |

Expected layout under `FIG_ROOT`:

```
<FIG_ROOT>/
  output/tables/                      # input CSVs (see schemas below)
  output/figures/                     # figures written here
  output/figures/final/               # final copies written here
  output/logs/                        # provenance JSON and label checks
  output/reports/                     # one markdown report per script
  intermediate/traj_groups_w8q.csv    # used by 11b only
  gate1v3_fitbit_ckm_w8q_min3.rds     # used by 11b only
```

## Input schemas

Column names must match exactly. Types are what `readr` infers.

### `output/tables/07_table_supp_sample_flow.csv` — Figure 1

| Column | Type | Meaning |
| --- | --- | --- |
| `stage` | character | Grouping label of the filtering stage |
| `step` | character | Step key; the script looks rows up by this key |
| `rule` | character | Description of the rule applied at the step |
| `n_before` | integer | Records before the step (`NA` on the first row) |
| `n_excluded` | integer | Records removed by the step (`NA` on the first row) |
| `n_after` | integer | Records remaining after the step |
| `events` | integer | Event count where applicable, otherwise `NA` |
| `source` | character | Short provenance tag |
| `note` | character | Free text |

The file must have exactly 12 rows and must contain these `step` keys:
`fitbit_any`, `valid_periods`, `has_ehr`, `prior_ehr`, `alive_adult`,
`landmark_before_cdr_end`, `has_fu_ehr`, `main_cohort`, `CKM_any_atrisk`,
`Cardiorenal_atrisk`.

```
stage,step,rule,n_before,n_excluded,n_after,events,source,note
Stage 0,fitbit_any,Initial population,NA,NA,100000,NA,step0,
Stage 1,valid_periods,At least three valid quarters,100000,20000,80000,NA,step1,
```

The script stops with an error if the chain does not balance
(`n_before - n_excluded == n_after`), if the rows are not continuous, or if any
count that would be printed is below 20.

### `output/tables/09b_g4_fitted_trajectories_long.csv` — trajectory lines

| Column | Type | Meaning |
| --- | --- | --- |
| `group` | integer/character | Trajectory group code (`1`–`4`) |
| `group_label` | character | Group name as stored upstream |
| `group_prior_pct` | numeric | Group size in per cent |
| `indicator` | character | `steps_mean`, `hr_mean` or `sleep_mean` |
| `indicator_label` | character | Human-readable indicator name |
| `unit` | character | Unit of the indicator |
| `period` | integer | Quarter index, `0`–`7` |
| `fitted_z` | numeric | Fitted value on the standardised scale |
| `fitted_value` | numeric | Fitted value on the original scale |
| `param_mean` | numeric | Mean used for the inverse transform |
| `param_sd` | numeric | SD used for the inverse transform |

96 rows (4 groups × 3 indicators × 8 quarters). `13_final_figures.R` reads
`group`, `indicator`, `period` and `fitted_value`.

```
group,group_label,group_prior_pct,indicator,indicator_label,unit,period,fitted_z,fitted_value,param_mean,param_sd
1,Group 1,25.00,steps_mean,Steps per day,steps/day,0,0.700,10500.000,8000.000,3500.000
```

### `output/tables/09c_observed_by_group_quarter.csv` — trajectory points

| Column | Type | Meaning |
| --- | --- | --- |
| `traj_group` | character | Trajectory group code (`"1"`–`"4"`) |
| `group_label` | character | Group name as stored upstream |
| `period` | integer | Quarter index, `0`–`7` |
| `indicator` | character | `steps_mean`, `hr_mean` or `sleep_mean` |
| `unit` | character | Unit of the indicator |
| `mean` | numeric | Observed mean in that group and quarter |
| `n_nonmissing` | integer | Non-missing observations behind the mean |
| `n_group_period` | integer | Individuals in the group; constant within a group |

96 rows. The script asserts that the plotted points equal the `mean` column in
all 96 group–quarter–indicator cells (maximum absolute difference 0).

```
traj_group,group_label,period,indicator,unit,mean,n_nonmissing,n_group_period
4,Group 4,0,steps_mean,steps/day,7800.00,2000,2000
```

### `output/tables/07d_table_main_hr_v33.csv` — main forest plot

| Column | Type | Meaning |
| --- | --- | --- |
| `outcome` | character | Outcome key; the figure uses `CKM_any` |
| `model` | character | `M1`, `M2` or `M3` |
| `group` | character | Trajectory group code (`"1"`–`"3"`; group 4 is the reference) |
| `label` | character | Group label as stored upstream |
| `hr` | numeric | Hazard ratio |
| `ci_low`, `ci_high` | numeric | 95% confidence limits |
| `p` | numeric | p value |

Additional columns may be present; they are ignored. The script takes the
`CKM_any` rows and requires exactly 9 (3 groups × 3 models).

```
outcome,model,group,label,hr,ci_low,ci_high,p
CKM_any,M1,1,Group 1,1.18,0.89,1.57,0.254
```

### `output/tables/07d_table_supp_sensitivity_v33.csv` — sensitivity forest plot

| Column | Type | Meaning |
| --- | --- | --- |
| `outcome` | character | Outcome key; must be `CKM_any` |
| `analysis` | character | `S1`, `S2`, `S3` or `S4` |
| `analysis_label` | character | Long description of the analysis |
| `n`, `events` | integer | Size and event count of the analysis |
| `hr_g1`, `lo_g1`, `hi_g1`, `p_g1` | numeric | Group 1 estimate and interval |
| `hr_g2`, `lo_g2`, `hi_g2`, `p_g2` | numeric | Group 2 estimate and interval |
| `hr_g3`, `lo_g3`, `hi_g3`, `p_g3` | numeric | Group 3 estimate and interval |

Four rows. The figure adds the primary model as a fifth row, taken from the
`M2` block of `07d_table_main_hr_v33.csv`, giving 15 plotted estimates.

```
outcome,analysis,analysis_label,n,events,hr_g1,lo_g1,hi_g1,p_g1,hr_g2,lo_g2,hi_g2,p_g2,hr_g3,lo_g3,hi_g3,p_g3
CKM_any,S1,Example analysis,4200,340,1.11,0.81,1.52,0.517,0.94,0.69,1.28,0.678,1.58,1.16,2.15,0.003
```

### `output/tables/11b_ckm_cuminc_curves.csv` — cumulative-incidence curves

| Column | Type | Meaning |
| --- | --- | --- |
| `group` | character | Trajectory group code |
| `time` | numeric | Years since the time origin |
| `cuminc` | numeric | Cumulative incidence (proportion, not per cent) |
| `lower`, `upper` | numeric | 95% confidence limits |

Step-function values, one row per event time per group, starting at time 0 with
incidence 0. Values must be non-decreasing within a group. `11b` writes this
file; `13_final_figures.R` reads it back.

```
group,time,cuminc,lower,upper
4,0,0,0,0
4,0.25,0.0021,0.0010,0.0040
```

### `output/tables/11b_ckm_risk_table.csv` — numbers at risk

| Column | Type | Meaning |
| --- | --- | --- |
| `group` | character | Trajectory group code |
| `time` | integer | Year, `0`–`5` |
| `n_at_risk` | integer | Individuals still under observation at that year |

Exactly 24 rows (4 groups × 6 time points).

```
group,time,n_at_risk
4,0,1100
```

### `intermediate/traj_groups_w8q.csv` — group assignment (11b only)

| Column | Type | Meaning |
| --- | --- | --- |
| `person_id` | character | Individual identifier |
| `traj_group` | integer | Group code, `1`–`4`; `4` is the reference group |
| `max_post_prob` | numeric | Posterior probability of the assigned group |

```
person_id,traj_group,max_post_prob
P0001,3,0.912
```

### Input `.rds` (11b only)

Not distributed with this repository. A serialised list with two data frames:

* `analytic` — one row per identifier. Columns read: `person_id`,
  `main_cohort` (logical), `CKM_any_prevalent` (0/1), `CKM_any_time_yr`
  (numeric), `landmark` (Date), `death_after_lm` (0/1), `death_date` (Date),
  `fu_end_ehr` (Date), `fu_end_admin` (Date).
* `status_long` — long format. Columns read: `person_id`, `outcome`
  (character), `censoring` (character), `event_date` (Date), `time_yr`
  (numeric), `event` (0/1).

The script asserts that the reconstructed follow-up time equals
`CKM_any_time_yr` for every row and that the event count matches its anchor.

## Running

Run from the repository root, with `FIG_ROOT` pointing at the directory that
holds the inputs:

```bash
export FIG_ROOT=/path/to/workdir

Rscript R/11a_figure_sample_flow.R
Rscript R/11b_figure_ckm_cuminc.R
Rscript R/13_final_figures.R
```

PowerShell:

```powershell
$env:FIG_ROOT = "D:\work"
& "C:\Program Files\R\R-4.5.3\bin\x64\Rscript.exe" R\13_final_figures.R
```

`13_final_figures.R` accepts an optional subset of figure names, so a single
figure can be redrawn in a few seconds:

```bash
Rscript R/13_final_figures.R Figure2_trajectory_G4
Rscript R/13_final_figures.R Figure4_forest_main FigureS3_forest_sensitivity
Rscript R/13_final_figures.R Figure3_ckm_cumulative_incidence
```

With no argument all four figures are written. The argument is matched against
the output file stem.

`11b` needs the `.rds`; point `FIG_INPUT_RDS` at a local copy:

```bash
export FIG_ROOT=/path/to/workdir
export FIG_INPUT_RDS=/path/to/gate1v3_fitbit_ckm_w8q_min3.rds
Rscript R/11b_figure_ckm_cuminc.R
```

## Outputs

| Script | Files | Size |
| --- | --- | --- |
| `11a` | `output/figures/11a_fig_sample_flow.pdf`, `.tiff`; copies to `output/figures/final/Figure1_sample_flow.*` | 174 × 240 mm, 600 dpi |
| `11b` | `output/figures/11b_fig_ckm_cuminc.pdf`, `.tiff` | 174 × 150 mm, 600 dpi |
| `13` | `output/figures/final/Figure2_trajectory_G4.pdf`, `.tiff` | 174 × 122 mm, 600 dpi |
| `13` | `output/figures/final/Figure3_ckm_cumulative_incidence.pdf`, `.tiff` | 174 × 150 mm, 600 dpi |
| `13` | `output/figures/final/Figure4_forest_main.pdf`, `.tiff` | 152 × 116 mm, 600 dpi |
| `13` | `output/figures/final/FigureS3_forest_sensitivity.pdf`, `.tiff` | 172 × 128 mm, 600 dpi |

Side products:

* `output/tables/11a_sample_flow_boxes.csv` — the box table that the flow
  diagram is drawn from
* `output/tables/11b_ckm_cuminc_curves.csv`,
  `output/tables/11b_ckm_cuminc_1_3_5y.csv`,
  `output/tables/11b_ckm_risk_table.csv`,
  `output/tables/11b_ckm_event_status.csv`
* `output/logs/11a_provenance.json`, `output/logs/11b_provenance.json`,
  `output/logs/13_provenance.json`, `output/logs/13_figure_labels.csv`
* `output/reports/step11a_sample_flow.md`, `output/reports/step11b_ckm_cuminc.md`,
  `output/reports/step13_final_figures.md`

All PDFs use the Arial family with embedded fonts. Figures carry no title, no
caption and no generator credit inside the image.

## Notes

* The scripts are drawing code. They do not fit models and do not recompute
  estimates: `13_final_figures.R` reads hazard ratios and cumulative-incidence
  values from stored CSVs, and `11b` writes its curves to CSV before drawing
  them, so the plotted values can be checked against the file.
* `R/07d_summary_v33.R` is included unchanged, except that no path edits were
  needed because it already uses relative paths. Its forest-plot section is
  coupled to the table-building section of the same script (it reads thirteen
  intermediate CSVs and also writes tables and a report), so it cannot be split
  without rewriting it. It is kept as the provenance of
  `07d_table_main_hr_v33.csv` and `07d_table_supp_sensitivity_v33.csv`, which
  are the direct inputs of the final forest plots. Running it end to end needs
  the upstream CSVs that are not part of this repository; the final figures
  themselves are reproduced by `13_final_figures.R` from the two tables above.
* Group colours are bound to the group code, not to the display order:
  group 4 (reference) `#CC79A7`, group 1 `#D55E00`, group 2 `#0072B2`,
  group 3 `#009E73`. Legends and the risk table are shown in the order
  G1 → G2 → G3 → G4.
* Every script stops with an error instead of drawing a figure whose input
  violates the expected row counts, the expected group sizes or the expected
  chain of counts. Those checks encode the published numbers; if you reuse the
  code with different data, edit the `EXPECTED_*` constants at the top of each
  script.
