# =============================================================================
# Step 11b: cumulative incidence of CKM disease by trajectory group (Figure B)
#   Population : CKM_any at-risk set of the main cohort (n = 4,648, 408 events).
#   Estimator  : Aalen-Johansen (survival::survfit on a multi-state Surv object),
#                so that death is a competing event and not a censoring event.
#   Time       : years from the landmark. For each person
#                  t = min(event date, death date, last EHR visit, CDR cut-off)
#                where the last two reproduce the existing censoring rule. The
#                reconstructed t is asserted to equal the stored
#                CKM_any_time_yr for every person (max absolute difference 0).
#   Status     : 1 = incident CKM disease, 2 = death before any CKM event and
#                before the censoring date (competing event), 0 = censored.
#                A death that falls after the last EHR visit is not observed as
#                a competing event, because the existing rule censors at the
#                last EHR visit first; such deaths are counted and reported.
#   95% CI     : from the Aalen-Johansen variance estimate of survfit
#                (conf.type = "log-log").
#   Every number is written to CSV by this script and the figure is drawn from
#   that CSV. No model is fitted, no existing table or figure is changed.
#   No console output (project rule).
# =============================================================================
suppressPackageStartupMessages({
  library(survival)
  library(readr)
  library(jsonlite)
})
SCRIPT_DIR <- Sys.getenv("FIG_SCRIPT_DIR", unset = "R")
source(file.path(SCRIPT_DIR, "_traj_group_contract.R"))

# All paths are relative to a single root. Override the root with the environment
# variable FIG_ROOT, or edit the default below; the layout under the root is
# output/tables (inputs), output/figures, output/logs and output/reports.
# The input RDS is not distributed with this repository; point FIG_INPUT_RDS at
# a local copy.
ROOT <- Sys.getenv("FIG_ROOT", unset = ".")
INPUT_RDS <- Sys.getenv("FIG_INPUT_RDS",
                        unset = file.path(ROOT, "gate1v3_fitbit_ckm_w8q_min3.rds"))
OUT_FIG <- file.path(ROOT, "output/figures")
OUT_TAB <- file.path(ROOT, "output/tables")
OUT_LOG <- file.path(ROOT, "output/logs")
OUT_REP <- file.path(ROOT, "output/reports")
invisible(lapply(c(OUT_FIG, OUT_TAB, OUT_LOG, OUT_REP), dir.create,
                 recursive = TRUE, showWarnings = FALSE))

W_MM <- 174
H_MM <- 150
EXPECTED_N <- 4648L
EXPECTED_EVENTS <- 408L
DAY <- 1 / 365.25
RISK_TIMES <- 0:5
CI_TIMES <- c(1, 3, 5)
GROUP_ORDER <- as.character(TRAJ_LEVELS)                  # 4, 1, 2, 3
GROUP_COL <- c(`4` = "#CC79A7", `1` = "#D55E00",
               `2` = "#0072B2", `3` = "#009E73")
GROUP_EN <- c(`4` = "Group 4 (reference)", `1` = "Group 1",
              `2` = "Group 2", `3` = "Group 3")
GROUP_ROW_LAB <- c(`4` = "Group 4 (ref.)", `1` = "Group 1",
                   `2` = "Group 2", `3` = "Group 3")

# --------------------------------------------------------------- build data
x <- readRDS(INPUT_RDS)
a <- x[["analytic"]]
sl <- x[["status_long"]]
d <- a[a[["main_cohort"]] %in% TRUE & a[["CKM_any_prevalent"]] %in% 0, ]
stopifnot(nrow(d) == EXPECTED_N)
d[["person_id"]] <- as.character(d[["person_id"]])

# the stored analytic time equals the status_long time for the EHR censoring rule
s <- sl[sl[["outcome"]] == "CKM_any" & sl[["censoring"]] == "ehr", ]
stopifnot(!anyDuplicated(s[["person_id"]]))
nb <- nrow(d)
d <- merge(d, s[, c("person_id", "event_date", "time_yr", "event")],
           by = "person_id", all.x = TRUE, sort = FALSE, validate = "1:1")
stopifnot(nrow(d) == nb, !anyNA(d[["time_yr"]]),
          sum(d[["event"]]) == EXPECTED_EVENTS,
          max(abs(d[["time_yr"]] - d[["CKM_any_time_yr"]])) == 0)

grp <- as.data.frame(readr::read_csv(TRAJ_GROUP_FILE, show_col_types = FALSE))
grp <- grp[, c(TRAJ_ID_COL, TRAJ_GROUP_COL)]
names(grp) <- c("person_id", "traj_group")
grp[["person_id"]] <- as.character(grp[["person_id"]])
stopifnot(!anyDuplicated(grp[["person_id"]]))
d <- merge(d, grp, by = "person_id", all.x = TRUE, sort = FALSE, validate = "1:1")
stopifnot(nrow(d) == nb, !anyNA(d[["traj_group"]]))
d[["traj_group"]] <- as_traj_group(d[["traj_group"]])

d[["t_event"]] <- ifelse(d[["event"]] %in% 1,
                         as.numeric(d[["event_date"]] - d[["landmark"]]) / 365.25,
                         NA_real_)
d[["t_death"]] <- ifelse(d[["death_after_lm"]] %in% 1 & !is.na(d[["death_date"]]),
                         as.numeric(d[["death_date"]] - d[["landmark"]]) / 365.25,
                         NA_real_)
d[["t_cens"]] <- as.numeric(pmin(d[["fu_end_ehr"]], d[["fu_end_admin"]]) -
                              d[["landmark"]]) / 365.25
d[["t"]] <- pmin(d[["t_event"]], d[["t_death"]], d[["t_cens"]], na.rm = TRUE)
stopifnot(!anyNA(d[["t"]]), all(d[["t"]] > 0),
          max(abs(d[["t"]] - d[["CKM_any_time_yr"]])) == 0)

d[["status"]] <- ifelse(!is.na(d[["t_event"]]) & d[["t_event"]] <= d[["t"]] + 1e-9, 1L,
                        ifelse(!is.na(d[["t_death"]]) & d[["t_death"]] <= d[["t"]] + 1e-9,
                               2L, 0L))
stopifnot(sum(d[["status"]] %in% 1) == EXPECTED_EVENTS,
          all(d[["death_after_lm"]][d[["status"]] %in% 2] %in% 1))

n_death_total <- sum(d[["death_after_lm"]] %in% 1)
n_death_competing <- sum(d[["status"]] %in% 2)
n_death_after_event <- sum(d[["status"]] %in% 1 & d[["death_after_lm"]] %in% 1)
n_death_after_cens <- n_death_total - n_death_competing - n_death_after_event
stopifnot(n_death_after_cens >= 0L)

status_tab <- data.frame(
  category = c("incident CKM disease (status 1)",
               "death as competing event (status 2)",
               "censored at last EHR visit or CDR cut-off (status 0)",
               "deaths after the landmark, total",
               "  death before any CKM event and before censoring",
               "  death after a CKM event",
               "  death after the censoring date (censored first)"),
  n = c(sum(d[["status"]] %in% 1), n_death_competing, sum(d[["status"]] %in% 0),
        n_death_total, n_death_competing, n_death_after_event, n_death_after_cens),
  stringsAsFactors = FALSE
)
readr::write_excel_csv(status_tab, file.path(OUT_TAB, "11b_ckm_event_status.csv"))

# --------------------------------------------------------------- fit
d[["st"]] <- factor(d[["status"]], levels = 0:2,
                    labels = c("censored", "ckm", "death"))
fit <- survival::survfit(survival::Surv(t, st) ~ traj_group, data = d,
                         conf.type = "log-log")
stopifnot(inherits(fit, "survfitms"), !is.null(fit[["pstate"]]),
          identical(colnames(fit[["pstate"]]), c("(s0)", "ckm", "death")))

strata_n <- fit[["strata"]]
strata_lab <- rep(sub("^traj_group=", "", names(strata_n)), strata_n)
stopifnot(length(strata_lab) == length(fit[["time"]]))
cum <- data.frame(group = strata_lab,
                  time = as.numeric(fit[["time"]]),
                  cuminc = as.numeric(fit[["pstate"]][, "ckm"]),
                  lower = as.numeric(fit[["lower"]][, "ckm"]),
                  upper = as.numeric(fit[["upper"]][, "ckm"]),
                  stringsAsFactors = FALSE)
cum <- cum[order(match(cum[["group"]], GROUP_ORDER), cum[["time"]]), ]
# explicit origin for the step function
cum <- rbind(
  data.frame(group = GROUP_ORDER, time = 0, cuminc = 0, lower = 0, upper = 0,
             stringsAsFactors = FALSE),
  cum
)
cum <- cum[order(match(cum[["group"]], GROUP_ORDER), cum[["time"]]), ]
stopifnot(!anyNA(cum[, c("cuminc", "lower", "upper")]),
          all(cum[["lower"]] <= cum[["cuminc"]] + 1e-12),
          all(cum[["cuminc"]] <= cum[["upper"]] + 1e-12))
for (g in GROUP_ORDER) {
  v <- cum[["cuminc"]][cum[["group"]] == g]
  stopifnot(all(diff(v) >= -1e-12))
}
readr::write_excel_csv(cum, file.path(OUT_TAB, "11b_ckm_cuminc_curves.csv"))

# ------------------------------------------------- estimates at 1, 3 and 5 years
sm <- summary(fit, times = CI_TIMES, extend = TRUE)
# in summary.survfit the strata element is already expanded to one entry per row
sm_lab <- sub("^traj_group=", "", as.character(sm[["strata"]]))
stopifnot(length(sm_lab) == nrow(sm[["pstate"]]))
ci <- data.frame(group = sm_lab, time = as.numeric(sm[["time"]]),
                 cuminc = as.numeric(sm[["pstate"]][, "ckm"]),
                 lower = as.numeric(sm[["lower"]][, "ckm"]),
                 upper = as.numeric(sm[["upper"]][, "ckm"]),
                 stringsAsFactors = FALSE)
ci <- ci[order(match(ci[["group"]], GROUP_ORDER), ci[["time"]]), ]
stopifnot(nrow(ci) == 12L, !anyNA(ci))
readr::write_excel_csv(ci, file.path(OUT_TAB, "11b_ckm_cuminc_1_3_5y.csv"))

# ------------------------------------------------------------------ risk table
risk <- do.call(rbind, lapply(GROUP_ORDER, function(g) {
  tt <- d[["t"]][as.character(d[["traj_group"]]) == g]
  data.frame(group = g, time = RISK_TIMES,
             n_at_risk = vapply(RISK_TIMES, function(k) sum(tt >= k), integer(1)),
             stringsAsFactors = FALSE)
}))
risk <- risk[order(match(risk[["group"]], GROUP_ORDER), risk[["time"]]), ]
group_sizes <- as.integer(table(as.character(d[["traj_group"]]))[
  match(GROUP_ORDER, names(table(as.character(d[["traj_group"]]))))])
stopifnot(all(risk[["n_at_risk"]][risk[["time"]] == 0] == group_sizes),
          sum(group_sizes) == EXPECTED_N)
readr::write_excel_csv(risk, file.path(OUT_TAB, "11b_ckm_risk_table.csv"))

max_t_group <- vapply(GROUP_ORDER, function(g) max(d[["t"]][as.character(d[["traj_group"]]) == g]),
                      numeric(1))
max_t <- max(max_t_group)
stopifnot(all(max_t_group > 5))

# ------------------------------------------------------------------- figure
cv <- as.data.frame(readr::read_csv(file.path(OUT_TAB, "11b_ckm_cuminc_curves.csv"),
                                    show_col_types = FALSE))
rk <- as.data.frame(readr::read_csv(file.path(OUT_TAB, "11b_ckm_risk_table.csv"),
                                    show_col_types = FALSE))
stopifnot(nrow(cv) > 0L, nrow(rk) == 4L * length(RISK_TIMES))
y_max <- ceiling(max(cv[["cuminc"]]) * 100 / 5) * 5 + 5
if (y_max < 4) y_max <- 4

draw_cif <- function() {
  layout(matrix(1:2, ncol = 1), heights = c(3.5, 1.4))
  op <- par(mar = c(2.9, 6.6, 1.0, 1.0), las = 1, family = "Arial",
            cex.axis = 0.75, cex.lab = 0.85, mgp = c(2.2, 0.5, 0))
  on.exit(par(op))
  plot(NA, xlim = c(0, max_t), ylim = c(0, y_max), axes = FALSE,
       xlab = "", ylab = "Cumulative incidence of CKM disease (%)")
  for (g in rev(GROUP_ORDER)) {
    s <- cv[cv[["group"]] == g, ]
    lines(c(s[["time"]], max_t), c(s[["cuminc"]], tail(s[["cuminc"]], 1L)) * 100,
          type = "s", col = GROUP_COL[[g]], lwd = 1.15)
  }
  axis(1, at = 0:floor(max_t), tcl = -0.25, cex.axis = 0.75)
  axis(2, at = seq(0, y_max, by = 5), tcl = -0.25, cex.axis = 0.75)
  box(col = "grey35", lwd = 0.5)
  mtext("Years since the landmark", side = 1, line = 1.55, cex = 0.85,
        family = "Arial")
  legend("topleft", legend = unname(GROUP_EN[GROUP_ORDER]),
         col = unname(GROUP_COL[GROUP_ORDER]), lty = 1, lwd = 1.15,
         bty = "n", cex = 0.72, seg.len = 1.6, inset = c(0.02, 0.02))

  par(mar = c(0.1, 6.6, 0, 1.0))
  plot.new()
  plot.window(xlim = c(0, max_t), ylim = c(0.2, 5.2))
  text(0, 4.95, "Number at risk", adj = c(0, 0.5), cex = 0.78, font = 2,
       xpd = NA)
  segments(0, 4.72, max_t, 4.72, col = "grey70", lwd = 0.4)
  for (k in seq_along(GROUP_ORDER)) {
    g <- GROUP_ORDER[k]
    yy <- 4.05 - (k - 1) * 0.85
    text(-0.02 * max_t, yy, unname(GROUP_ROW_LAB[[g]]), adj = c(1, 0.5),
         cex = 0.72, xpd = NA)
    r <- rk[rk[["group"]] == g, ]
    for (j in seq_along(RISK_TIMES)) {
      text(RISK_TIMES[j], yy, formatC(r[["n_at_risk"]][j], format = "d", big.mark = ","),
           adj = c(if (j == 1L) 0 else 0.5, 0.5), cex = 0.72)
    }
  }
  invisible(NULL)
}

grDevices::cairo_pdf(file.path(OUT_FIG, "11b_fig_ckm_cuminc.pdf"),
                     width = W_MM / 25.4, height = H_MM / 25.4,
                     family = "Arial", onefile = TRUE)
draw_cif()
invisible(grDevices::dev.off())

ragg::agg_tiff(file.path(OUT_FIG, "11b_fig_ckm_cuminc.tiff"),
               width = W_MM / 25.4, height = H_MM / 25.4, units = "in",
               res = 600, background = "white")
draw_cif()
invisible(grDevices::dev.off())

# ------------------------------------------------------------------- report
pct <- function(x) sprintf("%.2f", 100 * x)
rep_lines <- c(
  "# Step 11b - cumulative incidence of CKM disease by trajectory group",
  "",
  paste0("Generated by `R/11b_figure_ckm_cuminc.R` on ",
         format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ". Version v3.3."),
  "",
  paste0("Population: CKM disease at-risk set of the main cohort, n = ",
         format(EXPECTED_N, big.mark = ","), "; events = ", EXPECTED_EVENTS, "."),
  "Estimator: Aalen-Johansen. Death is handled as a competing event.",
  "",
  "## Definition of the endpoint",
  "",
  "- Time origin: the Fitbit landmark.",
  "- Follow-up time: minimum of the CKM event date, the death date, the last EHR visit and the CDR cut-off, in years.",
  "- Status 1: incident CKM disease.",
  "- Status 2: death before any CKM event and before the censoring date (competing event).",
  "- Status 0: censored at the last EHR visit or at the CDR cut-off, whichever is earlier (existing rule).",
  "- A death that occurs after the last EHR visit cannot be observed as a competing event, because the existing rule censors the person at the last EHR visit first. Such deaths are counted below.",
  "- 95% confidence intervals come from the Aalen-Johansen variance estimate (`survfit`, `conf.type = \"log-log\"`).",
  "",
  "## Event and censoring counts",
  "",
  "| Category | n |", "| --- | --- |",
  paste0("| ", status_tab[["category"]], " | ",
         formatC(status_tab[["n"]], format = "d", big.mark = ","), " |"),
  "",
  "## Maximum follow-up by group",
  "",
  "| Group | n | Maximum follow-up (years) |", "| --- | --- | --- |",
  vapply(seq_along(GROUP_ORDER), function(k) {
    g <- GROUP_ORDER[k]
    paste0("| ", GROUP_EN[[g]], " | ", formatC(group_sizes[k], format = "d", big.mark = ","),
           " | ", sprintf("%.2f", max_t_group[k]), " |")
  }, character(1)),
  "",
  "## Cumulative incidence at 1, 3 and 5 years",
  "",
  "| Group | Year | Cumulative incidence (%) | 95% CI (%) |", "| --- | --- | --- | --- |",
  vapply(seq_len(nrow(ci)), function(i) {
    paste0("| ", GROUP_EN[[ci[["group"]][i]]], " | ", ci[["time"]][i], " | ",
           pct(ci[["cuminc"]][i]), " | ", pct(ci[["lower"]][i]), " to ",
           pct(ci[["upper"]][i]), " |")
  }, character(1)),
  "",
  "## Files",
  "",
  "- `output/figures/11b_fig_ckm_cuminc.pdf`, `.tiff` (600 dpi)",
  "- `output/tables/11b_ckm_cuminc_curves.csv` (step-function values used by the figure)",
  "- `output/tables/11b_ckm_cuminc_1_3_5y.csv`",
  "- `output/tables/11b_ckm_risk_table.csv`",
  "- `output/tables/11b_ckm_event_status.csv`"
)
writeLines(rep_lines, file.path(OUT_REP, "step11b_ckm_cuminc.md"), useBytes = TRUE)

jsonlite::write_json(list(
  script = "R/11b_figure_ckm_cuminc.R",
  run_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  version = "v3.3",
  input_file = INPUT_RDS,
  population = "main_cohort == TRUE and CKM_any_prevalent == 0",
  n = EXPECTED_N,
  events = sum(d[["status"]] %in% 1),
  competing_events_death = n_death_competing,
  censored = sum(d[["status"]] %in% 0),
  deaths_total = n_death_total,
  deaths_after_event = n_death_after_event,
  deaths_after_censoring = n_death_after_cens,
  time_definition = "min(event date, death date, last EHR visit, CDR cut-off) minus landmark, in years; identical to the stored CKM_any_time_yr for all persons",
  estimator = "Aalen-Johansen, survival::survfit multi-state, conf.type = log-log",
  max_follow_up_by_group = as.list(stats::setNames(round(max_t_group, 4), GROUP_ORDER)),
  risk_times = RISK_TIMES,
  ci_times = CI_TIMES,
  width_mm = W_MM,
  height_mm = H_MM,
  outputs = c("output/figures/11b_fig_ckm_cuminc.pdf",
              "output/figures/11b_fig_ckm_cuminc.tiff",
              "output/tables/11b_ckm_cuminc_curves.csv",
              "output/tables/11b_ckm_cuminc_1_3_5y.csv",
              "output/tables/11b_ckm_risk_table.csv",
              "output/tables/11b_ckm_event_status.csv",
              "output/reports/step11b_ckm_cuminc.md"),
  r_version = R.version.string,
  git_commit = NULL
), file.path(OUT_LOG, "11b_provenance.json"), auto_unbox = TRUE, pretty = TRUE)
