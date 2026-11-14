# =============================================================================
# Step 13: final versions of Figures 2, 3, 4 and S3
#   Final numbering (2026-11-14): Figure 3 is the cumulative-incidence curve and
#   Figure 4 is the main forest plot; the two stems were swapped on that date and
#   the script now writes the swapped names.
#   Redraw only. Every number is read from a stored CSV written by an earlier
#   run: 09b fitted trajectories, 09c observed means by quarter, 07d main hazard
#   ratios, 07d sensitivity hazard ratios, 11b cumulative-incidence curves and
#   risk table. No model is fitted, no estimate is recomputed, no value changed.
#   Changes against the previous versions:
#     Figure 2  fitted lines kept, observed means from 09c added as points,
#               legend carries the group name with n (%), no in-figure title,
#               x axis "Quarter since Fitbit start".
#     Figure 4  x axis "Hazard ratio (95% CI)" with ticks 0.75/1.0/1.5/2.0,
#               "body size" replaced by "BMI", English group names on the y axis.
#     Figure S3 same axis change and English group names; analysis codes removed.
#     Figure 3  English group names, x axis truncated at 5 years, y axis
#               "Cumulative incidence of CKM conditions (%)".
#   Shared palette: group 4 (reference) pink, 1 vermillion, 2 blue, 3 green.
#   Arial, PDF with embedded fonts, TIFF 600 dpi. Files are written straight
#   into output/figures/final/ and overwrite the previous finals.
#   No console output (project rule).
# =============================================================================
suppressPackageStartupMessages({
  library(readr)
  library(ggplot2)
  library(jsonlite)
})

# All paths are relative to a single root. Override the root with the environment
# variable FIG_ROOT, or edit the default below; the layout under the root is
# output/tables (inputs), output/figures, output/logs and output/reports.
ROOT <- Sys.getenv("FIG_ROOT", unset = ".")
OUT_FIG <- file.path(ROOT, "output/figures")
OUT_TAB <- file.path(ROOT, "output/tables")
OUT_LOG <- file.path(ROOT, "output/logs")
OUT_REP <- file.path(ROOT, "output/reports")
FINAL <- file.path(OUT_FIG, "final")
invisible(lapply(c(OUT_FIG, OUT_TAB, OUT_LOG, OUT_REP, FINAL), dir.create,
                 recursive = TRUE, showWarnings = FALSE))

# optional command-line subset: Rscript R/13_final_figures.R Figure2_trajectory_G4
ARG_FIGS <- commandArgs(trailingOnly = TRUE)
WANT <- function(stem) {
  length(ARG_FIGS) == 0L ||
    any(vapply(ARG_FIGS, function(a) grepl(a, stem, fixed = TRUE), logical(1)))
}

# ------------------------------------------------------------- shared contract
# Display order is G1 to G4; the colour stays bound to the group code, so the
# palette does not move when the display order changes. G4 is the reference.
GROUP_ORDER <- c("1", "2", "3", "4")
GROUP_NAME <- c(
  `4` = "Moderate activity\u2013low heart rate",
  `1` = "High activity\u2013short sleep",
  `2` = "Moderately high activity\u2013intermediate heart rate",
  `3` = "Low activity\u2013high heart rate"
)
GROUP_COL <- c(`4` = "#CC79A7", `1` = "#D55E00", `2` = "#0072B2", `3` = "#009E73")
EXPECTED_G1 <- 1553L
EXPECTED_G2 <- 2032L
EXPECTED_G3 <- 2112L
EXPECTED_G4 <- 2073L
EXPECTED_N <- 7770L

# -------------------------------------------------------------- Figure 2 data
MODEL_FIT <- file.path(OUT_TAB, "09b_g4_fitted_trajectories_long.csv")
OBS_FILE <- file.path(OUT_TAB, "09c_observed_by_group_quarter.csv")
fit <- as.data.frame(readr::read_csv(MODEL_FIT, show_col_types = FALSE))
obs <- as.data.frame(readr::read_csv(OBS_FILE, show_col_types = FALSE))
stopifnot(nrow(fit) == 96L, nrow(obs) == 96L)
fit[["grp"]] <- as.character(fit[["group"]])
obs[["grp"]] <- as.character(obs[["traj_group"]])

d2 <- merge(fit[, c("grp", "indicator", "period", "fitted_value")],
            obs[, c("grp", "indicator", "period", "mean", "n_nonmissing",
                    "n_group_period")],
            by = c("grp", "indicator", "period"), validate = "1:1")
stopifnot(nrow(d2) == 96L, !anyNA(d2[, c("fitted_value", "mean")]))

grp_n <- tapply(obs[["n_group_period"]], obs[["grp"]], function(x) unique(x))
stopifnot(all(lengths(grp_n) == 1L))
grp_n <- vapply(GROUP_ORDER, function(g) as.integer(grp_n[[g]]), integer(1))
stopifnot(identical(unname(grp_n),
                    c(EXPECTED_G1, EXPECTED_G2, EXPECTED_G3, EXPECTED_G4)),
          sum(grp_n) == EXPECTED_N)
grp_pct <- round(100 * grp_n / EXPECTED_N, 1)
stopifnot(identical(unname(grp_pct), c(20.0, 26.2, 27.2, 26.7)))

fmt_n <- function(x) formatC(x, format = "d", big.mark = ",")
leg2 <- vapply(GROUP_ORDER, function(g) {
  nm <- unname(GROUP_NAME[[g]])
  if (g == "4") {
    sprintf("%s (reference; n = %s, %.1f%%)", nm, fmt_n(grp_n[[g]]), grp_pct[[g]])
  } else {
    sprintf("%s (n = %s, %.1f%%)", nm, fmt_n(grp_n[[g]]), grp_pct[[g]])
  }
}, character(1))
stopifnot(identical(unname(leg2), c(
  "High activity\u2013short sleep (n = 1,553, 20.0%)",
  "Moderately high activity\u2013intermediate heart rate (n = 2,032, 26.2%)",
  "Low activity\u2013high heart rate (n = 2,112, 27.2%)",
  "Moderate activity\u2013low heart rate (reference; n = 2,073, 26.7%)")))

d2[["glab"]] <- factor(d2[["grp"]], levels = GROUP_ORDER)
IND_LEVELS <- c("steps_mean", "hr_mean", "sleep_mean")
IND_TITLE <- c(steps_mean = "Steps per day (steps/day)",
               hr_mean = "Mean heart rate (bpm)",
               sleep_mean = "Sleep duration (min/day)")
d2[["ind"]] <- factor(unname(IND_TITLE[d2[["indicator"]]]), levels = unname(IND_TITLE))

p2 <- ggplot(d2, aes(x = period)) +
  geom_line(aes(y = fitted_value, colour = glab), linewidth = 0.8) +
  geom_point(aes(y = mean, colour = glab), size = 1.4, shape = 16) +
  facet_wrap(~ind, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = GROUP_COL, breaks = GROUP_ORDER, labels = leg2,
                      name = NULL, drop = FALSE,
                      guide = guide_legend(ncol = 1, byrow = TRUE)) +
  scale_x_continuous(breaks = 0:7, limits = c(0, 7),
                     expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = "Quarter since Fitbit start", y = NULL) +
  theme_bw(base_size = 9, base_family = "Arial") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        panel.border = element_rect(colour = "grey35", linewidth = 0.4),
        strip.background = element_rect(fill = "grey96", colour = "grey70",
                                        linewidth = 0.4),
        strip.text = element_text(size = 8.5, family = "Arial"),
        axis.text = element_text(colour = "black", size = 8, family = "Arial"),
        axis.title = element_text(size = 8.5, family = "Arial"),
        legend.text = element_text(size = 7.0, family = "Arial"),
        legend.key.width = unit(7, "mm"),
        legend.spacing.y = unit(0.4, "mm"),
        legend.position = "bottom",
        legend.margin = margin(1, 2, 1, 2),
        plot.margin = margin(5, 8, 4, 6))

# -------------------------------------------------------------- Figure 4 data
MAIN_HR <- file.path(OUT_TAB, "07d_table_main_hr_v33.csv")
SENS_HR <- file.path(OUT_TAB, "07d_table_supp_sensitivity_v33.csv")
hr <- as.data.frame(readr::read_csv(MAIN_HR, show_col_types = FALSE))
hr <- hr[hr[["outcome"]] == "CKM_any", ]
stopifnot(nrow(hr) == 9L, all(hr[["model"]] %in% c("M1", "M2", "M3")),
          all(hr[["group"]] %in% c("1", "2", "3")),
          !anyNA(hr[, c("hr", "ci_low", "ci_high")]))

MODEL_EN <- c(M1 = "M1: demographic and temporal",
              M2 = "M2: primary model",
              M3 = "M3: plus BMI, blood pressure, medication")
MODEL_COL <- c(M1 = "#E69F00", M2 = "#0072B2", M3 = "#009E73")
MODEL_SHAPE <- c(M1 = 17, M2 = 16, M3 = 15)
Y_LAB <- c(`1` = "High activity\u2013\nshort sleep",
           `2` = "Moderately high activity\u2013\nintermediate heart rate",
           `3` = "Low activity\u2013\nhigh heart rate")

d3 <- hr
d3[["model_f"]] <- factor(d3[["model"]], levels = names(MODEL_EN))
d3[["y"]] <- as.numeric(d3[["group"]])
d3[["y_plot"]] <- d3[["y"]] + unname(c(M1 = -0.22, M2 = 0, M3 = 0.22)[d3[["model"]]])
stopifnot(nrow(d3) == 9L, !anyNA(d3[["y_plot"]]))

x_scale <- scale_x_log10(breaks = c(0.75, 1, 1.5, 2),
                         labels = c("0.75", "1.0", "1.5", "2.0"),
                         limits = c(0.6, 2.4),
                         expand = expansion(mult = c(0.03, 0.03)))
dir_note <- list(
  annotate("text", x = 0.615, y = -Inf, label = "\u2190 lower risk", hjust = 0,
           vjust = -0.9, size = 2.6, family = "Arial", colour = "grey30"),
  annotate("text", x = 2.35, y = -Inf, label = "higher risk \u2192", hjust = 1,
           vjust = -0.9, size = 2.6, family = "Arial", colour = "grey30"))
theme_forest <- function() {
  theme_bw(base_size = 9, base_family = "Arial") +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3),
          panel.border = element_rect(colour = "grey35", linewidth = 0.4),
          axis.text = element_text(colour = "black", size = 8.2, family = "Arial"),
          axis.title = element_text(size = 9, family = "Arial"),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.text = element_text(size = 7.8, family = "Arial"),
          legend.key = element_blank(),
          legend.position = "bottom",
          plot.margin = margin(6, 8, 4, 10))
}

p3 <- ggplot(d3, aes(x = hr, y = y_plot, colour = model_f, shape = model_f)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey30", linewidth = 0.4) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0, linewidth = 0.5,
                orientation = "y") +
  geom_point(size = 2.1) +
  x_scale + dir_note +
  scale_y_continuous(breaks = 1:3, labels = unname(Y_LAB), limits = c(0.4, 3.6),
                     expand = c(0, 0)) +
  scale_colour_manual(values = MODEL_COL, name = NULL, labels = MODEL_EN, drop = FALSE) +
  scale_shape_manual(values = MODEL_SHAPE, name = NULL, labels = MODEL_EN, drop = FALSE) +
  guides(colour = guide_legend(nrow = 3, byrow = TRUE),
         shape = guide_legend(nrow = 3, byrow = TRUE)) +
  labs(x = "Hazard ratio (95% CI)") +
  theme_forest()

# ------------------------------------------------------------- Figure S3 data
sens <- as.data.frame(readr::read_csv(SENS_HR, show_col_types = FALSE))
stopifnot(nrow(sens) == 4L, all(sens[["outcome"]] == "CKM_any"))
main_m2 <- hr[hr[["model"]] == "M2", ]
stopifnot(nrow(main_m2) == 3L)
long_g <- do.call(rbind, lapply(seq_len(nrow(sens)), function(i) {
  data.frame(analysis = sens[["analysis"]][i],
             group = c("1", "2", "3"),
             hr = c(sens[["hr_g1"]][i], sens[["hr_g2"]][i], sens[["hr_g3"]][i]),
             ci_low = c(sens[["lo_g1"]][i], sens[["lo_g2"]][i], sens[["lo_g3"]][i]),
             ci_high = c(sens[["hi_g1"]][i], sens[["hi_g2"]][i], sens[["hi_g3"]][i]),
             stringsAsFactors = FALSE)
}))
d4 <- rbind(data.frame(analysis = "Main",
                       group = as.character(main_m2[["group"]]),
                       hr = main_m2[["hr"]], ci_low = main_m2[["ci_low"]],
                       ci_high = main_m2[["ci_high"]], stringsAsFactors = FALSE),
            long_g)
stopifnot(nrow(d4) == 15L, !anyNA(d4[, c("hr", "ci_low", "ci_high")]))
ANALYSIS_EN <- c(Main = "Primary model",
                 S1 = "Excluding events in the first 6 months",
                 S2 = "Restricted to baseline disease-free",
                 S3 = "Restricted to posterior probability \u2265 0.70",
                 S4 = "Additionally adjusted for COVID-19 fraction")
d4[["analysis_f"]] <- factor(d4[["analysis"]], levels = names(ANALYSIS_EN),
                             labels = unname(ANALYSIS_EN))
d4[["grp_f"]] <- factor(d4[["group"]], levels = c("1", "2", "3"))
d4[["y"]] <- as.numeric(d4[["analysis_f"]])
d4[["y_plot"]] <- d4[["y"]] +
  unname(c(`1` = -0.22, `2` = 0, `3` = 0.22)[d4[["group"]]])

pS3 <- ggplot(d4, aes(x = hr, y = y_plot, colour = grp_f, shape = grp_f)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey30", linewidth = 0.4) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0, linewidth = 0.5,
                orientation = "y") +
  geom_point(size = 2.0) +
  x_scale + dir_note +
  scale_y_continuous(breaks = 1:5, labels = unname(ANALYSIS_EN),
                     limits = c(0.4, 5.6), expand = c(0, 0)) +
  scale_colour_manual(values = GROUP_COL[c("1", "2", "3")], name = NULL,
                      labels = unname(GROUP_NAME[c("1", "2", "3")]), drop = FALSE) +
  scale_shape_manual(values = c(`1` = 16, `2` = 17, `3` = 15), name = NULL,
                     labels = unname(GROUP_NAME[c("1", "2", "3")]), drop = FALSE) +
  guides(colour = guide_legend(nrow = 3, byrow = TRUE),
         shape = guide_legend(nrow = 3, byrow = TRUE)) +
  labs(x = "Hazard ratio (95% CI)") +
  theme_forest()

# ---------------------------------------------------------------- save helper
save_pdf_tiff <- function(p, stem, w_mm, h_mm) {
  grDevices::cairo_pdf(file.path(FINAL, paste0(stem, ".pdf")),
                       width = w_mm / 25.4, height = h_mm / 25.4,
                       family = "Arial", onefile = TRUE)
  print(p)
  invisible(grDevices::dev.off())
  ragg::agg_tiff(file.path(FINAL, paste0(stem, ".tiff")),
                 width = w_mm / 25.4, height = h_mm / 25.4, units = "in",
                 res = 600, background = "white")
  print(p)
  invisible(grDevices::dev.off())
}

if (WANT("Figure2_trajectory_G4")) save_pdf_tiff(p2, "Figure2_trajectory_G4", 174, 122)
if (WANT("Figure4_forest_main")) save_pdf_tiff(p3, "Figure4_forest_main", 152, 116)
if (WANT("FigureS3_forest_sensitivity")) {
  save_pdf_tiff(pS3, "FigureS3_forest_sensitivity", 172, 128)
}

# -------------------------------------------------------------- Figure 3 data
CURVES <- file.path(OUT_TAB, "11b_ckm_cuminc_curves.csv")
RISK <- file.path(OUT_TAB, "11b_ckm_risk_table.csv")
cv <- as.data.frame(readr::read_csv(CURVES, show_col_types = FALSE))
rk <- as.data.frame(readr::read_csv(RISK, show_col_types = FALSE))
cv[["group"]] <- as.character(cv[["group"]])
rk[["group"]] <- as.character(rk[["group"]])
stopifnot(nrow(rk) == 24L, all(GROUP_ORDER %in% cv[["group"]]),
          max(cv[["time"]]) >= 5)
RISK_TIMES <- 0:5
X_MAX <- 5
y_max <- ceiling(max(cv[["cuminc"]][cv[["time"]] <= X_MAX]) * 100 / 5) * 5 + 5
if (y_max < 5) y_max <- 5
stopifnot(nrow(rk[rk[["time"]] %in% RISK_TIMES, ]) == 24L)

GROUP_EN <- c(`4` = "Moderate activity\u2013low heart rate (reference)",
              `1` = "High activity\u2013short sleep",
              `2` = "Moderately high activity\u2013intermediate heart rate",
              `3` = "Low activity\u2013high heart rate")
GROUP_ROW <- c(`4` = "Moderate activity\u2013\nlow heart rate (ref.)",
               `1` = "High activity\u2013\nshort sleep",
               `2` = "Moderately high activity\u2013\nintermediate heart rate",
               `3` = "Low activity\u2013\nhigh heart rate")

draw_cif <- function() {
  layout(matrix(1:2, ncol = 1), heights = c(3.5, 1.6))
  op <- par(mar = c(2.9, 7.8, 1.0, 1.0), las = 1, family = "Arial",
            cex.axis = 0.75, cex.lab = 0.85, mgp = c(2.2, 0.5, 0))
  on.exit(par(op))
  plot(NA, xlim = c(0, X_MAX), ylim = c(0, y_max), axes = FALSE, xlab = "",
       ylab = "Cumulative incidence of CKM conditions (%)")
  for (g in rev(GROUP_ORDER)) {
    s <- cv[cv[["group"]] == g & cv[["time"]] <= X_MAX, ]
    stopifnot(nrow(s) > 1L, !is.unsorted(s[["time"]]))
    lines(c(s[["time"]], X_MAX), c(s[["cuminc"]], tail(s[["cuminc"]], 1L)) * 100,
          type = "s", col = GROUP_COL[[g]], lwd = 1.15)
  }
  axis(1, at = 0:X_MAX, tcl = -0.25, cex.axis = 0.75)
  axis(2, at = seq(0, y_max, by = 5), tcl = -0.25, cex.axis = 0.75)
  box(col = "grey35", lwd = 0.5)
  mtext("Years since the landmark", side = 1, line = 1.55, cex = 0.85,
        family = "Arial")
  legend("topleft", legend = unname(GROUP_EN[GROUP_ORDER]),
         col = unname(GROUP_COL[GROUP_ORDER]), lty = 1, lwd = 1.15,
         bty = "n", cex = 0.62, seg.len = 1.5, inset = c(0.01, 0.02))

  par(mar = c(0.1, 7.8, 0, 1.0))
  plot.new()
  plot.window(xlim = c(0, X_MAX), ylim = c(0.2, 5.4))
  text(0, 5.05, "Number at risk", adj = c(0, 0.5), cex = 0.78, font = 2, xpd = NA)
  segments(0, 4.80, X_MAX, 4.80, col = "grey70", lwd = 0.4)
  for (k in seq_along(GROUP_ORDER)) {
    g <- GROUP_ORDER[k]
    yy <- 4.05 - (k - 1) * 0.85
    text(-0.02 * X_MAX, yy, unname(GROUP_ROW[[g]]), adj = c(1, 0.5),
         cex = 0.62, xpd = NA)
    r <- rk[rk[["group"]] == g & rk[["time"]] %in% RISK_TIMES, ]
    r <- r[order(r[["time"]]), ]
    stopifnot(nrow(r) == length(RISK_TIMES))
    for (j in seq_along(RISK_TIMES)) {
      text(RISK_TIMES[j], yy, formatC(r[["n_at_risk"]][j], format = "d", big.mark = ","),
           adj = c(if (j == 1L) 0 else 0.5, 0.5), cex = 0.66)
    }
  }
  invisible(NULL)
}

if (WANT("Figure3_ckm_cumulative_incidence")) {
  grDevices::cairo_pdf(file.path(FINAL, "Figure3_ckm_cumulative_incidence.pdf"),
                       width = 174 / 25.4, height = 150 / 25.4,
                       family = "Arial", onefile = TRUE)
  draw_cif()
  invisible(grDevices::dev.off())
  ragg::agg_tiff(file.path(FINAL, "Figure3_ckm_cumulative_incidence.tiff"),
                 width = 174 / 25.4, height = 150 / 25.4, units = "in",
                 res = 600, background = "white")
  draw_cif()
  invisible(grDevices::dev.off())
}

# ------------------------------------------------------------ label check file
chk <- rbind(
  data.frame(figure = "Figure 2", item = "legend", value = leg2,
             stringsAsFactors = FALSE),
  data.frame(figure = "Figure 2", item = "x axis",
             value = "Quarter since Fitbit start", stringsAsFactors = FALSE),
  data.frame(figure = "Figure 2", item = "facet",
             value = unname(IND_TITLE), stringsAsFactors = FALSE),
  data.frame(figure = "Figure 4", item = "x axis",
             value = "Hazard ratio (95% CI)", stringsAsFactors = FALSE),
  data.frame(figure = "Figure 4", item = "ticks",
             value = "0.75, 1.0, 1.5, 2.0 (log scale)", stringsAsFactors = FALSE),
  data.frame(figure = "Figure 4", item = "y axis", value = gsub("\n", " / ", Y_LAB),
             stringsAsFactors = FALSE),
  data.frame(figure = "Figure 4", item = "legend", value = unname(MODEL_EN),
             stringsAsFactors = FALSE),
  data.frame(figure = "Figure S3", item = "x axis",
             value = "Hazard ratio (95% CI)", stringsAsFactors = FALSE),
  data.frame(figure = "Figure S3", item = "y axis", value = unname(ANALYSIS_EN),
             stringsAsFactors = FALSE),
  data.frame(figure = "Figure S3", item = "legend",
             value = unname(GROUP_NAME[c("1", "2", "3")]), stringsAsFactors = FALSE),
  data.frame(figure = "Figure 3", item = "y axis",
             value = "Cumulative incidence of CKM conditions (%)",
             stringsAsFactors = FALSE),
  data.frame(figure = "Figure 3", item = "x axis",
             value = "Years since the landmark; 0 to 5", stringsAsFactors = FALSE),
  data.frame(figure = "Figure 3", item = "legend", value = unname(GROUP_EN),
             stringsAsFactors = FALSE)
)
readr::write_excel_csv(chk, file.path(OUT_LOG, "13_figure_labels.csv"))

# observed points versus the 09c table, cell by cell
ov <- data.frame(grp = d2[["grp"]], period = d2[["period"]],
                 indicator = d2[["indicator"]], plotted = d2[["mean"]],
                 stringsAsFactors = FALSE)
ov <- merge(ov, obs[, c("grp", "period", "indicator", "mean")],
            by = c("grp", "period", "indicator"), validate = "1:1")
stopifnot(nrow(ov) == 96L, max(abs(ov[["plotted"]] - ov[["mean"]])) == 0)

# ------------------------------------------------------------------- report
lines <- c(
  "# Step 13 - final versions of Figures 2, 3, 4 and S3",
  "",
  paste0("Generated by `R/13_final_figures.R` on ",
         format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "."),
  "",
  "Redraw only: every number is read from a stored CSV, no model is fitted and",
  "no estimate is recomputed.",
  "",
  "## Sources",
  "",
  "- Figure 2: `output/tables/09b_g4_fitted_trajectories_long.csv` (lines) and",
  "  `output/tables/09c_observed_by_group_quarter.csv` (points).",
  "- Figure 4: `output/tables/07d_table_main_hr_v33.csv` (CKM_any, M1 to M3).",
  "- Figure S3: `output/tables/07d_table_main_hr_v33.csv` (M2) and",
  "  `output/tables/07d_table_supp_sensitivity_v33.csv`.",
  "- Figure 3: `output/tables/11b_ckm_cuminc_curves.csv` and",
  "  `output/tables/11b_ckm_risk_table.csv`.",
  "",
  "## Group labels",
  "",
  paste0("- ", GROUP_ORDER, " = ", unname(GROUP_NAME[GROUP_ORDER])),
  "",
  "Display order is G1, G2, G3, G4 in Figure 2 (legend) and Figure 3 (legend and",
  "risk table). Colours stay bound to the group code: G4 keeps the reference",
  "annotation and its pink colour regardless of position.",
  "",
  paste0("- Figures written in this run: ",
         if (length(ARG_FIGS) == 0L) "all" else paste(ARG_FIGS, collapse = ", "), "."),
  "",
  "## Figure 2 legend",
  "",
  paste0("- ", leg2),
  "",
  "## Checks",
  "",
  paste0("- Group sizes: ", paste(sprintf("%s = %s (%.1f%%)", GROUP_ORDER,
                                         fmt_n(grp_n), grp_pct), collapse = "; "),
         "; total ", fmt_n(sum(grp_n)), "."),
  "- Observed points equal the 09c values in all 96 group-quarter-indicator cells",
  "  (maximum absolute difference 0).",
  "- Figure 3 is truncated at 5 years; the risk table is read for 0 to 5 years.",
  "",
  "## Files (written directly into `output/figures/final/`)",
  "",
  "- `Figure2_trajectory_G4.pdf` / `.tiff`",
  "- `Figure4_forest_main.pdf` / `.tiff`",
  "- `FigureS3_forest_sensitivity.pdf` / `.tiff`",
  "- `Figure3_ckm_cumulative_incidence.pdf` / `.tiff`",
  "- `output/logs/13_figure_labels.csv` (text and axis check)"
)
writeLines(lines, file.path(OUT_REP, "step13_final_figures.md"), useBytes = TRUE)

jsonlite::write_json(list(
  script = "R/13_final_figures.R",
  run_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  model_fitted = FALSE,
  inputs = c(MODEL_FIT, OBS_FILE, MAIN_HR, SENS_HR, CURVES, RISK),
  group_sizes = as.list(stats::setNames(grp_n, GROUP_ORDER)),
  group_percent = as.list(stats::setNames(grp_pct, GROUP_ORDER)),
  observed_points_equal_09c = TRUE,
  figure2 = list(width_mm = 174, height_mm = 122, dpi = 600,
                 x_axis = "Quarter since Fitbit start"),
  figure4 = list(width_mm = 152, height_mm = 116, dpi = 600,
                 x_axis = "Hazard ratio (95% CI)",
                 ticks = c(0.75, 1, 1.5, 2)),
  figure_s3 = list(width_mm = 172, height_mm = 128, dpi = 600,
                   x_axis = "Hazard ratio (95% CI)",
                   ticks = c(0.75, 1, 1.5, 2)),
  figure3 = list(width_mm = 174, height_mm = 150, dpi = 600,
                 x_max_years = 5, y_axis = "Cumulative incidence of CKM conditions (%)"),
  r_version = R.version.string,
  git_commit = NULL
), file.path(OUT_LOG, "13_provenance.json"), auto_unbox = TRUE, pretty = TRUE)
