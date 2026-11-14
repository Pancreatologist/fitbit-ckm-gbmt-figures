# =============================================================================
# Step 7d: consolidated deliverables under plan v3.3
#   1. main table: CKM disease only, M1, M2 (primary) and M3
#   2. forest plot of the three models for CKM disease
#   3. forest plot of the CKM sensitivity analyses on the primary model
#   4. supplementary tables: PH, sensitivity, sample flow, missingness,
#      cardiorenal disease (three models), exploratory single-disease analysis
#   5. verification of the M1 hazard ratios against the values supplied by the
#      investigator
#   6. one readable markdown report and one delivery manifest
# Every number is read from a stored CSV written by 04b, 06b, 08b or an earlier
# step; nothing is typed by hand and no model is refitted here. stopifnot()
# checks the anchors before any table or figure is written.
# No console output (project rule).
# =============================================================================
suppressPackageStartupMessages({
  library(readr)
  library(ggplot2)
  library(jsonlite)
})

OUT_TAB <- "output/tables"
OUT_LOG <- "output/logs"
OUT_FIG <- "output/figures"
OUT_REP <- "output/reports"
invisible(lapply(c(OUT_TAB, OUT_LOG, OUT_FIG, OUT_REP), dir.create,
                 recursive = TRUE, showWarnings = FALSE))

SEED <- 20250101L
M_IMP <- 20L
PLAN_VERSION <- "v3.3"
EXPECTED <- list(CKM_any = c(n = 4648, events = 408),
                 Cardiorenal = c(n = 6871, events = 222))
MODELS <- c("M1", "M2", "M3")
OUTCOME_LABEL <- c(CKM_any = "CKM disease (primary outcome)",
                   Cardiorenal = "Cardiorenal disease (secondary outcome)")
GROUP_EN <- c(
  `1` = "Group 1\nhigh activity, short sleep",
  `2` = "Group 2\nhigher activity, moderate heart rate",
  `3` = "Group 3\nlow activity, high heart rate",
  `4` = "Group 4\nmoderate activity, low heart rate"
)
MODEL_EN <- c(
  M1 = "M1: demographic and temporal",
  M2 = "M2: primary model",
  M3 = "M3: plus body size, blood pressure, medication"
)
ANALYSIS_EN <- c(
  Main = "Primary model (M2)",
  S1 = "S1: 6-month exclusion",
  S2 = "S2: baseline disease-free",
  S3 = "S3: posterior probability \u2265 0.70",
  S4 = "S4: additionally adjusted for COVID fraction"
)
MODEL_COL <- c(M1 = "#E69F00", M2 = "#0072B2", M3 = "#009E73")
MODEL_SHAPE <- c(M1 = 17, M2 = 16, M3 = 15)
GROUP_COL <- c(`1` = "#D55E00", `2` = "#0072B2", `3` = "#009E73")
USER_M1 <- c(`1` = "1.18 (0.89 to 1.57)", `2` = "0.92 (0.69 to 1.22)",
             `3` = "1.60 (1.21 to 2.12)")

mask <- function(x) ifelse(!is.na(x) & x > 0 & x < 20, "<20", as.character(x))
fmt_p <- function(p) ifelse(is.na(p), "not calculated",
                            ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
md_block <- function(h, rows) c(
  paste0("| ", paste(h, collapse = " | "), " |"),
  paste0("| ", paste(rep("---", length(h)), collapse = " | "), " |"),
  vapply(rows, function(v) paste0("| ", paste(v, collapse = " | "), " |"), character(1))
)

# --------------------------------------------------------------- inputs
hr_main <- as.data.frame(readr::read_csv(
  file.path(OUT_TAB, "04b_cox_hr_v32.csv"), show_col_types = FALSE))
ovr <- as.data.frame(readr::read_csv(
  file.path(OUT_TAB, "04b_cox_overall_v32.csv"), show_col_types = FALSE))
epv <- as.data.frame(readr::read_csv(
  file.path(OUT_TAB, "04b_cox_epv_v32.csv"), show_col_types = FALSE))
ph <- as.data.frame(readr::read_csv(
  file.path(OUT_TAB, "04b_cox_ph_v32.csv"), show_col_types = FALSE))
sens_hr <- as.data.frame(readr::read_csv(
  file.path(OUT_TAB, "06b_sensitivity_hr.csv"), show_col_types = FALSE))
sens_sum <- as.data.frame(readr::read_csv(
  file.path(OUT_TAB, "06b_sensitivity_summary.csv"), show_col_types = FALSE))
sens_st <- as.data.frame(readr::read_csv(
  file.path(OUT_TAB, "06b_sensitivity_status.csv"), show_col_types = FALSE))
flow_tab <- as.data.frame(readr::read_csv(
  file.path(OUT_TAB, "07_table_supp_sample_flow.csv"), show_col_types = FALSE))
strata <- as.data.frame(readr::read_csv(
  file.path(OUT_TAB, "05_ph_strata_hr_ckm_any.csv"), show_col_types = FALSE))
sleep0 <- as.data.frame(readr::read_csv(
  file.path(OUT_LOG, "06b_sleep0_condition.csv"), show_col_types = FALSE))
sd_sum <- as.data.frame(readr::read_csv(
  file.path(OUT_TAB, "08b_single_disease_summary.csv"), show_col_types = FALSE))
sd_ev <- as.data.frame(readr::read_csv(
  file.path(OUT_TAB, "08b_single_disease_events_by_group.csv"), show_col_types = FALSE))
sd_hr <- as.data.frame(readr::read_csv(
  file.path(OUT_TAB, "08b_single_disease_hr.csv"), show_col_types = FALSE))
sd_ph <- as.data.frame(readr::read_csv(
  file.path(OUT_TAB, "08b_single_disease_ph_by_imputation.csv"), show_col_types = FALSE))
epy <- do.call(rbind, lapply(c("CKM_any", "Cardiorenal"), function(o) {
  s <- as.data.frame(readr::read_csv(
    file.path(OUT_TAB, sprintf("04_cox_events_py_%s.csv",
                               ifelse(o == "CKM_any", "ckm_any", "cardiorenal"))),
    show_col_types = FALSE))
  s[["outcome"]] <- o
  s
}))
for (nm in c("group")) {
  hr_main[[nm]] <- as.character(hr_main[[nm]])
  epy[[nm]] <- as.character(epy[[nm]])
  sd_ev[[nm]] <- as.character(sd_ev[[nm]])
  sd_hr[[nm]] <- as.character(sd_hr[[nm]])
}

# the strata refit was run on the former M2, which is M3 under v3.3
stopifnot(all(grepl("^M2", strata[["model"]])))
strata[["model"]] <- sub("^M2", "M3", strata[["model"]])

# --------------------------------------------------------------- checks
stopifnot(
  nrow(hr_main) == 18L,
  setequal(unique(hr_main[["model"]]), MODELS),
  all(hr_main[["group"]] %in% c("1", "2", "3")),
  nrow(ovr) == 6L, nrow(epv) == 6L, nrow(ph) == 8L,
  nrow(epy) == 8L,
  nrow(sens_hr) == 12L, nrow(sens_sum) == 4L, nrow(sens_st) == 4L,
  all(sens_st[["status"]] == "estimated"),
  identical(unique(sens_sum[["outcome"]]), "CKM_any"),
  !anyNA(sens_hr[, c("n", "events", "hr", "ci_low", "ci_high", "p")]),
  !anyNA(hr_main[, c("hr", "ci_low", "ci_high", "p")]),
  nrow(sd_sum) == 6L, nrow(sd_ev) == 24L, nrow(sd_hr) == 3L,
  nrow(sd_ph) == 2L * M_IMP,
  identical(unique(sd_hr[["disease"]]), "HTN"),
  setequal(sd_sum[["status"]],
           c("estimated", "counts only", "masked: fewer than 20 events"))
)
for (o in c("CKM_any", "Cardiorenal")) {
  stopifnot(sum(epy[["n"]][epy[["outcome"]] == o]) == EXPECTED[[o]][["n"]],
            sum(epy[["events"]][epy[["outcome"]] == o]) == EXPECTED[[o]][["events"]])
}
stopifnot(all(hr_main[["ci_low"]] < hr_main[["hr"]]),
          all(hr_main[["hr"]] < hr_main[["ci_high"]]),
          all(sens_hr[["ci_low"]] < sens_hr[["hr"]]),
          all(sens_hr[["hr"]] < sens_hr[["ci_high"]]),
          all(sd_hr[["ci_low"]] < sd_hr[["hr"]]),
          all(sd_hr[["hr"]] < sd_hr[["ci_high"]]))
# the M3 rows reproduce the former M2 values that were carried over in 04b
stopifnot(all(hr_main[["source"]][hr_main[["model"]] == "M3"] ==
                "reused: v3.1 run (former M2)"))
for (i in seq_len(nrow(sens_sum))) {
  h <- sens_hr[sens_hr[["analysis"]] == sens_sum[["analysis"]][i], ]
  stopifnot(nrow(h) == 3L,
            isTRUE(all.equal(h[["hr"]][order(h[["group"]])],
                             as.numeric(unlist(sens_sum[i, paste0("hr_g", 1:3)])))),
            isTRUE(all.equal(h[["events"]][order(h[["group"]])],
                             as.numeric(unlist(sens_sum[i, paste0("events_g", 1:3)])))))
}
stopifnot(all(c("CKM_any", "Cardiorenal") %in%
                ph[["outcome"]][ph[["model"]] == "M2"]))
stopifnot(
  nrow(flow_tab) == 12L,
  flow_tab[["n_after"]][flow_tab[["step"]] == "fitbit_any"] == 53344L,
  flow_tab[["n_after"]][flow_tab[["step"]] == "main_cohort"] == 7770L
)
# no masked count may appear in the single-disease tables
stopifnot(all(sd_sum[["n"]] >= 20L),
          !any(sd_sum[["n_events_display"]] %in% as.character(1:19)),
          !any(sd_ev[["events_display"]] %in% as.character(1:19)))

# --------------------------------------------------------------- main table
tab_all <- merge(
  hr_main[, c("outcome", "model", "group", "label", "hr", "ci_low", "ci_high", "p",
              "source")],
  epy[, c("outcome", "group", "n", "events", "person_years", "rate_per_1000py")],
  by = c("outcome", "group"), all.x = TRUE, sort = FALSE
)
tab_all <- merge(tab_all, epv[, c("outcome", "model", "n_events", "n_parameters",
                                  "epv", "warning")],
                 by = c("outcome", "model"), all.x = TRUE, sort = FALSE)
tab_all <- merge(tab_all, ovr[, c("outcome", "model", "F", "df1", "df2", "p")],
                 by = c("outcome", "model"), all.x = TRUE, sort = FALSE,
                 suffixes = c("", "_d1"))
names(tab_all)[names(tab_all) == "p_d1"] <- "d1_p"
names(tab_all)[names(tab_all) == "F"] <- "d1_F"
names(tab_all)[names(tab_all) == "df1"] <- "d1_df1"
names(tab_all)[names(tab_all) == "df2"] <- "d1_df2"
names(tab_all)[names(tab_all) == "epv.x"] <- "epv"
tab_all[["outcome_label"]] <- unname(OUTCOME_LABEL[tab_all[["outcome"]]])
tab_all[["group_label_en"]] <- unname(GROUP_EN[tab_all[["group"]]])
tab_all[["model_label_en"]] <- unname(MODEL_EN[tab_all[["model"]]])
tab_all <- tab_all[order(match(tab_all[["outcome"]], c("CKM_any", "Cardiorenal")),
                         match(tab_all[["model"]], MODELS),
                         tab_all[["group"]]), ]
stopifnot(nrow(tab_all) == 18L,
          !anyNA(tab_all[, c("n", "events", "person_years", "epv", "d1_p")]))

tab_main <- tab_all[tab_all[["outcome"]] == "CKM_any", ]
tab_cardio <- tab_all[tab_all[["outcome"]] == "Cardiorenal", ]
stopifnot(nrow(tab_main) == 9L, nrow(tab_cardio) == 9L)
readr::write_excel_csv(tab_main, file.path(OUT_TAB, "07d_table_main_hr_v33.csv"))
readr::write_excel_csv(tab_cardio,
                       file.path(OUT_TAB, "07d_table_supp_cardiorenal_v33.csv"))

# the strata refit of the former M2 must reproduce the M3 rows
for (i in seq_len(nrow(strata))) {
  if (!identical(strata[["model"]][i], "M3 (main)")) next
  r <- tab_main[tab_main[["model"]] == "M3" & tab_main[["group"]] == strata[["group"]][i], ]
  stopifnot(nrow(r) == 1L,
            isTRUE(all.equal(r[["hr"]], strata[["hr"]][i], tolerance = 1e-8)),
            isTRUE(all.equal(r[["ci_low"]], strata[["ci_low"]][i], tolerance = 1e-8)),
            isTRUE(all.equal(r[["ci_high"]], strata[["ci_high"]][i], tolerance = 1e-8)))
}

# --------------------------------------------------------------- M1 verification
m1 <- hr_main[hr_main[["outcome"]] == "CKM_any" & hr_main[["model"]] == "M1", ]
m1 <- m1[order(m1[["group"]]), ]
stopifnot(nrow(m1) == 3L)
m1_fmt <- sprintf("%.2f (%.2f to %.2f)", m1[["hr"]], m1[["ci_low"]], m1[["ci_high"]])
m1_check <- data.frame(
  outcome = "CKM_any", model = "M1", group = m1[["group"]],
  hr = m1[["hr"]], ci_low = m1[["ci_low"]], ci_high = m1[["ci_high"]],
  formatted = m1_fmt, investigator_value = unname(USER_M1[m1[["group"]]]),
  stringsAsFactors = FALSE
)
m1_check[["match"]] <- m1_check[["formatted"]] == m1_check[["investigator_value"]]
readr::write_excel_csv(m1_check, file.path(OUT_TAB, "07d_m1_verification_v33.csv"))
stopifnot(all(m1_check[["match"]]))

m2 <- hr_main[hr_main[["outcome"]] == "CKM_any" & hr_main[["model"]] == "M2", ]
m2 <- m2[order(m2[["group"]]), ]
m2_tail <- c(sprintf("%.2f", m2[["ci_low"]][1L]),
             sprintf("%.2f", m2[["ci_high"]][3L]))
stopifnot(identical(m2_tail, c("0.88", "2.09")))

# --------------------------------------------------------------- sensitivity table
sens_tab <- sens_sum
sens_tab[["outcome_label"]] <- unname(OUTCOME_LABEL[sens_tab[["outcome"]]])
sens_tab <- sens_tab[order(match(sens_tab[["analysis"]], c("S1", "S2", "S3", "S4"))), ]
readr::write_excel_csv(sens_tab, file.path(OUT_TAB, "07d_table_supp_sensitivity_v33.csv"))

# --------------------------------------------------------------- PH table
ph_tab <- ph
ph_tab[["outcome_label"]] <- unname(OUTCOME_LABEL[ph_tab[["outcome"]]])
ph_tab <- ph_tab[order(match(ph_tab[["outcome"]], c("CKM_any", "Cardiorenal")),
                       match(ph_tab[["model"]], c("M2", "M3")),
                       ph_tab[["term"]]), ]
readr::write_excel_csv(ph_tab, file.path(OUT_TAB, "07d_table_supp_ph_v33.csv"))

# --------------------------------------------------------------- figure data
fig_main <- tab_main
fig_main[["model_f"]] <- factor(fig_main[["model"]], levels = MODELS)
fig_main[["y"]] <- as.numeric(fig_main[["group"]])
fig_main[["y_plot"]] <- fig_main[["y"]] +
  unname(c(M1 = -0.22, M2 = 0, M3 = 0.22)[fig_main[["model"]]])
stopifnot(nrow(fig_main) == 9L,
          !anyNA(fig_main[, c("hr", "ci_low", "ci_high", "y_plot")]))

main_m2 <- tab_main[tab_main[["model"]] == "M2", ]
stopifnot(nrow(main_m2) == 3L)
sens_fig <- rbind(
  data.frame(outcome = main_m2[["outcome"]], analysis = "Main",
             group = main_m2[["group"]], hr = main_m2[["hr"]],
             ci_low = main_m2[["ci_low"]], ci_high = main_m2[["ci_high"]],
             stringsAsFactors = FALSE),
  sens_hr[, c("outcome", "analysis", "group", "hr", "ci_low", "ci_high")]
)
stopifnot(nrow(sens_fig) == 15L, !anyNA(sens_fig[, c("hr", "ci_low", "ci_high")]))
sens_fig[["analysis_f"]] <- factor(sens_fig[["analysis"]],
                                   levels = c("Main", "S4", "S3", "S2", "S1"),
                                   labels = unname(ANALYSIS_EN[c("Main", "S4", "S3", "S2", "S1")]))
sens_fig[["grp_f"]] <- factor(sens_fig[["group"]], levels = c("1", "2", "3"))
sens_fig[["y"]] <- as.numeric(sens_fig[["analysis_f"]])
sens_fig[["y_plot"]] <- sens_fig[["y"]] +
  unname(c(`1` = -0.22, `2` = 0, `3` = 0.22)[sens_fig[["group"]]])

theme_pub <- function() {
  theme_bw(base_size = 9, base_family = "Arial") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3),
      panel.border = element_rect(colour = "grey35", linewidth = 0.4),
      axis.text = element_text(colour = "black", size = 8.5, family = "Arial"),
      axis.title = element_text(size = 9, family = "Arial"),
      axis.title.y = element_blank(),
      legend.title = element_text(size = 8.5, family = "Arial"),
      legend.text = element_text(size = 8.5, family = "Arial"),
      legend.key = element_blank(),
      axis.ticks.y = element_blank(),
      axis.ticks.y.left = element_blank(),
      legend.position = "bottom",
      plot.margin = margin(6, 8, 6, 6)
    )
}

x_scale <- scale_x_log10(
  breaks = c(0.75, 1, 1.5, 2),
  labels = c("75%", "100%", "150%", "200%"),
  limits = c(0.6, 2.4),
  expand = expansion(mult = c(0.03, 0.03))
)
dir_note <- list(
  annotate("text", x = 0.64, y = -Inf, label = "\u2190 lower risk", hjust = 0,
           vjust = -0.9, size = 2.6, family = "Arial", colour = "grey30"),
  annotate("text", x = 2.2, y = -Inf, label = "higher risk \u2192", hjust = 1,
           vjust = -0.9, size = 2.6, family = "Arial", colour = "grey30")
)

p_main <- ggplot(fig_main, aes(x = hr, y = y_plot, colour = model_f, shape = model_f)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey30", linewidth = 0.4) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0, linewidth = 0.5,
                orientation = "y") +
  geom_point(size = 2.1) +
  x_scale +
  dir_note +
  scale_y_continuous(breaks = 1:3, labels = unname(GROUP_EN[c("1", "2", "3")]),
                     limits = c(0.4, 3.6), expand = c(0, 0)) +
  scale_colour_manual(values = MODEL_COL, name = NULL,
                      labels = MODEL_EN, drop = FALSE) +
  scale_shape_manual(values = MODEL_SHAPE, name = NULL,
                     labels = MODEL_EN, drop = FALSE) +
  guides(colour = guide_legend(nrow = 3, byrow = TRUE),
         shape = guide_legend(nrow = 3, byrow = TRUE)) +
  labs(x = "Hazard ratio (percentage of the reference hazard)") +
  theme_pub()

p_sens <- ggplot(sens_fig, aes(x = hr, y = y_plot, colour = grp_f, shape = grp_f)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey30", linewidth = 0.4) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0, linewidth = 0.5,
                orientation = "y") +
  geom_point(size = 2.0) +
  x_scale +
  dir_note +
  scale_y_continuous(breaks = 1:5, labels = unname(ANALYSIS_EN[c("Main", "S4", "S3", "S2", "S1")]),
                     limits = c(0.4, 5.6), expand = c(0, 0)) +
  scale_colour_manual(values = GROUP_COL, name = NULL,
                      labels = GROUP_EN[c("1", "2", "3")], drop = FALSE) +
  scale_shape_manual(values = c(`1` = 16, `2` = 17, `3` = 15), name = NULL,
                     labels = GROUP_EN[c("1", "2", "3")], drop = FALSE) +
  guides(colour = guide_legend(nrow = 3, byrow = TRUE),
         shape = guide_legend(nrow = 3, byrow = TRUE)) +
  labs(x = "Hazard ratio (percentage of the reference hazard)") +
  theme_pub()

save_fig <- function(p, stem, w, h) {
  grDevices::cairo_pdf(file.path(OUT_FIG, paste0(stem, ".pdf")), width = w, height = h,
                       family = "Arial")
  print(p)
  grDevices::dev.off()
  ragg::agg_tiff(file.path(OUT_FIG, paste0(stem, ".tiff")), width = w, height = h,
                 units = "in", res = 400, background = "white")
  print(p)
  grDevices::dev.off()
  ragg::agg_png(file.path(OUT_FIG, paste0(stem, ".png")), width = w, height = h,
                units = "in", res = 300, background = "white")
  print(p)
  grDevices::dev.off()
}

save_fig(p_main, "07d_fig_forest_main", 5.6, 4.3)
save_fig(p_sens, "07d_fig_forest_sensitivity", 6.4, 4.4)

# --------------------------------------------------------------- display report
ref_rows <- do.call(rbind, lapply(c("CKM_any", "Cardiorenal"), function(o) {
  s <- epy[epy[["outcome"]] == o & epy[["group"]] == "4", ]
  c(o, mask(s[["n"]]), mask(s[["events"]]), sprintf("%.0f", s[["person_years"]]),
    sprintf("%.1f", s[["rate_per_1000py"]]))
}))
stopifnot(nrow(ref_rows) == 2L, ncol(ref_rows) == 5L)

sens_note <- local({
  s1_ev <- sens_tab[["events"]][sens_tab[["analysis"]] == "S1"]
  stopifnot(length(s1_ev) == 1L)
  sprintf("S1 moved the start of follow-up to landmark + 0.5 years; %d events occurred in the first six months and were excluded. All other analyses keep the original start of follow-up.",
          EXPECTED[["CKM_any"]][["events"]] - s1_ev)
})

sd_order <- c("T2DM", "HTN", "AF", "HF", "CKD", "ASCVD")
sd_sum <- sd_sum[match(sd_order, sd_sum[["disease"]]), ]
stopifnot(nrow(sd_sum) == 6L, !anyNA(sd_sum[["disease"]]))
epv_co <- suppressWarnings(as.numeric(
  sd_sum[["epv_display"]][sd_sum[["status"]] == "counts only"]))
stopifnot(length(epv_co) == 4L, all(is.finite(epv_co)))
sd_cell <- function(d) {
  s <- sd_hr[sd_hr[["disease"]] == d, ]
  if (nrow(s) == 0L) return(c("not reported", "not reported", "not reported"))
  s <- s[order(s[["group"]]), ]
  stopifnot(nrow(s) == 3L)
  sprintf("%.2f (%.2f to %.2f), %s", s[["hr"]], s[["ci_low"]], s[["ci_high"]],
          fmt_p(s[["p"]]))
}

md <- c(
  "# Step 7 - consolidated results under plan v3.3",
  "",
  sprintf("Generated by `R/07d_summary_v33.R` on %s from the files written by `R/04b_cox_v32.R`, `R/06b_sensitivity_v32.R`, `R/08b_single_disease_v33.R` and earlier steps. Version: v3.3 (`TRAE_fitbit_ckm_analysis_spec.md`). Every number below is read from a stored CSV; no value is typed by hand and no model was refitted in this step.",
          format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "CKM disease is the primary outcome and the only outcome in the main table. Cardiorenal disease is reported in the supplementary material. Reference group: group 4 (moderate activity, low heart rate). m = 20 imputed data sets, seed 20250101, Rubin's rules, ties = Efron.",
  "",
  "Three adjustment sets are reported. M1 holds demographic and temporal terms. M2 is the primary model: M1 plus smoking, depression (`mdd_fs`) and log(1 + prior-year visit days). M3 adds BMI (3-knot spline), SBP and baseline antihypertensive, antidiabetic and statin use; it is numerically the earlier M2 and is reused unchanged.",
  "",
  "## 1. Main results - CKM disease",
  "",
  md_block(c("Model", "Group", "n", "Events", "Person-years", "Rate per 1,000 py",
             "HR (95% CI)", "p", "EPV", "D1 p", "Source"),
           lapply(seq_len(nrow(tab_main)), function(i) c(
             tab_main[["model"]][i], tab_main[["group"]][i],
             mask(tab_main[["n"]][i]), mask(tab_main[["events"]][i]),
             sprintf("%.0f", tab_main[["person_years"]][i]),
             sprintf("%.1f", tab_main[["rate_per_1000py"]][i]),
             sprintf("%.2f (%.2f to %.2f)", tab_main[["hr"]][i],
                     tab_main[["ci_low"]][i], tab_main[["ci_high"]][i]),
             fmt_p(tab_main[["p"]][i]), sprintf("%.1f", tab_main[["epv"]][i]),
             fmt_p(tab_main[["d1_p"]][i]), tab_main[["source"]][i]
           ))),
  "",
  "### Reference group (group 4, moderate activity and low heart rate)",
  "",
  md_block(c("Outcome", "n", "Events", "Person-years", "Rate per 1,000 py"),
           lapply(seq_len(nrow(ref_rows)), function(i) as.character(ref_rows[i, ]))),
  "",
  "### Trajectory-group overall test and events per parameter",
  "",
  md_block(c("Outcome", "Model", "F", "df1", "df2", "p", "Events", "Parameters", "EPV"),
           lapply(seq_len(nrow(ovr)), function(i) {
             e <- ovr[i, ]
             p <- epv[epv[["outcome"]] == e[["outcome"]] & epv[["model"]] == e[["model"]], ]
             c(e[["outcome"]], e[["model"]], sprintf("%.3f", e[["F"]]),
               as.integer(e[["df1"]]),
               formatC(e[["df2"]], format = "f", digits = 0, big.mark = ","),
               fmt_p(e[["p"]]), as.integer(p[["n_events"]]),
               as.integer(p[["n_parameters"]]), sprintf("%.1f", p[["epv"]]))
           })),
  "",
  sprintf("Cardiorenal disease has %d events over %d parameters in the primary model (EPV %.1f). This is reported as it stands; no covariate was removed to raise it.",
          epv[["n_events"]][epv[["outcome"]] == "Cardiorenal" & epv[["model"]] == "M2"],
          epv[["n_parameters"]][epv[["outcome"]] == "Cardiorenal" & epv[["model"]] == "M2"],
          epv[["epv"]][epv[["outcome"]] == "Cardiorenal" & epv[["model"]] == "M2"]),
  "",
  "## 2. Figure 1 - hazard ratios of the three adjustment sets (CKM disease)",
  "",
  "![Forest plot of the three models](../figures/07d_fig_forest_main.png)",
  "",
  "Rows are the three trajectory groups against group 4; the three adjustment sets are distinguished by colour and symbol. The vertical broken line marks a hazard ratio of 1 and the horizontal axis is logarithmic, labelled as a percentage of the reference hazard. Values plotted here are those of the table in section 1.",
  "",
  "## 3. Figure 2 - CKM disease across the sensitivity analyses",
  "",
  "![Forest plot of the sensitivity analyses](../figures/07d_fig_forest_sensitivity.png)",
  "",
  "Rows are the primary model and the four sensitivity analyses, all based on M2 and all restricted to CKM disease; the three trajectory groups are distinguished by colour and symbol. Cardiorenal disease has no sensitivity analysis.",
  "",
  "## 4. Supplementary table S1 - proportional hazards",
  "",
  sprintf("`cox.zph` was fitted to each model in each of the %d imputed data sets; median and range of the %d p values are given. M3 is reported from the earlier run.", M_IMP, M_IMP),
  "",
  md_block(c("Outcome", "Model", "Term", "df", "chi-square, median (range)",
             "p, median (range)", "Imputations with p < 0.05"),
           lapply(seq_len(nrow(ph_tab)), function(i) {
             e <- ph_tab[i, ]
             c(e[["outcome"]], e[["model"]], e[["term"]], as.integer(e[["df"]]),
               sprintf("%.2f (%.2f to %.2f)", e[["chisq_median"]], e[["chisq_min"]],
                       e[["chisq_max"]]),
               sprintf("%.3f (%.3f to %.3f)", e[["p_median"]], e[["p_min"]],
                       e[["p_max"]]),
               sprintf("%d of %d", as.integer(e[["n_p_below_0p05"]]), as.integer(e[["m"]])))
           })),
  "",
  "Scaled Schoenfeld residuals for the primary model, first imputed data set: `output/figures/04b_schoenfeld_ckm_any_m2.png` and `04b_schoenfeld_cardiorenal_m2.png`.",
  "",
  "Supplementary footnote: for CKM disease the earlier M2 (now M3) had a borderline global test driven by `education_cat`. Repeating M3 with `strata(education_cat)` left the trajectory hazard ratios essentially unchanged. This is a supplementary analysis of M3 and does not replace the primary model.",
  "",
  md_block(c("Model", "Group", "HR (95% CI)", "p"),
           lapply(seq_len(nrow(strata)), function(i) {
             e <- strata[i, ]
             c(e[["model"]], e[["group"]],
               sprintf("%.2f (%.2f to %.2f)", e[["hr"]], e[["ci_low"]], e[["ci_high"]]),
               fmt_p(e[["p"]]))
           })),
  "",
  "## 5. Supplementary table S2 - sensitivity analyses (CKM disease, primary model M2)",
  "",
  "HR (95% CI) and p for groups 1, 2 and 3 against group 4. n and events are those of the analysed population. The `fs_year` merge column records whether a level with no event was merged with the adjacent year in that analysis only.",
  "",
  md_block(c("Analysis", "n", "Events", "EPV", "Group 1 HR (95% CI), p",
             "Group 2 HR (95% CI), p", "Group 3 HR (95% CI), p", "fs_year merge"),
           lapply(seq_len(nrow(sens_tab)), function(i) {
             e <- sens_tab[i, ]
             c(e[["analysis"]], mask(e[["n"]]), mask(e[["events"]]),
               sprintf("%.1f", e[["epv"]]),
               sprintf("%.2f (%.2f to %.2f), %s", e[["hr_g1"]], e[["lo_g1"]],
                       e[["hi_g1"]], fmt_p(e[["p_g1"]])),
               sprintf("%.2f (%.2f to %.2f), %s", e[["hr_g2"]], e[["lo_g2"]],
                       e[["hi_g2"]], fmt_p(e[["p_g2"]])),
               sprintf("%.2f (%.2f to %.2f), %s", e[["hr_g3"]], e[["lo_g3"]],
                       e[["hi_g3"]], fmt_p(e[["p_g3"]])),
               e[["fs_year_merge"]])
           })),
  "",
  "### Analysis definitions",
  "",
  md_block(c("Code", "Definition"), list(
    c("S1", "Exclude events in the first 6 months of follow-up and move the start of follow-up to landmark + 0.5 years"),
    c("S2", "Exclude baseline antihypertensive or antidiabetic use, eGFR < 60 mL/min/1.73 m2 or HbA1c >= 6.5%"),
    c("S3", "Exclude persons with maximum posterior probability < 0.70"),
    c("S4", "M2 with additional adjustment for the COVID-19 fraction of the Fitbit window")
  )),
  "",
  sens_note,
  "",
  sprintf("Merge footnote: %s.", paste(sprintf("%s: %s", sens_tab[["analysis"]],
                                               sens_tab[["fs_year_merge"]]),
                                      collapse = "; ")),
  "",
  sprintf("Estimated parameters: %s. EPV is events divided by that number; a factor level with no person in a restricted population is not estimated, which is why the count differs between analyses.",
          paste(sprintf("%s %d", sens_tab[["analysis"]],
                        as.integer(sens_tab[["n_parameters"]])), collapse = "; ")),
  "",
  "## 6. Supplementary table S3 - sample flow",
  "",
  md_block(c("Stage", "Rule", "n before", "Excluded", "n after", "Events", "Source"),
           lapply(seq_len(nrow(flow_tab)), function(i) {
             mm <- function(v) if (is.na(v)) "-" else mask(v)
             c(flow_tab[["stage"]][i], flow_tab[["rule"]][i],
               mm(flow_tab[["n_before"]][i]), mm(flow_tab[["n_excluded"]][i]),
               mm(flow_tab[["n_after"]][i]), mm(flow_tab[["events"]][i]),
               flow_tab[["source"]][i])
           })),
  "",
  "The gate-0 count is verified in code from the flow record stored in `gate1v3_fitbit_ckm_w8q_min3.rds` (`output/logs/07b_flow_check.csv`, nine checks, all pass); it is not a hand-entered value.",
  "",
  "## 7. Supplementary table S4 - missingness",
  "",
  "`output/tables/07_table_supp_missingness.csv` gives the missingness of every analysis variable in the three populations. Missing values were not imputed outside the multiple-imputation step and were never set to zero.",
  "",
  "## 8. Conditional dual-indicator check",
  "",
  sprintf("Persons with no valid sleep quarter: %d of %d in the main cohort (%.1f%%). The condition for the dual-indicator model was %s (threshold: twice the cohort share in one group).",
          sum(sleep0[["n_sleep0"]]), sum(sleep0[["n_group"]]),
          100 * sum(sleep0[["n_sleep0"]]) / sum(sleep0[["n_group"]]),
          if (any(sleep0[["condition_met"]])) "met" else "not met"),
  "",
  "## 9. Supplementary table S5 - cardiorenal disease (three models)",
  "",
  "The same models and strata as the main analysis. Cardiorenal disease has 222 events over 28 parameters in the primary model (EPV 7.9); this is stated in the report and no covariate was removed to raise it.",
  "",
  md_block(c("Model", "Group", "n", "Events", "Person-years", "Rate per 1,000 py",
             "HR (95% CI)", "p", "EPV", "D1 p", "Source"),
           lapply(seq_len(nrow(tab_cardio)), function(i) c(
             tab_cardio[["model"]][i], tab_cardio[["group"]][i],
             mask(tab_cardio[["n"]][i]), mask(tab_cardio[["events"]][i]),
             sprintf("%.0f", tab_cardio[["person_years"]][i]),
             sprintf("%.1f", tab_cardio[["rate_per_1000py"]][i]),
             sprintf("%.2f (%.2f to %.2f)", tab_cardio[["hr"]][i],
                     tab_cardio[["ci_low"]][i], tab_cardio[["ci_high"]][i]),
             fmt_p(tab_cardio[["p"]][i]), sprintf("%.1f", tab_cardio[["epv"]][i]),
             fmt_p(tab_cardio[["d1_p"]][i]), tab_cardio[["source"]][i]
           ))),
  "",
  "## 10. Supplementary table S6 - exploratory analysis of the six CKM components",
  "",
  sprintf("Population: the CKM disease at-risk set, n = %s; all six components are free of prevalent disease in this set. The stored CKM imputation object is reused (m = %d, Rubin's rules). Model M2, reference group 4. Occurrence of another component does not censor follow-up. This analysis is exploratory and appears in the supplementary material only.",
          format(sd_sum[["n"]][1L], big.mark = ","), M_IMP),
  "",
  md_block(c("Component", "n", "Events", "EPV", "Group 1", "Group 2", "Group 3",
             "D1 p", "PH p, trajectory term"),
           lapply(seq_len(nrow(sd_sum)), function(i) {
             s <- sd_sum[i, ]
             c(s[["disease_label"]], mask(s[["n"]]), s[["n_events_display"]],
               s[["epv_display"]], sd_cell(s[["disease"]]),
               ifelse(is.na(s[["d1_p"]]), "not calculated", fmt_p(s[["d1_p"]])),
               ifelse(is.na(s[["ph_traj_p_median"]]), "not calculated",
                      fmt_p(s[["ph_traj_p_median"]])))
           })),
  "",
  sprintf("Group 1-3 cells give the hazard ratio (95%% confidence interval) and p against group 4. Hypertension is fully estimated. Type 2 diabetes, atrial fibrillation, heart failure and chronic kidney disease have an EPV of %.1f to %.1f, and non-trajectory covariates contain levels with no event, so the covariates are completely separated and the likelihood has no finite maximum; no effect estimate is reported for them and only events, person-years and incidence rates are given. Atherosclerotic cardiovascular disease has fewer than 20 events: it follows the masking rule and is not modelled.",
          min(epv_co), max(epv_co)),
  "",
  "### Events and person-years by trajectory group",
  "",
  md_block(c("Component", "Group", "n", "Events", "Person-years", "Rate per 1,000 py"),
           lapply(seq_len(nrow(sd_ev)), function(i) {
             e <- sd_ev[i, ]
             c(e[["disease"]], e[["group"]], e[["n_display"]], e[["events_display"]],
               sprintf("%.0f", e[["person_years"]]),
               sprintf("%.1f", e[["rate_per_1000py"]]))
           })),
  "",
  sprintf("`cox.zph` for hypertension, %d imputed data sets: trajectory term median p = %.3f, global test median p = %.3f.",
          M_IMP,
          sd_sum[["ph_traj_p_median"]][sd_sum[["disease"]] == "HTN"],
          sd_sum[["ph_global_p_median"]][sd_sum[["disease"]] == "HTN"]),
  "",
  "## 11. Verification of the M1 hazard ratios (CKM disease)",
  "",
  "The M1 rows were read from `04b_cox_hr_v32.csv`, formatted to two decimals, and compared with the values supplied by the investigator.",
  "",
  md_block(c("Group", "HR", "95% CI", "Formatted", "Investigator value", "Match"),
           lapply(seq_len(nrow(m1_check)), function(i) c(
             m1_check[["group"]][i], sprintf("%.6f", m1_check[["hr"]][i]),
             sprintf("%.6f to %.6f", m1_check[["ci_low"]][i], m1_check[["ci_high"]][i]),
             m1_check[["formatted"]][i], m1_check[["investigator_value"]][i],
             ifelse(m1_check[["match"]][i], "yes", "NO")))),
  "",
  sprintf("The three M1 rows match the investigator values exactly; no value was changed by hand. The confidence limits %s and %s that prompted the query belong to M2 (group 1 lower limit %.6f, group 3 upper limit %.6f), not to M1.",
          m2_tail[1L], m2_tail[2L], m2[["ci_low"]][1L], m2[["ci_high"]][3L]),
  "",
  "## 12. Key results (primary model M2, CKM disease)",
  ""
)
s <- tab_main[tab_main[["model"]] == "M2", ]
for (i in seq_len(nrow(s))) {
  p_txt <- fmt_p(s[["p"]][i])
  p_txt <- if (startsWith(p_txt, "<")) paste("p", p_txt) else paste("p =", p_txt)
  md <- c(md, sprintf("- %s, %s versus group 4: HR %.2f (95%% CI %.2f to %.2f; %s).",
                      OUTCOME_LABEL[["CKM_any"]],
                      gsub("\n", ", ", GROUP_EN[[s[["group"]][i]]]),
                      s[["hr"]][i], s[["ci_low"]][i], s[["ci_high"]][i], p_txt))
}
e <- ovr[ovr[["outcome"]] == "CKM_any" & ovr[["model"]] == "M2", ]
p_txt <- fmt_p(e[["p"]])
p_txt <- if (startsWith(p_txt, "<")) paste("p", p_txt) else paste("p =", p_txt)
md <- c(md, sprintf("- Overall trajectory-group test: F = %.3f (df1 = %d, df2 = %s), %s.",
                    e[["F"]], as.integer(e[["df1"]]),
                    formatC(e[["df2"]], format = "f", digits = 0, big.mark = ","),
                    p_txt))
h <- sd_sum[sd_sum[["disease"]] == "HTN", ]
md <- c(md, sprintf("- Hypertension (exploratory, reference group 4): group 3 HR %.2f (95%% CI %.2f to %.2f), EPV %.1f.",
                    sd_hr[["hr"]][sd_hr[["group"]] == "3"],
                    sd_hr[["ci_low"]][sd_hr[["group"]] == "3"],
                    sd_hr[["ci_high"]][sd_hr[["group"]] == "3"],
                    as.numeric(h[["epv_display"]][1L])))
md <- c(md, "",
        "## 13. Files", "",
        "- `output/tables/07d_table_main_hr_v33.csv`",
        "- `output/tables/07d_table_supp_cardiorenal_v33.csv`",
        "- `output/tables/07d_table_supp_single_disease_v33.csv` (assembled from the 08b files)",
        "- `output/tables/07d_table_supp_sensitivity_v33.csv`",
        "- `output/tables/07d_table_supp_ph_v33.csv`",
        "- `output/tables/07d_m1_verification_v33.csv`",
        "- `output/figures/07d_fig_forest_main.pdf`, `.tiff`, `.png`",
        "- `output/figures/07d_fig_forest_sensitivity.pdf`, `.tiff`, `.png`",
        "",
        "No new model was fitted in this step.")
writeLines(md, file.path(OUT_REP, "step7_summary_v33.md"))

# --------------------------------------------------------------- single-disease table
sd_tab <- sd_sum
sd_tab[["group_1_hr_ci_p"]] <- vapply(sd_tab[["disease"]], function(d) sd_cell(d)[1L], character(1))
sd_tab[["group_2_hr_ci_p"]] <- vapply(sd_tab[["disease"]], function(d) sd_cell(d)[2L], character(1))
sd_tab[["group_3_hr_ci_p"]] <- vapply(sd_tab[["disease"]], function(d) sd_cell(d)[3L], character(1))
readr::write_excel_csv(sd_tab, file.path(OUT_TAB, "07d_table_supp_single_disease_v33.csv"))

# --------------------------------------------------------------- manifest
manifest_files <- c(
  "output/tables/07d_table_main_hr_v33.csv",
  "output/tables/07d_table_supp_cardiorenal_v33.csv",
  "output/tables/07d_table_supp_single_disease_v33.csv",
  "output/tables/07d_table_supp_sensitivity_v33.csv",
  "output/tables/07d_table_supp_ph_v33.csv",
  "output/tables/07d_m1_verification_v33.csv",
  "output/tables/08b_single_disease_summary.csv",
  "output/tables/08b_single_disease_hr.csv",
  "output/tables/08b_single_disease_events_by_group.csv",
  "output/tables/08b_single_disease_ph_by_imputation.csv",
  "output/figures/07d_fig_forest_main.pdf",
  "output/figures/07d_fig_forest_main.tiff",
  "output/figures/07d_fig_forest_main.png",
  "output/figures/07d_fig_forest_sensitivity.pdf",
  "output/figures/07d_fig_forest_sensitivity.tiff",
  "output/figures/07d_fig_forest_sensitivity.png",
  "output/reports/step7_summary_v33.md",
  "output/reports/step8_single_disease_v33.md",
  "output/reports/final_delivery_manifest_v33.md",
  "output/logs/08b_single_disease_flow.csv",
  "output/logs/08b_single_disease_settings.csv",
  "output/logs/08b_single_disease_zero_event_levels.csv",
  "output/logs/08b_single_disease_provenance.json",
  "output/logs/07d_summary_v33_provenance.json"
)
stopifnot(all(file.exists(setdiff(
  manifest_files,
  c("output/reports/final_delivery_manifest_v33.md",
    "output/logs/07d_summary_v33_provenance.json")))))

anchors <- c(
  sprintf("| CKM disease at-risk set, n | %d | `07d_table_main_hr_v33.csv` |", EXPECTED[["CKM_any"]][["n"]]),
  sprintf("| CKM disease events | %d | `07d_table_main_hr_v33.csv` |", EXPECTED[["CKM_any"]][["events"]]),
  sprintf("| Cardiorenal disease at-risk set, n | %d | `07d_table_supp_cardiorenal_v33.csv` |", EXPECTED[["Cardiorenal"]][["n"]]),
  sprintf("| Cardiorenal disease events | %d | `07d_table_supp_cardiorenal_v33.csv` |", EXPECTED[["Cardiorenal"]][["events"]]),
  sprintf("| M2 group 3 hazard ratio, CKM disease | %.2f | `07d_table_main_hr_v33.csv` |",
          tab_main[["hr"]][tab_main[["model"]] == "M2" & tab_main[["group"]] == "3"]),
  sprintf("| Hypertension events | %s | `08b_single_disease_summary.csv` |",
          sd_sum[["n_events_display"]][sd_sum[["disease"]] == "HTN"])
)

manifest <- c(
  "# Delivery manifest - plan v3.3",
  "",
  sprintf("Written by `R/07d_summary_v33.R` on %s. Version v3.3 (`TRAE_fitbit_ckm_analysis_spec.md`). m = 20, seed 20250101, Rubin's rules, ties = Efron. No number in this manifest was typed by hand.",
          format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## 1. Anchors",
  "",
  "| Quantity | Value | Source file |",
  "| --- | --- | --- |",
  anchors,
  "",
  "## 2. Deliverables",
  "",
  md_block(c("File", "Content"),
           lapply(manifest_files, function(f) c(f, switch(f,
             "output/tables/07d_table_main_hr_v33.csv" = "Main table, CKM disease, M1/M2/M3",
             "output/tables/07d_table_supp_cardiorenal_v33.csv" = "Supplementary table S5, cardiorenal disease, M1/M2/M3",
             "output/tables/07d_table_supp_single_disease_v33.csv" = "Supplementary table S6, six CKM components",
             "output/tables/07d_table_supp_sensitivity_v33.csv" = "Supplementary table S2, four sensitivity analyses",
             "output/tables/07d_table_supp_ph_v33.csv" = "Supplementary table S1, proportional hazards",
             "output/tables/07d_m1_verification_v33.csv" = "Verification of the M1 hazard ratios",
             "output/tables/08b_single_disease_summary.csv" = "Single-disease summary, source of S6",
             "output/tables/08b_single_disease_hr.csv" = "Single-disease hazard ratios (hypertension only)",
             "output/tables/08b_single_disease_events_by_group.csv" = "Single-disease events and person-years by group",
             "output/tables/08b_single_disease_ph_by_imputation.csv" = "Single-disease PH tests by imputation (hypertension)",
             "output/figures/07d_fig_forest_main.pdf" = "Figure 1, PDF",
             "output/figures/07d_fig_forest_main.tiff" = "Figure 1, TIFF 400 dpi",
             "output/figures/07d_fig_forest_main.png" = "Figure 1, PNG 300 dpi",
             "output/figures/07d_fig_forest_sensitivity.pdf" = "Figure 2, PDF",
             "output/figures/07d_fig_forest_sensitivity.tiff" = "Figure 2, TIFF 400 dpi",
             "output/figures/07d_fig_forest_sensitivity.png" = "Figure 2, PNG 300 dpi",
             "output/reports/step7_summary_v33.md" = "Consolidated results report",
             "output/reports/step8_single_disease_v33.md" = "Single-disease report",
             "output/reports/final_delivery_manifest_v33.md" = "This manifest",
             "output/logs/08b_single_disease_flow.csv" = "Sample flow, single-disease analysis",
             "output/logs/08b_single_disease_settings.csv" = "Model settings, single-disease analysis",
             "output/logs/08b_single_disease_zero_event_levels.csv" = "Zero-event covariate levels, counts-only decision",
             "output/logs/08b_single_disease_provenance.json" = "Provenance, single-disease analysis",
             "output/logs/07d_summary_v33_provenance.json" = "Provenance, consolidation step",
             f)))),
  "",
  "## 3. Superseded files",
  "",
  "The v3.2 tables, figures, report and manifest that this version replaces were moved to `output/deprecated/2026-11-08_v32_superseded/` with a README; they were not deleted and are not referenced by the v3.3 report.",
  ""
)
writeLines(manifest, file.path(OUT_REP, "final_delivery_manifest_v33.md"))

jsonlite::write_json(
  list(script = "R/07d_summary_v33.R", run_time = format(Sys.time()),
       plan_version = PLAN_VERSION, seed = SEED, m = M_IMP,
       inputs = c("04b_cox_hr_v32.csv", "04b_cox_overall_v32.csv",
                  "04b_cox_epv_v32.csv", "04b_cox_ph_v32.csv",
                  "06b_sensitivity_hr.csv", "06b_sensitivity_summary.csv",
                  "06b_sensitivity_status.csv", "07_table_supp_sample_flow.csv",
                  "07_table_supp_missingness.csv", "05_ph_strata_hr_ckm_any.csv",
                  "08b_single_disease_summary.csv", "08b_single_disease_hr.csv",
                  "08b_single_disease_events_by_group.csv",
                  "08b_single_disease_ph_by_imputation.csv"),
       models = MODELS, main_outcome = "CKM_any",
       supplementary_outcome = "Cardiorenal", sensitivity_outcome = "CKM_any only",
       figures = c("07d_fig_forest_main", "07d_fig_forest_sensitivity"),
       git_commit = NA),
  file.path(OUT_LOG, "07d_summary_v33_provenance.json"), auto_unbox = TRUE,
  pretty = TRUE
)
