# =============================================================================
# Step 11a: sample flow diagram (Figure 1), final version
#   All counts are read from output/tables/07_table_supp_sample_flow.csv; no
#   count is typed into this script. The seven population boxes follow the main
#   chain from 53,344 to 7,770. The "prior EHR" (4,620) and "alive and adult"
#   (9) exclusions are merged into one box (4,629), so no count below 20 appears.
#   At the foot of the diagram the chain branches into the two outcome cohorts:
#   CKM at-risk (4,648; 3,122 excluded at baseline) and cardiorenal at-risk
#   (6,871; 899 excluded at baseline), also read from the same flow table.
#   Box text is written in plain language: internal step names such as gate0,
#   "main cohort" or "analytic file" never appear in the figure.
#   174 mm wide, Arial, PDF with embedded fonts, 600 dpi TIFF.
#   No model is fitted and no table is changed. No console output.
# =============================================================================
suppressPackageStartupMessages({
  library(readr)
  library(jsonlite)
})

# All paths are relative to a single root. Override the root with the environment
# variable FIG_ROOT, or edit the default below; the layout under the root is
# output/tables (inputs), output/figures, output/logs and output/reports.
ROOT <- Sys.getenv("FIG_ROOT", unset = ".")
FLOW_FILE <- file.path(ROOT, "output/tables/07_table_supp_sample_flow.csv")
OUT_FIG <- file.path(ROOT, "output/figures")
OUT_TAB <- file.path(ROOT, "output/tables")
OUT_LOG <- file.path(ROOT, "output/logs")
OUT_REP <- file.path(ROOT, "output/reports")
FINAL <- file.path(OUT_FIG, "final")
invisible(lapply(c(OUT_FIG, OUT_TAB, OUT_LOG, OUT_REP, FINAL), dir.create,
                 recursive = TRUE, showWarnings = FALSE))

W_MM <- 174
H_MM <- 240
EXPECTED_FIRST <- 53344L
EXPECTED_LAST <- 7770L
EXPECTED_MERGED_EXCL <- 4629L
EXPECTED_CKM_N <- 4648L
EXPECTED_CKM_EXCL <- 3122L
EXPECTED_CR_N <- 6871L
EXPECTED_CR_EXCL <- 899L
MAIN_CHAIN <- c("fitbit_any", "valid_periods", "has_ehr", "prior_ehr",
                "alive_adult", "landmark_before_cdr_end", "has_fu_ehr",
                "main_cohort")

f_all <- as.data.frame(readr::read_csv(FLOW_FILE, show_col_types = FALSE))
stopifnot(nrow(f_all) == 12L, all(MAIN_CHAIN %in% f_all[["step"]]),
          !anyDuplicated(f_all[["step"]]))
f <- f_all[match(MAIN_CHAIN, f_all[["step"]]), ]
stopifnot(nrow(f) == 8L)

# gate0 row has no n_before in the stored table; every later row must balance
stopifnot(is.na(f[["n_before"]][1L]), is.na(f[["n_excluded"]][1L]),
          f[["n_after"]][1L] == EXPECTED_FIRST)
for (i in 2:8) {
  stopifnot(!is.na(f[["n_before"]][i]), !is.na(f[["n_excluded"]][i]),
            f[["n_before"]][i] - f[["n_excluded"]][i] == f[["n_after"]][i])
}
stopifnot(f[["n_after"]][8L] == EXPECTED_LAST,
          all(f[["n_after"]][-8L] == f[["n_before"]][-1L]))

# merge the "prior EHR" and "alive and adult" steps into one exclusion box
stopifnot(f[["n_excluded"]][4L] + f[["n_excluded"]][5L] == EXPECTED_MERGED_EXCL,
          f[["n_after"]][5L] == f[["n_before"]][6L])

# ------------------------------------------------------- plain-language labels
# Display wording only; every number still comes from the flow table above.
POP_DISPLAY <- c(
  fitbit_any = "Fitbit participants",
  valid_periods = "At least three valid quarters of Fitbit data",
  has_ehr = "At least one healthcare record",
  prior_ehr = "Alive and adult, with a healthcare record before follow-up started",
  landmark_before_cdr_end = "Follow-up start within the data coverage window",
  has_fu_ehr = "At least one healthcare visit during follow-up",
  main_cohort = "Analysis cohort"
)
EXC_DISPLAY <- c(
  valid_periods = "Fewer than three valid quarters of Fitbit data",
  has_ehr = "No healthcare record",
  prior_ehr = paste("No healthcare record before follow-up started,",
                    "or not alive and adult at that date"),
  landmark_before_cdr_end = "Follow-up start after the end of data coverage",
  has_fu_ehr = "No healthcare visit during follow-up",
  main_cohort = paste("Less than one year of healthcare records",
                      "before follow-up start")
)
stopifnot(all(MAIN_CHAIN[c(1:4, 6:8)] %in% names(POP_DISPLAY)),
          all(MAIN_CHAIN[c(2:4, 6:8)] %in% names(EXC_DISPLAY)))

pop_n <- f[["n_after"]][c(1:3, 5:8)]           # step 4 population absorbed
pop_step <- MAIN_CHAIN[c(1:3, 5:8)]
pop_text <- unname(POP_DISPLAY[c("fitbit_any", "valid_periods", "has_ehr",
                                 "prior_ehr", "landmark_before_cdr_end",
                                 "has_fu_ehr", "main_cohort")])
exc_n <- c(f[["n_excluded"]][2:3], EXPECTED_MERGED_EXCL, f[["n_excluded"]][6:8])
exc_step <- MAIN_CHAIN[c(2:4, 6:8)]
exc_text <- unname(EXC_DISPLAY[c("valid_periods", "has_ehr", "prior_ehr",
                                 "landmark_before_cdr_end", "has_fu_ehr",
                                 "main_cohort")])
stopifnot(length(pop_n) == 7L, length(pop_text) == 7L,
          length(exc_n) == 6L, length(exc_text) == 6L)

# ------------------------------------------------------------- outcome cohorts
oc <- f_all[f_all[["step"]] %in% c("CKM_any_atrisk", "Cardiorenal_atrisk"), ]
stopifnot(nrow(oc) == 2L)
ckm <- oc[oc[["step"]] == "CKM_any_atrisk", ]
cr <- oc[oc[["step"]] == "Cardiorenal_atrisk", ]
stopifnot(nrow(ckm) == 1L, nrow(cr) == 1L,
          ckm[["n_after"]] == EXPECTED_CKM_N,
          ckm[["n_excluded"]] == EXPECTED_CKM_EXCL,
          cr[["n_after"]] == EXPECTED_CR_N,
          cr[["n_excluded"]] == EXPECTED_CR_EXCL,
          ckm[["n_before"]] == EXPECTED_LAST, cr[["n_before"]] == EXPECTED_LAST)

out_title <- c("CKM conditions: at-risk cohort",
               "Cardiorenal disease: at-risk cohort")
out_n <- c(as.integer(ckm[["n_after"]]), as.integer(cr[["n_after"]]))
out_exc <- c(as.integer(ckm[["n_excluded"]]), as.integer(cr[["n_excluded"]]))

# every number that will be printed must be at least 20
all_counts <- c(pop_n, exc_n, out_n, out_exc)
stopifnot(all(all_counts >= 20L))

boxes <- data.frame(
  type = c(rep("population", 7L), rep("exclusion", 6L), rep("outcome", 2L)),
  order = c(1:7, 2:7, 1:2),
  step = c(pop_step, exc_step, c("CKM_any_atrisk", "Cardiorenal_atrisk")),
  rule = c(pop_text, exc_text, out_title),
  n = c(pop_n, exc_n, out_n),
  n_excluded = c(rep(NA_integer_, 7L), exc_n, out_exc),
  stringsAsFactors = FALSE
)
readr::write_excel_csv(boxes, file.path(OUT_TAB, "11a_sample_flow_boxes.csv"))

# ---------------------------------------------------------------- draw figure
b <- as.data.frame(readr::read_csv(file.path(OUT_TAB, "11a_sample_flow_boxes.csv"),
                                   show_col_types = FALSE))
pop <- b[b[["type"]] == "population", ]
exc <- b[b[["type"]] == "exclusion", ]
out <- b[b[["type"]] == "outcome", ]
stopifnot(nrow(pop) == 7L, nrow(exc) == 6L, nrow(out) == 2L)

BOX_H <- 0.052
TOP <- 0.985
LAST_BOTTOM <- 0.300
GAP <- (TOP - LAST_BOTTOM - 7 * BOX_H) / 6
XL <- c(0.030, 0.470)
XR <- c(0.540, 0.970)
MID_L <- mean(XL)
MID_R <- mean(XR)
BUS <- 0.245
OUT_TOP <- 0.190
OUT_BOT <- 0.075
COL_INK <- "grey15"
COL_BORDER <- "#0072B2"
COL_FILL <- "#DDE9F4"
COL_EXC_BORDER <- "#999999"
COL_EXC_FILL <- "#F4F4F4"

y_top <- TOP - (seq_len(7L) - 1L) * (BOX_H + GAP)
y_bot <- y_top - BOX_H
y_mid <- (y_top + y_bot) / 2

fmt_n <- function(x) formatC(x, format = "d", big.mark = ",")

draw_block <- function(x, y_centre, lines, cex, lh) {
  n <- length(lines)
  ys <- y_centre + ((n - 1) / 2 - (seq_len(n) - 1)) * lh
  for (k in seq_len(n)) text(x, ys[k], lines[k], cex = cex[k], col = COL_INK)
  invisible(NULL)
}

draw_flow <- function() {
  op <- par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i", family = "Arial")
  on.exit(par(op))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  for (i in seq_len(7L)) {
    rect(XL[1], y_bot[i], XL[2], y_top[i], col = COL_FILL,
         border = COL_BORDER, lwd = 0.7)
    rl <- strwrap(pop[["rule"]][i], width = 54)
    draw_block(MID_L, y_mid[i], c(rl, paste0("n = ", fmt_n(pop[["n"]][i]))),
               cex = c(rep(0.58, length(rl)), 0.74), lh = 0.0132)
    if (i < 7L) {
      arrows(MID_L, y_bot[i] - 0.002, MID_L, y_top[i + 1L] + 0.002,
             length = 0.045, col = "grey25", lwd = 0.7)
    }
  }
  for (j in seq_len(6L)) {
    yy <- (y_bot[j] + y_top[j + 1L]) / 2
    rect(XR[1], yy - 0.036, XR[2], yy + 0.036, col = COL_EXC_FILL,
         border = COL_EXC_BORDER, lwd = 0.6)
    rl <- strwrap(exc[["rule"]][j], width = 50)
    draw_block(MID_R, yy,
               c(rl, paste0("n excluded = ", fmt_n(exc[["n"]][j]))),
               cex = c(rep(0.54, length(rl)), 0.66), lh = 0.0125)
    arrows(XL[2] + 0.004, yy, XR[1] - 0.004, yy, length = 0.045,
           col = "grey25", lwd = 0.7)
  }
  # branch to the two outcome cohorts
  arrows(MID_L, y_bot[7L] - 0.002, MID_L, BUS + 0.002, length = 0.03,
         col = "grey25", lwd = 0.7)
  segments(MID_L, BUS, MID_R, BUS, col = "grey25", lwd = 0.7)
  for (k in seq_len(2L)) {
    xx <- c(MID_L, MID_R)[k]
    arrows(xx, BUS - 0.001, xx, OUT_TOP + 0.002, length = 0.045,
           col = "grey25", lwd = 0.7)
    xl <- if (k == 1L) XL else XR
    rect(xl[1], OUT_BOT, xl[2], OUT_TOP, col = COL_FILL,
         border = COL_BORDER, lwd = 0.7)
    draw_block(mean(xl), (OUT_TOP + OUT_BOT) / 2,
               c(out[["rule"]][k], paste0("n = ", fmt_n(out[["n"]][k])),
                 paste0("Excluded at baseline: ",
                        fmt_n(out[["n_excluded"]][k]))),
               cex = c(0.62, 0.74, 0.54), lh = 0.021)
  }
  invisible(NULL)
}

grDevices::cairo_pdf(file.path(OUT_FIG, "11a_fig_sample_flow.pdf"),
                     width = W_MM / 25.4, height = H_MM / 25.4,
                     family = "Arial", onefile = TRUE)
draw_flow()
invisible(grDevices::dev.off())

ragg::agg_tiff(file.path(OUT_FIG, "11a_fig_sample_flow.tiff"),
               width = W_MM / 25.4, height = H_MM / 25.4, units = "in",
               res = 600, background = "white")
draw_flow()
invisible(grDevices::dev.off())

# ------------------------------------------------------------- plain check
bad_words <- c("gate0", "gate", "main cohort", "analytic file")
shown_text <- tolower(c(pop[["rule"]], exc[["rule"]], out[["rule"]]))
stopifnot(!any(vapply(bad_words, function(w) {
  any(grepl(w, shown_text, fixed = TRUE))
}, logical(1))))

# ------------------------------------------------------------------- report
lines <- c(
  "# Step 11a - sample flow diagram (final)",
  "",
  paste0("Generated by `R/11a_figure_sample_flow.R` on ",
         format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "."),
  "",
  "All counts are read from `output/tables/07_table_supp_sample_flow.csv`;",
  "the script asserts that every step balances and that the chain is continuous.",
  "",
  "| Box | Population | n | Excluded at this step |",
  "| --- | --- | --- | --- |",
  vapply(seq_len(7L), function(i) {
    paste0("| ", i, " | ", pop[["rule"]][i], " | ", fmt_n(pop[["n"]][i]), " | ",
           if (i == 1L) "-" else fmt_n(exc[["n"]][i - 1L]), " |")
  }, character(1)),
  "",
  "| Outcome cohort | n | Excluded at baseline |",
  "| --- | --- | --- |",
  vapply(seq_len(2L), function(k) {
    paste0("| ", out[["rule"]][k], " | ", fmt_n(out[["n"]][k]), " | ",
           fmt_n(out[["n_excluded"]][k]), " |")
  }, character(1)),
  "",
  paste0("Every visible count is at least 20; the smallest is ",
         fmt_n(min(all_counts)), "."),
  "",
  "Output: `output/figures/11a_fig_sample_flow.pdf` and `.tiff`, copied to",
  "`output/figures/final/Figure1_sample_flow.*` by the same run."
)
writeLines(lines, file.path(OUT_REP, "step11a_sample_flow.md"), useBytes = TRUE)

jsonlite::write_json(list(
  script = "R/11a_figure_sample_flow.R",
  run_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  input_file = FLOW_FILE,
  boxes = 7L,
  exclusions = 6L,
  outcome_boxes = 2L,
  merged_exclusion = list(steps = c("prior_ehr", "alive_adult"),
                          n = EXPECTED_MERGED_EXCL),
  first_n = EXPECTED_FIRST,
  last_n = EXPECTED_LAST,
  ckm_at_risk = EXPECTED_CKM_N,
  cardiorenal_at_risk = EXPECTED_CR_N,
  minimum_visible_count = min(all_counts),
  masking = "no count below 20 is displayed or written; the two small-count steps are merged",
  display_wording = "plain language; internal step names are not shown",
  width_mm = W_MM,
  height_mm = H_MM,
  outputs = c("output/figures/11a_fig_sample_flow.pdf",
              "output/figures/11a_fig_sample_flow.tiff",
              "output/figures/final/Figure1_sample_flow.pdf",
              "output/figures/final/Figure1_sample_flow.tiff",
              "output/tables/11a_sample_flow_boxes.csv",
              "output/reports/step11a_sample_flow.md"),
  r_version = R.version.string,
  git_commit = NULL
), file.path(OUT_LOG, "11a_provenance.json"), auto_unbox = TRUE, pretty = TRUE)

invisible(file.copy(file.path(OUT_FIG, "11a_fig_sample_flow.pdf"),
                    file.path(OUT_FIG, "final", "Figure1_sample_flow.pdf"), overwrite = TRUE))
invisible(file.copy(file.path(OUT_FIG, "11a_fig_sample_flow.tiff"),
                    file.path(OUT_FIG, "final", "Figure1_sample_flow.tiff"), overwrite = TRUE))
