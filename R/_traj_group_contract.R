# =============================================================================
# Shared trajectory-group contract: single source of truth
# Sourced by R/01c_write_groups.R, R/02_descriptive.R and R/03_mice.R.
# Definitions only, no side effects.
#
# On disk intermediate/traj_groups_w8q.csv keeps an integer code 1..4, because a
# CSV cannot store factor levels. Every reader must pass that column through
# as_traj_group(), which gives a factor with levels c("4", "1", "2", "3"):
# group 4 first, i.e. the reference group of every downstream model.
# TRAJ_LABELS holds the user-supplied Chinese group names; do not translate them.
# =============================================================================

TRAJ_G <- 4L            # number of trajectory groups, fixed by the user
TRAJ_REF <- 4L          # reference group, first factor level

TRAJ_GROUP_FILE <- file.path(Sys.getenv("FIG_ROOT", unset = "."),
                             "intermediate/traj_groups_w8q.csv")
TRAJ_GROUP_COL <- "traj_group"
TRAJ_ID_COL <- "person_id"

# level order: reference group first, then 1, 2, 3
TRAJ_LEVELS <- c(TRAJ_REF, setdiff(seq_len(TRAJ_G), TRAJ_REF))

TRAJ_LABELS <- c(
  `1` = "高活动–短睡眠",
  `2` = "较高活动–中等心率",
  `3` = "低活动–高心率",
  `4` = "中等活动–低心率"
)

# integer / character codes -> factor with levels c("4","1","2","3")
as_traj_group <- function(x) {
  xi <- suppressWarnings(as.integer(as.character(x)))
  if (any(is.na(xi))) stop("traj_group contains missing or non-integer values")
  if (!all(xi %in% seq_len(TRAJ_G))) {
    stop("traj_group values outside 1..", TRAJ_G, ": ",
         paste(sort(unique(xi[!(xi %in% seq_len(TRAJ_G))])), collapse = ","))
  }
  factor(as.character(xi), levels = as.character(TRAJ_LEVELS))
}

# factor or codes -> Chinese group label (character)
traj_label <- function(x) {
  unname(TRAJ_LABELS[as.character(x)])
}
