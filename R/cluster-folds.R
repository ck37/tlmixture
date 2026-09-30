#' Cluster-grouped CV-TMLE folds
#'
#' Builds CV-TMLE folds that keep every cluster within a single fold, using
#' rsample::group_vfold_cv(). The grouping column is only used to assign the
#' folds: the returned splits are rebuilt on \code{data} itself, so the
#' grouping never appears in the analysis or assessment data (where it would
#' otherwise be treated as a confounder).
#'
#' group_vfold_cv() requires strata to be constant within each group, so the
#' outcome itself cannot be used. Instead folds are stratified on each
#' cluster's mean outcome (its event rate for a binary outcome), which keeps
#' the outcome balanced across folds as far as whole clusters allow. If
#' stratification is not possible (e.g. too few clusters), folds are grouped
#' without it.
#'
#' @param data The analysis data frame, without cluster id columns.
#' @param groups Grouping vector, one value per row, from
#'   cluster_fold_groups().
#' @param v Number of folds.
#' @param strata_values Optional numeric vector, one value per row (the
#'   outcome), summarized to a per-cluster mean for stratification.
#'
#' @return An rsample rset whose splits index rows of \code{data}.
#' @noRd
cluster_vfold_cv = function(data, groups, v, strata_values = NULL) {
  if (length(groups) != nrow(data)) {
    stop("The cluster grouping has ", length(groups), " values but data has ",
         nrow(data), " rows.")
  }
  if (length(unique(groups)) < v) {
    stop("There are ", length(unique(groups)), " clusters to form folds on, ",
         "fewer than folds_cvtmle = ", v, ".")
  }
  group_df = data.frame(.cluster_group = groups)

  grouped = NULL
  if (!is.null(strata_values)) {
    # Per-cluster mean outcome, constant within each cluster as rsample needs.
    group_df$.cluster_strata = stats::ave(as.numeric(strata_values), groups,
                                          FUN = function(x) mean(x, na.rm = TRUE))
    grouped = tryCatch(
      suppressWarnings(
        rsample::group_vfold_cv(group_df, group = ".cluster_group", v = v,
                                strata = ".cluster_strata")),
      error = function(e) NULL)
  }
  if (is.null(grouped)) {
    grouped = rsample::group_vfold_cv(group_df, group = ".cluster_group", v = v)
  }

  splits = lapply(grouped$splits, function(split) {
    rsample::make_splits(list(analysis = split$in_id,
                              assessment = rsample::complement(split)),
                         data = data)
  })
  rsample::manual_rset(splits, ids = grouped$id)
}

#' Cross-fitting fold assignment, keeping clusters together when ids are given
#'
#' Without \code{id} this is exactly the previous assignment,
#' \code{sample(rep(seq(cv_folds), length.out = n))}, so results are unchanged.
#' With \code{id}, whole clusters are assigned to folds, and the number of
#' folds is capped at the number of clusters.
#'
#' @return A list with \code{folds} (one fold per row) and \code{cv_folds}.
#' @noRd
crossfit_folds = function(n, cv_folds, id = NULL) {
  if (is.null(id)) {
    return(list(folds = sample(rep(seq(cv_folds), length.out = n)),
                cv_folds = cv_folds))
  }
  clusters = unique(id)
  cv_folds = min(cv_folds, length(clusters))
  cluster_folds = sample(rep(seq(cv_folds), length.out = length(clusters)))
  list(folds = cluster_folds[match(id, clusters)], cv_folds = cv_folds)
}

#' SuperLearner id and number of CV folds, given optional cluster ids
#'
#' SuperLearner keeps each id within one internal fold, but needs at least as
#' many distinct ids as folds, and cannot also stratify (callers turn
#' stratification off when \code{id} is returned). Without clusters this
#' returns \code{id = NULL} and \code{V} unchanged.
#' @noRd
sl_cluster_args = function(id, V) {
  if (is.null(id)) return(list(id = NULL, V = V))
  V_use = cluster_sl_folds(id, V)
  if (V_use < 2L) return(list(id = NULL, V = 2L))
  list(id = id, V = V_use)
}

#' TRUE if calling \code{fn} with an argument named \code{arg} is safe: it has
#' that argument or \code{...}.
#' @noRd
accepts_arg = function(fn, arg) {
  args = names(formals(fn))
  arg %in% args || "..." %in% args
}
