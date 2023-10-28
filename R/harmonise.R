#' @export
nn_harmonise_data <- function (exposure_dat, outcome_dat, action = 2)
{
  stopifnot(all(action %in% 1:3))
  check_required_columns(exposure_dat, "exposure")
  check_required_columns(outcome_dat, "outcome")
  res.tab <- merge(outcome_dat, exposure_dat, by = "SNP")
  ncombinations <- length(unique(res.tab$id.outcome))
  if (length(action) == 1) {
    action <- rep(action, ncombinations)
  }
  else if (length(action) != ncombinations) {
    stop("Action argument must be of length 1 (where the same action will be used for all outcomes), or number of unique id.outcome values (where each outcome will use a different action value)")
  }
  res.tab <- harmonise_cleanup_variables(res.tab)
  if (nrow(res.tab) == 0) {
    return(res.tab)
  }
  d <- data.frame(id.outcome = unique(res.tab$id.outcome),
                  action = action)
  res.tab <- merge(res.tab, d, by = "id.outcome")
  combs <- subset(res.tab, !duplicated(paste(id.exposure, id.outcome)),
                  select = c(id.exposure, id.outcome))
  fix.tab <- list()
  mr_cols <- c("beta.exposure", "beta.outcome", "se.exposure",
               "se.outcome")
  for (i in 1:nrow(combs)) {
    x <- subset(res.tab, id.exposure == combs$id.exposure[i] &
                  id.outcome == combs$id.outcome[i])
    #message("Harmonising ", x$exposure[1], " (", x$id.exposure[1],
    #        ") and ", x$outcome[1], " (", x$id.outcome[1], ")")
    x <- harmonise(x, 0.08, x$action[1])
    attr(x, "log")[["candidate_variants"]] <- sum(exposure_dat$id.exposure ==
                                                    x$id.exposure[1])
    attr(x, "log")[["variants_absent_from_reference"]] <- sum(exposure_dat$id.exposure ==
                                                                x$id.exposure[1]) - nrow(x)
    x$mr_keep[apply(x[, mr_cols], 1, function(y) any(is.na(y)))] <- FALSE
    attr(x, "log")[["total_variants"]] <- nrow(x)
    attr(x, "log")[["total_variants_for_mr"]] <- sum(x$mr_keep)
    attr(x, "log")[["proxy_variants"]] <- ifelse(is.null(x$proxy.outcome),
                                                 0, sum(x$proxy.outcome, na.rm = TRUE))
    fix.tab[[i]] <- x
  }
  jlog <- plyr::rbind.fill(lapply(fix.tab, function(x) attr(x,
                                                            "log")))
  fix.tab <- plyr::rbind.fill(fix.tab)
  attr(fix.tab, "log") <- jlog
  if (!"samplesize.outcome" %in% names(fix.tab)) {
    fix.tab$samplesize.outcome <- NA
  }
  return(fix.tab)
}
