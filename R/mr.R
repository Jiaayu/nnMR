nn_mr <- function (dat, parameters = default_parameters(), method_list = subset(mr_method_list(),
                                                                                use_by_default)$obj)
{
  mr_tab <- plyr::ddply(dat, c("id.exposure", "id.outcome"),
                        function(x1) {
                          x <- subset(x1, mr_keep)
                          if (nrow(x) == 0) {
                            message("No SNPs available for MR analysis of '",
                                    x1$id.exposure[1], "' on '", x1$id.outcome[1],
                                    "'")
                            return(NULL)
                          }
                          else {
                            #message("Analysing '", x1$id.exposure[1], "' on '",
                            #        x1$id.outcome[1], "'")
                            a <- 1
                          }
                          res <- lapply(method_list, function(meth) {
                            get(meth)(x$beta.exposure, x$beta.outcome, x$se.exposure,
                                      x$se.outcome, parameters)
                          })
                          methl <- mr_method_list()
                          mr_tab <- data.frame(outcome = x$outcome[1], exposure = x$exposure[1],
                                               method = methl$name[match(method_list, methl$obj)],
                                               nsnp = sapply(res, function(x) x$nsnp), b = sapply(res,
                                                                                                  function(x) x$b), se = sapply(res, function(x) x$se),
                                               pval = sapply(res, function(x) x$pval))
                          mr_tab <- subset(mr_tab, !(is.na(b) & is.na(se) &
                                                       is.na(pval)))
                          return(mr_tab)
                        })
  return(mr_tab)
}


check_required_columns <- function (dat, type = "exposure")
{
  required_columns <- c("SNP", paste0(c("id.", "", "beta.",
                                        "se.", "effect_allele.", "other_allele."), type))
  index <- required_columns %in% names(dat)
  if (!all(index)) {
    stop("The following required columns are missing from ",
         type, ": ", paste(required_columns[!index], collapse = ", "))
  }
  return(NULL)
}


harmonise_cleanup_variables <- function (res.tab)
{
  res.tab$beta.exposure <- as.numeric(res.tab$beta.exposure)
  res.tab$beta.outcome <- as.numeric(res.tab$beta.outcome)
  res.tab$eaf.exposure <- as.numeric(res.tab$eaf.exposure)
  res.tab$eaf.outcome[res.tab$eaf.outcome == "NR"] <- NA
  res.tab$eaf.outcome[res.tab$eaf.outcome == "NR "] <- NA
  res.tab$eaf.outcome <- as.numeric(res.tab$eaf.outcome)
  res.tab$effect_allele.exposure <- toupper(res.tab$effect_allele.exposure)
  res.tab$other_allele.exposure <- toupper(res.tab$other_allele.exposure)
  res.tab$effect_allele.outcome <- toupper(res.tab$effect_allele.outcome)
  res.tab$other_allele.outcome <- toupper(res.tab$other_allele.outcome)
  res.tab$other_allele.outcome[res.tab$other_allele.outcome ==
                                 ""] <- NA
  return(res.tab)
}

harmonise <- function (dat, tolerance, action)
{
  dat$orig_SNP <- dat$SNP
  SNP_index <- sapply(1:length(dat$SNP), function(i) sum(dat$SNP[1:i] ==
                                                           dat$SNP[i]))
  dat$SNP <- paste0(dat$SNP, "_", SNP_index)
  SNP <- dat$SNP
  A1 <- dat$effect_allele.exposure
  A2 <- dat$other_allele.exposure
  B1 <- dat$effect_allele.outcome
  B2 <- dat$other_allele.outcome
  betaA <- dat$beta.exposure
  betaB <- dat$beta.outcome
  fA <- dat$eaf.exposure
  fB <- dat$eaf.outcome
  dat <- subset(dat, select = -c(effect_allele.exposure, other_allele.exposure,
                                 effect_allele.outcome, other_allele.outcome, beta.exposure,
                                 beta.outcome, eaf.exposure, eaf.outcome))
  i22 <- !is.na(A1) & !is.na(A2) & !is.na(B1) & !is.na(B2)
  i21 <- !is.na(A1) & !is.na(A2) & !is.na(B1) & is.na(B2)
  i12 <- !is.na(A1) & is.na(A2) & !is.na(B1) & !is.na(B2)
  i11 <- !is.na(A1) & is.na(A2) & !is.na(B1) & is.na(B2)
  d22 <- harmonise_22(SNP[i22], A1[i22], A2[i22], B1[i22],
                      B2[i22], betaA[i22], betaB[i22], fA[i22], fB[i22], tolerance,
                      action)
  d21 <- harmonise_21(SNP[i21], A1[i21], A2[i21], B1[i21],
                      betaA[i21], betaB[i21], fA[i21], fB[i21], tolerance,
                      action)
  d12 <- harmonise_12(SNP[i12], A1[i12], B1[i12], B2[i12],
                      betaA[i12], betaB[i12], fA[i12], fB[i12], tolerance,
                      action)
  d11 <- harmonise_11(SNP[i11], A1[i11], B1[i11], betaA[i11],
                      betaB[i11], fA[i11], fB[i11], tolerance, action)
  jlog <- plyr::rbind.fill(as.data.frame(attr(d22, "log"),
                                         stringsAsFactors = FALSE), as.data.frame(attr(d21, "log"),
                                                                                  stringsAsFactors = FALSE), as.data.frame(attr(d12, "log"),
                                                                                                                           stringsAsFactors = FALSE), as.data.frame(attr(d11, "log"),
                                                                                                                                                                    stringsAsFactors = FALSE))
  jlog <- cbind(data.frame(id.exposure = dat$id.exposure[1],
                           id.outcome = dat$id.outcome[1], stringsAsFactors = FALSE),
                jlog)
  d <- rbind(d21, d22, d12, d11)
  d <- merge(d, dat, by = "SNP", all.x = TRUE)
  d$SNP <- d$orig_SNP
  d <- subset(d, select = -orig_SNP)
  d <- d[order(d$id.outcome), ]
  d$mr_keep <- TRUE
  if (action == 3) {
    d$mr_keep[d$palindromic | d$remove | d$ambiguous] <- FALSE
    if (any(d$palindromic)) {
      #message("Removing the following SNPs for being palindromic:\n",
      #        paste(d$SNP[d$palindromic], collapse = ", "))
      xxxxxx <- 1
    }
    if (any(d$remove)) {
      #message("Removing the following SNPs for incompatible alleles:\n",
      #        paste(d$SNP[d$remove], collapse = ", "))
      xxxxxx <- 1
    }
    jlog[["incompatible_alleles"]] <- sum(d$remove)
    if (any(d$ambiguous & !d$palindromic)) {
      #message("Removing the following SNPs for having incompatible allele frequencies:\n",
      #        paste(d$SNP[d$ambiguous], collapse = ", "))
      xxxx <- 1
    }
    jlog[["ambiguous_alleles"]] <- sum(d$ambiguous)
  }
  if (action == 2) {
    d$mr_keep[d$remove | d$ambiguous] <- FALSE
    if (any(d$remove)) {
      #message("Removing the following SNPs for incompatible alleles:\n",
      #        paste(d$SNP[d$remove], collapse = ", "))
      xxxx <- 1
    }
    jlog[["incompatible_alleles"]] <- sum(d$remove)
    if (any(d$ambiguous)) {
      #message("Removing the following SNPs for being palindromic with intermediate allele frequencies:\n",
      #        paste(d$SNP[d$ambiguous], collapse = ", "))
      xxxx <- 1
    }
    jlog[["ambiguous_alleles"]] <- sum(d$ambiguous)
  }
  if (action == 1) {
    d$mr_keep[d$remove] <- FALSE
    if (any(d$remove)) {
      #message("Removing the following SNPs for incompatible alleles:\n",
      #        paste(d$SNP[d$remove], collapse = ", "))
      xxxx <- 1
    }
    jlog[["incompatible_alleles"]] <- sum(d$remove)
  }
  attr(d, "log") <- jlog
  return(d)
}


harmonise_22 <- function (SNP, A1, A2, B1, B2, betaA, betaB, fA, fB, tolerance,
                          action)
{
  if (length(SNP) == 0)
    return(data.frame())
  jlog <- list()
  jlog[["alleles"]] <- "2-2"
  indel_index <- nchar(A1) > 1 | nchar(A2) > 1 | A1 == "D" |
    A1 == "I"
  temp <- recode_indels_22(A1[indel_index], A2[indel_index],
                           B1[indel_index], B2[indel_index])
  A1[indel_index] <- temp$A1
  A2[indel_index] <- temp$A2
  B1[indel_index] <- temp$B1
  B2[indel_index] <- temp$B2
  status1 <- (A1 == B1) & (A2 == B2)
  to_swap <- (A1 == B2) & (A2 == B1)
  jlog[["switched_alleles"]] <- sum(to_swap)
  Btemp <- B1[to_swap]
  B1[to_swap] <- B2[to_swap]
  B2[to_swap] <- Btemp
  betaB[to_swap] <- betaB[to_swap] * -1
  fB[to_swap] <- 1 - fB[to_swap]
  status1 <- (A1 == B1) & (A2 == B2)
  palindromic <- check_palindromic(A1, A2)
  i <- !palindromic & !status1
  B1[i] <- flip_alleles(B1[i])
  B2[i] <- flip_alleles(B2[i])
  status1 <- (A1 == B1) & (A2 == B2)
  jlog[["flipped_alleles_basic"]] <- sum(i)
  i <- !palindromic & !status1
  to_swap <- (A1 == B2) & (A2 == B1)
  Btemp <- B1[to_swap]
  B1[to_swap] <- B2[to_swap]
  B2[to_swap] <- Btemp
  betaB[to_swap] <- betaB[to_swap] * -1
  fB[to_swap] <- 1 - fB[to_swap]
  status1 <- (A1 == B1) & (A2 == B2)
  remove <- !status1
  remove[indel_index][!temp$keep] <- TRUE
  minf <- 0.5 - tolerance
  maxf <- 0.5 + tolerance
  tempfA <- fA
  tempfB <- fB
  tempfA[is.na(tempfA)] <- 0.5
  tempfB[is.na(tempfB)] <- 0.5
  ambiguousA <- tempfA > minf & tempfA < maxf
  ambiguousB <- tempfB > minf & tempfB < maxf
  if (action == 2) {
    status2 <- ((tempfA < 0.5 & tempfB > 0.5) | (tempfA >
                                                   0.5 & tempfB < 0.5)) & palindromic
    to_swap <- status2 & !remove
    betaB[to_swap] <- betaB[to_swap] * -1
    fB[to_swap] <- 1 - fB[to_swap]
    if (!is.null(jlog)) {
      jlog[["flipped_alleles_palindrome"]] <- sum(to_swap)
    }
  }
  else {
    if (!is.null(jlog)) {
      jlog[["flipped_alleles_palindrome"]] <- 0
    }
  }
  d <- data.frame(SNP = SNP, effect_allele.exposure = A1, other_allele.exposure = A2,
                  effect_allele.outcome = B1, other_allele.outcome = B2,
                  beta.exposure = betaA, beta.outcome = betaB, eaf.exposure = fA,
                  eaf.outcome = fB, remove = remove, palindromic = palindromic,
                  ambiguous = (ambiguousA | ambiguousB) & palindromic)
  attr(d, "log") <- jlog
  return(d)
}

recode_indels_22 <- function (A1, A2, B1, B2)
{
  ncA1 <- nchar(A1)
  ncA2 <- nchar(A2)
  ncB1 <- nchar(B1)
  ncB2 <- nchar(B2)
  i1 <- ncA1 > ncA2 & B1 == "I" & B2 == "D"
  B1[i1] <- A1[i1]
  B2[i1] <- A2[i1]
  i1 <- ncA1 < ncA2 & B1 == "I" & B2 == "D"
  B1[i1] <- A2[i1]
  B2[i1] <- A1[i1]
  i1 <- ncA1 > ncA2 & B1 == "D" & B2 == "I"
  B1[i1] <- A2[i1]
  B2[i1] <- A1[i1]
  i1 <- ncA1 < ncA2 & B1 == "D" & B2 == "I"
  B1[i1] <- A1[i1]
  B2[i1] <- A2[i1]
  i1 <- ncB1 > ncB2 & A1 == "I" & A2 == "D"
  A1[i1] <- B1[i1]
  A2[i1] <- B2[i1]
  i1 <- ncB1 < ncB2 & A1 == "I" & A2 == "D"
  A2[i1] <- B1[i1]
  A1[i1] <- B2[i1]
  i1 <- ncB1 > ncB2 & A1 == "D" & A2 == "I"
  A2[i1] <- B1[i1]
  A1[i1] <- B2[i1]
  i1 <- ncB1 < ncB2 & A1 == "D" & A2 == "I"
  A1[i1] <- B1[i1]
  A2[i1] <- B2[i1]
  keep <- rep(TRUE, length(A1))
  keep[(ncA1 > 1 & ncA1 == ncA2 & (B1 == "D" | B1 == "I"))] <- FALSE
  keep[(ncB1 > 1 & ncB1 == ncB2 & (A1 == "D" | A1 == "I"))] <- FALSE
  keep[A1 == A2] <- FALSE
  keep[B1 == B2] <- FALSE
  return(data.frame(A1 = A1, A2 = A2, B1 = B1, B2 = B2, keep = keep,
                    stringsAsFactors = FALSE))
}


check_palindromic <- function (A1, A2)
{
  (A1 == "T" & A2 == "A") | (A1 == "A" & A2 == "T") | (A1 ==
                                                         "G" & A2 == "C") | (A1 == "C" & A2 == "G")
}

flip_alleles <- function (x) {
  x <- toupper(x)
  x <- gsub("C", "g", x)
  x <- gsub("G", "c", x)
  x <- gsub("A", "t", x)
  x <- gsub("T", "a", x)
  return(toupper(x))
}


harmonise_21 <- function (SNP, A1, A2, B1, betaA, betaB, fA, fB, tolerance, action)
{
  if (length(SNP) == 0)
    return(data.frame())
  jlog <- list()
  jlog[["alleles"]] <- "2-1"
  n <- length(A1)
  B2 <- rep(NA, n)
  ambiguous <- rep(FALSE, n)
  palindromic <- check_palindromic(A1, A2)
  remove <- palindromic
  indel_index <- nchar(A1) > 1 | nchar(A2) > 1 | A1 == "D" |
    A1 == "I"
  temp <- recode_indels_21(A1[indel_index], A2[indel_index],
                           B1[indel_index])
  A1[indel_index] <- temp$A1
  A2[indel_index] <- temp$A2
  B1[indel_index] <- temp$B1
  B2[indel_index] <- temp$B2
  remove[indel_index][!temp$keep] <- TRUE
  status1 <- A1 == B1
  minf <- 0.5 - tolerance
  maxf <- 0.5 + tolerance
  tempfA <- fA
  tempfB <- fB
  tempfA[is.na(tempfA)] <- 0.5
  tempfB[is.na(tempfB)] <- 0.5
  freq_similar1 <- (tempfA < minf & tempfB < minf) | (tempfA >
                                                        maxf & tempfB > maxf)
  ambiguous[status1 & !freq_similar1] <- TRUE
  B2[status1] <- A2[status1]
  to_swap <- A2 == B1
  jlog[["switched_alleles"]] <- sum(to_swap)
  freq_similar2 <- (tempfA < minf & tempfB > maxf) | (tempfA >
                                                        maxf & tempfB < minf)
  ambiguous[to_swap & !freq_similar2] <- TRUE
  B2[to_swap] <- B1[to_swap]
  B1[to_swap] <- A1[to_swap]
  betaB[to_swap] <- betaB[to_swap] * -1
  fB[to_swap] <- 1 - fB[to_swap]
  to_flip <- A1 != B1 & A2 != B1
  jlog[["flipped_alleles_no_oa"]] <- sum(to_flip)
  ambiguous[to_flip] <- TRUE
  B1[to_flip] <- flip_alleles(B1[to_flip])
  status1 <- A1 == B1
  B2[status1] <- A2[status1]
  to_swap <- A2 == B1
  B2[to_swap] <- B1[to_swap]
  B1[to_swap] <- A1[to_swap]
  betaB[to_swap] <- betaB[to_swap] * -1
  fB[to_swap] <- 1 - fB[to_swap]
  d <- data.frame(SNP = SNP, effect_allele.exposure = A1, other_allele.exposure = A2,
                  effect_allele.outcome = B1, other_allele.outcome = B2,
                  beta.exposure = betaA, beta.outcome = betaB, eaf.exposure = fA,
                  eaf.outcome = fB, remove = remove, palindromic = palindromic,
                  ambiguous = ambiguous | palindromic)
  attr(d, "log") <- jlog
  return(d)
}

harmonise_12 <- function (SNP, A1, B1, B2, betaA, betaB, fA, fB, tolerance, action)
{
  if (length(SNP) == 0)
    return(data.frame())
  jlog <- list()
  jlog[["alleles"]] <- "1-2"
  n <- length(A1)
  A2 <- rep(NA, n)
  ambiguous <- rep(FALSE, n)
  palindromic <- check_palindromic(B1, B2)
  remove <- palindromic
  indel_index <- nchar(B1) > 1 | nchar(B2) > 1 | B1 == "D" |
    B1 == "I"
  temp <- recode_indels_21(A1[indel_index], B1[indel_index],
                           B2[indel_index])
  A1[indel_index] <- temp$A1
  A2[indel_index] <- temp$A2
  B1[indel_index] <- temp$B1
  B2[indel_index] <- temp$B2
  remove[indel_index][!temp$keep] <- TRUE
  status1 <- A1 == B1
  minf <- 0.5 - tolerance
  maxf <- 0.5 + tolerance
  tempfA <- fA
  tempfB <- fB
  tempfA[is.na(tempfA)] <- 0.5
  tempfB[is.na(tempfB)] <- 0.5
  freq_similar1 <- (tempfA < minf & tempfB < minf) | (tempfA >
                                                        maxf & tempfB > maxf)
  ambiguous[status1 & !freq_similar1] <- TRUE
  A2[status1] <- B2[status1]
  to_swap <- A1 == B2
  jlog[["switched_alleles"]] <- sum(to_swap)
  freq_similar2 <- (tempfA < minf & tempfB > maxf) | (tempfA >
                                                        maxf & tempfB < minf)
  ambiguous[to_swap & !freq_similar2] <- TRUE
  A2[to_swap] <- A1[to_swap]
  A1[to_swap] <- B1[to_swap]
  betaA[to_swap] <- betaA[to_swap] * -1
  fA[to_swap] <- 1 - fA[to_swap]
  to_flip <- A1 != B1 & A1 != B2
  jlog[["flipped_alleles_no_oa"]] <- sum(to_flip)
  ambiguous[to_flip] <- TRUE
  A1[to_flip] <- flip_alleles(A1[to_flip])
  status1 <- A1 == B1
  A2[status1] <- B2[status1]
  to_swap <- B2 == A1
  B2[to_swap] <- B1[to_swap]
  B1[to_swap] <- A1[to_swap]
  betaB[to_swap] <- betaB[to_swap] * -1
  fB[to_swap] <- 1 - fB[to_swap]
  d <- data.frame(SNP = SNP, effect_allele.exposure = A1, other_allele.exposure = A2,
                  effect_allele.outcome = B1, other_allele.outcome = B2,
                  beta.exposure = betaA, beta.outcome = betaB, eaf.exposure = fA,
                  eaf.outcome = fB, remove = remove, palindromic = palindromic,
                  ambiguous = ambiguous | palindromic)
  attr(d, "log") <- jlog
  return(d)
}

harmonise_11 <- function (SNP, A1, B1, betaA, betaB, fA, fB, tolerance, action)
{
  if (length(SNP) == 0)
    return(data.frame())
  jlog <- list()
  jlog[["alleles"]] <- "1-1"
  n <- length(A1)
  A2 <- rep(NA, n)
  B2 <- rep(NA, n)
  ambiguous <- rep(FALSE, n)
  palindromic <- FALSE
  status1 <- A1 == B1
  remove <- !status1
  minf <- 0.5 - tolerance
  maxf <- 0.5 + tolerance
  tempfA <- fA
  tempfB <- fB
  tempfA[is.na(tempfA)] <- 0.5
  tempfB[is.na(tempfB)] <- 0.5
  freq_similar1 <- (tempfA < minf & tempfB < minf) | (tempfA >
                                                        maxf & tempfB > maxf)
  ambiguous[status1 & !freq_similar1] <- TRUE
  d <- data.frame(SNP = SNP, effect_allele.exposure = A1, other_allele.exposure = A2,
                  effect_allele.outcome = B1, other_allele.outcome = B2,
                  beta.exposure = betaA, beta.outcome = betaB, eaf.exposure = fA,
                  eaf.outcome = fB, remove = remove, palindromic = palindromic,
                  ambiguous = ambiguous | palindromic)
  attr(d, "log") <- jlog
  return(d)
}
