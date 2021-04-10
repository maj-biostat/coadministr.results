## -----------------------------------------------------------------------------
##
## Script name: process_dat.R
##
## Purpose of script: extract summaries from the simulation outputs
##
## Author: MAJ
##
## Date Created: 2021-03-01
##
## -----------------------------------------------------------------------------
##
## Process the simulation results as stored in the *.rds files
## found under dat. Invoke each of + process_inferiority, +
## process_noninferiority + process_combined_rules manually to process the rules
## of interest. The `operating_characteristics.Rmd` file is dependent on the
## results from this processing.
##
## -----------------------------------------------------------------------------

library(data.table)
library(ggplot2)
library(parallel)

# globals
g_dec_thresh_inf <- 0.995
g_dec_thresh_ninf <- 0.995

plot1 <- function(dfig, fnfig = NULL,
                  width = 6, height = 4,
                  cohort = FALSE){
  p <- ggplot(data = dfig,
              aes(x = es, y = mu_pwr,
                  group = brand1, col = brand1)) +
    geom_point() +
    geom_line() +
    theme_bw() +
    scale_x_continuous("Effect size (prob scale)") +
    scale_y_continuous("Power")  +
    scale_color_discrete("Brand") +
    theme(legend.position="bottom")

  if(cohort){
    p <- p + facet_wrap(~paste0("Cohort ", cohort1))
  }

  ggsave(fnfig, plot = p,
         width = width, height = height,
         units = "in")
  return(p)
}

plot2 <- function(dfig, fnfig = NULL,
                  width = 6, height = 4,
                  cohort = FALSE){
  p <- ggplot(data = dfig,
              aes(x = es, y = mu_n,
                  group = brand1, col = brand1)) +
    geom_point() +
    geom_line() +
    theme_bw() +
    scale_x_continuous("Effect size (prob scale)") +
    scale_y_continuous("Total enrolled")  +
    scale_color_discrete("Brand") +
    theme(legend.position="bottom")

  if(cohort){
    p <- p + facet_wrap(~paste0("Cohort ", cohort1))
  }

  ggsave(fnfig, plot = p,
         width = width, height = height,
         units = "in")
  return(p)
}

extract1 <- function(res = NULL){

  # fns <- list.files(path = "dat", pattern = "^0")
  # r <- readRDS(file.path("dat", fns[10]))
  # res <- r$res
  # z <- res[[1]]
  # w <- z$lr[[1]]

  stopifnot(!is.null(res))
  lpar <- res[[1]]$lpar
  # For each simulation z (that contains the dataset, parameters and
  # results summary) data.table with row for each arm coad, brand, cohort.
  isim <- 1
  dd <- do.call(rbind, lapply(res, function(z) {
    stopifnot(!is.na(z$lr))
    for(idx in 1:ncol(z$lpar$beta)){
      stopifnot(identical(lpar$beta[, idx, with = FALSE],
                          z$lpar$beta[, idx, with = FALSE]))
    }
    # Maximum follow up time for each participant
    tbl <- z$d[, .(curr_time = max(futime)), keyby = id]
    dd <- do.call(rbind, lapply(z$lr, function(w){
      # message()
      bb <- z$lpar$beta[
        timept == 1,.(es, pae, brand, cohort, lab)][
          order(lab, brand, cohort)]

      d <- data.table(sim_id = isim,
                      look = w$look,
                      enrolled = lpar$looks[w$look],
                      coad = w$dtbl$arm,
                      brand = w$dtbl$brand,
                      cohort = w$dtbl$cohort,
                      y = w$dtbl$y,
                      n = w$dtbl$n,
                      # empirical prob
                      p_obs = w$dtbl$observed,
                      p_tru = bb$pae,
                      # posterior mean prob
                      p_hat = w$p,
                      es = bb$es,
                      curr_time = as.numeric(tbl$curr_time[lpar$looks[w$look]])
      )
      d
    })) # w
    isim <<- isim + 1
    dd
  })) # z
  list(dd = dd, lpar = lpar)
}

extract2 <- function(res = NULL){
  # i = 7; r <- readRDS(fns[i]); res = r$res; z <- res[[1]]; w <- z$lr[[1]]
  # z <- res[[2]]
  stopifnot(!is.null(res))
  lpar <- res[[1]]$lpar
  isim <- 1
  # For each simulation z (that contains the dataset, parameters and
  # results summary) data.table with row for each arm coad, brand.
  isim <- 1
  dd <-  do.call(rbind, lapply(res, function(z) {
    stopifnot(!is.na(z$lr))
    for(idx in 1:ncol(z$lpar$beta)){
      stopifnot(identical(lpar$beta[, idx, with = FALSE],
                          z$lpar$beta[, idx, with = FALSE]))
    }
    # Maximum follow up time for each participant
    tbl <- z$d[, .(curr_time = max(futime)), keyby = id]

    # Results by coad and brand (avg out cohort)
    dd <- do.call(rbind, lapply(z$lr, function(w){
      # message()
      dtbl <- w$dtbl[, .(y = sum(y), n = sum(n)), keyby = .(arm, brand)]
      dtbl[, observed := y/n]
      d <- data.table(sim_id = isim,
                      look = w$look,
                      enrolled = lpar$looks[w$look],
                      coad = dtbl$arm,
                      brand = dtbl$brand,
                      y = dtbl$y,
                      n = dtbl$n,
                      # empirical prob
                      p_obs = dtbl$observed,
                      # posterior mean prob
                      p_hat = c(w$p_az[1], w$p_pf[2], w$p_az[2], w$p_pf[2]),
                      curr_time = as.numeric(tbl$curr_time[lpar$looks[w$look]])
      )
      d
    })) # w
    isim <<- isim + 1
    dd
  })) # z
  list(dd = dd, lpar = lpar)
}

extract3 <- function(res = NULL){

  # fns <- list.files(path = "dat", pattern = "^0")
  # r <- readRDS(file.path("dat", fns[3]))
  # res <- r$res
  # z <- res[[1]]
  # w <- z$lr[[1]]

  stopifnot(!is.null(res))
  lpar <- res[[1]]$lpar
  # Comparisons (for brand and cohort)
  isim <- 1
  dd <-  do.call(rbind, lapply(res, function(z) {
    stopifnot(!is.na(z$lr))
    for(idx in 1:ncol(z$lpar$beta)){
      stopifnot(identical(lpar$beta[, idx, with = FALSE],
                          z$lpar$beta[, idx, with = FALSE]))
    }
    # Maximum follow up time for each participant
    tbl <- z$d[, .(curr_time = max(futime)), keyby = id]
    # Results by coad and brand (avg out cohort)
    dd <- do.call(rbind, lapply(z$lr, function(w){

      # the mean effect size for each brand
      es <- z$lpar$beta[
        timept == 1 & arm == 1, .(es = mean(es)), keyby = .(brand, cohort)]

      dtbl <- cbind(w$dtbl[1:4, ], w$dtbl[5:8, ])
      names(dtbl) <- paste0(names(w$dtbl), rep(1:2, each = 6))
      d <- data.table(sim_id = isim,
                      look = w$look,
                      enrolled = lpar$looks[w$look])
      d <- cbind(d, dtbl)
      d[, es := es$es]
      d[, del_ni := z$lpar$delta_ni]
      d[, del_obs := observed1 - observed2]
      d[, del_hat := w$delta]
      d[, pr_inf := w$pr_inf]
      d[, pr_ninf := w$pr_ninf]
      d[, curr_time := as.numeric(tbl$curr_time[lpar$looks[w$look]])]
      d
    })) # w
    isim <<- isim + 1
    dd
  })) # z
  list(dd = dd, lpar = lpar)
}

extract4 <- function(res = NULL){

  # fns <- list.files(path = "dat", pattern = "^0")
  # r <- readRDS(file.path("dat", fns[2]))
  # res <- r$res
  # z <- res[[1]]
  # w <- z$lr[[1]]

  stopifnot(!is.null(res))
  lpar <- res[[1]]$lpar
  # Comparisons (brand)
  isim <- 1
  dd <-  do.call(rbind, lapply(res, function(z) {
    stopifnot(!is.na(z$lr))
    for(idx in 1:ncol(z$lpar$beta)){
      stopifnot(identical(lpar$beta[, idx, with = FALSE],
                          z$lpar$beta[, idx, with = FALSE]))
    }
    # Maximum follow up time for each participant
    tbl <- z$d[, .(curr_time = max(futime)), keyby = id]
    # Results by coad and brand (avg out cohort)
    dd <- do.call(rbind, lapply(z$lr, function(w){

      # the mean effect size for each brand
      es <- z$lpar$beta[
        timept == 1 & arm == 1, .(es = mean(es)), keyby = brand]

      dtbl <- w$dtbl[, .(y = sum(y), n = sum(n)), keyby = .(arm, brand)]
      dtbl[, observed := y/n]
      dtbl <- cbind(dtbl[1:2, ], dtbl[3:4, ])
      names(dtbl) <- paste0(names(dtbl), rep(1:2, each = 5))
      d <- data.table(sim_id = isim,
                      look = w$look,
                      enrolled = lpar$looks[w$look])
      d <- cbind(d, dtbl)
      d[, es := es$es]
      d[, del_ni := z$lpar$delta_ni]
      d[, del_obs := observed1 - observed2]
      d[, del_hat := c(w$delta_az, w$delta_pf)]
      d[, pr_inf := c(w$pr_inf_az, w$pr_inf_pf)]
      d[, pr_ninf := c(w$pr_ninf_az, w$pr_ninf_pf)]
      d[, curr_time := as.numeric(tbl$curr_time[lpar$looks[w$look]])]
      d
    })) # w
    isim <<- isim + 1
    dd
  })) # z
  list(dd = dd, lpar = lpar)
}

extract_all <- function(dir = "dat"){
  if(.Platform$OS.type == "unix") {
    ncores <- parallel::detectCores()-2
  } else {
    ncores <- 1
  }
  fns <- list.files(path = dir, pattern = "^0")
  i = 1;
  extractor <- function(ff = NULL){
    dd <- do.call(
      rbind,
      mclapply(X=seq_along(fns),
               mc.cores = ncores,
               FUN=function(i) {
                 r <- readRDS(file.path(dir, fns[i]))
                 l <- ff(r$res)
                 d <- l$dd
                 d$fn <- fns[i]
                 d$pars_id <- i
                 d
               }))
    dd
  }

  # data by coad, brand and cohort, estimates etc.
  message("Extract 1...")
  dd1 <- extractor(extract1)
  message("Extract 2...")
  dd2 <- extractor(extract2)
  message("Extract 3...")
  dd3 <- extractor(extract3)
  message("Extract 4...")
  dd4 <- extractor(extract4)

  list(dd1 = dd1, dd2 = dd2, dd3 = dd3, dd4 = dd4)
}


# Need to check for anomolous sample size inferiority at es = 0
process <- function(fn = "l.rds"){
  l <- extract_all()
  saveRDS(l, fn)
  l
}

ref_grid <- function(l, id = 1){

  if(id == 1){
    # for ensuring complete set of parsid
    es <- sort(unique(l$dd4$es))
    dref <- data.table(
      expand.grid(pars_id = sort(unique(l$dd1$pars_id)),
                  brand1 = c("AZ", "PF"))
    )
    dref[, es := rep(es, length = .N)]
  } else if(id == 2){
    # for ensuring complete set of parsid
    es <- sort(unique(l$dd4$es))
    dref <- data.table(
      expand.grid(pars_id = sort(unique(l$dd1$pars_id)),
                  brand1 = c("AZ", "PF"),
                  cohort1 = 1:2)
    )
    dref[, es := rep(es, length = .N)]
  }

  return(dref)
}

decisions3 <- function(l){
  # The power for each dose is described.
  dd <- copy(l$dd3)
  suppressWarnings(dd[, inf := NULL])
  suppressWarnings(dd[, ninf := NULL])
  suppressWarnings(dd[, fn := NULL])

  # were there any inferiority decisions at either of the covid doses?
  dd[, inf := as.numeric(pr_inf > g_dec_thresh_inf)]
  # cumulative decisions of inf over the duration of trial
  dd[, inf := cumsum(inf), keyby = .(sim_id, pars_id, brand1, cohort1)]
  dd[inf >= 1, inf := 1]

  # any noninferiority decisions at either of the covid doses?
  dd[, ninf := as.numeric(pr_ninf > g_dec_thresh_ninf)]
  # cumulative decisions of ninf over the duration of trial
  dd[, ninf := cumsum(ninf), keyby = .(sim_id, pars_id, brand1, cohort1)]
  dd[ninf >= 1, ninf := 1]

  dd[, anystop := as.numeric(inf | ninf), by = .(sim_id, pars_id, brand1, cohort1) ]
  dd[, anystop := cumsum(anystop), by = .(sim_id, pars_id, brand1, cohort1) ]
  dd[anystop >= 1, anystop := 1]

  dd[anystop == 1, decision := NA_character_ ]
  dd[anystop == 1 & inf == 1, decision := "Inferior" ]
  dd[anystop == 1 & ninf == 1 & is.na(decision), decision := "Noninferior" ]

  dd
}

decisions4 <- function(l){
  # The power for each dose is described.
  dd <- copy(l$dd4)
  suppressWarnings(dd[, inf := NULL])
  suppressWarnings(dd[, ninf := NULL])
  suppressWarnings(dd[, fn := NULL])

  # were there any inferiority decisions at either of the covid doses?
  dd[, inf := as.numeric(pr_inf > g_dec_thresh_inf)]
  # cumulative decisions of inf over the duration of trial
  dd[, inf := cumsum(inf), keyby = .(sim_id, pars_id, brand1)]
  dd[inf >= 1, inf := 1]

  # any noninferiority decisions at either of the covid doses?
  dd[, ninf := as.numeric(pr_ninf > g_dec_thresh_ninf)]
  # cumulative decisions of ninf over the duration of trial
  dd[, ninf := cumsum(ninf), keyby = .(sim_id, pars_id, brand1)]
  dd[ninf >= 1, ninf := 1]

  dd[, anystop := as.numeric(inf | ninf), by = .(sim_id, pars_id, brand1) ]
  dd[, anystop := cumsum(anystop), by = .(sim_id, pars_id, brand1) ]
  dd[anystop >= 1, anystop := 1]

  dd[anystop == 1, decision := NA_character_ ]
  dd[anystop == 1 & inf == 1, decision := "Inferior" ]
  dd[anystop == 1 & ninf == 1 & is.na(decision), decision := "Noninferior" ]

  dd
}

# Inferiority by brand:
# Power by es
# Total enrolled participants by es
inferiority1 <- function(l, ff = NULL){

  # stopifnot(!is.null(ff))

  dref <- ref_grid(l, 1)
  sort(unique(l$dd1$enrolled))

  # What/.when decisions made?
  dd <- decisions4(l)

  # ddtt1 <- dd[inf == 1 & enrolled == 500, ]
  # ddtt2 <- l$dd1[sim_id == 1 & pars_id == 10 & enrolled == 500 , ]
  #
  # ddtt1[sim_id == 1 & enrolled == 500, ]
  # ddtt2


  # power = proportion of simulations (by brand) that we would declare
  # inferiority at the nominated decision threshold
  dtmp1 <- dd[inf == 1, .(mu_pwr = length(unique(sim_id))/max(dd$sim_id) ),
              keyby = .(pars_id, brand1, es)]
  dtmp1 <- merge(dref, dtmp1, by = c("pars_id", "brand1", "es"), all.x = TRUE)
  dtmp1[is.na(mu_pwr), mu_pwr := 0]
  # sample size
  dtmp2 <- dd[inf == 1, .(mu_n = mean(enrolled),
                          sd_n = sd(enrolled)),
              keyby = .(pars_id, brand1, es)]
  dtmp2 <- merge(dref, dtmp2, by = c("pars_id", "brand1", "es"), all.x = TRUE)

  # Power
  p1 <- plot1(dtmp1, fnfig = "fig/pwr_inf1.png")#; print(p1)
  # Total participants enrolled when we hit an inferiority decision
  p2 <- plot2(dtmp2, fnfig = "fig/n_inf1.png")#; print(p2)

  # Summary tables
  dtmp <- merge(dtmp1, dtmp2, by = c("pars_id", "brand1","es"))
  k1 <- knitr::kable(dtmp[
    brand1 == "AZ", .(es, mu_pwr,
                      mu_n = round(mu_n),
                      sd_n = round(sd_n))][
                        order(es)], digits = 3)

  k2 <- knitr::kable(dtmp[
    brand1 == "PF", .(es, mu_pwr,
                      mu_n = round(mu_n),
                      sd_n = round(sd_n))][
                        order(es)], digits = 3)

  list(p1 = p1, p2 = p2, dtmp1 = dtmp1, dtmp2 = dtmp2, k1 = k1, k2 = k2)
}

# Noninferiority by brand:
# Power by es
# Total enrolled participants by es
non_inferiority1 <- function(l){
  dref <- ref_grid(l, 1)
  sort(unique(l$dd1$enrolled))

  # What/.when decisions made?
  dd <- decisions4(l)

  # power = proportion of simulations (by brand) that we would declare
  # inferiority at the nominated decision threshold
  dtmp1 <- dd[ninf == 1, .(mu_pwr = length(unique(sim_id))/max(dd$sim_id) ),
              keyby = .(pars_id, brand1, es)]
  dtmp1 <- merge(dref, dtmp1, by = c("pars_id", "brand1", "es"), all.x = TRUE)
  dtmp1[is.na(mu_pwr), mu_pwr := 0]
  # Total participants enrolled
  dtmp2 <- dd[ninf == 1, .(mu_n = mean(enrolled),
                          sd_n = sd(enrolled)),
              keyby = .(pars_id, brand1, es)]
  dtmp2 <- merge(dref, dtmp2, by = c("pars_id", "brand1", "es"), all.x = TRUE)

  # Power
  p1 <- plot1(dtmp1, fnfig = "fig/pwr_ninf1.png");  #print(p1)
  # Total participants enrolled when we hit an inferiority decision
  p2 <- plot2(dtmp2, fnfig = "fig/n_ninf1.png");  #print(p2)
  # Summary tables
  dtmp <- merge(dtmp1, dtmp2, by = c("pars_id", "brand1","es"))
  k1 <- knitr::kable(dtmp[
    brand1 == "AZ", .(es, mu_pwr,
                      mu_n = round(mu_n),
                      sd_n = round(sd_n))][
                        order(es)], digits = 3)
  k2 <- knitr::kable(dtmp[
    brand1 == "PF", .(es, mu_pwr,
                      mu_n = round(mu_n),
                      sd_n = round(sd_n))][
                        order(es)], digits = 3)
  list(p1 = p1, p2 = p2, dtmp1 = dtmp1, dtmp2 = dtmp2, k1 = k1, k2 = k2)
}

# Inferiority or Noninferiority by brand:
# Power by es
# Total enrolled participants by es
combined1 <- function(l){
  dref <- ref_grid(l, 1)
  sort(unique(l$dd1$enrolled))

  dd <- decisions4(l)

  # compute power as the proportion of times that we declared inferiority
  # over the course of a trial out of the simulations for that specific
  # configuration.
  dtmp1 <- dd[anystop == 1,
              .(mu_pwr = length(unique(sim_id))/max(dd$sim_id) ),
              keyby = .(pars_id, brand1, es)]
  dtmp1 <- merge(dref, dtmp1, by = c("pars_id", "brand1", "es"), all.x = TRUE)
  dtmp1[is.na(mu_pwr), mu_pwr := 0]
  # Total participants enrolled
  dtmp2 <- dd[anystop == 1, .(mu_n = mean(enrolled),
                              sd_n = sd(enrolled)),
              keyby = .(pars_id, brand1, es)]
  dtmp2 <- merge(dref, dtmp2, by = c("pars_id", "brand1", "es"), all.x = TRUE)

  # Power
  p1 <- plot1(dtmp1, fnfig = "fig/pwr_cmb1.png"); #print(p1)
  # Total participants enrolled when we hit an inferiority decision
  p2 <- plot2(dtmp2, fnfig = "fig/n_cmb1.png"); #print(p2)
  # Summary tables
  dtmp <- merge(dtmp1, dtmp2, by = c("pars_id", "brand1","es"))
  k1 <- knitr::kable(dtmp[
    brand1 == "AZ", .(es, mu_pwr,
                      mu_n = round(mu_n),
                      sd_n = round(sd_n))][
                        order(es)], digits = 3)
  k2 <- knitr::kable(dtmp[
    brand1 == "PF", .(es, mu_pwr,
                      mu_n = round(mu_n),
                      sd_n = round(sd_n))][
                        order(es)], digits = 3)
  list(p1 = p1, p2 = p2, dtmp1 = dtmp1, dtmp2 = dtmp2, k1 = k1, k2 = k2)
}

# Inferiority by brand and cohort:
# Power by es
# Total enrolled participants by es
inferiority2 <- function(l){

  dref <- ref_grid(l, 2)
  sort(unique(l$dd1$enrolled))

  # What/.when decisions made?
  dd <- decisions3(l)

  # power = proportion of simulations (by brand) that we would declare
  # inferiority at the nominated decision threshold
  dtmp1 <- dd[inf == 1, .(mu_pwr = length(unique(sim_id))/max(dd$sim_id) ),
              keyby = .(pars_id, brand1, cohort1, es)]
  dtmp1 <- merge(dref, dtmp1, by = c("pars_id", "brand1", "cohort1","es"), all.x = TRUE)
  dtmp1[is.na(mu_pwr), mu_pwr := 0]
  # sample size
  dtmp2 <- dd[inf == 1, .(mu_n = mean(enrolled),
                          sd_n = sd(enrolled)),
              keyby = .(pars_id, brand1, cohort1, es)]
  dtmp2 <- merge(dref, dtmp2, by = c("pars_id", "brand1", "cohort1", "es"), all.x = TRUE)

  # Power
  p1 <- plot1(dtmp1, fnfig = "fig/pwr_inf2.png",
              cohort = TRUE); print(p1)
  # Total participants enrolled when we hit an inferiority decision
  p2 <- plot2(dtmp2, fnfig = "fig/n_inf2.png",
              cohort = TRUE); print(p2)

  # Summary tables
  dtmp <- merge(dtmp1, dtmp2, by = c("pars_id", "brand1", "cohort1", "es"))
  k1 <- knitr::kable(dtmp[
    brand1 == "AZ" & cohort1 == 1, .(es, mu_pwr,
                      mu_n = round(mu_n),
                      sd_n = round(sd_n))][
                        order(es)], digits = 3)

  k2 <- knitr::kable(dtmp[
    brand1 == "PF" & cohort1 == 1, .(es, mu_pwr,
                      mu_n = round(mu_n),
                      sd_n = round(sd_n))][
                        order(es)], digits = 3)

  list(p1 = p1, p2 = p2, dtmp1 = dtmp1, dtmp2 = dtmp2, k1 = k1, k2 = k2)
}

non_inferiority2 <- function(l){
  dref <- ref_grid(l, 2)
  sort(unique(l$dd1$enrolled))

  # What/.when decisions made?
  dd <- decisions3(l)

  # power = proportion of simulations (by brand) that we would declare
  # inferiority at the nominated decision threshold
  dtmp1 <- dd[ninf == 1, .(mu_pwr = length(unique(sim_id))/max(dd$sim_id) ),
              keyby = .(pars_id, brand1, cohort1, es)]
  dtmp1 <- merge(dref, dtmp1, by = c("pars_id", "brand1", "cohort1", "es"), all.x = TRUE)
  dtmp1[is.na(mu_pwr), mu_pwr := 0]
  # Total participants enrolled
  dtmp2 <- dd[ninf == 1, .(mu_n = mean(enrolled),
                           sd_n = sd(enrolled)),
              keyby = .(pars_id, brand1, cohort1, es)]
  dtmp2 <- merge(dref, dtmp2, by = c("pars_id", "brand1", "cohort1", "es"), all.x = TRUE)

  # Power
  p1 <- plot1(dtmp1, fnfig = "fig/pwr_ninf2.png",
              cohort = TRUE);  print(p1)
  # Total participants enrolled when we hit an inferiority decision
  p2 <- plot2(dtmp2, fnfig = "fig/n_ninf2.png",
              cohort = TRUE);  print(p2)
  # Summary tables
  dtmp <- merge(dtmp1, dtmp2, by = c("pars_id", "brand1", "cohort1", "es"))
  k1 <- knitr::kable(dtmp[
    brand1 == "AZ" & cohort1 == 1, .(es, mu_pwr,
                      mu_n = round(mu_n),
                      sd_n = round(sd_n))][
                        order(es)], digits = 3)
  k2 <- knitr::kable(dtmp[
    brand1 == "PF" & cohort1 == 1, .(es, mu_pwr,
                      mu_n = round(mu_n),
                      sd_n = round(sd_n))][
                        order(es)], digits = 3)
  list(p1 = p1, p2 = p2, dtmp1 = dtmp1, dtmp2 = dtmp2, k1 = k1, k2 = k2)
}

combined2 <- function(l){
  dref <- ref_grid(l, 2)
  sort(unique(l$dd1$enrolled))

  dd <- decisions3(l)

  # compute power as the proportion of times that we declared inferiority
  # over the course of a trial out of the simulations for that specific
  # configuration.
  dtmp1 <- dd[anystop == 1,
              .(mu_pwr = length(unique(sim_id))/max(dd$sim_id) ),
              keyby = .(pars_id, brand1, cohort1, es)]
  dtmp1 <- merge(dref, dtmp1, by = c("pars_id", "brand1", "cohort1", "es"), all.x = TRUE)
  dtmp1[is.na(mu_pwr), mu_pwr := 0]
  # Total participants enrolled
  dtmp2 <- dd[anystop == 1, .(mu_n = mean(enrolled),
                              sd_n = sd(enrolled)),
              keyby = .(pars_id, brand1, cohort1, es)]
  dtmp2 <- merge(dref, dtmp2, by = c("pars_id", "brand1", "cohort1", "es"), all.x = TRUE)

  # Power
  p1 <- plot1(dtmp1, fnfig = "fig/pwr_cmb2.png",
              cohort = TRUE); #print(p1)
  # Total participants enrolled when we hit an inferiority decision
  p2 <- plot2(dtmp2, fnfig = "fig/n_cmb2.png",
              cohort = TRUE); #print(p2)
  # Summary tables
  dtmp <- merge(dtmp1, dtmp2, by = c("pars_id", "brand1", "cohort1", "es"))
  k1 <- knitr::kable(dtmp[
    brand1 == "AZ" & cohort1 == 1, .(es, mu_pwr,
                      mu_n = round(mu_n),
                      sd_n = round(sd_n))][
                        order(es)], digits = 3)
  k2 <- knitr::kable(dtmp[
    brand1 == "PF" & cohort1 == 1, .(es, mu_pwr,
                      mu_n = round(mu_n),
                      sd_n = round(sd_n))][
                        order(es)], digits = 3)
  list(p1 = p1, p2 = p2, dtmp1 = dtmp1, dtmp2 = dtmp2, k1 = k1, k2 = k2)
}


cumprobstop2 <- function(l){
  dref <- ref_grid(l, 2)
  sort(unique(l$dd1$enrolled))

  dd <- decisions3(l)
  dtmp1 <- dd[anystop == 1, .(pr_stop = .N/max(dd$sim_id) ),
              by = .(pars_id, enrolled, es)][order(pars_id, enrolled)]

  dtmp1[, c_pr_stop := cumsum(pr_stop), by = pars_id]

}


main <- function(){
  l <- process()
  li <- inferiority2(l)
  ln <- non_inferiority2(l)
  lc <- combined2(l)

  li$k1; li$k2
  ln$k1; ln$k2
  lc$k1; lc$k2
}


