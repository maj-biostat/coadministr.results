## -----------------------------------------------------------------------------
##
## Script name: run_sim.R
##
## Purpose of script: entry point for running a new simulation.
##
## Author: MAJ
##
## Date Created: 2021-03-01
##
## -----------------------------------------------------------------------------
##
## Update the contents of `run1` as required and then either call manually
## or add a call to `run1` and source the script from the commandline with
## Rscript.
##
## -----------------------------------------------------------------------------

TRIALDESIGN <- "DESIGN1"
RUNSIM <- T
tg_stan <- new.env()

library(colormap)
library(matrixStats)
library(data.table)
library(digest)
library(parallel)
library(coadministr)


if(TRIALDESIGN == "DESIGN1"){
  library(rstan)
  # options(mc.cores = parallel::detectCores()-1)
  rstan_options(auto_write = TRUE)
  Sys.setenv(LOCAL_CPPFLAGS = '-march=native')
  # tg_stan$model_code <- rstan::stan_model(file = "mod_rand_int2.stan",
  #                                         auto_write = TRUE)
}

# if(TRIALDESIGN == "DESIGN2"){
#   library(rjags)
#   library(runjags)
# }

options(mc.cores = 1)


dat_ex <- function(){

  glab_brand <- c("AZ", "PF")
  glab_trt <- rbind(c("CVD+FLU", "CVD+PBO"),
                            c("PBO", "FLU"))

  lb <- calc_betas()
  set.seed(1)
  d <- coadministr::get_data(J = 1000,
                beta = lb$b[[10]],
                sd_id = 0.5,
                n_per_yr = 2000,
                p_miss = 0.0)

  # dfig1 <- melt(dfig1, measure.vars = c("observed", "modelled"))
  # dfig2 <- d[timept == 1, .(y=sum(y), n=.N, mu = mean(y)), by = .(arm)][order(arm)]
  # dfig2[, arm := glab_trt[1, arm]]

  p <- ggplot(d, aes(x = timept,
                     y = pae,
                     col = glab_trt[1, arm])) +
    geom_jitter(pch = 20, width = 0.01, height = 0, alpha = 0.6) +
    theme_bw() +
    scale_y_continuous("Probability of reaction", limits = c(0, 1)) +
    scale_x_continuous("Timepoint", breaks = c(1, 2)) +
    scale_color_discrete("Treatment") +
    theme(legend.position = "bottom") +
    facet_wrap(paste0(glab_brand[brand]) ~
                 paste0("Enrol CVD dose ",
                        cohort, " (cohort ", cohort, ")") )
  ggsave("fig/dat_ex2.png", plot = p,
         width = 6, height = 4, units = "in")


}


calc_betas <- function(){

  # For all effect sizes look at equal and differential coad effects
  # by brand
  # pr_fx_coad <- data.table(
  #   expand.grid(c(0, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1,
  #                 0.125, 0.15, 0.20, 0.225, 0.25, 0.3 ),
  #               c(-0.1, 0.0, 0.1))
  #   )
  # pr_fx_coad[, fx_az := Var1]
  # pr_fx_coad[, fx_pf := Var1 + Var2]
  # pr_fx_coad[fx_pf < 0 , fx_pf := 0]
  # pr_fx_coad[, Var1 := NULL]; pr_fx_coad[, Var2 := NULL]

  pr_fx_coad <- data.table(
    expand.grid(c(0, 0.05, 0.075, 0.1, 0.125, 0.15,
                  0.20, 0.225, 0.25, 0.3 ),
                c(0.0))
  )
  pr_fx_coad[, fx_az := Var1]
  pr_fx_coad[, fx_pf := Var1 + Var2]
  pr_fx_coad[fx_pf < 0 , fx_pf := 0]
  pr_fx_coad[, Var1 := NULL]; pr_fx_coad[, Var2 := NULL]

  b <- lapply(1:nrow(pr_fx_coad), function(i){
    # currently uses the same increase in Pr AE for coadministration
    # with flu regardless of brand:
    pr_ae <- coadministr::pr_adverse_evt(
      pae_pbo = 0.02,
      pae_flu = 0.05,
      pae_coad_az = pr_fx_coad[i, fx_az],
      # will always be > 0
      pae_coad_pf = pr_fx_coad[i, fx_pf]
    )
    pr_ae
  })

  list(b = b, pr_fx_coad = pr_fx_coad)
}

run1 <- function(){

  lb <- calc_betas()
  message("Number of scenarios ", length(lb$b))

  if(.Platform$OS.type == "unix") {
    ncores <- parallel::detectCores()-2
  } else {
    ncores <- 1
  }
  message("Number of CPU cores used ", ncores)

  nsim <- 50
  checkit1 <- function(l){
    for(i in 1:length(l)){
      if(any(is.na(l[[i]]))){
        message("NA content at ", i)
      }
    }
  }

  i <- 10
  for(i in 1:nrow(lb$pr_fx_coad)){

    message("Starting simulation at ", Sys.time()); x = 1
    res <- parallel::mclapply(X=1:nsim, mc.cores = ncores  , FUN=function(x) {
      # set.seed(x) # helps with debugging
      message(Sys.time(), " rep ", x);
      ll <- tryCatch({
        coadministr::trial(
          scenario = i,
          looks = seq(500, 1000, by = 100),
          n_per_yr = 2000,
          delta_ni = 0.15,
          beta = lb$b[[i]],
          sd_id = 0.5,
          p_miss = 0.15,
          Ncoad = 2, Nbrand = 2, Ncohort = 2,
          a0 = 1, b0 = 1, eta0 = 0.01
        )
      },
      error=function(e) {
        message(" ERROR " , e);
        return(NA)
      })
      ll
    })

    message(" completed at ", Sys.time())
    fname <- paste0(stringr::str_pad(i, width = 3, pad = "0"),
                    "_results_", lb$pr_fx_coad[i, fx_az],"_", lb$pr_fx_coad[i, fx_pf], ".rds")
    checkit1(res)
    saveRDS(list(res = res),
            file.path("dat", fname))

  }
}




if(RUNSIM){
  if(TRIALDESIGN == "DESIGN1") run1()
  if(TRIALDESIGN == "DESIGN2") run2()
}
