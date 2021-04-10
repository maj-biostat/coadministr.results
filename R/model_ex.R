library(coadministr.stanc)
library(data.table)
library(rstan)
library(pwr)
rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')
options(mc.cores = 1)


run_stan_model <- function(){
  set.seed(1)
  dtbl <- data.table(arm = rep(c("CVD+FLU", "CVD+PBO"), each = 4),
                     brand = rep(c("AZ", "AZ", "PF", "PF"), len = 8),
                     cohort = rep(1:2, len = 8),
                     y = c( 70,  60,  50,  70,  50,  40,  30,  50),
                     n = c(100, 100, 100, 100, 100, 100, 100, 100))
  dtbl[, observed := y/n]
  # Matrix of form:
  # AZ1 PF1
  # AZ2 PF2
  yflu <- array(dtbl[arm == "CVD+FLU", y], dim = c(2, 2))
  ypbo <- array(dtbl[arm == "CVD+PBO", y], dim = c(2, 2))
  nflu <- array(dtbl[arm == "CVD+FLU", n], dim = c(2, 2))
  npbo <- array(dtbl[arm == "CVD+PBO", n], dim = c(2, 2))

  ld = list(
    yflu = yflu, ypbo = ypbo, nflu = nflu, npbo = npbo,
    Ncoad = 2, Nbrand = 2, Ncohort = 2,
    a0 = 1, b0 = 1, eta0 = 0.01
  )

  l1 <- rstan::sampling(coadministr.stanc::stanc_betabinhier(), data = ld,
                    chains = 1, thin = 1,
                    iter = 2000, warmup = 1000, refresh = 0)

  m <- rstan::extract(l1, pars = c("pflu", "ppbo"))
  post <- data.table(cbind(m$pflu[, , 1], m$pflu[, , 2],
                           m$ppbo[, , 1], m$ppbo[, , 2])
  )
  colnames(post) <- c("p_cvdflu_az_1", "p_cvdflu_az_2",
                      "p_cvdflu_pf_1", "p_cvdflu_pf_2",
                      "p_cvdpbo_az_1", "p_cvdpbo_az_2",
                      "p_cvdpbo_pf_1", "p_cvdpbo_pf_2")

  # compare the means with the true values
  mu <- colMeans(post)
  dtbl[, p_hat := mu]
  dtbl[, p_hat_name := names(mu)]
  dtbl

  pwr.2p.test(h = ES.h(p1 = 0.10, p2 = 0.05), sig.level = 0.05, power = .80)



  # bind together az cohort and brand for flu vs pbo
  post1 <- post[, .(p_cvdflu_az = c(p_cvdflu_az_1, p_cvdflu_az_2),
                    p_cvdpbo_az = c(p_cvdpbo_az_1, p_cvdpbo_az_2))]
  post2 <- post[, .(p_cvdflu_pf = c(p_cvdflu_pf_1, p_cvdflu_pf_2),
                    p_cvdpbo_pf = c(p_cvdpbo_pf_1, p_cvdpbo_pf_2))]

  p <- colMeans(post1)
  p.out <- pwr.2p.test(h = ES.h(p1 = p[1], p2 = p[2]), sig.level = 0.05, power = .80)
  plot(p.out)

  # posterior means cohort and brand
  llrr <- list()
  llrr$p         <- colMeans(post)
  # averages out the cohort variability
  llrr$p_az      <- colMeans(post1)
  llrr$p_pf      <- colMeans(post2)

  # differences in prob of ae for cmb vs mono across brand and cohort
  llrr$delta <- c(mean(post[[1]] - post[[5]]),
                  mean(post[[2]] - post[[6]]),
                  mean(post[[3]] - post[[7]]),
                  mean(post[[4]] - post[[8]]))

  # difference in prob of ae for cmb vs mono across brand
  llrr$delta_az <- mean(post1[[1]] - post1[[2]])
  llrr$delta_pf <- mean(post2[[1]] - post2[[2]])


  # INFERIORITY

  # pr ae rate in cmb is higher than mono across brand and cohort
  llrr$pr_inf <- c(mean(post[[1]] - post[[5]] > 0),
                   mean(post[[2]] - post[[6]] > 0),
                   mean(post[[3]] - post[[7]] > 0),
                   mean(post[[4]] - post[[8]] > 0))

  # pr ae rate in cmb is higher than mono across brand
  llrr$pr_inf_az <- mean(post1[[1]] - post1[[2]] > 0)
  llrr$pr_inf_pf <- mean(post2[[1]] - post2[[2]] > 0)


  hist(post1[[1]] - post1[[2]])


  # NON INFERIORITY

  llrr$pr_ninf   <- c(mean(post[[1]] - post[[5]] < 0.2),
                      mean(post[[2]] - post[[6]] < 0.2),
                      mean(post[[3]] - post[[7]] < 0.2),
                      mean(post[[4]] - post[[8]] < 0.2))

  llrr$pr_ninf_az <- mean(post1[[1]] - post1[[2]] < 0.2)
  llrr$pr_ninf_pf <- mean(post2[[1]] - post2[[2]] < 0.2)
}







