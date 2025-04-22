prior.samp.reduced <- function(chains, nprov = 6, ab = TRUE, upsilon = TRUE) {
  l <- lapply(1:chains, function(x) {
    tl <- exp(rnorm(1, 5.67, 0.01))
    tdiff <- rgamma(1, 112, scale = 0.226) 
    th <- tl + tdiff
    s_upsilon <- exp(rnorm(35, -2.5, 0.1))
    sm1 <- 0:4
    k <- -1.05*sm1^2 + 4.22*sm1 + 4.08
    nups <- nprov*35
    if (ab) {nups <-  nups - 5}
    lst <- list(
      'rho' = rgamma(5, k, scale = 0.045),
      'HA' = rgamma(1, 5.4, scale = 0.134),
      'TL' = rnorm(1, 284, 2),
      #'TL' = tl,
      'HL' = -rgamma(1, 3.6, scale = 2.253),
      'TH' = rnorm(1, 305, 2),
      #'TH' = th,
      'HH' = rgamma(1, 7.6, scale = 3.12),
      'rho_mult' = exp(rnorm(nprov*5, 0, 0.15)),
      'HL_mult' = exp(rnorm(nprov, 0, 0.15)),
      'HA_mult' = exp(rnorm(nprov, 0, 0.15)),
      'HH_mult' = exp(rnorm(nprov, 0, 0.15)),
      'TL_mult' = exp(rnorm(nprov, 0, 0.005)),
      'TH_mult' = exp(rnorm(nprov, 0, 0.005)),
      's_eps' = exp(rnorm(nprov*5, -1.5, 0.1)),
      's_upsilon' = s_upsilon,
      'upsilon' = exp(rnorm(nups, 0, rep(s_upsilon, 6)))
    )
    if (!upsilon) {lst <- lst[1:(length(lst)-2)]}
    return(lst)
  })
  return(l)
}
