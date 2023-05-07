library(dclone)
library(splines)
library(tidyverse)
library(car)
count <- dplyr::count


## Calculate TA from other parameters
ta.fun <- function(HL, HH, TL, TH) {
  c3 <- 0.0001987
  den <- c3*log(-HL/HH) + (HL/TL) - (HH/TH)
  tphi <- (HL - HH)/den
  return(tphi)
}

tempvec <- seq(0, 40, by = 0.1)

## Dev time function
calc_pred3 <- function(temp, rho25, HA, TL, HL, TH, HH, TA) {
  c3 <- 0.0001987
  tK = temp + 273
  num = rho25 * tK/TA * exp(HA/c3* (1/TA-1/tK))
  den1 = exp(HL/c3 *(1/TL-1/tK))
  den2 = exp(HH/c3 *(1/TH-1/tK))
  tau = 1/(num/(1 + den1 + den2))
  return(tau)
}

## Obtain either parametric curves or bias-reduced splines from parameter vectors
get_curves <- function(data, temp = tempvec, spline = FALSE) {
  lst <- lapply(1:nrow(data), function(r) {
    row <- data[r,]
    ta <- with(row, ta.fun(HL, HH, TL, TH))
    if (spline) {
      temps <- seq(5, 35, by = 5)
      times <- with(row, calc_pred3(temps, rho, HA, TL, HL, TH, HH, ta))
      ups <- unlist(row[,grep('upsilon.', names(row))])
      sigma <- row$s_upsilon
      u.ups <- pnorm(ups, 0, 1)
      c.ups <- qlnorm(u.ups, 0, sigma)
      times <- times*c.ups
      rates <- 1/times
      if (length(temp) != length(temps)) {
        i.spline <- interpSpline(temps, rates)
        df <- as.data.frame(predict(i.spline, temp))
        names(df) <- c('temp', 'rate')
        df$rate <- pmax(df$rate, 0)
        df$time <- 1/df$rate
      }
      else {df <- data.frame('temp' = temps,
                             'rate' = pmax(rates, 0),
                             'time' = 1/pmax(rates, 0))}
    }
    else {
      times <- with(row, calc_pred3(temp, rho, HA, TL, HL, TH, HH, ta))
      rates <- 1/times
      df <- data.frame('temp' = temp, 'rate' = rates, 'time' = times)
    }
    df$index <- r
    return(df)
  })
  df2 <- bind_rows(lst)
  return(df2)
}

scale.est <- function(x) {IQR(x)/1.349}

powertrans <- function(y, lambda = NULL) {
  if (is.null(lambda)) {lambda <- powerTransform(y)$lambda}
  if (lambda == 0) {
    return(log(y))
  }
  else {
    return((y^lambda - 1)/lambda)
  }
}

qmatch <- function(prior, post) {
  l.prior <- log(prior)
  mu.pr <- mean(l.prior)
  sig.pr <- sd(l.prior)
  mu.pt <- median(post)
  sig.pt <- scale.est(post)
  
  z <- (post - mu.pt)/sig.pt
  x <- z*sig.pr + mu.pr
  return(exp(x))
}


