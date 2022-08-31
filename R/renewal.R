## constant-strength intervention + susceptible depletion
renewal_strength <- function(R0=2.5,
                             Pt=function(t) 1,
                             N=1e6,
                             genfun=function(x) pgamma(x, 5, 1),
                             I0=1,
                             tmax=100,
                             genmax=14) {
  gen <- genfun(0:genmax+0.5)
  gen <- diff(gen)
  gen <- gen/sum(gen)
  
  tvec <- seq(0, tmax)
  
  Ivec <- rep(0, length(tvec))
  Ivec[1] <- I0
  
  Svec <- rep(0, length(tvec))
  Svec[1] <- N - I0
  
  for (i in 2:length(tvec)) {
    Ivec[i] <- R0 * Svec[i-1] * Pt(i) * sum(Ivec[max(1, i-genmax):(i-1)] * gen[min(i-1, genmax):1])/N
    Svec[i] <- Svec[i-1] - Ivec[i]
  }
  
  data.frame(
    tvec=tvec,
    Ivec=Ivec,
    Svec=Svec
  )
}

## constant-speed intervention + susceptible depletion
renewal_speed <- function(R0=2.5,
                          ht=function(t) 0,
                          N=1e6,
                          genfun=function(x) pgamma(x, 5, 1),
                          I0=1,
                          tmax=100,
                          genmax=14) {
  gen <- genfun(0:genmax+0.5)
  gen <- diff(gen)
  gen <- gen/sum(gen)
  
  tvec <- seq(0, tmax)
  
  Ivec <- rep(0, length(tvec))
  Ivec[1] <- I0
  
  Svec <- rep(0, length(tvec))
  Svec[1] <- N - I0
  
  for (i in 2:length(tvec)) {
    Ivec[i] <- R0 * Svec[i-1] * sum(rev(exp(-cumsum(rev(sapply(max(1, i-genmax):(i-1), ht))))) * Ivec[max(1, i-genmax):(i-1)] * gen[min(i-1, genmax):1])/N
    Svec[i] <- Svec[i-1] - Ivec[i]
  }
  
  data.frame(
    tvec=tvec,
    Ivec=Ivec,
    Svec=Svec
  )
}

## constant-strength and speed intervention + susceptible depletion
renewal_general <- function(R0=2.5,
                            Pt=function(t) 1,
                            ht=function(t) 0,
                            N=1e6,
                            genfun=function(x) pgamma(x, 5, 1),
                          I0=1,
                          tmax=100,
                          genmax=14) {
  gen <- genfun(0:genmax+0.5)
  gen <- diff(gen)
  gen <- gen/sum(gen)
  
  tvec <- seq(0, tmax)
  
  Ivec <- rep(0, length(tvec))
  Ivec[1] <- I0
  
  Svec <- rep(0, length(tvec))
  Svec[1] <- N - I0
  
  for (i in 2:length(tvec)) {
    Ivec[i] <- R0 * Svec[i-1] * Pt(i) * sum(rev(exp(-cumsum(rev(sapply(max(1, i-genmax):(i-1), ht))))) * Ivec[max(1, i-genmax):(i-1)] * gen[min(i-1, genmax):1])/N
    Svec[i] <- Svec[i-1] - Ivec[i]
  }
  
  data.frame(
    tvec=tvec,
    Ivec=Ivec,
    Svec=Svec
  )
}
