# Your Name (s1801322, JianqiaoSun)
# Place your function definitions and associated code documentation my_code.R

#' Simulated confidence intervals
#' @param m number of simulations
#' @param n number of observations used for each sample
#' @param lambda true parameter value
#' @param alpha The nominal (desired) error probability for the
#'   confidence interval.
#'   Default: 0.05, to get 95% confidence intervals
#' @param type 1, 2, or 3, selecting one of the three parameterisation
#'   alternatives, theta=lambda, theta=sqrt(lambda), or theta=log(lambda).
#'   Default: 1
#'
#' @return A data frame with m rows and two columns, named Lower
#'   and Upper containing the confidence intervals, one for each row.

multi_pois_CI <- function(m, n, lambda, alpha = 0.05, type = 1) {
  Lower = c()
  Upper = c()
  for (i in 1:m) {
    simul_data = rpois(n, lambda)
    Lower[i] = pois_CI(simul_data, alpha, type)[1]
    Upper[i] = pois_CI(simul_data, alpha, type)[2]
  }
  return(data.frame(Lower, Upper))
}

#' Simulated confidence intervals, a tidy version containing all 3 types
#' @param m number of simulations
#' @param n number of observations used for each sample
#' @param lambda true parameter value
#' @param alpha The nominal (desired) error probability for the
#'   confidence interval.
#'   Default: 0.05, to get 95% confidence intervals
#'
#' @return A data frame with 3*m rows and two columns, named Lower,
#'   Upper, and Type containing the confidence intervals from three calls to multi_pois_CI.

tidy_multi_pois_CI <- function(m, n, lambda, alpha = 0.05) {
  Lower = c()
  Upper = c()
  Type = c()
  for (i in 1:m) {
    Type[3*i-2] = 1
    Type[3*i-1] = 2
    Type[3*i] = 3
    simul_data = rpois(n, lambda)
    Lower[3*i-2] = pois_CI(simul_data, alpha, 1)[1]
    Upper[3*i-2] = pois_CI(simul_data, alpha, 1)[2]
    Lower[3*i-1] = pois_CI(simul_data, alpha, 2)[1]
    Upper[3*i-1] = pois_CI(simul_data, alpha, 2)[2]
    Lower[3*i] = pois_CI(simul_data, alpha, 3)[1]
    Upper[3*i] = pois_CI(simul_data, alpha, 3)[2]
  }
  return(data.frame(Lower, Upper, Type))
}

#' Logarithm of p(N,Y)
#' @param N a vector of length J
#' @param Y a data frame of size J×2
#'
#' @return The log-probability. -Inf if input is not valid.

log_prob_NY <- function(N, Y, xi, a, b) {
  J = length(N)
  logp = 0
  if(ncol(Y)==2 & nrow(Y)==J) {
    for (j in 1:J) {
      logp = logp + log(dgeom(N[j],xi)) + lchoose(N[j],Y[j,1]) + lchoose(N[j],Y[j,2])
    }
  } else {
    return(-Inf)
  }
  if(logp==-Inf) {
    return(-Inf)
  } else {
    return(logp+lbeta(a+sum(Y), b+2*sum(N)-sum(Y)) - lbeta(a,b))
  }
}

#' Importance sampling
#' @param K number of samples to generate
#' @param Y a data frame of size J×2
#'
#' @return A data frame with K rows and J+2 columns, each row contains
#'   J vector for the prior distributions of N, Log_Weights for
#'   re-normalised log-importance-weights, Phi for detection probability

arch_importance <- function(K, Y, xi, a, b) {
  J = nrow(Y)
  df = as.data.frame(matrix(nrow=K,ncol=J+2))
  for (j in 1:J) {
    names(df)[j] = paste0("N",j)
  }
  names(df)[J+1] = "Log_Weights"
  names(df)[J+2] = "Phi"

  for (k in 1:K) {
    N = rgeom(J, xi)
    df[k,1:J] = N
    df[k,J+1] = try(log_prob_NY(N,Y,xi,a,b) - sum(log(dgeom(N,xi))))
    df[k,J+2] = try(rbeta(1,a+sum(Y), b+2*sum(N)-sum(Y)))
  }
  df$Log_Weights = try(df$Log_Weights-max(df$Log_Weights))
  return(df)
}
