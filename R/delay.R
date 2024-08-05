# mkuhn, 2021-04-06
# delay distribution functions


#' Delayed Exponential Distribution
#'
#' @description
#' Density, distribution function, quantile function and random generation for the delayed exponential distribution with `rate`-parameter.
#'
#' @details
#' Additional arguments are forwarded via `...` to the underlying functions of the exponential distribution in the `stats`-package.
#' The numerical arguments other than `n` are recycled to the length of the result. Only the first elements of the logical arguments are used.
#'
#' @param x A numeric vector of values for which to get the density.
#' @param q A numeric vector of quantile values.
#' @param p A numeric vector of probabilities.
#' @param n integer. Number of random observations requested.
#' @param delay numeric. The delay, must be non-negative.
#' @param rate numeric. The event rate, must be non-negative.
#' @param ... further arguments are passed on to the underlying non-delayed function, e.g., [stats::dexp()]
#' @return `dexp_delayed` gives the density, `pexp_delayed` gives the distribution function, `qexp_delayed` gives the quantile function,
#' and `rexp_delayed` generates a pseudo-random sample from the delayed exponential distribution.
#'
#' The length of the result is determined by `n` for `rexp_delayed`, and is the maximum of the lengths of the numerical arguments for the other functions.
#' @keywords distribution
#' @name DelayedExponential
NULL

#' @rdname DelayedExponential
#' @export
dexp_delayed <- function(x, delay, rate = 1, ...) stats::dexp(x = x - delay, rate = rate, ...)
#' @rdname DelayedExponential
#' @export
pexp_delayed <- function(q, delay, rate = 1, ...) stats::pexp(q = q - delay, rate = rate, ...)
#' @rdname DelayedExponential
#' @export
qexp_delayed <- function(p, delay, rate = 1, ...) delay + stats::qexp(p = p, rate = rate, ...)
#' @rdname DelayedExponential
#' @export
rexp_delayed <- function(n, delay, rate = 1) delay + stats::rexp(n = n, rate = rate)


#' Delayed Weibull Distribution
#'
#' @description
#' Density, distribution function, quantile function and random generation for the delayed Weibull distribution with parameters
#' as in the Weibull distribution functions in R's stats-package, namely:
#' * `delay`
#' * `shape`
#' * `scale` (inverse of rate)
#'
#' @details
#' Additional arguments are forwarded via `...` to the underlying functions of the exponential distribution in the stats-package.
#'
#' The numerical arguments other than `n` are recycled to the length of the result. Only the first elements of the logical arguments are used.
#'
#' @param x A numeric vector of values for which to get the density.
#' @param q A numeric vector of quantile values.
#' @param p A numeric vector of probabilities.
#' @param n integer. Number of random observations requested.
#' @param delay numeric. The delay, must be non-negative.
#' @param shape numeric. Shape parameter, must be positive.
#' @param scale numeric. Scale parameter (inverse of rate), must be positive.
#' @param ... further arguments are passed on to the underlying non-delayed function, e.g., [stats::dweibull()]
#' @return `dweib_delayed` gives the density, `pweib_delayed` gives the distribution function, `qweib_delayed` gives the quantile function,
#' and `rweib_delayed` generates a pseudo-random sample from the delayed Weibull distribution.
#'
#' The length of the result is determined by `n` for `rweib_delayed`, and is the maximum of the lengths of the numerical arguments for the other functions.
#' @keywords distribution
#' @name DelayedWeibull
NULL

#' @rdname DelayedWeibull
#' @export
dweib_delayed <- function(x, delay, shape, scale = 1, ...) stats::dweibull(x = x - delay, shape = shape, scale = scale, ...)
#' @rdname DelayedWeibull
#' @export
pweib_delayed <- function(q, delay, shape, scale = 1, ...) stats::pweibull(q = q - delay, shape = shape, scale = scale, ...)
#' @rdname DelayedWeibull
#' @export
qweib_delayed <- function(p, delay, shape, scale = 1, ...) delay + stats::qweibull(p = p, shape = shape, scale = scale, ...)
#' @rdname DelayedWeibull
#' @export
rweib_delayed <- function(n, delay, shape, scale = 1) delay + stats::rweibull(n = n, shape = shape, scale = scale)


#' Get delay distribution function
#' @param distribution character(1). delay distribution.
#' @param type character(1). type of function, cdf: cumulative distribution function, density or random function
#' @param twoGroup logical(1). Do we have two groups?
#' @param bind character. Names of parameters that are bind between the two groups.
#' @return selected distribution function or parameter names
getDist <- function(distribution = c("exponential", "weibull"), type = c("cdf", "prob", "density", "random", "param"),
                    twoGroup = FALSE, bind = NULL) {
  distribution <- match.arg(distribution)
  type <- match.arg(type)

  switch(type,
         # cumulative distribution function
         prob = ,
         cdf     = c(pexp_delayed, pweib_delayed),
         # density function
         density = c(dexp_delayed, dweib_delayed),
         # random draw function
         random  = c(rexp_delayed, rweib_delayed),
         param   = {
           par_list <- list(exponential = c("delay", "rate"),
                            weibull = c("delay", "shape", "scale"))
           if (twoGroup) {
             # bind only parameters from chosen distribution
             myPars <- par_list[[1L + (distribution == 'weibull')]]
             bind <- intersect(myPars, bind) #intersect: enforces original order
             par_gr <- purrr::map(par_list, ~ setdiff(., bind))
             # bind parameters first
             purrr::map(par_gr, ~ c(bind,
                                    paste(rep(., 2L), rep(c("x", "y"), each = length(.)),
                                          sep = ".")))
           } else par_list

         },
         stop(glue("Unknown attribute of distribution {distribution}."), call. = FALSE)
  )[[1L + (distribution == 'weibull')]]
}

