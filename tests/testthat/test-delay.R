# mkuhn, 2022-03-16
# tests for the delayed distribution functions


test_that('density of delayed distributions', {
  tPoints <- seq.int(from = -1, to = 11, length.out = 27)

  rateVals <- c(.06, .32, .821, 1.14, 1.78, 2.19, 5.116, 11.2, rlnorm(n=3, meanlog = 1, sdlog = 3))
  shapeVals <- c(.0001, .0052, .14, .87, 1.26, 3.9, 14.8, 21, rpois(n=3, lambda = 8))

  # delayed exponential with delay=0 coincides with stats::dexp
  purrr::walk(.x = rateVals, .f = ~expect_identical(dexp_delayed(x = tPoints, delay = 0L, rate = .x),
                                                    stats::dexp(x = tPoints, rate = .x)))

  # delayed Weibull with delay=0 coincides with stats::dweibull
  purrr::walk2(.x = rateVals**-1, .y = shapeVals,
              .f = ~expect_identical(dweib_delayed(x = tPoints, delay = 0L, scale = .x, shape = .y),
                                     stats::dweibull(x = tPoints, scale = .x, shape = .y)))

  # delayed exponential coincides with delayed weibull with shape = 1 fixed
  delayTimes <- c(0, 1, 4, 7, 11, 13)
  delayTimes <- c(delayTimes, rpois(n=length(rateVals) - length(delayTimes), lambda = 9))
  purrr::walk2(.x = delayTimes, .y = rateVals,
               .f = ~expect_equal(dexp_delayed(x = tPoints, delay = .x, rate = .y),
                                  dweib_delayed(x = tPoints, delay = .x, shape = 1L, scale = .y**-1)))
})

test_that('distribution function of delayed distributions', {
  qPoints <- seq.int(from = 0, to = 1, length.out = 27)

  rateVals <- c(.1, .52, 1.4, 1.58, 3.9, 11.2)
  shapeVals <- c(.001, .14, .84, 1.6, 5.9, 13.8)
  # CDF of delayed exponential with delay=0 coincides with stats::pexp
  purrr::walk(.x = rateVals, .f = ~expect_identical(pexp_delayed(q = qPoints, delay = 0L, rate = .x),
                                                           stats::pexp(q = qPoints, rate = .x)))

  # CDF of delayed Weibull with delay=0 coincides with stats::pweibull
  purrr::walk2(.x = rateVals**-1, .y = shapeVals,
               .f = ~expect_identical(pweib_delayed(q = qPoints, delay = 0L, scale = .x, shape = .y),
                                      stats::pweibull(q = qPoints, scale = .x, shape = .y)))

  # delayed exponential coincides with delayed weibull with shape = 1 fixed
  delayTimes <- c(0, 1, 4, 7, 11, 13)
  purrr::walk2(.x = delayTimes, .y = rateVals,
               .f = ~expect_equal(pexp_delayed(q = qPoints, delay = .x, rate = .y),
                                  pweib_delayed(q = qPoints, delay = .x, shape = 1L, scale = .y**-1)))
})



