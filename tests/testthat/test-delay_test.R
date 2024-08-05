# test cases for significance tests for delay

test_that('Structure of test objects.', {
  set.seed(123)

  x <- rexp_delayed(n = 131L, delay = 5, rate = .1)
  y <- rexp_delayed(n = 111L, delay = 5, rate = .3)

  ted_d <- test_diff(x = x, y = y, distribution = 'expon', param = 'delay',
                     R = 19, type = 'bootstrap')

  expect_s3_class(ted_d, class = 'incubate_test')
  expect_identical(names(ted_d), c('t_obs', 'testDist', 'R', 'chisq_df_hat', 'param', 'P'))
  expect_type(ted_d$P$bootstrap, type = 'double')
  expect_lte( ted_d$P$bootstrap, expected = 1L)
  expect_gte( ted_d$P$bootstrap, expected = 0L)
})

test_that('GOF-test on single-group exponentials', {
  testthat::skip_on_cran()
  testthat::skip(message = 'Too long to run every time!')
  future::plan(future.callr::callr, workers = 3L)

  # GOF-tests on true exponential data with varying sample size, delay and rate
  # results in a matrix of dimension #scenarios x #replications
  fitting_expos <- future.apply::future_replicate(n = 467L, simplify = FALSE, expr = {
    scenarios <- expand.grid(n = c(10, 25, 50), delay = c(0, 5, 15), rate = c(.01, .2, .4, 1, 1.5, 4))
    # fit exponential models with varying n, delay and rate
    purrr::pmap(.l = scenarios,
                .f = ~ delay_model(x = rexp_delayed(n = ..1, delay = ..2, rate = ..3),
                                   distribution = 'exponential'))
  }) %>% purrr::transpose() # get a list of scenarios, each containing its models of replicated data

  # list: for each scenario, the vector of Moran's GOF-test p-value
  GOF_pvals <- list(
    moran = purrr::map(.x = fitting_expos, .f = ~ purrr::map_dbl(., function(fit) test_GOF(fit, method = 'moran')$p.value)),
    pearson = purrr::map(.x = fitting_expos, .f = ~ purrr::map_dbl(., function(fit) test_GOF(fit, method = 'pearson')$p.value))
    #ad = purrr::map(.x = fitting_expos, .f = ~ purrr::map_dbl(., function(fit) test_GOF(fit, method = 'ad')$p.value))
  )

  # # to visualize the GOF P-values:
  # as_tibble(GOF_pvals[['moran']], .name_repair = 'unique') %>%
  #   tidyr::pivot_longer(cols = everything()) %>%
  #   ggplot(aes(x=value, col = name)) +
  #   geom_freqpoly(binwidth = .15) + xlim(0, 1)

  # expect uniform P-values for GOF under valid H0
  # go over the three types of GOF-tests, within test type, check each scenario
  # use purrr::flatten(GOF_pvals) as data argument to walk to test *all* in one go
  purrr::walk(GOF_pvals[['pearson']], .f = ~expect_gt(object = mean(.), expected = 0.25))
  purrr::walk(GOF_pvals[['pearson']], .f = ~expect_lt(object = mean(.), expected = 0.75))
  # purrr::walk(GOF_pvals[['ad']], .f = ~expect_gt(object = mean(.), expected = 0.25))
  # purrr::walk(GOF_pvals[['ad']], .f = ~expect_lt(object = mean(.), expected = 0.75))

  purrr::walk(GOF_pvals[['moran']], .f = ~expect_gt(object = mean(.), expected = 0.35))
  purrr::walk(GOF_pvals[['moran']], .f = ~expect_lt(object = mean(.), expected = 0.65))
  # Moran's test shows indeed P-values on average close to 0.5
  purrr::walk(GOF_pvals[['moran']], .f = ~expect_equal(object = mean(.), expected = 0.5, tolerance = .25, info = 'moran'))
  # Moran's test does not give P-values very close to 0, less than expected under uniformity
  GOF_moran_pvals_KSpval <- purrr::map_dbl(GOF_pvals[['moran']], .f = ~suppressWarnings(ks.test(.[which(. > .085)],
                                                                                                y = 'punif', min = 0.085)$p.value))
  # not too many small P-values
  expect_lte(length(which(GOF_moran_pvals_KSpval < .05)) / length(GOF_moran_pvals_KSpval),
             expected = .3)
  # some high P-values
  expect_gte(length(which(GOF_moran_pvals_KSpval > .6)) / length(GOF_moran_pvals_KSpval),
             expected = .15)

  future::plan(future::sequential)
})


test_that("Test difference in delay for two exponential fits", {
  testthat::skip_on_cran()
  testthat::skip(message = 'Too long to run every time!')

  future::plan(future.callr::callr, workers = 3L)

  set.seed(12345)
  x <- rexp_delayed(13L, delay = 11, rate = .05)
  y <- rexp_delayed(17L, delay = 11, rate = .08)

  # increasing effect
  te_diff_delays <- purrr::map(purrr::set_names(c(0, 9, 19)),
                               ~ test_diff(x = x + .x, y = y, param = "delay", type = 'bootstrap', R = 399))

  te_diff_delays_P_bs <- purrr::map_dbl(te_diff_delays, ~ purrr::chuck(., "P", "bootstrap"))

  # null model (no effect) has a high P-value
  expect_gt(te_diff_delays_P_bs[['0']], expected = .1)

  # the bigger the effect (=difference in delay) the smaller the P-value
  expect_true(all(diff(te_diff_delays_P_bs) < 0L))
  # negative correlation betw effect size and P-values
  expect_lt(cor(x = as.integer(names(te_diff_delays)),
                y = te_diff_delays_P_bs),
            expected = -.67)

  #test effect of sample size: increase n and power goes up.
  set.seed(123456)
  # data with difference in delay by 2.5 time units
  #+but different sample sizes
  xs <- purrr::map(purrr::set_names(c(9, 20, 32, 37)),
                   ~ rexp_delayed(., delay = 6.5, rate = .07))

  ys <- purrr::map(purrr::set_names(c(10, 19, 30, 38)),
                   ~ rexp_delayed(., delay = 9, rate = .07))

  te_diff_delays_n <- purrr::map2(.x = xs, .y = ys,
                                  .f = ~ suppressWarnings(test_diff(x = .x, y = .y, param = "delay", type = 'bootstrap', R = 399)))

  expect_lt(cor(x = as.integer(names(te_diff_delays_n)),
                y = purrr::map_dbl(te_diff_delays_n, ~ purrr::chuck(., "P", "bootstrap"))),
            expected = -.67)

  future::plan(future::sequential)
})



test_that("Test difference in delay when H0 is true (no difference in delay)", {

  testthat::skip_on_cran()
  testthat::skip(message = 'Too long to run every time!')

  future::plan(future.callr::callr, workers = 3L)

  set.seed(20210506)

  testres_P_H0 <- future.apply::future_vapply(X = seq_len(21L), FUN.VALUE = double(1L),
                                              FUN = function(dummy) {
                                                x <- rexp_delayed(13, delay = 4, rate = .07)
                                                y <- rexp_delayed(11, delay = 4, rate = .1)

                                                # return P-value of bootstrap test
                                                Pval <- NA_real_
                                                try(Pval <- purrr::chuck(test_diff(x = x, y = y, param = "delay", distribution = "exp", R = 301),
                                                                         "P", "bootstrap"),
                                                    silent = TRUE)

                                                Pval
                                              }, future.seed = TRUE)

  testres_P_H0 <- testres_P_H0[is.finite(testres_P_H0)]

  # KS-test does not reject H0: uniform distribution
  expect_gt(suppressWarnings(stats::ks.test(x = testres_P_H0, y = "punif")$p.value), expected = .2)

  future::plan(future::sequential)

})

