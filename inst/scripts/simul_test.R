#!/usr/bin/env Rscript
# Evaluate test for differences in delayed exponential or Weibull setting
#
# Test delay parameter


# init -----
library('incubate')
# minimal version check:
#+ 0.7.6 for GOF-Pvalues for restricted & unrestricted model: e.g. gof_mo0 (was gof_mo) and gof_mo1 (new)
#+ 0.9.8 for names for P-values have changed: boot => bootstrap, gof_mo0 => moran, etc
#+ 1.1.9.9000 script is developed as part of the incubate package (not separate as part of the MS)
#+ 1.1.9.9014 ties='density' as default now also for tests
#+ 1.1.9.9016 avoid attributes, use transform() for Pearson/AD GOF tests
stopifnot( packageVersion('incubate') >= '1.1.9.9016' )
cat('incubate package version: ', toString(packageVersion('incubate')), '\n')

library('dplyr', warn.conflicts = FALSE)
library('purrr')
library('tidyr', warn.conflicts = FALSE)
library('tibble')

suppressPackageStartupMessages(library('R.utils'))

TODAY <- Sys.Date()


# command line arguments -----
cmdArgs <- R.utils::commandArgs(trailingOnly=TRUE,
                                asValues = TRUE,
                                excludeReserved = FALSE, excludeEnvVars = TRUE,
                                defaults = list(resultsDir = getwd(), dist='exponential', scenario='MS', slice=0, seed=as.integer(TODAY),
                                                chnkSize=0, workers=5, R=150, mcnrep=100))


if (any(c('help', 'h') %in% names(cmdArgs))){
  cat('Run Monte-Carlo simulations with delayed Exponential or Weibull data in a two group setting.\n')
  cat('A test for difference in delay (and sometimes delay+rate) is performed.\n')
  cat('Sample size, delay, scale and scale ratio (between the two groups) and shape use different fixed values (see code in this script).\n')
  cat('Flexible command line parameter options are:\n')
  cat('  --help\t print this help\n')
  cat('  --resultsDir=\t specify the directory where to put the result files. Defaults to the directory where Rscript is executed.\n')
  cat('  --dist=\t specify distribution that governs the data generation\n')
  cat('  --scenario=\t with respect to the delay in both groups, choose a scenario for the simulation:\n\t\t\tDELAYEQ = no difference in delay,\n\t\t\tDELAYGT = 2nd group y with bigger delay.\n\t\t\tMS = only relevant scenarios shown in manuscript (default)\n\t\t\tALL = all cases\n')
  cat('  --slice=\t if given, pick only this number of scenarios for the simulation\n')
  cat('  --seed=\t if given, set random seed at the start of the script. Default is date-dependent. \n')
  cat('  --chnkSize=\t chunk size to write out results having processed so many scenarios. Default is no chunking (=0).\n')
  cat('  --workers=\t number of parallel computations using `future.callr` and `future.apply`. The first level of parallelization is across the MC-replications for each simulation setting.\n')
  cat('  --R=\t\t number of samples within parametric bootstrap test: it determines the resolution for our P-value, e.g.,\n\t\t R=100 will allow for P-values at per-cent resolution\n')
  cat('  --mcnrep=\t size of Monte-Carlo study: it is the number of replicated bootstrap data sets on which statistical tests are done.\n')
  cat('  --print\t show scenarios to simulate and exit.\n')
  quit(save = 'no')
}

myResultsDir <- cmdArgs[["resultsDir"]]
stopifnot( is.character(myResultsDir), dir.exists(myResultsDir),
           # check read & write permission (first octal information)
           (file.mode(myResultsDir) %>% as.character() %>% substr(1,1) %>% as.octmode() & 6) == '6')

myDist <- cmdArgs[["dist"]]
stopifnot( is.character(myDist), length(myDist) == 1L )
myDist <- match.arg(arg = myDist, choices = c("exponential", "weibull"))
isExpon <- isTRUE(myDist == 'exponential')
stopifnot( isExpon || isTRUE(myDist == 'weibull'))

myWorkers <- cmdArgs[["workers"]]
stopifnot( is.numeric(myWorkers), length(myWorkers) == 1L, myWorkers >= 1L )

myChnkSize <- cmdArgs[["chnkSize"]]
stopifnot( is.numeric(myChnkSize), length(myChnkSize) == 1L )

myR <- cmdArgs[["R"]]
stopifnot( is.numeric(myR), length(myR) == 1L, myR >= 1L )

myMCNrep <- cmdArgs[["mcnrep"]]
stopifnot( is.numeric(myMCNrep), length(myMCNrep) == 1L, myMCNrep >= 1L )

mySlice <- cmdArgs[["slice"]]
stopifnot( ! is.null(mySlice), is.numeric(mySlice), length(mySlice) == 1L )

mySeed <- cmdArgs[['seed']]
stopifnot( is.numeric(mySeed), length(mySeed) == 1L, mySeed >= 0L )

myScenario <- cmdArgs[["scenario"]]
stopifnot( ! is.null(myScenario), is.character(myScenario), length(myScenario) == 1L )
myScenario <- match.arg(arg = toupper(myScenario), choices = c('DELAYEQ', 'DELAYGT', 'MS', 'ALL'))

myPrint <- isTRUE(any(c('print', 'p') %in% names(cmdArgs)))





# set up simulation setting -----

if (mySeed > 0L) set.seed(mySeed)

# shape values according to distribution

simSetting <- tidyr::expand_grid(n_x = 8L, #c(8, 10, 12),
                                 delay_x = 5,
                                 delay_y = c(5, 7, 9, 11, 13, 15), #, 20, 100, 1000),
                                 scale_x = c(5, 10), #c(1, 2, 5),
                                 scale_ratio = c(2, 1, .5),
                                 # effectively filter for distribution
                                 shape = if (isExpon) 1 else c(.5, 2))

# avoid duplicates:
# by convention, group y is not less delayed than group x
simSetting <- simSetting %>%
  dplyr::filter(delay_y >= delay_x - 1e-11) %>%
  # use equally sized groups
  dplyr::mutate(n_y = n_x, .after = n_x)


#filter only relevant scale_ratio combinations
# for exponential:
# [G1] scale_x = 10 & scale_ratio = 1  (rate_x = .1, rate_ratio = 1)
# [G2] scale_x =  5 & scale_ratio = 1  (rate_x = .2, rate_ratio = 1)
# [G3] scale_x =  5 & scale_ratio =  2 [i.e., scale_y = 10]  (rate_x = .2, rate_ratio = .5 [i.e., rate_y = .1])
# [G4] scale_x = 10 & scale_ratio = .5 [i.e., scale_y =  5]  (rate_x = .1, rate_ratio = 2  [i.e., rate_y = .2])
# for weibull:
# dd=0, k=.5|2, scale_x = 10 & scale_ratio = 1
# dd=5, k=.5|2, scale_x = 10 & scale_ratio = 1


# this contains the cases which are needed in the manuscript
simFilterMS <- if (isExpon) {
  simSetting %>%
    dplyr::filter(dplyr::near(scale_ratio, 1) | dplyr::near(scale_x, 5) & dplyr::near(scale_ratio, 2) |
                    dplyr::near(scale_x, 10) & dplyr::near(scale_ratio, .5)) } else
      tibble::tibble(scale_x = 10, scale_ratio = 1)

# filter for target scenario!
# use all capital letters (see definition of myScenario)
switch (myScenario,
        DELAYEQ = {
          simSetting <- simSetting %>%
            dplyr::filter(dplyr::near(delay_x, delay_y))
                          # # for both distributions:
                          # dplyr::near(scale_x, 10), dplyr::near(scale_ratio, 1))
        },

        DELAYGT = {
          simSetting <- simSetting %>%
            dplyr::filter(delay_y > delay_x + 1e-11)
            #dplyr::inner_join(simFilter, by = names(simFilter))
        },

        MS = {
          simSetting <- simSetting %>%
            dplyr::inner_join(simFilterMS, by = names(simFilterMS))
        },

        ALL = { invisible(NULL) }, # keep everything! (also unused border cases)
        stop("Unknown target!")
)



if (mySlice > 0L) {
  simSetting <- simSetting %>%
    dplyr::slice(seq_len(length.out = ceiling(mySlice)))
}

if (myPrint) {
  print(knitr::kable(simSetting, format = 'pipe', digits = 2))
  cat('\nIn total,', NROW(simSetting), 'simulation scenarios.\n')
  cat('Each scenario is covered by ', myMCNrep, 'data replications.\n')
  cat('A bootstrap test has R=', myR, 'parametric bootstrap samples.\n')
  cat('Seed set initially is: ', if (mySeed>0) mySeed else '-not set-', '\n')
  cat('Results directory is set to ', myResultsDir, '\n')
  quit(save = 'no')
}

# set up parallel computing ----
if (myWorkers > 1L){
  library('future.callr')
  library('future.apply')

  future::plan(strategy = future.callr::callr, workers = myWorkers)
  # two level future
  # future::plan(list(
  #   tweak(future.callr::callr, workers = 2L),
  #   tweak(multicore, workers = 4L)
  # ))
}


# functions -----

#' Monte-Carlo simulation for a single simulation setting.
#' Uses parallel computation (future_replicate) to go through the (=nrep) MC-simulations.
#' Each bootstrap test is also future-aware (and would pick up a nested future-plan setting)
#' @param xx numeric. parameters that specify the simulation model for both groups
#' @return dataframe. P-values in the different Monte-Carlo runs.
doMCSim <- function(xx){
  n_x <- xx[[1]]
  n_y <- xx[[2]]
  delay_x <- xx[[3]]
  delay_y <- xx[[4]]

  scale_x <- xx[[5]]
  scale_ratio <- xx[[6]]
  shape <- xx[[7]]


  testList <- future.apply::future_replicate(n = myMCNrep, expr = {
    x <- y <- 1 #dummy init
    # generate data
    if (isExpon) {
      stopifnot( dplyr::near(shape, 1L) )
      x <- rexp_delayed(n = n_x, delay = delay_x, rate = 1/scale_x)
      y <- rexp_delayed(n = n_y, delay = delay_y, rate = 1/(scale_x * scale_ratio))
    } else {
      # weibull
      x <- rweib_delayed(n = n_x, delay = delay_x, scale = scale_x, shape = shape)
      y <- rweib_delayed(n = n_y, delay = delay_y, scale = scale_x * scale_ratio, shape = shape)
    }

    te_diff <- NULL
    # test_diff also uses parallel computations depending on future-settings
    try(expr = {
      te_diff <- test_diff(x = x, y = y, distribution = 'expon', param = 'delay', R = myR)
      # get bootstrap P-value for combined test: delay+rate (if the scale (=1/rate for exponential) is indeed different)
      if (scale_ratio != 1){
        te_diff2 <- test_diff(x = x, y = y, distribution = 'expon', param = c('delay', 'rate'), R = myR, type = 'bootstrap')
        # store P-value of delay+rate in original test_diff-object
        te_diff$P$bootstrap2 <- purrr::pluck(te_diff2, 'P', 'bootstrap', .default = NA_real_)
      }#fi
    }, silent = TRUE)
    te_diff
  }, simplify = FALSE, future.seed = TRUE)

  # drop NULLs
  testList <- purrr::compact(testList)

  # long dataframe of P-values
  tibble::tibble(
    run = seq(length(testList)),
    R_nominal = purrr::map_dbl(testList, "R"),
    R_effective = purrr::map_dbl(testList, list("testDist", length), .default = NA_real_),
    boot  = purrr::map_dbl(testList, list("P", "bootstrap"), .default = NA_real_),
    boot2 = purrr::map_dbl(testList, list("P", "bootstrap2"), .default = NA_real_),
    gof_moran0 = purrr::map_dbl(testList, list("P", "moran"), .default = NA_real_),
    gof_moran1 = purrr::map_dbl(testList, list("P", "moran1"), .default = NA_real_),
    gof_pearson0 = purrr::map_dbl(testList, list("P", "pearson"), .default = NA_real_),
    gof_pearson1 = purrr::map_dbl(testList, list("P", "pearson1"), .default = NA_real_),
    gof_ad0 = purrr::map_dbl(testList, list("P", "ad"), .default = NA_real_),
    gof_ad1 = purrr::map_dbl(testList, list("P", "ad1"), .default = NA_real_),
    lr = purrr::map_dbl(testList, list("P", "lr"), .default = NA_real_),
    lr_pp = purrr::map_dbl(testList, list("P", "lr_pp"), .default = NA_real_) ) %>%
    # make long format, all output in column 'P'
    tidyr::pivot_longer(cols = !run, names_to = "method", values_to = "value",
                        # drop NAs (for instance, missing boot2 P-value)
                        values_drop_na = TRUE)

}


#' run MC-simulations for each scenario sequentially (row-by-row)
#' @return dataframe resulting P-values as list column and add meta information
applyMCSims <- function(simSetDF){
  simSetDF %>%
    mutate(., P = apply(as.matrix(.), MARGIN = 1L, FUN = doMCSim))
}


# run & save ----
addMetaData <- function(da, timeTag) {
  # add comment as text
  comment(da) <- list(seed = mySeed, R = myR, mcnrep = myMCNrep, workers = myWorkers, chnkSize = myChnkSize,
                      host = Sys.info()[["nodename"]],
                      rversion = R.version.string,
                      incubate = as.character(packageVersion('incubate')),
                      date = TODAY,
                      time = timeTag) %>%
    #paste(names(.), ., sep = '=', collapse = ',')
    deparse
  da
}

DATETIME_TAG <- format(Sys.time(), format = "%Y-%m-%d-%Hh%Mm%Ss")
rdsBaseName <- paste0("simRes_test_", DATETIME_TAG)
rdsName <- file.path(myResultsDir, paste0(rdsBaseName, ".rds"))

if (myChnkSize < 1L || NROW(simSetting) <= myChnkSize){
  # no chunking
  simSetting <- applyMCSims(simSetDF = simSetting) %>%
    addMetaData(timeTag = DATETIME_TAG)
  saveRDS(simSetting, file = rdsName)
} else {
  # work in chunks
  rowIdx <- seq_len(NROW(simSetting))
  chnkNbr <- (length(rowIdx) %/% myChnkSize)+1L
  stopifnot( chnkNbr > 1L, chnkNbr <= 999999L )
  # stripe over the scenarios
  rowIdxLst <- split(rowIdx, f = rep_len(x=seq_len(chnkNbr), length.out = length(rowIdx)))
  stopifnot( length(rowIdxLst) == chnkNbr )
  for (i in seq_along(rowIdxLst)){
    simSetting_chnk <- dplyr::slice(simSetting, rowIdxLst[[i]])
    simSetting_chnk <- applyMCSims(simSetDF = simSetting_chnk)
    cat("Write out chunk ", i, "\n")
    saveRDS(simSetting_chnk %>% addMetaData(timeTag = DATETIME_TAG),
            file = file.path(myResultsDir, paste0(rdsBaseName, "_", sprintf("%06d", i), ".rds")))
  }#rof

  # merge chunked output!
  chnkFileNames <- list.files(path = myResultsDir, pattern = paste0('^', rdsBaseName, '_[[:digit:]]+[.]rds$'),
             full.names = TRUE)
  if ( length(chnkFileNames) ){
    chnkFiles <- purrr::map(.x = chnkFileNames, .f = readRDS)

    saveRDS(chnkFiles %>%
              dplyr::bind_rows() %>%
              addMetaData(timeTag = DATETIME_TAG),
            file = rdsName)

    if ( file.exists(rdsName) && (! exists('infoRDS') || ! inherits( try(infoRDS(rdsName), silent = TRUE), "try-error")) ){
      message("Removing ", length(chnkFileNames), " intermediate chunked RDS-files!")
      try(file.remove(chnkFileNames))
    } #fi remove RDS-chunk-files

  } else {
    warning("Did not find chunked RDS-output.")
  }

} #esle chunking



# example for a visualization of test results
# simSetting %>%
#   filter(dplyr::near(n, 10), dplyr::near(delay_x, 5), dplyr::near(rate_x, .1)) %>%
#   unnest(cols = P) %>%
#   ggplot(mapping = aes(x = P, col = method)) +
#   geom_freqpoly(bins = 12) + xlim(0,1) +
#   coord_trans(y = "log1p") +
#   facet_grid(rows = vars(delay_y), cols = vars(rate_ratio), labeller = label_both) +
#   labs(x = "P-value", title = "Test Results under different group effects", subtitle = "**n = 10**, delay~x~ = 5")



# teardown ----

# output the latest warnings:
cat("\n+++\nThese are warnings from the script:\n+++\n")
warnings()

if (myWorkers > 1L && isNamespaceLoaded('future')) future::plan(strategy = future::sequential)

cat('\n\n~fine~\n')

