
#' Extract the parameters for the specified group.
#'
#' The parameters of the requested group are named using the canonical parameter names of the distribution.
#'
#' For a one-group setting or when `group=NULL` it simply returns the given parameter.
#' This is an internal helper function
#' used in [coef.incubate_fit()], [bsDataStep()] and in the factory method [objFunFactory()] below.
#' @param par named parameters (as simple vector or as list)
#' @param group character. Which group to extract parameters for?
#' @param twoGroup flag. Is it a two-group setting?
#' @param oNames character. Original parameter names from distribution.
#' @param bind character. Which parameters are bind together in a two-group setting?
#' @return named vector of parameters from the relevant group
getPars <- function(par, group = "x", twoGroup, oNames, bind) {
  if ( ! twoGroup || is.null(group) ) return(par)

  stopifnot( is.character(group), nzchar(group) )

  # extract all group parameters and restore original name (e.g. remove ".x")
  par.gr <- purrr::set_names(par[grepl(pattern = paste0(".", group), x = names(par), fixed = TRUE)],
                              nm = setdiff(oNames, bind))

  # restore original order
  # contract: bind is intersected and has right order
  # contract: bind comes first in par
  c(par.gr, par[bind])[oNames]
}


#' Factory method for objective function, either according to maximum product of spacings estimation ('MPSE')
#' or according to standard maximum likelihood estimation ('MLE0').
#'
#' Given the observed data this factory method produces an MPSE objective function implementation
#' which is the negative of the MPSE-criterion H or the negative log-likelihood for MLE.
#'
#' @details
#' From the observations, negative or infinite values are discarded. In any case, the objective function is to be minimized.
#'
#' @param x numeric. observations
#' @param y numeric. observations in second group.
#' @param method character(1). Specifies the method for which to build the objective function. Default value is `MPSE`. `MLE0` is the standard MLE-method, calculating the likelihood function as the product of density values
#' @param distribution character(1). delayed distribution family
#' @param bind character. parameter names that are bind together (i.e. equated) between both groups
#' @param ties character. How to handle ties within data of a group.
#' @param verbose integer flag. How much verbosity in output? The higher the more output. Default value is 0 which is no output.
#' @return the objective function (e.g., the negative MPSE criterion) for given choice of model parameters or `NULL` upon errors
objFunFactory <- function(x, y=NULL, method = c('MPSE', 'MLE0'), distribution = c("exponential", "weibull"), bind=NULL,
                             ties=c('density', 'equidist', 'random', 'error'), verbose = 0L) {

  # setup ----
  stopifnot( is.numeric(x), length(x) > 0L, is.null(y) || is.numeric(y) && length(y) > 0L )
  method <- match.arg(method)
  distribution <- match.arg(distribution)
  stopifnot( is.null(bind) || is.character(bind) && length(bind) >= 1L )
  ties <- match.arg(ties)


  #standard ('original') names of distribution
  oNames <- getDist(distribution, type = "param", twoGroup = FALSE, bind = NULL)
  #bind: intersect with oNames enforces the canonical order of dist-parameters!
  bind <- intersect(oNames, bind)


  # data preparation ----

  # Break ties in case of ties='break'
  # @param obs: data vector
  # @return sorted, cleaned up data vector or NULL in case of trouble
  preprocess <- function(obs) {

    if ( is.null(obs) || ! is.numeric(obs)) return(NULL)

    ind_neg <- which(obs < 0L)
    if (length(ind_neg)){
      warning("Negative values in data", deparse(substitute(obs)), "! These are dropped.", call. = FALSE)
      obs <- obs[-ind_neg]
    }# fi

    # drop NA and +/-Inf & sort
    obs <- sort(obs[is.finite(obs)])

    if (!length(obs)) {
      warning("No valid data! Only non-negative and finite real values are valid.", call. = FALSE)
      return(invisible(NULL))
    }# fi

    if (is.null(obs)) return(NULL)

    # tie break
    # || ties == 'groupedML') # groupedML not implemented yet
    if ( startsWith(method, 'MLE') || ties == 'density' ) return(obs)

    diffobs <- diff(obs)
    stopifnot( all(diffobs >= 0L) ) # i.e. sorted obs

    tiesDiffInd <- which(diffobs == 0L) # < .Machine$double.xmin

    if (length(tiesDiffInd)){
      #rl <- rle(diff(tiesDiffInd))
      if (verbose > 0L){
        #length(which(rl$values > 1L))+1L,
        cat(glue('{length(tiesDiffInd) + sum(diff(tiesDiffInd)>1L) + 1L} tied observations ',
                 'in {sum(diff(tiesDiffInd)>1L) + 1L} group(s) within data vector.\n'))
      }

      if ( ties == 'error' ) stop('Ties within data are not allowed!', call. = FALSE)

      roundOffPrecision <- estimRoundingError(obs, maxObs = 1000L)
      if (verbose > 0L){
        cat(glue("Round-off error has magnitude {roundOffPrecision}\n"))
      }


      # rounding radius can't be wider than smallest observed diff.
      # plogis to mitigate the effect of sample size: the larger the sample the more we can 'trust' the observed minimal diff
      # obs[1L] = min(obs) = diff of minimal obs with 0
      rr <- .5 * min(stats::plogis(q = length(obs), scale = 11) * diffobs[which(diffobs > 0L)],
                     # rounding precision here
                     roundOffPrecision, obs[[1L]], na.rm = TRUE)

      ## modify tied observations per group of ties
      startInd <- endInd <- 1L
      repeat {
        #proceed to end of tie-group
        while (endInd < length(tiesDiffInd) && tiesDiffInd[endInd+1L] == tiesDiffInd[endInd] + 1L) {endInd <- endInd+1L}
        #include adjacent index to complete tie-group
        obsInd <- c(tiesDiffInd[startInd:endInd], tiesDiffInd[endInd]+1L)
        stopifnot( stats::sd(obs[obsInd]) == 0L ) #check: tie-group
        obs[obsInd] <- obs[obsInd] + if (ties == 'random') {
          # sort ensures that data after tie-break is still sorted from small to large
          sort(stats::runif(n = length(obsInd), min = -rr, max = +rr)) } else {
            stopifnot( ties == 'equidist' )
            # use evenly spaced observations to break tie as proposed by Cheng (1989) on Moran test statistic
            #+They first use the ties = 'density' approach for initial estimation of parameters for Moran's statistic
            seq.int(from = -rr, to = +rr, length.out = length(obsInd))
          }
        startInd <- endInd <- endInd+1L
        if ( startInd > length(tiesDiffInd) ) break
      } #repeat

      if (verbose > 1L && length(obs) < 50L ){
        cat(glue("New data: {paste(obs, collapse = ', ')}\n"))
      }
    } #fi tiesdiff

    # we have broken all ties
    stopifnot( !any(diff(obs)==0L) )

    obs
  } #fn preprocess

  # overwrite the data vectors with pre-processed data
  if (is.null({x <- preprocess(obs = x)})) return(invisible(NULL))
  y <- preprocess(obs = y)

  # do we have two groups after pre-processing?
  twoGroup <- !is.null(y) && is.numeric(y) && length(y)


  # checks ------------------------------------------------------------------

  if (!twoGroup && !is.null(bind)) warning("bind= has a given non-null argument but it is ignored as we have only a single group!",
                                           call. = FALSE)

  if (startsWith(method, 'MLE') && twoGroup ) {
    warning('MLE fitting is currently only supported for single group setting!', call. = FALSE)
    return(invisible(NULL))
  } # MLE


  # check that there is enough data (here we also look at bind= if twoGroup)
  if (!twoGroup && length(x) < length(oNames) || twoGroup && length(x) + length(y) < 2L * length(oNames) - length(bind) && min(length(x), length(y)) < length(oNames) - length(bind)) {
    warning('Too few valid observations provided!', call. = FALSE)
    return(invisible(NULL))
  }


  # optimization arguments -----
  par_names <- getDist(distribution = distribution, type = "param",
                       twoGroup = twoGroup, bind = bind)


  # get start values for a group of observations
  # @return list with par and upper limit for delay
  getParSetting.gr <- function(obs){
    # contract: obs is sorted!

    DELAY_MIN <- 1e-9

    parV <- switch (EXPR = distribution,
                    # min(obs) = obs[1L]
                    exponential = c( max(DELAY_MIN, obs[[1L]] - 2/length(obs)),
                                     mean(obs - obs[[1L]] + 2/length(obs))**-1L ),
                    weibull = {
                      # start values from 'Weibull plot'
                      #+using the empirical distribution function
                      ## in MASS::fitdistr they simplify:
                      # lx <- log(x)
                      # m <- mean(lx)
                      # v <- var(lx)
                      # shape <- 1.2/sqrt(v)
                      # scale <- exp(m + 0.572/shape)
                      # take out extreme values for robustness (against possible outliers)
                      #+ when at least 20 observations
                      obs_f <- obs[max(1L, floor(length(obs)*.02)):ceiling(length(obs)*.95)] #assume sorted data
                      start_y <- log(-log(1-(seq_along(obs_f)-.3)/(length(obs_f)+.4)))
                      # cf. lm.fit(x = cbind(1, log(obs)), y = start_y))$coefficients
                      start_shape <- stats::cor(log(obs_f), start_y) * stats::sd(start_y) / stats::sd(log(obs_f))
                      start_scale <- exp(mean(log(obs_f) - mean(start_y) / start_shape))


                      c( max(DELAY_MIN, obs[[1L]] - 2/length(obs)),
                         start_shape,
                         start_scale )
                    },
                    stop("Unknown distribution provided!", call. = FALSE)
    )

    list(
      par = parV,
      delay_upper = max(DELAY_MIN, min(obs) - .1/length(obs), min(obs)*.99999)
    )
  }# fn getParSetting.gr


  # parameter bounds: lower & upper
  lowerVec <- upperVec <- purrr::set_names(rep_len(NA_real_, length(par_names)),
                                           nm = par_names)

  PAR_LOW <- sqrt(.Machine$double.eps)
  PAR_BOUNDS <- list(delay = c(lower = 0, upper = NA_real_),
                     rate  = c(lower = PAR_LOW, upper = +Inf),
                     shape = c(lower = PAR_LOW, upper = +Inf),
                     scale = c(lower = PAR_LOW, upper = +Inf))


  # alas, purrr::iwalk did not work for me here
  for (na in names(PAR_BOUNDS)) {
    idx <- startsWith(par_names, prefix = na)
    if (any(idx)) {
      lowerVec[idx] <- purrr::chuck(PAR_BOUNDS, na, 'lower')
      upperVec[idx] <- purrr::chuck(PAR_BOUNDS, na, 'upper')
    } #fi
  } #rof

  par0_x <- getParSetting.gr(x)
  parV <-
    if (! twoGroup) {
      upperVec['delay'] <- par0_x[['delay_upper']]
      par0_x[['par']]

    } else {
      stopifnot( twoGroup )

      # all parameters are bound
      if ( length(bind) == length(oNames) ) {
        # treat x and y as a single group for start value heuristic
        par0_xy <- getParSetting.gr(c(x,y))

        upperVec['delay'] <- par0_xy[['delay_upper']]
        par0_xy[['par']]


      } else { #twoGroup, not all params bound!

        par0_y <- getParSetting.gr(y)

        start_x <- par0_x[['par']]
        start_y <- par0_y[['par']]

        # set upper bound for delay parameter(s)!
        if ('delay' %in% bind) {
          upperVec['delay'] <- min(par0_x[['delay_upper']], par0_y[['delay_upper']])
        } else {
          upperVec['delay.x'] <- par0_x[['delay_upper']]
          upperVec['delay.y'] <- par0_y[['delay_upper']]
        } # fi


        if (distribution == 'exponential') {

          if (is.null(bind)) c(start_x, start_y) else
            if (identical(bind, 'delay')) c(min(start_x[[1L]], start_y[[1L]]), start_x[-1L], start_y[-1L]) else
              # harmonic mean of rates (go towards higher variability setting)
              if (identical(bind, 'rate')) c(2L * start_x[[2L]] * start_y[[2L]] / (start_x[[2L]] + start_y[[2L]]),
                                             start_x[-2L], start_y[-2L])


        } else {
          stopifnot( distribution == 'weibull' )

          switch( EXPR = paste(bind, collapse = '+'),
                  delay = {
                    c(min(start_x[[1L]], start_y[[1L]]), start_x[-1L], start_y[-1L])
                  },
                  shape = {
                    # geometric mean of shapes
                    c(sqrt(start_x[[2L]] * start_y[[2L]]), start_x[-2L], start_y[-2L])
                  },
                  scale = {
                    # arithmetic mean of scales (corresponds to harmonic mean of rates)
                    c((start_x[[3L]] + start_y[[3L]])/2L, start_x[-3L], start_y[-3L])
                  },
                  `delay+shape` = {
                    c(min(start_x[[1L]], start_y[[1L]]), sqrt(start_x[[2L]] * start_y[[2L]]), start_x[[3]], start_y[[3]])
                  },
                  `delay+scale` = {
                    c(min(start_x[[1L]], start_y[[1L]]), (start_x[[3L]] + start_y[[3L]])/2L, start_x[[2]], start_y[[2]])
                  },
                  `shape+scale` = {
                    c(sqrt(start_x[[2L]] * start_y[[2L]]), (start_x[[3L]] + start_y[[3L]])/2L, start_x[[1L]], start_y[[1L]])
                  },
                  # default: bind=NULL
                  {
                    stopifnot( is.null(bind) )
                    c(start_x, start_y)
                  })

        } # weibull
      } #twoGroup, not all params bound!
    } # twoGrp

  stopifnot( ! any(is.na(upperVec), is.na(lowerVec)) )

  optim_args <- list(
    par = parV,
    method = "L-BFGS-B",
    lower = lowerVec,
    upper = upperVec,
    # QQQ something like exp(trunc(log(par))) where par is the start parameters
    control = list(parscale = pmin.int(1e11, pmax.int(sqrt(.Machine$double.eps), parV)))
  )




  # objective function ----

  # log spacings:
  # calculate the differences in EDF (for given parameters in group) of adjacent observations on log scale
  # @return n+1 cumulative diffs on log-scale
  getCumDiffs <- function(pars, group) {
    # get() would (by default) start looking in the execution environment of getCumDiffs and would find the object in its parent
    # or: env = rlang::env_parent() # parent of caller env is (normally!?) the function-env
    # or: env = rlang::fn_env(getCumDiffs) # referring directly to the function environment
    obs <- rlang::env_get(env = rlang::env_parent(rlang::current_env(), n=1L), nm = group, inherit = FALSE)
    pars.gr <- getPars(pars, group = group, twoGroup = twoGroup, oNames = oNames, bind = bind)

    # calculate spacings
    # contract: data is sorted!
    cumDiffs <- diff(c(0L,
                       purrr::exec(getDist(distribution, type = "cdf"), !!! c(list(q=obs), pars.gr)),
                       1L))


    # use densFun for ties
    # we check difference of obs directly (not cumDiffs)
    #+because cumDiffs can be 0 even if obs are different, in particular for non-suitable parameters!
    ind_t <- which(diff(obs) == 0L)
    if ( length(ind_t) ){
      stopifnot( ties == 'density' ) # other tie-strategies have already dealt with ties in *preprocess*
      # increase index by 1 to get from diff(obs)-indices to cumDiffs-indices
      cumDiffs[1L+ind_t] <- purrr::exec(getDist(distribution, type = "dens"), !!! c(list(x = obs[ind_t]), pars.gr))
    } #fi

    # respect the machine's numerical lower limit
    cumDiffs[which(cumDiffs < .Machine$double.xmin)] <- .Machine$double.xmin

    log(cumDiffs)

  }# fn getCumDiffs


  # objective function like negative mean log-spacings for MPSE
  # Estimate parameters by minimizing this function.
  # `pars` the parameter vector.
  # `aggregated` a logical flag. For two group case, `FALSE` returns individual mean log cum-diffs per group
  objFun <- function(pars, aggregated = TRUE) {
    stopifnot( length(par_names) == length(pars) )
    pars <- purrr::set_names(pars, par_names)

    switch(method,
           MPSE = {
             - if (! twoGroup) {
               mean(getCumDiffs(pars, group = "x"))
             } else {
               #twoGroup:
               #the approach to first merge x and y and then do the cumDiffs, log and mean does *not* work out
               #because the parameters should be optimized within group.
               #merged data lead to frequent non-convergence or visually bad fits
               res <- c(mean(getCumDiffs(pars, group = "x")), mean(getCumDiffs(pars, group = "y")))

               if (aggregated) stats::weighted.mean(res, w = c(length(x), length(y))) else res
             }
           },
           MLE0 = {
             # MLE0 is currently only implemented for single group situation
             stopifnot( ! twoGroup )
             nObs <- length(x)
             xc <- x - pars[['delay']]
             - if (distribution == 'exponential')
               nObs * (log(pars[['rate']]) - pars[['rate']] * mean(xc)) else
                 nObs * (log(pars[['shape']]) - pars[['shape']] * log(pars[['scale']]) + (pars[['shape']]-1L) * mean(log(xc)) - mean(xc**pars[['shape']])/pars[['scale']]**pars[['shape']])
           },
           stop(glue('Objective function for method {method} is not implemented!'), call. = FALSE)
    )
  }

  # attach analytical solution for MLE
  if ( method == 'MLE0' && ! twoGroup && distribution == 'exponential' ){
    attr(objFun, which = "opt") <- list(par = c(delay = x[[1L]], rate = 1L/(mean(x) - x[[1L]])),
                                        value = length(x) * ( log(mean(x) - x[[1L]]) + 1L ),
                                        convergence = 0L,
                                        message = "analytic solution for standard MLE ('MLE0')",
                                        counts = 1L)
  }

  objFun
}



#' Fit optimal parameters according to the objective function (either MPSE or MLE0).
#'
#' The objective function carries the given data in its environment and it is to be minimized.
#' R's standard routine `stats::optim` does the numerical optimization, using numerical derivatives.
#' or the analytical solution is returned directly if available.
#' @param objFun objective function to be minimized
#' @param optim_args list of own arguments for optimization. If `NULL` it uses the default optim arguments associated to the objective function.
#' @param verbose integer that indicates the level of verboseness. Default 0 is quiet.
#' @return optimization object including a named parameter vector or `NULL` in case of errors during optimization
delay_fit <- function(objFun, optim_args = NULL, verbose = 0) {

  if (is.null(objFun)) return(invisible(NULL))
  stopifnot( is.function(objFun) )
  objFunEnv <- rlang::fn_env(objFun)

  # check if there is already a solution provided by the objective function
  optObj <- attr(objFun, which = "opt", exact = TRUE)

  if ( is.list(optObj) && all( c('par', 'value', 'convergence') %in% names(optObj)) ){
    if (verbose > 0L) cat("Using provided (analytical) solution to objective function.\n")
  } else {
    optObj <- NULL #start from scratch
    # numeric optimization
    if (verbose > 0L) message("Start with numeric optimziation of objective function.")

    if (is.null(optim_args)) optim_args <- rlang::env_get(objFunEnv, nm = "optim_args")
    stopifnot( is.list(optim_args), 'par' %in% names(optim_args) )
    # set objective function (overwrite entry 'fn' if it is already present)
    optim_args[["fn"]] <- objFun

    # optim: 1st attempt
    try({optObj <- purrr::exec(stats::optim, !!! optim_args)}, silent = TRUE)

    if (is.null(optObj)){
      if (verbose > 0L) warning(glue("{rlang::env_get(env = objFunEnv, nm = 'method')}-optimization failed during model fit!"),
                                call. = FALSE)
    } else if ( isTRUE(optObj$convergence > 0L) ){
      # do a 2nd attempt of optim in case it did not converge in the first place
      if (verbose > 1L) message("No proper convergence during 1st optimization in delay fit. Re-try with different parameter scaling.")

      # Use parameter values of non-converged fit as new start values (and adapt parscale accordingly)
      #+The objFun is to be minimized,  smaller is better!
      if ( isTRUE(is.numeric(optObj$par) && all(is.finite(optObj$par)) && optObj$value < objFun(optim_args$par)) ){
        optim_args[['par']] <- optObj$par  # purrr::assign_in(where = "par", value = optObj$par)

        if ( isTRUE("parscale" %in% names(optim_args[["control"]])) ){
          newparsc <- abs(optim_args[['par']])
          newparsc[which(newparsc < 1e-7)] <- 1e-7
          newparsc[which(newparsc > 1e8)] <- 1e8
          optim_args[['control']][['parscale']] <- newparsc
        }

        # optim: 2nd attempt
        optObj <- NULL
        try({optObj <- purrr::exec(stats::optim, !!! optim_args)}, silent = TRUE)

        if ( is.null(optObj) || isTRUE(optObj$convergence > 0L && verbose > 0L) ) warning("No proper convergence after re-try.",
                                                                                          call. = FALSE)
      }## fi rescaling for 2nd attempt
    }## fi 2nd attempt necessary?

    # set names to parameter vector
    if (! is.null(optObj)){
      stopifnot( 'par' %in% names(optObj) )
      # set canonical names for parameters
      optObj$par <- purrr::set_names(optObj$par, rlang::env_get(env = objFunEnv, nm = "par_names"))
      # save optim_args in optimization object (but w/o objective function)
      optim_args$fn <- NULL
      optObj <- append(optObj, values = list(optim_args = optim_args))
    }
  } #esle numeric optimization

  optObj
}



#' Fit a delayed Exponential or Weibull model to one or two given sample(s).
#'
#' Maximum product spacing is used to fit the parameters.
#' Numerical optimization is done by `stats::optim`.
#' @param x numeric. observations of 1st group. Can also be a list of data from two groups.
#' @param y numeric. observations from 2nd group
#' @param distribution character. Which delayed distribution is assumed? Exponential or Weibull.
#' @param method character. Which method to fit the model? 'MPSE' = maximum product of spacings estimation *or* 'MLE0' = standard maximum likelihood estimation
#' @param bind character. parameter names that are bind together in 2-group situation.
#' @param ties character. How to handle ties.
#' @param optim_args list. optimization arguments to use. Use `NULL` to use the data-dependent default values.
#' @param verbose integer. level of verboseness. Default 0 is quiet.
#' @return `incubate_fit` the delay-model fit object. Or `NULL` if optimization failed (e.g. too few observations).
#' @export
delay_model <- function(x = stop('Specify observations!', call. = FALSE), y = NULL,
                        distribution = c("exponential", "weibull"), method = c('MPSE', 'MLE0'),
                        bind=NULL, ties=c('density', 'equidist', 'random', 'error'),
                        optim_args=NULL, verbose = 0) {


  # setup -------------------------------------------------------------------

  # unpack x if it is a list of two vectors
  if (is.list(x)){
    stopifnot( length(x) == 2L )
    y <- x[[2L]]
    x <- x[[1L]]
  }

  # enforce that the first argument x= is properly instantiated
  stopifnot( !is.null(x) && is.numeric(x) && length(x) )

  distribution <- match.arg(distribution)
  method <- toupper(method)
  if (length(method) == 1L && method == 'MSE') {
    message("The method name 'MPSE' is prefered over 'MSE'!")
    method <- 'MPSE'
  }
  method <- match.arg(method)
  ties <- match.arg(ties)

  if (is.logical(verbose)) verbose <- as.numeric(verbose)
  if ( is.null(verbose) || ! is.numeric(verbose) || ! is.finite(verbose) ) verbose <- 0
  verbose <- verbose[[1L]]


  # objective function ------------------------------------------------------

  objFun <- objFunFactory(x = x, y = y, method = method, distribution = distribution,
                          bind = bind, ties = ties, verbose = verbose)
  objFunEnv <- rlang::fn_env(objFun)

  # optimise objective function
  optObj <- delay_fit(objFun, optim_args = optim_args, verbose = verbose)

  if (is.null(optObj)) return(invisible(NULL))

  twoGroup <- rlang::env_get(env = objFunEnv, nm = "twoGroup")
  # overwrite data with  pre-processed data
  x <- rlang::env_get(env = objFunEnv, nm = "x")
  y <- rlang::env_get(env = objFunEnv, nm = "y", default = NULL)

  structure(
    list(
      data = if (twoGroup) list(x = x, y = y) else x,
      distribution = distribution,
      method = method,
      twoGroup = twoGroup,
      bind = bind,
      ties = ties,
      objFun = objFun,
      par = optObj$par,
      val = optObj$value, ##objFun(optObj$par),
      optimizer = purrr::compact(optObj[c('convergence', 'message', 'counts', 'optim_args')])),
    class = "incubate_fit")
}

#' @export
print.incubate_fit <- function(x, ...){
  coe <- coef(x)
  cat(glue::glue_data(x, .sep = "\n",
                      "Fit a delayed {distribution} through {c('Maximum Product of Spacings Estimation (MPSE)', 'standard Maximum Likelihood Estimation (MLE0)')[[1L+(method=='MLE0')]]} for {c('a single group', 'two independent groups')[[1L+twoGroup]]}.",
                      "Data: {if (twoGroup) paste(lengths(data), collapse = ' and ') else length(data)} observations, ranging from {paste(signif(range(data), 4), collapse = ' to ')}",
                      "Fitted coefficients: {paste(paste('\n  ', names(coe)), signif(coe,5L), sep = ': ', collapse = ' ')}\n\n")
  )
}

#' Coefficients of a delay-model fit.
#' @param object object that is a `incubate_fit`
#' @param group character string to request the canonical parameter for one group
#' @param ... further arguments, currently not used.
#' @return named coefficient vector
#' @export
coef.incubate_fit <- function(object, group = NULL, ...){
  #stopifnot( inherits(object, "incubate_fit") )

  getPars(object[["par"]], group = group, twoGroup = object[["twoGroup"]],
          # use original parameter names of distribution
          oNames = getDist(object[["distribution"]], type = "param", twoGroup = FALSE, bind = NULL),
          # contract: bind was intersected with parameter names and, hence, has right order
          bind = object[["bind"]])
}

#' @export
summary.incubate_fit <- function(object, ...){
  print(object)
}

#' Refit an `incubate_fit`-object with specified optimization arguments.
#' If more things need to be changed use `delay_model`.
#' @param object `incubate_fit`-object
#' @param optim_args optimization arguments
#' @param verbose integer flag. Requested verbosity during `delay_fit`
#' @param ... further arguments, currently not used.
#' @return The updated fitted object of class `incubate_fit`
#' @export
update.incubate_fit <- function(object, optim_args, verbose = 0, ...){

  stopifnot( all(c("objFun", "twoGroup", "data", "par", "val", "optimizer") %in% names(object)) )

  ## fit model with given optim_args
  optObj <- delay_fit(object[["objFun"]], optim_args = optim_args, verbose = verbose)

  if (is.null(optObj)) return(invisible(NULL))

  # update all relevant fields in the list
  object[c("par", "val", "optimizer")] <- list(par = optObj$par, val = optObj$value,
                                               # drop NULLs from list (e.g. if optim_args is not present)
                                               optimizer = purrr::compact(optObj[c("convergence", "message", "counts", "optim_args")]))

  object
}

#' @export
plot.incubate_fit <- function(x, y, title, subtitle, ...){
  stopifnot( inherits(x, "incubate_fit") )
  # parameter y comes from the plot-generic. y is not used here.

  rlang::check_installed(pkg = 'ggplot2', reason = 'to get plots', version = '3.3')

  cumFun <- getDist(x[["distribution"]], type = "cdf")

  p <- grNames <- NULL

  # catch the one-group case!
  if ( isTRUE(x[["twoGroup"]]) ){
    stopifnot( is.list(x[["data"]]) )

    grNames <- names(x[["data"]])

    p <- ggplot2::ggplot(data = tibble::enframe(unlist(x[["data"]]), name = "group"),
                         mapping = ggplot2::aes_(x = ~value, col = ~substr(x = group, 1L, 1L))) +
      # add estimated delay model(s)
      purrr::map(.x = grNames,
                 .f = ~ ggplot2::geom_function(mapping = ggplot2::aes(col = .x), inherit.aes = FALSE,
                                               fun = cumFun,
                                               args = coef(x, group = .x), linetype = "dashed"))
  } else {
    grNames <- "x"
    p <- ggplot2::ggplot(data = tibble::tibble(value=x[["data"]]),
                         mapping = ggplot2::aes_(x = ~value)) +
      # add estimated delay model
      ggplot2::geom_function(inherit.aes = FALSE,
                             fun = cumFun,
                             args = coef(x, group = grNames), linetype = "dashed")
  }


  if (missing(title)) title <- glue::glue_data(x, "Fitted {distribution} delay {c('model', 'models')[1L+twoGroup]}")
  if (missing(subtitle)) subtitle <- paste(purrr::map_chr(grNames,
                                                          ~ paste(names(coef(x, group = .)), signif(coef(x, group = .), 4),
                                                                  sep = ": ", collapse = ", ")),
                                           collapse = " - ")


  p +
    ggplot2::stat_ecdf(pad=TRUE) +
    ggplot2::xlim(0L, NA) +
    ggplot2::labs(x = 'Time', y = 'Cumulative prop. of events',
                  col = 'Group',
                  title = title, subtitle = subtitle) +
    ggplot2::scale_y_reverse()
}


#' @export
simulate.incubate_fit <- function(object, nsim = 1, seed = NULL, ...){
  stopifnot(inherits(object, 'incubate_fit'))

  ranFun <- getDist(object$distribution, type = "r")
  nObs <- if (isTRUE(object$twoGroup)) lengths(object$data) else length(object$data)

  # arguments to the random function generation
  ranFunArgsX <- as.list(c(n=nObs[[1L]], coef(object, group = "x")))
  ranFunArgsY <- if (isTRUE(object$twoGroup)) as.list(c(n=nObs[[2L]], coef(object, group = "y")))

  simExpr <- if (isTRUE(object$twoGroup))
    expression(list(x=rlang::exec(ranFun, !!! ranFunArgsX),
                    y=rlang::exec(ranFun, !!! ranFunArgsY))) else
                      expression(rlang::exec(ranFun, !!! ranFunArgsX))

  if (nsim > 1000L){
    future.apply::future_replicate(n = nsim, expr = eval(simExpr), simplify = FALSE, future.seed = TRUE)
  } else {
    replicate(n = nsim, expr = eval(simExpr), simplify = FALSE)
  }
}

#' Generate bootstrap distribution of model parameters to fitted incubate model.
#'
#' Bootstrap data are here estimated coefficients from models fitted to bootstrap samples.
#' The bootstrap data is used to make bootstrap inference in the second step.
#' It is an internal function, the main entry point is [confint.incubate_fit()].
#' @param object an `incubate_fit`-object
#' @param bs_data character. Which type of bootstrap method to generate data?
#' @param R integer. Number of bootstrapped model coefficient estimates
#' @param useBoot flag. Do you want to use the boot-package? Default value is `FALSE`.
#' @param smd_factor numeric. smooth-delay factor: influence the amount of smoothing. 0 means no smoothing at all. Default is 0.25 (as was optimal in simulation for log-quantile together with log-delay-shift = 5)
#' @return bootstrap data, either as matrix or of class `boot` (depending on the `useBoot`-flag)
bsDataStep <- function(object, bs_data = c('parametric', 'ordinary'), R, useBoot = FALSE, smd_factor = 0.25) {
  bs_data <- match.arg(bs_data)
  stopifnot(is.numeric(R), length(R) == 1L, R > 1)
  R <- ceiling(R)
  useBoot <- isTRUE(useBoot)
  stopifnot( is.numeric(smd_factor), length(smd_factor) == 1L, smd_factor >= 0L )
  smooth_delay <- isTRUE(smd_factor > 0L)
  ranFun <- getDist(object$distribution, type = "r")
  dFun <- getDist(object$distribution, type = "d")
  twoGroup <- isTRUE(object$twoGroup)
  nObs <- if (twoGroup) lengths(object$data) else length(object$data)

  if (smooth_delay && bs_data != 'parametric') {
    smooth_delay <- FALSE
    smd_factor <- 0L
    # how could smooth_delay work also for ordinary bootstrap?!
    warning('Smoothing of delay is only implemented for parametric bootstrap!', call. = FALSE)
  }


  # smooth delay: sample delay values according to objective function (where delay is varied and other parameters are kept fixed) in the vicinity of the estimated delay
  # This reflects the certainty we have in the delay estimation. Low variability in data will lead to a quickly deteriorating objective function.
  # return vector of length R with delay candidate values
  getSMDCandidates <- function(group = 'x'){
    obs <- if (twoGroup) object$data[[group]] else object$data
    obs1 <- obs[[1L]]
    del_coef <- coef(object, group = group)[['delay']]

    # avoid smoothing if 1st observation or estimted delay is too close to zer0
    if ( min(obs1, del_coef) < sqrt(.Machine$double.eps) ) return(rep_len(del_coef, length.out = R))

    stopifnot( is.function(object$objFun) )

    groupIdx <- 1L + (twoGroup && group == 'y')
    coefVect <- coef(object) # objective function expects parameters for all involved groups
    del_ind <- grep('delay', names(coefVect)) # indices of coefficients that involve delay, e.g. 'delay' or 'delay.x' or 'delay.y'
    # in case of a delay per group ('delay.x' and 'delay.y') use the right one
    if (length(del_ind) > 1L) del_ind <- del_ind[[groupIdx]]


    # look at differences of first observations
    obs_d <- diff(obs[seq_len(min(23L, nObs[[groupIdx]]))])
    obs_d <- obs_d[is.finite(obs_d) & obs_d > 0L] #get rid of ties
    obs_d <- if (! length(obs_d)) .0001 else min(obs_d)

    # candidate region for delay parameters
    #+min(..) ensures that we are not too close at obs1, otherwise for MLE we have only a single point
    #+ del_coef - (obs1 - del_coef) = 2 * del_coef - obs1
    del_interv <- c(low = max(0L, min(del_coef - (obs1 - del_coef), del_coef - obs_d,
                                      obs1 - .0001, obs1 * .9999, na.rm = TRUE)),
                    high = obs1)

    #+areas for delay with high values of objective function are more likely to be sampled
    #+candidate region: symmetric around coef_del as midpoint, up to smallest observed value
    #+candidate region becomes finer sampled the broader the interval is
    #+point estimate for delay is part of sample (if lower bound is not cut to be 0, via max in from= argument)
    delayCandDF <- tibble(delay = seq.int(from = del_interv[['low']], to = del_interv[['high']],
                                          # uneven number of grid points (hence, MPSE-estimate for delay will be one of the grid points)
                                          # grid step width at most 0.005
                                          length.out = max(997L, 2L * min(ceiling(R/2), 100L*ceiling(diff(del_interv)))+1L)),
                          # fixing the parameter estimates other than delay
                          objVal = purrr::map_dbl(.x = .data[["delay"]],
                                                  .f = ~ object$objFun(pars = replace(coefVect, del_ind, .x),
                                                                       aggregated = FALSE)[[groupIdx]]))

      # we like to drop last entry (delay = 1st observation) as objective function tends to explode
      # but we have to keep last entry if it corresponds to the delay estimate (e.g., as is the case for MLE0-fitting)
    if (delayCandDF$delay[NROW(delayCandDF)] > del_coef) delayCandDF <- delayCandDF[-NROW(delayCandDF),, drop = FALSE]
    # relative change to optimal value, will be negative as objective function is minimized
    delayCandDF$objValInv <- (object$val - delayCandDF$objVal) / (object$val+.01)
    # shift upwards into non-negative area
    delayCandDF$objValInv <- delayCandDF$objValInv - min(delayCandDF$objValInv, na.rm = TRUE)
    # scale to be between 0 and 1
    delayCandDF$objValInv <- (delayCandDF$objValInv / (max(delayCandDF$objValInv, na.rm = TRUE) + .01))**(1L/(smd_factor+.01))
    delayCandDF$cumSum0 <- cumsum(delayCandDF$objValInv)
    # scale cumSum0 to 1.
    delayCandDF$cumSum <- delayCandDF$cumSum0 / max(delayCandDF$cumSum0)
    # lag-1: have it start with 0 and end with a single 1 (the last objValInv is most often 0 as)
    delayCandDF$cumSum <- c(0L, delayCandDF$cumSum[-NROW(delayCandDF)])

    # rightmost.closed = TRUE for the unlikely/impossible?! case that we draw a 1 by runif
    delayCandDF$delay[findInterval(x = stats::runif(R), vec = delayCandDF$cumSum, rightmost.closed = TRUE)]
  }

  delayCandX <- if (smooth_delay) getSMDCandidates(group = 'x')
  delayCandY <- if (smooth_delay && twoGroup) getSMDCandidates(group = 'y')

  if (useBoot) {
    stopifnot(!twoGroup) # for the time being only single group calls are supported!
    boot::boot(data = object$data,
               statistic = function(d, i) coef(delay_model(x=d[i], distribution = object$distribution,
                                                           ties = object$ties,
                                                           method = object$method, bind = object$bind)),
               sim = bs_data, mle = coef(object), R = R,
               ran.gen = function(d, coe){ # ran.gen function is only used for parametric bootstrap
                 if (smooth_delay){
                   coe[['delay']] <- delayCandX[sample.int(n = R, size = 1L)]
                 }
                 rlang::exec(ranFun, !!! as.list(c(n=nObs[[1L]], coe)))
               })

  } else {
    # own implementation: we inline data generation (simulate) and model fitting in one function
    # get coefficients from bootstrapped data
    #+(either by ordinary bootstrap of data or by parametric bootstrap)
    coefFun <- switch(bs_data,
                      ordinary = function(dummy) {
                        # draw bootstrap samples from the data
                        x <- (if (twoGroup) object$data$x else object$data)[sample.int(n = nObs[[1L]], replace = TRUE)]
                        y <- if (twoGroup) object$data$y[sample.int(n = nObs[[2L]], replace = TRUE)]

                        coef(delay_model(x=x, y=y, distribution = object$distribution,
                                         ties = object$ties,
                                         method = object$method, bind = object$bind))
                      },
                      parametric = {
                        # generate data from the fitted model
                        # for performance reasons, we 'inline' the simulate code, cf. test_diff

                        # arguments to the random function generation
                        ranFunArgsX <- as.list(c(n=nObs[[1L]], coef(object, group = "x")))
                        ranFunArgsY <- if (twoGroup) as.list(c(n=nObs[[2L]], coef(object, group = "y")))

                        function(ind) {
                          if (smooth_delay){
                            #+smooth delay according to how sure are we about the delay-estimate:
                            #+the more sure the smaller is the smoothing
                            ranFunArgsX[['delay']] <- delayCandX[ind]
                            if (twoGroup) ranFunArgsY[['delay']] <- delayCandY[ind]
                          }

                          # cf simulate (but inlined here for performance reasons)
                          x <- rlang::exec(ranFun, !!! ranFunArgsX)
                          y <- if (twoGroup) rlang::exec(ranFun, !!! ranFunArgsY)

                          coef(delay_model(x=x, y=y, distribution = object$distribution,
                                           ties = object$ties,
                                           method = object$method, bind = object$bind))
                        }
                      },
                      stop('Unkown bootstrap data generation type!')
    )

    future.apply::future_vapply(X = seq_len(R), FUN.VALUE = double(length(coef(object))),
                                FUN = coefFun, future.seed = TRUE)

    # more clear and shorter but less efficient!
    # future.apply::future_vapply(simulate(object, nsim = R), FUN.VALUE = double(length(cf)),
    #  FUN = \(d) coef(delay_model(x=d, distribution = object$distribution, ties = object$ties, method = object$method, bind = object$bind)))

  }
}

#' Confidence intervals for parameters of incubate-model fits.
#'
#' Bias-corrected bootstrap confidence limits (either quantile-based or normal-approximation based) are generated.
#' Optionally, there are also variants that use a log-transformation first.
#' At least R=1000 bootstrap replications are recommended. Default are quantile-based confidence intervals that internally use a log-transformation.
#' @param object object of class `incubate_fit`
#' @param parm character. Which parameters to get confidence interval for?
#' @param level numeric. Which is the requested confidence level for the interval? Default value is 0.95
#' @param R number of bootstrap replications. Used only if not `bs_data`-object is provided.
#' @param bs_data character or bootstrap data object. If character, it specifies which type of bootstrap is requested and the bootstrap data will be generated. Data can also be provided here directly. If missing it uses parametric bootstrap.
#' @param bs_infer character. Which type of bootstrap inference is requested to generate the confidence interval?
#' @param useBoot logical. Delegate bootstrap confint calculation to the `boot`-package?
#' @param ... further arguments, currently not used.
#' @return A matrix (or vector) with columns giving lower and upper confidence limits for each parameter.
#' @export
confint.incubate_fit <- function(object, parm, level = 0.95, R = 199L,
                                 bs_data, bs_infer = c('logquantile', 'lognormal', 'quantile', 'quantile0', 'normal', 'normal0'),
                                 useBoot=FALSE, ...){
  stopifnot(inherits(object, 'incubate_fit'))
  stopifnot(is.numeric(level), length(level) == 1L, level < 1L, level > 0L)
  stopifnot(is.numeric(R), length(R) == 1L, R > 0)
  if (missing(bs_data)) bs_data <- 'parametric'
  if (is.vector(bs_data) && is.character(bs_data)) bs_data <- match.arg(bs_data[[1L]], choices = c('parametric', 'ordinary'))
  bs_infer <- match.arg(bs_infer)
  logTransform <- isTRUE(startsWith(bs_infer, 'log'))

  twoGroup <- isTRUE(object$twoGroup)
  nObs <- if (twoGroup) lengths(object$data) else length(object$data)

  useBoot <- isTRUE(useBoot) || inherits(bs_data, 'boot')

  genBootstrapData <- is.character(bs_data) && length(bs_data == 1L) && ! is.na(bs_data) && nzchar(bs_data)
  stopifnot( genBootstrapData || useBoot && inherits(bs_data, 'boot') || is.matrix(bs_data) )


  # check if we can really use boot
  if ( useBoot &&
       (! requireNamespace("boot", quietly = TRUE) || twoGroup || ! bs_infer %in% c('normal', 'lognormal', 'quantile', 'logquantile', 'quantile0')) ) {
    warning('Using own implementation as package', sQuote('boot'), 'is not available or scenario not implemented.',
            call. = FALSE)
    useBoot <- FALSE
  }

  cf <- coef(object)
  pnames <- names(cf)
  stopifnot( is.numeric(cf), is.character(pnames), nzchar(pnames), length(cf) == length(pnames) )

  if (missing(parm)) parm <- pnames else
    if (is.numeric(parm)) parm <- pnames[parm]
  parm <- intersect(pnames, parm) # in any case

  if (is.null(parm) || ! length(parm) || any(! nzchar(parm))) {
    warning('Invalid parameter name given in argument parm=', call. = FALSE)
    return(invisible(NULL))
  }

  stopifnot( is.character(parm), length(parm) >= 1L )

  a <- (1L - level) / 2L
  a <- c(a, 1L - a)

  # if not already provided get bootstrap data (i.e. coefficients) from fitted model to bootstrapped observations
  if (genBootstrapData) {
    bs_data <- bsDataStep(object = object, bs_data = bs_data, R = R, useBoot = useBoot)
    if (R < 999) warning('Be cautious with the confidence interval(s) because the number of bootstrap samples R is rather low (R<999).',
                       call. = FALSE)
  }
  stopifnot( ! is.vector(bs_data) && ! is.character(bs_data) )
  # set R according to the provided bs_data (in particular important when both R & bs_data object are given)
  R <- if (useBoot) bs_data[['R']] else NCOL(bs_data)


  # logShift: needed only when log-transformation is requested. Start with standard value for all parameters
  logshift <- purrr::set_names(rep_len(.0001, length.out=length(pnames)), nm = pnames)
  # for delay, the transformation needs to be independent of the scale of delay, so we subtract the minimum and add a shift
  #+use fixed logshift_delay = 5 (which performed well in simulation at single group, exponential distribution, together with smd=0.25)
  if (logTransform){
    LOGSHIFT_DELAY <- 5
    for (i in which(startsWith(pnames, 'delay'))){
      logshift[i] <- -min(if (useBoot) bs_data$t[,i] else bs_data[i,], na.rm = TRUE) + LOGSHIFT_DELAY
      # using low quantiles would make it less dependent on R but then we needed to check that x-logshift remains positive (for log)
      #stats::quantile(..i.., probs = c(0, 0.001), na.rm = TRUE, names = FALSE) # catch when diff() > LOGSHIFT_DELAY
    }#rof
  }#fi

  # do bootstrap inference on bootstrap data
  ci <- if (useBoot) {
    stopifnot( inherits(bs_data, 'boot') )

    # 'perc' just takes the quantiles,
    #+'basic' uses quantiles of the difference to the observed value (bias-correction)
    ci_type <- switch(bs_infer,
                         quantile0 = 'perc',
                         quantile =,
                         logquantile = 'basic',
                         normal =,
                         lognormal = 'norm',
                         stop('This boot.ci-type is not supported!'))

    matrix(unlist(
      purrr::map(seq_len(length.out = length(coef(object))), .f = ~ {
        # the output of boot.ci can have different CIs as named matrix list entries
        ci_bo <- {if (logTransform)
          boot::boot.ci(bs_data, index = ., conf = level, type = ci_type,
                        h = function(t) log(t + logshift[[.]]), hdot = function(t) 1/(t + logshift[[.]]),
                        hinv = function(t) exp(t) - logshift[[.]]) else
            boot::boot.ci(bs_data, index = ., conf = level, type = ci_type)}[[switch(ci_type,
                                                                                       norm = 'normal',
                                                                                       perc = 'percent',
                                                                                       ci_type)]]
        # depending on the CI-type: normal yields 3 columns, perc and others give 5 columns
        stopifnot( is.matrix(ci_bo), NCOL(ci_bo) > 2L )
        # the last two columns are always the lower and upper bound
        ci_bo[, c(NCOL(ci_bo)-1L, NCOL(ci_bo))] })),
      ncol = 2L, byrow = TRUE)
  } else {

    stopifnot( is.matrix(bs_data) )

    # bootstrapped confidence limits
    # bias-correction for parametric bootstrap only!?
    #delayH_mle_bias <- mean(delay_mle_bs) - delayH_mle
    switch(bs_infer,
           quantile0 = {
             t(apply(bs_data, 1L, stats::quantile, probs = a, na.rm = TRUE))
           },
           quantile = {
             # bias-corrected quantile-based CI
             # see Davison, p28
             # vector - matrix: vector is expanded column-wise, and the row-dimension fits (=number of coefs)
             2L * cf - t(apply(bs_data, 1L, stats::quantile, probs = rev(a), na.rm = TRUE))

           },
           logquantile = local({
             # #bs_min <- apply(bs_data, 1L, min) - .15
             # bs_min <- purrr::set_names(rep.int(-.001, length(cf)), nm = names(cf))
             # # for delay, the transformation should be independent of the scale of delay
             # if ('delay' %in% names(bs_min)) bs_min['delay'] <- min(bs_data['delay',], na.rm = TRUE) - .1

             ## bias-corrected normal-based CI after log-transformation
             -logshift + exp(
               2L * log(cf + logshift) - log(t(apply(bs_data, 1L, stats::quantile, probs = rev(a), na.rm = TRUE))+logshift)
             )
           }),
           normal0 = {
             t(c(1L, 1L) %o% .rowMeans(bs_data, m = length(cf), n = R) + stats::qnorm(a) %o% apply(bs_data, 1L, stats::sd))
           },
           normal = {
             ## bias-corrected normal-based CI
             ## ci_delay_mle <- delayH_mle - delayH_mle_bias + c(-1, 1) * qnorm(.975) * delayH_mle_sd
             t(c(1L, 1L) %o% (2L * cf - .rowMeans(bs_data, m = length(cf), n = R)) + stats::qnorm(a) %o% apply(bs_data, 1L, stats::sd))
           },
           lognormal = local({
             # #bs_min <- apply(bs_data, 1L, min) - .15
             # bs_min <- purrr::set_names(rep.int(-.001, length(cf)), nm = names(cf))
             # # for delay, the transformation should be independent of the scale of delay
             # if ('delay' %in% names(bs_min)) bs_min['delay'] <- min(bs_data['delay',], na.rm = TRUE) - .1

             bs_data_h <- log(bs_data + logshift)
             ## bias-corrected normal-based CI after log-transformation
             -logshift + exp(
               t(c(1L, 1L) %o% (2L * log(cf + logshift) - .rowMeans(bs_data_h, m = length(cf), n = R)) + stats::qnorm(a) %o% apply(bs_data_h, 1L, stats::sd)))
           }),
           stop('This type of bootstrap confidence interval is not supported!')
    )
  } #esle useBoot

  # ensure formatted row and column names
  rownames(ci) <- pnames
  colnames(ci) <- paste0(format(a*100, trim = TRUE, nsmall = 1L), '%')

  # enforce parameter bounds also for CI
  # all parameters are non-negative!
  ci[which(ci<0L)] <- 0L


  ci[parm, , drop = FALSE]
}

#' Transform observed data to unit interval
#'
#' The transformation is the probability integral transform. It uses the cumulative distribution function with the estimated parameters of the model fit.
#' All available data in the model fit is transformed.
#'
#' @note
#' This S3-method implementation is quite different from its default method that allows for non-standard evaluation on data frames, primarily for interactive use.
#' But the name `transform` just fits so nicely to the intended purpose that it is re-used for the probability integral transform.
#'
#' @param _data a fitted model object of class `incubate_fit`
#' @param ... currently ignored
#' @return The transformed data, either a vector (for single group) or a list with entries x and y (in two group scenario)
#' @export
transform.incubate_fit <- function(`_data`, ...){
  stopifnot(inherits(`_data`, "incubate_fit"))

  cdfFun <- getDist(`_data`$distribution, type = "cdf")

  twoGroup <- `_data`$twoGroup
  x <- if (twoGroup) `_data`$data$x else `_data`$data

  tr <- purrr::exec(cdfFun, !!! c(list(q=x), coef(`_data`, group = 'x')))
  if (twoGroup) tr <- list(x = tr, y = purrr::exec(cdfFun, !!! c(list(q=`_data`$data$y), coef(`_data`, group = 'y'))))

  tr
}
