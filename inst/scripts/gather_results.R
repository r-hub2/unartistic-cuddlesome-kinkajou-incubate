#!/usr/bin/env Rscript
# mkuhn, 2022-07-07
# gather individual simulation results


# setup --------------------------------------------------------------------

library('rlang')
suppressPackageStartupMessages(library('purrr'))
suppressPackageStartupMessages(library('dplyr'))
library('tidyr')
library('glue')
suppressPackageStartupMessages(library('R.utils'))



cmdArgs <- R.utils::commandArgs(trailingOnly=TRUE,
                                asValues = TRUE,
                                excludeReserved = FALSE, excludeEnvVars = TRUE,
                                defaults = list(resultsDir=file.path(getwd(), 'results'),
                                                resultsTag='MS',
                                                type = 'test'))

if (any(c('help', 'h') %in% names(cmdArgs))){
  cat('Gather Monte-Carlo simulation results and save it as a common list.\n')
  cat('Parameter options are:\n')
  cat('  --help\t print this help\n')
  cat('  --resultsDir=\t specify the directory where to find and to put the result files. Defaults to directory "results".\n')
  cat('  --resultsTag=\tspecify a name suffix for results file. Default is "MS".\n')
  cat('  --type=\t what type of results to gather? "test" (default) or "confint"\n')
  cat('  --removeTemp\t flag to clean temporary results file after they have been saved.\n')
  quit(save = 'no')
}

myRemoveTemp <- isTRUE(any(c('r', 'removeTemp') %in% names(cmdArgs)))

myResultsDir <- cmdArgs[["resultsDir"]]
stopifnot( is.character(myResultsDir), dir.exists(myResultsDir) )

myResultsTag <- cmdArgs[["resultsTag"]]
stopifnot( is.character(myResultsTag), length(myResultsTag) == 1L, nzchar(myResultsTag) )

myType <- cmdArgs[["type"]]
stopifnot( is.character(myType), length(myType) == 1L, nzchar(myType) )
myType <- match.arg(arg = tolower(myType), choices = c("test", "confint"))

simResFileNames <- list.files(myResultsDir,
                              pattern = paste0('simRes_',myType,'_[234]\\d+.+[.]rds$'),
                              full.names=TRUE)

if (! length(simResFileNames)) {
	cat('No temporary result files found!\n')
	q(save='no')
}


# run-namespace
RUN_NS <- tidyr::crossing(L1=LETTERS, L2=LETTERS) %>%
  dplyr::transmute(L=paste0(L1, L2)) %>%
  dplyr::pull(L)


# previous results --------------------------------------------------------

# look at previous results (already saved)
resData <- NULL
indOffset <- 0L
RES_FILEN <- file.path(myResultsDir, paste0('simRes_', myType, '_', myResultsTag,'.rds'))

if (file.exists(RES_FILEN)){
	resData <- readRDS(RES_FILEN)

	if ( ! is.list(resData) || is.null(names(resData))){
		cat('\nResults file is not a list! We start over from scratch now!\n')
		resData <- NULL
	} else {
	  cat(glue('There have been {length(resData)} entries saved already!'), '\n')
	  # run numbers only for tests
	  if (myType == 'test'){
	    indOffset <- max(0L,
	                     which(RUN_NS %in% purrr::map_chr(resData, ~ { substr(.x[['run']][[1L]], start = 1L, stop = nchar(RUN_NS[[1L]])) })),
	                     na.rm = TRUE)
	    cat(glue('Offset index is {indOffset}.'), '\n')
	  }
	}
}


# read in temporary results data --------------------------------------------------

# read in temporary result files, process the results
# @return list name (from metadata) and data (unnested)
readResultFile <- function(rdsFN, ind) {
	rdsF <- readRDS(rdsFN)
	rdsFC <- comment(rdsF)
	mdList <- eval(parse(text = rdsFC))
	stopifnot( is.list(mdList), all( c('host', 'time') %in% names(mdList)) )

	resName <- paste(mdList[['host']], mdList[['time']], sep='||')

	unnestVar <- c('P', 'ci_res')[[1L+(myType == 'confint')]]
	resData <- rdsF %>%
		tidyr::unnest(cols = all_of(unnestVar))

	if (myType == 'test'){
	  stopifnot( indOffset + ind <= length(RUN_NS) )
	  resData <- resData %>%
	    dplyr::mutate(run = paste0(RUN_NS[[indOffset + ind]], run))
	}

	# pass on comment to unnested dataframe
	comment(resData) <- rdsFC

	rlang::list2(!!resName := resData)
}


resCandidates <- purrr::imap(.x = simResFileNames, .f = readResultFile) %>%
  # drop list level from map
  purrr::flatten()



# save all results data ---------------------------------------------------


if (is.null(resData)){
	cat('\nStart with fresh results from scratch!\n')
	saveRDS(resCandidates, file = RES_FILEN)
} else {
	resDuplicates <- intersect(names(resData), names(resCandidates))
	if ( length(resDuplicates) ){
		cat(glue("These temporary result dataframes are already stored in the {myResultsTag}-results file:\n  * ",
		         "{paste(resDuplicates, collapse = '\n  * ')}", .trim = FALSE), "\n\n")
		cat('Please clean up temporary results files that are already saved in result list, first!\n')
		q(save='no')
	}
	cat('\nAdd result candidates to existing result list!\n')
	saveRDS(c(resData, resCandidates), file = RES_FILEN)
}

if (myRemoveTemp){
  cat('\nAbout to remove temporary result files!\n')
  file.remove(simResFileNames)
}

cat('\n~fine~\n')
