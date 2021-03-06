#' @rawNamespace useDynLib(agTrendTMB, .registration=TRUE); useDynLib(agTrendTMB_TMBExports)
#' @importFrom stats coef cov.wt dnorm lm model.matrix optim pnorm qnorm sd 
#' @keywords internal
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL


.onAttach <- function(library, pkgname)
{
  info <-utils::packageDescription(pkgname)
  package <- info$Package
  version <- info$Version
  date <- info$Date
  packageStartupMessage(
      paste(package, version, paste("(",date, ")", sep=""), "\n")
  )
}
