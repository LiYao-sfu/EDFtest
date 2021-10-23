#' EDFtest: Goodness-of-fit Tests based on Empirical Distribution Function
#'
#' \code{EDFtest} provides software for the calculation of goodness-of-fit test statistics
#' and their P-values. The three statistics computed are the Empirical Distribution
#' function statistics called Cramér-von Mises, Anderson-Darling, and Watson statistic.
#'
#' @name EDFtest-package
#' @aliases EDFtest-package
#' @docType package
#' @author
#'
#' \strong{Maintainer}: Li Yao \email{yaoliy@sfu.ca}
#'
#' Authors:
#' \itemize{
#'   \item Richard Lockhart
#'   \item Li Yao
#' }
#'
#' @references Stephens, M.A. (1974). EDF Statistics for Goodness of Fit and Some Comparisons. Journal of the American Statistical Association, Vol. 69, 730-737. Available at \url{https://doi.org/10.2307/2286009}.
#' @references Stephens, M.A. (1976). Asymptotic results for goodness-of-fit statistics with unknown parameters. Annals of Statistics, Vol. 4, 357-369. Available at \url{https://doi.org/10.1214/aos/1176343411}
#'
#' @seealso
#' Useful links:
#' \itemize{
#'   \item \url{https://github.com/LiYao-sfu/EDFtest}
#'   \item Report bugs at \url{https://github.com/LiYao-sfu/EDFtest/issues}
#' }
#'
#' @keywords package
#'
#' @import stats
#' @importFrom CompQuadForm imhof
#' @importFrom rmutil rlaplace
NULL
