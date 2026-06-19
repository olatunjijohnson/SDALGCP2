#' Simulated aggregated disease-count data
#'
#' A small, self-contained example dataset used throughout the help pages and
#' vignettes. It is simulated from the model the package fits: an 8x8 lattice of
#' regions, a spatially structured covariate, a latent Gaussian spatial field
#' with exponential covariance, and Poisson counts with a population offset. The
#' true fixed effects are \code{(Intercept) = -6} and \code{x1 = 0.6}; the latent
#' field has variance \eqn{\sigma^2 = 0.3} and exponential scale \eqn{\phi = 4}.
#'
#' @format An \code{\link[sf]{sf}} object of 64 \code{POLYGON} regions with columns:
#' \describe{
#'   \item{region}{integer region identifier (1-64).}
#'   \item{cases}{observed disease count in the region.}
#'   \item{x1}{a standardised, spatially structured covariate.}
#'   \item{pop}{population at risk (the offset; use \code{offset(log(pop))}).}
#'   \item{geometry}{the region polygon.}
#' }
#' @source Simulated; see \code{data-raw/sdalgcp_data.R} in the package sources.
#' @seealso \code{\link{liver}} for a real disease-count example.
#' @examples
#' data(sdalgcp_data)
#' summary(sdalgcp_data$cases)
#' plot(sdalgcp_data["cases"])
"sdalgcp_data"

#' Primary biliary cirrhosis incidence in North East England
#'
#' A real aggregated disease-count dataset: incident primary biliary cirrhosis
#' (a chronic liver disease) cases by Lower-layer Super Output Area (LSOA) in the
#' Newcastle and Gateshead area of North East England, with population and area
#' deprivation covariates. This is the case study of Johnson et al. (2019) and a
#' realistic test bed for the spatial model: \code{cases ~ deprivation +
#' offset(log(pop))}.
#'
#' @format An \code{\link[sf]{sf}} object of 545 LSOA polygons
#'   (British National Grid, EPSG:27700) with columns:
#' \describe{
#'   \item{lsoa}{LSOA 2004 census code.}
#'   \item{cases}{observed incident case count in the LSOA.}
#'   \item{pop}{population at risk (the offset; use \code{offset(log(pop))}).}
#'   \item{IMD}{Index of Multiple Deprivation score (higher = more deprived).}
#'   \item{Income}{income-deprivation score.}
#'   \item{Employment}{employment-deprivation score.}
#'   \item{geometry}{the LSOA polygon.}
#' }
#' @source Johnson, O., Diggle, P. and Giorgi, E. (2019), "A spatially discrete
#'   approximation to log-Gaussian Cox processes for modelling aggregated disease
#'   count data", \emph{Statistics in Medicine}, 38(24), 4871-4884.
#'   \doi{10.1002/sim.8339}. Population and area-deprivation covariates are from
#'   the 2004 English indices of deprivation (Lower-layer Super Output Area level).
#'   See \code{data-raw/liver.R} in the package sources.
#' @seealso \code{\link{sdalgcp_data}} for a small simulated example.
#' @examples
#' data(liver)
#' summary(liver$cases)
#' plot(liver["IMD"])
"liver"
