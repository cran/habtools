#' Calculate rugosity, fractal dimension, and height for a DEM
#'
#' @param data A dem of class RasterLayer.
#' @param lvec Scales to use for calculation.
#' @param method_fd method for the calculation of rugosity and fractal dimension.
#' Can be "hvar", "sd", "cubes", or "area". Defaults to "hvar".
#' @param method_rg Method to be used for the rugosity calculation. Defaults to "area".
#' @param parallel Logical. Use parallel processing? Defaults to FALSE.
#' @param ncores Number of cores to use if parallel = TRUE.
#' @param ... Additional arguments see [fd()].
#'
#' @seealso [fd()]
#' @seealso [rg()]
#' @seealso [hr()]
#'
#' @details Uses area method for rugosity and hvar method for fractal dimension calculations as default.
#' @return A dataframe with the three complexity metrics.
#' @export
#'
#' @examples
#'
#' dem <- dem_sample(horseshoe, L = 1)
#' rdh(dem, lvec = c(0.125, 0.25, 0.5, 1))
rdh <- function(data, lvec, method_fd = "hvar", method_rg = "area",
                parallel = FALSE, ncores = (parallel::detectCores()-1),  ...){

  message(paste0("fd calculation using ", method_fd, " method."))

  if (method_fd %in% c("sd", "hvar")) {
    ddata <- fd(data, method = method_fd, lvec = lvec, keep_data = TRUE,
                parallel = parallel, ncores = ncores, ...)
    d <- ddata$D
  } else if (method_fd %in% c("cubes", "area")) {
    ddata <- fd(data, method = method_fd, lvec = lvec, keep_data = TRUE,
                 ...)
    d <- ddata$D
  } else {
    stop("Invalid method_fd specification.")
  }

  message(paste0("rg calculation using ", method_rg, " method."))

  if (method_rg == "area") {
    r <- rg(data, method = "area")
  } else if (method_rg == "hvar") {
    if (method_fd == "area") {
      r <- rg(data, method = "hvar", parallel = parallel, ncores = ncores)
    } else {
      message(paste0("L0 is set to ", min(ddata$lvec), "."))
      r <- rg(data, method = "hvar", L0 = min(ddata$lvec), parallel = parallel, ncores = ncores)
    }
  } else {
    stop("Invalid method_rg specification.")
  }

  h <- hr(data)

  data.frame(R = r, D = d, H = h)
}
