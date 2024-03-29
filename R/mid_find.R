#' Find midpoint of a DEM
#'
#' @param data A DEM in RasterLayer format.
#'
#' @return A data frame with x and y midpoints.
#' @export
#'
#' @examples
#' mid_find(horseshoe)
mid_find <- function(data){
  data.frame(
    x_mid = mean(raster::extent(data)[c(1,2)],),
    y_mid = mean(raster::extent(data)[c(3,4)])
  )
}
