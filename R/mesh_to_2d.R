#' Transform 3D mesh into 2D outline
#'
#' @description
#' Turns a 3D Mesh file into an xy data frame.
#'
#' @param mesh A mesh3d object.
#' @param L0 (Optional) The desired DEM resolution in same units at the 3D mesh.
#' @param plot logical. Plot the output?
#' @param silent logical. Defaults to not showing warnings.
#'
#' @return A data frame.
#' @export
#'
#' @details
#' The function uses the vertices of the mesh object and projects them on the XY plane.
#' Then, only points that define the perimeter of the shape are maintained.
#'
#' @importFrom graphics polygon
#'
#' @examples
#' mcap_2d <- mesh_to_2d(mcap, plot = TRUE)
#'
#' geometry::polyarea(mcap_2d$x, mcap_2d$y) # area
#' planar(mcap)
#'
#' perimeter(mcap_2d) # perimeter
#' circularity(mcap_2d) # circularity
#' fd_boxes(mcap_2d) # fractal dimension
#'

mesh_to_2d <- function(mesh, L0 = NULL, plot = FALSE, silent = TRUE){

  res <- Rvcg::vcgMeshres(mesh)$res[[1]]
  if (missing(L0)){
    L0 <- res
    if (!silent){
      message(paste("L0 is set to mesh resolution (", L0, ")", sep = ""))
    }
  }
  if (L0 < res) {
    if (!silent){
      warning("L0 is smaller than mesh resolution")
    }
  }

  if (L0 > res) {
    mesh <- Rvcg::vcgQEdecim(mesh, edgeLength = L0, silent = TRUE)
  }

  # get normals of faces
  m <- mesh
  n <- Rvcg::vcgFaceNormals(mesh)
  # remove faces facing downward
  t <- which(n[3,]<0)
  m$it <- m$it[,-t]
  m$vb[3,] <- 0
  m <- Rvcg::vcgQEdecim(m, edgeLength = L0, silent = TRUE)
  #m <- Rvcg::vcgUniformRemesh(m, voxelSize = L0)
  x <- m$vb[1,]
  y <- m$vb[2,]
  dt <- data.frame(x=x, y=y)
  poly <- concaveman::concaveman(as.matrix(dt), concavity = 1, length_threshold = L0)
  poly <- data.frame(poly)
  names(poly) <- c("x", "y")
  if (plot) {
    plot(poly, asp=1)
    polygon(poly)
  }
  return(poly)
}
