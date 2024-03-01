## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 4,
  fig.height = 4
)

## ----setup, message = FALSE---------------------------------------------------
library(habtools)
library(rgl)
options(rgl.printRglwidget = TRUE)

## ---- eval=FALSE--------------------------------------------------------------
#  plot3d(mcap)

## -----------------------------------------------------------------------------
resvec <- Rvcg::vcgMeshres(mcap)[[2]] # vector of resolutions
hist(resvec)
summary(resvec)

## -----------------------------------------------------------------------------
mcap_uniform <- Rvcg::vcgUniformRemesh(mcap, silent = TRUE, multiSample = TRUE, voxelSize = min(resvec), mergeClost = TRUE)
Rvcg::vcgMeshres(mcap_uniform)[[1]]
summary(Rvcg::vcgMeshres(mcap_uniform)[[2]])

## -----------------------------------------------------------------------------
# fractal dimension
fd(mcap_uniform, method = "cubes", plot = TRUE, diagnose = TRUE)

## -----------------------------------------------------------------------------
# rugosity
rg(mcap_uniform)

# height range 
hr(mcap_uniform)


## -----------------------------------------------------------------------------
planar(mcap)
surface_area(mcap)

## -----------------------------------------------------------------------------
# convexity
convexity(mcap)

# packing
packing(mcap)

# sphericity
sphericity(mcap)

# second moment of area
sma(mcap)

# second moment of volume
smv(mcap)

# mechanical shape factor
csf(mcap)


## ---- message=FALSE-----------------------------------------------------------
dem <- mesh_to_dem(mcap_uniform, res = 0.002, fill = FALSE)
raster::plot(dem)
rg(dem, method = "area") 

## -----------------------------------------------------------------------------
pts <- mesh_to_2d(mcap_uniform)
plot(pts, asp=1)

# perimeter
perimeter(pts)

# circularity
circularity(pts)

# fractal dimension
fd(pts, method = "boxes", keep_data = FALSE, plot = TRUE)


