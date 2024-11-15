## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message=FALSE-----------------------------------------------------
library(habtools)
library(raster)

## ----fig.width=4, fig.height=4------------------------------------------------
# simulate fractal terrain
surf <- sim_dem(L = 64, smoothness = 0.5)

# fractal dimension using height variation method
fd(surf, method="hvar", plot = TRUE, diagnose = TRUE)

## ----fig.width=4, fig.height=4------------------------------------------------
# height variation
fd(surf, method = "hvar")

# standard deviation
fd(surf, method = "sd")

# area
fd(surf, method = "area")

# cube counting
fd(surf, method = "cubes")

## ----fig.width=4, fig.height=4------------------------------------------------
surf_z <- surf
values(surf_z) <- values(surf_z) * ((64-1) / hr(surf_z))
hr(surf_z)

# height variation
fd(surf_z, method = "hvar")

# standard deviation
fd(surf_z, method = "sd")

# area
fd(surf_z, method = "area")

# cube counting
fd(surf_z, method = "cubes")


## ----fig.width=4, fig.height=4------------------------------------------------
# cube counting
fd(mcap, method = "cubes", plot=TRUE, diagnose=TRUE)


## ----fig.width=4, fig.height=4------------------------------------------------
# project coral as xy coordinates
mcap_2d <- mesh_to_2d(mcap)

# box counting
fd(mcap_2d, method = "boxes", plot=TRUE, diagnose=TRUE)


## ----fig.width=4, fig.height=4------------------------------------------------
dem1 <- dem_crop(horseshoe, x0 = -470.8104, y0 = 1270.625, L = 2, plot = TRUE)
drop1 <- detect_drop(dem1, d = 0.1)
plot(drop1)
# This DEM does not have many drops
fd(dem1, method = "hvar", lvec = c(1, 0.5, 0.25, 0.125), diagnose=TRUE)

dem2 <- dem_crop(horseshoe, x0 = -466.8104, y0 = 1266.625, L = 2, plot = TRUE)
drop2 <- detect_drop(dem2, d = 0.1)
plot(drop2)
fd(dem2, method = "hvar", lvec = c(1, 0.5, 0.25, 0.125), diagnose=TRUE)

