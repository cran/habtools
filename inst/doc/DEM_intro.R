## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 4,
  fig.height = 4
)

## ----setup, echo = TRUE, message = FALSE--------------------------------------
library(habtools)
library(raster)
library(dplyr)
library(ggplot2)

## -----------------------------------------------------------------------------
# Load a raster
dem <- horseshoe
res(dem)
plot(dem)

## -----------------------------------------------------------------------------
dem1 <- dem_crop(horseshoe, x0 = -466, y0 = 1269, L = 2, plot = TRUE)
plot(dem1)

## -----------------------------------------------------------------------------
hr(dem1)

## ----eval=FALSE, include=FALSE------------------------------------------------
#  # Height variation method
#  rg(dem1, method = "hvar", L0 = 0.05, parallel = FALSE) # Parallel = TRUE enables parallel processing using multiple cores to speed up the calculations using the height variation method. Only use this if you have a powerful computer with at least four cores.
#  
#  # Area method
#  rg(dem1, method = "area", L0 = 0.05)

## ---- eval=FALSE--------------------------------------------------------------
#  # Height variation method
#  rg(dem1, method = "hvar", L0 = 0.05, parallel = FALSE) # Parallel = TRUE enables parallel processing using multiple cores to speed up the calculations using the height variation method. Only use this if you have a powerful computer with at least four cores.
#  #> [1] 1.6123
#  
#  # Area method
#  rg(dem1, method = "area", L0 = 0.05)
#  #> [1] 1.619947

## -----------------------------------------------------------------------------
# Height variation method
fd(dem1, method = "hvar", lvec = c(0.25, 0.5, 1, 2), plot = TRUE, diagnose = TRUE)

## -----------------------------------------------------------------------------
# Area method
fd(dem1, method = "area", lvec = c(0.03125, 0.0625, 0.125, 0.25), diagnose = TRUE)

## ----eval=FALSE, echo=FALSE---------------------------------------------------
#  rdh(dem1, lvec = c(0.125, 0.25, 0.5, 1, 2), method_fd = "hvar", method_rg = "hvar")
#  rdh(dem1, lvec = c( 0.125, 0.25, 0.5, 1, 2), method_fd = "hvar")
#  rdh(dem1, lvec = c(0.03125, 0.0625, 0.125, 0.25), method_fd = "area")

## ----eval=FALSE---------------------------------------------------------------
#  rdh(dem1, lvec = c(0.125, 0.25, 0.5, 1, 2), method_fd = "hvar", method_rg = "hvar")
#  #> fd calculation using hvar method.
#  #> rg calculation using hvar method.
#  #> L0 is set to 0.125.
#  #>          R        D         H
#  #> 1 1.552032 2.300395 0.9766939
#  rdh(dem1, lvec = c( 0.125, 0.25, 0.5, 1, 2), method_fd = "hvar")
#  #> fd calculation using hvar method.
#  #> rg calculation using area method.
#  #> L0 is set to the resolution of the raster: 0.01.
#  #>         R        D         H
#  #> 1 2.12501 2.300395 0.9766939
#  rdh(dem1, lvec = c(0.03125, 0.0625, 0.125, 0.25), method_fd = "area")
#  #> fd calculation using area method.
#  #> rg calculation using area method.
#  #> L0 is set to the resolution of the raster: 0.01.
#  #>         R        D         H
#  #> 1 2.12501 2.231573 0.9766939

## -----------------------------------------------------------------------------
dem_list <- dem_split(dem, size = 2)
# calculate one metric for all squares
sapply(dem_list, hr)

# calculate multiple metrics
data_rdh <- suppressMessages(lapply(dem_list, rdh, method_fd = "hvar", lvec = c(0.25, 0.5, 1, 2))) %>%
  bind_rows()


## -----------------------------------------------------------------------------
ggplot(data_rdh) +
  geom_point(aes(x = R, y = H, color = D, size = D)) +
  theme_classic()

## -----------------------------------------------------------------------------
dem <- dem_sample(horseshoe, L=2, plot=TRUE)
plot(dem)

