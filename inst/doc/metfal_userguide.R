## ----results = "asis",echo=FALSE,warning = FALSE,message=FALSE-----------
library(metfal)
simdpath = system.file("extdata", "simdata.csv", package = "metfal")
simd = read.csv(simdpath,header=FALSE)
colnames(simd) = NULL
pander::pandoc.table(simd[1:7, 1:9],split.table = Inf,style="rmarkdown")

## ----results = "markup",warning = FALSE----------------------------------
## loading metfal package
library(metfal)

## finding the path of simdata.csv
simdpath = system.file("extdata", "simdata.csv", package = "metfal")

## Read the metabolomics dataset
simd = read.metfile(simdpath,nrowheader=4,ncolheader=7, metcol=1,smplrow=4, resmpl="")

## ---- results="asis",echo=FALSE------------------------------------------
pander::pandoc.table(simd[1:3, 1:4])

## ---- results="markup"---------------------------------------------------
metf = makefeature(data=simd,wsize=10,ssize=0.5,defswitch=0.2)

## ---- results="markup"---------------------------------------------------
falsig = falsignal(l=metf, cutoffvar=0.03,cutoffswt=1,cutofflong=NULL)

## ---- results="asis", eval = FALSE---------------------------------------
#  write.metfile(simd,falsig,type=2,trans=FALSE)

## ---- results="asis",eval = FALSE----------------------------------------
#  plot_mrate(falsig,type="1")

## ----eval=FALSE----------------------------------------------------------
#  plot_qc(falsig)

## ----fig.width=4, fig.height=2,echo=FALSE--------------------------------
library(png)
library(grid)
img <- readPNG(system.file("extdata", "heat.png", package = "metfal"))
grid.raster(img)

## ----fig.width=4, fig.height=2,echo=FALSE--------------------------------
img <- readPNG(system.file("extdata", "var1.png", package = "metfal"))
grid.raster(img)

## ----fig.width=4, fig.height=2,echo=FALSE--------------------------------
img <- readPNG(system.file("extdata", "mm1.png", package = "metfal"))
grid.raster(img)

## ----fig.width=4, fig.height=2,echo=FALSE--------------------------------
img <- readPNG(system.file("extdata", "lb1.png", package = "metfal"))
grid.raster(img)

## ----fig.width=4, fig.height=2,echo=FALSE--------------------------------
img <- readPNG(system.file("extdata", "swi1.png", package = "metfal"))
grid.raster(img)

## ----eval = FALSE--------------------------------------------------------
#  ### find the path of simulated metabolomics data
#  simdpath = system.file("extdata", "simdata.csv", package = "metfal")
#  
#  ### read in the simulated data
#  simd = read.metfile(simdpath,nrowheader=4,ncolheader=7, metcol=1,smplrow=4)
#  
#  ### Generate missing rate matrix and extract 4 features
#  metf = makefeature(data=simd,wsize=10,ssize=0.5,defswitch=0.2)
#  
#  ### find false signal
#  falsig = falsignal(l=metf)
#  
#  ### plot qc figures
#  plot_mrate(falsig)
#  plot_qc(falsig)
#  
#  ### output false and good metabolites csv files
#  write.metfile(simd,falsig)

