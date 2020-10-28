# RuminatoR - an R package to classify feeding behavior of ruminants

`RuminatoR` allows to classify different behaviors of cows and other
ruminants based on pressure signals measured in the headcollar.

This is a reimplementation based on Scheidegger (2008).


## Installation

1. Install [R](https://cloud.r-project.org/) and [R-Studio](https://www.rstudio.com/products/RStudio/) (or any other editor).

2. Install `devtools` (type in the R command line)
```
install.packages("devtools")
```

3. Install RuminatoR (type in the R command line)
```
library(devtools)
install_github("scheidan/RuminatoR")
```

## Usage

... TODO ...

```R
library(RuminatoR)
library(ggplot2)  # optional, only for visualization

## --------------------------
## 1) train classifier

rf.classifier <- train(data.training,
                       min.amplitude=30, min.dt=6)

## ---------------------------------
## 2) apply classifier to new data

dat.class <- classify(rf.classifier, new.data,
                      max.dt.peak = 30, min.n.peaks = 10)

head(dat.class)
table(dat.class$activity)

## -----------
## 3) visualize result

## color line
ggplot(dat.class, aes(x=time, y=pressure, color=activity)) +
    geom_path(aes(group=1))

## color background
ggplot(dat.class) +
    ## background
    geom_rect(aes(NULL, NULL,
                  xmin=time-0.5, xmax=time+0.5,
                  ymin=-Inf, ymax=Inf,
                  fill=activity)) +
    ## add pressure
    geom_path(aes(x=time, y=pressure))
```

## References

Scheidegger, A. (2008). Klassifikation des Fressverhaltens von Kühen. Degree thesis, Zurich University of Applied Sciences, Zurich.
