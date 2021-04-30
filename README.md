# RuminatoR - an R package to classify feeding behavior of ruminants

`RuminatoR` allows to classify different behaviors of cows and other
ruminants based on pressure signals measured in the headcollar.

This is a reimplementation based on Scheidegger (2008).


## Installation

1. Install [R](https://cloud.r-project.org/) and [R-Studio](https://www.rstudio.com/products/RStudio/) (or any other editor).

2. Install `remotes` (type in the R command line)
```
install.packages("remotes")
```

3. Install `RuminatoR` (type in the R command line)
```
library(remotes)
install_github("scheidan/RuminatoR")
```

## Usage

You can get help in R with `?train` or `?classify`.

### Example

We assume the data are available as dataframes in the following structure:
```R
head(trainings.data)
                 time pressure   activity
1 2008-09-04 20:34:50   1076.8 ruminating
2 2008-09-04 20:34:50   1043.6 ruminating
3 2008-09-04 20:34:50   1029.1 ruminating

head(new.data)
                 time pressure
1 2008-09-04 06:59:59   1085.4
2 2008-09-04 06:59:59   1069.1
3 2008-09-04 06:59:59   1049.8
```

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
table(dat.class$activity.predicted)

## -----------
## 3) visualize result

## color line
ggplot(dat.class, aes(x=time, y=pressure, color=activity.predicted)) +
    geom_path(aes(group=1))

## color background
ggplot(dat.class) +
    ## background
    geom_rect(aes(NULL, NULL,
                  xmin=time-0.5, xmax=time+0.5,
                  ymin=-Inf, ymax=Inf,
                  fill=activity.predicted)) +
    ## add pressure
    geom_path(aes(x=time, y=pressure))
```

## References

Scheidegger, A. (2008). Klassifikation des Fressverhaltens von Kühen. Degree thesis, Zurich University of Applied Sciences, Zurich.
