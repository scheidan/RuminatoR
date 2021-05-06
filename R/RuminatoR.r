## -------------------------------------------------------
##
## RuminatoR
##
## Andreas Scheidegger -- andreas.scheidegger@eawag.ch
## -------------------------------------------------------


##' @import randomForest
##' @importFrom stats mad median na.omit predict runmed sd
NULL

##' RuminatoR - analyze the feeding behavior of ruminates
##'
##' .. TODO ..
##' 
##' @name Ruminator
##' @author Andreas Scheidegger
##' @docType package
NULL


##' @title Identify peaks in time series
##'
##' @param pressure time series of pressure data
##' @param min.amplitude  minimal amplitude to count as peak
##' @param min.dt minimal time difference between peaks
##' @return a vector of the indicies of the peaks
identify.peaks <- function(pressure, min.amplitude=30, min.dt=6){
    ## smooth data wth running median
    pressure.smooth <- runmed(x=pressure, k=101)

    ## find peaks with >minimal.amplitude
    n <- length(pressure.smooth)
    tmp <- data.frame('tn1'=c(NA,pressure[1:(n-1)]),
                      't'  =pressure,
                      'tp1'=c(pressure[2:n],NA))

    peak <- which((tmp$t >= tmp$tn1) & (tmp$t > tmp$tp1) & # local maximum
                  (tmp$t - pressure.smooth > min.amplitude) ) # check amplitude

    ## make sure peaks have enough distance
    i <- 1
    while (length(peak)>i) {
	if (peak[i+1]-peak[i] < min.dt){ # if time difference too small, take larger peak          
            if (pressure[peak[i]] > pressure[peak[i+1]]) {
                peak <- peak[-(i+1)] 
            } else { 
                peak <- peak[-i] }
            i <- i-1
        }        
	i <- i+1
    }

    peak
}

##' @title Compute the statistics for classification
##' 
##' @param data data.frame wite a columns pressure
##' @param prediction if FALSE, the obeserved activities are added to the output 
##' @param min.amplitude  minimal amplitude to count as peak
##' @param min.dt minimal time difference between peaks
##' @return a list with the statistics of the peaks, and the time points
compute.statistics <- function(data, prediction=FALSE,
                               min.amplitude=30, min.dt=6){
    ## identify indices of peaks
    peaks <- identify.peaks(data$pressure, min.amplitude=min.amplitude, min.dt=min.dt)
    dt.peak <- diff(peaks)            # time difference between peaks

    ## compute statistics for classification
    stats <- matrix(NA, nrow=length(peaks), ncol=13)

    colnames(stats) <- c('Mean.dt.peak','Med.dt.peak','Mad.dt.peak','sd.dt.peak',
                         'Mad.dt.peak.z','sd.dt.peak.z','Mad.dt.peak.v','sd.dt.peak.v',
                         'maxdt.peak','Mad.Peak','sd.Peak','Mad.pressure','sd.pressure')

    for (i in 31:(length(dt.peak)-31)) {

        ## -- statistics about the last / next 30 peaks --
        ## Mean and median of time difference of peaks
        stats[i,1] <- mean(dt.peak[(i-20):(i+20)])
        stats[i,2] <- median(dt.peak[(i-20):(i+20)])

        ## MAD and sd of time difference of peaks
        stats[i,3] <- mad(dt.peak[(i-20):(i+20)])
        stats[i,4] <- sd(dt.peak[(i-20):(i+20)])

        ## MAD and sd of future time difference, Mad.dt.peak.z
        stats[i,5] <- mad(dt.peak[i:(i+20)])
        stats[i,6] <- sd(dt.peak[i:(i+20)])

        ## MAD and sd of past time difference, Mad.dt.peak.v
        stats[i,7] <- mad(dt.peak[(i-20):i])
        stats[i,8] <- sd(dt.peak[(i-20):i])

        ## maximum time difference of peaks
        stats[i,9] <- max(dt.peak[(i-30):(i+30)])

        ## MAD and sd of peak pressure
        stats[i,10] <- mad(data$pressure[peaks[(i-20):(i+20)]])
        stats[i,11] <- sd(data$pressure[peaks[(i-20):(i+20)]])

        ## -- statistics about the last / next 30 peaks --
        ## MAD and sd of pressures +/- 10 sek, Mad.pressure
        stats[i,12] <- mad(data$pressure[(peaks[i]-100):(peaks[i]+100)])
        stats[i,13] <- sd(data$pressure[(peaks[i]-100):(peaks[i]+100)])
    }

    ## add observed activities and remove NAs
    nNa <- !is.na(stats[,1])
    ret.list <- list(X = stats[nNa,],
                     time = data$time[peaks][nNa] )
    if(!prediction) { ret.list[['activity']]  <-  data$activity[peaks][nNa] }

    return(ret.list)
}


##' @title Train Random Forest
##' 
##' @param data A \code{data.frame} wite a columns \code{time}, \code{pressure}, and \code{activity}
##' @param min.amplitude  minimal amplitude to count as peak
##' @param min.dt minimal time difference between peaks
##' @param ... arguments passed to \code{randomForest}
##' @return list with the trained randomForest and the parameters used.
##' @author Andreas Scheidegger
##' @export
train <- function(data, min.amplitude=30, min.dt=6, ...){

    if(!all(c("time", "pressure", "activity") %in% colnames(data))){
        stop("Training data must have columns 'time', 'pressure', and 'activity'!")
    }

    data$activity <- as.factor(data$activity)
    
    ## compute statistics
    dd <- compute.statistics(data, prediction=FALSE,
                             min.amplitude=min.amplitude , min.dt=min.dt)

    ## train
    list(classifier = randomForest::randomForest(x=dd$X, y=dd$activity, ...),
         parameter = c(min.amplitude=min.amplitude, min.dt=min.dt))

}


## return most frequent value of a vector
val.most <- function(a) {
    ta <- table(a)
    b <- which(ta == max(ta))

    ## if mulitble with same frequencey, choose one randomly
    if (length (b) > 1) {
        b <- sample(b, 1)
    }
    names(ta[b])
}



## Group peaks to blocks with same activity
##
## max.dt.peaks:  maximal time difference between two peaks in the same group
## min.n.peaks:  minimal number of peaks to count as group
group.peaks <- function(peaks,
                        max.dt.peaks = 30,
                        min.n.peaks = 10){
    
    groups <- rep(NA, length(peaks))

    groups[1] <- 1
    j <- 1
    for (i in 2:length(peaks)) {
	if(peaks[i]-peaks[i-1] < max.dt.peaks) {
            groups[i] <- groups[i-1]
        } else {     # if time difference larger, make new group
            j <- j+1 
            groups[i] <- j
        } 
 	
    }

    ## delete groups that are too shortgroups
    groups[groups %in% which(table(groups)<min.n.peaks)] <- NA

    groups
}


##' @title Classify pressure time series
##'
##' @param RF a trained random forest
##' @param newdata data.frame wite a columns \code{time}, and \code{pressure}
##' @param max.dt.peaks  maximal time difference between two peaks in the same group
##' @param min.n.peaks  minimal number of peaks to count as activity group
##' @return A \code{data.frame} with additional column of the classfied activities
##' @author Andreas Scheidegger
##' @export
classify <- function(RF, newdata, max.dt.peaks = 30, min.n.peaks = 10){

    if(!all(c("time", "pressure") %in% colnames(newdata))){
        stop("Data must have columns 'time' and 'pressure'!")
    }

    
    ## -- get peaks
    peaks <- identify.peaks(newdata$pressure)

    ## -- group peaks
    groups <- group.peaks(peaks, max.dt.peaks = max.dt.peaks,
                          min.n.peaks = min.n.peaks)
    
    ## -- compute statistics
    dd <- compute.statistics(newdata, prediction=TRUE,
                             min.amplitude = RF$parameter["min.amplitude"],
                             min.dt = RF$parameter["min.dt"])

    ## -- classify
    activity <- predict(RF$classifier, newdata=dd$X, type="response")
    
    ## -- assign predictions to groups

    ## peaks that belong to no group
    activity[is.na(groups)] <- NA

    ## choose the most frequent prediction for whole group
    for (i in na.omit(unique(groups))) {
        activity[groups == i] <- val.most(activity[groups==i])
    }

    ## -- map predictions to time series
    newdata$activity.predicted <- NA
    
    for (g in 1:max(groups, na.rm=TRUE)) {
        idx <- which(groups == g)
        if(length(idx) > 0) {
            t.start <- peaks[min(idx)]
            t.end <- peaks[max(idx)]
            newdata$activity.predicted[t.start:t.end] <- as.character(activity[min(idx)])
        }
    }

    newdata$activity.predicted <- as.factor(newdata$activity.predicted)
   
    ## return data frame with activity
    newdata
}
