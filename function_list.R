#############################################################################################
## function list
#############################################################################################

axis_limits <- function(species){
  if(species[1]=="CO2"){
    limits <- c(375,442)
  } else if(species[1]=="CH4"){
    limits <- c(1885,2120)
  }
  return(limits)
}

#############################################################################################

DateTime.DST <- function(year,mon,day,hour,min,sec) {
  
  if (mon<=9) {
    mon = as.character(paste("0",as.character(mon),sep=""))
  }              
  if (day<=9) {
    day = as.character(paste("0",as.character(day),sep=""))
  }
  if (hour<=9) {
    hour = as.character(paste("0",as.character(hour),sep=""))
  }
  if (min<=9) {
    min = as.character(paste("0",as.character(min),sep=""))
  }
  if (sec<=9) {
    sec = as.character(paste("0",as.character(sec),sep=""))
  }    
  
  dt <- as.character(paste(year,"-",mon,"-",day," ",hour,":",min,":",sec,sep=""))
  print(dt)
  dt <- as.POSIXct(dt,tz="ETC/GMT+1")
  
  return(dt)
}

#############################################################################################       

Daytime_mean <- function(timeseries) {
  library(lubridate)
  
  subset_daytime <- with(timeseries, timeseries[hour(sampling_datetime) >= 11 & hour(sampling_datetime) < 18,])
  subset_daytime$date <- as.character(subset_daytime$sampling_datetime)
  subset_daytime$date <- substr(subset_daytime$date,1,nchar(subset_daytime$date)-9)
  subset_daytime$ymd<- ymd(subset_daytime$date)
  timeseries <- aggregate(concentration ~ ymd, subset_daytime, mean)
  colnames(timeseries) <- c("sampling_datetime","concentration")
  timeseries$sampling_datetime <- as.POSIXct(timeseries$sampling_datetime)
  return(timeseries)
}

#############################################################################################

diff_Loess_REBS <- function(smoothed, REBS_df, timeseries){
  diff_LOESS_REBS <- data.frame(smoothed[1], REBS_df$fit, timeseries$sampling_datetime)
  colnames(diff_LOESS_REBS) <- c("smoothed", "REBS", "sampling_datetime")
  diff_LOESS_REBS$diff <- NA
  for (i in 1:nrow(diff_LOESS_REBS)) {
    diff_LOESS_REBS[i,4] = diff_LOESS_REBS[i,1]-diff_LOESS_REBS[i,2]
  }
  return(diff_LOESS_REBS)
}

#############################################################################################

ExtractDJF <- function(ccgvu_curve) {
  # species=get_species(ccgvu_curve)
  ccgvu_curve$date <- ccgvu_curve$sampling_datetime
  ccgvu_curve$season <- mkseas(ccgvu_curve, width="DJF")
  # DJF_curve <- paste("ccgvu_curve_",species[1],"_DJF", sep="")
  return(subset(ccgvu_curve, ccgvu_curve$season=="DJF")) # & !grepl("02-29", ccgvu_curve$date)))
}

#############################################################################################

Extract_winter <- function(data) {

  for (i in 1:nrow(data)) {
    if (month(data$sampling_datetime[i]) == 4 | month(data$sampling_datetime[i]) == 5 | 
        month(data$sampling_datetime[i]) == 6 | month(data$sampling_datetime[i]) == 7 | 
        month(data$sampling_datetime[i]) == 8 | month(data$sampling_datetime[i]) == 9 |
        month(data$sampling_datetime[i]) == 10) {
      data$summer[i] <- 0 
      data$no_anom[i] <- NA
      data$pos_SA[i] <- NA
      data$neg_SA[i] <- NA
    }
  }
  
  return(data)
  
}

#############################################################################################

Extract_summer <- function(data) {
  
  for (i in 1:nrow(data)) {
    if (month(data$sampling_datetime[i]) == 11 | month(data$sampling_datetime[i]) == 12 | 
        month(data$sampling_datetime[i]) == 1 | month(data$sampling_datetime[i]) == 2 | 
        month(data$sampling_datetime[i]) == 3) {
      data$winter[i] <- 0 
      data$no_anom[i] <- NA
      data$pos_SA[i] <- NA
      data$neg_SA[i] <- NA
    }
  }
  
  return(data)
  
}

#############################################################################################

find.freq <- function(x) {
  
  n <- length(x)
  spec <- spec.ar(c(x),plot=FALSE)
  print(max(spec$spec))
  if(max(spec$spec)>10) # Arbitrary threshold chosen by trial and error.
  {
    period <- round(1/spec$freq[which.max(spec$spec)])
    if(period==Inf) # Find next local maximum
    {
      j <- which(diff(spec$spec)>0)
      if(length(j)>0)
      {
        nextmax <- j[1] + which.max(spec$spec[j[1]:500])
        period <- round(1/spec$freq[nextmax])
      }
      else
        period <- 1
    }
  }
  else
    period <- 1
  return(period)
}

#############################################################################################

find_leap = function(x) {
  day(x) == 29 & month(x) == 2 
}

#############################################################################################

g_legend <- function(ggplot) {
    
  tmp <- ggplot_gtable(ggplot_build(ggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
    
}

#############################################################################################

get_species <- function(ccgvu_curve){
  if(mean(ccgvu_curve$orig, na.rm=TRUE)>=1400){
    species="CH4"
    units="ppb"
    print("species=CH4")
  } else {
    species="CO2"
    units="ppm"
    print("species=CO2")
  }
  spec_units <- list(species, units)
  return(spec_units)
}

#############################################################################################  
  
loess.sm <- function(data, x, y, span) {
  
  ss <- data
  dd <- data
  dd[[x]] <- decimal_date(dd[[x]])
  
  loessMod <- loess(dd[[y]] ~ dd[[x]], data=dd, span=span)
  s <- predict(loessMod)
  smoothed <- data.frame(ss[[x]], s)
  
  colnames(smoothed) = c(x,"smoothed")
  return(smoothed)
}

#############################################################################################

LoessSmooth <- function(timeseries, span){
  
  library(FedData)
  timeseries$date <- decimal_date(timeseries$sampling_datetime)
  loessMod <- loess(concentration ~ date, data=timeseries, span=span)
  smoothed <- predict(loessMod)
  
  if(timeseries$species[1]=="co2"){
    print("species=co2")
    plot(timeseries$concentration, x=timeseries$date, type="l", main="Loess Smoothing and Prediction",
         xlab="Date", ylab=expression(paste(CO[2]," (ppm)")))
    lines(smoothed, x=timeseries$date, col="red")
    smoothed <- data.frame(smoothed, timeseries$sampling_datetime)
    colnames(smoothed) <- c(paste("smoothed",(substr_right(span, 2)),sep=""), "sampling_datetime")
  } else if(timeseries$species[1]=="ch4"){
    print("species=ch4")
    plot(timeseries$concentration, x=timeseries$date, type="l", main="Loess Smoothing and Prediction",
         xlab="Date", ylab=expression(paste(CH[4]," (ppb)")))
    lines(smoothed, x=timeseries$date, col="red")
    smoothed <- data.frame(smoothed, timeseries$sampling_datetime)
    colnames(smoothed) <- c(paste("smoothed",(substr_right(span, 2)),sep=""), "sampling_datetime")
  }
  return(smoothed)
}

#############################################################################################

monthly_spei <- function(df, sp) {
  df$summer_SPEI <- NA
  for (i in 1:nrow(df)) {
    if (df$sampling_datetime[i] %in% sp$date) {
      df$summer_SPEI[i] <- sp$spei[which(sp$date == df$sampling_datetime[i])]
    }
  }
  df$summer_SPEI <- na_interpolation(df$summer_SPEI)
  for (i in 1:nrow(df)) {
    if (((month(df$sampling_datetime[i]) < 4) | (month(df$sampling_datetime[i]) > 10)) | 
        (year(df$sampling_datetime[i]) > 2018)) {
      df$summer_SPEI[i] <- NA
    }
  }
  return(df)
}

#############################################################################################

mmnorm <-
structure(function (data, minval = 0, maxval = 1) 
{
    d = dim(data)
    c = class(data)
    cnames = colnames(data)
    classes = data[, d[2]]
    data = data[, -d[2]]
    minvect = apply(data, 2, min)
    maxvect = apply(data, 2, max)
    rangevect = maxvect - minvect
    zdata = scale(data, center = minvect, scale = rangevect)
    newminvect = rep(minval, d[2] - 1)
    newmaxvect = rep(maxval, d[2] - 1)
    newrangevect = newmaxvect - newminvect
    zdata2 = scale(zdata, center = FALSE, scale = (1/newrangevect))
    zdata3 = zdata2 + newminvect
    zdata3 = cbind(zdata3, classes)
    if (c == "data.frame") 
        zdata3 = as.data.frame(zdata3)
    colnames(zdata3) = cnames
    return(zdata3)
}, source = c("function (data,minval=0,maxval=1) ", "{", "#This is a function to apply min-max normalization to a matrix or dataframe.", 
"#Min-max normalization subtracts the minimum of an attribute from each value", 
"#of the attribute and then divides the difference by the range of the attribute.", 
"#These new values are multiplied by the given range of the attribute", 
"#and finally added to the given minimum value of the attribute.", 
"#These operations transform the data into [minval,mxval].", 
"#Usually minval=0 and maxval=1.", "#Uses the scale function found in the R base package.", 
"#Input: data= The matrix or dataframe to be scaled", "", "", 
"#store all attributes of the original data", "d=dim(data)", 
"c=class(data)", "cnames=colnames(data)", "", "#remove classes from dataset", 
"classes=data[,d[2]]", "data=data[,-d[2]]", "", "minvect=apply(data,2,min)", 
"maxvect=apply(data,2,max)", "rangevect=maxvect-minvect", "zdata=scale(data,center=minvect,scale=rangevect)", 
"", "#remove attributes added by the function scale and turn resulting", 
"#vector back into a matrix with original dimensions", "#attributes(zdata)=NULL", 
"#zdata=matrix(zdata,dim(data)[1],dim(data)[2])", "", "newminvect=rep(minval,d[2]-1)", 
"newmaxvect=rep(maxval,d[2]-1)", "newrangevect=newmaxvect-newminvect", 
"zdata2=scale(zdata,center=FALSE,scale=(1/newrangevect))", "zdata3=zdata2+newminvect", 
"", "zdata3=cbind(zdata3,classes)", "", "if (c==\"data.frame\") zdata3=as.data.frame(zdata3)", 
"colnames(zdata3)=cnames", "return(zdata3)", "", "}"))

#############################################################################################

NA.loess.smooth <- function(data, x, y, span) {
  
  dd <- subset(data, !is.na(data[[y]]))
  ss <- subset(data, data[[x]] %in% dd[[x]])
  dd[[x]] <- decimal_date(dd[[x]])

  loessMod <- loess(dd[[y]] ~ dd[[x]], data=dd, span=span) 
  s <- predict(loessMod)
  
  smoothed <- data.frame(ss[[x]], s)
  colnames(smoothed) <- c(x,y)
  
  interp <- data.frame(data[[x]])
  interp[[y]] <- NA
  colnames(interp) <- c(x,y)
  
  k=1
  for (i in 1:nrow(interp)) {
    if (interp[i,x] %in% smoothed[[x]]) {
      for (j in k:nrow(smoothed)) {
        if (smoothed[j,x] == interp[i,x]) {
          interp[i,y] <- smoothed[j,y]
          k = j
          break
        }
      }
    }
  }
  i <- na_interpolation(interp[[y]])
  interp$i <- i
  
  dd <- interp
  dd[[x]] <- decimal_date(dd[[x]])
  
  loessMod <- loess(dd$i ~ dd[[x]], data=dd, span=span) 
  s <- predict(loessMod)
  
  smoothed <- data.frame(interp[[x]], s)
  colnames(smoothed) = c(x,"smoothed")
  
  return(smoothed)
  
}   

#############################################################################################

NaN.to.na <- function(data) {
  
  for (i in 1:nrow(data)) {
    if (is.nan(data$concentration[i])) {
      data$concentration[i] <- NA
    }
  }
  return(data)
}

#############################################################################################

negative.to.na <- function(data) {
  
  for (i in 1:nrow(data)) {
    if (!is.na(data$std_dev[i])) {
      if (data$std_dev[i] < 0) {
        data$std_dev[i] <- NA
      }
    }
  }
  return(data)
}

#############################################################################################

sd_nr <- function(data) {
  
  sd_nr <- setNames(data.frame(matrix(ncol = 1, nrow = nrow(data))), 
                         c("sd_nr"))
  
  for (i in 1:nrow(data)) {
    if (!is.na(data$median_conc[i])) {
      conc <- setNames(data.frame(matrix(ncol = 1, nrow = 504)), c("concentration"))
      conc$concentration <- data$concentration[(i-251):(i+252)]
      
      for (j in 1:nrow(conc)) {
        if (!is.na(conc$concentration[j])) {
          if (conc$concentration[j] > data$median_conc[i]) {
            conc$concentration[j] <- NA
          }
        }
      }

      normalized <- setNames(data.frame(matrix(ncol = 2, nrow = sum(!is.na(conc$concentration)))), 
                             c("negative", "positive"))
      
      normalized$negative <- subset(conc$concentration, !is.na(conc$concentration)) 
      for (j in 1:nrow(normalized)) {
        normalized$positive[j] <- normalized$negative[j] + (2*(data$median_conc[i]-normalized$negative[j]))
      }
      
      df_norm <- data.frame(res = c(normalized[,"negative"], normalized[,"positive"]))
      sd_nr$sd_nr[i] <- sd(df_norm$res)
      print(i)
    }
  }
  return(sd_nr)
}

#############################################################################################

taylor.diagram.x <- function (ref, model, add = FALSE, col = "red", pch = 19, pos.cor = TRUE, 
                    xlab = "", ylab = "", main = "Taylor Diagram", show.gamma = TRUE, 
                    ngamma = 3, gamma.col = 8, sd.arcs = 0, ref.sd = FALSE, 
                    sd.method = "sample", grad.corr.lines = c(0.2, 0.4, 0.6, 0.8, 0.9), 
                    pcex = 1, cex.axis = 1, normalize = FALSE, mar = c(5, 4, 6, 6),
                    text, ...) #the added parameter
{
  grad.corr.full <- c(0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99, 
                      1)
  R <- cor(ref, model, use = "pairwise")
  if (is.list(ref)) 
    ref <- unlist(ref)
  if (is.list(model)) 
    ref <- unlist(model)
  SD <- function(x, subn) {
    meanx <- mean(x, na.rm = TRUE)
    devx <- x - meanx
    ssd <- sqrt(sum(devx * devx, na.rm = TRUE)/(length(x[!is.na(x)]) - 
                                                  subn))
    return(ssd)
  }
  subn <- sd.method != "sample"
  sd.r <- SD(ref, subn)
  sd.f <- SD(model, subn)
  if (normalize) {
    sd.f <- sd.f/sd.r
    sd.r <- 1
  }
  maxsd <- 1.5 * max(sd.f, sd.r)
  oldpar <- par("mar", "xpd", "xaxs", "yaxs")
  if (!add) {
    par(mar = mar)
    if (pos.cor) {
      if (nchar(ylab) == 0) 
        ylab = "Standard deviation"
      plot(0, xlim = c(0, maxsd * 1.1), ylim = c(0, maxsd * 
                                                   1.1), xaxs = "i", yaxs = "i", axes = FALSE, main = main, 
           xlab = "", ylab = ylab, type = "n", cex = cex.axis, 
           ...)
      mtext(xlab, side = 1, line = 2.3)
      if (grad.corr.lines[1]) {
        for (gcl in grad.corr.lines) lines(c(0, maxsd * 
                                               gcl), c(0, maxsd * sqrt(1 - gcl^2)), lty = 3)
      }
      segments(c(0, 0), c(0, 0), c(0, maxsd), c(maxsd, 
                                                0))
      axis.ticks <- pretty(c(0, maxsd))
      axis.ticks <- axis.ticks[axis.ticks <= maxsd]
      axis(1, at = axis.ticks, cex.axis = cex.axis)
      axis(2, at = axis.ticks, cex.axis = cex.axis)
      if (sd.arcs[1]) {
        if (length(sd.arcs) == 1) 
          sd.arcs <- axis.ticks
        for (sdarc in sd.arcs) {
          xcurve <- cos(seq(0, pi/2, by = 0.03)) * sdarc
          ycurve <- sin(seq(0, pi/2, by = 0.03)) * sdarc
          lines(xcurve, ycurve, col = "blue", lty = 3)
        }
      }
      if (show.gamma[1]) {
        if (length(show.gamma) > 1) 
          gamma <- show.gamma
        else gamma <- pretty(c(0, maxsd), n = ngamma)[-1]
        if (gamma[length(gamma)] > maxsd) 
          gamma <- gamma[-length(gamma)]
        labelpos <- seq(45, 70, length.out = length(gamma))
        for (gindex in 1:length(gamma)) {
          xcurve <- cos(seq(0, pi, by = 0.03)) * gamma[gindex] + 
            sd.r
          endcurve <- which(xcurve < 0)
          endcurve <- ifelse(length(endcurve), min(endcurve) - 
                               1, 105)
          ycurve <- sin(seq(0, pi, by = 0.03)) * gamma[gindex]
          maxcurve <- xcurve * xcurve + ycurve * ycurve
          startcurve <- which(maxcurve > maxsd * maxsd)
          startcurve <- ifelse(length(startcurve), max(startcurve) + 
                                 1, 0)
          lines(xcurve[startcurve:endcurve], ycurve[startcurve:endcurve], 
                col = gamma.col)
          if (xcurve[labelpos[gindex]] > 0) 
            boxed.labels(xcurve[labelpos[gindex]], ycurve[labelpos[gindex]], 
                         gamma[gindex], border = FALSE)
        }
      }
      xcurve <- cos(seq(0, pi/2, by = 0.01)) * maxsd
      ycurve <- sin(seq(0, pi/2, by = 0.01)) * maxsd
      lines(xcurve, ycurve)
      bigtickangles <- acos(seq(0.1, 0.9, by = 0.1))
      medtickangles <- acos(seq(0.05, 0.95, by = 0.1))
      smltickangles <- acos(seq(0.91, 0.99, by = 0.01))
      segments(cos(bigtickangles) * maxsd, sin(bigtickangles) * 
                 maxsd, cos(bigtickangles) * 0.97 * maxsd, sin(bigtickangles) * 
                 0.97 * maxsd)
      par(xpd = TRUE)
      if (ref.sd) {
        xcurve <- cos(seq(0, pi/2, by = 0.01)) * sd.r
        ycurve <- sin(seq(0, pi/2, by = 0.01)) * sd.r
        lines(xcurve, ycurve)
      }
      points(sd.r, 0, cex = pcex)
      text(cos(c(bigtickangles, acos(c(0.95, 0.99)))) * 
             1.05 * maxsd, sin(c(bigtickangles, acos(c(0.95, 
                                                       0.99)))) * 1.05 * maxsd, c(seq(0.1, 0.9, by = 0.1), 
                                                                                  0.95, 0.99), cex = cex.axis)
      text(maxsd * 0.8, maxsd * 0.8, "Correlation", srt = 315, 
           cex = cex.axis)
      segments(cos(medtickangles) * maxsd, sin(medtickangles) * 
                 maxsd, cos(medtickangles) * 0.98 * maxsd, sin(medtickangles) * 
                 0.98 * maxsd)
      segments(cos(smltickangles) * maxsd, sin(smltickangles) * 
                 maxsd, cos(smltickangles) * 0.99 * maxsd, sin(smltickangles) * 
                 0.99 * maxsd)
    }
    else {
      x <- ref
      y <- model
      R <- cor(x, y, use = "pairwise.complete.obs")
      E <- mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE)
      xprime <- x - mean(x, na.rm = TRUE)
      yprime <- y - mean(y, na.rm = TRUE)
      sumofsquares <- (xprime - yprime)^2
      Eprime <- sqrt(sum(sumofsquares)/length(complete.cases(x)))
      E2 <- E^2 + Eprime^2
      if (add == FALSE) {
        maxray <- 1.5 * max(sd.f, sd.r)
        plot(c(-maxray, maxray), c(0, maxray), type = "n", 
             asp = 1, bty = "n", xaxt = "n", yaxt = "n", 
             xlim = c(-1.1 * maxray, 1.1 * maxray), xlab = xlab, 
             ylab = ylab, main = main, cex = cex.axis)
        discrete <- seq(180, 0, by = -1)
        listepoints <- NULL
        for (i in discrete) {
          listepoints <- cbind(listepoints, maxray * 
                                 cos(i * pi/180), maxray * sin(i * pi/180))
        }
        listepoints <- matrix(listepoints, 2, length(listepoints)/2)
        listepoints <- t(listepoints)
        lines(listepoints[, 1], listepoints[, 2])
        lines(c(-maxray, maxray), c(0, 0))
        lines(c(0, 0), c(0, maxray))
        for (i in grad.corr.lines) {
          lines(c(0, maxray * i), c(0, maxray * sqrt(1 - 
                                                       i^2)), lty = 3)
          lines(c(0, -maxray * i), c(0, maxray * sqrt(1 - 
                                                        i^2)), lty = 3)
        }
        for (i in grad.corr.full) {
          text(1.05 * maxray * i, 1.05 * maxray * sqrt(1 - 
                                                         i^2), i, cex = cex.axis, adj = cos(i)/2)
          text(-1.05 * maxray * i, 1.05 * maxray * sqrt(1 - 
                                                          i^2), -i, cex = cex.axis, adj = 1 - cos(i)/2)
        }
        seq.sd <- seq.int(0, 2 * maxray, by = (maxray/10))[-1]
        for (i in seq.sd) {
          xcircle <- sd.r + (cos(discrete * pi/180) * 
                               i)
          ycircle <- sin(discrete * pi/180) * i
          for (j in 1:length(xcircle)) {
            if ((xcircle[j]^2 + ycircle[j]^2) < (maxray^2)) {
              points(xcircle[j], ycircle[j], col = "darkgreen", 
                     pch = ".")
              if (j == 10) 
                text(xcircle[j], ycircle[j], signif(i, 
                                                    2), cex = cex.axis, col = "darkgreen", 
                     srt = 90)
            }
          }
        }
        seq.sd <- seq.int(0, maxray, length.out = 5)
        for (i in seq.sd) {
          xcircle <- cos(discrete * pi/180) * i
          ycircle <- sin(discrete * pi/180) * i
          if (i) 
            lines(xcircle, ycircle, lty = 3, col = "blue")
          text(min(xcircle), -0.06 * maxray, signif(i, 
                                                    2), cex = cex.axis, col = "blue")
          text(max(xcircle), -0.06 * maxray, signif(i, 
                                                    2), cex = cex.axis, col = "blue")
        }
        text(0, -0.14 * maxray, "Standard Deviation", 
             cex = cex.axis, col = "blue")
        text(0, -0.22 * maxray, "Centered RMS Difference", 
             cex = cex.axis, col = "darkgreen")
        points(sd.r, 0, pch = 22, bg = "darkgreen", cex = pcex)
        text(0, 1.2 * maxray, "Correlation Coefficient", 
             cex = cex.axis)
      }
      S <- (2 * (1 + R))/(sd.f + (1/sd.f))^2
    }
  }
  points(sd.f * R, sd.f * sin(acos(R)), pch = pch, col = col, 
         cex = pcex)
  text(sd.f * R, sd.f * sin(acos(R)),  #Added line
       labels=text, cex = pcex, pos=3) #pos argument can be adjusted as needed
  invisible(oldpar)
}

#############################################################################################

zero.to.na <- function(data) {
  
  for (i in 1:nrow(data)) {
    if (!is.na(data$concentration[i])) {
      if (data$concentration[i] < 1) {
        data$concentration[i] <- NA
      }
    }
  }
  return(data)
}
