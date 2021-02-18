options(java.parameters = "- Xmx1024m")
source("C:/Users/aresovsk/Documents/R/main.R")

#############################################################################


## 1. Read in the 2011-2020 data and subset dates
  OPE_co2_2011_2019 <- readRDS('OPE_2011-2019_co2_120m.rds')
  OPE_co2_2011_2019$sampling_datetime <- as.POSIXct(OPE_co2_2011_2019$sampling_datetime, format="%Y-%m-%d %H:%M:%OS", tz="GMT")
  OPE_co2_2011_2019 <- OPE_co2_2011_2019[OPE_co2_2011_2019$sampling_datetime < "2019-05-01",]
  OPE_co2_2011_2019 <- OPE_co2_2011_2019[,c(11,13)]
  OPE_co2_newL2 <- readRDS('OPE_2011-2020_co2_120m.rds')
  OPE_co2_newL2$sampling_datetime <- as.POSIXct(OPE_co2_newL2$sampling_datetime, format="%Y-%m-%d %H:%M:%OS", tz="GMT")
  OPE_co2_newL2 <- OPE_co2_newL2[OPE_co2_newL2$sampling_datetime >= "2019-05-01",]
  OPE_co2_newL2 <- OPE_co2_newL2[,c(1,2)]
  OPE_co2_2011_2020 <- rbind(OPE_co2_2011_2019, OPE_co2_newL2)
  rm(OPE_co2_2011_2019, OPE_co2_newL2)
  
## Fix missing dates
  hourly <- data.frame(seq(OPE_co2_2011_2020$sampling_datetime[1], OPE_co2_2011_2020$sampling_datetime[nrow(OPE_co2_2011_2020)], by="hour"))
  colnames(hourly) <- c("sampling_datetime")
  hourly$concentration <- NA
  df <- hourly
  for (i in 1:nrow(df)) {  
    print(i)
    if (df$sampling_datetime[i] %in% OPE_co2_2011_2020$sampling_datetime) {
      df$concentration[i] <- OPE_co2_2011_2020$concentration[which(df$sampling_datetime[i] == OPE_co2_2011_2020$sampling_datetime)]
    }
  }
  OPE_co2_2011_2020 <- df


## 2. Read in NRT data
  OPE_co2_NRT <- read.xlsx('C:/Users/aresovsk/ICOS_data/ICOS_NRT_growing_20200928/ICOS_ATC_NRT_OPE_2020-06-01_2020-09-27_120.0_CO2.xlsx',
                           sheetIndex = "hourly", header=TRUE)
  OPE_co2_NRT$sampling_datetime <- as.POSIXct(OPE_co2_NRT$sampling_datetime, format="%Y-%m-%d %H:%M:%OS", tz="GMT")
  OPE_co2_NRT <- OPE_co2_NRT[,c(1,2)]
  
  pastyear.hourly <- data.frame(seq(as.POSIXct((OPE_co2_2011_2020$sampling_datetime[nrow(OPE_co2_2011_2020)] +3600), tz="GMT"), 
                                    as.POSIXct(OPE_co2_NRT$sampling_datetime[nrow(OPE_co2_NRT)], tz="GMT"), by="hour"))
  colnames(pastyear.hourly) <- c("sampling_datetime")
  pastyear.hourly$concentration <- NA
  df <- pastyear.hourly
  for (i in 1:nrow(df)) {
    print(i)
    if (df$sampling_datetime[i] %in% OPE_co2_NRT$sampling_datetime) {
      df$concentration[i] <- OPE_co2_NRT$concentration[which(df$sampling_datetime[i] == OPE_co2_NRT$sampling_datetime)]
    }
  }
  OPE_co2_NRT <- df

  ## Concatenate
  OPE_co2_all <- rbind(OPE_co2_2011_2020, OPE_co2_NRT)
  ## Use only afternoon values
  OPE_co2_all <- subset(OPE_co2_all, (hour(as.POSIXct(OPE_co2_all$sampling_datetime)) >= 12) & 
                                     (hour(as.POSIXct(OPE_co2_all$sampling_datetime)) <= 17)) # or
  ## Use only nighttime values
  #OPE_co2_all <- subset(OPE_co2_all, (hour(as.POSIXct(OPE_co2_all$sampling_datetime)) >= 20) | 
  #                                   (hour(as.POSIXct(OPE_co2_all$sampling_datetime)) <= 5))


## 3. Aggregate hourly to daily
  OPE_co2_all$Date  <- as.Date(OPE_co2_all$sampling_datetime, "%Y-%m-%d", tz="GMT")
  OPE_co2_daily <- aggregate(concentration ~ Date, OPE_co2_all, mean, na.rm=TRUE, na.action=na.pass)
  OPE_co2_daily <- NaN.to.na(OPE_co2_daily)
  colnames(OPE_co2_daily) <- c("sampling_datetime", "concentration")
  OPE_co2_daily$sampling_datetime <- as.POSIXct(OPE_co2_daily$sampling_datetime, format="%Y-%m-%d", tz="GMT")
  OPE_co2_daily$sampling_datetime <- round(OPE_co2_daily$sampling_datetime, units = "days")
  OPE_co2_daily$sampling_datetime <- as.POSIXct(OPE_co2_daily$sampling_datetime, format="%Y-%m-%d", tz="GMT")
  OPE_co2_daily$concentration <- na_interpolation(OPE_co2_daily$concentration)   
  ## Remove leap days
  OPE_co2_daily <- subset(OPE_co2_daily, OPE_co2_daily$sampling_datetime %in% 
                          OPE_co2_daily$sampling_datetime[!find_leap(OPE_co2_daily$sampling_datetime)])
  

## 4. Extract the 90-day CCGvu baseline curve
  colNames <- c("func","poly","smooth","trend","detrend","smcycle","harm","res","smres","trres","ressm","gr")
  ccgvu_90d_OPE_co2 <- getCcgcrv(OPE_co2_daily, "sampling_datetime", "concentration", 
                             ccgcrvParameters="-all -short\ 90", merge=FALSE, colNames=colNames)
  ccgvu_30d_OPE_co2 <- getCcgcrv(OPE_co2_daily, "sampling_datetime", "concentration", 
                             ccgcrvParameters="-all -short\ 30", merge=FALSE, colNames=colNames)
  
  ## Round dates
  ccgvu_90d_OPE_co2$sampling_datetime <- round_date(ccgvu_90d_OPE_co2$sampling_datetime, unit = "day")
  ccgvu_30d_OPE_co2$sampling_datetime <- round_date(ccgvu_30d_OPE_co2$sampling_datetime, unit = "day")


## 5. Calculate delta-C values
  
  ccgvu_90d_OPE_co2$dC <- NA
  for (i in 1:nrow(ccgvu_90d_OPE_co2)) {
    ccgvu_90d_OPE_co2$dC[i] <- ccgvu_90d_OPE_co2$smcycle[i] - ccgvu_90d_OPE_co2$harm[i]
  } 
  ccgvu_30d_OPE_co2$dC <- NA
  for (i in 1:nrow(ccgvu_30d_OPE_co2)) {
    ccgvu_30d_OPE_co2$dC[i] <- ccgvu_30d_OPE_co2$smcycle[i] - ccgvu_30d_OPE_co2$harm[i]
  } 


## 6. Calculate season-adjusted sigma values 
  ccgvu_90d_OPE_co2$sd <- seasonal_sigma(ccgvu_90d_OPE_co2, ndays=(nrow(ccgvu_90d_OPE_co2)-1))
  ccgvu_30d_OPE_co2$sd <- seasonal_sigma(ccgvu_30d_OPE_co2, ndays=(nrow(ccgvu_30d_OPE_co2)-1))
  
#  # Use 2*sd ? 
  for (i in 1:nrow(ccgvu_90d_OPE_co2)) {
    ccgvu_90d_OPE_co2$sd[i] <- ccgvu_90d_OPE_co2$sd[i] * 2
  }
  for (i in 1:nrow(ccgvu_30d_OPE_co2)) {
    ccgvu_30d_OPE_co2$sd[i] <- ccgvu_30d_OPE_co2$sd[i] * 2
  }
  
  ccgvu_90d_OPE_co2$plus.sigma <- NA
  ccgvu_90d_OPE_co2$minus.sigma <- NA
  ccgvu_30d_OPE_co2$plus.sigma <- NA
  ccgvu_30d_OPE_co2$minus.sigma <- NA
  for (i in 1:nrow(ccgvu_90d_OPE_co2)) {
    ccgvu_90d_OPE_co2$plus.sigma[i]  <- ccgvu_90d_OPE_co2$harm[i] + ccgvu_90d_OPE_co2$trend[i] + ccgvu_90d_OPE_co2$sd[i]
    ccgvu_90d_OPE_co2$minus.sigma[i] <- ccgvu_90d_OPE_co2$harm[i] + ccgvu_90d_OPE_co2$trend[i] - ccgvu_90d_OPE_co2$sd[i]
    ccgvu_30d_OPE_co2$plus.sigma[i] <- ccgvu_30d_OPE_co2$harm[i] + ccgvu_30d_OPE_co2$trend[i] + ccgvu_30d_OPE_co2$sd[i]
    ccgvu_30d_OPE_co2$minus.sigma[i] <- ccgvu_30d_OPE_co2$harm[i] + ccgvu_30d_OPE_co2$trend[i] - ccgvu_30d_OPE_co2$sd[i]
  }

  ## Make smoothed curves
  OPE_90d_smoothed_co2 <- loess.sm(ccgvu_90d_OPE_co2, "sampling_datetime", "orig", span=(90/nrow(ccgvu_90d_OPE_co2)))
  OPE_30d_smoothed_co2 <- loess.sm(ccgvu_30d_OPE_co2, "sampling_datetime", "orig", span=(30/nrow(ccgvu_30d_OPE_co2)))
  
  
## 7. Identify seasonal anomalies 
  OPE_90d_smoothed_co2$no_anom <- NA
  OPE_90d_smoothed_co2$pos_SA <- NA
  OPE_90d_smoothed_co2$neg_SA <- NA
  OPE_90d_smoothed_co2$SSE <- NA
  for (i in 1:nrow(OPE_90d_smoothed_co2)) {
    if ((OPE_90d_smoothed_co2$smoothed[i] >= ccgvu_90d_OPE_co2$plus.sigma[i])) {
      OPE_90d_smoothed_co2$pos_SA[i] <- OPE_90d_smoothed_co2$smoothed[i] - ccgvu_90d_OPE_co2$plus.sigma[i]
      OPE_90d_smoothed_co2$SSE[i] <- OPE_90d_smoothed_co2$smoothed[i]
    } else if ((OPE_90d_smoothed_co2$smoothed[i] <= ccgvu_90d_OPE_co2$minus.sigma[i])) {
      OPE_90d_smoothed_co2$neg_SA[i] <- OPE_90d_smoothed_co2$smoothed[i] - ccgvu_90d_OPE_co2$minus.sigma[i]
      OPE_90d_smoothed_co2$SSE[i] <- OPE_90d_smoothed_co2$smoothed[i]
    } else {
      OPE_90d_smoothed_co2$no_anom[i] <- 0
    }
  }
  OPE_30d_smoothed_co2$no_anom <- NA
  OPE_30d_smoothed_co2$pos_SA <- NA
  OPE_30d_smoothed_co2$neg_SA <- NA
  OPE_30d_smoothed_co2$SSE <- NA
  for (i in 2:(nrow(OPE_30d_smoothed_co2)-1)) {
    if (((OPE_30d_smoothed_co2$smoothed[i] >= ccgvu_30d_OPE_co2$plus.sigma[i]) & (OPE_30d_smoothed_co2$smoothed[i-1] >= ccgvu_30d_OPE_co2$plus.sigma[i-1])) |
       ((OPE_30d_smoothed_co2$smoothed[i] >= ccgvu_30d_OPE_co2$plus.sigma[i]) & (OPE_30d_smoothed_co2$smoothed[i+1] >= ccgvu_30d_OPE_co2$plus.sigma[i+1]))) {
      OPE_30d_smoothed_co2$pos_SA[i] <- OPE_30d_smoothed_co2$smoothed[i] - ccgvu_30d_OPE_co2$plus.sigma[i] 
      OPE_30d_smoothed_co2$SSE[i] <- OPE_30d_smoothed_co2$smoothed[i]
    } else if (((OPE_30d_smoothed_co2$smoothed[i] <= ccgvu_30d_OPE_co2$minus.sigma[i]) & (OPE_30d_smoothed_co2$smoothed[i-1] <= ccgvu_30d_OPE_co2$minus.sigma[i-1])) |
              ((OPE_30d_smoothed_co2$smoothed[i] <= ccgvu_30d_OPE_co2$minus.sigma[i]) & (OPE_30d_smoothed_co2$smoothed[i+1] <= ccgvu_30d_OPE_co2$minus.sigma[i+1]))) {
      OPE_30d_smoothed_co2$neg_SA[i] <- OPE_30d_smoothed_co2$smoothed[i] - ccgvu_30d_OPE_co2$minus.sigma[i]
      OPE_30d_smoothed_co2$SSE[i] <- OPE_30d_smoothed_co2$smoothed[i]
    } else {
      OPE_30d_smoothed_co2$no_anom[i] <- 0
    }
  }
  
  
  # Make the plots look a bit nicer:
  for (i in 2:nrow(OPE_90d_smoothed_co2)) {
    if ((is.na(OPE_90d_smoothed_co2$pos_SA[i]) & is.na(OPE_90d_smoothed_co2$no_anom[i+1]) & is.na(OPE_90d_smoothed_co2$neg_SA[i+1])) | 
        (is.na(OPE_90d_smoothed_co2$pos_SA[i]) & is.na(OPE_90d_smoothed_co2$no_anom[i-1]) & is.na(OPE_90d_smoothed_co2$neg_SA[i-1]))) {
      OPE_90d_smoothed_co2$pos_SA[i] = 0
    } else if ((is.na(OPE_90d_smoothed_co2$neg_SA[i]) & is.na(OPE_90d_smoothed_co2$no_anom[i+1]) & is.na(OPE_90d_smoothed_co2$pos_SA[i+1])) |  
               (is.na(OPE_90d_smoothed_co2$neg_SA[i]) & is.na(OPE_90d_smoothed_co2$no_anom[i-1]) & is.na(OPE_90d_smoothed_co2$pos_SA[i-1]))) {
        OPE_90d_smoothed_co2$neg_SA[i] = 0 
    } 
  }
  for (i in 2:nrow(OPE_30d_smoothed_co2)) {
    if ((is.na(OPE_30d_smoothed_co2$pos_SA[i]) & is.na(OPE_30d_smoothed_co2$no_anom[i+1]) & is.na(OPE_30d_smoothed_co2$neg_SA[i+1])) | 
        (is.na(OPE_30d_smoothed_co2$pos_SA[i]) & is.na(OPE_30d_smoothed_co2$no_anom[i-1]) & is.na(OPE_30d_smoothed_co2$neg_SA[i-1]))) {
      OPE_30d_smoothed_co2$pos_SA[i] = 0
    } else if ((is.na(OPE_30d_smoothed_co2$neg_SA[i]) & is.na(OPE_30d_smoothed_co2$no_anom[i+1]) & is.na(OPE_30d_smoothed_co2$pos_SA[i+1])) |  
               (is.na(OPE_30d_smoothed_co2$neg_SA[i]) & is.na(OPE_30d_smoothed_co2$no_anom[i-1]) & is.na(OPE_30d_smoothed_co2$pos_SA[i-1]))) {
        OPE_30d_smoothed_co2$neg_SA[i] = 0 
    } 
  }   


## 9: Anomaly plots
  
  # Seasonal anomaly plot
  # Show only growing season
  OPE_90d_smco2.copy <- OPE_90d_smoothed_co2
  OPE_90d_smco2.copy$winter <- NA
  OPE_90d_smco2.copy <- Extract_summer(OPE_90d_smco2.copy)

  # Set x-axis parameters
  x=seq(as.Date("2012-01-01"), as.Date("2021-01-01"), by="years",)
  labels=date_format("%Y")(x)
  breaks=as.POSIXct(sort(x))

  ggplot() +
    geom_line(data=OPE_90d_smco2.copy, aes(x=sampling_datetime, y=no_anom)) +
    geom_line(data=OPE_90d_smco2.copy, aes(x=sampling_datetime, y=pos_SA, color="coral1")) +
    geom_line(data=OPE_90d_smco2.copy, aes(x=sampling_datetime, y=neg_SA, color="turquoise3")) +
    geom_line(data=OPE_90d_smco2.copy, aes(x=sampling_datetime, y=winter), linetype="dotted") +
    scale_x_datetime(limits=range(breaks), labels = labels, breaks = breaks) +
    scale_y_continuous(limits=c(-3.0, 4.6), sec.axis = sec_axis( ~ .*1/2)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "transparent", color = NA),
          plot.background = element_rect(fill = "transparent",color = NA),
          legend.title=element_blank(),legend.position="none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(hjust = 0.98, vjust = -14.1, size = 9, margin=margin(-7,0,-7,0))) + labs(title="OPE")
  
  # Synoptic anomaly plot
  # Show only extended winter season
  OPE_30d_smco2.copy <- OPE_30d_smoothed_co2
  OPE_30d_smco2.copy$summer <- NA
  OPE_30d_smco2.copy <- Extract_winter(OPE_30d_smco2.copy)
  
  # Set x-axis parameters
  x=seq(as.Date("2012-01-01"), as.Date("2021-01-01"), by="years",)
  labels=date_format("%Y")(x)
  breaks=as.POSIXct(sort(x))

  ggplot() +
    geom_line(data=OPE_30d_smco2.copy, aes(x=sampling_datetime, y=no_anom)) +
    geom_line(data=OPE_30d_smco2.copy, aes(x=sampling_datetime, y=pos_SA, color="coral1")) +
    geom_line(data=OPE_30d_smco2.copy, aes(x=sampling_datetime, y=neg_SA, color="turquoise3")) +
    geom_line(data=OPE_30d_smco2.copy, aes(x=sampling_datetime, y=summer), linetype="dotted") +
    scale_x_datetime(limits=range(breaks), labels = labels, breaks = breaks) +
    scale_y_continuous(limits=c(-3.0, 4.6), sec.axis = sec_axis( ~ .*1/2)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "transparent", color = NA),
          plot.background = element_rect(fill = "transparent",color = NA),
          legend.title=element_blank(),legend.position="none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(hjust = 0.98, vjust = -14.1, size = 9, margin=margin(-7,0,-7,0))) + labs(title="OPE")
  
