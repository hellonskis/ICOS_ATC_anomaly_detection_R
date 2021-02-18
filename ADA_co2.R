options(java.parameters = "- Xmx1024m")
source("C:/Users/aresovsk/Documents/R/main.R")

#############################################################################

site = "OPE"

## 1. Read in the 2011-2020 data and subset dates
  hist_df <- readRDS('OPE_2011-2019_co2_120m.rds')
  hist_df$sampling_datetime <- as.POSIXct(hist_df$sampling_datetime, format="%Y-%m-%d %H:%M:%OS", tz="GMT")
  hist_df <- hist_df[,c(11,13)]
  hist_df <- hist_df[hist_df$sampling_datetime < "2019-05-01",]
  
  L2_df <- readRDS('OPE_2011-2020_co2_120m.rds')
  L2_df$sampling_datetime <- as.POSIXct(L2_df$sampling_datetime, format="%Y-%m-%d %H:%M:%OS", tz="GMT")
  L2_df <- L2_df[L2_df$sampling_datetime >= "2019-05-01",]
  L2_df <- L2_df[,c(1,2)]

  hist_co2 <- rbind(hist_df, L2_df)
  rm(hist_df, L2_df)
  
  ## Fix missing dates
  hourly <- data.frame(seq(hist_co2$sampling_datetime[1], hist_co2$sampling_datetime[nrow(hist_co2)], by="hour"))
  colnames(hourly) <- c("sampling_datetime")
  hourly$concentration <- NA
  df <- hourly
  for (i in 1:nrow(df)) {  
    print(i)
    if (df$sampling_datetime[i] %in% hist_co2$sampling_datetime) {
      df$concentration[i] <- hist_co2$concentration[which(df$sampling_datetime[i] == hist_co2$sampling_datetime)]
    }
  }
  hist_co2 <- df


## 2. Read in NRT data
  NRT_df <- read.xlsx('ICOS_ATC_NRT_OPE_2020-06-01_2020-09-27_120.0_CO2.xlsx',
                           sheetIndex = "hourly", header=TRUE)
  NRT_df$sampling_datetime <- as.POSIXct(NRT_df$sampling_datetime, format="%Y-%m-%d %H:%M:%OS", tz="GMT")
  NRT_df <- NRT_df[,c(1,2)]
  
  ## Fix missing dates
  hourly <- data.frame(seq(as.POSIXct((hist_co2$sampling_datetime[nrow(hist_co2)] +3600), tz="GMT"), 
                                    as.POSIXct(NRT_df$sampling_datetime[nrow(NRT_df)], tz="GMT"), by="hour"))
  colnames(pastyear.hourly) <- c("sampling_datetime")
  pastyear.hourly$concentration <- NA
  df <- pastyear.hourly
  for (i in 1:nrow(df)) {
    print(i)
    if (df$sampling_datetime[i] %in% NRT_df$sampling_datetime) {
      df$concentration[i] <- NRT_df$concentration[which(df$sampling_datetime[i] == NRT_df$sampling_datetime)]
    }
  }
  NRT_df <- df

  ## Concatenate
  df_co2 <- rbind(hist_co2, NRT_df)
  ## Use only afternoon values (non-mountain sites)
  df_co2 <- subset(df_co2, (hour(as.POSIXct(df_co2$sampling_datetime)) >= 12) & 
                           (hour(as.POSIXct(df_co2$sampling_datetime)) <= 17)) # or
  ## Use only nighttime values (mountain sites)
  #df_co2 <- subset(df_co2, (hour(as.POSIXct(df_co2$sampling_datetime)) >= 20) | 
  #                         (hour(as.POSIXct(df_co2$sampling_datetime)) <= 5))


## 3. Aggregate hourly to daily
  df_co2$Date  <- as.Date(df_co2$sampling_datetime, "%Y-%m-%d", tz="GMT")
  df_co2_daily <- aggregate(concentration ~ Date, df_co2, mean, na.rm=TRUE, na.action=na.pass)
  df_co2_daily <- NaN.to.na(df_co2_daily)
  colnames(df_co2_daily) <- c("sampling_datetime", "concentration")
  df_co2_daily$sampling_datetime <- as.POSIXct(df_co2_daily$sampling_datetime, format="%Y-%m-%d", tz="GMT")
  df_co2_daily$sampling_datetime <- round(df_co2_daily$sampling_datetime, units = "days")
  df_co2_daily$sampling_datetime <- as.POSIXct(df_co2_daily$sampling_datetime, format="%Y-%m-%d", tz="GMT")
  df_co2_daily$concentration <- na_interpolation(df_co2_daily$concentration)   
  ## Remove leap days
  df_co2_daily <- subset(df_co2_daily, df_co2_daily$sampling_datetime %in% 
                         df_co2_daily$sampling_datetime[!find_leap(df_co2_daily$sampling_datetime)])
  

## 4. Extract the 90-day CCGvu baseline curve
  colNames <- c("func","poly","smooth","trend","detrend","smcycle","harm","res","smres","trres","ressm","gr")
  ccgvu_90d_co2 <- getCcgcrv(df_co2_daily, "sampling_datetime", "concentration", 
                             ccgcrvParameters="-all -short\ 90", merge=FALSE, colNames=colNames)
  ccgvu_30d_co2 <- getCcgcrv(df_co2_daily, "sampling_datetime", "concentration", 
                             ccgcrvParameters="-all -short\ 30", merge=FALSE, colNames=colNames)
  
  ## Round dates
  ccgvu_90d_co2$sampling_datetime <- round_date(ccgvu_90d_co2$sampling_datetime, unit = "day")
  ccgvu_30d_co2$sampling_datetime <- round_date(ccgvu_30d_co2$sampling_datetime, unit = "day")


## 5. Calculate delta-C values
  
  ccgvu_90d_co2$dC <- NA
  for (i in 1:nrow(ccgvu_90d_o2)) {
    ccgvu_90d_co2$dC[i] <- ccgvu_90d_co2$smcycle[i] - ccgvu_90d_co2$harm[i]
  } 
  ccgvu_30d_co2$dC <- NA
  for (i in 1:nrow(ccgvu_30d_co2)) {
    ccgvu_30d_co2$dC[i] <- ccgvu_30d_co2$smcycle[i] - ccgvu_30d_co2$harm[i]
  } 


## 6. Calculate season-adjusted sigma values 
  ccgvu_90d_co2$sd <- seasonal_sigma(ccgvu_90d_co2, ndays=(nrow(ccgvu_90d_co2)-1))
  ccgvu_30d_co2$sd <- seasonal_sigma(ccgvu_30d_co2, ndays=(nrow(ccgvu_30d_co2)-1))
  
#  # Use 2*sd ? 
  for (i in 1:nrow(ccgvu_90d_co2)) {
    ccgvu_90d_co2$sd[i] <- ccgvu_90d_co2$sd[i] * 2
  }
  for (i in 1:nrow(ccgvu_30d_co2)) {
    ccgvu_30d_co2$sd[i] <- ccgvu_30d_co2$sd[i] * 2
  }
  
  ccgvu_90d_co2$plus.sigma <- NA
  ccgvu_90d_o2$minus.sigma <- NA
  ccgvu_30d_co2$plus.sigma <- NA
  ccgvu_30d_co2$minus.sigma <- NA
  for (i in 1:nrow(ccgvu_90d_co2)) {
    ccgvu_90d_co2$plus.sigma[i]  <- ccgvu_90d_co2$harm[i] + ccgvu_90d_co2$trend[i] + ccgvu_90d_co2$sd[i]
    ccgvu_90d_co2$minus.sigma[i] <- ccgvu_90d_co2$harm[i] + ccgvu_90d_co2$trend[i] - ccgvu_90d_co2$sd[i]
    ccgvu_30d_co2$plus.sigma[i] <- ccgvu_30d_co2$harm[i] + ccgvu_30d_co2$trend[i] + ccgvu_30d_co2$sd[i]
    ccgvu_30d_co2$minus.sigma[i] <- ccgvu_30d_co2$harm[i] + ccgvu_30d_co2$trend[i] - ccgvu_30d_co2$sd[i]
  }

  ## Make smoothed curves
  df_90d_smoothed_co2 <- loess.sm(ccgvu_90d_co2, "sampling_datetime", "orig", span=(90/nrow(ccgvu_90d_co2)))
  df_30d_smoothed_co2 <- loess.sm(ccgvu_30d_co2, "sampling_datetime", "orig", span=(30/nrow(ccgvu_30d_co2)))
  
  
## 7. Identify seasonal anomalies 
  df_90d_smoothed_co2$no_anom <- NA
  df_90d_smoothed_co2$pos_SA <- NA
  df_90d_smoothed_co2$neg_SA <- NA
  df_90d_smoothed_co2$SSE <- NA
  for (i in 1:nrow(df_90d_smoothed_co2)) {
    if ((df_90d_smoothed_co2$smoothed[i] >= ccgvu_90d_co2$plus.sigma[i])) {
      df_90d_smoothed_co2$pos_SA[i] <- df_90d_smoothed_co2$smoothed[i] - ccgvu_90d_co2$plus.sigma[i]
      df_90d_smoothed_co2$SSE[i] <- df_90d_smoothed_co2$smoothed[i]
    } else if ((df_90d_smoothed_co2$smoothed[i] <= ccgvu_90d_co2$minus.sigma[i])) {
      df_90d_smoothed_co2$neg_SA[i] <- df_90d_smoothed_co2$smoothed[i] - ccgvu_90d_co2$minus.sigma[i]
      df_90d_smoothed_co2$SSE[i] <- df_90d_smoothed_co2$smoothed[i]
    } else {
      df_90d_smoothed_co2$no_anom[i] <- 0
    }
  }
  df_30d_smoothed_co2$no_anom <- NA
  df_30d_smoothed_co2$pos_SA <- NA
  df_30d_smoothed_co2$neg_SA <- NA
  df_30d_smoothed_co2$SSE <- NA
  for (i in 2:(nrow(df_30d_smoothed_co2)-1)) {
    if (((df_30d_smoothed_co2$smoothed[i] >= ccgvu_30d_co2$plus.sigma[i]) & (df_30d_smoothed_co2$smoothed[i-1] >= ccgvu_30d_co2$plus.sigma[i-1])) |
       ((df_30d_smoothed_co2$smoothed[i] >= ccgvu_30d_co2$plus.sigma[i]) & (df_30d_smoothed_co2$smoothed[i+1] >= ccgvu_30d_co2$plus.sigma[i+1]))) {
      df_30d_smoothed_co2$pos_SA[i] <- df_30d_smoothed_co2$smoothed[i] - ccgvu_30d_co2$plus.sigma[i] 
      df_30d_smoothed_co2$SSE[i] <- df_30d_smoothed_co2$smoothed[i]
    } else if (((df_30d_smoothed_co2$smoothed[i] <= ccgvu_30d_co2$minus.sigma[i]) & (df_30d_smoothed_co2$smoothed[i-1] <= ccgvu_30d_co2$minus.sigma[i-1])) |
              ((df_30d_smoothed_co2$smoothed[i] <= ccgvu_30d_co2$minus.sigma[i]) & (df_30d_smoothed_co2$smoothed[i+1] <= ccgvu_30d_co2$minus.sigma[i+1]))) {
      df_30d_smoothed_co2$neg_SA[i] <- df_30d_smoothed_co2$smoothed[i] - ccgvu_30d_co2$minus.sigma[i]
      df_30d_smoothed_co2$SSE[i] <- df_30d_smoothed_co2$smoothed[i]
    } else {
      df_30d_smoothed_co2$no_anom[i] <- 0
    }
  }
  
  
  # Make the plots look a bit nicer:
  for (i in 2:nrow(df_90d_smoothed_co2)) {
    if ((is.na(df_90d_smoothed_co2$pos_SA[i]) & is.na(df_90d_smoothed_co2$no_anom[i+1]) & is.na(df_90d_smoothed_co2$neg_SA[i+1])) | 
        (is.na(df_90d_smoothed_co2$pos_SA[i]) & is.na(df_90d_smoothed_co2$no_anom[i-1]) & is.na(df_90d_smoothed_co2$neg_SA[i-1]))) {
      df_90d_smoothed_co2$pos_SA[i] = 0
    } else if ((is.na(df_90d_smoothed_co2$neg_SA[i]) & is.na(df_90d_smoothed_co2$no_anom[i+1]) & is.na(df_90d_smoothed_co2$pos_SA[i+1])) |  
               (is.na(df_90d_smoothed_co2$neg_SA[i]) & is.na(df_90d_smoothed_co2$no_anom[i-1]) & is.na(df_90d_smoothed_co2$pos_SA[i-1]))) {
        df_90d_smoothed_co2$neg_SA[i] = 0 
    } 
  }
  for (i in 2:nrow(df_30d_smoothed_co2)) {
    if ((is.na(df_30d_smoothed_co2$pos_SA[i]) & is.na(df_30d_smoothed_co2$no_anom[i+1]) & is.na(df_30d_smoothed_co2$neg_SA[i+1])) | 
        (is.na(df_30d_smoothed_co2$pos_SA[i]) & is.na(df_30d_smoothed_co2$no_anom[i-1]) & is.na(df_30d_smoothed_co2$neg_SA[i-1]))) {
      df_30d_smoothed_co2$pos_SA[i] = 0
    } else if ((is.na(df_30d_smoothed_co2$neg_SA[i]) & is.na(df_30d_smoothed_co2$no_anom[i+1]) & is.na(df_30d_smoothed_co2$pos_SA[i+1])) |  
               (is.na(df_30d_smoothed_co2$neg_SA[i]) & is.na(df_30d_smoothed_co2$no_anom[i-1]) & is.na(df_30d_smoothed_co2$pos_SA[i-1]))) {
        df_30d_smoothed_co2$neg_SA[i] = 0 
    } 
  }   


## 9: Anomaly plots
  
  # Seasonal anomaly plot
  # Show only growing season
  df_90d_smco2.copy <- df_90d_smoothed_co2
  df_90d_smco2.copy$winter <- NA
  df_90d_smco2.copy <- Extract_summer(df_90d_smco2.copy)

  # Set x-axis parameters
  x=seq(as.Date("2012-01-01"), as.Date("2021-01-01"), by="years",)
  labels=date_format("%Y")(x)
  breaks=as.POSIXct(sort(x))

  ggplot() +
    geom_line(data=df_90d_smco2.copy, aes(x=sampling_datetime, y=no_anom)) +
    geom_line(data=df_90d_smco2.copy, aes(x=sampling_datetime, y=pos_SA, color="coral1")) +
    geom_line(data=df_90d_smco2.copy, aes(x=sampling_datetime, y=neg_SA, color="turquoise3")) +
    geom_line(data=df_90d_smco2.copy, aes(x=sampling_datetime, y=winter), linetype="dotted") +
    scale_x_datetime(limits=range(breaks), labels = labels, breaks = breaks) +
    scale_y_continuous(limits=c(-3.0, 4.6), sec.axis = sec_axis( ~ .*1/2)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "transparent", color = NA),
          plot.background = element_rect(fill = "transparent",color = NA),
          legend.title=element_blank(),legend.position="none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(hjust = 0.98, vjust = -14.1, size = 9, margin=margin(-7,0,-7,0))) + labs(title=paste(site))
  
  # Synoptic anomaly plot
  # Show only extended winter season
  df_30d_smco2.copy <- df_30d_smoothed_co2
  df_30d_smco2.copy$summer <- NA
  df_30d_smco2.copy <- Extract_winter(df_30d_smco2.copy)
  
  # Set x-axis parameters
  x=seq(as.Date("2012-01-01"), as.Date("2021-01-01"), by="years",)
  labels=date_format("%Y")(x)
  breaks=as.POSIXct(sort(x))

  ggplot() +
    geom_line(data=df_30d_smco2.copy, aes(x=sampling_datetime, y=no_anom)) +
    geom_line(data=df_30d_smco2.copy, aes(x=sampling_datetime, y=pos_SA, color="coral1")) +
    geom_line(data=df_30d_smco2.copy, aes(x=sampling_datetime, y=neg_SA, color="turquoise3")) +
    geom_line(data=df_30d_smco2.copy, aes(x=sampling_datetime, y=summer), linetype="dotted") +
    scale_x_datetime(limits=range(breaks), labels = labels, breaks = breaks) +
    scale_y_continuous(limits=c(-3.0, 4.6), sec.axis = sec_axis( ~ .*1/2)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "transparent", color = NA),
          plot.background = element_rect(fill = "transparent",color = NA),
          legend.title=element_blank(),legend.position="none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(hjust = 0.98, vjust = -14.1, size = 9, margin=margin(-7,0,-7,0))) + labs(title=paste(site))
  
