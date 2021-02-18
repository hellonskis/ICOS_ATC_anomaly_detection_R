CCGCRV_PROG_PATH <- paste("C:/Users/aresovsk/Downloads/ccgcrv", ifelse(.Platform$OS.type == "windows", ".exe", ""), sep = "")

# -############################################################################-
#' Apply ccg curve fitting routines to time series data
#'
#' See \url{ftp://ftp.cmdl.noaa.gov/pub/john/ccgcrv/help/ccgvu.html} for details.
#'
#' @param datatable \code{data.frame} containing the source data.
#' @param datetimeColumnName String. Name of the \code{datatable}'s column containing
#' the datetime data, in a format readable by \code{\link{decimal_date}} function.
#' @param valueColumnName String. Name of the \code{datatable}'s column containing
#' the to-be-fitted data.
#' @param ccgcrvParameters String or string vector. Parameters to be passed to the
#' \code{ccgcrv} script. Default is \code{-all}
#' @param merge Logical. If \code{TRUE}, the ccg curve fitting value will be added
#' to the input \code{datatable}.
#' @param colNames String vector. When custom \code{ccgcrvParameters} are set,
#' user should give the name of the columns outputted by CCGCRV.
#'
#' @return \code{data.table} with curve fitting timeseries
#' @export
#'
#' @examples
#' \dontrun{
#' datatable <- getFlaskData("AMS", "co2")
#' getCcgcrv(datatable, "sampling_datetime", "concentration")
#' getCcgcrv(datatable, "sampling_datetime", "concentration", calculatedValues = c("smooth", "trend"))
#' }
getCcgcrv <- function(datatable, datetimeColumnName, valueColumnName, ccgcrvParameters, calculatedValues, merge = FALSE, colNames) {
  library(lubridate)

  datatable <- copy(datatable)

  if (nrow(datatable) == 0) {
    return(NULL)
  }

  ## Temporarilly rename datetime and value columns
  setnames(datatable, datetimeColumnName, "datetime")
  setnames(datatable, valueColumnName,    "value")
  
  ## Re-order requested columns. CCGCRV output columns order is hard-coded, it
  ## does not depend on input parameter order...
  orderedCcgcrvColumns <- c("func","poly","smooth","trend","detrend","smcycle","harm","res","smres","trres","ressm","gr")
  if (!missing(calculatedValues)) {
    calculatedValues <- as.character(sort(factor(calculatedValues, levels = orderedCcgcrvColumns)))
  }

  ## Format source data.frame
  # Keep only datetime and value columns
  if (ncol(datatable) > 2) { 
    inputDt <- datatable[!is.na(value), .(datetime, value)]} else { 
    inputDt <- data.table(datatable[!is.na(datatable$value),])}
  #inputDt <- data.table(datatable[!is.na(datatable$value),])
  # Convert datetime to decimal date
  inputDt <- inputDt[, datetime := decimal_date(datetime)]
  # Sort by datetime
  inputDt <- inputDt[order(datetime)]

  ## Get curve fitting data.frame
  sourceFileName <- paste("input_", runif(1), ".txt", sep = "")
  write.table(inputDt, file = sourceFileName, col.names = FALSE, row.names = FALSE)
  if (missing(ccgcrvParameters)) {
    if (!missing(calculatedValues)) {
      ccgcrvParameters <- paste("-orig -", paste(calculatedValues, collapse = " -"), sep = "")
      colNames <- c("decimal_date","orig", calculatedValues)
    } else {
      ccgcrvParameters <- "-all"
      colNames <- c("decimal_date","orig", orderedCcgcrvColumns)
    }
  } else {
    colNames <- c("decimal_date","orig", as.character(sort(factor(colNames, levels = orderedCcgcrvColumns)))) ## Re-order user-defined colNames to match CCGCRV hard-coded colNames
  }
  command        <- paste(CCGCRV_PROG_PATH, paste(ccgcrvParameters), sourceFileName)
  outputString   <- system(command, intern = TRUE)
  file.remove(sourceFileName)
  
  ## Convert the ccgcrv output to data.frame
  ccgcrvDt <- fread(input     = paste(outputString, collapse = "\n"),
                    col.names = colNames)
  ccgcrvDt <- ccgcrvDt[, datetime := date_decimal(decimal_date, tz = "GMT")]
  ccgcrvDt <- ccgcrvDt[, datetime := numericToPosix(round(as.numeric(datetime)))] # Remove small difference due to double decimal conversion
  ccgcrvDt <- ccgcrvDt[, decimal_date := NULL]

  ## Merge ccgcrv output an dinput data.frame
  if (merge) {
    ccgcrvDt <- merge(ccgcrvDt, datatable, by = "datetime", all = TRUE)
  }

  ## Give datetime and value columns their original name
  setnames(ccgcrvDt, "datetime", datetimeColumnName)

  return(ccgcrvDt)
}


################################################################################################################

Ccgcrv_pm_baseline <- function(data, span) {
  
  afternoon <- subset(data, (hour(data$sampling_datetime) >= 14) & (hour(data$sampling_datetime) <= 17))
  ## Remove values with std_dev > 0.5 ppm
  for (i in 1:nrow(afternoon)) {
    if (!is.na(afternoon$std_dev[i])) {
      if (afternoon$std_dev[i] > 0.5) {
        afternoon$concentration[i] = NA
      }
    }
  }
  
  ## Aggregate to daily
  afternoon$Date  <- as.Date(afternoon$sampling_datetime, "%Y-%m-%d", tz="GMT")
  ccg <- aggregate(concentration ~ Date, afternoon, mean, na.rm=TRUE, na.action=na.pass)
  ccg <- NaN.to.na(ccg)
  colnames(ccg) <- c("sampling_datetime", "concentration")
  ccg$sampling_datetime <- as.POSIXct(ccg$sampling_datetime, format="%Y-%m-%d", tz="GMT")
  ccg$sampling_datetime <- round(ccg$sampling_datetime, units = "days")
  ccg$sampling_datetime <- as.POSIXct(ccg$sampling_datetime, format="%Y-%m-%d", tz="GMT")
  ccg$concentration <- na_interpolation(ccg$concentration)
  
  colNames <- c("func","poly","smooth","trend","detrend","smcycle","harm","res","smres","trres","ressm","gr")
  ccgvu <- getCcgcrv(ccg, "sampling_datetime", "concentration", ccgcrvParameters=paste("-all -short\ ",span), merge=FALSE, colNames=colNames) 
  ccgvu$sampling_datetime <- round_date(ccg$sampling_datetime, unit = "day")
  
  return(ccgvu)
  
}


################################################################################################################

Ccgcrv_am_baseline <- function(data, span) {
  
  morning <- subset(data, (hour(data$sampling_datetime) >= 2) & (hour(data$sampling_datetime) <= 5))
  ## Remove values with std_dev > 0.5 ppm
  for (i in 1:nrow(morning)) {
    if (!is.na(morning$std_dev[i])) {
      if (morning$std_dev[i] > 0.5) {
        morning$concentration[i] = NA
      }
    }
  }
  
  ## Aggregate to daily
  morning$Date  <- as.Date(morning$sampling_datetime, "%Y-%m-%d", tz="GMT")
  ccg <- aggregate(concentration ~ Date, morning, mean, na.rm=TRUE, na.action=na.pass)
  ccg <- NaN.to.na(ccg)
  colnames(ccg) <- c("sampling_datetime", "concentration")
  ccg$sampling_datetime <- as.POSIXct(ccg$sampling_datetime, format="%Y-%m-%d", tz="GMT")
  ccg$sampling_datetime <- round(ccg$sampling_datetime, units = "days")
  ccg$sampling_datetime <- as.POSIXct(ccg$sampling_datetime, format="%Y-%m-%d", tz="GMT")
  ccg$concentration <- na_interpolation(ccg$concentration)
  
  colNames <- c("func","poly","smooth","trend","detrend","smcycle","harm","res","smres","trres","ressm","gr")
  ccgvu <- getCcgcrv(ccg, "sampling_datetime", "concentration", ccgcrvParameters=paste("-all -short\ ",span), merge=FALSE, colNames=colNames) 
  ccgvu$sampling_datetime <- round_date(ccg$sampling_datetime, unit = "day")
  
  return(ccgvu)
  
}
  

  
 


  
   

  
