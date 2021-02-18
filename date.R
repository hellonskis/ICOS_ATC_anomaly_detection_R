
getTimeAgo <- function(lengthOfTime, startDate = Sys.Date()) {
  timeAgo <- seq(startDate, length = 2, by = paste("-", lengthOfTime, sep = ""))[2]
  timeAgo <- as.POSIXct(as.character.Date(timeAgo), tz = "GMT")
  return(timeAgo)
}

getToday <- function() {
  today <- as.POSIXct(as.character.Date(Sys.Date()), tz = "GMT")
  return(today)
}

# -############################################################################-
#' Get a date older than any recorded data
#'
#' This function gives a numerical equivalent to ICOS' \code{FD} (First Date) tag.
#'
#' @return \code{POSIXct}
#' @export
#'
#' @examples getFirstAvailableDate()
getFirstAvailableDate <- function() {
  return(as.POSIXct("1950-01-01 00:00:00", tz = "GMT"))
}

# -############################################################################-
#' Converts numeric dates to POSIXct dates
#'
#' @param numericDates Integer vector. Number of seconds since 1970-01-01 00:00:00
#'
#' @return Vector of \code{POSIXct}  dates
#' @export
#'
#' @examples numericToPosix(c(1073520000, 1249466400))
numericToPosix <- function(numericDates) {
  posixDates <- as.POSIXct(x      = numericDates,
                           origin = as.POSIXct("1970-01-01", tz = "GMT"),
                           tz     = "GMT")
  return(posixDates)
}

getLastMomentOfDay <- function(datetime) {
  return(as.POSIXct(paste(format(datetime, "%Y-%m-%d"), "23:59:59"), tz = "GMT"))
}


# -############################################################################-
#' Get the "light-related" daily period night/dawn/day/dusk for each local hour
#'
#' @param dayNightOnly Boolean. If \code{TRUE}, return only day and night hours.
#'
#' @return \code{data.table}
#' @export
#'
#' @examples
#' getDailyPeriodsDt()
#' getDailyPeriodsDt(dayNightOnly = TRUE)
getDailyPeriodsDt <- function(dayNightOnly = FALSE) {
  dailyPeriodsDt <- data.table(
    localHour  = c(seq(0,4),        seq(5,11),      seq(12,16),    seq(17,23)),
    periodName = c(rep("night", 5), rep("dawn", 7), rep("day", 5), rep("dusk", 7))
    )

  if (dayNightOnly) {
    dailyPeriodsDt <- dailyPeriodsDt[periodName %in% c("night", "day")]
  }
  return(dailyPeriodsDt)
}

# -############################################################################-
#' Add the "light-related" daily period night/dawn/day/dusk to input \code{datatable}
#'
#' @param datatable \code{data.table} containing at least a column with local hours
#' as numeric.
#' @param localHourColumn  String. Name of the column containing the local hours
#' as numeric.
#' @param dayNightOnly Boolean. If \code{TRUE}, return only day and night hours.
#'
#' @return \code{data.table}. Same as input, with an additional \code{periodName}
#' column.
#' @export
#'
#' @examples
#' datatable <- data.table(value = rnorm(24), localHour = seq(0,23))
#' addDailyPeriods(datatable)
#' addDailyPeriods(datatable, dayNightOnly = TRUE)
addDailyPeriods <- function(datatable, localHourColumn = "localHour", dayNightOnly = FALSE) {
  datatable <- merge(x     = datatable,
                     by.x  = localHourColumn,
                     all.x = TRUE,
                     y     = getDailyPeriodsDt(dayNightOnly),
                     by.y  = "localHour")
  return(datatable)
}

# -############################################################################-
#' Title
#'
#' /!\ For a use with data.table, use getSeason() function instead.
#'
#' @param dataframe
#' @param site
#' @param gmtColumnName
#' @param localDatetimeColumn
#'
#' @return
#' @export
#'
#' @examples
addLocalDatetime <- function(dataframe, site, gmtColumnName = "sampling_datetime", localDatetimeColumn = "local_datetime") {
  timeshiftSecond <- getSiteTimezone(site)
  dataframe[, localDatetimeColumn] <- dataframe[, gmtColumnName] + timeshiftSecond
  return(dataframe)
}

# -############################################################################-
#' Add a \code{localDatetime} column to input \code{datatable}
#'
#' @param datatable \code{data.table}
#' @param site
#' @param gmtColumnName
#' @param addDailyPeriod Logical.
#' @param dayNightOnly Logical.
#' @param addLocalHour Logical.
#'
#' @return
#' @export
#'
#' @examples
addLocalTime <- function(datatable, site, gmtColumnName, addDailyPeriod = FALSE, dayNightOnly = FALSE, addLocalHour = FALSE) {
  timeshiftSecond <- getSiteTimezone(site)
  datatable <- datatable[, localDatetime := get(gmtColumnName) + timeshiftSecond]
  datatable <- datatable[, localHour := hour(localDatetime)]
  if (addDailyPeriod) {
    datatable <- addDailyPeriods(datatable, localHourColumn = "localHour", dayNightOnly = dayNightOnly)
  }
  if (!addLocalHour) {
    datatable <- datatable[, localHour := NULL]
  }
  return(datatable)
}

addSunriset <- function(datatable, site, timezone = c("local", "gmt"), localDatetimeColumn = "local_datetime") {
  latLonDt        <- getSiteLatLon(site)
  geoposition     <- matrix(c(latLonDt$lon, latLonDt$lat), nrow = 1)
  timeshiftSecond <- getSiteTimezone(site)

  ## Sunrise
  sunriseDt <- data.table(maptools::sunriset(
    crds        = geoposition,
    dateTime    = datatable[, get(localDatetimeColumn)],
    direction   = "sunrise",
    POSIXct.out = TRUE))
  if ("gmt" %in% timezone) {
    datatable <- datatable[, sunriseGMT := sunriseDt$time]
  }
  if ("local" %in% timezone) {
    datatable <- datatable[, sunriseLocal := sunriseDt$time + timeshiftSecond]
  }

  ## Sunset
  sunsetDt <- data.table(maptools::sunriset(
    crds        = geoposition,
    dateTime    = datatable[, get(localDatetimeColumn)],
    direction   = "sunset",
    POSIXct.out = TRUE))
  if ("gmt" %in% timezone) {
    datatable <- datatable[, sunsetGMT := sunsetDt$time]
  }
  if ("local" %in% timezone) {
    datatable <- datatable[, sunsetLocal := sunsetDt$time + timeshiftSecond]
  }

  return(datatable)
}

# -############################################################################-
#' Add a \code{season} column to the input \code{dataframe}
#'
#' /!\ For a use with data.table, use getSeason() function instead.
#'
#' @param dataframe \code{data.frame} containing at least a \code{POSIXct} column
#' @param datetimeColumn String. Name of the column containing the \code{POSIXct}
#' data from which the season should be extracted
#'
#' @return Same as \code{dataframe}, with an additinal \code{season} column
#' @export
#'
#' @examples
#' dataframe <- data.frame(datetime = c(as.POSIXct("2016-01-01"), as.POSIXct("2016-06-30")))
#' dataframe <- addSeason(dataframe, "datetime")
addSeason <- function(dataframe, datetimeColumn){
  datetime <- dataframe[, datetimeColumn]
  datetime <- 100*month(datetime) + day(datetime)
  ## input Seasons upper limits in the form MMDD in the "break =" option:
  cuts <- base::cut(datetime, breaks = c(0,321,621,0921,1221,1231))
  # rename the resulting groups (could've been done within cut(...levels=) if "Winter" wasn't double
  levels(cuts) <- c("Winter", "Spring", "Summer", "Autumn", "Winter")
  dataframe$season <- cuts
  return(dataframe)
}

# -############################################################################-
#' Get the season of the input \code{datetimes}
#'
#' @param datetimes
#'
#' @return
#' @export
#'
#' @examples
getSeason <- function(datetimes) {
  ## Convert dates to a MMDD-formatted numeric value
  datetimes <- 100*month(datetimes) + day(datetimes)
  ## "Cut" by season
  seasonsBoundaries <- c(0,321,621,0921,1221,1231)
  cuts <- base::cut(datetimes, breaks = seasonsBoundaries)
  # Rename the resulting groups (could've been done within cut(...levels=) if "Winter" wasn't double
  levels(cuts) <- c("Winter", "Spring", "Summer", "Autumn", "Winter")
  return(cuts)
}

# -############################################################################-
#' Add rows with NA values where no data is available for given datetime
#'
#' This function introduces NA in timeseries where data is missing to avoid that
#' \link{\code{ggplot()}} draws a line between the limits of the "data hole".
#'
#' @param datatable \code{data.table} containing data to be filled
#' @param expectedDatetimes Vector of \code{POSIXct}. All datetimes where a value
#' is expected.
#' @param columnGroupsDt \code{data.table} containing only the grouping columns.
#' For each line of this data.table, and for each \code{expectedDatetimes}, a value
#' (either a NA or the value already existing in \code{datatable}) is expected in output.
#' @param datetimeColumn String. Name of the column containing the datetimes
#'
#' @return \code{data.table}
#' @export
#'
#' @examples \dontrun{"NO EXAMPLE"}
fillWithMissingDatetime <- function(datatable, expectedDatetimes, columnGroupsDt, datetimeColumn) {
  expectedDatetimesDt <- columnGroupsDt[rep(seq(1, nrow(columnGroupsDt)), length(expectedDatetimes))]
  expectedDatetimesDt <- expectedDatetimesDt[, (datetimeColumn) := sort(rep(expectedDatetimes, nrow(columnGroupsDt)))]
  filledDatatable     <- merge(datatable, expectedDatetimesDt, all.y = TRUE, by = intersect(names(datatable), names(expectedDatetimesDt)))
  setcolorder(filledDatatable, names(datatable))
  return(filledDatatable)
}