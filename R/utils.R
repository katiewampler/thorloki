utils::globalVariables(c("site", "time", "value"))


#' Format date for using in file names
#'
#' @param date a date in the form YYYY-MM-DD, default is today's date
#' @param yeardigit either 2 or 4 to indicate number of year digits desired
#' @importFrom lubridate year month day
#' @importFrom stringr str_pad
#' @return a string of the date with underscores between the year, month, and day
#' @export
#'
#' @examples
#' file_date("2023-02-13")

file_date <- function(date = Sys.Date(), yeardigit=4){
  stopifnot(yeardigit %in% c(2,4))
  if(yeardigit == 2){
    year <- stringr::str_sub(as.character(lubridate::year(date)), start=-2, end=-1)
  }else{year <- lubridate::year(date)}
  month <- stringr::str_pad(lubridate::month(date), width= 2, side="left", "0")
  day <- stringr::str_pad(lubridate::day(date), width= 2, side="left", "0")
  new_date <- paste(year, month, day, sep="_")
  return(new_date)
}

#' Convert date to one used for samples
#'
#' @param date a date in the form YYYY-MM-DD, default is today's date
#' @param yeardigit either 2 or 4 to indicate number of year digits desired
#'
#' @return a character in the form of YYMMDD or YYYYMMDD used for logging samples
#' @export
#' @importFrom lubridate year month day
#' @importFrom stringr str_pad
#' @examples
#' samp_date("2023-02-13")

samp_date <- function(date=Sys.Date(), yeardigit=2){
  stopifnot(yeardigit %in% c(2,4))
  if(yeardigit == 2){
    year <- stringr::str_sub(as.character(lubridate::year(date)), start=-2, end=-1)
  }else{year <- lubridate::year(date)}
  month <- stringr::str_pad(lubridate::month(date), width= 2, side="left", "0")
  day <- stringr::str_pad(lubridate::day(date), width= 2, side="left", "0")
  new_date <- paste(year, month, day, sep="")
  return(new_date)
}

#' Calculate percent change
#'
#' @param old the original value
#' @param new the new value
#'
#' @return the percent change between the old and new values as a percent (%)
#' @export
#'
#' @examples
#' per_change(1,2)
#'
per_change <- function(old, new){
  stopifnot(is.numeric(c(old, new)))
  change <- (new-old)/old * 100
  return(change)
}

#' Split up a string into pieces
#'
#' A wrapper for str_split from stringr that will return a vector if
#' the string has length 1
#'
#' @importFrom stringr str_split_fixed
#' @param string Input vector. Either a character vector, or something coercible to one.
#' @param pattern Pattern to look for.
#' @param n Maximum number of pieces to return. Default (Inf) uses all possible split positions.
#' @param piece A vector of numbers, if it's not NA it will return that piece of each vector
#' @return if the string is length 1 or 1 piece from each string is requested, it will
#' return a vector. If multiple pieces of multiple strings is requested
#' it will create a character matrix.
#' @export
#'
#' @examples
#' string <- "hello/world/its/katie"
#' str_split2(string, "/")
#' str_split2(string, "/", piece=c(2:4))
#' string2 <- c("loki_is_here", "thor_is_here_too")
#' str_split2(string2,"_")
#' str_split2(string2,"_", piece=c(1,3))
#'
str_split2 <- function(string, pattern, n=Inf, piece=NA){
  val <- stringr::str_split_fixed(string, pattern=pattern, n=n)
  if(sum(is.na(piece))==0){
    val <- val[,piece]
  }
  if(length(string) == 1){
    val <- as.vector(val)}

  return(val)
}


