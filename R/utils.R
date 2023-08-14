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


#' Get season from date
#'
#' Spring is defined as months March through May, Summer is months June through
#' August, Fall is September through November, and Winter is Decemeber through
#' February.
#'
#' @param Date a vector of dates
#'
#' @return a vector of seasons
#' @export
#'
#' @examples
#' dates <- c("2023-02-13", "2023-07-16")
#' as.Date(dates)
#' season(dates)
season <- function(Date){
  month <- data.frame(dates = Date, month = month(Date))
  month$season <- NA
  month$season[month$month %in% 3:5] <- "Spring"
  month$season[month$month %in% 6:8] <- "Summer"
  month$season[month$month %in% 9:11] <- "Fall"
  month$season[month$month %in% c(12,1,2)] <- "Winter"

  return(month$season)
}

#' Katie's personal ggplot theme
#'
#' Keeps with the green theme branding of my work. Uses Barlow Semi Condensed font
#' and used the green colors from google slides.
#'
#' @details
#' Make sure the font is installed from here: (https://fonts.google.com/specimen/Barlow+Semi+Condensed)
#' Then Just once you'll need to run the following code:
#' font_add("Barlow Semi Condensed", "BarlowSemiCondensed-Regular.ttf")
#' font_add("Barlow Semi Condensed ExtraBold", "BarlowSemiCondensed-ExtraBold.ttf")
#' font_add("Barlow Semi Condensed SemiBold", "BarlowSemiCondensed-SemiBold.ttf")
#'
#' @importFrom showtext showtext_auto
#' @importFrom sysfonts font_add
#' @importFrom ggthemes theme_clean
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_point geom_smooth facet_wrap
#' @export
#' @examples
#' x <- sample(1:250, 50)
#' y <- sample(c("A","B", "C"), 50, replace=T)
#' z <- sample(c("LOW","MED","HIGH"), 50, replace=T)
#' w <- sample(1:250, 50)
#' df <- data.frame(x=x, y=y, z=z, w=w)
#'
#' ggplot(df, aes(x=y, y=x, fill=z)) + geom_boxplot() + theme_green()
#' ggplot(df, aes(x=x, y=w, color=y)) + geom_point() + geom_smooth(method="lm") + theme_green()
#'
theme_green <- function(){
  #set up fonts
  sysfonts::font_add("Barlow Semi Condensed", "BarlowSemiCondensed-Regular.ttf")
  sysfonts::font_add("Barlow Semi Condensed ExtraBold", "BarlowSemiCondensed-ExtraBold.ttf")
  sysfonts::font_add("Barlow Semi Condensed SemiBold", "BarlowSemiCondensed-SemiBold.ttf")
  showtext::showtext_auto()
  font <- "Barlow Semi Condensed"
  bold <- "Barlow Semi Condensed SemiBold"
  exbold <- "Barlow Semi Condensed ExtraBold"

  ggthemes::theme_clean() %+replace%
  ggplot2::theme(
    #change background colors
    plot.background = element_rect(fill = "#d9ead3", colour = "#d9ead3"),
    legend.background = element_rect(fill = "#d9ead3", colour = "#d9ead3"),
    panel.background = element_rect(fill = "#eef6eb", colour = "#d9ead3"),
    legend.key = element_rect(fill = "#d9ead3", colour = "#d9ead3"),
    strip.background =element_rect(fill="gray60"),
    axis.line.x.bottom = element_line(colour = "#274e13"),
    axis.line.y.left = element_line(colour ="#274e13"),


    #change text
    axis.title = element_text(family=exbold, size=40, color="#38761dff"),
    axis.text =  element_text(family=font, size=30, color="#274e13"),
    legend.title = element_text(family=exbold, size=40, color="#38761dff"),
    legend.text = element_text(family=font, size=30, color="#274e13"),
    strip.text = element_text(family=bold, size=40, color="#274e13")
  )
}


#' ggplot theme for publication figures
#'
#' Used to help ensure consistent formatting across publication figures. Based on
#' ggclean with a few modifications.
#'
#' @importFrom ggthemes theme_clean
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_point geom_smooth facet_wrap
#' @export
#'
#' @examples
#' x <- sample(1:250, 50)
#' y <- sample(c("A","B", "C"), 50, replace=T)
#' z <- sample(c("LOW","MED","HIGH"), 50, replace=T)
#' w <- sample(1:250, 50)
#' df <- data.frame(x=x, y=y, z=z, w=w)
#'
#' ggplot(df, aes(x=y, y=x, fill=z)) + geom_boxplot() + theme_pub()
#' ggplot(df, aes(x=x, y=w, color=y)) + geom_point() + geom_smooth(method="lm") + theme_pub()
theme_pub <- function(){
  ggthemes::theme_clean() %+replace%
    ggplot2::theme(
      #change background colors
      plot.background = element_rect(fill = "white", colour = "white"),
      legend.background = element_rect(fill = "white", colour = "white"),

      #change text
      axis.title = element_text(face="bold", size=20),
      legend.title = element_text(face="bold", size=20),
      axis.text = element_text(size=16),
      legend.text = element_text(size=16),
      strip.text = element_text(size=20))
}
