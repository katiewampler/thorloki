% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{str_split2}
\alias{str_split2}
\title{Split up a string into pieces}
\usage{
str_split2(string, pattern, n = Inf, piece = NA)
}
\arguments{
\item{string}{Input vector. Either a character vector, or something coercible to one.}

\item{pattern}{Pattern to look for.}

\item{n}{Maximum number of pieces to return. Default (Inf) uses all possible split positions.}

\item{piece}{A vector of numbers, if it's not NA it will return that piece of each vector}
}
\value{
if the string is length 1 or 1 piece from each string is requested, it will
return a vector. If multiple pieces of multiple strings is requested
it will create a character matrix.
}
\description{
A wrapper for str_split from stringr that will return a vector if
the string has length 1
}
\examples{
string <- "hello/world/its/katie"
str_split2(string, "/")
str_split2(string, "/", piece=c(2:4))
string2 <- c("loki_is_here", "thor_is_here_too")
str_split2(string2,"_")
str_split2(string2,"_", piece=c(1,3))

}
