% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ssn-functions.R
\name{ds_select}
\alias{ds_select}
\title{Perform double selection using a p-value limit to determine which variable to include}
\usage{
ds_select(focal = "DNBR", response = "DOC_mgL", p_limit = 0.1, df, covar)
}
\arguments{
\item{focal}{a character of the column name of the variable you're interested in}

\item{response}{a character of the column name of the response variable}

\item{p_limit}{variable with p-value below this limit will be included in the final model}

\item{df}{the dataframe with the data to fit in it}

\item{covar}{a vector of column names of the variables to include as potential covariates}
}
\value{
a formula for the final model
}
\description{
Perform double selection using a p-value limit to determine which variable to include
}
