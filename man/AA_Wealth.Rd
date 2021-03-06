% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AAUtilities.r
\name{AA_Wealth}
\alias{AA_Wealth}
\title{Future Wealth}
\usage{
AA_Wealth(p, r, sd, t, amt = 1)
}
\arguments{
\item{p}{Probability of doing worse than result}

\item{r}{Mean return}

\item{sd}{Standard deviation of returns}

\item{t}{Period in years}

\item{amt}{Iniitial investment. Default = 1}
}
\value{
value Future wealth
}
\description{
Calculates the wealth for which there is only a p chance of doing worse than
over t years given a mean of r and a standard deviation of sd. Assumes a
lognormal distribution.  Specify inputs in percentage formats (e.g. 8%=8, not
0.08).
}
\examples{
AA_Wealth(.5,8,12,10)

}

