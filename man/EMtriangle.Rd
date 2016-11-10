\name{EMtriangle}
\alias{EMtriangle}
\title{
Endler and Mielke triangle colour space
}
\description{
Plots a triangle colour space for trichromatic Endler and Mielke (2005) colour vision models.
}
\usage{
EMtriangle(ylim = c(-0.8, 0.8), xlim = c(-0.8, 0.8), vnames=TRUE)
}
\arguments{
  \item{xlim, ylim}{x and y limits.}
  \item{vnames}{
Whether vector names should be plotted.}
}
\details{
The original model is available for tetrachromatic animals only. Doris Gomez derived a trichromatic version which is available in softaware AVICOL (Gomez, 2006) and implemented here.
}
\references{
Endler, J. A., and P. Mielke. 2005. Comparing entire colour patterns as birds see them. Biol J Linn Soc 86:405-431.

Gomez, D. 2006. AVICOL, a program to analyse spectrometric data. Last update october 2011. Free executable available at:\cr
\url{http://sites.google.com/site/avicolprogram/}  or from the author at \email{dodogomez@yahoo.fr}
}
\author{
Felipe M. Gawryszewski \email{f.gawry@gmail.com}
}
\seealso{
\code{\link{EMtetrahedron}},  \code{\link{EMmodel}}
}
\examples{
EMtriangle()
points(x=0,y=0, pch=16, col="red")
}