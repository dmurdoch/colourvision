\name{EMtetrahedron}
\alias{EMtetrahedron}
\title{
Endler and Mielke (2005) tetrahedron colour space}
\description{
Plots Endler and Mielke (2005) tetrahedron colour space for tetrachromatic animals.}
\usage{
EMtetrahedron(x, y, z, s.col = "red", f.col = "black",
              vnames = c("u","s","m","l"), type = "p",
              radius = 0.01, add = F, xlab = "",
              ylab = "", zlab = "", box = F, axes = F,
              ylim = c(-0.75, 0.75), xlim = c(-0.75, 0.75),
              zlim = c(-0.75, 0.75), aspect = T,
              mar = c(1, 1, 1, 1), vectors=FALSE, ...)
}
\arguments{
  \item{x}{x coordinate of points to be plotted}
  \item{y}{y coordinate of points to be plotted}
  \item{z}{z coordinate of points to be plotted}
  \item{s.col}{Colour to be used for plotted items. See \code{\link{plot3d}}.}
  \item{f.col}{Colour of tetrahedron lines. See \code{\link{plot3d}}.}
  \item{vnames}{Vector names.}
  \item{type}{see \code{\link{plot3d}} function.}
  \item{radius}{see \code{\link{plot3d}} function.}
  \item{add}{see \code{\link{plot3d}} function.}
  \item{xlab}{see \code{\link{plot3d}} function.}
  \item{ylab}{see \code{\link{plot3d}} function.}
  \item{zlab}{see \code{\link{plot3d}} function.}
  \item{box}{see \code{\link{plot3d}} function.}
  \item{axes}{see \code{\link{plot3d}} function.}
  \item{ylim}{see \code{\link{plot3d}} function.}
  \item{xlim}{see \code{\link{plot3d}} function.}
  \item{zlim}{see \code{\link{plot3d}} function.}
  \item{aspect}{see \code{\link{plot3d}} function.}
  \item{mar}{see \code{\link{plot3d}} function.}
  \item{vectors}{Whether vectors representing direction of photoreceptor outputs should be plotted.}
  \item{...}{Other arguments passed to function \code{\link{plot3d}}.}
}
\references{
Endler, J. A., and P. Mielke. 2005. Comparing entire colour patterns as birds see them. Biol J Linn Soc 86:405-431.
}
\author{
Felipe M. Gawryszewski \email{f.gawry@gmail.com}
}
\seealso{
\code{\link{EMtriangle}}, \code{\link{EMmodel}}
}
\examples{
\dontrun{
library(rgl)
EMtetrahedron(x=0,y=0,z=0)
}
}