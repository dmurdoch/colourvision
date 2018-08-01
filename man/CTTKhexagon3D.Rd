\name{CTTKhexagon3D}
\alias{CTTKhexagon3D}
\title{
Chittka (1992) colour space for tetrachromatic animals.
}
\description{
Plots a hexagonal trapezohedron representing Chittka (1992) colour space for tetrachromatic animals (Thery and Casas, 2002).}
\usage{
CTTKhexagon3D(x, y, z, s.col = "red", f.col = "black",
              vnames = c("E1","E2","E3","E4"), type = "p",
              radius = 0.01, add = F, xlab = "", ylab = "", zlab = "",
              box = F, axes = F, ylim = c(-1, 1), xlim = c(-1, 1),
              zlim = c(-1,1), aspect = T, vectors=F, ...)
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
  \item{vectors}{Whether vectors representing direction of photoreceptor outputs should be plotted.}
  \item{...}{Other arguments passed to function \code{\link{plot3d}}.}
}
\references{
Chittka, L. 1992. The colour hexagon: a chromaticity diagram based on photoreceptor excitations as a generalized representation of colour opponency. J Comp Physiol A 170, 533-543.

Thery, M., and J. Casas. 2002. Predator and prey views of spider camouflage. Nature 415, 133-133.
}
\author{
Felipe M. Gawryszewski \email{f.gawry@gmail.com}
}
\seealso{
\code{\link{CTTKmodel}}, \code{\link{CTTKhexagon}}
}
\examples{
\dontrun{
library(rgl)
CTTKhexagon3D(x=0,y=0,z=0)
}
}