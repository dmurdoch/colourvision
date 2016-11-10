\name{CTTKhexagon3D}
\alias{CTTKhexagon3D}
\title{
Chittka (1992) colour space for tetrachromatic animals.
}
\description{
Plots a hexagonal trapezohedron representing Chittka (1992) colour space for tetrachromatic animals.}
\usage{
CTTKhexagon3D(x, y, z, type="p", s.col="red",
              radius=0.01, f.col="grey", vnames=TRUE)
}
\arguments{
  \item{x}{
x coordinate of points to be plotted
}
  \item{y}{
y coordinate of points to be plotted
}
  \item{z}{
z coordinate of points to be plotted
}
  \item{type}{
Plot points "p" or spheres "s". See \code{plot3d}.}
  \item{s.col}{
Color to be used for plotted items. See \code{plot3d}.}
  \item{radius}{
Radius of spheres. See \code{plot3d}.}
  \item{f.col}{
Colour of polygon lines. See \code{plot3d}.}
  \item{vnames}{
Whether vector names should be plotted.}
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
CTTKhexagon3D(x=0,y=0,z=0)}