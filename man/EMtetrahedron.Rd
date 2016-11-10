\name{EMtetrahedron}
\alias{EMtetrahedron}
\title{
Endler and Mielke tetrahedron colour space}
\description{
Plots Endler and Mielke (2005) tetrahedron colour space for tetrachromatic animals.}
\usage{
EMtetrahedron(x, y, z, type="p", s.col="red",
              radius=0.01, f.col="black", vnames=TRUE)
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
Colour of tetrahedron lines. See \code{plot3d}.}
  \item{vnames}{
Whether vector names should be plotted.}
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
EMtetrahedron(x=0,y=0,z=0)
}