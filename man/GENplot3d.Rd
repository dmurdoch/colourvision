\name{GENplot3d}
\alias{GENplot3d}
\title{
Generic model colour space 3D plot
}
\description{
Plots models based on the GENmodel( ) function for tetrachromatic animals.
}
\usage{
GENplot3d(model, col.names=c("X1","X2","X3"),
          vectors=TRUE, vnames=TRUE, vsize="auto",
          xlab="x", ylab="y", zlab="z",
          xlim="auto", ylim="auto", zlim="auto", asp=1, ...)
}
\arguments{
  \item{model}{Output of a colour vision model.}
  \item{col.names}{Column names to be plotted.}
  \item{vectors}{Whether vectors representing direction of photoreceptor outputs should be plotted.}
  \item{vnames}{Whether vector names should be plotted.}
  \item{vsize}{Length of vectors to be plotted. Default calculates length automatically.}
  \item{xlab}{see \code{\link{par3d}} function.}
  \item{ylab}{see \code{\link{par3d}} function.}
  \item{zlab}{see \code{\link{par3d}} function.}
  \item{xlim}{x-axis range. Default calculates range automatically. See \code{\link{par3d}} function.}
  \item{ylim}{y-axis range. Default calculates range automatically. See \code{\link{par3d}} function.}
  \item{zlim}{z-axis range. Default calculates range automatically. See \code{\link{par3d}} function.}
  \item{asp}{see \code{\link{plot3d}} function.}
  \item{...}{Other arguments passed to \code{\link{plot3d}} function.}
}
\author{
Felipe M. Gawryszewski \email{f.gawry@gmail.com}
}
\seealso{
\code{\link{CTTKhexagon}}, \code{\link{CTTKhexagon3D}}, \code{\link{EMtriangle}}, \code{\link{EMtetrahedron}}, \code{\link{RNLplot}}, \code{\link{RNLplot3d}}, \code{\link{GENplot}}, \code{\link{plot.colourvision}}, \code{\link{plot3d.colourvision}}
}