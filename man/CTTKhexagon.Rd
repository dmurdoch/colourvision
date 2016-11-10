\name{CTTKhexagon}
\alias{CTTKhexagon}
\title{
Chittka (1992) colour hexagon}
\description{
Plots Chittka (1992) colour hexagon for trichromatic animals.
}
\usage{
CTTKhexagon(vnames=TRUE)
}
\arguments{
  \item{vnames}{
Whether vector names should be plotted.
}

}
\references{
Chittka, L. 1992. The colour hexagon: a chromaticity diagram based on photoreceptor excitations as a generalized representation of colour opponency. J Comp Physiol A 170:533-543.
}
\author{
Felipe M. Gawryszewski \email{f.gawry@gmail.com}
}
\seealso{
\code{\link{CTTKmodel}}, \code{\link{CTTKhexagon3D}}
}
\examples{
##Honeybee photoreceptor sensitivity curves
data("bee")

##Grey background:
## with 7 percent reflectance from 300 to 700nm:
Rb <- data.frame(300:700, rep(7, length(300:700)))

## Read CIE D65 standard illuminant already converted to quantum flux:
data("D65")

##Reflectance data
## with a sigmoid spectrum and midpoint at 500nm and 550 nm
R1<-logistic(x=seq(300,700,1), x0=500, L=50, k=0.04)
R2<-logistic(x=seq(300,700,1), x0=550, L=50, k=0.04)
R<-cbind(R1, R2[,2])

## Run colour vision model:
model<-CTTKmodel(photo="tri", R=R, I=D65, Rb=Rb,
C=bee)

##plot data in the colour space
CTTKhexagon()
points(x=model$x, y=model$y)
}