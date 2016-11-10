\name{EMmodel}
\alias{EMmodel}
\title{
Endler and Mielke colour vision model
}
\description{
Endler and Mielke trichromatic and tetrachromatic colour vision models.
}
\usage{
EMmodel(photo = c("tri", "tetra"), R, I, Rb, C,
        interpolate=TRUE, nm=seq(300,700,1))
}
\arguments{
  \item{photo}{
Trichromatic (\code{"tri"}) or tetrachromatic (\code{"tetra"}) colour vision model.
}
  \item{R}{
Reflectance of observed objects. A data frame with first column corresponding to wavelength values and following columns with reflectance values. \code{R} must be in the same scale as \code{Rb} (percentage or proportion).
}
  \item{I}{
Irradiance spectrum. A data frame with two columns only: first column corresponding to wavelength values and second column with irradiance values. Irradiance values must be in quantum flux units.
}
  \item{Rb}{
Background reflectance. A data frame with two columns only: first column corresponding to wavelength values and second column with reflectance values. \code{Rb} must be in the same scale as \code{R} (percentage or proportion).
}
  \item{C}{
Photoreceptor sensitivity curves, from lowest to longest lambda-max. A data frame: first column corresponding to wavelength values and following columns with photoreceptor sensitivity values (see function \code{\link{photor}}). 
}
  \item{interpolate}{
Whether data files should be interpolated before further calculations. See \code{\link{approx}}.}
  \item{nm}{
A sequence of numeric values specifying where interpolation is to take place. See \code{\link{approx}}.}

}
\details{
The original model is available for tetrachromatic animals only. Doris Gomez derived a trichromatic version which is available in software AVICOL (Gomez, 2006) and implemented here. First, relative quantum catches are log-transformed:
\deqn{f_{i} = \ln{q_i}}{fi = ln(qi)}
where \eqn{q_i}{qi} is the relative quantum catch of photoreceptor type i, given by \code{\link{Qr}}. The model uses only relative output values, so that photoreceptor outputs (u,s,m for thrichromats and u,s,m,l for tetrachromats) are given by:
\deqn{\{u,s,m,l\} = \frac{f_i}{\sum_{i=1}^n f_i}}{{u,s,m,l} = fi'/sum(f1+f2+f3...)}
Then, for thrichromatic vision, coordinates in the colour space are found by (Gomez 2006):
\deqn{x = \frac{\sqrt{3}}{3}(m-u)}{x = (sqrt(3)/3)*(m-u)}
\deqn{y = \frac{2}{3}(s-\frac{1}{2}(u+m))}{y = (2/3)*(s-0.5*(u+m))}
And for tetrachromatic vision (Endler and Mielke 2005):
\deqn{x = \sqrt{\frac{3}{2}}(\frac{1-2s-m-u}{2})}{x = sqrt(3/2)*((1-2s-m-u)/2)}
\deqn{y = \frac{-1+3m+u}{2\sqrt{2}}}{y = (-1+3m+u)/2*sqrt(2)}
\deqn{z = u-\frac{1}{4}}{z = u-0.25}
}

\value{
for trichromats:
\item{f1,f2,f3}{Photoreceptor quantum-catches values after the von Kries transformation (see function \code{\link{Qr}} and the log-transformation.}
\item{u,s,m}{Photoreceptor output values.}
\item{x,y}{Coordinates in the animal colour space. See function \code{\link{EMtriangle}}}
\item{deltaSo}{Euclidian distance to the origin of the colour space. It represents the conspicouness of the stimulus (\code{R}) in relation to the background (\code{Rb}).}
for tetrachromats:
\item{f1,f2,f3,f4}{Photoreceptor quantum-catches values after the von Kries transformation (see function \code{\link{Qr}}.}
\item{u,s,m,l}{Photoreceptor output values.}
\item{x,y,z}{Coordinates in the animal colour space. See function \code{\link{EMtetrahedron}}}
\item{deltaSo}{Euclidian distance to the origin of the colour space. It represents the conspicouness of the stimulus (\code{R}) in relation to the background (\code{Rb}).}
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
\code{\link{EMtriangle}}, \code{\link{EMtetrahedron}}, \code{\link{photor}}, \code{\link{CTTKmodel}}, \code{\link{RNLmodel}}
}
\examples{
##Photoreceptor sensitivity curves
##with lambda max at 350nm, 450nm and 550nm:
C<-photor(lambda.max=c(350,450,550))

##Gray background
##with 7 percent reflectance from 300 to 700nm:
Rb <- data.frame(300:700, rep(7, length(300:700)))

## Read CIE D65 standard illuminant
data("D65")

##Reflectance data
## with a sigmoid spectrum and midpoint at 500nm and 550 nm
R1<-logistic(x=seq(300,700,1), x0=500, L=50, k=0.04)
R2<-logistic(x=seq(300,700,1), x0=550, L=50, k=0.04)
R<-cbind(R1, R2[,2])
R[,2]<-R[,2]+10
R[,3]<-R[,3]+10

## Run model 
EMmodel(photo="tri",
       R=R, I=D65, Rb=Rb, C=C)
}