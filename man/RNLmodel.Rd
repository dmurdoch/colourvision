\name{RNLmodel}
\alias{RNLmodel}
\title{
Receptor Noise Limited Models
}
\description{
Receptor noise limited colour vision models for dichromatic, trichromatic and tetrachromatic animals (Vorobyev & Osorio 1998; Vorobyev et al. 1998).
}
\usage{
RNLmodel(model=c("linear","log"), photo=c("di","tri","tetra"),
         R, Rb, I, C, noise=FALSE, v, n, e,
         interpolate=TRUE, nm=seq(300,700,1))
}
\arguments{
\item{model}{Use linear (\code{"linear"}; Vorobyev & Osorio, 1998), or log-linear  (\code{"log"}; Vorobyev et al. 1998) RNLmodel version. \code{model="log"} is preferred.

\code{model="linear"} can be used only if comparison is to be made between two very similar colours.
}
  \item{photo}{
Dichromatic (\code{"di"}), trichromatic (\code{"tri"}) or tetrachromatic (\code{"tetra"}) colour vision model.
}
  \item{R}{
Reflectance of observed objects. A data frame with first column corresponding to wavelength values and following columns with reflectance values. \code{R} must be in the same scale as \code{Rb} (percentage or proportion).
}
  \item{Rb}{
Background reflectance. A data frame with two columns only: first column corresponding to wavelength values and second column with reflectance values. \code{Rb} must be in the same scale as \code{R} (percentage or proportion).
}
  \item{I}{
Irradiance spectrum. A data frame with two columns only: first column corresponding to wavelength values and second column with irradiance values. Irradiance values must be in quantum flux units.
}
  \item{C}{
Photoreceptor sensitivity curves, from lowest to longest lambda-max. A data frame: first column corresponding to wavelength values and following columns with photoreceptor sensitivity values (see function \code{\link{photor}}). 
}
  \item{noise}{
Logical. Whether receptor noise is provided (\code{noise = TRUE}) or calculated from photoreceptor relative abundances (\code{noise = FALSE}).
}
  \item{e}{
Receptor noise of each photoreceptor type. It is used when \code{noise = TRUE}
}
  \item{n}{
Relative number of each photoreceptor type in the retina. Usually increases with lambda-max. It is used to calculate \code{e} when \code{noise = FALSE}.
}
  \item{v}{
Noise-to-signal ratio of a single photoreceptor. It is used to calculate \code{e} when \code{noise = FALSE}.
}
  \item{interpolate}{
Whether data files should be interpolated before further calculations. See \code{\link{approx}}.}
  \item{nm}{
A sequence of numeric values specifying where interpolation is to take place. See \code{\link{approx}}.}
}
\value{
for dichromatic animals:
\item{e1,e2}{Noise of photoreceptors channels.}
\item{f1,f2}{photoreceptor outputs.}
\item{x}{coordinate in the colour space.}
\item{deltaSo}{Euclidian distance to the origin of the colour space. It represents the conspicouness of the stimulus (\code{R}) in relation to the background (\code{Rb}).}

for trichromatic animals:
\item{e1,e2,e3}{Noise of photoreceptors channels.}
\item{f1,f2,f3}{photoreceptor outputs.}
\item{x,y}{coordinates in the colour space.}
\item{deltaSo}{Euclidian distance to the origin of the colour space. It represents the conspicouness of the stimulus (\code{R}) in relation to the background (\code{Rb}).}

for tetrachromatic animals:
\item{e1,e2,e3,e4}{Noise of photoreceptors channels.}
\item{f1,f2,f3,f4}{photoreceptor outputs.}
\item{x,y,z}{coordinates in the colour space.}
\item{deltaSo}{Euclidian distance to the origin of the colour space. It represents the conspicouness of the stimulus (\code{R}) in relation to the background.chromaticity distance between reflectance  and the background (\code{Rb}).}
}
\details{
The receptor noise limited model was originally developed to calculate \eqn{\Delta S}{deltaS} between two reflectance curves directly, without finding colour locus coordinates (e.g. \code{x,y}; Vorobyev and Osorio 1998). This function uses later formulae for finding colour loci in a chromaticity diagram (Hempel de Ibarra et al. 2001; Renoult et al. 2015). To compare two reflectance spectra one should calculate the Euclidean distance between them (e.g. \eqn{\Delta S = \sqrt{(x_1-x_2)^2+(y_1-y_2)^2}}{deltaS=sqrt((x1-x2)^2+(y1-y2)^2)}).

In lack of a direct measurement, receptor noise (\eqn{e_i}{ei}) can be estimated by the relative abundance of photoreceptor types in the retina, and a measurement of a single photoreceptor noise-to-signal ratio:
\deqn{e_i=\frac{\nu}{\sqrt{\eta _i}}}{ei=v/sqrt(n)}
where \eqn{\nu}{v} is the noise-to-signal ratio of a single photoreceptor, and \eqn{\eta}{n} is the relative abundance of photoreceptor i in the retina. Alternativelly, noise may be dependent of the intensity, but this possibility is not implement in \code{colourvision} yet. Noise dependent of intensity usually holds for low light conditions only (Vorobyev et al. 1998). 
}
\references{
Hempel de Ibarra, N., M. Giurfa, and M. Vorobyev. 2001. Detection of coloured patterns by honeybees through chromatic and achromatic cues. J Comp Physiol A 187:215-224.

Renoult, J. P., A. Kelber, and H. M. Schaefer. 2015. Colour spaces in ecology and evolutionary biology. Biol Rev Camb Philos Soc, doi: 10.1111/brv.12230

Vorobyev, M., and D. Osorio. 1998. Receptor noise as a determinant of colour thresholds. Proceedings of the Royal Society B 265:351-358.

Vorobyev, M., D. Osorio, A. T. D. Bennett, N. J. Marshall, and I. C. Cuthill. 1998. Tetrachromacy, oil droplets and bird plumage colours. J Comp Physiol A 183:621-633.
}

\author{
Felipe M. Gawryszewski \email{f.gawry@gmail.com}
}
\seealso{
\code{\link{photor}}, \code{\link{CTTKmodel}}, \code{\link{EMmodel}}
}
\examples{
## Photoreceptor sensitivity spectra
##with lambda max at 350nm, 450nm and 550nm:
C<-photor(lambda.max=c(350,450,550))

##Grey background
##with 7 percent reflectance from 300 to 700nm:
Rb <- data.frame(300:700, rep(7, length(300:700)))

## Read CIE D65 standard illuminant:
data("D65")

##Reflectance data
## with a sigmoid spectrum and midpoint at 500nm and 550 nm
R1<-logistic(x=seq(300,700,1), x0=500, L=50, k=0.04)
R2<-logistic(x=seq(300,700,1), x0=550, L=50, k=0.04)
R<-cbind(R1, R2[,2])
R[,2]<-R[,2]+10
R[,3]<-R[,3]+10

## Run model 
RNLmodel(photo="tri", model="log",
       R=R, I=D65, Rb=Rb,
       C=C,
       noise=TRUE, e = c(0.13, 0.06, 0.12))
}
