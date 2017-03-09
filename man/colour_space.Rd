\name{colour_space}
\alias{colour_space}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
N-dimensional colour spaces
}
\description{
Generates a colour space based on any number of photorecptor types and finds a colour loci for a given photoreceptor output.}
\usage{
colour_space(n, length, q)
}
\arguments{
  \item{n}{Number of photoreceptor types. Function accepts any number of photoreceptor types >=3. For instance, trichromatic: \code{n=3}; tetrachromatic: \code{n=4}, etc.
}
  \item{length}{Vector length used to construct the colour space. Typically \code{length=1}}
  \item{q}{Photoreceptor output values}
}
\details{Based on Pike (2012), but with a fixed vector length, instead of fixed edge length. This function is used internally in colour vision models.
}
\value{
A list with the following dimentions:
\item{coordinates}{Colour loci of photoreceptor output values \code{q} in the colour space}
\item{vector_matrix}{Coordinates of colour space vertices}
}
\references{
Pike, T.W. 2012. Generalised chromaticity diagrams for animals with n-chromatic colour vision. Journal of Insect Behavior 255: 277-286.}
\author{
Felipe M. Gawryszewski \email{f.gawry@gmail.com}
}
\seealso{
\code{\link{Qr}}, \code{\link{CTTKmodel}}, \code{\link{EMmodel}}
}

\examples{
#A trichromatic colour space based on Endler and Mielke (2005)
tri<-colour_space(n=3, length=0.75, q=c(0.5,0.2,0.3))

#showing:
#(1) Limits of the colour space (triangle)
plot(-1000, ylim=c(-1,1), xlim=c(-1,1), asp=1, ylab="X2", xlab="X1")
polygon(tri$vector_matrix, lty=2)

#(2) Vectors (length=0.75) used to build the colour space (arrows)
arrows(x0=0,y0=0, x1=tri$vector_matrix[1,1], y1=tri$vector_matrix[1,2], col="red")
arrows(x0=0,y0=0, x1=tri$vector_matrix[2,1], y1=tri$vector_matrix[2,2], col="red")
arrows(x0=0,y0=0, x1=tri$vector_matrix[3,1], y1=tri$vector_matrix[3,2], col="red")

#(3) Colour loci of given photoreceptor outputs
points(x=tri$coordinates[[1]], y=tri$coordinates[[2]], pch=21, col="blue", bg="blue")
}