\name{matrixrank}
\alias{matrixrank}
\title{
Rank of a Matrix
}
\description{
This function returns the rank of an input matrix
}
\usage{
  matrixrank(x, tolerance = 1e-08)
}
\arguments{
  \item{x}{
The input matrix. Make sure this is a matrix object
}
  \item{tolerance}{
The tolerance of the numerical zeros. All singular values whose absolute value is smaller than the tolerance will be treated as zero.
}
}
\details{
This function identifies the singular values, and count the number of such values significant away from 0. The result is the rank of the input matrix.
}
\value{
The rank of the input matrix.
}
\author{
Lingsong Zhang (lingsong@purdue.edu)
}
\seealso{
See Also in \code{\link{svd}}.
}
\examples{
#generate a random matrix
x<-matrix(rnorm(100), nrow=20);

#calculate the rank of this matrix
y<-matrixrank(x);
y
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{svd}
\keyword{matrix}% __ONLY ONE__ keyword per line
