\name{overallmean}
\alias{overallmean}
\title{
	Global mean matrix
}
\description{
This function returns the global mean matrix of an input matrix
}
\usage{
  overallmean(x)
}
\arguments{
  \item{x}{
The input matrix. Make sure this is a matrix object
}
}
\details{
This function calculates the mean of all cells within the input matrix. If the matrix has dimension m x n, this function will return a m x n matrix, where each cell has the same result, as the average of the mn elements in the corresponding input matrix.
}
\value{
The overall mean matrix
}
\author{
Yao Wang (wang1150@purdue.edu);
Lingsong Zhang (lingsong@purdue.edu)
}
\seealso{
See Also in \code{\link{svd}}, \code{\link{apply}}, \code{\link{columnmean}}, \code{\link{doublemean}}, \code{\link{rowmean}}.
}
\examples{
#generate a random matrix
x<-matrix(rnorm(100), nrow=20);

#calculate the row mean matrix
y<-overallmean(x);
y
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{svd}
\keyword{mean}
\keyword{matrix}% __ONLY ONE__ keyword per line
