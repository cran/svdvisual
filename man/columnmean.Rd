\name{columnmean}
\alias{columnmean}
\title{
Column mean matrix
}
\description{
This function returns the column mean matrix of an input matrix
}
\usage{
  columnmean(x)
}
\arguments{
  \item{x}{
The input matrix. Make sure this is a matrix object
}
}
\details{
This function calculates the mean of each column within the input matrix. If the matrix has dimension m x n, this function will return a m x n matrix, where each column has the same result, the average of the n elements in the corresponding column in the input matrix.
}
\value{
The column mean matrix
}
\author{
Yao Wang (wang1150@purdue.edu);
Lingsong Zhang (lingsong@purdue.edu)
}
\seealso{
See Also in \code{\link{svd}}, \code{\link{apply}}, \code{\link{rowmean}}, \code{\link{doublemean}}, \code{\link{overallmean}}.
}
\examples{
#generate a random matrix
x<-matrix(rnorm(100), nrow=20);

#calculate the column mean matrix
y<-columnmean(x);
y
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{svd}
\keyword{mean}
\keyword{matrix}% __ONLY ONE__ keyword per line
