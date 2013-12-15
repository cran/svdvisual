\name{rowmean}
\alias{rowmean}
\title{
Row mean matrix
}
\description{
This function returns the row mean matrix of an input matrix
}
\usage{
  rowmean(x)
}
\arguments{
  \item{x}{
The input matrix. Make sure this is a matrix object
}
}
\details{
This function calculates the mean of each row within the input matrix. If the matrix has dimension m x n, this function will return a m x n matrix, where each row has the same result, and the average of the n elements in the corresponding row in the input matrix.
}
\value{
The row mean matrix
}
\author{
Yao Wang (wang1150@purdue.edu);
Lingsong Zhang (lingsong@purdue.edu)
}
\seealso{
See Also in \code{\link{svd}}, \code{\link{apply}}, \code{\link{columnmean}}, \code{\link{doublemean}}, \code{\link{overallmean}}.
}
\examples{
#generate a random matrix
x<-matrix(rnorm(100), nrow=20);

#calculate the row mean matrix
y<-rowmean(x);
y
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{svd}
\keyword{mean}
\keyword{matrix}% __ONLY ONE__ keyword per line
