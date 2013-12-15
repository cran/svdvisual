\name{doublemean}
\alias{doublemean}
\title{
Double mean matrix
}
\description{
This function returns the double mean matrix of an input matrix
}
\usage{
  doublemean(x)
}
\arguments{
  \item{x}{
The input matrix. Make sure this is a matrix object
}
}
\details{
This function calculates the mean of each row, each column and the overall mean within the input matrix. If the matrix has dimension m x n, this function will return a m x n matrix, where each cell is the sum of the row mean and the column mean, subtract by the overall mean.

Let X be the \eqn{m \times n} input matrix, \eqn{x_{ij}} is the i, j cell within it. Let D be the resulting matrix of this program, and \eqn{d_{ij}} as the corresponding cell. Let
\deqn{
	r_i =\frac{1}{n}\sum_{j=1}^n x_{ij}, 
	c_j=\frac{1}{m}\sum_{i=1}^m x_{ij},
	o=\frac{1}{mn}\sum_{i=1}^m\sum_{j=1}^n x_{ij}.
	}
Then
\deqn{
	d_{ij}=r_i+c_j-o;
	}
}
\value{
The double mean matrix
}
\author{
Yao Wang (wang1150@purdue.edu);
Lingsong Zhang (lingsong@purdue.edu)
}
\seealso{
See Also in \code{\link{svd}}, \code{\link{apply}}, \code{\link{columnmean}}, \code{\link{rowmean}}, \code{\link{overallmean}}.
}
\examples{
#generate a random matrix
x<-matrix(rnorm(100), nrow=20);

#calculate the row mean matrix
y<-doublemean(x);
y
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{svd}
\keyword{mean}
\keyword{matrix}% __ONLY ONE__ keyword per line
