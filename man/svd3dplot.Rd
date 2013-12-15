\name{svd3dplot}
\alias{svd3dplot}
\title{
The SVD three dimensional plots: surface plot and/or image plot for SVD decomposition
}
\description{
This function provides surface or image plot for singular value decomposition method. The plot includes several subplots for the following components: original data, cumulative approximation matrix, residual matrix, and several rank 1 SVD components. 
}
\usage{
svd3dplot(data, ncomp = 3, irow=F, icol=F, isurface = T, iimage = F, 
      xlab = "Column", ylab = "Row", zlab = "", ...)
}
\arguments{
  \item{data}{
The input data matrix
}
  \item{ncomp}{
The number of components to calculate. The default value is 3. If the specified number is larger than the rank of the matrix, this option will automatically be set to the rank of the matrix minus one. When either irow and icol is specified as \code{TRUE}, the output will include a mean matrix, and ncomp-1 SVD component of the demeaned matrix. 
}
  \item{irow}{
  	a logical number. If irow and icol both are \code{TRUE}, this program will calculate a \code{\link{doublemean}}. The resulting SVD will be based on the demeaned matrix (i.e., removing the double mean). If only irow is \code{TRUE}, this program will calculate a \code{\link{rowmean}}, and the resulting SVD will be based on the corresponding demeaned matrix (i.e., removing the row mean).
  }
  \item{icol}{
  	a logical number. If irow and icol both are \code{TRUE}, this program will calculate a \code{\link{doublemean}}. The resulting SVD will be based on the demeaned matrix (i.e., removing the double mean). If only icol is \code{TRUE}, this program will calculate a \code{\link{columnmean}}, and the resulting SVD will be based on the corresponding demeaned matrix (i.e., removing the column mean).  	
  }
  \item{isurface}{
Whether the surface plot will be generated. The default value is \code{TRUE}.
}
  \item{iimage}{
Whether the image plot will be generated. The default value is \code{FALSE}.
}
  \item{xlab}{
The \code{xlab} option for the plots. The default value is \code{Column}.
}
  \item{ylab}{
The \code{ylab} option for the plots. The default value is \code{Row}.
}
  \item{zlab}{
The \code{zlab} option for the plots. The default value is empty.
}
  \item{\dots}{
other related ploting options for \code{\link{wireframe}} or \code{\link{levelplot}} in the trellis plot.
}
}
\value{
The code will generate either a surface plot or an image plot for the individual SVD components, the original data, the approximation data and the residual data. 
}
\references{
See detailed explanation of this visualizaiton in \\
Zhang, L., Marron, J. S., Shen, H. and Zhu, Z. (2007), Singular Value Decomposition and Its Visualization, Journal of Computational and Graphical Statistics. 
}
\author{
Lingsong Zhang (lingsong@purdue.edu)
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See Also as \code{\link{svd}}, \code{\link{wireframe}}, and \code{\link{levelplot}}.
}
\examples{
#generate a random sample
#generate a random matrix
x<-matrix(rnorm(100), nrow=20);

#generate a surface plot
svd3dplot(x);

#generate an image plot
svd3dplot(x, iimage=TRUE);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{svd}
\keyword{visualization}% __ONLY ONE__ keyword per line
