svd3dplot <-
function(data, ncomp=3, isurface=T, iimage=F, xlab='Column', ylab='Row', zlab='', ...){
# this is the svd three dimensional plot


data.rank<-matrixrank(data);
if (ncomp>=data.rank){
	print('The specified rank is larger than the rank of the input matrix');
	print('The parameter will be reset as the rank of the input matrix -1');
	ncomp=data.rank-1;
}
data.dim<-dim(data);
nrow=data.dim[1];
ncol=data.dim[2];
ncell=nrow*ncol;

data.svd<-svd(data, nu=ncomp, nv=ncomp);
umat=data.svd$u;
vmat=data.svd$v;
svec=data.svd$d[1:ncomp];

app=umat%*%diag(svec)%*%t(vmat);
res=data-app;
rowmat=seq(1:nrow)%*% t(rep(1, ncol));
colmat=rep(1, nrow) %*% t(seq(1, ncol));

#generate ploting matrix
plotmat<-data.frame(matvalue=as.vector(data), rows=as.vector(rowmat), columns=as.vector(colmat), label=rep("Original data", ncell));
for (i in 1:ncomp){
	tempmat=svec[i]*(umat[, i]%*% t(vmat[, i]));
	plotmat<-rbind(plotmat, data.frame(matvalue=as.vector(tempmat), rows=as.vector(rowmat), columns=as.vector(colmat), label=rep(paste("SVD", i), ncell)));
}
	plotmat<-rbind(plotmat, data.frame(matvalue=as.vector(app), rows=as.vector(rowmat), columns=as.vector(colmat), label=rep("Approximation", ncell)));
	plotmat<-rbind(plotmat, data.frame(matvalue=as.vector(res), rows=as.vector(rowmat), columns=as.vector(colmat), label=rep("Residual", ncell)));

#require('lattice')

if (isurface==iimage){
	print('Can not simultaneously generate image plots and surface plots');
	print('Will only show the image plot');
	isurface=F;
	iimage=T;
}

nplot=ncomp+3;
nplotrow=ceiling(nplot/3);
localcondindex=c(1, nplot-1, nplot, 2:(nplot-2))

if (isurface){
	print(wireframe(matvalue~columns+rows|label, data=plotmat, xlab=xlab, ylab=ylab, zlab=zlab, main='SVD surface plot',  index.cond=list(localcondindex), layout=c(3, nplotrow), ...));
}

if (iimage){
	print(levelplot(matvalue~columns+rows|label, data=plotmat, xlab=xlab, ylab=ylab, zlab=zlab, main='SVD image plot', index.cond=list(localcondindex), layout=c(3, nplotrow), ...));
}
}
