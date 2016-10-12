#' @title Inhomogeneous anisotropic spatio-temporal \eqn{K}-function
#' @description This function compute an estimation of the inhomogeneous anisotropic spatio-temporal \eqn{K}-function
#' @param xyt Coordinates and times \eqn{(x,y,t)} of the point pattern.
#' @param s.region Two-column matrix specifying polygonal region containing all data locations. If s.region is missing, the Ripley-Rasson estimate of the spatial domain is considered.
#' @param t.region Vector containing the minimum and maximum values of the time interval. If t.region is missing, the range of \eqn{xyt[,3]} is considered.
#' @param lambda Vector of values of the spatio-temporal intensity function evaluated at the points \eqn{(x,y,t)} in \eqn{W x T}. If lambda is missing, the estimate of the anisotropic spatio-temporal \eqn{K}-function is computed as for the homogeneous case, i.e. considering \eqn{n/|W x T|} as an estimate of the spatio-temporal intensity.
#' @param ds Vector of distances \eqn{u} at which \eqn{\hat{K}_{\phi}(r,t)} is computed.
#' @param dt Vector of times \eqn{v} at which \eqn{\hat{K}_{\phi}(r,t)} is computed.
#' @param ang Vector of angles in radians at which \eqn{\hat{K}_{\phi}(r,t)} is computed.  If ang is missing, the function HASTKfunct is evaluated in the vector common angles.
#' @param correction A character vector specifying the edge correction(s) to be applied among "border", "modified.border", "translate" and "none". The default is "border".
#' @return A list containing:
#' \itemize{
#' \item \code{astkl}: Array containing \code{nds} by \code{ndt} matrices whose elements are \eqn{\hat{K}_{\phi}(r,t)} evaluated in each ang-vector.
#' \item \code{ds}: If \code{dist} is missing, a vector of distances \code{u} at which \eqn{\hat{K}_{\phi}(u,v)} is computed.
#' \item \code{dt}: If \code{times} is missing, a vector of distances \code{v} at which \eqn{\hat{K}_{\phi}(u,v)} is computed.
#' \item \code{lambda}: Value of the estimation of \eqn{\hat{\rho}^{2}}.
#' \item \code{ang}: If \code{ang} is missing, the function si evalueten on the vector \eqn{(pi/6,pi/4,pi/3,pi/2,2pi/3,3pi/4,5pi/6,pi)}.
#' }
#' @author Francisco J. Rodriguez-Cortes <cortesf@@uji.es> \url{https://fjrodriguezcortes.wordpress.com}
#' @references Gabriel, E., Rowlingson, B., Diggle P J. (2013). \code{stpp}: an R package for plotting, simulating and analyzing Spatio-Temporal Point Patterns. Journal of Statistical Software 53, 1-29.
#' @references Illian, J. B., Penttinen, A., Stoyan, H. and Stoyan, D. (2008). Statistical Analysis and Modelling of Spatial Point Patterns. John Wiley and Sons, London.
#' @references Gonzalez, J. A., Rodriguez-Cortes, F. J., Cronie, O., Mateu, J. (2016). Spatio-temporal point process statistics: a review. Spatial Statistics. Accepted.
#' @references Ohser, J. and D. Stoyan (1981). On the second-order and orientation analysis of planar stationary point processes. Biometrical Journal 23, 523-533.
#' @examples
#' ## Not run:
#' #################
#'
#' # Realisations of the homogeneous spatio-temporal Poisson processes
#' stp <- rpp(100)$xyt
#' # Generated spatio-temporal point pattern
#' plot(stp)
#' 
#' # Estimation of the anisotropic homogeneous spatio-temporal K-functions
#' out <- HASTKfunct(stp)
#' 
#' z1 <- out$astkf[,,1]
#' z2 <- out$astkf[,,2]
#' z3 <- out$astkf[,,3]
#' z4 <- out$astkf[,,4]
#' 
#' # Plot
#' par(mfrow=c(2,2))
#' persp(out$ds,out$dt,z1,theta=-45,phi=30,zlim=range(z1,na.rm=TRUE),ticktype="detailed",xlab="\n r = distance",ylab="\n t = time",zlab="",nticks=6,cex.axis=1.3,cex.lab=1.7)
#' persp(out$ds,out$dt,z2,theta=-45,phi=30,zlim=range(z2,na.rm=TRUE),ticktype="detailed",xlab="\n r = distance",ylab="\n t = time",zlab="",nticks=6,cex.axis=1.3,cex.lab=1.7)
#' persp(out$ds,out$dt,z3,theta=-45,phi=30,zlim=range(z3,na.rm=TRUE),ticktype="detailed",xlab="\n r = distance",ylab="\n t = time",zlab="",nticks=6,cex.axis=1.3,cex.lab=1.7)
#' persp(out$ds,out$dt,z4,theta=-45,phi=30,zlim=range(z4,na.rm=TRUE),ticktype="detailed",xlab="\n r = distance",ylab="\n t = time",zlab="",nticks=6,cex.axis=1.3,cex.lab=1.7)
#' 
#' ## End(Not run)
#' #################
HASTKfunct <- function(xyt,s.region,t.region,lambda,ds,dt,ang,correction="border") {
  
  correc=c("none","border","modified.border","translate")
  id <- match(correction, correc, nomatch = NA)
  if (any(nbg <- is.na(id))) {
    mess <- paste("unrecognised correction method:", paste(dQuote(correction[nbg]),collapse = ", "))
    stop(mess, call. = FALSE)
  }
  id=unique(id)	
  correc2=rep(0,4)
  correc2[id]=1	
  
  if (missing(s.region)) {
    x <- xyt[,1]
    y <- xyt[,2]
    W <- ripras(x,y)
    poly <- W$bdry
    X <- poly[[1]]$x
    Y <- poly[[1]]$y
    s.region <- cbind(X,Y)
    }
  
  bsw <- owin(poly=list(x=s.region[,1],y=s.region[,2]))
  
  if (missing(ds)){
    rect <- as.rectangle(bsw)
    maxd <- min(diff(rect$xrange),diff(rect$yrange))/4
    ds <- seq(0, maxd, len=20)
    ds <- sort(ds)
    }
  if(ds[1]==0) {ds=ds[-1]}
  
  if (missing(t.region)){
    xr = range(xyt[,3],na.rm=TRUE)
    xw = diff(xr)
    t.region <- c(xr[1]-0.01*xw,xr[2]+0.01*xw)
    }
  
  bsupt <- max(t.region)
  binft <- min(t.region)

  if (missing(dt)){
    maxt <- (bsupt-binft)/4
    dt <- seq(0, maxt, len=20)
    dt <- sort(dt)
  }
  if(dt[1]==0) {dt=dt[-1]}
  
  if (missing(ang)) {ang <-c(pi/6,pi/4,pi/3,pi/2,(2*pi)/3,(3*pi)/4,(5*pi)/6,pi)}

  pts <- xyt[,1:2]
  xytimes <- xyt[,3]
  ptsx <- pts[, 1]
  ptsy <- pts[, 2]
  ptst <- xytimes
  npt <- length(ptsx)
  nds <- length(ds)
  area <- area(bsw)
  ndt <- length(dt)
  lgt <- (bsupt - binft)
  ang <- sort(ang)
  vol <- area*lgt
  
  wbi = array(0,dim=c(npt,nds,ndt))
  wbimod = array(0,dim=c(npt,nds,ndt))
  wt = array(0,dim=c(npt,npt))
  
  pppxy = ppp(x=ptsx,y=ptsy,window=bsw)

  if(missing(lambda))
  {
    misl <- 1
    lambda <- rep(npt/vol,npt)
  }
  else misl <- 0
  if (length(lambda)==1) lambda <- rep(lambda,npt)
 
  #  correction=="border" and "modified border"
  
  if(any(correction=="border")|any(correction=="modified.border"))
  {
    bi=bdist.points(pppxy)
    bj=.bdist.times(ptst,t.region)
    
    for(i in 1:nds) 
    { 
      for(j in 1:ndt)
      {
        wbi[,i,j] = (bi>ds[i])*(bj>dt[j])/sum((bi>ds[i])*(bj>dt[j])/lambda)
        wbimod[,i,j] = (bi>ds[i])*(bj>dt[j])/(eroded.areas(bsw,ds[i])*.eroded.areat(t.region,dt[j]))
      } }
    wbi[is.na(wbi)]=0
  }
  
  # correction=="translate"
  
  if(any(correction=="translate"))
  {
    wtt = .overlap.tint(xytimes,t.region)
    wts = edge.Trans(pppxy)
    wt = wtt*wts
    wt=1/wt
  }

astkl <- sapply(ang, function(ang) astk.ang(ptsx=ptsx,ptsy=ptsy,ptst=ptst,npt=npt,lambda=lambda,ag=ang,ds=dt,nds=ndt,dt=dt,ndt=ndt,vol=vol,wbi=wbi,wbimod=wbimod,wt=wt,correc2=correc2)$astkf,simplify="array")
  
invisible(return(list(astkf=astkl,ds=ds,dt=dt,lambda=lambda,ang=ang,correction=correction)))
}

  # sub-functions

astk.ang <- function(ptsx,ptsy,ptst,npt,lambda,ag,ds,nds,dt,ndt,vol,wbi,wbimod,wt,correc2){

astkf <- array(0,dim=c(nds, ndt))
 
storage.mode(astkf) <- "double"
 
 astk <- .Fortran("astkfunct",as.double(ptsx),as.double(ptsy),as.double(ptst),as.integer(npt),as.double(lambda), 
                 as.double(ag),as.double(ds),as.integer(nds),as.double(dt),as.integer(ndt),as.double(vol),
                 as.double(wbi),as.double(wbimod),as.double(wt),as.integer(correc2),(astkf))
  
 astkf <- astk[[16]]
  
 invisible(return(list(astkf=astkf)))
 }

.overlap.tint=function(times,t.region)
{
  if (missing(t.region)) t.region=range(times)
  ntimes=length(times)
  
  a = diff(range(t.region))
  
  wt=matrix(a,ncol=ntimes,nrow=ntimes)
  for(i in 1:ntimes)
  { for(j in 1:ntimes)
  {
    if (i!=j)
    {
      b = a-abs(times[i]-times[j])
      wt[i,j]=a/b
    }}}
  invisible(return(wt))
}

.bdist.times=function(times, t.region)
{
  if (missing(t.region)) t.region=range(times)
  ntimes=length(times)
  a=min(t.region)
  b=max(t.region)
  
  bj=NULL
  for(j in 1:ntimes)
    bj=c(bj,min(c(abs(times[j]-a),abs(times[j]-b))))
  
  invisible(return(bj))
}

.eroded.areat=function(t.region,dist)
{
  a = diff(range(t.region))
  b = a-dist
  invisible(return(b))
}