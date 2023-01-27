library(raster)
library(dplyr)
library(Matrix)
library(Rcpp)
library(R.matlab)
library(INLA)
library(pracma)
library(plot.matrix)
library(sf)
library(spData)

##################

norm_vec <- function(x) (sqrt(sum(x^2)))

CrossProduct3D <- function(x, y, i=1:3) {
  # Project inputs into 3D, since the cross product only makes sense in 3D.
  To3D <- function(x) head(c(x, rep(0, 3)), 3)
  x <- To3D(x)
  y <- To3D(y)
  
  # Indices should be treated cyclically (i.e., index 4 is "really" index 1, and
  # so on).  Index3D() lets us do that using R's convention of 1-based (rather
  # than 0-based) arrays.
  Index3D <- function(i) (i - 1) %% 3 + 1
  
  # The i'th component of the cross product is:
  # (x[i + 1] * y[i + 2]) - (x[i + 2] * y[i + 1])
  # as long as we treat the indices cyclically.
  return (x[Index3D(i + 1)] * y[Index3D(i + 2)] -
            x[Index3D(i + 2)] * y[Index3D(i + 1)])
}

triangleArea<-function(p1,p2,p3){
  a<-norm_vec(p1-p2)
  b<-norm_vec(p1-p3)
  c<-norm_vec(p2-p3)
  
  s<-(a+b+c)/2
  A <- sqrt(abs(s*(s-a)*(s-b)*(s-c)))
  
  return(A)
}

PrecisionMaker <- function(vLoc,tv,tt=NULL) {
  
  library(Matrix)
  
  #Calculate triangel area
  meshArea<-rep(NA,nrow(tv))
  
  for (i in 1:nrow(tv)) {
    p1<-vLoc[tv[i,1],]
    p2<-vLoc[tv[i,2],]
    p3<-vLoc[tv[i,3],]
    meshArea[i]<-triangleArea(p1,p2,p3)
  }
  
  #Calculation locations of centroids
  
  meshCLoc<-(vLoc[tv[,1],]+vLoc[tv[,2],]+vLoc[tv[,3],])/3
  
  #Calculate locations of centroids of each edge
  
  meshELoc_tmp<-list()
  
  meshELoc_tmp[[1]]<-as.matrix((vLoc[tv[,2],]+vLoc[tv[,3],])/2)
  meshELoc_tmp[[2]]<-as.matrix((vLoc[tv[,1],]+vLoc[tv[,3],])/2)
  meshELoc_tmp[[3]]<-as.matrix((vLoc[tv[,1],]+vLoc[tv[,2],])/2)
  
  meshELoc=matrix(NA,3*nrow(meshELoc_tmp[[1]]),ncol(meshELoc_tmp[[1]]))
  
  meshELoc[3*c(1:nrow(meshELoc_tmp[[1]]))-2,]=as.matrix(meshELoc_tmp[[1]])
  meshELoc[3*c(1:nrow(meshELoc_tmp[[1]]))-1,]=as.matrix(meshELoc_tmp[[2]])
  meshELoc[3*c(1:nrow(meshELoc_tmp[[1]])),]=as.matrix(meshELoc_tmp[[3]])
  
  
  #Calculate vertex to triangles map
  
  triangleName<-rownames(tv)
  
  ## Option 1: might not be ideal if index is wrong
  # max(tv)
  
  ## Option 2
  
  meshVT<-list()
  for (i in sort(unique(as.numeric(as.matrix(tv))))) {
    templist<-list(which(apply(tv, 1, function(x) {i %in% x})==TRUE))
    names(templist)<-i
    meshVT<-c(meshVT,templist)
  }
  
  # Get weighting based on meshVT
  
  Lam<-list()
  for (ii in 1:nrow(vLoc)) {
    lam<-meshArea[meshVT[[ii]]]
    lam<-lam/sum(lam)
    templist<-list(sparseMatrix(i=c(1:length(lam)),j=c(1:length(lam)),x=lam))
    names(templist)<-ii
    Lam<-c(Lam,templist)
  }
  
  # Create new coordinate system
  
  # Seems could only be used for globe coordinate (3D)
  
  Weights<-list()
  
  rVal<-NULL
  rIdx<-NULL
  vIdx<-NULL
  startIdx<-1
  
  for (ii in 1:nrow(vLoc)) {
    v1<-c(0,0,0)
    v2<-c(0,0,0)
    if(abs(vLoc[ii,3])!=1){
      v1[1]<-vLoc[ii,2]*(-1)
      v1[2]<-vLoc[ii,1]
    }
    else{v1[1]<-1}
    v1<-v1/norm_vec(v1)
    v2<-CrossProduct3D(as.numeric(vLoc[ii,]),v1)
    v2<-v2/norm_vec(v2)
    
    # Project coordinates just after this
    vOld<-meshCLoc[meshVT[[ii]],]
    
    for (iii in 1:length(meshVT[[ii]])) {
      vOld[iii,]<-vOld[iii,]-vLoc[ii,]
    }
    
    xNew<-as.matrix(vOld)%*%v1
    yNew<-as.matrix(vOld)%*%v2
    
    # Solve weighted least squares
    # Weights are a matrix (sparse matrix class) with dim t(A)%*%Lam[[ii]]
    
    A<-cbind(rep(1,length(xNew)),xNew,yNew)
    A<-as.matrix(A)
    WW<-solve(t(A)%*%Lam[[ii]]%*%A,t(A)%*%Lam[[ii]])
    
    Weights[[ii]]<-WW[1,]
    
    rVal<-c(rVal,WW[1,])
    rIdx<-c(rIdx,meshVT[[ii]])
    vIdx<-c(vIdx,ii*rep(1,length(meshVT[[ii]])))
    startIdx<-c(startIdx,startIdx[length(startIdx)]+length(meshVT[[ii]]))
  }
  
  AtoV<-sparseMatrix(i=vIdx,j=rIdx,x=rVal)
  
  output<-list(meshArea,meshCLoc,meshELoc,meshVT,Lam,Weights,AtoV,rVal,vIdx,rIdx,startIdx)
  
  names(output)<-c("meshArea","meshCLoc","meshELoc","meshVT","Lam","Weights","AtoV","rVal","vIdx","rIdx","startIdx")
  
  return(output)
}

sphericalComplex <- function(l,phi,theta) {
  m<-c(-l:l)
  
  CC<-sqrt((2*l+1)/(4*pi)*factorial(l-m)/factorial(l+m))
  
  if(l==0){
    P<-(as.matrix(legendre(l,cos(theta))))
  }
  
  else{
    P<-t(as.matrix(legendre(l,cos(theta))))
    mpos<-c(1:l)
    CC2<-((-1)^mpos*factorial(l-mpos)/factorial(l+mpos))
    
    if(length(CC2)==1){
      Ptemp<-(P[,-1]*(CC2))
      P<-cbind(Ptemp,P)
    }
    else{
      Ptemp<-(P[,-1]%*%diag(CC2))
      P<-cbind(Ptemp[,ncol(Ptemp):1],P)
    }
  }
  
  #Yhar<-P%*%exp(1i*kron(m,phi))%*%diag(CC)
  if(length(CC)==1){
    Yhar<-P*exp(1i*kron(t(m),phi))*(CC)
  }
  else{
    Yhar<-P*exp(1i*kron(t(m),phi))%*%diag(CC)
  }
  
  if(l==0){
    YderPhi=matrix(0,nrow(P),1)
  }
  else{
    lNext<-l+1
    mNext<-c(-lNext:lNext)
    CNext<-sqrt((2*lNext+1)/(4*pi)*factorial(lNext-mNext)/factorial(lNext+mNext))
    PNext<-t(as.matrix(legendre(lNext,cos(theta))))
    mposNext<-c(1:lNext)
    C2Next<-((-1)^mposNext)*factorial(lNext-mposNext)/factorial(lNext+mposNext)
    if(length(C2Next)==1){
      PNtemp<-PNext[,2:ncol(PNext)]*(C2Next)
    }
    else{
      PNtemp<-PNext[,2:ncol(PNext)]%*%diag(C2Next)
    }
    PNext<-cbind(PNtemp[,ncol(PNtemp):1],PNext)
    
    YderPhi=-0.5*(PNext[,3:ncol(PNext)]+PNext[,1:(ncol(PNext)-2)]%*%diag((l-m+1)*(l-m+2)))%*%diag(1/m)
    YderPhi[,l+1]=0
    YderPhi=YderPhi*exp(1i*kron(t(m),phi))%*%diag(CC)%*%diag(1i*m)
  }
  if(length(CC)==1){
    YderTheta=matrix(0,nrow(P),1)*exp(1i*kron(t(m),phi))*CC
  }
  else{
    P2=cbind(P[,2:ncol(P)],rep(0,nrow(P)))
    P0=cbind(rep(0,nrow(P)),P[,1:(ncol(P)-1)])
    
    Pder= 0.5*P2-0.5*P0%*%diag((l+m)*(l-m+1))
    
    YderTheta=Pder*exp(1i*kron(t(m),phi))%*%diag(CC) 
  }
  
  result<-list(Yhar,YderPhi,YderTheta)
  
  return(result)
}


sphericalReal <- function(l,phi,theta) {
  YharTmp<-sphericalComplex(l,phi,theta)[[1]]
  YderPhiTmp<-sphericalComplex(l,phi,theta)[[2]]
  YderThetaTmp<-sphericalComplex(l,phi,theta)[[3]]
  if(l==0){
    Yhar=Real(YharTmp)
    YderPhi=Real(YderPhiTmp)
    YderTheta=Real(YderThetaTmp)
  }
  else{
    Yhar=YharTmp
    YderPhi=YderPhiTmp
    YderTheta=YderThetaTmp
    m<-c(1:l)
    negM<-rev(m)
    posM<-l+1+m
    
    if(length(m)==1){
      Yhar[,negM]<-1/(1i*sqrt(2))*(YharTmp[,negM]-YharTmp[,posM]*((-1)^m))
      Yhar[,posM]<-1/sqrt(2)*(YharTmp[,posM]+YharTmp[,negM]*((-1)^m))
      
      YderPhi[,negM]<-1/(1i*sqrt(2))*(YderPhiTmp[,negM]-YderPhiTmp[,posM]*((-1)^m))
      YderPhi[,posM]<-1/sqrt(2)*(YderPhiTmp[,posM]+YderPhiTmp[,negM]*((-1)^m))
      
      YderTheta[,negM]<-1/(1i*sqrt(2))*(YderThetaTmp[,negM]-YderThetaTmp[,posM]*((-1)^m))
      YderTheta[,posM]<-1/sqrt(2)*(YderThetaTmp[,posM]+YderThetaTmp[,negM]*((-1)^m))
    }
    else{
      Yhar[,negM]<-1/(1i*sqrt(2))*(YharTmp[,negM]-YharTmp[,posM]%*%diag((-1)^m))
      Yhar[,posM]<-1/sqrt(2)*(YharTmp[,posM]+YharTmp[,negM]%*%diag((-1)^m))
      
      YderPhi[,negM]<-1/(1i*sqrt(2))*(YderPhiTmp[,negM]-YderPhiTmp[,posM]%*%diag((-1)^m))
      YderPhi[,posM]<-1/sqrt(2)*(YderPhiTmp[,posM]+YderPhiTmp[,negM]%*%diag((-1)^m))
      
      YderTheta[,negM]<-1/(1i*sqrt(2))*(YderThetaTmp[,negM]-YderThetaTmp[,posM]%*%diag((-1)^m))
      YderTheta[,posM]<-1/sqrt(2)*(YderThetaTmp[,posM]+YderThetaTmp[,negM]%*%diag((-1)^m))
    }
    Yhar<-Real(Yhar)
    YderPhi=Real(YderPhi)
    YderTheta=Real(YderTheta)
  }
  
  result=list(Yhar,YderPhi,YderTheta)
  
  return(result)
}

sphericalR3 <- function(l,phi,theta) {
  # Yhar, YderPhi, and YderTheta comes from complete function
  Yhar<-sphericalReal(l,phi,theta)[[1]]
  YderPhi<-sphericalReal(l,phi,theta)[[2]]
  YderTheta<-sphericalReal(l,phi,theta)[[3]]
  
  vx<-matrix(NA,nrow = nrow(YderPhi),ncol = 2*ncol(YderPhi))
  vy<-matrix(NA,nrow = nrow(YderPhi),ncol = 2*ncol(YderPhi))
  vz<-matrix(NA,nrow = nrow(YderPhi),ncol = 2*ncol(YderPhi))
  
  vecTheta<-YderTheta
  vecPhi<-YderPhi
  
  len<-ncol(YderPhi)
  
  vx[,1:len]<-kron(t(rep(1,len)),cos(theta)*cos(phi))*vecTheta+kron(t(rep(1,len)),-sin(phi))*vecPhi
  
  vy[,1:len]<-kron(t(rep(1,len)),cos(theta)*sin(phi))*vecTheta+kron(t(rep(1,len)),cos(phi))*vecPhi
  
  vz[,1:len]<-kron(t(rep(1,len)),-sin(theta))*vecTheta
  
  vecTheta<-(-YderPhi)
  vecPhi<-YderTheta
  
  # Not finished for exact expression of vx, vy, and vz
  vx[,(len+1):ncol(vx)]<-kron(t(rep(1,len)),cos(theta)*cos(phi))*vecTheta+kron(t(rep(1,len)),-sin(phi))*vecPhi
  
  vy[,(len+1):ncol(vy)]<-kron(t(rep(1,len)),cos(theta)*sin(phi))*vecTheta+kron(t(rep(1,len)),cos(phi))*vecPhi
  
  vz[,(len+1):ncol(vz)]<-kron(t(rep(1,len)),-sin(theta))*vecTheta
  
  result=list(vx,vy,vz)
  
  return(result)
}

calcBasis <- function(rhoMaxL, sigmaMaxL, vectorMaxL,meshCLoc,meshELoc) {
  c2s<-cart2sph(meshCLoc)
  
  az=c2s[,1]
  el=c2s[,2]
  
  c2sE<-cart2sph(meshELoc)
  
  azE=c2sE[,1]
  elE=c2sE[,2]
  
  basisRho=NULL
  basisSigma=NULL
  
  for (i in 0:rhoMaxL) {
    basisRho=cbind(basisRho,sphericalReal(i,az,pi/2-el)[[1]])
  }
  for (i in 0:sigmaMaxL) {
    basisSigma=cbind(basisSigma,sphericalReal(i,az,pi/2-el)[[1]])
  }
  
  if(vectorMaxL==0){
    basisDVx=matrix(0,3*nrow(basisRho),1)
    basisDVy=matrix(0,3*nrow(basisRho),1)
    basisDVz=matrix(0,3*nrow(basisRho),1)
  }
  else{
    basisDVx=NULL
    basisDVy=NULL
    basisDVz=NULL
    for (i in 1:vectorMaxL) {
      DVres<-sphericalR3(i,azE,pi/2-elE)
      basisDVx=cbind(basisDVx,DVres[[1]])
      basisDVy=cbind(basisDVy,DVres[[2]])
      basisDVz=cbind(basisDVz,DVres[[3]])
    }
  }
  
  idxRho=c(1:ncol(basisRho))
  idxSigma=idxRho[length(idxRho)]+c(1:ncol(basisSigma))
  idxDV=idxSigma[length(idxSigma)]+c(1:ncol(basisDVx))
  
  result=list(basisRho,basisSigma,basisDVx,basisDVy,basisDVz,idxRho,idxSigma,idxDV)
  
  names(result)<-c("basisRho","basisSigma","basisDVx","basisDVy","basisDVz","idxRho","idxSigma","idxDV")
  return(result)
}

##################

# Here rgeneric is being applied, some functions defined above need to be defined again inside rgeneric

# The core part for rgeneric is 
# initial: initialize parameters (labeled as "thetas" inside rgeneric)
# Q: define the precision matrix. Here graph works the same as Q, both defines the precision matrix
# mu: defines the meam values. I set all to 0
# log.prior: prior setting for all "thetas"

'inla.rgeneric.nonstationary.test.fixed.sigma' <- function(
  cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
          "log.prior", "quit"),
  theta = NULL) {
  
  
  require(Matrix)
  require(Rcpp)
  require(R.matlab)
  require(INLA)
  require(pracma)
  
  makeDV<-function(meshArea,invert){
    LEN<-length(meshArea)
    if(invert==0){
      RET<-sparseMatrix(i=c(1:LEN),j=c(1:LEN),x=meshArea)
    }
    else{
      RET<-sparseMatrix(i=c(1:LEN),j=c(1:LEN),x=1/meshArea)
    }
    return(RET)
  }
  
  #View(as.matrix(makeDV(restest$meshArea,invert = 0)))
  #View(as.matrix(makeDV(restest$meshArea,invert = 1)))
  
  ## makeDRho2
  
  makeDRho2 <- function(basisRho,PAR,invert) {
    if(ncol(basisRho)==1){
      RHO<-basisRho*PAR
      RHO<-as.numeric(RHO)
      LEN<-length(RHO)
    }
    else{
      RHO<-basisRho%*%PAR
      RHO<-as.numeric(RHO)
      LEN<-length(RHO)
    }
    if(invert==0){
      RES<-sparseMatrix(i=c(1:LEN),j=c(1:LEN),x=exp(RHO)^2)
    }
    else{
      RES<-sparseMatrix(i=c(1:LEN),j=c(1:LEN),x=exp(-RHO)^2)
    }
    return(RES)
  }
  
  #makeDRho2(c(0,0),0,0)
  
  ## makeDSigma
  
  makeDSigma <- function(basisSigma,PAR,invert) {
    if(ncol(basisSigma)==1){
      SIGMA<-basisSigma*PAR
      SIGMA<-as.numeric(SIGMA)
      LEN<-length(SIGMA)
    }
    else{
      SIGMA<-basisSigma%*%PAR
      SIGMA<-as.numeric(SIGMA)
      LEN<-length(SIGMA)
    }
    if(invert==0){
      RES<-sparseMatrix(i=c(1:LEN),j=c(1:LEN),x=exp(SIGMA))
    }
    else{
      RES<-sparseMatrix(i=c(1:LEN),j=c(1:LEN),x=exp(-SIGMA))
    }
    return(RES)
  }
  
  #makeDSigma(c(0,0),0,0)
  
  ## AH_maker
  
  ## A certain parameter needDiff in C++ code is set as 0, might need to change it in some other situations.
  AH_maker <- function(basisDVx,basisDVy,basisDVz,PAR,meshVLoc,meshTV,meshTT,meshStartIdx,meshRidx,meshRVal) {
    if(ncol(basisDVx)==1){
      VX<-basisDVx*PAR
      VX<-as.numeric(VX)
    }
    else{
      VX<-basisDVx%*%PAR
      VX<-as.numeric(VX)
    }
    if(ncol(basisDVy)==1){
      VY<-basisDVx*PAR
      VY<-as.numeric(VY)
    }
    else{
      VY<-basisDVy%*%PAR
      VY<-as.numeric(VY)
    }
    if(ncol(basisDVz)==1){
      VZ<-basisDVz*PAR
      VZ<-as.numeric(VZ)
    }
    else{
      VZ<-basisDVz%*%PAR
      VZ<-as.numeric(VZ)
    }
    
    NT<-dim(meshTV)[1]
    NV<-dim(meshVLoc)[1]
    NW<-length(meshRidx)
    
    sourceCpp("rAH_finalDer.cpp")
    
    AH_final<-makeAH(NT,NV,NW,as.numeric(t(as.matrix(meshVLoc))),as.numeric(t(as.matrix(meshTV))),as.numeric(t(as.matrix(meshTT))),VX,VY,VZ,0,VX,VY,VZ,meshStartIdx,meshRidx,meshRVal)
    
    colnames(AH_final)=c("I","J","V")
    
    AH_omit<-na.omit(AH_final)
    
    SPmat<-sparseMatrix(i=AH_omit$I,j=AH_omit$J,x=AH_omit$V)
    
    return(SPmat)
  }
  
  makeQ <- function(meshArea,basisRho,basisSigma,basisDVx,basisDVy,basisDVz,meshVLoc,meshTV,meshTT,PAR,par_idxRho,par_idxSigma,par_AH,meshStartIdx,meshRidx,meshRVal) {
    
    maxid<-max(c(par_idxRho,par_idxSigma,par_AH))
    
    newid<-rep(NA,maxid)
    
    SAHseq=1
    for (SAH in c(par_idxRho,par_AH)) {
      newid[SAH]=SAHseq
      SAHseq=SAHseq+1
    }
    
    par_idxRho<-newid[par_idxRho]
    par_AH<-newid[par_AH]
    
    DV<-makeDV(meshArea,0)
    DVInv<-makeDV(meshArea,1)
    
    DRho2<-makeDRho2(basisRho,PAR[par_idxRho],0)
    DRho2Inv<-makeDRho2(basisRho,PAR[par_idxRho],1)
    
    DSigmaInv<-makeDSigma(basisSigma,rep(0,length(par_idxSigma)),1)
    
    ### A slight modifications are made here to have all the parameters for basisSigma are fixed to be 0
    
    ### Original code is here
    ###  DSigmaInv<-makeDSigma(basisSigma,PAR[par_idxSigma],1)
    
    AH<-AH_maker(basisDVx,basisDVy,basisDVz,PAR[par_AH],meshVLoc,meshTV,meshTT,meshStartIdx,meshRidx,meshRVal)
    
    
    # print(DV)
    # print(DVInv)
    # print(DRho2)
    # print(DRho2Inv)
    # print(DSigmaInv)
    
    Aop<-(DVInv%*%DRho2)^(1/2)%*%(DV%*%DRho2Inv-AH)%*%(DSigmaInv/sqrt(4*pi))
    Q<-t(Aop)%*%Aop
    
    test_invQ<-inla.qinv(Q)
    
    diag_invQ<-diag(test_invQ)
    
    sd_diag_Q <- as(diag(sqrt(diag_invQ)), "sparseMatrix")  
    
    new_Q<-sd_diag_Q%*%Q%*%sd_diag_Q
    
    return(new_Q)
  }
  
  #Internal function
  
  interpret.theta <- function() {
    return(theta)
  }
  
  graph = function() {
    return(Q())
  }
  
  Q <- function() {
    
    PAR <- interpret.theta()
    
    # Final function for making Q matrix
    
    res<-makeQ(restest$meshArea,(nonstat_basis_rho[[1]]),(nonstat_basis_rho[[2]]),(nonstat_basis_rho[[3]]),(nonstat_basis_rho[[4]]),(nonstat_basis_rho[[5]]),VLOC,TV,TT,PAR,nonstat_basis_rho[[6]],nonstat_basis_rho[[7]],nonstat_basis_rho[[8]],restest$startIdx,restest$rIdx,restest$rVal)
    
    return(res)
  }
  
  mu <- function()
  {
    return(numeric(0))
  }
  
  log.norm.const <- function() {
    return(numeric(0))
  }
  
  log.prior <- function() {
    param = interpret.theta()
    
    #res <- sum(dgamma(param,0.00005,1000,log = TRUE))+sum(log(param))
    # Smaller variances are being applied here
    res <- sum(dnorm(param,0,1,log = TRUE))
    
    return(res)
  }
  
  initial <- function() {
    return(init_vect)
  }
  
  quit <- function() {
    return(invisible())
  }
  
  if (length(theta) == 0) {
    theta = initial()
  }
  
  res <- do.call(match.arg(cmd), args = list())
  return(res)
}

addlandsea_from_Optim <- function(meshCLoc,method="R") {
  c2s<-cart2sph(meshCLoc)
  
  az<-c2s[,1]
  el<-c2s[,2]
  
  lon<-az*360/(2*pi)
  lat<-el*360/(2*pi)
  
  if(method=="matlab"){
    xv=lon*0
    for (i in 1:8) {
      xv=xv+inpolygon(lon, lat, na.omit(as.numeric(land_sea_store$S[,,1][,i]$X)), na.omit(as.numeric(land_sea_store$S[,,1][,i]$Y)))
    }
    lMat=as.numeric(as.logical(xv))
    gridLandSeaE=rbind(rbind(lMat,lMat),lMat)
    lMatE<-as.numeric(gridLandSeaE)
    
    result<-list(lMat=lMat,lMatE=lMatE)
  }
  else{
    points <- data.frame(lon, lat)
    sf_use_s2(FALSE)
    pts <- st_as_sf(points, coords=1:2, crs=4326)
    lMat <- !is.na(as.numeric(st_intersects(pts, world)))
    lMat<-as.numeric(lMat)
    gridLandSeaE=rbind(rbind(lMat,lMat),lMat)
    lMatE<-as.numeric(gridLandSeaE)
    result<-list(lMat=lMat,lMatE=lMatE)
  }
  return(result)
}

addLandSeaCov <- function(basisinfo,lMat) {
  basisRho<-cbind(basisinfo$basisRho,lMat)
  currMax<-max(c(basisinfo$idxRho,basisinfo$idxSigma,basisinfo$idxDV))
  idxRho<-c(basisinfo$idxRho,(currMax+1))
  result<-list(basisRho=basisRho,idxRho=idxRho)
  return(result)
}

mat_get_distance <- function(precisionmaker,cBorder) {
  meshCLoc<-precisionmaker$meshCLoc
  dGrid<-rep(NA,nrow(meshCLoc))
  idx=floor(seq(1,nrow(cBorder[[1]][[1]]),length.out = 10000))
  for (i in 1:length(dGrid)) {
    dGrid[i]<-min(na.omit(as.numeric(pdist2(meshCLoc[i,],cBorder[[1]][[1]][idx,]))))
  }
  return(dGrid)
}

addBorder <- function(border_distance,basisRho,idxRho,basisSigma,idxSigma,maxinfo,bsize=0.03) {
  bIdx<- as.numeric(border_distance<=bsize)
  basisRho<-cbind(basisRho,bIdx)
  basisSigma<-cbind(basisSigma,bIdx)
  idxRho<-c(idxRho,(maxinfo+1))
  idxSigma<-c(idxSigma,(maxinfo+2))
  result<-list(basisRho=basisRho,idxRho=idxRho,basisSigma=basisSigma,idxSigma=idxSigma)
  return(result)
}

#### This is particularly set for situation that sigma set to 0 (like 4,0,4)

addLandSeaCov_NS <- function(basisinfo,als_optim) {
  # Rho
  tmpBas1<-basisinfo$basisRho
  tmpBas2<-tmpBas1
  
  for (i in 1:ncol(tmpBas1)) {
    tmpBas1[,i]<-tmpBas1[,i]*als_optim[[1]]
    tmpBas2[,i]<-tmpBas2[,i]*(1-als_optim[[1]])
  }
  
  res_basisRho<-cbind(tmpBas1,tmpBas2)
  
  #idxRho
  cMax<-max(basisinfo$idxDV)+1
  res_idxRho<-c(basisinfo$idxRho,c(cMax:(cMax+ncol(tmpBas1)-1)))
  
  #DVx
  tmpBas1<-basisinfo$basisDVx
  tmpBas2<-tmpBas1
  
  for (i in 1:ncol(tmpBas1)) {
    tmpBas1[,i]<-tmpBas1[,i]*als_optim[[2]]
    tmpBas2[,i]<-tmpBas2[,i]*(1-als_optim[[2]])
  }
  
  res_basisDVx<-cbind(tmpBas1,tmpBas2)
  
  #DVy
  tmpBas1<-basisinfo$basisDVy
  tmpBas2<-tmpBas1
  
  for (i in 1:ncol(tmpBas1)) {
    tmpBas1[,i]<-tmpBas1[,i]*als_optim[[2]]
    tmpBas2[,i]<-tmpBas2[,i]*(1-als_optim[[2]])
  }
  
  res_basisDVy<-cbind(tmpBas1,tmpBas2)
  
  #DVz
  tmpBas1<-basisinfo$basisDVz
  tmpBas2<-tmpBas1
  
  for (i in 1:ncol(tmpBas1)) {
    tmpBas1[,i]<-tmpBas1[,i]*als_optim[[2]]
    tmpBas2[,i]<-tmpBas2[,i]*(1-als_optim[[2]])
  }
  
  res_basisDVz<-cbind(tmpBas1,tmpBas2)
  
  #idxDv
  cMax<-max(res_idxRho)+1
  res_idxDv<-c(basisinfo$idxDV,c(cMax:(cMax+ncol(tmpBas1)-1)))
  
  result<-list(basisRho=res_basisRho,basisDVx=res_basisDVx,basisDVy=res_basisDVy,basisDVz=res_basisDVz,idxRho=res_idxRho,idxDV=res_idxDv)
  
  return(result)
}



makeDV<-function(meshArea,invert){
  LEN<-length(meshArea)
  if(invert==0){
    RET<-sparseMatrix(i=c(1:LEN),j=c(1:LEN),x=meshArea)
  }
  else{
    RET<-sparseMatrix(i=c(1:LEN),j=c(1:LEN),x=1/meshArea)
  }
  return(RET)
}

#View(as.matrix(makeDV(restest$meshArea,invert = 0)))
#View(as.matrix(makeDV(restest$meshArea,invert = 1)))

## makeDRho2

makeDRho2 <- function(basisRho,PAR,invert) {
  if(ncol(basisRho)==1){
    RHO<-basisRho*PAR
    RHO<-as.numeric(RHO)
    LEN<-length(RHO)
  }
  else{
    RHO<-basisRho%*%PAR
    RHO<-as.numeric(RHO)
    LEN<-length(RHO)
  }
  if(invert==0){
    RES<-sparseMatrix(i=c(1:LEN),j=c(1:LEN),x=exp(RHO)^2)
  }
  else{
    RES<-sparseMatrix(i=c(1:LEN),j=c(1:LEN),x=exp(-RHO)^2)
  }
  return(RES)
}

#makeDRho2(c(0,0),0,0)

## makeDSigma

makeDSigma <- function(basisSigma,PAR,invert) {
  if(ncol(basisSigma)==1){
    SIGMA<-basisSigma*PAR
    SIGMA<-as.numeric(SIGMA)
    LEN<-length(SIGMA)
  }
  else{
    SIGMA<-basisSigma%*%PAR
    SIGMA<-as.numeric(SIGMA)
    LEN<-length(SIGMA)
  }
  if(invert==0){
    RES<-sparseMatrix(i=c(1:LEN),j=c(1:LEN),x=exp(SIGMA))
  }
  else{
    RES<-sparseMatrix(i=c(1:LEN),j=c(1:LEN),x=exp(-SIGMA))
  }
  return(RES)
}

#makeDSigma(c(0,0),0,0)

## AH_maker

## A certain parameter needDiff in C++ code is set as 0, might need to change it in some other situations.
AH_maker <- function(basisDVx,basisDVy,basisDVz,PAR,meshVLoc,meshTV,meshTT,meshStartIdx,meshRidx,meshRVal) {
  if(ncol(basisDVx)==1){
    VX<-basisDVx*PAR
    VX<-as.numeric(VX)
  }
  else{
    VX<-basisDVx%*%PAR
    VX<-as.numeric(VX)
  }
  if(ncol(basisDVy)==1){
    VY<-basisDVx*PAR
    VY<-as.numeric(VY)
  }
  else{
    VY<-basisDVy%*%PAR
    VY<-as.numeric(VY)
  }
  if(ncol(basisDVz)==1){
    VZ<-basisDVz*PAR
    VZ<-as.numeric(VZ)
  }
  else{
    VZ<-basisDVz%*%PAR
    VZ<-as.numeric(VZ)
  }
  
  NT<-dim(meshTV)[1]
  NV<-dim(meshVLoc)[1]
  NW<-length(meshRidx)
  
  sourceCpp("rAH_finalDer.cpp")
  
  AH_final<-makeAH(NT,NV,NW,as.numeric(t(as.matrix(meshVLoc))),as.numeric(t(as.matrix(meshTV))),as.numeric(t(as.matrix(meshTT))),VX,VY,VZ,0,VX,VY,VZ,meshStartIdx,meshRidx,meshRVal)
  
  colnames(AH_final)=c("I","J","V")
  
  AH_omit<-na.omit(AH_final)
  
  SPmat<-sparseMatrix(i=AH_omit$I,j=AH_omit$J,x=AH_omit$V)
  
  return(SPmat)
}

makeQ <- function(meshArea,basisRho,basisSigma,basisDVx,basisDVy,basisDVz,meshVLoc,meshTV,meshTT,PAR,par_idxRho,par_idxSigma,par_AH,meshStartIdx,meshRidx,meshRVal) {
  
  maxid<-max(c(par_idxRho,par_idxSigma,par_AH))
  
  newid<-rep(NA,maxid)
  
  SAHseq=1
  for (SAH in c(par_idxRho,par_AH)) {
    newid[SAH]=SAHseq
    SAHseq=SAHseq+1
  }
  
  par_idxRho<-newid[par_idxRho]
  par_AH<-newid[par_AH]
  
  DV<-makeDV(meshArea,0)
  DVInv<-makeDV(meshArea,1)
  
  DRho2<-makeDRho2(basisRho,PAR[par_idxRho],0)
  DRho2Inv<-makeDRho2(basisRho,PAR[par_idxRho],1)
  
  DSigmaInv<-makeDSigma(basisSigma,rep(0,length(par_idxSigma)),1)
  
  ### A slight modifications are made here to have all the parameters for basisSigma are fixed to be 0
  
  ### Original code is here
  ###  DSigmaInv<-makeDSigma(basisSigma,PAR[par_idxSigma],1)
  
  AH<-AH_maker(basisDVx,basisDVy,basisDVz,PAR[par_AH],meshVLoc,meshTV,meshTT,meshStartIdx,meshRidx,meshRVal)
  
  
  # print(DV)
  # print(DVInv)
  # print(DRho2)
  # print(DRho2Inv)
  # print(DSigmaInv)
  
  Aop<-(DVInv%*%DRho2)^(1/2)%*%(DV%*%DRho2Inv-AH)%*%(DSigmaInv/sqrt(4*pi))
  Q<-t(Aop)%*%Aop
  
  test_invQ<-inla.qinv(Q)
  
  diag_invQ<-diag(test_invQ)
  
  sd_diag_Q <- as(diag(sqrt(diag_invQ)), "sparseMatrix")  
  
  new_Q<-sd_diag_Q%*%Q%*%sd_diag_Q
  
  return(new_Q)
}


makeAmat <- function(MESH, LOC,TV) {
  A.est<-inla.spde.make.A(mesh = MESH,loc = LOC)
  
  A_list<-list()
  for (i in 1:nrow(A.est)) {
    A_list[[i]]<-which(A.est[i,]!=0)
  }
  
  # A_list<-NULL
  # for (i in 1:nrow(A.est)) {
  #   searchtemp<-sort(which(A.est[i,]!=0))
  #   currchar<-NULL
  #   for(j in 1:length(searchtemp)){
  #     currchar<-paste(currchar,searchtemp[j],sep="-")
  #   }
  #   currchar<-substring(currchar,2)
  #   A_list<-c(A_list,currchar)
  # }
  
  searchindex=NULL
  for (i in 1:nrow(TV)) {
    indextemp=sort(TV[i,])
    indexchar=paste(indextemp[1],indextemp[2],indextemp[3],sep = "-")
    searchindex=c(searchindex,indexchar)
  }
  
  targetmat<-matrix(0,nrow(LOC),nrow(TV))
  
  # for (i in 1:nrow(targetmat)) {
  #   locres<-which(grepl(A_list[i],searchindex))[1]
  #   targetmat[i,locres]=1
  # }
  
  for (i in 1:nrow(targetmat)) {
    locres<-ifelse(length(A_list[[i]])==3,which(grepl(A_list[[i]][1],searchindex) & grepl(A_list[[i]][2],searchindex) & grepl(A_list[[i]][3],searchindex) ), ifelse(length(A_list[[i]])==2, which(grepl(A_list[[i]][1],searchindex) & grepl(A_list[[i]][2],searchindex)  ), which(grepl(A_list[[i]][1],searchindex)) ))[1]
    targetmat[i,locres]=1
  }
  
  targetmat<-as(targetmat, "sparseMatrix")  
  
  return(targetmat)
}

# Changed expandA functions to real "sparseMatrix"
expandA <- function(Amat,repl) {
  dimreal<-dim(Amat)
  dimrepl<-dim(Amat)*repl
  #targetmat<-as(matrix(0,dimrepl[1],dimrepl[2]), "sparseMatrix")  
  targetmat<-sparseMatrix(dimrepl[1],dimrepl[2])
  for (i in 1:repl) {
    targetmat[(dimreal[1]*(i-1)+1):(dimreal[1]*i),(dimreal[2]*(i-1)+1):(dimreal[2]*i)]<-Amat
  }
  return(targetmat)
}

#####################

### Pardiso not working anymore, change to taucs

inla.setOption(smtp="taucs")

library(ggplot2)
library(sf)
library("rnaturalearth")
library("rnaturalearthdata")

world <- ne_countries(scale = "medium", returnclass = "sf")

locmesh = inla.mesh.create(globe = 10)
restest<-PrecisionMaker(locmesh$loc,locmesh$graph$tv,locmesh$graph$tt)
meshClonlat<-inla.mesh.map(restest$meshCLoc,"longlat",inverse = FALSE)
meshCdf<-as.data.frame(meshClonlat)

### Try to simulate the data

nonstat_basis_rho<-calcBasis(1,0,1,as.matrix(restest$meshCLoc),as.matrix(restest$meshELoc))

border_store<-readMat("border.mat")

als_optim<-addlandsea_from_Optim(restest$meshCLoc)

extendbasis<-addLandSeaCov_NS(nonstat_basis_rho,als_optim)

stat_distance<-mat_get_distance(restest,border_store$cBorder)

maxinfo<-max(c(extendbasis$idxRho,nonstat_basis_rho$idxSigma,extendbasis$idxDV))

extendRho<-addBorder(stat_distance,extendbasis$basisRho,extendbasis$idxRho,nonstat_basis_rho$basisSigma,nonstat_basis_rho$idxSigma,maxinfo)

combine_basis<-list(basisRho=extendRho$basisRho,basisSigma=nonstat_basis_rho$basisSigma,basisDVx=extendbasis$basisDVx,basisDVy=extendbasis$basisDVy,basisDVz=extendbasis$basisDVz,idxRho=extendRho$idxRho,idxSigma=nonstat_basis_rho$idxSigma,idxDV=extendbasis$idxDV)

maxid<-max(c(combine_basis$idxRho,combine_basis$idxSigma,combine_basis$idxDV))

set.seed(1234)
partest=runif(21)
partest[1:4]=rnorm(4,1,0.5)
partest[11:14]=rnorm(4,-1,0.5)
partest[21]=0
partest[5:10]=rnorm(6,1,0.5)
partest[15:20]=rnorm(6,-1,0.5)

saveTrue<-partest

test_precmat<-makeQ(restest$meshArea,(combine_basis[[1]]),(combine_basis[[2]]),(combine_basis[[3]]),(combine_basis[[4]]),(combine_basis[[5]]),locmesh$loc,locmesh$graph$tv,locmesh$graph$tt,partest,combine_basis[[6]],combine_basis[[7]],combine_basis[[8]],restest$startIdx,restest$rIdx,restest$rVal)

test_invQ<-inla.qinv(test_precmat)

diag_invQ<-diag(test_invQ)

sd_diag_Q <- as(diag(sqrt(diag_invQ)), "sparseMatrix")  

new_Q<-sd_diag_Q%*%test_precmat%*%sd_diag_Q

test_new_Qinv<-inla.qinv(new_Q)

diag_newQ_inv<-diag(test_new_Qinv)

summary(diag_newQ_inv)

testQ_Chol<-Cholesky(new_Q, LDL = FALSE, perm = TRUE)


etest=rep(0,dim(new_Q)[1])
#etest[5531]=1
etest[280]=1
covtest<-solve(new_Q,etest)
covtest<-as.numeric(covtest)

diag_invQ<-diag(test_new_Qinv)
#corrtest<-covtest/sqrt(diag_invQ[5531]*diag_invQ)
corrtest<-covtest/sqrt(diag_invQ[280]*diag_invQ)

summary(corrtest)

precpdata<-cbind(meshCdf,corrtest)

etest=etest+1

precpdata<-cbind(precpdata,etest)

precpdata<-as.data.frame(precpdata)

colnames(precpdata)<-c("lon","lat","PRECT","SIZE")

### Set index for replicates here
index<<-as.numeric(Sys.getenv("SGE_TASK_ID"))
replicate_level<-seq(1,100,by=1)

z<-NULL

y<-NULL

set.seed(321)

for (i in 1:replicate_level[index]) {
  ztemp <- rnorm(nrow(new_Q))
  
  ytemp <- as.vector(solve(testQ_Chol, solve(testQ_Chol,
                                             ztemp,
                                             system = 'Lt'
  ), system = 'Pt'))
  
  z<-rbind(z,ztemp)
  
  y<-rbind(y,ytemp)
}

##### Important! Add one more random effects for the y as nugget effects

##### Here only use the same effects for land and ocean as the test. Could change it into different settings if necessary

y_rand<-matrix(NA,nrow(y),ncol(y))

set.seed(1230)
for (i in 1:nrow(y_rand)) {
  for (j in 1:ncol(y_rand)) {
    y_rand[i,j]<-rnorm(1,0,sqrt(0.05))
  }
}

y<-y+y_rand

# Now, start INLA to estimate non-stationary effect with land/ocean

meshtest<-inla.mesh.create(globe = 10)

restest<-PrecisionMaker(meshtest$loc,meshtest$graph$tv,meshtest$graph$tt)

#stat_basis_rho<-calcBasis(0,0,0,as.matrix(restest$meshCLoc),as.matrix(restest$meshELoc))

####################

nonstat_basis_rho<-calcBasis(1,0,1,as.matrix(restest$meshCLoc),as.matrix(restest$meshELoc))

border_store<-readMat("border.mat")

als_optim<-addlandsea_from_Optim(restest$meshCLoc)

extendbasis<-addLandSeaCov_NS(nonstat_basis_rho,als_optim)

stat_distance<-mat_get_distance(restest,border_store$cBorder)

maxinfo<-max(c(extendbasis$idxRho,extendbasis$idxSigma,extendbasis$idxDV))

extendRho<-addBorder(stat_distance,extendbasis$basisRho,extendbasis$idxRho,nonstat_basis_rho$basisSigma,nonstat_basis_rho$idxSigma,maxinfo)

combine_basis<-list(basisRho=extendRho$basisRho,basisSigma=nonstat_basis_rho$basisSigma,basisDVx=extendbasis$basisDVx,basisDVy=extendbasis$basisDVy,basisDVz=extendbasis$basisDVz,idxRho=extendRho$idxRho,idxSigma=nonstat_basis_rho$idxSigma,idxDV=extendbasis$idxDV)

savedegree<-length(c(combine_basis$idxRho,combine_basis$idxDV))

set.seed(1234)
partest=runif(savedegree)*0.1-0.05

maxid<-max(c(combine_basis$idxRho,combine_basis$idxSigma,combine_basis$idxDV))


library(plot.matrix)

library("INLA")

stat.model.test <- inla.rgeneric.define(inla.rgeneric.nonstationary.test.fixed.sigma, init_vect=partest, restest=restest, nonstat_basis_rho=combine_basis,VLOC=meshtest$loc,TV=meshtest$graph$tv,TT=meshtest$graph$tt)

#testMat<-readMat("data.mat")

sf_use_s2(FALSE) 
pts <- st_as_sf(meshCdf, coords=1:2, crs=4326)
landsea<-!is.na(as.numeric(st_intersects(pts, world)))
landidx<-which(landsea)
seaidx<-which(!landsea)

globe_prect_y<-NULL

for (i in 1:replicate_level[index]) {
  globy_land<-as.numeric(y[i,])
  globy_sea<-as.numeric(y[i,])
  
  globy_land[seaidx]=NA
  globy_sea[landidx]=NA
  
  globe_temp<-cbind(globy_land,globy_sea)
  globe_prect_y<-rbind(globe_prect_y,globe_temp)
}

testdf<-data.frame(idx=rep(c(1:(length(globe_prect_y)/replicate_level[index])),replicate_level[index]),y=globe_prect_y)

idx.rep<-NULL

for (i in 1:replicate_level[index]) {
  idx.rep<-c(idx.rep,rep(i,dim(meshCdf)[1]))
}

loc_5rep<-NULL

for (i in 1:replicate_level[index]) {
  loc_5rep<-rbind(loc_5rep,as.matrix(restest$meshCLoc))
}

#locinfo<-readMat("locloc.mat")

A.est<-makeAmat(meshtest,restest$meshCLoc,meshtest$graph$tv)

A.est<-expandA(A.est,replicate_level[index])

s.index<-inla.spde.make.index(name = "S",n.spde = nrow(meshtest$graph$tv),n.repl = replicate_level[index])

globetest<-inla.stack(data=list(y=globe_prect_y),A=list(A.est,1),effects=list(s.index,Intercept=rep(1,nrow(globe_prect_y))))

formula<-as.formula("y~-1+Intercept+f(S,model=stat.model.test,replicate=S.repl)")


#### Don't know why but impossible to use auto
savetime<-system.time(mres_non <- inla(formula, data = inla.stack.data(globetest), 
                                       family = c("gaussian","gaussian"),control.predictor = 
                                         list(A=inla.stack.A(globetest),
                                              compute=TRUE),control.inla = 
                                         list(int.strategy = "eb"),verbose = TRUE,control.compute = list(dic=TRUE)))

savePost<-mres_non$summary.hyperpar[3:nrow(mres_non$summary.hyperpar),1]

saveHyper<-mres_non$summary.hyperpar

SaveMargHyper<-mres_non$marginals.hyperpar

saveRMSE_hyper<-sqrt((saveTrue-savePost)^2)

trueVal<-as.numeric(t(y))

estVal<-mres_non$summary.fitted.values[1:length(trueVal),1]

saveRMSE_est<-sqrt((trueVal-estVal)^2)

saveres<-list(saveRMSE_est=saveRMSE_est,saveRMSE_hyper=saveRMSE_hyper,saveTrue=saveTrue,savePost=savePost,saveHyper=saveHyper,SaveMargHyper=SaveMargHyper)

saveRDS(saveres,paste0("rep_",replicate_level[index],"_degree_1","_NS_LS_gaus",".RData"))