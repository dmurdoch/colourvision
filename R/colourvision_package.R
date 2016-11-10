spec.denoise <- function(specfiles, spar = 0.7){
  result<-NULL
  name<-names(specfiles)
  for (i in 2:ncol(specfiles)){
    temp<-smooth.spline(x = specfiles[,1], y = specfiles[,i], spar = spar)$y
    result<-cbind(result,temp)
  }
  result<-data.frame(specfiles[,1], result)
  names(result)<-name
  return(result)
}

energytoflux<-function (datum) {
  r<-datum[,2]*datum[,1]*8.357922e-05
  r<-data.frame(datum[,1],r)
  names(r)<-names(datum)
  return(r)
}

logistic<-function(x=seq(300,700,1), x0, L, k) {
  y<-L/(1+exp(-k*(x-x0)))
  r<-data.frame(x,y)
  colnames(r)<-c("Wavelength","R")
  return(r)
}

photor <- function (lambda.max,
               lambda = seq(300,700,1),
               beta.band = FALSE) {
  A = 69.7
  B = 28
  b = 0.922
  C = -14.9
  c = 1.104
  D = 0.674
  
  nphotor<-length(lambda.max)
  r <- matrix(ncol=nphotor, nrow = length(lambda))
  
  for (k in 1:nphotor) {
  
  a = 0.8795 + 0.0459*exp(-1*(((lambda.max[[k]]-300)^2)/11940))
  
    for(i in 1:length(lambda)) {
      
      x <- lambda.max[[k]]/lambda[[i]]
      r[i,k] <- 1/( exp(A*(a-x)) + exp(B*(b-x)) + exp(C*(c-x)) + D)
    }
  
    r[,k]<-r[,k]/max(r[,k])
  
  #beta band
  if (beta.band == TRUE) {
    beta <- vector(length = length(lambda))
    for(i in 1:length(lambda)) {
      lambda.max.beta = 189 + 0.315*lambda.max[[k]]
      b = -40.5 + 0.195*lambda.max[[k]]
      A.beta = 0.26
      beta[[i]] <- A.beta*exp(-((lambda[[i]]-lambda.max.beta)/b)^2)
    }
  }
  
  if (beta.band == TRUE) {r[,k] <- r[,k] + beta}
  
  }
  
  r <- cbind(lambda, r)
  r<-data.frame(r)
  
  colnames(r)<-c("Wavelength", paste("lambda.max", lambda.max, sep=""))
  return(r)
  
}


Q <- function (R,I,C,interpolate,nm) {
  
  if(interpolate == TRUE) {
  
  I <- data.frame(nm, approx(x = I[,1], y = I[,2], xout = nm, method="linear")$y)
  R <- data.frame(nm, approx(x = R[,1], y = R[,2], xout = nm, method="linear")$y)
  C <- data.frame(nm, approx(x = C[,1], y = C[,2], xout = nm, method="linear")$y)
  }
  
  if (interpolate == FALSE) {
    test<-c(length(R[,1])==length(I[,1]), length(R[,1])==length(C[,1]), length(C[,1])==length(I[,1]))
    ifelse(any(test==FALSE), yes=stop("Uneven number of rows. Use 'interpolate=TRUE'.", call.=FALSE), no="")
    ifelse(any(c(R[,1]==I[,1],R[,1]==C[,1],I[,1]==C[,1]))==FALSE, yes=stop("Different wavelenght values of model parameters. Use 'interpolate=TRUE'.", call. = FALSE), no="")
    a<-R[,1]
    b<-R[2:length(R[,1]),1]
    ifelse(sd(b-a[1:length(a)-1], na.rm = T)!=0, yes=stop("Uneven wavelenght intervals. Use 'interpolate=TRUE'.", call. = FALSE), no="")
  }
  r <- sum(I[,2]*R[,2]*C[,2])
  return(r)
}

Qr <- function(R, I, Rb, C, interpolate, nm) {
  Qr <- Q(I=I,R=R,C=C,interpolate=interpolate,nm=nm)
  QEr <- Q(I=I,R=Rb,C=C,interpolate=interpolate,nm=nm)
  r<-Qr/QEr
  return(r)
}

EMmodel <- function (photo=c("tri","tetra"),
                     R,
                     I,
                     Rb,
                     C,
                     interpolate=TRUE,
                     nm=seq(300,700,1))
{
  nphoto=ncol(C)-1
  ifelse(photo=="tetra" && nphoto!=4, yes=warning("Tetrachromatic model but C argument does not have four sensitivity curves.", call.=FALSE), no="")
  ifelse(photo=="tri" && nphoto!=3, yes=warning("Trichromatic model but C argument does not have three sensitivity curves.", call.=FALSE), no="")
  ifelse(any(ncol(I)>2), yes=warning("I argument with more than two columns. Only the first two will be used.", call. = FALSE), no="")
  ifelse(any(ncol(Rb)>2), yes=warning("Rb argument with more than two columns. Only the first two will be used.", call. = FALSE), no="")
  maxR<-apply(data.frame(R[,2:ncol(R)]), 2, max)
  ifelse( any(maxR<=1)&&max(Rb[,2])>1 || any(maxR>1)&&max(Rb[,2])<=1, yes=warning("There seems to be a problem with input files. 'R' and 'Rb' must be in the same scale. Both must be either in percentage (0-100%) or proportion (0-1).", call. = FALSE), no="")
  
  internal<- function (photo=c("tri","tetra"),
                       R,
                       I,
                       Rb,
                       C,
                       interpolate,
                       nm)
  {
    
    log<-"log"
    S1 <- Qr(I=I, R=R, Rb=Rb, C=C[,c(1,2)], interpolate=interpolate, nm=nm)
    if (log=="log") {S1 <- log(S1)}
    
    S2 <- Qr(I=I, R=R, Rb=Rb, C=C[,c(1,3)], interpolate=interpolate, nm=nm)
    if (log=="log") {S2 <- log(S2)}
    
    S3 <- Qr(I=I, R=R, Rb=Rb, C=C[,c(1,4)], interpolate=interpolate, nm=nm)
    if (log=="log") {S3 <- log(S3)}
    
    if(photo=="tetra") {S4 <- Qr(I=I, R=R, Rb=Rb, C=C[,c(1,5)], interpolate=interpolate, nm=nm)
    if (log=="log") {S4 <- log(S4)}
    }
    
    if(photo=="tri") {S <-c(S1,S2,S3)}
    if(photo=="tetra") {S <-c(S1,S2,S3,S4)}
    
    u <- S[[1]]/sum(S)
    s <- S[[2]]/sum(S)
    m <- S[[3]]/sum(S)
    if(photo=="tetra") {l <- S[[4]]/sum(S)}
    
    if(photo=="tri") {
      x <- (2/3)*((sqrt(3)/2)*(m-u) )
      y <- (2/3)*(s-0.5*(u+m))
    }
    
    if(photo=="tetra") {
      x <- ((1-2*s-m-u)/2)*sqrt(3/2)
      y <- (-1+3*m+u)/(2*sqrt(2))
      z <- u - (1/4)
    }
    
    if(photo=="tetra") {
      
      deltaSo<-sqrt(x^2+y^2+z^2)
      
      r<-c(S1, S2, S3, S4, u, s, m, l, x, y, z, deltaSo)
      r<-as.vector(r)
      names(r)<-c("f1", "f2", "f3", "f4", "u", "s", "m", "l", "x", "y", "z", "deltaSo")  
    }
    
    if(photo=="tri") {
      
      deltaSo<-sqrt(x^2+y^2)
      
      r<-c(S1, S2, S3, u, s, m, x, y, deltaSo)
      r<-as.vector(r)
      names(r)<-c("f1", "f2", "f3", "u", "s", "m", "x", "y", "deltaSo") 
    }
    
    return(r) 
  }
  
  n.spectra<-ncol(R)-1
  R.list<-vector(length=n.spectra, "list")
  for (i in 1:n.spectra) {R.list[[i]]<-data.frame(R[,1],R[,i+1])}
  
  r<-sapply(X=R.list, FUN=internal, photo=photo, I=I,
            Rb=Rb, C=C, interpolate=interpolate, nm=nm)
  r<-as.data.frame(t(r))
  rownames(r)<-names(R)[2:ncol(R)]
  
  if(photo=="tri") {
    test <- r[,c("f1","f2","f3")]
    ifelse (any(test<0),
            yes=warning("Log-transformation might be generating photoreceptor outputs < 0.", call.=FALSE), no="")
  }
  if(photo=="tetra") {
    test <- r[,c("f1","f2","f3","f4")]
    ifelse (any(test<0),
            yes=warning("Log-transformation might be generating photoreceptor outputs < 0.",
                        call.=FALSE), no="")
  }
  
  attr(r, "model name") <- "Endler and Mielke model"
  attr(r, "number photor. types") <- photo
  return(r)
  
}


noise_e <- function (v,n) {
  e<-v/sqrt(n)
  return(e)
}

RNLmodel <- function(model = c("linear", "log"),
                   photo=c("di", "tri","tetra"), R, Rb, I, C,
                   noise=FALSE, v, n, e,
                   interpolate=TRUE,
                   nm=seq(300,700,1)) {
  dependent <- FALSE
  nphoto=ncol(C)-1
  ifelse(photo=="tetra" && nphoto!=4, yes=warning("Tetrachromatic model but 'C' argument has a number of sensitivity curves different than four.", call.=FALSE), no="")
  ifelse(photo=="tri" && nphoto!=3, yes=warning("Trichromatic model but 'C' argument has a number of sensitivity curves different than three.", call.=FALSE), no="")
  ifelse(photo=="di" && nphoto!=2, yes=warning("Dichromatic model but 'C' argument has a number of sensitivity curves different than two.", call.=FALSE), no="")
  ifelse(photo=="tetra" && noise==T && length(e)!=4, yes=warning("Tetrachromatic model but 'e' argument has a number of parameters different than four.", call.=FALSE), no="")
  ifelse(photo=="tri" && noise==T && length(e)!=3, yes=warning("Trichromatic model but 'e' argument has a number of parameters different than three.", call.=FALSE), no="")
  ifelse(photo=="di" && noise==T && length(e)!=2, yes=warning("Dichromatic model but 'e' argument has a number of parameters different than two.", call.=FALSE), no="")
  ifelse(photo=="tetra" && noise==F && length(n)!=4, yes=warning("Tetrachromatic model but 'n' argument has a number of parameters different than four.", call.=FALSE), no="")
  ifelse(photo=="tri" && noise==F && length(n)!=3, yes=warning("Trichromatic model but 'n' argument does has a number of parameters different than three.", call.=FALSE), no="")
  ifelse(photo=="di" && noise==F && length(n)!=2, yes=warning("Dichromatic model but 'n' argument does has a number of parameters different than two.", call.=FALSE), no="")
  ifelse(any(ncol(I)>2), yes=warning("'I' argument with more than two columns. Only the first two will be used.", call. = FALSE), no="")
  ifelse(any(ncol(Rb)>2), yes=warning("'Rb' argument with more than two columns. Only the first two will be used.", call. = FALSE), no="")
  maxR<-apply(data.frame(R[,2:ncol(R)]), 2, max)
  ifelse( any(maxR<=1)&&max(Rb[,2])>1 || any(maxR>1)&&max(Rb[,2])<=1, yes=warning("There seems to be a problem with input files. 'R' and 'Rb' must be in the same scale. Both must be either in percentage (0-100%) or proportion (0-1).", call. = FALSE), no="")
  
  internal <-function(model,
                      photo, R, Rb, I, C,
                      noise, dependent, v, n, e,
                      interpolate,
                      nm) {
  
    #Sr values for colour patch 1 (R)
    S1.1 <- Qr(I=I, R=R, Rb=Rb, C=C[,c(1,2)], interpolate=interpolate, nm=nm)
    if (model=="linear") {S1.1 <- S1.1}
    if (model=="log") {S1.1 <- log10(S1.1)}
    
    S2.1 <- Qr(I=I, R=R, Rb=Rb, C=C[,c(1,3)], interpolate=interpolate, nm=nm)
    if (model=="linear") {S2.1 <- S2.1}
    if (model=="log") {S2.1 <- log10(S2.1)}
    
    if(photo=="tri"||photo=="tetra") {
      S3.1 <- Qr(I=I, R=R, Rb=Rb, C=C[,c(1,4)], interpolate=interpolate, nm=nm)
      if (model=="linear") {S3.1 <- S3.1}
      if (model=="log") {S3.1 <- log10(S3.1)}
    }
    
    if(photo=="tetra") {
      S4.1 <- Qr(I=I, R=R, Rb=Rb, C=C[,c(1,5)], interpolate=interpolate, nm=nm)
      if (model=="linear") {S4.1 <- S4.1}
      if (model=="log") {S4.1 <- log10(S4.1)}
    }
    
    if (noise==FALSE) {
      if (dependent==FALSE) {
      noiseS1<-noise_e(v=v,n=n[[1]])
      noiseS2<-noise_e(v=v,n=n[[2]])
      if(photo=="tri"||photo=="tetra") {noiseS3<-noise_e(v=v,n=n[[3]])}
      if(photo=="tetra") {noiseS4<-noise_e(v=v,n=n[[4]])}
      }
      if (dependent==TRUE) {
        if(interpolate == TRUE) {ifelse((nm[[2]]-nm[[1]])!=1,
                                       yes=stop("When 'dependent = TRUE' interpolation interval must be 1nm.", call. = FALSE), no="")}
        if(interpolate == FALSE) {ifelse((R[2,1]-R[1,1])!=1,
                                        yes=stop("When 'dependent = TRUE' wavelength interval must be 1nm. Consider using 'interpolate = TRUE'.", call. = FALSE),no="")}
        noiseS1 <- sqrt((v^2/n[[1]]) + (1/Q(R=R,I=I,C=C[,c(1,2)],interpolate=interpolate,nm=nm)))
        noiseS2 <- sqrt((v^2/n[[2]]) +  (1/Q(R=R,I=I,C=C[,c(1,3)],interpolate=interpolate,nm=nm)))
        if(photo=="tri"||photo=="tetra") {noiseS3 <- sqrt((v^2/n[[3]]) + (1/Q(R=R,I=I,C=C[,c(1,4)],interpolate=interpolate,nm=nm)))}
        if(photo=="tetra") {noiseS4 <- sqrt((v^2/n[[4]]) + (1/Q(R=R,I=I,C=C[,c(1,5)],interpolate=interpolate,nm=nm)))}
      }
    }
    
    if (noise==TRUE) {  
      noiseS1<-e[[1]]
      noiseS2<-e[[2]]
      if(photo=="tri"||photo=="tetra") {noiseS3<-e[[3]]}
      if(photo=="tetra") {noiseS4<-e[[4]]}    
    }
    
    if(photo=="tetra"){
      
      #Renoult et al. 2015
      A1<-sqrt( 1/( (noiseS3^2)+(noiseS4^2) ) )
      B1<-sqrt( ( (noiseS3^2)+(noiseS4^2) ) / ( (noiseS2*noiseS3)^2 + (noiseS2*noiseS4)^2 + (noiseS3*noiseS4)^2 )  )
      C1<-sqrt( ((noiseS2*noiseS3)^2 + (noiseS2*noiseS4)^2 + (noiseS3*noiseS4)^2) / 
                  ( (noiseS2*noiseS3*noiseS4)^2 +
                      (noiseS1*noiseS3*noiseS4)^2 +
                      (noiseS1*noiseS2*noiseS4)^2 +
                      (noiseS1*noiseS2*noiseS3)^2 )  )
      
      a1<- (noiseS3^2) / ( (noiseS3^2) + (noiseS4^2) )
      b1<- (noiseS4^2) / ( (noiseS3^2)+ (noiseS4^2) )
      
      a2<- ( (noiseS2*noiseS3)^2 )/( (noiseS2*noiseS3)^2 + (noiseS2*noiseS4)^2 + (noiseS3*noiseS4)^2 )
      b2<- ( (noiseS2*noiseS4)^2 )/( (noiseS2*noiseS3)^2 + (noiseS2*noiseS4)^2 + (noiseS3*noiseS4)^2 )
      c2<- ( (noiseS3*noiseS4)^2 )/( (noiseS2*noiseS3)^2 + (noiseS2*noiseS4)^2 + (noiseS3*noiseS4)^2 )
      
  
      x1<-A1*(S4.1-S3.1)
      y1<-B1*(S2.1-(a1*S4.1+b1*S3.1))
      z1<-C1*(S1.1-(a2*S4.1+b2*S3.1+c2*S2.1))
  
      delta_e<-sqrt((x1)^2+(y1)^2+(z1)^2)
      
      r<-c(noiseS1,noiseS2,noiseS3,noiseS4,S1.1,S2.1,S3.1,S4.1,x1,y1,z1,delta_e)
      r<-as.vector(r)
      names(r)<-c("e1","e2","e3","e4","f1","f2","f3","f4","x","y","z","deltaSo")
    }
    
    if(photo=="tri"){
      
      #Hempel de Ibarra, Giurfa and Vorobyev (2001)
      A1<-sqrt(1/(noiseS2^2+noiseS3^2))
      B1<-sqrt((noiseS2^2+noiseS3^2)/(((noiseS1*noiseS2)^2)+((noiseS1*noiseS3)^2)+((noiseS2*noiseS3)^2)))
      a1<-noiseS2^2/(noiseS2^2+noiseS3^2)
      b1<-noiseS3^2/(noiseS2^2+noiseS3^2)
      X1 <- A1*(S3.1-S2.1)
      Y1 <- B1*(S1.1-(a1*S3.1+b1*S2.1))
  
      delta_e <- sqrt(X1^2+Y1^2)
      
      r<-c(noiseS1,noiseS2,noiseS3,S1.1,S2.1,S3.1,X1,Y1,delta_e)
      r<-as.vector(r)
      names(r)<-c("e1","e2","e3","f1","f2","f3","x","y","deltaSo")
    }
    
    if(photo=="di"){
    
      #Hempel de Ibarra, Giurfa and Vorobyev (2001)
      X1 <- sqrt(1/(noiseS1^2+noiseS2^2))*(S2.1-S1.1)
      delta_e <- sqrt ( X1^2 )
      
      r<-c(noiseS1,noiseS2,S1.1,S2.1,X1,delta_e)
      r<-as.vector(r)
      names(r)<-c("e1","e2","f1","f2","x","deltaSo")
      
    }
  
  return(r)
      
  }
  
  
  n.spectra<-ncol(R)-1
  R.list<-vector(length=n.spectra, "list")
  for (i in 1:n.spectra) {R.list[[i]]<-data.frame(R[,1],R[,i+1])}
  
  r<-sapply(X=R.list, FUN=internal, model=model, photo=photo,
            I=I,
            Rb=Rb, C=C,
            noise=noise,
            dependent=dependent,
            v=v,
            n=n,
            e=e,
            interpolate=interpolate,
            nm=nm)
  r<-as.data.frame(t(r))
  rownames(r)<-names(R)[2:ncol(R)]
  
  if(model=="log" && photo=="di") {
    test <- r[,c("f1","f2")]
    ifelse (any(test<0),
            yes=warning("Photoreceptor output < 0.", call.=FALSE), no ="")
  }
  
  if(model=="log" && photo=="tri") {
    test <- r[,c("f1","f2","f3")]
    ifelse (any(test<0),
            yes=warning("Photoreceptor output < 0.", call.=FALSE), no ="")
  }
  if(model=="log" && photo=="tetra") {
    test <- r[,c("f1","f2","f3","f4")]
    ifelse (any(test<0),
            yes=warning("Photoreceptor output < 0.", call.=FALSE), no="")
  }
  
  attr(r, "model name") <- "Receptor noise limited model"
  attr(r, "photoreceptor function") <- model
  attr(r, "number photor. types") <- photo
  attr(r, "noise calculated") <- noise
  if (noise==FALSE) {
    attr(r, "dependent") <- dependent
    attr(r, "v") <- v
    attr(r, "n") <- n
    }
  return(r)
  
}


CTTKmodel <- function (photo=c("tri", "tetra"),
                  R,
                  I,
                  Rb,
                  C,
                  interpolate=TRUE,
                  nm=seq(300,700,1))
{

  nphoto=ncol(C)-1
  ifelse(photo=="tetra" && nphoto!=4, yes=warning("Tetrachromatic model but C argument does not have four sensitivity curves.", call.=FALSE), no="")
  ifelse(photo=="tri" && nphoto!=3, yes=warning("Trichromatic model but C argument does not have three sensitivity curves.", call.=FALSE), no="")
  ifelse(any(ncol(I)>2), yes=warning("I argument with more than two columns. Only the first two will be used.", call. = FALSE), no="")
  ifelse(any(ncol(Rb)>2), yes=warning("Rb argument with more than two columns. Only the first two will be used.", call. = FALSE), no="")
  maxR<-apply(data.frame(R[,2:ncol(R)]), 2, max)
  ifelse( any(maxR<=1)&&max(Rb[,2])>1 || any(maxR>1)&&max(Rb[,2])<=1, yes=warning("There seems to be a problem with input files. 'R' and 'Rb' must be in the same scale. Both must be either in percentage (0-100%) or proportion (0-1).", call. = FALSE), no="")
  
  internal <- function (photo,
                        R,
                        I,
                        Rb,
                        C,
                        interpolate,
                        nm) {
  
    Puv <- Qr(I=I, R=R, Rb=Rb, C=C[,c(1,2)], interpolate=interpolate, nm=nm)
    Pb <- Qr(I=I, R=R, Rb=Rb, C=C[,c(1,3)], interpolate=interpolate, nm=nm)
    Pg <- Qr(I=I, R=R, Rb=Rb, C=C[,c(1,4)], interpolate=interpolate, nm=nm)
    if(photo=="tetra") {
      Pr <- Qr(I=I, R=R, Rb=Rb, C=C[,c(1,5)], interpolate=interpolate, nm=nm)
    }
    
    Euv <- Puv/(1+Puv)
    Eb <- Pb/(1+Pb)
    Eg <- Pg/(1+Pg)
    if(photo=="tetra") {Er <- Pr/(1+Pr)}
  
    if(photo=="tri") { 
  
    E <- c(Euv, Eb, Eg)
    names(E) <- c("Euv", "Eb", "Eg")
      
    x <- (sqrt(3)/2)*(Eg-Euv)
    y <- Eb-0.5*(Eg+Euv)
    
    deltaSo<-sqrt(x^2+y^2) 
    
    E <- c(Puv, Pb, Pg, Euv, Eb, Eg, x, y, deltaSo)
    E <- as.vector(E)
    names(E) <- c("P1", "P2", "P3", "E1", "E2", "E3", "x", "y", "deltaSo")
    }
    
    if(photo=="tetra") {  
    
    E <- c(Euv, Eb, Eg, Er)
    E <- as.vector(E)
    names(E) <- c("Euv", "Eb", "Eg", "Er")
      
    x <- (sqrt(3)*sqrt(2))/3 * (E[["Eg"]]-E[["Er"]]) 
    y <- E[["Euv"]] - (1/3)*(E[["Eb"]]+E[["Eg"]]+E[["Er"]])  
    z <- ((2*sqrt(2))/3) * ( ( 0.5*(E[["Eg"]]+E[["Er"]]) ) - E[["Eb"]] )
  
    deltaSo<-sqrt(x^2+y^2+z^2)
    
    E <- c(Puv, Pb, Pg, Pr, Euv, Eb, Eg, Er, x, y, z, deltaSo)
    E <- as.vector(E)
    names(E) <- c("P1", "P2", "P3", "P4", "E1", "E2", "E3", "E4", "x", "y", "z", "deltaSo")
    }
  return(E)
  }
  
  n.spectra<-ncol(R)-1
  R.list<-vector(length=n.spectra, "list")
  for (i in 1:n.spectra) {R.list[[i]]<-data.frame(R[,1],R[,i+1])}
  
  r<-sapply(X=R.list, FUN=internal, photo=photo, I=I,
            Rb=Rb, C=C, interpolate=interpolate, nm=nm)
  r<-as.data.frame(t(r))
  rownames(r)<-names(R)[2:ncol(R)]
  
  attr(r, "model name") <- "Colour hexagon model"
  attr(r, "number photor. types") <- photo
  
  return(r)
  
}


CTTKhexagon <- function (vnames=TRUE) {

  graphics::plot(x=0,y=0, pch=16, bty="n",yaxt="n",xaxt="n", col="white", ylim=c(-1.2,1.2), xlim=c(-1.2,1.2), asp=1, ann=FALSE)
  graphics::polygon(x=c(0,0.86660254,0.86660254,0,-0.86660254,-0.86660254,0),
          y=c(1,0.5,-0.5,-1,-0.5,0.5,1))
  if (vnames==TRUE) {
    graphics::text(x=0,y=1,labels=expression(E[2]),pos=3)
    graphics::text(x=-0.86660254,y=-0.5,labels=expression(E[1]), pos=1)
    graphics::text(x=0.86660254,y=-0.5, labels=expression(E[3]), pos=4)
  }
}

CTTKhexagon3D <- function (x,y,z, type = "p", s.col="red", radius=0.01, f.col="grey", vnames=TRUE) {
  
  requireNamespace("rgl")
  
  rgl::plot3d(x=x,y=y,z=z, col = s.col, type = type, add=F,
         xlab = "", ylab="", zlab="",
         box=F, axes=F, radius=radius, ylim=c(-1,1), xlim=c(-1,1), zlim=c(-1,1),
         aspect = T)
  
  
  #Photoreceptor vector vertices
  E4<-c(-0.8164966, -0.3333333,  0.4714045)
  E3<-c( 0.8164966, -0.3333333,  0.4714045)
  E2<-c( 0.0000000, -0.3333333, -0.9428090)
  E1<-c( 0.0000000,  1.0000000,  0.0000000)
  
  #Hexagonal trapezohedron vertices
  x.vertex<-c(0.0000000,  0.0000000,  0.8164966,
              0.8164966, -0.8164966, -0.8164966, 
              0.0000000,  0.0000000,  0.0000000,
              0.8164966, 0.8164966, -0.8164966, -0.8164966,  0.0000000)
  y.vertex<-c(1.0000000, -0.3333333, -0.3333333,
              0.3333333,  0.3333333, -0.3333333,
              0.3333333, -1.0000000,  0.6666667,
              0.6666667, -0.6666667,  0.6666667, -0.6666667, -0.6666667)
  z.vertex<-c(0.0000000, -0.9428090,  0.4714045,
              -0.4714045, -0.4714045,  0.4714045,
              0.9428090,  0.0000000, -0.9428090,
              0.4714045, -0.4714045,  0.4714045, -0.4714045,  0.9428090)
  
  rgl::plot3d(x=x.vertex[c(1,10)],y=y.vertex[c(1,10)],z=z.vertex[c(1,10)], col = f.col, type = "l", add=T, lwd=1)
  rgl::plot3d(x=x.vertex[c(1,9)],y=y.vertex[c(1,9)],z=z.vertex[c(1,9)], col = f.col, type = "l", add=T, lwd=1)
  rgl::plot3d(x=x.vertex[c(1,12)],y=y.vertex[c(1,12)],z=z.vertex[c(1,12)], col = f.col, type = "l", add=T, lwd=1)
  rgl::plot3d(x=x.vertex[c(5,9)],y=y.vertex[c(5,9)],z=z.vertex[c(5,9)], col = f.col, type = "l", add=T, lwd=1)
  rgl::plot3d(x=x.vertex[c(5,12)],y=y.vertex[c(5,12)],z=z.vertex[c(5,12)], col = f.col, type = "l", add=T, lwd=1)
  rgl::plot3d(x=x.vertex[c(6,12)],y=y.vertex[c(6,12)],z=z.vertex[c(6,12)], col = f.col, type = "l", add=T, lwd=1)
  rgl::plot3d(x=x.vertex[c(7,12)],y=y.vertex[c(7,12)],z=z.vertex[c(7,12)], col = f.col, type = "l", add=T, lwd=1)
  rgl::plot3d(x=x.vertex[c(7,10)],y=y.vertex[c(7,10)],z=z.vertex[c(7,10)], col = f.col, type = "l", add=T, lwd=1)
  rgl::plot3d(x=x.vertex[c(7,14)],y=y.vertex[c(7,14)],z=z.vertex[c(7,14)], col = f.col, type = "l", add=T, lwd=1)
  rgl::plot3d(x=x.vertex[c(3,10)],y=y.vertex[c(3,10)],z=z.vertex[c(3,10)], col = f.col, type = "l", add=T, lwd=1)
  rgl::plot3d(x=x.vertex[c(3,14)],y=y.vertex[c(3,14)],z=z.vertex[c(3,14)], col = f.col, type = "l", add=T, lwd=1)
  rgl::plot3d(x=x.vertex[c(6,14)],y=y.vertex[c(6,14)],z=z.vertex[c(6,14)], col = f.col, type = "l", add=T, lwd=1)
  rgl::plot3d(x=x.vertex[c(8,14)],y=y.vertex[c(8,14)],z=z.vertex[c(8,14)], col = f.col, type = "l", add=T, lwd=1)
  rgl::plot3d(x=x.vertex[c(8,13)],y=y.vertex[c(8,13)],z=z.vertex[c(8,13)], col = f.col, type = "l", add=T, lwd=1)
  rgl::plot3d(x=x.vertex[c(8,11)],y=y.vertex[c(8,11)],z=z.vertex[c(8,11)], col = f.col, type = "l", add=T, lwd=1)
  rgl::plot3d(x=x.vertex[c(5,13)],y=y.vertex[c(5,13)],z=z.vertex[c(5,13)], col = f.col, type = "l", add=T, lwd=1)
  rgl::plot3d(x=x.vertex[c(2,9)],y=y.vertex[c(2,9)],z=z.vertex[c(2,9)], col = f.col, type = "l", add=T, lwd=1)
  rgl::plot3d(x=x.vertex[c(2,11)],y=y.vertex[c(2,11)],z=z.vertex[c(2,11)], col = f.col, type = "l", add=T, lwd=1)
  rgl::plot3d(x=x.vertex[c(4,11)],y=y.vertex[c(4,11)],z=z.vertex[c(4,11)], col = f.col, type = "l", add=T, lwd=1)
  rgl::plot3d(x=x.vertex[c(4,9)],y=y.vertex[c(4,9)],z=z.vertex[c(4,9)], col = f.col, type = "l", add=T, lwd=1)
  rgl::plot3d(x=x.vertex[c(4,10)],y=y.vertex[c(4,10)],z=z.vertex[c(4,10)], col = f.col, type = "l", add=T, lwd=1)
  rgl::plot3d(x=x.vertex[c(2,13)],y=y.vertex[c(2,13)],z=z.vertex[c(2,13)], col = f.col, type = "l", add=T, lwd=1)
  rgl::plot3d(x=x.vertex[c(6,13)],y=y.vertex[c(6,13)],z=z.vertex[c(6,13)], col = f.col, type = "l", add=T, lwd=1)
  rgl::plot3d(x=x.vertex[c(3,11)],y=y.vertex[c(3,11)],z=z.vertex[c(3,11)], col = f.col, type = "l", add=T, lwd=1)
  
  if (vnames==TRUE) {
    rgl::text3d(x=E1[[1]],y=E1[[2]],z=E1[[3]], texts="E1", cex=.75, adj=c(1,1))
    rgl::text3d(x=E2[[1]],y=E2[[2]],z=E2[[3]], texts="E2", cex=.75, adj=c(1,1))
    rgl::text3d(x=E3[[1]],y=E3[[2]],z=E3[[3]], texts="E3", cex=.75, adj=c(1,1))
    rgl::text3d(x=E4[[1]],y=E4[[2]],z=E4[[3]], texts="E4", cex=.75, adj=c(1,1))
  }
  
}



EMtriangle <- function (ylim=c(-0.8,0.8),
                        xlim=c(-0.8,0.8),
                        vnames=TRUE) {
  
  graphics::plot(x=0,y=0, pch=16, bty="n",yaxt="n",xaxt="n", col="white",
       ylim=ylim, xlim=xlim, asp=1, ann=FALSE)
  graphics::polygon(
    
    x=c((2/3)*((sqrt(3)/2)*(0-1)),
              (2/3)*((sqrt(3)/2)*(0-0)),
              (2/3)*((sqrt(3)/2)*(1-0))
      ),
    
    y=c((2/3)*(0-0.5*(1+0)),
      (2/3)*(1-0.5*(0+0)),
      (2/3)*(0-0.5*(0+1))
      )
  )
  if(vnames==TRUE) {
    graphics::text(x=(2/3)*((sqrt(3)/2)*(0-1)),
                 y=(2/3)*(0-0.5*(1+0)),
                 labels = "u",
                 pos=1)
    graphics::text(x=(2/3)*((sqrt(3)/2)*(0-0)),
                 y=(2/3)*(1-0.5*(0+0)),
                 labels = "s",
                 pos=3)
    graphics::text(x=((2/3)*((sqrt(3)/2)*(1-0))),
                 y=(2/3)*(0-0.5*(0+1)),
                 labels = "m",
                 pos=4)
  }
}


EMtetrahedron <- function (x,y,z, type = "p", s.col="red", radius=0.01, f.col="black", vnames=TRUE) {
  
  requireNamespace("rgl")
  rgl::rgl.viewpoint(  zoom = .75 )
  rgl::plot3d(x=x,y=y,z=z, col = s.col, type = type, add=F,
         xlab = "", ylab="", zlab="",
         box=F, axes=F, radius=radius, ylim=c(-0.75,0.75), xlim=c(-0.75,0.75), zlim=c(-0.75,0.75),
         aspect = T, mar=c(1,1,1,1))
  
  smu <- function (s,m,u,l) {
    
    x <- ((1-2*s-m-u)/2)*sqrt(3/2)
    y <- (-1+3*m+u)/(2*sqrt(2))
    z <- u - 1/4
    r<-c(x,y,z)
    names(r)<-c("x","y","z")
    return(r)
  }
  
  rgl::plot3d(x=c(smu(1,0,0)[["x"]],smu(0,0,0)[["x"]]),
         y=c(smu(1,0,0)[["y"]],smu(0,0,0)[["y"]]),
         z=c(smu(1,0,0)[["z"]],smu(0,0,0)[["z"]]),
         col = f.col, type = "l", add=T, lwd=1)
  rgl::plot3d(x=c(smu(0,0,0)[["x"]],smu(0,0,1)[["x"]]),
         y=c(smu(0,0,0)[["y"]],smu(0,0,1)[["y"]]),
         z=c(smu(0,0,0)[["z"]],smu(0,0,1)[["z"]]),
         col = f.col, type = "l", add=T, lwd=1)
  rgl::plot3d(x=c(smu(1,0,0)[["x"]],smu(0,0,1)[["x"]]),
         y=c(smu(1,0,0)[["y"]],smu(0,0,1)[["y"]]),
         z=c(smu(1,0,0)[["z"]],smu(0,0,1)[["z"]]),
         col = f.col, type = "l", add=T, lwd=1)
  rgl::plot3d(x=c(smu(0,1,0)[["x"]],smu(0,0,1)[["x"]]),
         y=c(smu(0,1,0)[["y"]],smu(0,0,1)[["y"]]),
         z=c(smu(0,1,0)[["z"]],smu(0,0,1)[["z"]]),
         col = f.col, type = "l", add=T, lwd=1)
  rgl::plot3d(x=c(smu(0,1,0)[["x"]],smu(0,0,0)[["x"]]),
         y=c(smu(0,1,0)[["y"]],smu(0,0,0)[["y"]]),
         z=c(smu(0,1,0)[["z"]],smu(0,0,0)[["z"]]),
         col = f.col, type = "l", add=T, lwd=1)
  rgl::plot3d(x=c(smu(0,1,0)[["x"]],smu(1,0,0)[["x"]]),
         y=c(smu(0,1,0)[["y"]],smu(1,0,0)[["y"]]),
         z=c(smu(0,1,0)[["z"]],smu(1,0,0)[["z"]]),
         col = f.col, type = "l", add=T, lwd=1)
  
  if (vnames==TRUE) {
    rgl::text3d(x=smu(1,0,0)[["x"]],y=smu(1,0,0)[["y"]],z=smu(1,0,0)[["z"]], texts="u", cex=1, adj=c(0,0))
    rgl::text3d(x=smu(0,1,0)[["x"]],y=smu(0,1,0)[["y"]],z=smu(0,1,0)[["z"]], texts="s", cex=1, adj=c(1,1))
    rgl::text3d(x=smu(0,0,1)[["x"]],y=smu(0,0,1)[["y"]],z=smu(0,0,1)[["z"]], texts="m", cex=1, adj=c(1,1))
    rgl::text3d(x=smu(0,0,0)[["x"]],y=smu(0,0,0)[["y"]],z=smu(0,0,0)[["z"]], texts="l", cex=1, adj=c(1,1))
  }
}