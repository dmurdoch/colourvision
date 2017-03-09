spec.denoise <- function(specfiles, spar = 0.7, ...){
  result<-NULL
  name<-names(specfiles)
  for (i in 2:ncol(specfiles)){
    temp<-smooth.spline(x = specfiles[,1], y = specfiles[,i], spar = spar, ...)$y
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
  A <- 69.7
  B <- 28
  b <- 0.922
  C <- -14.9
  c <- 1.104
  D <- 0.674
  
  nphotor<-length(lambda.max)
  r <- matrix(ncol=nphotor, nrow = length(lambda))
  
  for (k in 1:nphotor) {
  
  a <- 0.8795 + 0.0459*exp(-1*(((lambda.max[[k]]-300)^2)/11940))
  
    for(i in 1:length(lambda)) {
      
      x <- lambda.max[[k]]/lambda[[i]]
      r[i,k] <- 1/( exp(A*(a-x)) + exp(B*(b-x)) + exp(C*(c-x)) + D)
    }
  
  #beta band
  if (beta.band == TRUE) {
    beta <- vector(length = length(lambda))
    lambda.max.beta = 189 + 0.315*lambda.max[[k]]
    A.beta = 0.26
    for(j in 1:length(lambda)) {
      b.beta = -40.5 + 0.195*lambda.max[[k]]
      beta[[j]] <- A.beta*exp(-1*((lambda[[j]]-lambda.max.beta)/b.beta)^2)
    }
    r[,k] <- r[,k] + beta
    }
  }
  
  r <- cbind(lambda, r)
  r<-data.frame(r)
  
  colnames(r)<-c("Wavelength", paste("lambda.max", lambda.max, sep=""))
  return(r)

}


Q<-function (R, I, C, interpolate, nm) 
{
  if (interpolate == TRUE) {
    I <- data.frame(nm, approx(x = I[, 1], y = I[, 2], xout = nm, 
                               method = "linear")$y)
    R <- data.frame(nm, approx(x = R[, 1], y = R[, 2], xout = nm, 
                               method = "linear")$y)
    C <- data.frame(nm, approx(x = C[, 1], y = C[, 2], xout = nm, 
                               method = "linear")$y)
  }
  if (interpolate == FALSE) {
    test <- c(length(R[, 1]) == length(I[, 1]), length(R[, 
                                                         1]) == length(C[, 1]), length(C[, 1]) == length(I[, 
                                                                                                           1]))
    ifelse(any(test == FALSE), yes = stop("Uneven number of rows. Use 'interpolate=TRUE'.", 
                                          call. = FALSE), no = "")
    ifelse(any(c(R[, 1] == I[, 1], R[, 1] == C[, 1], I[, 
                                                       1] == C[, 1])) == FALSE, yes = stop("Different wavelenght values of model parameters. Use 'interpolate=TRUE'.", 
                                                                                           call. = FALSE), no = "")
    a <- R[, 1]
    b <- R[2:length(R[, 1]), 1]
    ifelse(sd(b - a[1:length(a) - 1], na.rm = T) != 0, yes = stop("Uneven wavelenght intervals. Use 'interpolate=TRUE'.", 
                                                                  call. = FALSE), no = "")
  }
  int<-abs(R[2, 1]-R[1, 1])
  r <- sum(I[, 2] * R[, 2] * C[, 2])*int
  return(r)
}


Qr <- function(R, I, Rb, C, interpolate, nm) {
  Qr <- Q(I=I,R=R,C=C,interpolate=interpolate,nm=nm)
  QEr <- Q(I=I,R=Rb,C=C,interpolate=interpolate,nm=nm)
  r<-Qr/QEr
  return(r)
}

colour_space<-function(n, length, q) {
  v1<-vector(length=n-1)
  v1[[n-1]]<- -1/(n-1)
  for (i in (length(v1)-1):1) {
    v1[[i]]<--sqrt(1-sum(v1[(i+1):length(v1)]^2))/(i)
  }
  
  #vector rotation
  v<-matrix(data=1, ncol=n-1, nrow=n-1)
  v[1,1]<--1
  for (i in 2:ncol(v)) {
    v[i,i]<-v[i-1,i-1]-1
  }
  
  for (i in 1:(ncol(v)-1)) {
    for (k in (i+1):nrow(v)) {
      v[k,i]<-0
    }
  }
  
  for (i in 1:nrow(v)) {
    v[i,]<-v[i,]*v1
  }
  v<-rbind(v1,v)
  v<-v*length
  
  col.names<-vector(length=ncol(v))
  for (i in 1:ncol(v)) {
    col.names[[i]]<-paste("X",i,sep="")
  }
  row.names<-vector(length=nrow(v))
  for (i in 1:nrow(v)) {
    row.names[[i]]<-paste("v",i,sep="")
  }
  rownames(v)<-row.names
  colnames(v)<-col.names
  
  X<-vector(length=(n-1))
  for (i in 1:(n-1)) {
    X[[n-i]]<-v[n+1-i,n-i] * (q[[n+1-i]] - sum(q[1:(n-i)])/(n-i) ) 
  }
  names(X)<-colnames(v)
  
  r<-vector("list", length=2)
  r[[1]]<-X
  r[[2]]<-v
  names(r)<-c("coordinates", "vector_matrix")
  class(r)<-"colourvision"
  return(r)
}

EMmodel <- function (photo=ncol(C)-1,
                     R,
                     I,
                     Rb,
                     C,
                     interpolate=TRUE,
                     nm=seq(300,700,1))
{
  photo1<-photo
  if(photo=="di"){photo1<-2}
  if(photo=="tri"){photo1<-3}
  if(photo=="tetra"){photo1<-4}
  if(photo=="penta"){photo1<-5}
  nphoto=ncol(C)-1
  
  ifelse(photo1 != nphoto, yes = warning("Argument 'C' has a number of sensitivity curves different than argument 'photo'.", 
                                         call. = FALSE), no = "")
  ifelse(any(ncol(I) > 2), yes = warning("'I' argument with more than two columns. Only the first two will be used.", 
                                         call. = FALSE), no = "")
  ifelse(any(ncol(Rb) > 2), yes = warning("'Rb' argument with more than two columns. Only the first two will be used.", 
                                          call. = FALSE), no = "")
  maxR <- apply(data.frame(R[, 2:ncol(R)]), 2, max)
  ifelse(any(maxR <= 1) && max(Rb[, 2]) > 1 || any(maxR > 1) && max(Rb[, 2]) <= 1,
         yes = warning("There seems to be a problem with input files. 'R' and 'Rb' must be in the same scale. Both must be either in percentage (0-100%) or proportion (0-1).", 
                       call. = FALSE), no = "")
  
  internal<- function (photo,
                       R,
                       I,
                       Rb,
                       C,
                       interpolate,
                       nm)
  {
    
    S<-vector(length=photo1)
    for (i in 1:photo1) {
      S[[i]]<-Qr(I=I, R=R, Rb=Rb, C=C[,c(1,1+i)], interpolate=interpolate, nm=nm)
    }
    S.log<-log(S)
    E<-S.log/sum(S.log)

    if(photo1==2) {
      x <- (3/4)*(E[[2]]-E[[1]])
      
      deltaSo<-abs(x)
      
      r<-c(S, E, x, deltaSo)
      r<-as.vector(r)
    }
    
    if(photo1==3) {
      
      x<-(-3*sqrt(3)/8)*(E[[2]]-E[[1]])
      y<-(3/4)*(E[[3]]-(E[[1]]+E[[2]])/2)
      
      deltaSo<-sqrt(x^2 + y^2)
      
      r<-c(S, E, x, y, deltaSo)
      r<-as.vector(r)
    }
    
    if(photo1==4) {
      x <- ((1-2*E[[2]]-E[[3]]-E[[1]])/2)*sqrt(3/2)
      y <- (-1+3*E[[3]]+E[[1]])/(2*sqrt(2))
      z <- E[[1]] - (1/4)

      deltaSo<-sqrt(x^2 + y^2 + z^2)
      
      r<-c(S, E, x, y, z, deltaSo)
      r<-as.vector(r)
    }
    
    if(photo1==3 || photo1>4) {
      
      r<-colour_space(n=photo1, length=0.75, q=E)      
      r<-c(S, E, r$coordinates, sqrt(sum(r$coordinates^2)))
      
    }
    return(r)
    
  }
  
  n.spectra<-ncol(R)-1
  R.list<-vector(length=n.spectra, "list")
  for (i in 1:n.spectra) {R.list[[i]]<-data.frame(R[,1],R[,i+1])}
  
  r<-sapply(X=R.list, FUN=internal, photo=photo, I=I,
            Rb=Rb, C=C, interpolate=interpolate, nm=nm)
  r<-as.data.frame(t(r))

  #names
  namesQr<-vector(length=photo1)
  namesE<-vector(length=photo1)
  namesX<-vector(length=photo1-1)
  for (i in 1:photo1) {
    namesQr[[i]]<-paste("Qr", i, sep="")
    namesE[[i]]<-paste("E", i, sep="")
  }
  for (i in 1:(photo1-1)) {
    namesX[[i]]<-paste("X", i, sep="")
  }
  
  r.names<-c(namesQr, namesE, namesX, "deltaS")
  
  colnames(r)<-r.names
  rownames(r)<-names(R)[2:ncol(R)]
  
  test <- log(r[,namesQr])
  ifelse (any(test<0),
            yes=warning("Log-transformation might be generating photoreceptor outputs < 0.", call.=FALSE), no="")
  ifelse (any(r[,"deltaS"]>0.75),
          yes=warning("deltaS > 0.75. Log-transformation might be generating photoreceptor outputs < 0.", call.=FALSE), no="")
  
  class(r)<-c("colourvision", "data.frame")
  attr(r, "model_name") <- "Endler and Mielke model"
  attr(r, "n_photor_types") <- photo1
  attr(r, "Rb") <- Rb
  attr(r, "I") <- I
  attr(r, "C") <- C
  attr(r, "Interpolate") <- interpolate
  attr(r, "nm") <- nm
  
  return(r)
  
}


noise_e<-function (noise, e, v, n) 
{
  dependent<-FALSE
  quantum<-NULL
  if (noise == TRUE) {
    if(dependent==FALSE) {
      r<-e
    }
    if(dependent==TRUE) {
      r<-sqrt ( (e^2)+(1/quantum) )
    }
  }
  if (noise == FALSE) {
    if(dependent==FALSE) {
      r <- v/sqrt(n)
    }
    if(dependent==TRUE) {
      r <- sqrt( ((v^2)/n)+(1/quantum) )
    }
  }
  return(r)
}


RNLmodel <- function (model = c("linear", "log"), photo=ncol(C)-1,
                      R1, R2=Rb, Rb, I, C, noise = FALSE, v=NA, n=NA, e=NA,
                      interpolate = TRUE, nm = seq(300, 700, 1)) 
{
  
  dependent <- FALSE
  nphoto = ncol(C) - 1
  photo1<-photo
  if (photo=="di") {photo1<-2}
  if (photo=="tri") {photo1<-3}
  if (photo=="tetra") {photo1<-4}
  if (photo=="penta") {photo1<-5}
  
  #warnings
  ifelse(photo1 != nphoto, yes = warning("Argument 'C' has a number of sensitivity curves different than argument 'photo'.", 
                                                   call. = FALSE), no = "")
  ifelse(photo1 != length(e) && noise == T, 
         yes = warning("Argument 'e' has a number of parameters different than 'photo'.", 
                       call. = FALSE), no = "")
  ifelse(photo1 != length(n) && noise == F, 
         yes = warning("Argument 'n' has a number of parameters different than 'photo'.", 
                       call. = FALSE), no = "")
  ifelse(any(ncol(I) > 2), yes = warning("'I' argument with more than two columns. Only the first two will be used.", 
                                         call. = FALSE), no = "")
  ifelse(any(ncol(Rb) > 2), yes = warning("'Rb' argument with more than two columns. Only the first two will be used.", 
                                          call. = FALSE), no = "")
  ifelse(any(ncol(R2) > 2), yes = warning("'R2' argument with more than two columns. Only the first two will be used.", 
                                          call. = FALSE), no = "")
  maxR <- apply(data.frame(R1[, 2:ncol(R1)]), 2, max)
  ifelse(any(maxR <= 1) && max(Rb[, 2]) > 1 || any(maxR > 1) && max(Rb[, 2]) <= 1,
         yes = warning("There seems to be a problem with input files. 'R1' and 'Rb' must be in the same scale. Both must be either in percentage (0-100%) or proportion (0-1).", 
                                            call. = FALSE), no = "")
  ifelse(any(maxR <= 1) && max(R2[, 2]) > 1 || any(maxR > 1) && max(R2[, 2]) <= 1,
         yes = warning("There seems to be a problem with input files. 'R1' and 'R2' must be in the same scale. Both must be either in percentage (0-100%) or proportion (0-1).", 
                       call. = FALSE), no = "")
  
  #internal function to be used with 'apply'
  internal <- function(model, photo, R1, R2=Rb, Rb, I, C, noise, 
                       v, n, e, interpolate, nm) {
    
    #relative photon catches
    S1<-vector(length=photo1)
    S2<-vector(length=photo1)
    for (i in 1:photo1) {
      S1[[i]] <- Qr(I = I, R = R1, Rb = Rb, C = C[, c(1, i+1)], interpolate = interpolate, nm = nm)
      S2[[i]] <- Qr(I = I, R = R2, Rb = Rb, C = C[, c(1, i+1)], interpolate = interpolate, nm = nm)
    }
    
    S1.Qr<-S1
    S2.Qr<-S2
    
    if (model == "log") {
      S1 <- log10(S1)
      S2 <- log10(S2)
    }
    
    #noise
    noise_values<-vector(length=photo1)
    if (dependent == FALSE) {
      for (i in 1:photo1) {
        noise_values[[i]]<-noise_e(noise = noise, e = e[[i]], v = v, n = n[[i]])
      }
    }
    
    #colour loci  
    if (photo1 == 4) {
      A1 <- sqrt(1/((noise_values[[3]]^2) + (noise_values[[4]]^2)))
      B1 <- sqrt(((noise_values[[3]]^2) + (noise_values[[4]]^2))/((noise_values[[2]] *noise_values[[3]])^2 + (noise_values[[2]] * noise_values[[4]])^2 + (noise_values[[3]] * noise_values[[4]])^2))
      C1 <- sqrt(((noise_values[[2]] * noise_values[[3]])^2 + (noise_values[[2]] * noise_values[[4]])^2 + 
                    (noise_values[[3]] * noise_values[[4]])^2)/((noise_values[[2]] * noise_values[[3]] * 
                                               noise_values[[4]])^2 + (noise_values[[1]] * noise_values[[3]] * noise_values[[4]])^2 + 
                                              (noise_values[[1]] * noise_values[[2]] * noise_values[[4]])^2 + (noise_values[[1]] * 
                                                                                   noise_values[[2]] * noise_values[[3]])^2))
      a1 <- (noise_values[[3]]^2)/((noise_values[[3]]^2) + (noise_values[[4]]^2))
      b1 <- (noise_values[[4]]^2)/((noise_values[[3]]^2) + (noise_values[[4]]^2))
      a2 <- ((noise_values[[2]] * noise_values[[3]])^2)/((noise_values[[2]] * noise_values[[3]])^2 + 
                                       (noise_values[[2]] * noise_values[[4]])^2 + (noise_values[[3]] * noise_values[[4]])^2)
      b2 <- ((noise_values[[2]] * noise_values[[4]])^2)/((noise_values[[2]] * noise_values[[3]])^2 + 
                                       (noise_values[[2]] * noise_values[[4]])^2 + (noise_values[[3]] * noise_values[[4]])^2)
      c2 <- ((noise_values[[3]] * noise_values[[4]])^2)/((noise_values[[2]] * noise_values[[3]])^2 + 
                                       (noise_values[[2]] * noise_values[[4]])^2 + (noise_values[[3]] * noise_values[[4]])^2)
      x1 <- A1 * (S1[[4]] - S1[[3]])
      y1 <- B1 * (S1[[2]] - (a1 * S1[[4]] + b1 * S1[[3]]))
      z1 <- C1 * (S1[[1]] - (a2 * S1[[4]] + b2 * S1[[3]] + c2 * S1[[2]]))

      x2 <- A1 * (S2[[4]] - S2[[3]])
      y2 <- B1 * (S2[[2]] - (a1 * S2[[4]] + b1 * S2[[3]]))
      z2 <- C1 * (S2[[1]] - (a2 * S2[[4]] + b2 * S2[[3]] + c2 * S2[[2]]))
      
      delta_e <- sqrt((x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2)
      
      r <- c(noise_values, S1.Qr, S2.Qr, S1, S2, x1, y1, z1, x2, y2, z2, delta_e)
      r <- as.vector(r)
    }
    
    if (photo1 == 3) {
      A1 <- sqrt(1/(noise_values[[2]]^2 + noise_values[[3]]^2))
      B1 <- sqrt((noise_values[[2]]^2 + noise_values[[3]]^2)/(((noise_values[[1]] * noise_values[[2]])^2) + 
                                            ((noise_values[[1]] * noise_values[[3]])^2) + ((noise_values[[2]] * noise_values[[3]])^2)))
      a1 <- noise_values[[2]]^2/(noise_values[[2]]^2 + noise_values[[3]]^2)
      b1 <- noise_values[[3]]^2/(noise_values[[2]]^2 + noise_values[[3]]^2)
      X1 <- A1 * (S1[[3]] - S1[[2]])
      Y1 <- B1 * (S1[[1]] - (a1 * S1[[3]] + b1 * S1[[2]]))
      X2 <- A1 * (S2[[3]] - S2[[2]])
      Y2 <- B1 * (S2[[1]] - (a1 * S2[[3]] + b1 * S2[[2]]))
      delta_e <- sqrt((X1-X2)^2 + (Y1-Y2)^2)
      r <- c(noise_values, S1.Qr, S2.Qr, S1, S2, X1, Y1, X2, Y2, delta_e)
      r <- as.vector(r)
    }
    
    if (photo1 == 2) {
      X1 <- sqrt(1/(noise_values[[1]]^2 + noise_values[[2]]^2)) * (S1[[2]] - S1[[1]])
      X2 <- sqrt(1/(noise_values[[1]]^2 + noise_values[[2]]^2)) * (S2[[2]] - S2[[1]])
      delta_e <- sqrt((X1-X2)^2)
      r <- c(noise_values, S1.Qr, S2.Qr, S1, S2, X1, X2, delta_e)
      r <- as.vector(r)
    }
    
    if (photo1>4) {
    noise_numerator<-combn(noise_values,(photo1-2))
      for (i in 2:nrow(noise_numerator)) {
        noise_numerator[1,]<-noise_numerator[1,]*noise_numerator[i,]
      }
    noise_numerator<-noise_numerator[1,]
    noise_numerator<-noise_numerator^2
    
    #q_numerator
    q1_numerator<-combn(S1,2)
    q2_numerator<-combn(S2,2)
    q_numerator<-q1_numerator-q2_numerator
    q_numerator<-q_numerator[1,]-q_numerator[2,]
    q_numerator<-q_numerator^2
    
    #numerator
    numerator<-sum(noise_numerator*q_numerator[photo1:1])
    
    #noise denominatior
    noise_denominator<-combn(noise_values,(photo1-1))
    for (i in 2:nrow(noise_denominator)) {
      noise_denominator[1,]<-noise_denominator[1,]*noise_denominator[i,]
    }
    noise_denominator<-sum(noise_denominator[1,]^2)
    
    #deltaS
    deltaS<-sqrt(numerator/noise_denominator)
    r<-c(noise_values, S1.Qr, S2.Qr, S1, S2, deltaS)
    }
    
    return(r)
  }
  
  #apply function for several spectra
  n.spectra <- ncol(R1) - 1
  R1.list <- vector(length = n.spectra, "list")
  for (i in 1:n.spectra) {
    R1.list[[i]] <- data.frame(R1[, 1], R1[, i + 1])
  }
  r <- sapply(X = R1.list, FUN = internal, model = model, photo = photo, 
              I = I, R2=R2, Rb = Rb, C = C, noise = noise,
              v = v, n = n, e = e, interpolate = interpolate, nm = nm)
  r <- as.data.frame(t(r))
  
  e_names<-vector(length=photo1)
  Qr_1names<-vector(length=photo1)
  Qr_2names<-vector(length=photo1)
  E_1names<-vector(length=photo1)
  E_2names<-vector(length=photo1)
  for (i in 1:photo1) {
    e_names[[i]]<-paste("e", i, sep="")
    Qr_1names[[i]]<-paste("Qr", i, "_R1", sep="")
    Qr_2names[[i]]<-paste("Qr", i, "_R2", sep="")
    E_1names[[i]]<-paste("E", i, "_R1", sep="")
    E_2names[[i]]<-paste("E", i, "_R2", sep="")
  }
  if (photo1<=4) {
    namesX1<-vector(length=(photo1-1))
    namesX2<-vector(length=(photo1-1))
    for (i in 1:(photo1-1)) {
      namesX1[[i]]<-paste("X", i, "_R1", sep="")
      namesX2[[i]]<-paste("X", i, "_R2", sep="")
    }
    colnames(r) <- c(e_names, Qr_1names, Qr_2names, E_1names, E_2names, namesX1, namesX2, "deltaS")
  }
  if (photo1>4) {
    colnames(r) <- c(e_names, Qr_1names, Qr_2names, E_1names, E_2names, "deltaS")
  }
    rownames(r) <- names(R1)[2:ncol(R1)]
  
  test1<- r[,c(E_1names, E_2names)]
  ifelse(any(test1 < 0), yes = warning("Photoreceptor output < 0.", 
                                        call. = FALSE), no = "")

  class(r)<-c("colourvision", "data.frame")
  attr(r, "model_name") <- "Receptor noise limited model"
  attr(r, "model_input_output") <- model
  attr(r, "n_photor_types") <- photo1
  attr(r, "R2") <- R2
  attr(r, "Rb") <- Rb
  attr(r, "I") <- I
  attr(r, "C") <- C
  attr(r, "noise calculated") <- noise
  attr(r, "v") <- v
  attr(r, "n") <- n
  attr(r, "Interpolate") <- interpolate
  attr(r, "nm") <- nm
  
    return(r)
}



CTTKmodel <- function (photo=ncol(C)-1,
                  R,
                  I,
                  Rb,
                  C,
                  interpolate=TRUE,
                  nm=seq(300,700,1))
{
  photo1<-photo
  if(photo=="di"){photo1<-2}
  if(photo=="tri"){photo1<-3}
  if(photo=="tetra"){photo1<-4}
  if(photo=="penta"){photo1<-5}
  nphoto=ncol(C)-1
  
  ifelse(photo1 != nphoto, yes = warning("Argument 'C' has a number of sensitivity curves different than argument 'photo'.", 
                                         call. = FALSE), no = "")
  ifelse(any(ncol(I) > 2), yes = warning("'I' argument with more than two columns. Only the first two will be used.", 
                                         call. = FALSE), no = "")
  ifelse(any(ncol(Rb) > 2), yes = warning("'Rb' argument with more than two columns. Only the first two will be used.", 
                                          call. = FALSE), no = "")
  maxR <- apply(data.frame(R[, 2:ncol(R)]), 2, max)
  ifelse(any(maxR <= 1) && max(Rb[, 2]) > 1 || any(maxR > 1) && max(Rb[, 2]) <= 1,
         yes = warning("There seems to be a problem with input files. 'R' and 'Rb' must be in the same scale. Both must be either in percentage (0-100%) or proportion (0-1).", 
                       call. = FALSE), no = "")

  internal <- function (photo,
                        R,
                        I,
                        Rb,
                        C,
                        interpolate,
                        nm) {
  
    P<-vector(length=photo1)
    for (i in 1:length(P)) {
      P[[i]]<-Qr(I=I, R=R, Rb=Rb, C=C[,c(1,i+1)], interpolate=interpolate, nm=nm)
    }

    E<-P/(1+P)

    if (photo1==2) {
      x <- E[[2]]-E[[1]]
      deltaSo<-abs(x)
      r <- c(P, E, x, deltaSo)
    }
    
    if (photo1==3) {
      x <- (sqrt(3)/2)*(E[[3]]-E[[1]])
      y <- E[[2]]-0.5*(E[[3]]+E[[1]])
      deltaSo<-sqrt(x^2+y^2) 
      r <- c(P, E, x, y, deltaSo)
    }
    
    if(photo1==4) {
      
      x <- (sqrt(3)*sqrt(2))/3 * (E[[3]]-E[[4]]) 
      y <- E[[1]] - (1/3)*(E[[2]]+E[[3]]+E[[4]])  
      z <- ((2*sqrt(2))/3) * ( ( 0.5*(E[[3]]+E[[4]]) ) - E[[2]] )
  
      deltaSo<-sqrt(x^2+y^2+z^2)
    
      r <- c(P, E, x, y, z, deltaSo)
    }

    if(photo1>4) {
      
      r<-colour_space(n=photo1, length=1, q=E)      
      r<-c(P, E, r$coordinates, sqrt(sum(r$coordinates^2)))
      
    }
    return(r)
  }
  
  n.spectra<-ncol(R)-1
  R.list<-vector(length=n.spectra, "list")
  for (i in 1:n.spectra) {R.list[[i]]<-data.frame(R[,1],R[,i+1])}
  
  r<-sapply(X=R.list, FUN=internal, photo=photo, I=I,
            Rb=Rb, C=C, interpolate=interpolate, nm=nm)
  r<-as.data.frame(t(r))
  
  #names
  namesQr<-vector(length=photo1)
  namesE<-vector(length=photo1)
  namesX<-vector(length=photo1-1)
  for (i in 1:photo1) {
    namesQr[[i]]<-paste("Qr", i, sep="")
    namesE[[i]]<-paste("E", i, sep="")
  }
  for (i in 1:(photo1-1)) {
    namesX[[i]]<-paste("X", i, sep="")
  }
  
  r.names<-c(namesQr, namesE, namesX, "deltaS")
  
  colnames(r)<-r.names
  rownames(r)<-names(R)[2:ncol(R)]
  
  class(r)<-c("colourvision", "data.frame")
  attr(r, "model_name") <- "Colour hexagon model"
  attr(r, "n_photor_types") <- photo1
  attr(r, "Rb") <- Rb
  attr(r, "I") <- I
  attr(r, "C") <- C
  attr(r, "Interpolate") <- interpolate
  attr(r, "nm") <- nm
  return(r)
  
}


CTTKhexagon <- function (x, y, vnames=c(expression(E[1]),expression(E[2]),expression(E[3])),
                         pch=16, bty="n",
                         yaxt="n",xaxt="n", col="black",
                         ylim=c(-1.2,1.2), xlim=c(-1.2,1.2),
                         asp=1, ann=FALSE, ...) {

  graphics::plot(x=x,y=y, pch=pch, bty=bty,yaxt=yaxt,xaxt=xaxt, col=col, ylim=ylim, xlim=xlim, asp=asp, ann=ann, ...)
  graphics::polygon(x=c(0,0.86660254,0.86660254,0,-0.86660254,-0.86660254,0),
          y=c(1,0.5,-0.5,-1,-0.5,0.5,1))
  graphics::text(x=0,y=1,labels=vnames[[2]],pos=3)
  graphics::text(x=-0.86660254,y=-0.5,labels=vnames[[1]], pos=1)
  graphics::text(x=0.86660254,y=-0.5, labels=vnames[[3]], pos=4)
}

CTTKhexagon3D<- function (x, y, z, s.col = "red", f.col = "black", vnames = c("E1","E2","E3","E4"),
                          type = "p", radius = 0.01, 
                          add = F, xlab = "", ylab = "", zlab = "",
                          box = F, axes = F, ylim = c(-1, 1), xlim = c(-1, 1),
                          zlim = c(-1,1), aspect = T, ...) 
{
  requireNamespace("rgl")
  rgl::plot3d(x = x, y = y, z = z, col = s.col, type = type, 
              add = add, xlab = xlab, ylab = ylab, zlab = zlab, box = box, axes = axes, 
              radius = radius, ylim = ylim, xlim = xlim, zlim = zlim, aspect = aspect, ...)
  E4 <- c(-0.8164966, -0.3333333, 0.4714045)
  E3 <- c(0.8164966, -0.3333333, 0.4714045)
  E2 <- c(0, -0.3333333, -0.942809)
  E1 <- c(0, 1, 0)
  x.vertex <- c(0, 0, 0.8164966, 0.8164966, -0.8164966, -0.8164966, 
                0, 0, 0, 0.8164966, 0.8164966, -0.8164966, -0.8164966, 
                0)
  y.vertex <- c(1, -0.3333333, -0.3333333, 0.3333333, 0.3333333, 
                -0.3333333, 0.3333333, -1, 0.6666667, 0.6666667, -0.6666667, 
                0.6666667, -0.6666667, -0.6666667)
  z.vertex <- c(0, -0.942809, 0.4714045, -0.4714045, -0.4714045, 
                0.4714045, 0.942809, 0, -0.942809, 0.4714045, -0.4714045, 
                0.4714045, -0.4714045, 0.942809)
  rgl::plot3d(x = x.vertex[c(1, 10)], y = y.vertex[c(1, 10)], 
              z = z.vertex[c(1, 10)], col = f.col, type = "l", add = T, 
              lwd = 1)
  rgl::plot3d(x = x.vertex[c(1, 9)], y = y.vertex[c(1, 9)], 
              z = z.vertex[c(1, 9)], col = f.col, type = "l", add = T, 
              lwd = 1)
  rgl::plot3d(x = x.vertex[c(1, 12)], y = y.vertex[c(1, 12)], 
              z = z.vertex[c(1, 12)], col = f.col, type = "l", add = T, 
              lwd = 1)
  rgl::plot3d(x = x.vertex[c(5, 9)], y = y.vertex[c(5, 9)], 
              z = z.vertex[c(5, 9)], col = f.col, type = "l", add = T, 
              lwd = 1)
  rgl::plot3d(x = x.vertex[c(5, 12)], y = y.vertex[c(5, 12)], 
              z = z.vertex[c(5, 12)], col = f.col, type = "l", add = T, 
              lwd = 1)
  rgl::plot3d(x = x.vertex[c(6, 12)], y = y.vertex[c(6, 12)], 
              z = z.vertex[c(6, 12)], col = f.col, type = "l", add = T, 
              lwd = 1)
  rgl::plot3d(x = x.vertex[c(7, 12)], y = y.vertex[c(7, 12)], 
              z = z.vertex[c(7, 12)], col = f.col, type = "l", add = T, 
              lwd = 1)
  rgl::plot3d(x = x.vertex[c(7, 10)], y = y.vertex[c(7, 10)], 
              z = z.vertex[c(7, 10)], col = f.col, type = "l", add = T, 
              lwd = 1)
  rgl::plot3d(x = x.vertex[c(7, 14)], y = y.vertex[c(7, 14)], 
              z = z.vertex[c(7, 14)], col = f.col, type = "l", add = T, 
              lwd = 1)
  rgl::plot3d(x = x.vertex[c(3, 10)], y = y.vertex[c(3, 10)], 
              z = z.vertex[c(3, 10)], col = f.col, type = "l", add = T, 
              lwd = 1)
  rgl::plot3d(x = x.vertex[c(3, 14)], y = y.vertex[c(3, 14)], 
              z = z.vertex[c(3, 14)], col = f.col, type = "l", add = T, 
              lwd = 1)
  rgl::plot3d(x = x.vertex[c(6, 14)], y = y.vertex[c(6, 14)], 
              z = z.vertex[c(6, 14)], col = f.col, type = "l", add = T, 
              lwd = 1)
  rgl::plot3d(x = x.vertex[c(8, 14)], y = y.vertex[c(8, 14)], 
              z = z.vertex[c(8, 14)], col = f.col, type = "l", add = T, 
              lwd = 1)
  rgl::plot3d(x = x.vertex[c(8, 13)], y = y.vertex[c(8, 13)], 
              z = z.vertex[c(8, 13)], col = f.col, type = "l", add = T, 
              lwd = 1)
  rgl::plot3d(x = x.vertex[c(8, 11)], y = y.vertex[c(8, 11)], 
              z = z.vertex[c(8, 11)], col = f.col, type = "l", add = T, 
              lwd = 1)
  rgl::plot3d(x = x.vertex[c(5, 13)], y = y.vertex[c(5, 13)], 
              z = z.vertex[c(5, 13)], col = f.col, type = "l", add = T, 
              lwd = 1)
  rgl::plot3d(x = x.vertex[c(2, 9)], y = y.vertex[c(2, 9)], 
              z = z.vertex[c(2, 9)], col = f.col, type = "l", add = T, 
              lwd = 1)
  rgl::plot3d(x = x.vertex[c(2, 11)], y = y.vertex[c(2, 11)], 
              z = z.vertex[c(2, 11)], col = f.col, type = "l", add = T, 
              lwd = 1)
  rgl::plot3d(x = x.vertex[c(4, 11)], y = y.vertex[c(4, 11)], 
              z = z.vertex[c(4, 11)], col = f.col, type = "l", add = T, 
              lwd = 1)
  rgl::plot3d(x = x.vertex[c(4, 9)], y = y.vertex[c(4, 9)], 
              z = z.vertex[c(4, 9)], col = f.col, type = "l", add = T, 
              lwd = 1)
  rgl::plot3d(x = x.vertex[c(4, 10)], y = y.vertex[c(4, 10)], 
              z = z.vertex[c(4, 10)], col = f.col, type = "l", add = T, 
              lwd = 1)
  rgl::plot3d(x = x.vertex[c(2, 13)], y = y.vertex[c(2, 13)], 
              z = z.vertex[c(2, 13)], col = f.col, type = "l", add = T, 
              lwd = 1)
  rgl::plot3d(x = x.vertex[c(6, 13)], y = y.vertex[c(6, 13)], 
              z = z.vertex[c(6, 13)], col = f.col, type = "l", add = T, 
              lwd = 1)
  rgl::plot3d(x = x.vertex[c(3, 11)], y = y.vertex[c(3, 11)], 
              z = z.vertex[c(3, 11)], col = f.col, type = "l", add = T, 
              lwd = 1)
  rgl::text3d(x = E1[[1]], y = E1[[2]], z = E1[[3]], texts = vnames[[1]], 
                cex = 0.75, adj = c(1, 1))
  rgl::text3d(x = E2[[1]], y = E2[[2]], z = E2[[3]], texts = vnames[[2]], 
                cex = 0.75, adj = c(1, 1))
  rgl::text3d(x = E3[[1]], y = E3[[2]], z = E3[[3]], texts = vnames[[3]], 
                cex = 0.75, adj = c(1, 1))
  rgl::text3d(x = E4[[1]], y = E4[[2]], z = E4[[3]], texts = vnames[[4]], 
                cex = 0.75, adj = c(1, 1))
}



EMtriangle <- function (x, y, vnames=c("u","s","m"),
                        ylim=c(-0.9,0.9),
                        xlim=c(-0.9,0.9),
                        pch=16, bty="n",yaxt="n",xaxt="n",
                        col="black", asp=1, ann=FALSE, ...) {
  
  graphics::plot(x=x,y=y, pch=pch, bty=bty,yaxt=yaxt,xaxt=xaxt, col=col,
       ylim=ylim, xlim=xlim, asp=asp, ann=ann, ...)

  #vertices
  u<-c(-0.6495191, -0.3750000)
  s<-c(0.6495191, -0.3750000)
  m<-c(0.00, 0.75)

    graphics::polygon(x=c(u[[1]],s[[1]],m[[1]]),
                      y=c(u[[2]],s[[2]],m[[2]]))
    
    graphics::text(x=u[[1]],
                 y=u[[2]],
                 labels = vnames[[1]],
                 pos=1)
    graphics::text(x=s[[1]],
                 y=s[[2]],
                 labels = vnames[[2]],
                 pos=1)
    graphics::text(x=m[[1]],
                 y=m[[2]],
                 labels = vnames[[3]],
                 pos=3)
}


EMtetrahedron <- function (x, y, z, s.col = "red", f.col = "black", vnames = c("u","s","m","l"), 
                           type = "p", radius = 0.01, 
                           add = F, xlab = "", ylab = "", zlab = "",
                           box = F, axes = F, ylim = c(-0.75, 0.75), xlim = c(-0.75, 0.75),
                           zlim = c(-0.75, 0.75), aspect = T, mar = c(1, 1, 1, 1), ...) 
{
  rgl::rgl.viewpoint(zoom = 0.75)
  rgl::plot3d(x = x, y = y, z = z, col = s.col, type = type, 
              add = add, xlab = xlab, ylab = ylab, zlab = zlab, box = box, axes = axes, 
              radius = radius, ylim = ylim, xlim = xlim,
              zlim = zlim, aspect = aspect, mar = mar, ...)
  
  smu <- function(s, m, u, l) {
    x <- ((1 - 2 * s - m - u)/2) * sqrt(3/2)
    y <- (-1 + 3 * m + u)/(2 * sqrt(2))
    z <- u - 1/4
    r <- c(x, y, z)
    names(r) <- c("x", "y", "z")
    return(r)
  }
  rgl::plot3d(x = c(smu(1, 0, 0, 0)[["x"]], smu(0, 0, 0, 1)[["x"]]), 
              y = c(smu(1, 0, 0, 0)[["y"]], smu(0, 0, 0, 1)[["y"]]),
              z = c(smu(1, 0, 0, 0)[["z"]], smu(0, 0, 0, 1)[["z"]]), col = f.col, 
              type = "l", add = T, lwd = 1)
  rgl::plot3d(x = c(smu(0, 0, 0, 1)[["x"]], smu(0, 0, 1, 0)[["x"]]), 
              y = c(smu(0, 0, 0, 1)[["y"]], smu(0, 0, 1, 0)[["y"]]),
              z = c(smu(0, 0, 0, 1)[["z"]], smu(0, 0, 1, 0)[["z"]]), col = f.col, 
              type = "l", add = T, lwd = 1)
  rgl::plot3d(x = c(smu(1, 0, 0, 0)[["x"]], smu(0, 0, 1, 0)[["x"]]), 
              y = c(smu(1, 0, 0, 0)[["y"]], smu(0, 0, 1, 0)[["y"]]),
              z = c(smu(1, 0, 0, 0)[["z"]], smu(0, 0, 1, 0)[["z"]]), col = f.col, 
              type = "l", add = T, lwd = 1)
  rgl::plot3d(x = c(smu(0, 1, 0, 0)[["x"]], smu(0, 0, 1, 0)[["x"]]), 
              y = c(smu(0, 1, 0, 0)[["y"]], smu(0, 0, 1, 0)[["y"]]),
              z = c(smu(0, 1, 0, 0)[["z"]], smu(0, 0, 1, 0)[["z"]]), col = f.col, 
              type = "l", add = T, lwd = 1)
  rgl::plot3d(x = c(smu(0, 1, 0, 0)[["x"]], smu(0, 0, 0, 1)[["x"]]), 
              y = c(smu(0, 1, 0, 0)[["y"]], smu(0, 0, 0, 1)[["y"]]),
              z = c(smu(0, 1, 0, 0)[["z"]], smu(0, 0, 0, 1)[["z"]]), col = f.col, 
              type = "l", add = T, lwd = 1)
  rgl::plot3d(x = c(smu(0, 1, 0, 0)[["x"]], smu(1, 0, 0, 0)[["x"]]), 
              y = c(smu(0, 1, 0, 0)[["y"]], smu(1, 0, 0, 0)[["y"]]),
              z = c(smu(0, 1, 0, 0)[["z"]], smu(1, 0, 0, 0)[["z"]]), col = f.col, 
              type = "l", add = T, lwd = 1)

  rgl::text3d(x = smu(1, 0, 0, 0)[["x"]], y = smu(1, 0, 0, 0)[["y"]], 
                z = smu(1, 0, 0, 0)[["z"]], texts = vnames[[1]], cex = 1, adj = c(0, 
                                                                       0))
  rgl::text3d(x = smu(0, 1, 0, 0)[["x"]], y = smu(0, 1, 0, 0)[["y"]], 
                z = smu(0, 1, 0, 0)[["z"]], texts = vnames[[2]], cex = 1, adj = c(1, 
                                                                       1))
  rgl::text3d(x = smu(0, 0, 1, 0)[["x"]], y = smu(0, 0, 1, 0)[["y"]], 
                z = smu(0, 0, 1, 0)[["z"]], texts = vnames[[3]], cex = 1, adj = c(1, 
                                                                       1))
  rgl::text3d(x = smu(0, 0, 0, 1)[["x"]], y = smu(0, 0, 0, 1)[["y"]], 
                z = smu(0, 0, 0, 1)[["z"]], texts = vnames[[4]], cex = 1, adj = c(1, 
                                                                       1))
}


RNLthres <-function (photo=ncol(C)-1,
                     Rb, I, C, noise=TRUE, v=NA, n=NA, e=NA,
                     interpolate=TRUE, nm=seq(300,700,1)) {
  
  
  dependent <- FALSE
  nphoto = ncol(C) - 1
  photo1<-photo
  if (photo=="di") {photo1<-2}
  if (photo=="tri") {photo1<-3}
  if (photo=="tetra") {photo1<-4}
  
  #warnings
  ifelse(photo1 != nphoto, yes = warning("Argument 'C' has a number of sensitivity curves different than argument 'photo'.", 
                                         call. = FALSE), no = "")
  ifelse(photo1 != length(e) && noise == T, 
         yes = warning("Argument 'e' has a number of parameters different than 'photo'.", 
                       call. = FALSE), no = "")
  ifelse(photo1 != length(n) && noise == F, 
         yes = warning("Argument 'n' has a number of parameters different than 'photo'.", 
                       call. = FALSE), no = "")
  ifelse(any(ncol(I) > 2), yes = warning("'I' argument with more than two columns. Only the first two will be used.", 
                                         call. = FALSE), no = "")
  ifelse(any(ncol(Rb) > 2), yes = warning("'Rb' argument with more than two columns. Only the first two will be used.", 
                                          call. = FALSE), no = "")

  #Rb photon catch
  SRb<-vector(length=photo1)
  for (i in 1:photo1) {
    SRb[[i]] <- Q(I = I, R = Rb, C = C[, c(1, i+1)], interpolate = interpolate, nm = nm)
  }  
  #k parameter (Vorobyev and Osorio 1998 eq. 2)
  k<-1/SRb
  
  #noise
  noise_values<-vector(length=photo1)
  if (dependent == FALSE) {
    for (i in 1:photo1) {
      noise_values[[i]]<-noise_e(noise = noise, e = e[[i]], v = v, n = n[[i]])
    }
  }
  

  #photoreceptors
  Cnew<-C
  if (interpolate == TRUE) {
    Cnew<-data.frame(nm=nm)
    for (i in 1:photo1) {
      temp<-approx(x = C[, 1], y = C[,i+1], xout = nm, method = "linear")$y
      Cnew<-cbind(Cnew, temp)
    }
  }
  
  #kR Vorobyev and Osorio eq. 6
  kR<-matrix(ncol=photo1,nrow=nrow(Cnew))
  for (i in 1:length(k)) {
    kR[,i]<-k[[i]]*Cnew[,i+1]
  }

  #threshold
  if (photo1>2) {
    T1<-combn(x=noise_values, m=(photo1-1))
    for (i in 2:nrow(T1)) {
      T1[1,]<-T1[1,]*T1[i,]
    }
  T1<-sum(T1[1,]^2)
  }
  if(photo1==2) {T1<-sum(noise_values^2)}
  
  cols<-1:photo1
  cols.comb<-combn(x=cols, m=2)
  T2.kR<-matrix(nrow=nrow(kR), ncol=ncol(cols.comb))
  for (i in 1:ncol(cols.comb)) {
    T2.kR[,i]<-kR[,cols.comb[1,i]]-kR[,cols.comb[2,i]]
  }
  T2.kR<-T2.kR^2
  
  if(photo1==2) {T2.e<-1}
  if(photo1==3) {T2.e<-noise_values[photo1:1]^2}
  if(photo1>=4) {
    T2.e<-combn(x=noise_values, m=(photo1-2))
    for (i in 2:nrow(T2.e)) {
      T2.e[1,]<-T2.e[1,]*T2.e[i,]
    }
    T2.e<-T2.e[1,]^2
    T2.e<-T2.e[length(T2.e):1]
  }

  T2<-T2.kR
  for (i in 1:ncol(T2)) {
    T2[,i]<-T2.kR[,i]*T2.e[[i]]
  }
  if(ncol(T2)>1) {
    T2<-rowSums(T2)
  }
  
  Thres <- sqrt(T1/T2)

  S<-log10(1/Thres)
  
  #results
  r<-data.frame(nm,Thres,S)
  names(r)<-c("nm","T","S")
  
  class(r)<-c("colourvision", "data.frame")
  attr(r, "model_name") <- "RNL Threshold"
  attr(r, "n_photor_types") <- photo1
  attr(r, "Rb") <- Rb
  attr(r, "I") <- I
  attr(r, "C") <- C
  attr(r, "noise calculated") <- noise
  attr(r, "v") <- v
  attr(r, "n") <- n
  attr(r, "Interpolate") <- interpolate
  attr(r, "nm") <- nm
  return (r)
  
}

plot.colourvision <- function (x, ...) {
  photo1<-attributes(x)$n_photor_types
  if(photo1==4) {stop("For a 3D plot use 'plot3d'.")}
  if(photo1>4) {stop("Plotting is not available for > 3-dimentions.")}
  model<-attributes(x)$model_name
  if (model=="Colour hexagon model") {
    if (photo1==3) {
      CTTKhexagon(x=x[,"X1"],y=x[,"X2"], ...)
    }
    if (photo1==2) {
      plot(x=x[,"X1"],y=rep(0, length(x[,"X1"])), ylim=c(-1,1),xlim=c(-1,1),
           ann=FALSE, axes = FALSE,
           panel.first=c(
             segments(x0=0,x1=1,y0=0,y1=0),
             segments(x0=0,x1=-1,y0=0,y1=0),
             segments(x0=1,x1=1,y0=-0.03,y1=0.03),
             segments(x0=-1,x1=-1,y0=-0.03,y1=0.03),
             points(x=0,y=0,pch=4), ...))
    }
  }
  if (model=="Endler and Mielke model") {
    if (photo1==3) {
      EMtriangle(x=x[,"X1"],y=x[,"X2"], ...)
    }
    if (photo1==2) {
      plot(x=x[,"X1"],y=rep(0, length(x[,"X1"])), ylim=c(-0.75,0.75),xlim=c(-0.75,0.75),
           ann=FALSE, axes = FALSE, xlab="x",
           panel.first=c(
             segments(x0=0,x1=0.75,y0=0,y1=0),
             segments(x0=0,x1=-0.75,y0=0,y1=0),
             segments(x0=0.75,x1=0.75,y0=-0.03,y1=0.03),
             segments(x0=-0.75,x1=-0.75,y0=-0.03,y1=0.03),
             points(x=0,y=0,pch=4)), ...)
    }     
  }
  
  if (model=="Receptor noise limited model") {
    if (photo1==3) {
      plot(x=x[,"X1_R1"],y=x[,"X2_R1"],
           ylab="y", xlab="x", ...)
    }
    if (photo1==2) {
      lim<-max(abs(c(max(x[,"X1_R1"]), min(x[,"X1_R1"]))))
      plot(x=x[,"X1_R1"],y=rep(0, length(x[,"X1_R1"])), ylim=c(-lim,lim),xlim=c(-lim,lim),
           ann=FALSE, axes = FALSE, xlab="x",
           panel.first=c(
             segments(x0=0,x1=lim,y0=0,y1=0),
             segments(x0=0,x1=-lim,y0=0,y1=0),
             segments(x0=lim,x1=lim,y0=-lim*0.03,y1=lim*0.03),
             segments(x0=-lim,x1=-lim,y0=-lim*0.03,y1=lim*0.03),
             points(x=0,y=0,pch=4)), ...)
    }   
  }
  
  if (model=="RNL Threshold") {
    plot(x$S~x$nm, xlab="Wavelength(nm)", ylab="Log Sensitivity", type="l")
  }
}

plot3d.colourvision <- function (x, ...) {
  photo1<-attributes(x)$n_photor_types
  if(photo1<=3) {stop("For a 2D plot use 'plot'.")}
  if(photo1>4) {stop("Plotting is not available for > 3-dimentions.")}
  model<-attributes(x)$model_name
  if(model=="RNL Threshold") {stop("For a colour threshold use 'plot'.")}
  if (model=="Colour hexagon model") {
    if (photo1==4) {
      CTTKhexagon3D(x=x[,"X1"],y=x[,"X2"], z=x[,"X3"], ...)
    }
  }
  if (model=="Endler and Mielke model") {
    if (photo1==4) {
      EMtetrahedron(x=x[,"X1"],y=x[,"X2"], z=x[,"X3"], ...)
    }
  }
  
  if (model=="Receptor noise limited model") {
    if (photo1==4) {
      plot3d(x=x[,"X1_R1"],y=x[,"X2_R1"], z=x[,"X3_R1"],
           ylab="y", xlab="x", zlab="z", ...)
    }
  }
}