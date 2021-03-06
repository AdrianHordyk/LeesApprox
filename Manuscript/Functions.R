
# Make data.frame of Operating Model parameters
MakeOM_DF <- function(OMlist, LinfCV=0.1, maxsd=2) {
  lout <- list()
  for (i in 1:length(OMlist)) {
    Name <- names(OMlist)[i]
    Stock <- OMlist[[i]]
    Linf <- mean(Stock@Linf)
    K <- mean(Stock@K)
    t0 <- mean(Stock@t0)
    M <- mean(Stock@M)
    maxage <- ceiling(-log(0.001)/M)
    
    L50 <- mean(Stock@L50)
    L95 <- L50 +  mean(Stock@L50_95)
    
    L5 <- mean(Stock@L5) * L50
    LFS <- mean(Stock@LFS) * L50
    Vmaxlen <- mean(Stock@Vmaxlen)
    sigmaR <- mean(Stock@Perr)
    steepness <- mean(Stock@h)
    
    alpha <- Stock@a
    beta <- Stock@b
    
    lout[[i]] <- data.frame(Stock=Name, Linf=Linf, K=K,
               t0=t0, M=M, maxage=maxage, L50=L50, L95=L95,
               L5=L5, LFS=LFS, Vmaxlen=Vmaxlen, sigmaR=sigmaR,
               steepness=steepness, alpha=alpha, beta=beta, LinfCV=LinfCV,
               maxsd=maxsd, maxL=Linf + LinfCV*Linf*maxsd)
  }
  do.call('rbind', lout)
}


EqGTG <- function(x, DF, ngtg=5, Fmulti=1, binwidth=1, R0=1E4) { # equilibrium GTG model
  
  for (i in 1:ncol(DF)) {
    assign(names(DF)[i], DF[x,i])
  }

  distGTG <- seq(from=-maxsd, to=maxsd, length.out = ngtg)
  rdist <- dnorm(distGTG, 0, 1)/sum(dnorm(distGTG, 0, 1))
  Linfgtg <- Linf + LinfCV*Linf*distGTG
  
  LenBins <- seq(0, to=max(Linfgtg), by=binwidth)
  Nbins <- length(LenBins) -1
  By <- LenBins[2] - LenBins[1]
  LenMids <- seq(from=By*0.5, by=By, length.out = Nbins)
  
  ages <- 0:maxage # in years

  LAA <- array(NA, dim=c(maxage+1, ngtg)) # Length-at-age by GTG 
  ind <- as.matrix(expand.grid(ages+1,1:ngtg))
  LAA[ind] <- (Linfgtg[ind[,2]] * (1-exp(-K *(ages[ind[,1]]-t0))))
  WAA <- alpha * LAA^beta  # Weight-at-age by GTG
  
  # Maturity
  
  L50gtg <- matrix(L50/Linf * Linfgtg, ncol=ngtg, nrow=maxage+1,byrow=TRUE)  # assume constant L50/Linf across GTGs
  L95gtg <- array(NA, dim=dim(L50gtg))
  L95gtg[ind] <- L50gtg[ind] + (L95-L50)
  MAA <- 1/(1 + exp(-log(19) * ((LAA - L50gtg)/(L95gtg-L50gtg)))) # Maturity-at-age by GTG
  
  # selectivity-at-length - fishery
  sl <- (LFS - L5) /((-log(0.05,2))^0.5)
  sr <- (Linf - LFS) / ((-log(Vmaxlen,2))^0.5) # selectivity parameters are constant for all years
  SAA <- matrix(dnormal(LAA, LFS, sl, sr), nrow=maxage+1, ncol=ngtg) # selectivty-at-age by GTG 
  SL <- dnormal(LenMids, LFS, sl, sr)
  
  M_array <- matrix(M, nrow=maxage+1, ncol=ngtg) 
  FAA <- array(NA, dim=dim(SAA))
  Fmort <- Fmulti * M 
  FAA <- SAA * Fmort # fishing mortality at age by GTG 
  ZAA <- FAA + M_array # Z-at-age by GTG 
  
  NAA <- CAA <- matrix(NA, nrow=maxage+1, ncol=ngtg)
  NAA[1,] <- R0 * rdist
  for (a in 1:maxage) {
    ageind <- a + 1 
    NAA[ageind,] <- NAA[ageind-1,] * exp(-ZAA[ageind-1,])
  }
  CAA <- FAA/ZAA * (1-exp(-ZAA*NAA)) # catch in numbers
  
  out <- list()
  df <- data.frame(Age=ages, 
                   N=as.vector(NAA),
                   Select=as.vector(SAA), 
                   Length=as.vector(LAA),
                   Weight=as.vector(WAA),
                   GTG=rep(1:ngtg, each=(maxage+1)),
                   Linf=rep(Linfgtg, each=(maxage+1)),
                   CAA=as.vector(CAA))
  out$df <- df
  out$LenBins <- LenBins
  out$LenMids <- LenMids
  out
}















addLines <- function(LenBins) {
  for (x in 1:length(LenBins)) {
    abline(v=LenBins[x], col="lightgray", lty=3)
  }
}

addBars <- function(probdf, col2) {

  temp <- probdf 
  for (x in 1:nrow(probdf)) {
    polygon(x=c(temp$Bin1[x], temp$Bin2[x], temp$Bin2[x], temp$Bin1[x]),
            y=c(0, 0, temp$Prob[x], temp$Prob[x]), col=col2, xpd=NA)
  }
}


Ftrend <- function(yr.st, yr.end, curF,
                   F.pat= c('stable', 'inc', 'dec'),
                   Fcv=0.3, plot=TRUE) {
  F.pat <- match.arg(F.pat)
  yrs <- yr.st:yr.end
  yr.mid <- ceiling(mean(yrs))
  yr.ind <- which(yrs == yr.mid)
  nyrs <- length(yrs)
  Ferr <- exp(rnorm(nyrs, 0, Fcv))
  if (F.pat == "inc") {
    Ftrend <- seq(from=0, to=curF, length.out=nyrs) * Ferr
    Ftrend <- Ftrend/(Ftrend[length(Ftrend)]/curF)

  } else if (F.pat == "dec") {
    Ftrend <- c(seq(from=0, to=2*curF, length.out=yr.ind),
                seq(from=2*curF, to=curF, length.out=nyrs-yr.ind)) * Ferr
    Ftrend <- Ftrend/(Ftrend[length(Ftrend)]/curF)
  } else if (F.pat == "stable") {
    Ftrend <- c(seq(from=0, to=curF, length.out=yr.ind),
                seq(from=curF, to=curF, length.out=nyrs-yr.ind)) * Ferr
  }

  if (plot) {
    plot(yrs, Ftrend, type="l", bty="l", las=1, xlab="Years", ylab='apical Fishing mortality')
  }

  Ftrend
}

BHSRR <- function(SBcurr, SB0, R0, steepness) {
  (4 * R0 * steepness * SBcurr)/(SB0/R0 * R0 * (1-steepness) + (5*steepness-1)*SBcurr)
}


GTGpopsim <- function(Linf, K, t0, M, L50, L95, LFS, L5, Vmaxlen, sigmaR, steepness,
                      annualF,alpha=1E-5, beta=3, LinfCV=0.1, ngtg=101, maxsd=2, binwidth=1,
                      R0=1000) {

  # Growth-type-groups
  distGTG <- seq(from=-maxsd, to=maxsd, length.out = ngtg)
  rdist <- dnorm(distGTG, 0, 1)/sum(dnorm(distGTG, 0, 1))
  Linfgtg <- Linf + LinfCV*Linf*distGTG

  LenBins <- seq(0, to=max(Linfgtg), by=binwidth)
  Nbins <- length(LenBins) -1
  By <- LenBins[2] - LenBins[1]
  LenMids <- seq(from=By*0.5, by=By, length.out = Nbins)

  nyrs <- length(annualF)

  maxage <- ceiling(-log(0.01)/M) # maxium age
  ages <- 1:maxage # in years

  maxAgeind <- maxage # maxage + 1
  ageVec <- 1:maxAgeind
  ind <- as.matrix(expand.grid(1:nyrs, ageVec,1:ngtg))
  LAA <- array(NA, dim=c(nyrs, maxAgeind, ngtg)) # Length-at-age by GTG and year
  LAA[ind] <- (Linfgtg[ind[,3]] * (1-exp(-K *(ages[ind[,2]]-t0))))
  WAA <- alpha * LAA^beta  # Weight-at-age by GTG and year

  # Maturity
  L50gtg <- L95gtg <- array(NA, dim=c(nyrs, maxAgeind, ngtg))
  L50gtg[ind] <- L50/Linf * Linfgtg  # assume constant L50/Linf but allow L50 to vary by year
  L95gtg[ind] <- L50gtg[ind] + (L95-L50)
  MAA <- 1/(1 + exp(-log(19) * ((LAA - L50gtg)/(L95gtg-L50gtg)))) # Maturity-at-age by GTG and year

  # selectivity-at-length - fishery
  sl <- (LFS - L5) /((-log(0.05,2))^0.5)
  sr <- (Linf - LFS) / ((-log(Vmaxlen,2))^0.5) # selectivity parameters are constant for all years
  SAA <- array(dnormal(LAA, LFS, sl, sr), dim=c(nyrs, maxAgeind, ngtg)) # selectivty-at-age by GTG and year
  SL <- dnormal(LenMids, LFS, sl, sr)

  M_array <- array(M, dim=c(nyrs, maxAgeind, ngtg))
  FAA <- array(NA, dim=c(nyrs, maxAgeind, ngtg))
  FAA[ind] <- SAA * annualF[ind[,1]] # fishing mortality at age by GTG and year
  ZAA <- FAA + M_array # Z-at-age by GTG and year


  # Unfished Year One
  Nunfished <- CAA <- array(NA, dim=c(nyrs, maxAgeind, ngtg))
  SB <- array(NA, dim=c(nyrs, maxage, ngtg))
  Nunfished[1,1,] <- rdist * R0 # distribute virgin recruitment
  Nunfished[1,2:maxAgeind,] <- matrix(Nunfished[1,1,], nrow=maxage-1, ncol=ngtg, byrow=TRUE) *
    exp(-apply(M_array[1,ageVec-1,], 2, cumsum))

  SB[1,,] <- Nunfished[1,,] * WAA[1,,] * MAA[1,,]
  SB0 <- sum(SB[1,,])
  SBcurr <- Rec <- rep(NA, nyrs+1)
  SBcurr[1] <- SB0
  Rec[1] <- R0
  
  Nfished <- Nunfished
  
  CAA[1,,] <- FAA[1,,]/ZAA[1,,] * (1-exp(-ZAA[1,,] * Nfished[1,,]))
  
  recmu <- -0.5 * (sigmaR)^2
  recdevs <- exp(rnorm(nyrs, recmu, sigmaR))
  for (yr in 2:nyrs) {
    # message("Year ", yr, " of ", nyrs)
    Rec[yr] <- R0 # BHSRR(SBcurr[yr-1], SB0, R0, steepness) # recruitment
    Nfished[yr,1,] <- Rec[yr] * recdevs[yr] * rdist
    Nfished[yr,2:maxAgeind,] <- Nfished[yr-1,1:(maxage-1),] * exp(-ZAA[yr-1,1:(maxage-1),])
    SB[yr,,] <- Nfished[yr,,] * WAA[yr,,] * MAA[yr,,]
    SBcurr[yr] <- sum(SB[yr,,])
    CAA[yr,,] <- FAA[yr,,]/ZAA[yr,,] * (1-exp(-ZAA[yr,,] * Nfished[yr,,]))
    
  }

  df <- data.frame(Yr=1:nyrs, Age=rep(ages, each=nyrs), N=as.vector(Nfished),
                   Select=as.vector(SAA), Length=as.vector(LAA),
                   Weight=as.vector(WAA),
                   GTG=rep(1:ngtg, each=(nyrs*(maxage))),
                   Linf=rep(Linfgtg, each=(nyrs*(maxage))),
                   CAA=as.vector(CAA))
  out <- list(df=df, LenBins=LenBins, LenMids=LenMids, N=Nfished)

  out
}


