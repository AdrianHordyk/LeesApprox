
# --- Code for Lee's Approximation Manuscript ----
# A. Hordyk
# December 2019

library(LeesApproxTMB)
library(dplyr)
library(Rcpp)
library(ggplot2)
library(cowplot)

DLMextra()
library(DLMextra)

setwd('Manuscript') 
source('Functions.r')
Rcpp::sourceCpp('src/LeesApprox.cpp') # load CPP version of Lee's approx


# ---- Operating Models ---- 
OMlist <- list('Queen triggerfish'=DLMextra::Queen_Triggerfish_STT_NOAA,
               'Stoplight parrotfish'=DLMextra::Stoplight_Parrotfish_STX_NOAA,
               'Yellowtail snapper'=DLMextra::Yellowtail_Snapper_PR_NOAA)

OM_DF <- MakeOM_DF(OMlist)

ngtg <- 7
binwidth <- 5

OM_DF <- data.frame(Stock="example",
                    Linf=100, K=0.13333, t0=0, M=0.2, maxage=ceiling(-log(0.01)/0.2),
                    L50=66, L95=70, L5=40, LFS=50, Vmaxlen=1, sigmaR=0, steepness=0.99,
                    alpha=1E-5, beta=3, LinfCV=0.1, maxsd=3, maxL=100+3*0.1*100)

run1 <- EqGTG(1, OM_DF, ngtg=ngtg, Fmulti=1, binwidth = binwidth)
DF <- run1[[1]]  
DF <- DF %>% group_by(Age) %>% mutate(N2=N/sum(N), N3=0) 



library(scatterplot3d)
x <- DF$Age
y <- DF$Length
z <- DF$N2

s <- scatterplot3d(x=range(DF$Age),y=c(0, max(DF$Length)),z=c(0,1), 
                   angle=45, grid=TRUE, box=FALSE,type="n",
                   xlab="Age", ylab="Length", zlab="Relative Frequency")

len <- DF$Length
sel <- DF$Select
len2 <- len[order(len)]
sel2 <- sel[order(len)]

s$points3d(rep(6, length(sel)), len2, sel2, type="l", col="slategray", lwd=2, lty=2)

library(colorspace)

colors <-  qualitative_hcl(ngtg, palette = "Dark 3")

for (x in 1:ngtg) {
  temp <- DF %>% filter(GTG==x)
  s$points3d(temp$Age, temp$Length, temp$N3, col=colors[x], type="l")
}

scaler <- seq(1, to=0.2, length.out=16) * 3

# age <- 15
# temp <- DF %>% filter(Age==age)


for (age in c(2,  6, 10, 15)) {
  temp <- DF %>% filter(Age==age)
  s$points3d(temp$Age, temp$Length, temp$N2*scaler[age], col='black', type="l", lwd=2, lty=2)
  # s$points3d(temp$Age, temp$Length, temp$N2/max(temp$N2), col=colors, type="h", pch=16)
  s$points3d(temp$Age, temp$Length, temp$N2*scaler[age], col=colors, type="h", pch=16,
             cex=1.2)
}

Figure1(ngtg, binwidth)


# 
age <- 15
temp <- DF %>% filter(Age==age)

b1 <- max(which(run1$LenBins < min(temp$Length))) + 1
b2 <- min(which(run1$LenBins > max(temp$Length)))

plot(c(run1$LenBins[b1], run1$LenBins[b2]), c(0, max(temp$N2)), type="n", pch=16,
     yaxs="i", xaxs="i", bty="l")
addLines(run1$LenBins)
lines(temp$Length, temp$N2, col='black', lwd=2, lty=2, xpd=NA)
points(temp$Length, temp$N2, col=colors, pch=16, cex=1.2, xpd=NA)





run2 <- EqGTG(1, OM_DF, ngtg=ngtg, Fmulti=1, binwidth = binwidth)
DF2 <- run2[[1]]  
DF2 <- DF2 %>% group_by(Age) %>% mutate(N2=N/sum(N), N3=0) 


# add fished lines
temp <- DF2 %>% filter(Age==age)
lines(temp2$Length, temp2$N2, col='darkgray', lwd=2, lty=2, xpd=NA)
points(temp2$Length, temp2$N2, col=colors, pch=16, cex=1.2, xpd=NA)





addBars <- function(probdf, col2) {
  yrs <- unique(probdf$Yr)
  temp <- probdf %>% dplyr::filter(Yr==yrs[2])
  for (x in 1:nrow(probdf)) {
    polygon(x=c(temp$Bin1[x], temp$Bin2[x], temp$Bin2[x], temp$Bin1[x]),
            y=c(0, 0, temp$Prob[x], temp$Prob[x]), col=col2)
  }
}




temp2 <- DF2 %>% filter(Age==age)
# s$points3d(temp$Age, temp$Length, temp$N2/max(temp$N2), col=colors, type="h", pch=16)
s$points3d(temp2$Age, temp2$Length, temp2$N2*scaler[age], col='darkgray', type="l",
           lty=3, lwd=2)













DF2 <- run1[[1]]
tt <- DF2 %>% filter(Yr==Age)

len <- unique(DF2$Length)
sel <- unique(DF2$Select)
len2 <- len[order(len)]
sel2 <- sel[order(len)]



t2 <- tt %>% group_by(Age) %>% mutate(N2=N/sum(N)) 
t2$N3 <- 0

t3 <- t2 %>% filter(Age %in% c(1,4, 6, 12))

library(scatterplot3d)
x <- c(0,t3$Age)
y <- c(0,t3$Length)
z <- c(0,t3$N2)

s <- scatterplot3d(x,y,z, angle=45, grid=TRUE, box=FALSE,type="h", pch=16,
                   xlab="Age", ylab="Length", zlab="Relative Frequency", color=c(0,rep(1:5, each=4)))

s$points3d(0:12, c(0,t2$Length[1:12]), rep(0, 13), type="l")
s$points3d(0:12, c(0,t2$Length[49:60]), rep(0, 13), type="l")

tt <- filter(t2, Age==12)
s$points3d(tt$Age, tt$Length, tt$N2, type="l")




x <- c(0,t3$Age)
y <- c(0,t3$Length)
z <- c(0,t3$N2)

s <- scatterplot3d(x,y,z, angle=45, grid=TRUE, box=FALSE,type="h", pch=16,
                   xlab="Age", ylab="Length", zlab="Relative Frequency", color=c(0,rep(1:5, each=4)))


s$points3d(0:12, c(0,t2$Length[1:12]), rep(0, 13), type="l")
s$points3d(0:12, c(0,t2$Length[49:60]), rep(0, 13), type="l")

tt <- filter(t2, Age==12)
s$points3d(tt$Age, tt$Length, tt$N2, type="l")




plot3d(x,y,z)

x <- t2$Age
y <- t2$Length
z <- t2$N2

plot3d(x,y,z, add=TRUE)


matplot(tt$Age, tt$Length, type="b")


devtools::install_github("AckerDWM/gg3D")


# standardize to sum to one in each age-cohort





segments3d(x[2:3],y[2:3],z[2:3],col=2,lwd=2)


theta=0; phi=20
ggplot(t2, aes(x=Age, y=Length, z=n3, color=as.factor(GTG))) +
  axes_3D(theta=theta, phi=phi) +
  stat_3D(theta=theta, phi=phi, geom="path") +
  axes_3D(theta=theta, phi=phi) +
  stat_3D(theta=theta, phi=phi) +
  axis_labs_3D(theta=theta, phi=phi, size=3, 
               hjust=c(1,1,1.2,1.2,1.2,1.2), 
               vjust=c(-.5,-.5,-.2,-.2,1.2,1.2)) +
  labs_3D(theta=theta, phi=phi, 
          hjust=c(1,0,0), vjust=c(1.5,1,-.2),
          labs=c("Age", "Length", "Z")) +
  theme_void()
 

t3 <- t2 %>% filter(Age==12)
plot(t3$Length, t3$N2, type="b")


lattice::wireframe(N ~ Age * Length, data=t2)




x <- 1:5/10
y <- 1:5
z <- x %o% y
z <- z + .2*z*runif(25) - .1*z

library(rgl)
persp3d(x, y, z, col="skyblue")
t2 <- t2 %>% arrange(Age, GTG, Length)



persp3d(t2$Age, t2$Length, t2$N2, col="skyblue")









# ---- Figure 1 ----

source('Figure_1.r')

Figure1()

# ---- Figure 2 ----
pars <- list()
pars$ngtg <- 7
pars$LinfCV <- 0.1
pars$Linf <- 100
pars$maxsd <- 2
pars$binwidth <- 5

pars$K <- 0.25
pars$t0 <- 0
pars$M <- 0.2
pars$maxage <- ceiling(-log(0.01)/pars$M)
pars$Ages <- 1:pars$maxage

pars$L5 <- 30
pars$LFS <- 60
pars$Vmax <- 0.6

pars$yr.st <- 1950
pars$yr.end <- 2017
pars$FVec <- c(0, 0.2, 0.6)


source('Figure_2.r')

Figure2(pars)

# ---- Figure 3 ----
source('Figure_3.r')

Figure3(pars)


# ---- Sensitivity Tests (Figure 4) ----
source('Figure_4.r')
Fig4DF <- Figure4()




# ---- Assessment Model (Figure 5) ----

AgeSampSize <- 250
LengthSampSize <- 250

Cobs <- rlnorm(length(annualF), 0, 0.1)
Iobs <- rlnorm(length(annualF), 0, 0.1)


SimPop <- GTGpopsim(Linf, K, t0, M, L50, L95, LFS, L5, Vmaxlen, sigmaR,
                    steepness, annualF,alpha, beta, LinfCV, ngtg=1001, maxsd,
                    binwidth, R0=1E5)

DF <- SimPop[[1]]
DF <- DF %>% group_by(Yr) %>% mutate(Catch=sum(CAA*Weight),
                                     TotalB=sum(N*Weight))

TSData <- DF %>% select(Yr, Catch, Index=TotalB) %>% distinct()
TSData$Index <- TSData$Index/mean(TSData$Index)

# Catch-at-age
DF$vulnN <- DF$N * DF$Select
CAA_DF <- DF %>% group_by(Yr, Age) %>% summarize(CAA=sum(CAA))
VAge_Comp <- CAA_DF %>% tidyr::pivot_wider(names_from='Yr',
                                           values_from="CAA")

AgeSamps <- sapply(2:length(annualF), function(i) 
  rmultinom(n=1, size=AgeSampSize, VAge_Comp[[i+1]]))

AgeSamps <- cbind(rep(0, max(DF$Age)), AgeSamps)
CAA <- t(AgeSamps)

# Catch-at-length
CAL <- matrix(0, nrow=length(annualF), ncol=length(SimPop$LenMids))
for (yr in 1:length(annualF)) {
  df <- DF %>% filter(Yr==yr)
  lenP <- rep(0, length(SimPop$LenMids))
  for (l in seq_along(SimPop$LenMids)) {
    ind <- df$Length >= SimPop$LenBins[l] & df$Length < SimPop$LenBins[l+1]
    lenP[l] <- sum(df$vulnN[ind])
  }
  lenP <- lenP/sum(lenP)
  if (!all(lenP == 0)) {
    CAL[yr,] <- t(rmultinom(n=1, size=LengthSampSize, prob=lenP))
  } 
}

# TO DO - Add Obs Error 
Data <- new("Data")
Data@CAL_bins <- SimPop$LenBins
Data@CAL_mids <- SimPop$LenMids
Data@CAL <- array(CAL, dim=c(1, length(annualF), length(SimPop$LenMids)))
Data@CAA <- array(CAA, dim=c(1, length(annualF), max(DF$Age)))
Data@Year <- unique(DF$Yr)
Data@Cat <- matrix(TSData$Catch * Cobs, nrow=1)
Data@Ind <- matrix(TSData$Index * Iobs, nrow=1)

Data@MaxAge <- max(DF$Age)
Data@Mort <- M
Data@vbt0 <- t0
Data@vbK <- K
Data@vbLinf <- Linf
Data@L50 <- L50 
Data@L95 <- L95 
Data@wla <- alpha
Data@CV_vbLinf <- LinfCV
Data@wlb <- beta
Data@steep <- steepness
Data@sigmaR <- sigmaR


CAA_multiplier <- 50
CAL_multiplier <- 0

# Fit Assessment Models 
ngtg_assess <- 3
Mod1 <- SCA_GTG(Data = Data, ngtg=ngtg_assess, 
                CAA_multiplier=CAA_multiplier,
                CAL_multiplier= CAL_multiplier)
Mod2 <- SCA_GTG(Data = Data, ngtg=ngtg_assess, use_LeesEffect = FALSE,
                CAA_multiplier=CAA_multiplier,
                CAL_multiplier= CAL_multiplier)
Mod3 <- SCA(Data=Data,
            CAA_multiplier=CAA_multiplier,
            CAL_multiplier= CAL_multiplier)

plot(annualF, type="l", ylim=c(0, max(annualF*1.5)))
lines(Mod1@FMort, col='blue')
lines(Mod2@FMort, col="green")
lines(Mod3@FMort, col="red")



plot(TSData$Index/max(TSData$Index), type="l", ylim=c(0, 1.5))
lines(Mod1@B_B0, col='blue')
lines(Mod2@B_B0, col="green")
lines(Mod3@B_B0, col="red")


MSEtool::compare_models(Mod1, Mod2, label = c("Lee's Effect", "No Effect"))


source("https://raw.githubusercontent.com/quang-huynh/LeesApprox/master/SCA_GTG_markdown.R")
plot(Mod1)

Mod1@info$data$CAA_n
Obs_C_at_age <- Mod1@Obs_C_at_age
C_at_age <- Mod1@C_at_age
info <- Mod1@info
ind_valid <- rowSums(Obs_C_at_age, na.rm = TRUE) > 0
plot_composition(info$Year[ind_valid], Obs_C_at_age[ind_valid, ], C_at_age[ind_valid, ], plot_type = "annual", ages = NULL, N = info$data$CAA_n[ind_valid])


Mod1@FMort 
Mod2@FMort


# TODO - add some obs error to catch and index 


system.time(
  GTG_3 <- SCA_GTG(Data = SimulatedData, truncate_CAL = FALSE, ngtg=11)
)

# Turn off Lee's Effect. Runtime of 1.1 seconds
system.time(
  GTG_3_noLee <- SCA_GTG(Data = SimulatedData, use_LeesEffect = FALSE, ngtg=11)
)

MSEtool::compare_models(GTG_3, GTG_3_noLee, label = c("Lee's Effect", "No Effect"))
GTG_3@FMort
GTG_3_noLee@FMort

