
Figure1 <- function(save=TRUE) {
  if (save) png('Figures/Figure1.png', units='mm',
                height=100, width=170, res=300)
  ngtg <- 11
  binwidth <- 5
  colors <-  qualitative_hcl(ngtg, palette = "Dark 3")
  
  OM_DF <- data.frame(Stock="example",
                      Linf=100, K=0.13333, t0=0, M=0.2, maxage=ceiling(-log(0.01)/0.2),
                      L50=66, L95=70, L5=45, LFS=50, Vmaxlen=1, sigmaR=0, steepness=0.99,
                      alpha=1E-5, beta=3, LinfCV=0.1, maxsd=3, maxL=100+3*0.1*100)
  
  run1 <- EqGTG(1, OM_DF, ngtg=ngtg, Fmulti=1, binwidth = binwidth)
  run2 <- EqGTG(1, OM_DF, ngtg=ngtg, Fmulti=0, binwidth = binwidth) # unfished
  DF <- run1[[1]]  
  DF <- DF %>% group_by(Age) %>% mutate(N2=N/sum(N), N3=0) 
  
  layout.matrix <- matrix(c(1,1,2,3), nrow=2, ncol=2)
  layout(mat = layout.matrix,
         heights = c(1, 1), # Heights of the two rows
         widths = c(1,0.5)) # Widths of the two columns
  
  # layout.show(3)
  par(oma=c(0,0,0,0))
  maxY <- ceiling(max(DF$Length)/5) * 5
  
  s <- scatterplot3d(x=range(DF$Age),y=c(0, maxY),z=c(0,1), 
                     y.margin.add=0.5,
                     angle=45, grid=TRUE, box=FALSE,type="n",
                     xlab="Age", ylab="", zlab="Relative Frequency",
                     mar= c(3,3,1,1), xpd=NA)
  text(x = 8, y = 1.5, 'Length', srt = 45)
  text(0.25, 5.5, 'a)', xpd=NA)
  
  Lens <- seq(0, 140, by=1)
  sl <- (OM_DF$LFS - OM_DF$L5) /((-log(0.05,2))^0.5)
  sr <- (OM_DF$Linf - OM_DF$LFS) / ((-log(OM_DF$Vmaxlen,2))^0.5) 
  sel <- dnormal(Lens, OM_DF$LFS, sl, sr) 
  s$points3d(rep(5, length(sel)), Lens, sel, type="l", col="slategray", lwd=2, lty=1)
  # for (l in seq_along(Lens)) {
  #   x <- c(5,5)
  #   y <- c(l,l)
  #   z <- sel[which(Lens==l)]
  #   z <- c(0, z)
  #   s$points3d(x,y,z, type="l", col="slategray", lwd=1, lty=2)
  # }
  
  for (x in 1:ngtg) {
    temp <- DF %>% filter(GTG==x)
    s$points3d(temp$Age, temp$Length, temp$N3, col=colors[x], type="l")
  }
  
  scaler <- seq(1, to=0.1, length.out=max(DF$Age)) * 3.5
  
  for (age in c(1,  5, 10, 15, 20,24)) {
    temp <- DF %>% filter(Age==age)
    s$points3d(temp$Age, temp$Length, temp$N2*scaler[age], col='black', type="l", lwd=2, lty=2)
    # s$points3d(temp$Age, temp$Length, temp$N2/max(temp$N2), col=colors, type="h", pch=16)
    s$points3d(temp$Age, temp$Length, temp$N2*scaler[age], col=colors, type="h", pch=16,
               cex=1)
  }
  
  # Age 15 plot 
  DF2 <- DF %>% filter(Age==15)
  DF3 <- run2$df %>% filter(Age==15)
  DF3 <- DF3 %>% mutate(N2=N/sum(N))
  
  b1 <- run1$LenBins[max(which(run1$LenBins < min(DF2$Length)))]
  b2 <- run1$LenBins[min(which(run1$LenBins > max(DF2$Length)))]
  yline <- 2.8
  ymax <- max(DF3$N2)*1.05
  par(mar=c(2,2,1,0))
  plot(c(b1, b2), c(0, ymax), type="n", yaxs="i", xlab="", ylab="",
       bty="n", las=2, axes=FALSE)
  axis(side=1, labels=FALSE)
  axis(side=2, las=2)
  mtext(side=2, "Relative Frequency", line=yline)
  text(b1, ymax*1.05, 'b)', xpd=NA)
  
  addLines(run1$LenBins)
  points(DF3$Length, DF3$N2, pch=16, lty=2, lwd=1, type='b', col="darkgray",xpd=NA)
  
  lines(DF2$Length, DF2$N2, pch=16, lty=2, lwd=1)
  points(DF2$Length, DF2$N2, pch=16, cex=1, col=colors, xpd=NA)
  
  # Histogram compared to high res model 
  run1a <- EqGTG(1, OM_DF, ngtg=1001, Fmulti=1, binwidth = binwidth)
  D4 <- run1a$df %>% filter(Age==15)
  D4 <- D4 %>% mutate(N2=N/sum(N))
  xout <- sort(c(run1$LenBins, DF2$Length))
  prob1 <- calcprobR(tempLens=DF2$Length, tempNs=DF2$N2, xout, LenBins=run1$LenBins)
  prob2 <- calcprobR(tempLens=D4$Length, tempNs=D4$N2, xout, LenBins=run1$LenBins)
  
  
  bind1 <- which(run1$LenBins == b1)
  bind2 <- which(run1$LenBins == b2)
  
  makeTransparent = function(..., alpha=0.5) {
    if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
    alpha = floor(255*alpha)  
    newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
    
    .makeTransparent = function(col, alpha) {
      rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
    }
    newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
    return(newColor)
  }
  
  par(mar=c(4,2,0,0))
  col1 <- makeTransparent('grey', alpha=0.8)
  col2 <- makeTransparent('blue', alpha=0.2)
  labs <- run1$LenMids[bind1:bind2]
  maxy <- max(c(prob1, prob2))
  x <- barplot(prob2[bind1:bind2], 
               ylim=c(0, maxy), col=col1, las=2)
  text(cex=1, x=x, y=-0.03, labs, xpd=TRUE, srt=90)
  barplot(prob1[bind1:bind2], add=TRUE, col=col2, axes=FALSE)
  mtext(side=1, "Length Classes", line=2.2, xpd=NA)
  mtext(side=2, "Relative Frequency", line=yline)
  text(x[1], 0.225, 'c)', xpd=NA)
  if(save) dev.off()
  
}


# Figure1 <- function(ngtg, binwidth, height=4, width=6, res=400, save=TRUE) {
#   if (save) png('Figures/Figure1.png', units='in', height=height, width=width,
#                 res=res)
# 
#   
#   ngtg_low <- ngtg
#   ngtg_high <- 1001
#   age <- 15
#   yr <- c(1,68)
# 
#   par(mfcol=c(2,3), oma=c(4,4,0,0), mar=c(1,1,1,1))
#   pch1 <- 21
#   pch2 <- 24
#   col1 <- 'black'
#   col2 <- 'darkgray'
#   xlab <- "Length class"
#   ylab <- "Relative frequency"
#   xline <- 2
#   yline <- 2
#   m.cex <- 1
#   yaxs <- "i"
#   xaxs <- "r"
#   xpos <- function(df) min(df$Length)
#   ypos <- function(maxy, ind=1) maxy$maxy[ind]
#   tex.cex <- 0.8
# 
#   # Life-history parameters #
#   M <- 0.2
#   Linf <- 100
#   LinfCV <- 0.1
#   K <- 0.13333
#   t0 <- 0
#   sigmaR <- 0
#   steepness <- 0.99
# 
#   L50 <- 66 # maturity
#   L95 <- 70
# 
#   # Selectivity
#   L5 <- 40
#   LFS <- 50
#   Vmaxlen <- 1
# 
#   maxsd <- 3
# 
#   # Trend in fishing mortality
#   curF <- 0.4
#   yr.st <- 1950
#   yr.end <- 2017
#   F.trend <- Ftrend(yr.st, yr.end, curF, 'stable', plot=FALSE)
#   nyrs <- length(F.trend)
# 
#   # Low Res GTG #
#   SimPop <- GTGpopsim(Linf, K, t0, M, L50, L95, LFS, L5, Vmaxlen, sigmaR,
#                       steepness, F.trend, ngtg=ngtg_low, binwidth = binwidth,
#                       maxsd=maxsd)
# 
#   LenMids <- SimPop$LenMids
#   LenBins <- SimPop$LenBins
#   t1 <- SimPop$df %>% filter(Yr %in% yr, Age==age)
#   LenComp <- t1 %>% group_by(Yr, Age) %>% mutate(Bin=cut(Length, LenBins)) %>%
#     group_by(Yr, Age, Bin) %>% summarize(P=sum(N)) %>% ungroup() %>%
#     group_by(Yr, Age) %>% mutate(Prob = P/sum(P), Bin = as.character(Bin))
# 
#   n1 <- gsub("\\(", "", x=LenComp$Bin) %>% strsplit(",") %>%
#     data.frame(stringsAsFactors=FALSE)
#   Bin1 <- n1[1,] %>% as.numeric()
#   Bin2 <- gsub("\\]", "", x=n1[2,]) %>% as.numeric()
# 
#   ProbDF <- data.frame(Yr=LenComp$Yr, Prob=LenComp$Prob, Bin1=Bin1, Bin2=Bin2)
#   maxy <- ProbDF %>% group_by(Yr) %>% summarize(maxy=max(Prob))
# 
#   unfished <- t1 %>% filter(Yr==yr[1])
#   fished <- t1 %>% filter(Yr==yr[2])
#   plot(range(unfished$Length), c(0,1), ylim=c(0, max(ProbDF$Prob)), type="n",
#        bty="l", xlab="", ylab="", axes=FALSE, yaxs=yaxs, xaxs=xaxs)
#   axis(side=1, labels=FALSE)
#   axis(side=2)
#   addLines(LenBins)
#   points(unfished$Length, unfished$N/max(unfished$N)*maxy$maxy[1], pch=pch1,
#          col="black", bg=col1, xpd=NA)
#   points(fished$Length, fished$N/max(fished$N)*maxy$maxy[2], pch=pch2,
#          col='black', bg=col2, xpd=NA)
#   text(xpos(fished), ypos(maxy), "a)", cex=tex.cex, xpd=NA)
# 
#   plot(range(unfished$Length), c(0,1), ylim=c(0, maxy$maxy[2]), type="n",
#        bty="l", xlab="", ylab="", axes=FALSE, yaxs=yaxs, xaxs=xaxs)
#   axis(side=1)
#   axis(side=2)
#   addLines(LenBins)
#   addBars(ProbDF, col2)
#   text(xpos(fished), ypos(maxy,2), "b)", cex=tex.cex, xpd=NA)
#   mtext(side=2, line=yline, ylab, outer=TRUE, cex=m.cex, xpd=NA)
# 
#   # High Res GTG #
#   SimPop <- GTGpopsim(Linf, K, t0, M, L50, L95, LFS, L5, Vmaxlen, sigmaR,
#                       steepness, F.trend, ngtg=ngtg_high, binwidth = binwidth,
#                       maxsd=maxsd)
# 
#   LenMids <- SimPop$LenMids
#   LenBins <- SimPop$LenBins
#   t1 <- SimPop$df %>% filter(Yr %in% yr, Age==age)
#   LenComp <- t1 %>% group_by(Yr, Age) %>% mutate(Bin=cut(Length, LenBins)) %>%
#     group_by(Yr, Age, Bin) %>% summarize(P=sum(N)) %>% ungroup() %>%
#     group_by(Yr, Age) %>% mutate(Prob = P/sum(P), Bin = as.character(Bin))
# 
#   n1 <- gsub("\\(", "", x=LenComp$Bin) %>% strsplit(",") %>%
#     data.frame(stringsAsFactors=FALSE)
#   Bin1 <- n1[1,] %>% as.numeric()
#   Bin2 <- gsub("\\]", "", x=n1[2,]) %>% as.numeric()
# 
#   ProbDF <- data.frame(Yr=LenComp$Yr, Prob=LenComp$Prob, Bin1=Bin1, Bin2=Bin2)
#   maxy <- ProbDF %>% group_by(Yr) %>% summarize(maxy=max(Prob))
# 
#   unfished <- t1 %>% filter(Yr==yr[1])
#   fished <- t1 %>% filter(Yr==yr[2])
#   plot(range(unfished$Length), c(0,1), ylim=c(0, max(ProbDF$Prob)), type="n",
#        bty="l", xlab="", ylab="", axes=FALSE, yaxs=yaxs, xaxs=xaxs)
#   axis(side=1, labels=FALSE)
#   axis(side=2, labels=FALSE)
#   addLines(LenBins)
#   points(unfished$Length, unfished$N/max(unfished$N)*maxy$maxy[1], pch=pch1,
#          col="black", bg=col1, xpd=NA)
#   points(fished$Length, fished$N/max(fished$N)*maxy$maxy[2], pch=pch2,
#          col='black', bg=col2, xpd=NA)
#   mtext(side=1, line=xline, xlab, outer=TRUE, cex=m.cex, xpd=NA)
#   text(xpos(fished), ypos(maxy), "c)", cex=tex.cex, xpd=NA)
# 
# 
#   plot(range(unfished$Length), c(0,1), ylim=c(0, maxy$maxy[2]), type="n",
#        bty="l", xlab="", ylab="", axes=FALSE, yaxs=yaxs, xaxs=xaxs)
#   axis(side=1)
#   axis(side=2, labels=FALSE)
#   addLines(LenBins)
#   addBars(ProbDF, col2)
#   text(xpos(fished), ypos(maxy,2), "d)", cex=tex.cex, xpd=NA)
# 
# 
#   # Approx Method Res GTG
#   SimPop <- GTGpopsim(Linf, K, t0, M, L50, L95, LFS, L5, Vmaxlen, sigmaR,
#                       steepness, F.trend, ngtg=ngtg_low, binwidth = binwidth,
#                       maxsd=maxsd)
#   LenMids <- SimPop$LenMids
#   LenBins <- SimPop$LenBins
# 
#   t1 <- SimPop$df %>% filter(Yr %in% yr, Age==age)
#   xout <- sort(c(LenBins, t1$Length %>% unique()))
# 
#   # unfished & unfished
#   Problist <- list()
#   for (i in seq_along(yr)) {
#     tempdf <- t1 %>% filter(Yr==yr[i])
#     yout <- linear_int(tempdf$Length, tempdf$N, xout)
#     Prob <- calcprob(x=tempdf$Length, y=tempdf$N, xout=xout, LenBins)
#     Problist[[i]] <- data.frame(Yr=yr[i], Prob=Prob,
#                                 Bin1=LenBins[1:(length(LenBins)-1)],
#                                 Bin2=LenBins[2:(length(LenBins))])
#   }
#   ProbDF <- do.call("rbind", Problist)
#   maxy <- ProbDF %>% group_by(Yr) %>% summarize(maxy=max(Prob))
# 
#   unfished <- t1 %>% filter(Yr==yr[1])
#   fished <- t1 %>% filter(Yr==yr[2])
#   plot(range(unfished$Length), c(0,1), ylim=c(0, max(ProbDF$Prob)), type="n",
#        bty="l", xlab="", ylab="", axes=FALSE, yaxs=yaxs, xaxs=xaxs)
#   axis(side=1, labels=FALSE)
#   axis(side=2, labels=FALSE)
#   addLines(LenBins)
#   lines(unfished$Length, unfished$N/max(unfished$N)*maxy$maxy[1], pch=pch1,
#         col="black", bg=col1, type="b", xpd=NA)
#   lines(fished$Length, fished$N/max(fished$N)*maxy$maxy[2], pch=pch2,
#         col='black', bg=col2, type="b", xpd=NA)
#   text(xpos(fished), ypos(maxy), "e)", cex=tex.cex, xpd=NA)
# 
#   plot(range(unfished$Length), c(0,1), ylim=c(0, maxy$maxy[2]), type="n",
#        bty="l", xlab="", ylab="", axes=FALSE, yaxs=yaxs, xaxs=xaxs)
#   axis(side=1)
#   axis(side=2, labels=FALSE)
#   addLines(LenBins)
#   addBars(ProbDF, col2)
#   text(xpos(fished), ypos(maxy,2), "f)", cex=tex.cex, xpd=NA)
# 
#   if(save) dev.off()
# }

# 
# Figure1 <- function(x, OM_DF, ngtg, binwidth, Fmulti=1, age=10,
#                     height=4, width=6, res=400, save=TRUE) {
#   if (save) png('Figures/Figure1.png', units='in', height=height, width=width,
#                 res=res)
#   
#   binwidth <- binwidth
#   ngtg_low <- ngtg
#   ngtg_high <- 1001
#   
#   
#   par(mfcol=c(2,3), oma=c(4,4,0,0), mar=c(1,1,1,1))
#   pch1 <- 21
#   pch2 <- 24
#   col1 <- 'black'
#   col2 <- 'darkgray'
#   xlab <- "Length class"
#   ylab <- "Relative frequency"
#   xline <- 2
#   yline <- 2
#   m.cex <- 1
#   yaxs <- "i"
#   xaxs <- "r"
#   xpos <- function(df) min(df$Length)
#   ypos <- function(maxy, ind=1) maxy$maxy[ind]
#   tex.cex <- 0.8
#   
#  
#   
#   # Low Res GTG #
#   SimPop <- EqGTG(x, OM_DF, ngtg=ngtg, Fmulti=0, binwidth = binwidth)
#   SimPop2 <- EqGTG(x, OM_DF, ngtg=ngtg, Fmulti=Fmulti, binwidth = binwidth)
#   
#   LenMids <- SimPop$LenMids
#   LenBins <- SimPop$LenBins
#   t1 <- SimPop$df %>% filter(Age==age)
#   LenComp <- t1 %>% group_by(Age) %>% mutate(Bin=cut(Length, LenBins)) %>%
#     group_by(Age, Bin) %>% summarize(P=sum(N)) %>% ungroup() %>%
#     group_by(Age) %>% mutate(Prob = P/sum(P), Bin = as.character(Bin))
#   
#   n1 <- gsub("\\(", "", x=LenComp$Bin) %>% strsplit(",") %>%
#     data.frame(stringsAsFactors=FALSE)
#   Bin1 <- n1[1,] %>% as.numeric()
#   Bin2 <- gsub("\\]", "", x=n1[2,]) %>% as.numeric()
#   
#   ProbDF <- data.frame(Prob=LenComp$Prob, Bin1=Bin1, Bin2=Bin2)
#   maxy <- ProbDF%>% summarize(maxy=max(Prob))
#   
#   unfished <- SimPop[[1]] %>% filter(Age==age)
#   fished <- SimPop2[[1]] %>% filter(Age==age)
# 
#   plot(c(min(Bin1), max(Bin2)), c(0,1), ylim=c(0, max(ProbDF$Prob)), type="n",
#        bty="l", xlab="", ylab="", axes=FALSE, yaxs=yaxs, xaxs=xaxs)
#   axis(side=1, labels=FALSE)
#   axis(side=2)
#   addLines(LenBins)
#   points(unfished$Length, unfished$N/sum(unfished$N), pch=pch1,
#          col="black", bg=col1, xpd=NA)
#   points(fished$Length, fished$N/sum(fished$N), pch=pch2,
#          col='black', bg=col2, xpd=NA)
#   text(xpos(fished), ypos(maxy), "a)", cex=tex.cex, xpd=NA)
#   
#   plot(c(min(Bin1), max(Bin2)), c(0,1), ylim=c(0, maxy$maxy[1]), type="n",
#        bty="l", xlab="", ylab="", axes=FALSE, yaxs=yaxs, xaxs=xaxs)
#   axis(side=1)
#   axis(side=2)
#   addLines(LenBins)
#   addBars(ProbDF, col2)
#   text(xpos(fished), ypos(maxy,2), "b)", cex=tex.cex, xpd=NA)
#   mtext(side=2, line=yline, ylab, outer=TRUE, cex=m.cex, xpd=NA)
#   
#   
#   # High Res GTG #
#   SimPop <- EqGTG(x, OM_DF, ngtg=ngtg_high, Fmulti=0, binwidth = binwidth)
#   SimPop2 <- EqGTG(x, OM_DF, ngtg=ngtg_high, Fmulti=Fmulti, binwidth = binwidth)
#   LenMids <- SimPop$LenMids
#   LenBins <- SimPop$LenBins
#   t1 <- SimPop$df %>% filter(Age==age)
#   LenComp <- t1 %>% group_by(Age) %>% mutate(Bin=cut(Length, LenBins)) %>%
#     group_by(Age, Bin) %>% summarize(P=sum(N)) %>% ungroup() %>%
#     group_by(Age) %>% mutate(Prob = P/sum(P), Bin = as.character(Bin))
#   
#   n1 <- gsub("\\(", "", x=LenComp$Bin) %>% strsplit(",") %>%
#     data.frame(stringsAsFactors=FALSE)
#   Bin1 <- n1[1,] %>% as.numeric()
#   Bin2 <- gsub("\\]", "", x=n1[2,]) %>% as.numeric()
#   
#   ProbDF <- data.frame(Prob=LenComp$Prob, Bin1=Bin1, Bin2=Bin2)
#   # ProbDF$Prob <- ProbDF$Prob/max(ProbDF$Prob)
#   
#   unfished <- SimPop[[1]] %>% filter(Age==age)
#   yvals <- unfished$N/sum(unfished$N)
#   fished <- SimPop2[[1]] %>% filter(Age==age)
#   
#   plot(c(min(Bin1), max(Bin2)), c(0,max(yvals)), ylim=c(0, max(yvals)), type="n",
#        bty="l", xlab="", ylab="", axes=FALSE, yaxs=yaxs, xaxs=xaxs)
#   axis(side=1, labels=FALSE)
#   axis(side=2)
#   addLines(LenBins)
#   
#   # yvals <- yvals/max(yvals) 
#   points(unfished$Length, yvals, pch=pch1,
#          col="black", bg=col1, xpd=NA)
#   yvals2<- fished$N/sum(fished$N)
#   # yvals2 <- yvals2/sum(yvals2) 
#   points(fished$Length, yvals2, pch=pch2,
#          col='black', bg=col2, xpd=NA)
#   text(xpos(fished), ypos(maxy), "a)", cex=tex.cex, xpd=NA)
#   
#   plot(c(min(Bin1), max(Bin2)), c(0,max(ProbDF$Prob)), 
#        ylim=c(0, max(ProbDF$Prob)), type="n",
#        bty="l", xlab="", ylab="", axes=FALSE, yaxs=yaxs, xaxs=xaxs)
#   axis(side=1)
#   axis(side=2)
#   addLines(LenBins)
#   addBars(ProbDF, col2)
#   text(xpos(fished), ypos(maxy,2), "b)", cex=tex.cex, xpd=NA)
#   mtext(side=2, line=yline, ylab, outer=TRUE, cex=m.cex, xpd=NA)
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   LenMids <- SimPop$LenMids
#   LenBins <- SimPop$LenBins
#   t1 <- SimPop$df %>% filter(Yr %in% yr, Age==age)
#   LenComp <- t1 %>% group_by(Yr, Age) %>% mutate(Bin=cut(Length, LenBins)) %>%
#     group_by(Yr, Age, Bin) %>% summarize(P=sum(N)) %>% ungroup() %>%
#     group_by(Yr, Age) %>% mutate(Prob = P/sum(P), Bin = as.character(Bin))
#   
#   n1 <- gsub("\\(", "", x=LenComp$Bin) %>% strsplit(",") %>%
#     data.frame(stringsAsFactors=FALSE)
#   Bin1 <- n1[1,] %>% as.numeric()
#   Bin2 <- gsub("\\]", "", x=n1[2,]) %>% as.numeric()
#   
#   ProbDF <- data.frame(Yr=LenComp$Yr, Prob=LenComp$Prob, Bin1=Bin1, Bin2=Bin2)
#   maxy <- ProbDF %>% group_by(Yr) %>% summarize(maxy=max(Prob))
#   
#   unfished <- t1 %>% filter(Yr==yr[1])
#   fished <- t1 %>% filter(Yr==yr[2])
#   plot(range(unfished$Length), c(0,1), ylim=c(0, max(ProbDF$Prob)), type="n",
#        bty="l", xlab="", ylab="", axes=FALSE, yaxs=yaxs, xaxs=xaxs)
#   axis(side=1, labels=FALSE)
#   axis(side=2, labels=FALSE)
#   addLines(LenBins)
#   points(unfished$Length, unfished$N/max(unfished$N)*maxy$maxy[1], pch=pch1,
#          col="black", bg=col1, xpd=NA)
#   points(fished$Length, fished$N/max(fished$N)*maxy$maxy[2], pch=pch2,
#          col='black', bg=col2, xpd=NA)
#   mtext(side=1, line=xline, xlab, outer=TRUE, cex=m.cex, xpd=NA)
#   text(xpos(fished), ypos(maxy), "c)", cex=tex.cex, xpd=NA)
#   
#   
#   plot(range(unfished$Length), c(0,1), ylim=c(0, maxy$maxy[2]), type="n",
#        bty="l", xlab="", ylab="", axes=FALSE, yaxs=yaxs, xaxs=xaxs)
#   axis(side=1)
#   axis(side=2, labels=FALSE)
#   addLines(LenBins)
#   addBars(ProbDF, col2)
#   text(xpos(fished), ypos(maxy,2), "d)", cex=tex.cex, xpd=NA)
#   
#   
#   # Approx Method Res GTG
#   SimPop <- GTGpopsim(Linf, K, t0, M, L50, L95, LFS, L5, Vmaxlen, sigmaR,
#                       steepness, F.trend, ngtg=ngtg_low, binwidth = binwidth,
#                       maxsd=maxsd)
#   LenMids <- SimPop$LenMids
#   LenBins <- SimPop$LenBins
#   
#   t1 <- SimPop$df %>% filter(Yr %in% yr, Age==age)
#   xout <- sort(c(LenBins, t1$Length %>% unique()))
#   
#   # unfished & unfished
#   Problist <- list()
#   for (i in seq_along(yr)) {
#     tempdf <- t1 %>% filter(Yr==yr[i])
#     yout <- linear_int(tempdf$Length, tempdf$N, xout)
#     Prob <- calcprob(x=tempdf$Length, y=tempdf$N, xout=xout, LenBins)
#     Problist[[i]] <- data.frame(Yr=yr[i], Prob=Prob,
#                                 Bin1=LenBins[1:(length(LenBins)-1)],
#                                 Bin2=LenBins[2:(length(LenBins))])
#   }
#   ProbDF <- do.call("rbind", Problist)
#   maxy <- ProbDF %>% group_by(Yr) %>% summarize(maxy=max(Prob))
#   
#   unfished <- t1 %>% filter(Yr==yr[1])
#   fished <- t1 %>% filter(Yr==yr[2])
#   plot(range(unfished$Length), c(0,1), ylim=c(0, max(ProbDF$Prob)), type="n",
#        bty="l", xlab="", ylab="", axes=FALSE, yaxs=yaxs, xaxs=xaxs)
#   axis(side=1, labels=FALSE)
#   axis(side=2, labels=FALSE)
#   addLines(LenBins)
#   lines(unfished$Length, unfished$N/max(unfished$N)*maxy$maxy[1], pch=pch1,
#         col="black", bg=col1, type="b", xpd=NA)
#   lines(fished$Length, fished$N/max(fished$N)*maxy$maxy[2], pch=pch2,
#         col='black', bg=col2, type="b", xpd=NA)
#   text(xpos(fished), ypos(maxy), "e)", cex=tex.cex, xpd=NA)
#   
#   plot(range(unfished$Length), c(0,1), ylim=c(0, maxy$maxy[2]), type="n",
#        bty="l", xlab="", ylab="", axes=FALSE, yaxs=yaxs, xaxs=xaxs)
#   axis(side=1)
#   axis(side=2, labels=FALSE)
#   addLines(LenBins)
#   addBars(ProbDF, col2)
#   text(xpos(fished), ypos(maxy,2), "f)", cex=tex.cex, xpd=NA)
#   
#   if(save) dev.off()
# }
# 
