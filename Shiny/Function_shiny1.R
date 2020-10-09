
# functions ------------------------------------------------------------------


get_RR = function(h2, K, a){
  T    = -qnorm(K, 0,1)
  z    = dnorm(T)
  i    = z/K
  Tr   = (T - a * i * h2) / (sqrt(1 - a * a * h2 * h2 * i * (i-T)))
  Kr   = 1 - pnorm(Tr)
  RR   = Kr / K
  Tr2   = (T -i*h2) / sqrt(1-(0.5*h2*h2*i)*(i-T))
  Kr2   = 1 - pnorm(Tr2)
  RR2   = Kr2 / K
  return(list("RRabs"=round(Kr,3), "RRrel"=round(RR,3), "K"=round(K,3), "RRabs2"=round(Kr2,3), "h2"=round(h2,3), "Tr"=round(Tr,3),"Kr"=round(Kr,3), "T"=round(T,3), "i"=round(i,3)))
}

est_h2 = function(K, Kr, a){
  T  = -qnorm(K, 0, 1)
  Tr = -qnorm(Kr, 0, 1)
  z  = dnorm(T)
  i  = z / K
  h2 = (T - Tr * sqrt(1-(1-T/i)*(T^2-Tr^2))) / (a * i+(i-T)*Tr^2)
  return(list("h2"=h2, "T"=T, "Tr"=Tr, "z"=z, "i"=i))
}

get_CDRR = function(h2x, h2y, Kx, Ky, rg, a){
  Tx = -qnorm(Kx, 0, 1) 
  Ty = -qnorm(Ky, 0, 1) 
  zx = dnorm(Tx) 
  zy = dnorm(Ty) 
  ix = zx / Kx 
  iy = zy / Ky
  #risk of disorder x in relatives of those with disorder x
  Txx   = (Tx - a * ix * h2x) / (sqrt(1 - a * a * h2x * h2x * ix * (ix-Tx)))
  Kxx   = 1 - pnorm(Txx)
  RRxx   = Kxx / Kx
  #risk of disorder y in relatives of those with disorder y
  Tyy   = (Ty - a * iy * h2y) / (sqrt(1 - a * a * h2y * h2y * iy * (iy-Ty)))
  Kyy   = 1 - pnorm(Tyy)
  RRyy   = Kyy / Ky
  # calculate the genetic covariance:
  covg = rg * sqrt(h2x * h2y)
  # calculate risk ratio for relatives
  Txy  = (Tx - a * covg * iy) / sqrt(1 - a^2 * covg^2 * iy * (iy-Ty))
  Kxy  = 1 - pnorm(Txy)
  RRxy = Kxy/Kx # relative risk of disease x given that parent has disease y
  # calculate risk ratio for relatives other way around
  Tyx  = (Ty - a * covg * ix) / sqrt(1 - a^2 * covg^2 * ix * (ix-Tx))
  Kyx  = 1 - pnorm(Tyx)
  RRyx = Kyx/Ky # relative risk of disease x given that parent has disease y
  return(list("Kx"=Kx,"Kxx"=Kxx,"RRxx"=RRxx,"Ky"=Ky, "Kyy"=Kyy,"RRyy"=RRyy, "Kxy"=Kxy, 
              "RRxy"=RRxy, "Kyx"=Kyx, "RRyx"=RRyx, "rg"=rg, "covg"=covg, "h2x"=h2x, "h2y"=h2y ))  
}



# Plot Functions ----------------------------------------------------------


mk_plot1 <- function(K,h2, xlab=NA, ylab=NA){
  
  Txp = -qnorm(K, 0,1)
  z = dnorm(Txp) 
  i = z / K[1]
  mean = -(Txp) ;sd=1
  lb=0; ub=4
  
  x <- seq(Txp-6,Txp+3,length=1000)*sd + mean
  hx <- dnorm(x,mean,sd)
  K <- round(K,2)
  K_prevalence <- c(as.expression(bquote(italic(K)[x]~'='~.(K))))
  h2=h2
  p2 <- plot(x, hx, type="n", xlab="", ylab="",
             main=c(as.expression((bquote(bold("Phenotype x:"~h^2~"="~ .(h2)))))), cex.main =2, axes=FALSE)
  i <- x >= lb & x <= ub
  l <- x <= lb & x <= ub
  lines(x, hx)
  polygon(c(lb,x[i],ub), c(0,hx[i],0), col="darkblue")
  polygon(c(-5,x[l],lb), c(0,hx[l],0), col="grey") 
  
  
  mtext(K_prevalence, side = 1)
  
  abline(h=2)
  abline(v=0, h=-2)
  
  return(p2)  
}  

mk_plot1.1par <- function(K, xlab=NA, ylab=NA){
  
  Txp = -qnorm(K, 0,1)
  z = dnorm(Txp) 
  i = z / K[1]
  mean = -(Txp) ;sd=1
  lb=0; ub=4
  
  x <- seq(Txp-6,Txp+3,length=1000)*sd + mean
  hx <- dnorm(x,mean,sd)
  K <- round(K,2)
  #K_prevalence <- c(as.expression(bquote(italic(K)~'='~.(K))))
  K_prevalence <- c(as.expression(bquote(italic(K)[relative(aR)]~'='~.(K))))
  
  p2 <- plot(x, hx, type="n", xlab="", ylab="",
             main="", axes=FALSE)
  i <- x >= lb & x <= ub
  l <- x <= lb & x <= ub
  lines(x, hx)
  polygon(c(lb,x[i],ub), c(0,hx[i],0), col="darkblue")
  polygon(c(-5,x[l],lb), c(0,hx[l],0), col="grey") 
  
  
  mtext(K_prevalence, side = 1)
  
  abline(h=2)
  abline(v=0, h=-2)
  
  return(p2)  
}  

mk_plot1.2par <- function(K, xlab=NA, ylab=NA){
  
  Txp = -qnorm(K, 0,1)
  z = dnorm(Txp) 
  i = z / K[1]
  mean = -(Txp) ;sd=1
  lb=0; ub=4
  
  x <- seq(Txp-6,Txp+3,length=1000)*sd + mean
  hx <- dnorm(x,mean,sd)
  K <- round(K,2)
  #K_prevalence <- c(as.expression(bquote(italic(K)~'='~.(K))))
  K_prevalence <- c(as.expression(bquote(italic(K)[2][par]~'='~.(K))))
  
  p2 <- plot(x, hx, type="n", xlab="", ylab="",
             main="", axes=FALSE)
  i <- x >= lb & x <= ub
  l <- x <= lb & x <= ub
  lines(x, hx)
  polygon(c(lb,x[i],ub), c(0,hx[i],0), col="darkblue")
  polygon(c(-5,x[l],lb), c(0,hx[l],0), col="grey") 
  
  
  mtext(K_prevalence, side = 1)
  
  abline(h=2)
  abline(v=0, h=-2)
  
  return(p2)  
}  

mk_plot1.1 <- function(K, h2){
  
  Txp = -qnorm(K, 0,1)
  z = dnorm(Txp) 
  i = z / K[1]
  mean = -(Txp) ;sd=1
  lb=0; ub=4
  
  x <- seq(Txp-6,Txp+3,length=1000)*sd + mean
  hx <- dnorm(x,mean,sd)
  K <- round(K,2)
  K_prevalence <- c(as.expression(bquote(italic(K)[y]~'='~.(K))))
  h2=h2
  p2 <- plot(x, hx, type="n", xlab="", ylab="",
             main=c(as.expression((bquote(bold("Phenotype y:"~h^2~"="~ .(h2)))))), cex.main =2, axes=FALSE)
  i <- x >= lb & x <= ub
  l <- x <= lb & x <= ub
  lines(x, hx)
  polygon(c(lb,x[i],ub), c(0,hx[i],0), col="orange")
  polygon(c(-5,x[l],lb), c(0,hx[l],0), col="grey") 
  
  
  mtext(K_prevalence, side = 1)
  
  abline(h=2)
  abline(v=0, h=-2)
  
  return(p2)  
}  

mk_plot1.1.1par <- function(K, xlab=NA, ylab=NA){
  
  Txp = -qnorm(K, 0,1)
  z = dnorm(Txp) 
  i = z / K[1]
  mean = -(Txp) ;sd=1
  lb=0; ub=4
  
  x <- seq(Txp-6,Txp+3,length=1000)*sd + mean
  hx <- dnorm(x,mean,sd)
  K <- round(K,2)
  #K_prevalence <- c(as.expression(bquote(italic(K)~'='~.(K))))
  K_prevalence <- c(as.expression(bquote(italic(K)[relative (aR)]~'='~.(K))))
  
  p2 <- plot(x, hx, type="n", xlab="", ylab="",
             main="", axes=FALSE)
  i <- x >= lb & x <= ub
  l <- x <= lb & x <= ub
  lines(x, hx)
  polygon(c(lb,x[i],ub), c(0,hx[i],0), col="orange")
  polygon(c(-5,x[l],lb), c(0,hx[l],0), col="grey") 
  
  mtext(K_prevalence, side = 1)
  
  abline(h=2)
  abline(v=0, h=-2)
  
  return(p2)  
}  

mk_plot1.1.2par <- function(K, xlab=NA, ylab=NA){
  
  Txp = -qnorm(K, 0,1)
  z = dnorm(Txp) 
  i = z / K[1]
  mean = -(Txp) ;sd=1
  lb=0; ub=4
  
  x <- seq(Txp-6,Txp+3,length=1000)*sd + mean
  hx <- dnorm(x,mean,sd)
  K <- round(K,2)
  #K_prevalence <- c(as.expression(bquote(italic(K)~'='~.(K))))
  K_prevalence <- c(as.expression(bquote(italic(K)[2][par]~'='~.(K))))
  
  p2 <- plot(x, hx, type="n", xlab="", ylab="",
             main="", axes=FALSE)
  i <- x >= lb & x <= ub
  l <- x <= lb & x <= ub
  lines(x, hx)
  polygon(c(lb,x[i],ub), c(0,hx[i],0), col="orange")
  polygon(c(-5,x[l],lb), c(0,hx[l],0), col="grey") 
  
  
  mtext(K_prevalence, side = 1)
  
  abline(h=2)
  abline(v=0, h=-2)
  
  return(p2)  
}  


mk_plot2 <- function(K, xlab="NA", ylab="NA"){
  Prevalence <- round(K*100,0)
  par(mar=c(5.1,0,4.1,0))
  
  result <- matrix(nrow = 1, ncol = 100)
  result <- sample(as.matrix(c(rep("b",Prevalence),rep("a",100-Prevalence))))
  
  cols <- c("green", "red")
  nr <- 10
  nc <- 10
  
  # create data.frame of positions and colors
  m <- matrix(cols[factor(result)], nr, nc)
  DF <- data.frame(row = c(row(m)), col = c(col(m)[, nc:1]), value = c(m), gender =  sample(c(rep("male",50),rep("female",50))),
                   stringsAsFactors = FALSE)    
  
  plot(col ~ row, DF, col = DF$value, asp = 1, cex=5,
       xlim = c(0, nr), ylim = c(0, nc),
       axes = FALSE, xlab = "", ylab = "", type = "n")
  
  library(png)
  #setwd("~/Desktop")
  man.grey <- readPNG("male.dark.grey.png")
  man.blue <- readPNG("male.dark.blue.png")
  female.grey <- readPNG("female.dark.grey.png")
  female.blue <- readPNG("female_dark.blue.png")
  
  #blue
  Rm <- subset(DF, value == "red" & gender=="male")
  with(Rm, rasterImage(man.blue, 
                       row-.465, col-.465, row+.465, col+.465, 
                       xlim = c(0, nr), ylim = c(0, nc),
                       xlab = "", ylab = ""))
  
  #blue
  Rf <- subset(DF, value == "red" & gender=="female")
  with(Rf, rasterImage(female.blue, 
                       row-.465, col-.465, row+.465, col+.465, 
                       xlim = c(0, nr), ylim = c(0, nc),
                       xlab = "", ylab = ""))
  
  
  #grey.man <- man
  Gm <- subset(DF, value == "green" & gender == "male")
  with(Gm, rasterImage(man.grey, 
                       row-.465, col-.465, row+.465, col+.465, 
                       xlim = c(0, nr), ylim = c(0, nc),
                       xlab = "", ylab = ""))
  
  
  
  #grey.female <- female
  #green.female[,,2] <- female[,,4] # fill in green dimension
  Gf <- subset(DF, value == "green" & gender == "female")
  with(Gf, rasterImage(female.grey, 
                       row-.465, col-.465, row+.465, col+.465, 
                       xlim = c(0, nr), ylim = c(0, nc),
                       xlab = "", ylab = ""))
  
  
  mtext("Sample drawn from \n general population", side = 1)
}


mk_plot2.2 <- function(K, xlab="NA", ylab="NA"){
  Prevalence <- round(K*100,0)
  par(mar=c(5.1,0,4.1,0))
  
  result <- matrix(nrow = 1, ncol = 100)
  result <- sample(as.matrix(c(rep("b",Prevalence),rep("a",100-Prevalence))))
  
  cols <- c("green", "red")
  nr <- 10
  nc <- 10
  
  # create data.frame of positions and colors
  m <- matrix(cols[factor(result)], nr, nc)
  DF <- data.frame(row = c(row(m)), col = c(col(m)[, nc:1]), value = c(m), gender =  sample(c(rep("male",50),rep("female",50))),
                   stringsAsFactors = FALSE)    
  
  plot(col ~ row, DF, col = DF$value, asp = 1, cex=5,
       xlim = c(0, nr), ylim = c(0, nc),
       axes = FALSE, xlab = "", ylab = "", type = "n")
  
  library(png)
  #setwd("~/Desktop")
  man.grey <- readPNG("male.dark.grey.png")
  man.blue <- readPNG("male.orange.png")
  female.grey <- readPNG("female.dark.grey.png")
  female.blue <- readPNG("female.orange.png")
  
  #blue
  Rm <- subset(DF, value == "red" & gender=="male")
  with(Rm, rasterImage(man.blue, 
                       row-.465, col-.465, row+.465, col+.465, 
                       xlim = c(0, nr), ylim = c(0, nc),
                       xlab = "", ylab = ""))
  
  #blue
  Rf <- subset(DF, value == "red" & gender=="female")
  with(Rf, rasterImage(female.blue, 
                       row-.465, col-.465, row+.465, col+.465, 
                       xlim = c(0, nr), ylim = c(0, nc),
                       xlab = "", ylab = ""))
  
  
  #grey.man <- man
  Gm <- subset(DF, value == "green" & gender == "male")
  with(Gm, rasterImage(man.grey, 
                       row-.465, col-.465, row+.465, col+.465, 
                       xlim = c(0, nr), ylim = c(0, nc),
                       xlab = "", ylab = ""))
  
  
  
  #grey.female <- female
  #green.female[,,2] <- female[,,4] # fill in green dimension
  Gf <- subset(DF, value == "green" & gender == "female")
  with(Gf, rasterImage(female.grey, 
                       row-.465, col-.465, row+.465, col+.465, 
                       xlim = c(0, nr), ylim = c(0, nc),
                       xlab = "", ylab = ""))
  
  
  mtext("Sample drawn from \n general population", side = 1)
}

mk_plot3 <- function(K, Kr){
  Prevalence <- round(Kr*100,0)
  par(mar=c(5.1,0,4.1,0))
  
  result <- matrix(nrow = 1, ncol = 100)
  result <- sample(as.matrix(c(rep("b",Prevalence),rep("a",100-Prevalence))))
  
  cols <- c("green", "red")
  nr <- 10
  nc <- 10
  
  # create data.frame of positions and colors
  m <- matrix(cols[factor(result)], nr, nc)
  DF <- data.frame(row = c(row(m)), col = c(col(m)[, nc:1]), value = c(m), gender =  sample(c(rep("male",50),rep("female",50))),
                   stringsAsFactors = FALSE)    
  
  
  rr <- round(Kr/K,1)
  #rr <- round(1/2.1,1)
  RR <- c(as.expression(bquote(italic(RR)[x][y]~'='~.(rr))))
  
  plot(col ~ row, DF, col = DF$value, asp = 1, cex=5,
       xlim = c(0, nr), ylim = c(0, nc),
       axes = FALSE, xlab = "", ylab = "", type = "n")
  
  library(png)
  #setwd("~/Desktop")
  man.grey <- readPNG("male.dark.grey.png")
  man.blue <- readPNG("male.dark.blue.png")
  female.grey <- readPNG("female.dark.grey.png")
  female.blue <- readPNG("female_dark.blue.png")
  
  #blue
  Rm <- subset(DF, value == "red" & gender=="male")
  with(Rm, rasterImage(man.blue, 
                       row-.465, col-.465, row+.465, col+.465, 
                       xlim = c(0, nr), ylim = c(0, nc),
                       xlab = "", ylab = ""))
  
  #blue
  Rf <- subset(DF, value == "red" & gender=="female")
  with(Rf, rasterImage(female.blue, 
                       row-.465, col-.465, row+.465, col+.465, 
                       xlim = c(0, nr), ylim = c(0, nc),
                       xlab = "", ylab = ""))
  
  
  #grey.man <- man
  Gm <- subset(DF, value == "green" & gender == "male")
  with(Gm, rasterImage(man.grey, 
                       row-.465, col-.465, row+.465, col+.465, 
                       xlim = c(0, nr), ylim = c(0, nc),
                       xlab = "", ylab = ""))
  
  
  
  #grey.female <- female
  #green.female[,,2] <- female[,,4] # fill in green dimension
  Gf <- subset(DF, value == "green" & gender == "female")
  with(Gf, rasterImage(female.grey, 
                       row-.465, col-.465, row+.465, col+.465, 
                       xlim = c(0, nr), ylim = c(0, nc),
                       xlab = "", ylab = ""))
  
  mtext(RR, side = 3)
  mtext("Sample drawn from those \n related (aR) to one affected person", side = 1)
}

mk_plot3.1 <- function(K, Kr){
  Prevalence <- round(Kr*100,0)
  par(mar=c(5.1,0,4.1,0))
  
  result <- matrix(nrow = 1, ncol = 100)
  result <- sample(as.matrix(c(rep("b",Prevalence),rep("a",100-Prevalence))))
  
  cols <- c("green", "red")
  nr <- 10
  nc <- 10
  
  # create data.frame of positions and colors
  m <- matrix(cols[factor(result)], nr, nc)
  DF <- data.frame(row = c(row(m)), col = c(col(m)[, nc:1]), value = c(m), gender =  sample(c(rep("male",50),rep("female",50))),
                   stringsAsFactors = FALSE)    
  
  
  rr <- round(Kr/K,1)
  #rr <- round(1/2.1,1)
  RR <- c(as.expression(bquote(italic(RR)[y][x]~'='~.(rr))))
  
  plot(col ~ row, DF, col = DF$value, asp = 1, cex=5,
       xlim = c(0, nr), ylim = c(0, nc),
       axes = FALSE, xlab = "", ylab = "", type = "n")
  
  library(png)
  #setwd("~/Desktop")
  man.grey <- readPNG("male.dark.grey.png")
  man.blue <- readPNG("male.orange.png")
  female.grey <- readPNG("female.dark.grey.png")
  female.blue <- readPNG("female.orange.png")
  
  #blue
  Rm <- subset(DF, value == "red" & gender=="male")
  with(Rm, rasterImage(man.blue, 
                       row-.465, col-.465, row+.465, col+.465, 
                       xlim = c(0, nr), ylim = c(0, nc),
                       xlab = "", ylab = ""))
  
  #blue
  Rf <- subset(DF, value == "red" & gender=="female")
  with(Rf, rasterImage(female.blue, 
                       row-.465, col-.465, row+.465, col+.465, 
                       xlim = c(0, nr), ylim = c(0, nc),
                       xlab = "", ylab = ""))
  
  
  #grey.man <- man
  Gm <- subset(DF, value == "green" & gender == "male")
  with(Gm, rasterImage(man.grey, 
                       row-.465, col-.465, row+.465, col+.465, 
                       xlim = c(0, nr), ylim = c(0, nc),
                       xlab = "", ylab = ""))
  
  
  
  #grey.female <- female
  #green.female[,,2] <- female[,,4] # fill in green dimension
  Gf <- subset(DF, value == "green" & gender == "female")
  with(Gf, rasterImage(female.grey, 
                       row-.465, col-.465, row+.465, col+.465, 
                       xlim = c(0, nr), ylim = c(0, nc),
                       xlab = "", ylab = ""))
  
  mtext(RR, side = 3)
  mtext("Sample drawn from those \n related (aR) to one affected person", side = 1)
}

mk_plot4 <- function(K,Kr2, xlab="NA", ylab="NA"){
  Prevalence <- round(Kr2*100,0)
  par(mar=c(5.1,0,4.1,0))
  
  result <- matrix(nrow = 1, ncol = 100)
  result <- sample(as.matrix(c(rep("b",Prevalence),rep("a",100-Prevalence))))
  
  cols <- c("green", "red")
  nr <- 10
  nc <- 10
  
  # create data.frame of positions and colors
  m <- matrix(cols[factor(result)], nr, nc)
  DF <- data.frame(row = c(row(m)), col = c(col(m)[, nc:1]), value = c(m), gender =  sample(c(rep("male",50),rep("female",50))),
                   stringsAsFactors = FALSE)    
  
  
  rr2 <- round(Kr2/K,0)
  RR <- c(as.expression(bquote(italic(RR)[2][par]~'='~.(rr2))))
  
  plot(col ~ row, DF, col = DF$value, asp = 1, cex=5,
       xlim = c(0, nr), ylim = c(0, nc),
       axes = FALSE, xlab = "", ylab = "", type = "n")
  
  library(png)
  #setwd("~/Desktop")
  man.grey <- readPNG("male.dark.grey.png")
  man.blue <- readPNG("male.dark.blue.png")
  female.grey <- readPNG("female.dark.grey.png")
  female.blue <- readPNG("female_dark.blue.png")
  
  #blue
  Rm <- subset(DF, value == "red" & gender=="male")
  with(Rm, rasterImage(man.blue, 
                       row-.465, col-.465, row+.465, col+.465, 
                       xlim = c(0, nr), ylim = c(0, nc),
                       xlab = "", ylab = ""))
  
  #blue
  Rf <- subset(DF, value == "red" & gender=="female")
  with(Rf, rasterImage(female.blue, 
                       row-.465, col-.465, row+.465, col+.465, 
                       xlim = c(0, nr), ylim = c(0, nc),
                       xlab = "", ylab = ""))
  
  
  #grey.man <- man
  Gm <- subset(DF, value == "green" & gender == "male")
  with(Gm, rasterImage(man.grey, 
                       row-.465, col-.465, row+.465, col+.465, 
                       xlim = c(0, nr), ylim = c(0, nc),
                       xlab = "", ylab = ""))
  
  
  
  #grey.female <- female
  #green.female[,,2] <- female[,,4] # fill in green dimension
  Gf <- subset(DF, value == "green" & gender == "female")
  with(Gf, rasterImage(female.grey, 
                       row-.465, col-.465, row+.465, col+.465, 
                       xlim = c(0, nr), ylim = c(0, nc),
                       xlab = "", ylab = ""))
  
  mtext(RR, side = 3)
  mtext("Sample drawn from \n those with two affected parents", side = 1)
}

mk_plot4.1 <- function(K,Kr2, xlab="NA", ylab="NA"){
  Prevalence <- round(Kr2*100,0)
  par(mar=c(5.1,0,4.1,0))
  
  result <- matrix(nrow = 1, ncol = 100)
  result <- sample(as.matrix(c(rep("b",Prevalence),rep("a",100-Prevalence))))
  
  cols <- c("green", "red")
  nr <- 10
  nc <- 10
  
  # create data.frame of positions and colors
  m <- matrix(cols[factor(result)], nr, nc)
  DF <- data.frame(row = c(row(m)), col = c(col(m)[, nc:1]), value = c(m), gender =  sample(c(rep("male",50),rep("female",50))),
                   stringsAsFactors = FALSE)    
  
  
  rr2 <- round(Kr2/K,0)
  RR <- c(as.expression(bquote(italic(RR)[2][par]~'='~.(rr2))))
  
  plot(col ~ row, DF, col = DF$value, asp = 1, cex=5,
       xlim = c(0, nr), ylim = c(0, nc),
       axes = FALSE, xlab = "", ylab = "", type = "n")
  
  library(png)
  #setwd("~/Desktop")
  man.grey <- readPNG("male.dark.grey.png")
  man.blue <- readPNG("male.orange.png")
  female.grey <- readPNG("female.dark.grey.png")
  female.blue <- readPNG("female.orange.png")
  
  #blue
  Rm <- subset(DF, value == "red" & gender=="male")
  with(Rm, rasterImage(man.blue, 
                       row-.465, col-.465, row+.465, col+.465, 
                       xlim = c(0, nr), ylim = c(0, nc),
                       xlab = "", ylab = ""))
  
  #blue
  Rf <- subset(DF, value == "red" & gender=="female")
  with(Rf, rasterImage(female.blue, 
                       row-.465, col-.465, row+.465, col+.465, 
                       xlim = c(0, nr), ylim = c(0, nc),
                       xlab = "", ylab = ""))
  
  
  #grey.man <- man
  Gm <- subset(DF, value == "green" & gender == "male")
  with(Gm, rasterImage(man.grey, 
                       row-.465, col-.465, row+.465, col+.465, 
                       xlim = c(0, nr), ylim = c(0, nc),
                       xlab = "", ylab = ""))
  
  
  
  #grey.female <- female
  #green.female[,,2] <- female[,,4] # fill in green dimension
  Gf <- subset(DF, value == "green" & gender == "female")
  with(Gf, rasterImage(female.grey, 
                       row-.465, col-.465, row+.465, col+.465, 
                       xlim = c(0, nr), ylim = c(0, nc),
                       xlab = "", ylab = ""))
  
  mtext(RR, side = 3)
  mtext("Sample drawn from \n those with two affected parents", side = 1)
}

mk_plot6<- function(Kx,Kxx,Ky, Kyy, Kxy, Kyx){
  df <- as.data.frame(rbind(Kx,Kxx,Ky, Kyy, Kxy, Kyx))
  df$names <- row.names(df)
  names <- df$names
  par(mar =c(5,5,5,2))
  p7 <-barplot(df$V1, horiz = F,beside = F,names.arg=names, cex.names = 1,las=1,
               col = viridis(6))
  abline(h=0)
  
}


mk_plot7<- function(RRxx, RRyy, RRxy, RRyx){
  df1 <- as.data.frame(rbind(RRxx, RRyy, RRxy, RRyx))
  df1$names <- row.names(df1)
  names <- df1$names
  par(mar =c(5,5,5,2))
  p8 <-barplot(df1$V1, horiz = F,beside = F,names.arg=names, cex.names = 1,las=1,
               col = viridis(4))
  abline(h=0)
}


mk_plot8<- function(h2x, h2y,rg, covg){
  df2 <- as.data.frame(rbind(h2x, h2y,rg, covg))
  df2$names <- row.names(df2)
  names <- df2$names
  par(mar =c(5,5,5,2))
  p8 <-barplot(df2$V1, horiz = F,beside = F,names.arg=names, cex.names = 1,las=1,
               col = viridis(4))
  abline(h=0)
}

mk_plot9 <- function(rg){
  samples = 100
  r = rg
  
  library('MASS')
  data = mvrnorm(n=samples, mu=c(0, 0), Sigma=matrix(c(1, rg, rg, 1), nrow=2), empirical=TRUE)
  x = data[, 1]  # standard normal (mu=0, sd=1)
  y = data[, 2]  # standard normal (mu=0, sd=1)
  par(mar = c(5.1, 4.1, 4.1, 2.1)) 
  plot(x,y, pch = 16, frame = FALSE, main = "rg between x and y", cex.main =2)
  abline(lm(y ~ x, data = mtcars), col = "red")
  return(p)  
}  

mk_plot10 <- function(K, xlab=NA, ylab=NA){
  
  Txp = -qnorm(K, 0,1)
  z = dnorm(Txp) 
  i = z / K[1]
  mean = -(Txp) ;sd=1
  lb=0; ub=4
  
  x <- seq(Txp-6,Txp+3,length=1000)*sd + mean
  hx <- dnorm(x,mean,sd)
  K <- round(K,3)
  K_prevalence <- c(as.expression(bquote(italic(K)[xy]~'='~.(K))))
  #h2=h2
  p2 <- plot(x, hx, type="n", xlab="", ylab="",
             main="Liability x given \n relative diagnosed with y", cex.main =2, axes=FALSE)
  i <- x >= lb & x <= ub
  l <- x <= lb & x <= ub
  lines(x, hx)
  polygon(c(lb,x[i],ub), c(0,hx[i],0), col="darkblue")
  polygon(c(-5,x[l],lb), c(0,hx[l],0), col="grey") 
  
  
  mtext(K_prevalence, side = 1)
  
  abline(h=2)
  abline(v=0, h=-2)
  
  return(p2)  
}  

mk_plot11 <- function(K, xlab=NA, ylab=NA){
  
  Txp = -qnorm(K, 0,1)
  z = dnorm(Txp) 
  i = z / K[1]
  mean = -(Txp) ;sd=1
  lb=0; ub=4
  
  x <- seq(Txp-6,Txp+3,length=1000)*sd + mean
  hx <- dnorm(x,mean,sd)
  K <- round(K,3)
  K_prevalence <- c(as.expression(bquote(italic(K)[yx]~'='~.(K))))
  #h2=h2
  p2 <- plot(x, hx, type="n", xlab="", ylab="",
             main="Liability y given \n relative diagnosed with x", cex.main =2, axes=FALSE)
  i <- x >= lb & x <= ub
  l <- x <= lb & x <= ub
  lines(x, hx)
  polygon(c(lb,x[i],ub), c(0,hx[i],0), col="orange")
  polygon(c(-5,x[l],lb), c(0,hx[l],0), col="grey") 
  
  
  mtext(K_prevalence, side = 1)
  
  abline(h=2)
  abline(v=0, h=-2)
  
  return(p2)  
}  



# plot_liability ----------------------------------------------------------



mk_plot12 <- function(h2x,Kx){
  Tx = -qnorm(Kx, 0,1)
  Ks <- c(Kx)  #Lifetime risk
  Heritability <- c(h2x)
  Threshold <-c(Tx)
  #layout(matrix(c(1,2),1, 2, byrow = TRUE))
  h2 <- c(as.expression(bquote(italic(h)^2~'='~.(Heritability))))
  K_pop <- c(as.expression(bquote(italic(K)~'='~.(Ks))))
  T0 = Threshold
  z = dnorm(T0) 
  #print(z)
  i1 = z / Ks # mean phenotypic liability of those with disease
  # print(i)
  mean = -(T0) ;sd=1
  lb=0; ub=4
  x <- seq(T0-6,T0+3,length=1000)*sd + mean
  hx <- dnorm(x,mean,sd)
  plot(x, hx, type="n", xlab="", ylab="",
       main="Liability Threshold \n Model", cex.main = 1.5, axes=FALSE)
  i <- x >= lb & x <= ub
  l <- x <= lb & x <= ub
  lines(x, hx)
  polygon(c(lb,x[i],ub), c(0,hx[i],0), col="deeppink")
  polygon(c(-5,x[l],lb), c(0,hx[l],0), col="lightgrey")
  i_true <- i[which(i==T)]
  max.true <- length(i_true)
  text(hx[max.true],(z/3),"K", cex = 2, adj=-0.2, col="black")
  mtext("T", side=1, at=hx[max.true],col = "black", cex =2)
  text(hx[max.true],z,"z" , col ="black", cex=2 )
}

mk_plot13 <- function(K){
  
  Txp = -qnorm(K, 0,1)
  z = dnorm(Txp) 
  i = z / K[1]
  mean = -(Txp) ;sd=1
  lb=0; ub=4
  
  x <- seq(Txp-6,Txp+3,length=1000)*sd + mean
  hx <- dnorm(x,mean,sd)
  K <- round(K,2)
  K_prevalence <- c(as.expression(bquote(italic(K)~'='~.(K))))
  h2=h2
  p2 <- plot(x, hx, type="n", xlab="", ylab="",
             main="Lifetime risk \n general population", cex.main =1.5, axes=FALSE)
  i <- x >= lb & x <= ub
  l <- x <= lb & x <= ub
  lines(x, hx)
  polygon(c(lb,x[i],ub), c(0,hx[i],0), col="darkslategray1")
  polygon(c(-5,x[l],lb), c(0,hx[l],0), col="grey") 
  
  
  mtext(K_prevalence, side = 1)
  
#  abline(h=2)
#  abline(v=0, h=-2)
  
  return(p2)  
}  

mk_plot14 <- function(K, xlab=NA, ylab=NA){
  
  Txp = -qnorm(K, 0,1)
  z = dnorm(Txp) 
  i = z / K[1]
  mean = -(Txp) ;sd=1
  lb=0; ub=4
  
  x <- seq(Txp-6,Txp+3,length=1000)*sd + mean
  hx <- dnorm(x,mean,sd)
  K <- round(K,3)
  #K_prevalence <- c(as.expression(bquote(italic(K)~'='~.(K))))
  K_prevalence <- c(as.expression(bquote(italic(K)[1]~'='~.(K))))
  
  p2 <- plot(x, hx, type="n", xlab="", ylab="",
             main="Lifetime risk individuals with \n one affected relative", cex.main =1.5, axes=FALSE)
  i <- x >= lb & x <= ub
  l <- x <= lb & x <= ub
  lines(x, hx)
  polygon(c(lb,x[i],ub), c(0,hx[i],0), col="darkolivegreen1")
  polygon(c(-5,x[l],lb), c(0,hx[l],0), col="grey") 
  
  
  mtext(K_prevalence, side = 1)
  
  # abline(h=2)
  # abline(v=0, h=-2)
  
  return(p2)  
}  





