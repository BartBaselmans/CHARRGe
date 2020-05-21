# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.

library(shiny)
library(ggplot2)
library(dplyr)
library(png)
library(grid)
library(viridis)    


get_CDRR = function(h2x, h2y, Kx, Ky, rg, a){

    # define the liability threshold model parameters
    Tx   = -qnorm(Kx, 0, 1) # liability threshold for disease x
    Ty   = -qnorm(Ky, 0, 1) # liability threhsold for disease y
    zx   = dnorm(Tx) # height of the PDF at t for trait x
    zy   = dnorm(Ty) # height of the PDF at t for trait y
    ix   = zx / Kx # mean liability of affected population
    iy   = zy / Ky # mean liability of affected population
    covg = rg * sqrt(h2x * h2y)
    
    # index has disorder x risk of disorder y for relative
    Tx_y = (Tx - a * covg * iy) / sqrt(1 - a * a * covg * covg * iy * (iy-Ty))
    Kx_y = 1 - pnorm(Tx_y)
    Zx_y = Kx_y/Kx # relative risk of disease x given that parent has disease y

    return(list("CD_RRabs"=Kx_y, "CD_RRratio"=Zx_y))
}

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

    return(list("RRabs"=round(Kr,3), "RRrel"=round(RR,3), "K"=round(K,3), "RRabs2"=round(Kr2,3), "h2"=round(h2,3)))
}

est_h2 = function(K, Kr, a){
    
    T  = -qnorm(K, 0, 1)
    Tr = -qnorm(Kr, 0, 1)
    z  = dnorm(T)
    i  = z / K
    
    h2 = (T - Tr * sqrt(1-(1-T/i)*(T^2-Tr^2))) / (a * i+(i-T)*Tr^2)
    
    return(h2)
    
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

mk_plot1 <- function(K, xlab=NA, ylab=NA){
    
    Txp = -qnorm(K, 0,1)
    z = dnorm(Txp) 
    i = z / K[1]
    mean = -(Txp) ;sd=1
    lb=0; ub=4
    
    x <- seq(Txp-6,Txp+3,length=1000)*sd + mean
    hx <- dnorm(x,mean,sd)
    K <- round(K,2)
    K_prevalence <- c(as.expression(bquote(italic(K)~'='~.(K))))
    
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


mk_plot3 <- function(K, Kr, xlab="NA", ylab="NA"){
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
  RR <- c(as.expression(bquote(italic(RR)[1][par]~'='~.(rr))))
  
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
  mtext("Sample drawn from those related to an affected person \n with coefficient of relationship as selected", side = 1)
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

    
rr2 <- round(Kr2/K,1)
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


# 
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


# Define UI for application that draws a histogram
ui <- navbarPage("CHARRGe",
    
    navbarMenu("Risk in relatives",
        tabPanel("Single-trait",
            fluidPage(
                fluidRow(
                    column(3, offset=0,
                           wellPanel(
                               id = "tPanel",style = "overflow-y:scroll; max-height: 600px",
                                h4("Trait parameters"),
                                sliderInput("st_K", 
                                            "Lifetime risk of disease (K):",
                                            min = 0, max = 0.2, value = 0.1, step = 0.01),
                    
                                sliderInput("st_h2", 
                                            HTML("heritability (h<sup>2</sup>):"),
                                            min = 0, max = 1, value = 0.5, step = 0.01),
                        
                                radioButtons("st_a",
                                             label = HTML("Coefficient of relationship (a<sub>R</sub>)"),
                                             choices = c("0.5 (first-degree)",
                                                         "0.25 (second-degree)",
                                                         "0.125 (third-degree)",
                                                         "0.06125 (fourth-degree)"),
                                             selected = "0.5 (first-degree)"),
                                br(),
                                actionButton("goButton1", "calculate")
                           )
                           
                     ),

            fluidRow(
                        column(4, offset=0,
                                plotOutput("plot1")
                         
                         ),
                        column(4, offset=0,
                               plotOutput("plot2")
                               
                        ),
                        column(4, offset=0,
                               plotOutput("plot3")
                               
                        ),
                        column(4, offset=0,
                               plotOutput("plot4")
                               
                        ),
                        column(4, offset=3,
                               plotOutput("plot5")
                               
                        ),
                        column(4, offset=0,
                               plotOutput("plot6")
                              
                    )
                )
            )
        )
    )
  ),

navbarMenu("Cross Disorder",
tabPanel("CDRR",
         fluidRow(
           column(5, 
                  wellPanel(
                    h4("Trait parameters"),
                    sliderInput("CDRR_h2x", 
                                HTML("Heritability of disease x (h<sup>2</sup><sub>x</sub>):"),
                                min = 0, max = 1, value = 0.7, step = 0.01),
                    
                    sliderInput("CDRR_h2y", 
                                HTML("Heritability of disease y (h<sup>2</sup><sub>y</sub>):"),
                                min = 0, max = 1, value = 0.35, step = 0.01),
                    

                    sliderInput("CDRR_Kx",
                                HTML("Lifetime risk of disease x (K<sub>x</sub>):"),
                                min = 0, max = 0.2, value = 0.01, step = 0.01),
                    
                    sliderInput("CDRR_Ky", 
                                HTML("Lifteime risk of disease y (K<sub>y</sub>):"),
                                min = 0, max = 0.2, value = 0.15, step = 0.01),
                    

                    sliderInput("CDRR_rg",
                                HTML("Genetic correlation (r<sub>g</sub>):"),
                                min = -1, max = 1, value = 0.35, step = 0.01),
        

                    selectInput("CDRR_a",
                                label = HTML("coefficient of relationship (a<sub>R</sub>)"),
                                choices = c("0.5 (first-degree)",
                                            "0.25 (second-degree)",
                                            "0.125 (third-degree)",
                                            "0.06125 (fourth-degree)"),
                                selected = "0.5 (first-degree)"),
                    br(),
                    actionButton("goButton3", "calculate"),
               
                    
                    HTML('<hr style="border-color: black;">'),
                    
                    h4(HTML("Lifetime risk of disorder x in relatives (aR) of those with disorder x (<i>K<sub>x,x</sub></i>):")),
                    textOutput("Kxx"),
                    h4(HTML("Increased Risk in relatives (aR) of those with disorder x (<i>R<sub>x,x</sub></i>):")),
                    textOutput("RRxx"),
                    
                    HTML('<hr style="border-color: black;">'),
                    h4(HTML("Lifetime risk of disorder y in relatives (aR) of those with disorder y (<i>K<sub>y,y</sub></i>):")),
                    textOutput("Kyy"),
                    h4(HTML("Increased Risk in relatives (aR) of those with disorder y (<i>R<sub>y,y</sub></i>):")),
                    textOutput("RRyy"),
                    
                    HTML('<hr style="border-color: black;">'),
                    h4(HTML("lifetime risk of disorder x in relatives (aR) of those with disorder y  (<i>K<sub>x,y</sub></i>):")),
                    textOutput("Kxy"),
                    h4(HTML("Increased Risk in relatives (aR) of disorder x of those with disorder y (<i>R<sub>x,y</sub></i>):")),
                    textOutput("RRxy"),
                    
                    HTML('<hr style="border-color: black;">'),
                    h4(HTML("lifetime risk of disorder y in relatives (aR) of those with disorder x  (<i>K<sub>y,x</sub></i>):")),
                    textOutput("Kyx"),
                    h4(HTML("Increased Risk in relatives (aR) of disorder y of those with disorder x (<i>R<sub>y,x</sub></i>):")),
                    textOutput("RRyx"),
                    
                    HTML('<hr style="border-color: black;">')
                  )
           ),
            fluidRow(
              column(3, offset=0,
                     plotOutput("plot9")
              ),
              column(3, offset=0,
                     plotOutput("plot7")
                         ),
                column(3, offset=1,
                       plotOutput("plot8")
                    )
           )
        )
      )
  ),


navbarMenu("Citation",
           tabPanel("CHARGGe",
                    HTML(" <br> If you use the CHARRGe calculator, please cite: <br> <b> Baselmans et al. (2020). Risk in relatives, 
                    heritability, SNP-based heritability and genetic correlations in psychiatric disorders: a review. <i>Biological Psychiatry </b>")
                 
        )
    )
                    
)
                       
# Define server logic required to draw a histogram
server <- function(input, output, session) {

    output$plot1 = renderPlot({
        input$goButton1
        isolate({
            a       = c(0.06125, 0.125, 0.25, 0.5)
            h2      = seq(0,1,0.01)
            marker  = c(input$st_h2, get_RR(h2=input$st_h2, K=input$st_K, a=as.numeric(gsub(" .*", "", input$st_a)))$K)
            mk_plot1(K=marker[2], xlab="", ylab="")
                    
        })        
    })


    output$plot2 = renderPlot({
        input$goButton1
        isolate({
            a       = c(0.06125, 0.125, 0.25, 0.5)
            h2      = seq(0,1,0.01)
            marker  = c(input$st_h2, get_RR(h2=input$st_h2, K=input$st_K, a=as.numeric(gsub(" .*", "", input$st_a)))$K)
            mk_plot2(K=marker[2], xlab="", ylab="")
        })
    })
    
    output$plot3 = renderPlot({
        input$goButton1
        isolate({
            a       = c(0.06125, 0.125, 0.25, 0.5)
            h2      = seq(0,1,0.01)
            marker  = c(input$st_h2, get_RR(h2=input$st_h2, K=input$st_K, a=as.numeric(gsub(" .*", "", input$st_a)))$RRabs)
            mk_plot1(K=marker[2], xlab="Threshold Model", ylab="")
            
        })        
    })
    
    output$plot4 = renderPlot({
        input$goButton1
        isolate({
            a       = c(0.06125, 0.125, 0.25, 0.5)
            h2      = seq(0,1,0.01)
            marker  = c(input$st_h2, get_RR(h2=input$st_h2, K=input$st_K, a=as.numeric(gsub(" .*", "", input$st_a)))$RRabs)
            marker1  = c(input$st_h2, get_RR(h2=input$st_h2, K=input$st_K, a=as.numeric(gsub(" .*", "", input$st_a)))$K)
            mk_plot3(K=marker1[2], Kr=marker[2], xlab="", ylab="")
        })
    })
    
  
    
    output$plot5 = renderPlot({
        input$goButton1
        isolate({
            a       = c(0.06125, 0.125, 0.25, 0.5)
            h2      = seq(0,1,0.01)
            marker  = c(input$st_h2, get_RR(h2=input$st_h2, K=input$st_K, a=as.numeric(gsub(" .*", "", input$st_a)))$RRabs2)
            mk_plot1(K=marker[2], xlab="Threshold Model", ylab="")
            
        })        
    })
    
    output$plot6 = renderPlot({
        input$goButton1
        isolate({
            a       = c(0.06125, 0.125, 0.25, 0.5)
            h2      = seq(0,1,0.01)
            marker  = c(input$st_h2, get_RR(h2=input$st_h2, K=input$st_K, a=as.numeric(gsub(" .*", "", input$st_a)))$RRabs2)
            marker1  = c(input$st_h2, get_RR(h2=input$st_h2, K=input$st_K, a=as.numeric(gsub(" .*", "", input$st_a)))$K)
            mk_plot4(K=marker1[2], Kr2=marker[2], xlab="Threshold Model", ylab="")
        })
    })


    output$Kxx = renderText({
      input$goButton3
      isolate({
        a = as.numeric(gsub(" .*", "", input$CDRR_a))
        signif(get_CDRR(h2x=input$CDRR_h2x,
                      h2y=input$CDRR_h2y,
                      Kx=input$CDRR_Kx,
                      Ky=input$CDRR_Ky,
                      rg=input$CDRR_rg,
                      a=a)$Kxx, 3)
      })
    })


    output$RRxx = renderText({
      input$goButton3
      isolate({
        a = as.numeric(gsub(" .*", "", input$CDRR_a))
        signif(get_CDRR(h2x=input$CDRR_h2x,
                        h2y=input$CDRR_h2y,
                        Kx=input$CDRR_Kx,
                        Ky=input$CDRR_Ky,
                        rg=input$CDRR_rg,
                        a=a)$RRxx, 3)
      })
    })
    
    
    output$Kyy = renderText({
      input$goButton3
      isolate({
        a = as.numeric(gsub(" .*", "", input$CDRR_a))
        signif(get_CDRR(h2x=input$CDRR_h2x,
                        h2y=input$CDRR_h2y,
                        Kx=input$CDRR_Kx,
                        Ky=input$CDRR_Ky,
                        rg=input$CDRR_rg,
                        a=a)$Kyy, 3)
      })
    })
    
    
    output$RRyy = renderText({
      input$goButton3
      isolate({
        a = as.numeric(gsub(" .*", "", input$CDRR_a))
        signif(get_CDRR(h2x=input$CDRR_h2x,
                        h2y=input$CDRR_h2y,
                        Kx=input$CDRR_Kx,
                        Ky=input$CDRR_Ky,
                        rg=input$CDRR_rg,
                        a=a)$RRyy, 3)
      })
    })
    
    output$Kxy = renderText({
      input$goButton3
      isolate({
        a = as.numeric(gsub(" .*", "", input$CDRR_a))
        signif(get_CDRR(h2x=input$CDRR_h2x,
                        h2y=input$CDRR_h2y,
                        Kx=input$CDRR_Kx,
                        Ky=input$CDRR_Ky,
                        rg=input$CDRR_rg,
                        a=a)$Kxy, 3)
      })
    })    
    
    output$RRxy = renderText({
      input$goButton3
      isolate({
        a = as.numeric(gsub(" .*", "", input$CDRR_a))
        signif(get_CDRR(h2x=input$CDRR_h2x,
                        h2y=input$CDRR_h2y,
                        Kx=input$CDRR_Kx,
                        Ky=input$CDRR_Ky,
                        rg=input$CDRR_rg,
                        a=a)$RRxy, 3)
      })
    })  
           
    output$Kyx = renderText({
      input$goButton3
      isolate({
        a = as.numeric(gsub(" .*", "", input$CDRR_a))
        signif(get_CDRR(h2x=input$CDRR_h2x,
                        h2y=input$CDRR_h2y,
                        Kx=input$CDRR_Kx,
                        Ky=input$CDRR_Ky,
                        rg=input$CDRR_rg,
                        a=a)$Kyx, 3)
      })
    })  
    
    
    output$RRyx = renderText({
      input$goButton3
      isolate({
        a = as.numeric(gsub(" .*", "", input$CDRR_a))
        signif(get_CDRR(h2x=input$CDRR_h2x,
                        h2y=input$CDRR_h2y,
                        Kx=input$CDRR_Kx,
                        Ky=input$CDRR_Ky,
                        rg=input$CDRR_rg,
                        a=a)$RRyx, 3)
      })
    })  
    

    output$plot7 = renderPlot({
      input$goButton3
      isolate({
        marker1  = c(get_CDRR(h2x=input$CDRR_h2x,h2y=input$CDRR_h2y,Kx=input$CDRR_Kx,Ky=input$CDRR_Ky,rg=input$CDRR_rg,a=as.numeric(gsub(" .*", "", input$CDRR_a)))$Kx)
        marker2  = c(get_CDRR(h2x=input$CDRR_h2x,h2y=input$CDRR_h2y,Kx=input$CDRR_Kx,Ky=input$CDRR_Ky,rg=input$CDRR_rg,a=as.numeric(gsub(" .*", "", input$CDRR_a)))$Kxx)
        marker3  = c(get_CDRR(h2x=input$CDRR_h2x,h2y=input$CDRR_h2y,Kx=input$CDRR_Kx,Ky=input$CDRR_Ky,rg=input$CDRR_rg,a=as.numeric(gsub(" .*", "", input$CDRR_a)))$Ky)
        marker4  = c(get_CDRR(h2x=input$CDRR_h2x,h2y=input$CDRR_h2y,Kx=input$CDRR_Kx,Ky=input$CDRR_Ky,rg=input$CDRR_rg,a=as.numeric(gsub(" .*", "", input$CDRR_a)))$Kyy)
        marker5  = c(get_CDRR(h2x=input$CDRR_h2x,h2y=input$CDRR_h2y,Kx=input$CDRR_Kx,Ky=input$CDRR_Ky,rg=input$CDRR_rg,a=as.numeric(gsub(" .*", "", input$CDRR_a)))$Kxy)
        marker6  = c(get_CDRR(h2x=input$CDRR_h2x,h2y=input$CDRR_h2y,Kx=input$CDRR_Kx,Ky=input$CDRR_Ky,rg=input$CDRR_rg,a=as.numeric(gsub(" .*", "", input$CDRR_a)))$Kyx)
        mk_plot6(Kx=marker1[1],Kxx=marker2[1],Ky=marker3[1],Kyy=marker4[1],Kxy=marker5[1], Kyx=marker6[1])
      })
    })    
  
    output$plot8 = renderPlot({
      input$goButton3
      isolate({
        marker1  = c(get_CDRR(h2x=input$CDRR_h2x,h2y=input$CDRR_h2y,Kx=input$CDRR_Kx,Ky=input$CDRR_Ky,rg=input$CDRR_rg,a=as.numeric(gsub(" .*", "", input$CDRR_a)))$RRxx)
        marker2  = c(get_CDRR(h2x=input$CDRR_h2x,h2y=input$CDRR_h2y,Kx=input$CDRR_Kx,Ky=input$CDRR_Ky,rg=input$CDRR_rg,a=as.numeric(gsub(" .*", "", input$CDRR_a)))$RRyy)
        marker3  = c(get_CDRR(h2x=input$CDRR_h2x,h2y=input$CDRR_h2y,Kx=input$CDRR_Kx,Ky=input$CDRR_Ky,rg=input$CDRR_rg,a=as.numeric(gsub(" .*", "", input$CDRR_a)))$RRxy)
        marker4  = c(get_CDRR(h2x=input$CDRR_h2x,h2y=input$CDRR_h2y,Kx=input$CDRR_Kx,Ky=input$CDRR_Ky,rg=input$CDRR_rg,a=as.numeric(gsub(" .*", "", input$CDRR_a)))$RRyx)
        mk_plot7(RRxx=marker1[1],RRyy=marker2[1],RRxy=marker3[1], RRyx=marker4[1])
      })
    })       

    output$plot9 = renderPlot({
      input$goButton3
      isolate({
        marker1  = c(get_CDRR(h2x=input$CDRR_h2x,h2y=input$CDRR_h2y,Kx=input$CDRR_Kx,Ky=input$CDRR_Ky,rg=input$CDRR_rg,a=as.numeric(gsub(" .*", "", input$CDRR_a)))$h2x)
        marker2  = c(get_CDRR(h2x=input$CDRR_h2x,h2y=input$CDRR_h2y,Kx=input$CDRR_Kx,Ky=input$CDRR_Ky,rg=input$CDRR_rg,a=as.numeric(gsub(" .*", "", input$CDRR_a)))$h2y)
        marker3  = c(get_CDRR(h2x=input$CDRR_h2x,h2y=input$CDRR_h2y,Kx=input$CDRR_Kx,Ky=input$CDRR_Ky,rg=input$CDRR_rg,a=as.numeric(gsub(" .*", "", input$CDRR_a)))$rg)
        marker4  = c(get_CDRR(h2x=input$CDRR_h2x,h2y=input$CDRR_h2y,Kx=input$CDRR_Kx,Ky=input$CDRR_Ky,rg=input$CDRR_rg,a=as.numeric(gsub(" .*", "", input$CDRR_a)))$covg)
        mk_plot8(h2x=marker1[1],h2y=marker2[1],rg=marker3[1],covg=marker4[1])
      })
    })        
    
}

# Run the application 
shinyApp(ui = ui, server = server)
