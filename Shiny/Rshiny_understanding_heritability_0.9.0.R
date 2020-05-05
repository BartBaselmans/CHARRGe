#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#

#

library(shiny)
library(ggplot2)
library(dplyr)
#library(png)
library(grid)

#For undirected graph
# library("igraph")
# library("tidyverse")
# library("corrr")
# library("ggraph")

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

est_rg = function(Kx, Ky, Kx_x, Ky_y, Kx_y, a){
    
    # Kx_y, risk of disease x for relatives of disease y
    
    # parameters for trait X and Y
    Tx  = -qnorm(Kx, 0, 1)
    Ty  = -qnorm(Ky, 0, 1)
    
    zx  = dnorm(Tx)
    zy  = dnorm(Ty)
    
    ix  = zx/Kx
    iy  = zy/Ky
    
    
    Tx_x = -qnorm(Kx_x, 0, 1)
    Ty_y = -qnorm(Ky_y, 0, 1)
    
    Tx_y = -qnorm(Kx_y, 0, 1)
    
    h2x = (Tx-Tx_x * sqrt(1-(1-Tx/ix)*(Tx^2-Tx_x^2))) / (a * (ix+(ix-Tx)*Tx_x^2))
    h2y = (Ty-Ty_y * sqrt(1-(1-Ty/iy)*(Ty^2-Ty_y^2))) / (a * (iy+(iy-Ty)*Ty_y^2))
    
    hxy = (Tx-Tx_y * sqrt(1-(1-Ty/iy)*(Tx^2-Tx_y^2))) / (a * (iy+(iy-Ty)*Tx_y^2))
    
    rgxy  = hxy / sqrt(h2x * h2y)

        return(list("h2x"=h2x, "h2y"=h2y, "hxy"=hxy,
                "rgxy"=rgxy))
    
}

mk_plot = function(df, ylim=NA, xlab="", ylab="", m=NA){
    
    if(! is.na(ylim)){
        df = df %>% filter(y >= ylim[1] & y <= ylim[2])
    }
    
    p1 = ggplot(d=df, aes(x=x, y=y, group=a, colour=a)) +
        geom_line(size=1.5) +
        labs(x = xlab) +
        labs(y = ylab) +
        scale_colour_hue(labels=c("fourth-degree", 
                                  "third-degree",
                                  "second-degree",
                                  "first-degree"),
                         guide=guide_legend(nrow=2)) +
        geom_point(x=m[1], y=m[2], size=4.5, stroke=1, col="black", shape=21) +
        geom_point(x=m[1], y=m[2], size=1.5, col="black") +
        theme_classic() +
        theme(panel.grid.major = element_line(colour = "grey", linetype=3, size=0.3),
              panel.grid.minor = element_line(colour = "grey", linetype=3, size=0.2),
              legend.position = "bottom",
              legend.title = element_blank())
    
    return(p1)
}


mk_plot1 <- function(K, xlab=NA, ylab=NA){
    
    Txp = -qnorm(K, 0,1)
    z = dnorm(Txp) 
    i = z / K[1]
    mean = -(Txp) ;sd=1
    lb=0; ub=4
    
    x <- seq(Txp-6,Txp+3,length=1000)*sd + mean
    hx <- dnorm(x,mean,sd)
    round(K,2)
    K_prevalence <- c(as.expression(bquote(italic(K)~'='~.(K))))
    
    p2 <- plot(x, hx, type="n", xlab="", ylab="",
               main="", axes=FALSE)
    i <- x >= lb & x <= ub
    l <- x <= lb & x <= ub
    lines(x, hx)
    polygon(c(lb,x[i],ub), c(0,hx[i],0), col="darkblue")
    polygon(c(-5,x[l],lb), c(0,hx[l],0), col="grey") 
    
    #mtext(K_prevalence, side = 3, line = -3.5,at = 4,font = 2, cex = 1, adj=1)
    mtext(K_prevalence, side = 1)
    
      
    #axis(1, at=seq(-5, 2, 1), pos=0, cex.axis =1)
    abline(h=2)
    abline(v=0, h=-2)
    
return(p2)  
}  


mk_plot2 <- function(K, xlab="NA", ylab="NA"){
#par(mar=c(0,2,2,2))
par(mar=c(5.1, 4.1, 4.1, 2.1))
Prevalence <- round(K*100,0)
result <- matrix(nrow = 1, ncol = 100)
result <- sample(as.matrix(c(rep("b",Prevalence),rep("a",100-Prevalence))))

#sum(result[1,]=="a")
# input parameters - nr * nc should equal length(x)
cols <- c("grey", "darkblue")
nr <- 10
nc <- 10

# create data.frame of positions and colors
m <- matrix(cols[factor(result)], nr, nc)
DF <- data.frame(row = c(row(m)), col = c(col(m)[, nc:1]), value = c(m), gender =  sample(c(rep("male",50),rep("female",50))),
                 stringsAsFactors = FALSE)


p3 <- plot(col ~ row, DF, col = DF$value, xlab= "",ylab = "", pch = 15, cex = 1.8, asp = 1,
     xlim = c(0, nr), ylim = c(0, nc),
     axes = FALSE)
mtext("Sample drawn from \n the population", side =1)

return(p3)
}


mk_plot3 <- function(K, Kr, xlab="NA", ylab="NA"){
    #par(mar=c(0,2,2,2))
    par(mar=c(5.1, 4.1, 4.1, 2.1))
    Prevalence <- round(Kr*100,0)
    result <- matrix(nrow = 1, ncol = 100)
    result <- sample(as.matrix(c(rep("b",Prevalence),rep("a",100-Prevalence))))
    
    #sum(result[1,]=="a")
    # input parameters - nr * nc should equal length(x)
    cols <- c("grey", "darkblue")
    nr <- 10
    nc <- 10
    
    # create data.frame of positions and colors
    m <- matrix(cols[factor(result)], nr, nc)
    DF <- data.frame(row = c(row(m)), col = c(col(m)[, nc:1]), value = c(m), gender =  sample(c(rep("male",50),rep("female",50))),
                     stringsAsFactors = FALSE)
    
    
    #t <- get_RR(h2, K, a)$RRrel
    rr <- round(Kr/K,3)
    RR <- c(as.expression(bquote(italic(RR)[1]~'='~.(rr))))
    
    p4 <- plot(col ~ row, DF, col = DF$value, xlab= "",ylab = "", pch = 15, cex = 1.8, asp = 1,
               xlim = c(0, nr), ylim = c(0, nc),
               axes = FALSE)
    mtext("Sample drawn from \n those with one affected parent", side = 1)
    mtext(RR, side = 3)
    return(p4)
}



mk_plot4 <- function(K,Kr2, xlab="NA", ylab="NA"){
    #par(mar=c(0,2,2,2))
    par(mar=c(5.1, 4.1, 4.1, 2.1))
    Prevalence <- round(Kr2*100,0)
    result <- matrix(nrow = 1, ncol = 100)
    result <- sample(as.matrix(c(rep("b",Prevalence),rep("a",100-Prevalence))))
    
    #sum(result[1,]=="a")
    # input parameters - nr * nc should equal length(x)
    cols <- c("grey", "darkblue")
    nr <- 10
    nc <- 10
    
    # create data.frame of positions and colors
    m <- matrix(cols[factor(result)], nr, nc)
    DF <- data.frame(row = c(row(m)), col = c(col(m)[, nc:1]), value = c(m), gender =  sample(c(rep("male",50),rep("female",50))),
                     stringsAsFactors = FALSE)
    
    rr2 <- round(Kr2/K,2)
    RR <- c(as.expression(bquote(italic(RR)[2][par]~'='~.(rr2))))
    
    p5 <- plot(col ~ row, DF, col = DF$value, xlab= "",ylab = "", pch = 15, cex = 1.8, asp = 1,
               xlim = c(0, nr), ylim = c(0, nc),
               axes = FALSE)
    mtext("Sample drawn from \n those with two affected parents", side = 1)
    mtext(RR, side = 3)
    return(p5)
}


mk_plot5 <- function(h2,K, Kr, Kr2, xlab="NA", ylab="NA"){

ProbDisease1<- function(K,h2,xrange){
    sapply(xrange,function(x) pnorm( (x-qnorm(1-K))/sqrt(1-h2) ))
  } 

ProbDisease2<- function(Kr,h2,xrange){
  sapply(xrange,function(x) pnorm( (x-qnorm(1-Kr))/sqrt(1-h2) ))
} 

ProbDisease3<- function(Kr2,h2,xrange){
  sapply(xrange,function(x) pnorm( (x-qnorm(1-Kr2))/sqrt(1-h2) ))
} 


l3<- c(as.expression((bquote(population))),
       as.expression((bquote(1~parent))),
       as.expression((bquote(2~parents))))
       
xrange <- seq(-3,+5,len=100)
pop<- ProbDisease1(K=K,h2=h2,xrange)
par1 <- ProbDisease2(K=Kr,h2=h2,xrange)
par2 <- ProbDisease3(K=Kr2,h2=h2,xrange)

#par(mar=c(5.1, 4.1, 4.1, 2.1))
p6 <- matplot(xrange,cbind(pop,par1,par2), mgp=c(3,1,0),
        frame.plot = FALSE, type="l",lty=1:3,col=c("black","red","green"),lwd=2, 
        xlab = "Genetic Liability",cex.axis=1.2, cex.lab =1.2,
        ylab = "Probability", xaxt="n", cex.lab =1.2)
        axis(1, at = seq(-3, 5, by = 1), cex.axis=1.2)
        legend("topleft", legend = l3 ,
               col=c("black","red","green"), lty = 1:3, lwd=3, cex = 1.2,bty = "n")
        return(p6)
}

 mk_plot6<- function(h2x, h2y, hxy, rgxy ){
  df <- as.data.frame(rbind(h2x,h2y,hxy,rgxy))
  df$names <- row.names(df)
  names <- df$names
  par(mar =c(5,5,5,2))
  p7 <-barplot(df$V1, horiz = F,beside = F,names.arg=names, cex.names = 1,las=1,
             ylim = c(-1,1),
             col = c("darkgrey","darkgreen",
                     "darkorange","darkblue"))
  abline(h=0)
  return(p7)

}
# 
# 
# 
# mk_plot7<- function(rg){
#   rg <- tibble(x = "h2x", y="h2y", r = round(rg,2))
#   graph_cors <- rg %>%
#     filter(abs(r) > .01) %>%
#     graph_from_data_frame(directed = FALSE)
#   
#   p8 <- ggraph(graph_cors) +
#     geom_edge_link(aes(edge_alpha = abs(r), label =r, edge_width = abs(r), color = r),show.legend = FALSE) +
#     guides(edge_alpha = "none", edge_width = "none") +
#     scale_edge_colour_gradientn(limits = c(-1, 1), colors = c("blue", "red")) +
#     geom_node_point(color = "grey95", size = 15) +
#     geom_node_text(aes(label = name), repel = F) +
#     theme_graph() +
#     labs(title = "")
#   return(p8)
# }  


# Define UI for application that draws a histogram
ui <- navbarPage("Understanding heritability",
    
    navbarMenu("Risk for relatives",
        tabPanel("Single-trait",
            fluidPage(
                fluidRow(
                    column(3, offset=0,
                           wellPanel(
                               id = "tPanel",style = "overflow-y:scroll; max-height: 600px",
                                h4("Trait parameters"),
                                sliderInput("st_K", 
                                            "Lifetime risk of disease (K):",
                                            min = 0, max = 0.2, value = 0.1, step = 0.001),
                    
                                sliderInput("st_h2", 
                                            HTML("heritability (h<sup>2</sup>):"),
                                            min = 0, max = 1, value = 0.5, step = 0.01),
                        
                                radioButtons("st_a",
                                             label = HTML("IBD between relatives (a<sub>R</sub>)"),
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
                         #h6(textOutput("st_RRabs"), align="center")
                         ),
                        column(4, offset=0,
                               plotOutput("plot2")
                               #h6(textOutput("st_RRabs"), align="center")
                        ),
                        column(4, offset=0,
                               plotOutput("plot3")
                               #h6(textOutput("st_RRabs"), align="center")
                        ),
                        column(4, offset=0,
                               plotOutput("plot4")
                               #h6(textOutput("st_RRabs"), align="center")
                        ),
                        column(4, offset=3,
                               plotOutput("plot5")
                               #h6(textOutput("st_RRabs"), align="center")
                        ),
                        column(4, offset=0,
                               plotOutput("plot6")
                               #h6(textOutput("st_RRabs"), align="center")
                    )
                )
            )
        )
    )
  ),
  
  navbarMenu("Probability",
             tabPanel("single trait",
                      fluidRow(
                        column(4, 
                               wellPanel(
                                 h4("Trait parameters"),
                                 sliderInput("pr_K", 
                                             HTML("Lifetime risk of disease (K):"),
                                             min = 0, max = 0.5, value = 0.1, step = 0.001),
                                 
                                 sliderInput("pr_h2", 
                                             HTML("heritability (h<sup>2</sup>):"),
                                             min = 0, max = 1, value = 0.5, step = 0.01),
                                 

                                 radioButtons("pr_a",
                                              label = HTML("IBD between relatives (a<sub>R</sub>)"),
                                              choices = c("0.5 (first-degree)",
                                                          "0.25 (second-degree)",
                                                          "0.125 (third-degree)",
                                                          "0.06125 (fourth-degree)"),
                                              selected = "0.5 (first-degree)"),

                                 br(),
                                 actionButton("goButton2", "calculate"),
                               )
                          ),
                        
                        fluidRow(
                          column(4, offset=0,
                                 plotOutput("plot7"),
                                 #h6(textOutput("st_RRabs"), align="center")
                          )
                      )
             )
  )
),
navbarMenu("h2 & rg",
tabPanel("calculations",
         fluidRow(
           column(5, 
                  wellPanel(
                    h4("Trait parameters"),
                    sliderInput("rg_Kx", 
                                HTML("Lifetime risk of disease (K<sub>x</sub>):"),
                                min = 0, max = 0.2, value = 0.1, step = 0.001),
                    
                    sliderInput("rg_Ky", 
                                HTML("Lifetime risk of disease (K<sub>y</sub>):"),
                                min = 0, max = 0.2, value = 0.05, step = 0.001),
                    
                    # sliderInput("rg_Kz",
                    #             HTML("Lifetime risk of disease (K<sub>y</sub>):"),
                    #             min = 0, max = 0.5, value = 0.1, step = 0.001),
                    
                    sliderInput("rg_Kxx",
                                HTML("Lifetime risk of disease x given relatives with disease x (K<sub>Rx,x</sub>):"),
                                min = 0, max = 0.5, value = 0.15, step = 0.001),
                    
                    sliderInput("rg_Kyy", 
                                HTML("Lifteime risk of disease y given relatives with disease y (K<sub>Ry,y</sub>):"),
                                min = 0, max = 0.5, value = 0.15, step = 0.001),
                    
                    # sliderInput("rg_Kzz",
                    #             HTML("Lifteime risk of disease y given relatives with disease z (K<sub>Rz,z</sub>):"),
                    #             min = 0, max = 0.5, value = 0.15, step = 0.001),
                    
                    sliderInput("rg_Kxy",
                                HTML("Lifetime risk of disease x given relatives with disease y (K<sub>Rx,y</sub>):"),
                                min = 0, max = 0.5, value = 0.05, step = 0.001),
                    
                    # sliderInput("rg_Kxz",
                    #             HTML("Lifetime risk of disease x given relatives with disease y (K<sub>Rx,z</sub>):"),
                    #             min = 0, max = 0.5, value = 0.125, step = 0.001),
                    # 
                    # sliderInput("rg_Kyz",
                    #             HTML("Lifetime risk of disease x given relatives with disease y (K<sub>Ry,z</sub>):"),
                    #             min = 0, max = 0.5, value = 0.125, step = 0.001),

                    selectInput("rg_a",
                                label = HTML("IBD between relatives (a<sub>R</sub>)"),
                                choices = c("0.5 (first-degree)",
                                            "0.25 (second-degree)",
                                            "0.125 (third-degree)",
                                            "0.06125 (fourth-degree)"),
                                selected = "0.5 (first-degree)"),
                    br(),
                    actionButton("goButton3", "calculate"),
                  
                    
                    HTML('<hr style="border-color: black;">'),
                    h4(HTML("heritability x (<i>h<sup>2</sup><sub>x</sub></i>):")),
                    textOutput("h2x"),
                    h4(HTML("heritability y (<i>h<sup>2</sup><sub>y</sub></i>):")),
                    textOutput("h2y"),
                    h4(HTML("coheritability (<i>h<sub>xy</sub></i>):")),
                    textOutput("cohestxy"),
                    h4(HTML("genetic correlation (<i>r<sub>gxy</sub></i>):")),
                    textOutput("rgestxy"),
                    HTML('<hr style="border-color: black;">')
                  )
           ),
           fluidRow(
             column(4, offset=0,
                    plotOutput("plot8")
                    #h6(textOutput("st_RRabs"), align="center")
             ),
               column(5, offset=0,
                      plotOutput("plot9")
                      #h6(textOutput("st_RRabs"), align="center")
               )
)
)
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
            param   = expand.grid(a,h2)
            colnames(param) = c("a","b")
            param$c = get_RR(h2=param$b, K=input$st_K, a=param$a)$K
            param$a = as.factor(param$a)
            marker  = c(input$st_h2, get_RR(h2=input$st_h2, K=input$st_K, a=as.numeric(gsub(" .*", "", input$st_a)))$K)
            mk_plot1(K=marker[2], xlab="", ylab="")
                    
        })        
    })


    output$plot2 = renderPlot({
        input$goButton1
        isolate({
            a       = c(0.06125, 0.125, 0.25, 0.5)
            h2      = seq(0,1,0.01)
            param   = expand.grid(a,h2)
            colnames(param) = c("a","x")
            param$y = get_RR(h2=param$x, K=input$st_K, a=param$a)$K
            param$a = as.factor(param$a)
            marker  = c(input$st_h2, get_RR(h2=input$st_h2, K=input$st_K, a=as.numeric(gsub(" .*", "", input$st_a)))$K)
            mk_plot2(K=marker[2], xlab="", ylab="")
        })
    })
    
    output$plot3 = renderPlot({
        input$goButton1
        isolate({
            a       = c(0.06125, 0.125, 0.25, 0.5)
            h2      = seq(0,1,0.01)
            param   = expand.grid(a,h2)
            colnames(param) = c("a","b")
            param$c = get_RR(h2=param$b, K=input$st_K, a=param$a)$K
            param$a = as.factor(param$a)
            marker  = c(input$st_h2, get_RR(h2=input$st_h2, K=input$st_K, a=as.numeric(gsub(" .*", "", input$st_a)))$RRabs)
            mk_plot1(K=marker[2], xlab="Threshold Model", ylab="BLA")
            
        })        
    })
    
    output$plot4 = renderPlot({
        input$goButton1
        isolate({
            a       = c(0.06125, 0.125, 0.25, 0.5)
            h2      = seq(0,1,0.01)
            param   = expand.grid(a,h2)
            colnames(param) = c("a","x")
            param$y = get_RR(h2=param$x, K=input$st_K, a=param$a)$K
            param$a = as.factor(param$a)
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
            param   = expand.grid(a,h2)
            colnames(param) = c("a","b")
            param$c = get_RR(h2=param$b, K=input$st_K, a=param$a)$K
            param$a = as.factor(param$a)
            marker  = c(input$st_h2, get_RR(h2=input$st_h2, K=input$st_K, a=as.numeric(gsub(" .*", "", input$st_a)))$RRabs2)
            mk_plot1(K=marker[2], xlab="Threshold Model", ylab="BLA")
            
        })        
    })
    
    output$plot6 = renderPlot({
        input$goButton1
        isolate({
            a       = c(0.06125, 0.125, 0.25, 0.5)
            h2      = seq(0,1,0.01)
            param   = expand.grid(a,h2)
            colnames(param) = c("a","x")
            param$y = get_RR(h2=param$x, K=input$st_K, a=param$a)$K
            param$a = as.factor(param$a)
            marker  = c(input$st_h2, get_RR(h2=input$st_h2, K=input$st_K, a=as.numeric(gsub(" .*", "", input$st_a)))$RRabs2)
            marker1  = c(input$st_h2, get_RR(h2=input$st_h2, K=input$st_K, a=as.numeric(gsub(" .*", "", input$st_a)))$K)
            mk_plot4(K=marker1[2], Kr2=marker[2], xlab="Threshold Model", ylab="BLA")
        })
    })
    
    output$plot7 = renderPlot({
      input$goButton2
      isolate({
        a       = c(0.06125, 0.125, 0.25, 0.5)
        h2      = seq(0,1,0.01)
        param   = expand.grid(a,h2)
        colnames(param) = c("a","x")
        param$y = get_RR(h2=param$x, K=input$pr_K, a=param$a)$K
        param$a = as.factor(param$a)
        marker  = c(input$st_h2, get_RR(h2=input$pr_h2, K=input$pr_K, a=as.numeric(gsub(" .*", "", input$pr_a)))$K)
        marker1  = c(input$st_h2, get_RR(h2=input$pr_h2, K=input$pr_K, a=as.numeric(gsub(" .*", "", input$pr_a)))$RRabs)
        marker2  = c(input$st_h2, get_RR(h2=input$pr_h2, K=input$pr_K, a=as.numeric(gsub(" .*", "", input$pr_a)))$RRabs2)
        marker3  = c(input$st_h2, get_RR(h2=input$pr_h2, K=input$pr_K, a=as.numeric(gsub(" .*", "", input$pr_a)))$h2)
        mk_plot5(h2=marker3[2],K=marker[2], Kr=marker1[2], Kr2=marker2[2], xlab="Threshold Model", ylab="BLA")
      })
    })

    output$h2x = renderText({
      input$goButton3
      isolate({
        a = as.numeric(gsub(" .*", "", input$rg_a))
        signif(est_rg(Kx=input$rg_Kx,
                      Ky=input$rg_Ky,
                      Kx_x=input$rg_Kxx,
                      Ky_y=input$rg_Kyy,
                      Kx_y=input$rg_Kxy,
                      a=a)$h2x, 3)
      })
    })

    output$plot8 = renderPlot({
      input$goButton3
      isolate({
        # a       = c(0.06125, 0.125, 0.25, 0.5)
        # h2      = seq(0,1,0.01)
        # param   = expand.grid(a,h2)
        # colnames(param) = c("a","x")
        # param$y = get_RR(h2=param$x, K=input$pr_K, a=param$a)$K
        # param$a = as.factor(param$a)
        marker1  = c(est_rg(Kx=input$rg_Kx, Ky=input$rg_Ky, Kx_x=input$rg_Kxx, Ky_y=input$rg_Kyy, Kx_y=input$rg_Kxy, a=as.numeric(gsub(" .*", "", input$rg_a)))$h2x)
        marker2  = c(est_rg(Kx=input$rg_Kx, Ky=input$rg_Ky, Kx_x=input$rg_Kxx, Ky_y=input$rg_Kyy, Kx_y=input$rg_Kxy, a=as.numeric(gsub(" .*", "", input$rg_a)))$h2y)
        marker3  = c(est_rg(Kx=input$rg_Kx, Ky=input$rg_Ky, Kx_x=input$rg_Kxx, Ky_y=input$rg_Kyy, Kx_y=input$rg_Kxy, a=as.numeric(gsub(" .*", "", input$rg_a)))$hxy)
        marker4  = c(est_rg(Kx=input$rg_Kx, Ky=input$rg_Ky, Kx_x=input$rg_Kxx, Ky_y=input$rg_Kyy, Kx_y=input$rg_Kxy, a=as.numeric(gsub(" .*", "", input$rg_a)))$rgxy)
        mk_plot6(h2x=marker1[1],h2y=marker2[1],hxy=marker3[1], rgxy=marker4[1])
      })
    })    
    
    
     
    output$h2y = renderText({
      input$goButton3
      isolate({
        a = as.numeric(gsub(" .*", "", input$rg_a))
        signif(est_rg(Kx=input$rg_Kx,
                      Ky=input$rg_Ky,
                      Kx_x=input$rg_Kxx,
                      Ky_y=input$rg_Kyy,
                      Kx_y=input$rg_Kxy,
                      a=a)$h2y, 3)
      })
    })
    
    
    output$cohestxy = renderText({
      input$goButton3
      isolate({
        a = as.numeric(gsub(" .*", "", input$rg_a))
        signif(est_rg(Kx=input$rg_Kx,
                      Ky=input$rg_Ky,
                      Kx_x=input$rg_Kxx,
                      Ky_y=input$rg_Kyy,
                      Kx_y=input$rg_Kxy,
                      a=a)$hxy, 3)
      })
    })
  
    output$rgestxy = renderText({
      input$goButton3
      isolate({
        a = as.numeric(gsub(" .*", "", input$rg_a))
        signif(est_rg(Kx=input$rg_Kx,
                      Ky=input$rg_Ky,
                      Kx_x=input$rg_Kxx,
                      Ky_y=input$rg_Kyy,
                      Kx_y=input$rg_Kxy,
                      a=a)$rgxy, 3)
      })
    })    

    observeEvent(input$rg_Kx,{
      updateSliderInput(session, "rg_Kxx", min = input$rg_Kx)
      })  
    
    observeEvent(input$rg_Ky,{
      updateSliderInput(session, "rg_Kyy", min = input$rg_Ky)
      
    })  
    observeEvent(input$rg_Ky,{
      updateSliderInput(session, "rg_Kxy", min = input$rg_Ky)
      
    })  
      
}

# Run the application 
shinyApp(ui = ui, server = server)
