#Source the functions needed to run the app
source("Function_shiny1.R")

library(shiny)
library(ggplot2)
library(dplyr)
library(png)
library(grid)
library(viridis)    
library(shinyhelper)
library(htmlTable)

# Define UI for application that draws a histogram
# ui_h2 -------------------------------------------------------------------
ui <- navbarPage("CHARRGe",
                 tabPanel("Heritability",
                          fluidRow(
                            column(3, offset=0,
                                   wellPanel(
                                     id = "tPanel",style = "overflow-y:scroll; max-height: 600px",
                                     h6("For theory behind the model, click on the green help icon:"),
                                     #h5("Trait parameters"),
                                     sliderInput("h2_K", 
                                                  HTML("Lifetime risk of disease (K):"),
                                                  min = 0, max = 0.35, value = 0.15, step = 0.01) %>%
                                    helper(icon = "question",
                                                  colour = "green",
                                                  type = "markdown",
                                                  content = "Heritability"),
                                     sliderInput("h2_K1", 
                                                  HTML("Lifetime risk of disease with diagnosed relative (K<sub>1</sub>):"),
                                                  min = 0, max = 0.35, value = 0.2, step = 0.01),
                                     radioButtons("h2_a",
                                                  label = HTML("Coefficient of relationship (a<sub>R</sub>)"),
                                                  choices = c("0.5 (first-degree)",
                                                              "0.25 (second-degree)",
                                                              "0.125 (third-degree)",
                                                              "0.06125 (fourth-degree)"),
                                                  selected = "0.5 (first-degree)"),
                                   )
                            ),
# output_h2 ---------------------------------------------------------------
fluidRow(
                             column(8, offset=0,
                                     textOutput("text1"),
                                     br(),
                                     textOutput("text2"),
                              ),
                             column(4, offset = 0,
                                     uiOutput('Table1')
                                     ),
                             column(2, offset = 0,
                                     br(),
                             actionButton("button1", "Show result")
                                     ),
                             column(1, offset = 0,
                                    br(),
                                    uiOutput("Table2")
                                    ),
                             br(),
                             br(),
                              column(5, offset=0,
                                     plotOutput("plot1")
                              )
                            ),
                          ),
                        ),
# ui risk in relatives ----------------------------------------------------
               tabPanel("Risk in Relatives",
                fluidRow(
                    column(3, offset=0,
                           wellPanel(
                              id = "tPanel",style = "overflow-y:scroll; max-height: 600px",
                              h6("For theory behind the model, click on the green help icon:"),
                              sliderInput("st_Kx", 
                                            HTML("Lifetime risk of disease (K<sub>x</sub>):"),
                                            min = 0, max = 0.2, value = 0.15, step = 0.01) %>%
                              helper(icon = "question",
                                            colour = "green",
                                            type = "markdown",
                                            content = "RiskinRelatives"),
                              sliderInput("st_h2x", 
                                            HTML("heritability (h<sup>2</sup><sub>x</sub>):"),
                                            min = 0, max = 1, value = 0.5, step = 0.01),
                              sliderInput("st_Ky", 
                                           HTML("Lifetime risk of disease (K<sub>y</sub>):"),
                                           min = 0, max = 0.2, value = 0.01, step = 0.01),
                              sliderInput("st_h2y", 
                                           HTML("heritability (h<sup>2</sup><sub>y</sub>):"),
                                           min = 0, max = 1, value = 0.7, step = 0.01),
                              radioButtons("st_a",
                                           label = HTML("Coefficient of relationship (a<sub>R</sub>)"),
                                           choices = c("0.5 (first-degree)",
                                                         "0.25 (second-degree)",
                                                         "0.125 (third-degree)",
                                                         "0.06125 (fourth-degree)"),
                                             selected = "0.5 (first-degree)"),
                              br(),
                           )
                        ),
# Output_risk_in_relatives ------------------------------------------------
fluidRow(
                            column(8, offset=0,
                                    textOutput("text3"),
                                    br(),
                                    textOutput("text4")
                            ),
                            br(),
                            br(),
                            column(3, offset=0,
                                    uiOutput('Table3')
                            ),
                            column(2, offset = 0,
                                    br(),
                                    actionButton("button10", "Show")
                            ),
                            column(2, offset = 0,
                            uiOutput("Table4")
                            ),
                            br(),
                            br(),
                            br(),
                            column(8, offset=0,
                                    plotOutput("plot2")
                            )
                          ),
                        ),
                      ),
# ui_CDRR -----------------------------------------------------------------
#navbarMenu("Cross Disorder",
              tabPanel("Cross Disorder Risk in Relatives",
                fluidRow(
                    column(3, 
                          wellPanel(
                            h6("For theory behind the model, click on the green help icon:"),
                            sliderInput("CDRR_h2x", 
                                        HTML("Heritability of disease x (h<sup>2</sup><sub>x</sub>):"),
                                        min = 0, max = 1, value = 0.7, step = 0.01) %>%
                            helper(icon = "question",
                                        colour = "green",
                                        type = "markdown",
                                        content = "CDRR"),
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
                                        min = -1, max = 1, value = 0.7, step = 0.01),
                            selectInput("CDRR_a",
                                        label = HTML("coefficient of relationship (a<sub>R</sub>)"),
                                        choices = c("0.5 (first-degree)",
                                                    "0.25 (second-degree)",
                                                    "0.125 (third-degree)",
                                                    "0.06125 (fourth-degree)"),
                                        selected = "0.5 (first-degree)"),
                                    )
                                  ),
# Output_CDRR -------------------------------------------------------------
fluidRow(
                            column(8, offset=0,
                                   textOutput("text5"),
                                   br(),
                                   textOutput("text6"),
                            ),   
                            br(),
                            br(),
                            column(3, offset=0,
                                    uiOutput('Table5')
                            ),
          column(2, 
          br(),
          actionButton("button11", "Show"),
        ),
           column(2, offset = 0,
           br(),
           uiOutput("Table6"),
           br(),
           br(),
         ),
          column(8, offset=0,
          plotOutput("plot3")
        )
      )
    )
),



# Citation ----------------------------------------------------------------


#navbarMenu("Citation",
           tabPanel("About CHARRGe",
                    HTML(" <br> If you use the CHARRGe calculator, please cite: <br> 
                    <b> Bart M.L. Baselmans, Loïc Yengo, Wouter van Rheenen, Naomi R. Wray,
Risk in Relatives, Heritability, SNP-Based Heritability, and Genetic Correlations in Psychiatric Disorders: A Review,
Biological Psychiatry,2020 <br>

<br>ISSN 0006-3223,<br>
<br>https://doi.org/10.1016/j.biopsych.2020.05.034.<br>
<br>(http://www.sciencedirect.com/science/article/pii/S0006322320316693)<br>

<br>Abstract: The genetic contribution to psychiatric disorders is observed through the increased rates of disorders in the relatives of those diagnosed with disorders. These increased rates are observed to be nonspecific; for example, children of those with schizophrenia have increased rates of schizophrenia but also a broad range of other psychiatric diagnoses. While many factors contribute to risk, epidemiological evidence suggests that the genetic contribution carries the highest risk burden. The patterns of inheritance are consistent with a polygenic architecture of many contributing risk loci. The genetic studies of the past decade have provided empirical evidence identifying thousands of DNA variants associated with psychiatric disorders. Here, we describe how these latest results are consistent with observations from epidemiology. We provide an R tool (CHARRGe) to calculate genetic parameters from epidemiological parameters and vice versa. We discuss how the single nucleotide polymorphism–based estimates of heritability and genetic correlation relate to those estimated from family records.<br>

<br>Keywords: Family register data; Genetic correlation; GWAS; Heritability; Psychiatric genetics; Risk in relatives
 </b>")
                 
        )
    )
                

# Server ------------------------------------------------------------------
# Server_h2 ---------------------------------------------------------------
# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
#Text1
 output$text1 <- renderText({ 
     "The CHARRGe calculator can be used to get familiar with the Liability Threshold Model. It covers different
     aspects of the model and is divided in three sections (1) heritability (2) risk in relatives, and (3) cross disorder risk in relatives."
})
#Text2 
 output$text2 <- renderText({ 
   "On this page, you can set the lifetime risk of the general population (K), the life time risk of individuals with one affected 
   relative, and the coefficient of relationship (e.g. parent). With only these three parameter, it is possible to estimate the heritability of liability. 
   The table shows the correpsonding threshold for chosen K, and K1. Furthermore, it provides the parameter z, which is defined as the height 
   of the standard normal curve at threshold (T) and i, the mean phenotypic liability of cases drawn from the population, which is derived as: z/K. Now try to estimate the 
   heritability yourself! To check your answer, click on the 'Show result' button. Note, each time that you change the sliders, you have to click 
   on this button to see the 'updated' results. If the heritability goes higher than 1, then the ratio of K1/K selected is higher than is realistic. 
   For theory behind the liability threshold model, click on the green question mark. 
   The theory used assumes that only genetic factors contribute to shared risk between relatives."
 })

#Set condition so that K1 can not be lower than K 
 observeEvent(input$h2_K,
              {
              updateSliderInput(session, "h2_K1", min = input$h2_K)
              })   

#Obtaining parameters
  output$Table1 <- renderUI({ 
        K <- input$h2_K
        K1 <- input$h2_K1
        T0 <- est_h2(K=input$h2_K,Kr=input$h2_K1,a=as.numeric(gsub(" .*", "", input$CDRR_a)))$T 
        T1 <- est_h2(K=input$h2_K,Kr=input$h2_K1,a=as.numeric(gsub(" .*", "", input$CDRR_a)))$Tr
        z <- est_h2(K=input$h2_K,Kr=input$h2_K1,a=as.numeric(gsub(" .*", "", input$CDRR_a)))$z
        i <- est_h2(K=input$h2_K,Kr=input$h2_K1,a=as.numeric(gsub(" .*", "", input$CDRR_a)))$i
        #make vector of variables
        data <- round(c(K,K1,T0,T1,z,i),2)
        #table aspects
        rowNames <- c("Parameters")
        colNames <- c("K","K1","T","T1","z","i")
        #Actual making the table
          HTML(
                htmlTable(matrix(data,
                ncol=6, byrow = TRUE),
                header =  colNames,
                rnames = rowNames,
                n.rgroup = c(1),
                n.cgroup = c(1),
                caption=paste0("Parameters from the liability threshold model that are needed to calculate the heritability"))
              )
            })    
 
#Reactive table output (h2)
  df1 <- eventReactive(input$button1, 
              {
              data.frame(heritability=est_h2(K=input$h2_K,Kr=input$h2_K1,a=as.numeric(gsub(" .*", "", input$CDRR_a)))$h2)
              })
#Generate Table  
  output$Table2 <- renderTable({
    df1()
  })
 
#Figures  
  output$plot1 = renderPlot({
          layout(matrix(c(1,2,3), 1,3, byrow = TRUE))
                  mk_plot12(h2x=0.7,Kx=0.15) 
                  mk_plot13(K=input$h2_K)
                  mk_plot14(K=input$h2_K1)
          }, height = 300, width = 900)  
  
# Server_risk_relatives ---------------------------------------------------
#Text3
  output$text3 <- renderText({ 
    "This section calculates the expected risk in relatives, when lifetime risk and heritability are known, and assuming the liability threshold model of common polygenic disease. You can choose the 
    lifetime risk of the general population (N = 100 people icons) for two disorders (Kx and Ky) as well as their corresponding heritability estimates (h2x and h2y)
    using the sliders. By presenting two diseases side by side allows comparisons based on different combinations of parameters. The graphical output displays the liability model for (1) the general population, (2) the population drawn from those
    having one affected relative, and (3) population drawn from those having both parents affected." 
})      
#Text4
  output$text4 <- renderText({     
    "In the 'Heritability' section, it was shown that h2 could be estimated with the lifetime risk K and K1. If K1 is however, not known,
    it is possible to estimate it with only K, h2, and the coefficient of relationship. In the table, the threshold (Tx and Ty) correpsonding to lifetime risk Kx and Ky 
    are provided as well as their corresponding heritability esimates (h2x and h2y). With these parameters you should be able to estimate threshold T1x and T1y as well 
    as their corresponding lifetime risk parameters. To check your answer, click on the 'Show result' button. Note, each time that you change the sliders, you have to click 
   on this button to see the 'updated' results. When heritability is set to 1, the maximum risk in relatives expected is displayed. For theory behind the liability threshold model, click on the green question mark. 
   The theory used assumes that only genetic factors contribute to shared risk between relatives."
})

#Obtaining parameters for reactive table 3
  output$Table3 <- renderUI({ 
        h2x <- input$st_h2x
        h2y <- input$st_h2y
        Kx <- input$st_Kx
        Ky <- input$st_Ky
        Tx <- round(get_RR(h2=input$st_h2x, K=input$st_Kx, a=as.numeric(gsub(" .*", "", input$st_a)))$T,2)
        Ty <- round(get_RR(h2=input$st_h2y, K=input$st_Ky, a=as.numeric(gsub(" .*", "", input$st_a)))$T,2)
        ix <- round(get_RR(h2=input$st_h2x, K=input$st_Kx, a=as.numeric(gsub(" .*", "", input$st_a)))$i,2)
        iy <- round(get_RR(h2=input$st_h2y, K=input$st_Ky, a=as.numeric(gsub(" .*", "", input$st_a)))$i,2)
        #make matrix of obtained data
        data <- matrix(c(h2x,h2y,Kx, Ky, Tx, Ty,ix,iy),ncol = 4, byrow = F,
                          dimnames = list(c("trait x", "trait y"),
                          c("h2", "K", "T","i")))
        #Table aspects
        rowNames <- c("x", "y")
        colNames <- c("h2","K","T","i")
        #Actual making the table  
        htmlTable(data,
              ctable = c("solid", "double"),
              caption = "Parameters general population")
  }) 
  
#Reactive table output
      df2 <- eventReactive(input$button10, {
              data.frame("trait"=c("trait x","trait y"),
              "T1"=c(round(get_RR(h2=input$st_h2x, K=input$st_Kx, a=as.numeric(gsub(" .*", "", input$st_a)))$Tr,2),
                      round(get_RR(h2=input$st_h2y, K=input$st_Ky, a=as.numeric(gsub(" .*", "", input$st_a)))$Tr,2)),
              "K1"= c(round(get_RR(h2=input$st_h2x, K=input$st_Kx, a=as.numeric(gsub(" .*", "", input$st_a)))$Kr,2),
                      round(get_RR(h2=input$st_h2y, K=input$st_Ky, a=as.numeric(gsub(" .*", "", input$st_a)))$Kr,2))
                        )
                      })
#Call the table
      output$Table4 <-  renderTable({
                        df2()
                          })

#Figure2
#Get relevant parameters (obtained from get_RR function, which is sourced)    
  output$plot2 = renderPlot({
          marker1  = c(input$st_h2x, get_RR(h2=input$st_h2x, K=input$st_Kx, a=as.numeric(gsub(" .*", "", input$st_a)))$K)
          marker2  = c(input$st_h2y, get_RR(h2=input$st_h2y, K=input$st_Ky, a=as.numeric(gsub(" .*", "", input$st_a)))$K)
          marker3  = c(input$st_h2x, get_RR(h2=input$st_h2x, K=input$st_Kx, a=as.numeric(gsub(" .*", "", input$st_a)))$RRabs)
          marker4  = c(input$st_h2y, get_RR(h2=input$st_h2y, K=input$st_Ky, a=as.numeric(gsub(" .*", "", input$st_a)))$RRabs)
          marker5  = c(input$st_h2x, get_RR(h2=input$st_h2x, K=input$st_Kx, a=as.numeric(gsub(" .*", "", input$st_a)))$RRabs2)
          marker6  = c(input$st_h2y, get_RR(h2=input$st_h2y, K=input$st_Ky, a=as.numeric(gsub(" .*", "", input$st_a)))$RRabs2)
#Layout figure  
          layout(matrix(c(1,2,3,4,
                          5,6,7,8,
                          9,10,11,12), 3,4, byrow = TRUE))
#Fill layout figure (based on different functions)  
          mk_plot1(K=marker1[2], h2=marker1[1])
          mk_plot2(K=marker1[2])
          mk_plot1.1(K=marker2[2], h2=marker2[1])
          mk_plot2.2(K=marker2[2])
          mk_plot1.1par(K=marker3[2])
          mk_plot3(K=marker1[2],Kr=marker3[2])
          mk_plot1.1.1par(K=marker4[2])
          mk_plot3.1(K=marker2[2],Kr=marker4[2])
          mk_plot1.2par(K=marker5[2])
          mk_plot4(K=marker1[2],Kr2=marker5[2])
          mk_plot1.1.2par(K=marker6[2])
          mk_plot4.1(K=marker2[2],Kr2=marker6[2])
    }, height = 800, width = 900)

# Server_Cross_Disorder_risk_relatives ---------------------------------------------------      

  output$text5 <- renderText({ 
    "Just as epidemiological studies can investigate the increased risk of disorder x in relatives of those with disorder x, 
    they can also collect the data to estimate the increased risk of disorder y in relatives of those with disorder x, 
    and vice versa. Here, the top row (figure) is a visual representation of the input variables, such as the life time prevalence of 
    phenotype x and y as well as their heritability estimates. Additionally, the genetic correlation between x and y are provided. 
    The second row shows the visualization of the cross disorder risk ratio of disorder x when a relative is diagnosed with y,
    as well as the other way around."
    })
  
  output$text6 <- renderText({ 
    "Note, that a genetic correlation provides you information about the amount of overlap in genetic
    effects influencing two traits. However, if those traits have low heritability estimates, a high genetic correlation might not
    be that informative in terms of explaining the observed relation between the two traits. In addition to the genetic correlation,
    the co-heritability is a useful parameter. This parameter, is the proportion of the observed phenotypic correlation that is explained 
    by genetic effects and can be estimated using the heritability of two traits and the genetic correlation between them. 
    Try it yourself. Click on the green questionmark for more methodological background."
    })  

#Obtain relevant parameters (obtained from sliders)  
 output$Table5 <- renderUI({
        h2x <- input$CDRR_h2x
        h2y <- input$CDRR_h2y
        rg <- input$CDRR_rg
#Table properties
        rowNames <- c("heritability")
        colNames <- c("h2x","h2y","rg")
#making data matrix
        data <- matrix(c(h2x,h2y,rg),ncol = 3, byrow = F,
                        dimnames = list(c("estimates"),
                        c("h2x", "h2y", "rg")
                        )
                      )
#Calling the actual table                       
        htmlTable(data,
                ctable = c("solid", "double"),
                caption = "Chosen Parameters set by sliders")
        })

#Reactive Table output
  df3 <- eventReactive(input$button11, {
          data.frame(
            "Coheritability"=c(round(get_CDRR(h2x=input$CDRR_h2x,h2y=input$CDRR_h2y,Kx=input$CDRR_Kx,Ky=input$CDRR_Ky,rg=input$CDRR_rg,a=as.numeric(gsub(" .*", "", input$CDRR_a)))$covg,2))
               )
            })
#Calling the reactive Table
  output$Table6 <- renderTable({
  df3()
  })

##Figure3
#Get relevant parameters (obtained from get_RR function, which is sourced)    
output$plot3 = renderPlot({
        marker1  = get_CDRR(h2x=input$CDRR_h2x,h2y=input$CDRR_h2y,Kx=input$CDRR_Kx,Ky=input$CDRR_Ky,rg=input$CDRR_rg,a=as.numeric(gsub(" .*", "", input$CDRR_a)))$Kx
        marker2  = get_CDRR(h2x=input$CDRR_h2x,h2y=input$CDRR_h2y,Kx=input$CDRR_Kx,Ky=input$CDRR_Ky,rg=input$CDRR_rg,a=as.numeric(gsub(" .*", "", input$CDRR_a)))$h2x
        marker3  = get_CDRR(h2x=input$CDRR_h2x,h2y=input$CDRR_h2y,Kx=input$CDRR_Kx,Ky=input$CDRR_Ky,rg=input$CDRR_rg,a=as.numeric(gsub(" .*", "", input$CDRR_a)))$Ky
        marker4  = get_CDRR(h2x=input$CDRR_h2x,h2y=input$CDRR_h2y,Kx=input$CDRR_Kx,Ky=input$CDRR_Ky,rg=input$CDRR_rg,a=as.numeric(gsub(" .*", "", input$CDRR_a)))$h2y
        marker5  = get_CDRR(h2x=input$CDRR_h2x,h2y=input$CDRR_h2y,Kx=input$CDRR_Kx,Ky=input$CDRR_Ky,rg=input$CDRR_rg,a=as.numeric(gsub(" .*", "", input$CDRR_a)))$Kxy
        marker6  = get_CDRR(h2x=input$CDRR_h2x,h2y=input$CDRR_h2y,Kx=input$CDRR_Kx,Ky=input$CDRR_Ky,rg=input$CDRR_rg,a=as.numeric(gsub(" .*", "", input$CDRR_a)))$Kyx
        layout(matrix(c(1,2,3,4,
                        5,6,7,8
                        ), 2,4, byrow = TRUE))
        mk_plot1(K=marker1[1],h2=marker2[1])
        mk_plot9(rg=input$CDRR_rg)
        mk_plot1.1(K=marker3[1], h2=marker4[1])
        plot.new()
        mk_plot10(K=marker5[1])
        mk_plot3(K=marker1[1],Kr=marker5[1])
        mk_plot11(K=marker6[1])
        mk_plot3.1(K=marker3[1],Kr=marker6[1])
        }, height = 600, width = 900)


# Questionmarks ---------------------------------------------------       
#Activate the help questionmarks    
    observe_helpers(help_dir = getwd(),
                    withMathJax = TRUE)
    

}
# Run the application 
shinyApp(ui = ui, server = server)
