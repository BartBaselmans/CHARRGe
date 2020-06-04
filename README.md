# Risk in relatives, heritability, SNP-based heritability and genetic correlations in psychiatric disorders: a review

On this Github, there are two supplementary files that complement the content discussed in **Risk in relatives, heritability, SNP-based heritability and genetic correlations in psychiatric disorders: a review**, published in Biological Psychiatry.

# R-Markdown

The R-markdown R-code  can be used to replicate the figures used in the manuscript. Furthermore it provides R-code how you can calculate heritability, genetic correlatio and (cross-disorder) risk in relatives using the 'classic' **liability threshold model**. It discusses the basic principles of GREML and LD score regression and provides R-code to simulate genetic data to estimate heritability and genetic correlation.

The R-libraries needed to knit the document are:
 ```  
 #install.packages("rmarkdown") 
library(rmarkdown)   
library(data.table)
library(knitr)
 ```  
# Shiny
The R-libraries needed to run this app are:
 ```  
 #install.packages("rmarkdown") 
library(shiny)
library(ggplot2)
library(dplyr)
library(png)
library(grid)
library(viridis)  
 ```  


The Shiny-app can be used for teaching purposes and consists of two pages. \
**Risk for relatives**:
You can change the life time risk of disease (*K*), heritability, and the genetic relationship between relatives (IBD) and estimate the risk in relatives having one or two affected parents.

**heritability, genetic correlation and cross-disorder risk ratio**:
By changing the genetic parameters lifetime risk of disease X and Y, the lifetime risk of disease (x or y) given relatives with disease (x or y) as well as the lifetime risk of disease x given relatives with disease y, you can derive the heritability (x and y), co-heritability as well as the genetic correlation between the two disorders.

