# Risk in relatives, heritability, SNP-based heritability and genetic correlations in psychiatric disorders: a review

On this Github, there are two supplementary files that complement the content discussed in [INSERT LINK TO PAPER], published in Biological Psychiatry.

# R-Markdown

The R-markdown R-code to replicate the figures used in the manuscript. Furthermore it provides R-code how you can calculate heritability, genetic correlatio and risk in relatives using the 'classic' **liability threshold model**. It discusses the basic principles of GREML and LD score regression and provides R-code to simulate genetic data to estimate heritability and genetic correlation.

The R-libraries needed to knit the document are:
 ```  
 #install.packages("rmarkdown") 
library(rmarkdown)   
library(data.table)
library(knitr)
 ```  
# Shiny

The Shiny-app can be used for teaching purposes and consists of three pages
**Risk for relatives**
You can change the life time risk of disease (*K*), heritability, and the genetic relationship between relatives (IBD)

**Probability**
The probability to be diagnosed of a complex disorder is very non-linear and depends on heritability and lifetime risk of disease. You can vary those parameters to see their influences on probability

**heritability and genetic correlation**
By changing the genetic parameters lifetime risk of disease X and Y, the lifetime risk of disease (x or y) given relatives with disease (x or y) as well as the lifetime risk of disease x given relatives with disease y you can derive the heritability (x and y), co-heritability as well as the genetic correlation.

