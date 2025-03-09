## Script to render Rmarkdown
# Author: Tyler Pollard
# Date: 22 Aug 2024

library(rmarkdown)

render(input = "HurricaneAnalysis.Rmd", 
       output_file = "README.md")

