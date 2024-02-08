
# MASCARA

This repo contains the code to reproduce the MASCARA manuscript.

Some figures are generated with additional scripts, the png files are inlcuded here.
The scripts used to create them are contained in order in /run_first.

Dependencies managed with renv.

To get started, clone this repo, go to the /R directory and open a fresh Rstudio session by opening MASCARA_Tutorial.Rmd (or any of the .R or .Rmd files).
Doing this will initialise renv for this project.

Once open run this snippet:
```
if(!"renv" %in% installed.packages()[,1]){
  install.packages("renv")
}

renv::restore(prompt = FALSE)
```

This should ensure all required packages are installed.

The gtf file Oryza_sativa.IRGSP-1.0.50.gtf is required for annotation of Haider data, this is not necessary for running Paper_markdown.Rmd but can be obtained from: https://rapdb.dna.affrc.go.jp/download/irgsp1.html.
