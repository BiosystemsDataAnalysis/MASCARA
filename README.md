
# MASCARA

This repo contains the code to reproduce the MASCARA manuscript.

Some figures are generated with additional scripts, the png files are inlcuded here.
The scripts used to create them are contained in order in /run_first.

Dependencies managed with renv, open R and set working directory to ~/MASCARA/R, rum the following snippet to restore the environment:

`if(!"renv" %in% installed.packages()[,1]){
  install.packages("renv")
}

renv::restore(prompt = FALSE)`

This should enable all scripts to be run from their respective folders.

The gtf file Oryza_sativa.IRGSP-1.0.50.gtf is required for annotation of Haider data, this is not necessary for running Paper_markdown.Rmd. This can be obtained from: https://rapdb.dna.affrc.go.jp/download/irgsp1.html.
