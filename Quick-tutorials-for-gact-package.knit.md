---
title: 'Brief Introduction to R package gact'
author: "Palle Duun Rohde, Izel Fourie Sørensen, & Peter Sørensen"
date: "2022-10-29"
bibliography: qg2021.bib
biblio-style: apalike
link-citations: yes
output:
  bookdown::pdf_document2:
    dev: png
    includes:
      in_header: preamble.tex
  html_document:
    includes:
      in_header: mathjax_header.html
---




# Introduction
The practical is based on the R package `gact` (Rohde et al.2022)). This package provides an infrastructure for working with large-scale genomic association data linked to different types of genomic features.

The most recent version of `gact` can be obtained from github:


```r
library(devtools)
devtools::install_github("psoerensen/gact")
```

## Load packages used {.unlisted .unnumbered} 


```r
library(gact)
library(qgg)

library(org.Hs.eg.db)
library(reactome.db)

library(corrplot)
```

## Download and install GDT database {.unlisted .unnumbered} 
The function `gact()` dowload and install the GDT database: 





