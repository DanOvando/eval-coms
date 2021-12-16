# A History and Evaluation of Catch-Only Stock Assessment Models 

## Abstract

Understanding the status of fish stocks is a critical step in ensuring the ecological and economic sustainability of marine ecosystems. However, at least half of global catch and a vast majority of global fisheries lack formal stock assessments, largely due to a lack of sufficient data. Catch data, loosely referring to any catch records be it inclusive of discards or not, are the only type of fishery data available across a wide range of fisheries at a global scale. This has given rise to a long list of so-called “catch-only” models, intended to estimate aspects of stock status based primarily on characteristics of a fishery's catch history. In this paper, we review the history, performance and potential of “catch-only” models to estimate stock biomass status. While individual catch-only models often report good performance, repeated efforts to examine the performance of these models have consistently found them to be imprecise and biased when applied to new data-limited fisheries. We demonstrate that a large reason for this is the simple lack of information on stock status contained in the shape of a catch history alone. Off-the-shelf use of catch-only models can lead to poor and biased estimates of stock status, potentially hindering efforts at effective management.

## Reproducing Results


This project is set up with [`renv`](https://rstudio.github.io/renv/articles/renv.html) to manage package dependencies. Inside R (and with your working directory set correctly) run `renv::restore()`. Follow all prompts. This will install the correct versions of all the packages needed to replicate our results. Packages are installed in a stand-alone project library for this paper, and will not affect your installed R packages anywhere else. 

From there, `make-eval-coms.R` reproduces the results from Ovando _et al_ (2021)


 Ovando, D., Free, C.M., Jensen, O.P., Hilborn, R., 2021. A history and evaluation of catch-only stock assessment models. Fish and Fisheries n/a. https://doi.org/10.1111/faf.12637

  
