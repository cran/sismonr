## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  echo = TRUE,
  fig.align = "center",
  eval = F
)

## ---- eval = F-----------------------------------------------------------
#  install.packages("sismonr")

## ------------------------------------------------------------------------
#  library(sismonr)

## ------------------------------------------------------------------------
#  myinsilicosystem = createInSilicoSystem(G = 10, PC.p = 0.7, ploidy = 2)

## ---- eval = F-----------------------------------------------------------
#  ?insilicosystemargs

## ------------------------------------------------------------------------
#  class(myinsilicosystem)
#  names(myinsilicosystem)

## ------------------------------------------------------------------------
#  myinsilicosystem$genes

## ------------------------------------------------------------------------
#  ## system with only protein-coding genes, all regulators of transcription (PC.TC.p),
#  ## and all regulations are activations (positive regulation - TC.pos.p)
#  myinsilicosystem2 = createInSilicoSystem(G = 15, PC.p = 1, PC.TC.p = 1, TC.pos.p = 1)
#  myinsilicosystem2$genes
#  
#  ## Changing the function used to sample transcription rates for the genes
#  myinsilicosystem3 = createInSilicoSystem(G = 10,
#                                           basal_transcription_rate_samplingfct = function(x){runif(x, 0.1, 0.8)})
#  myinsilicosystem3$genes

## ------------------------------------------------------------------------
#  myinsilicosystem$edg

## ------------------------------------------------------------------------
#  myinsilicosystem$complexes
#  myinsilicosystem$complexeskinetics

## ------------------------------------------------------------------------
#  myinsilicosystem$complexesTargetReaction

## ---- fig.width = 6, fig.height = 6--------------------------------------
#  plotGRN(myinsilicosystem)

## ---- fig.width = 6, fig.height = 6--------------------------------------
#  plotGRN(myinsilicosystem, edgeType = "TC")

## ------------------------------------------------------------------------
#  names(myinsilicosystem$mosystem)
#  myinsilicosystem$mosystem$TCRN_edg

## ------------------------------------------------------------------------
#  myinsilicosystem$mosystem$TLRN_edg
#  myinsilicosystem$mosystem$RDRN_edg
#  myinsilicosystem$mosystem$PDRN_edg
#  myinsilicosystem$mosystem$PTMRN_edg

## ------------------------------------------------------------------------
#  mypop = createInSilicoPopulation(3, myinsilicosystem, ngenevariants = 4)

## ---- eval = F-----------------------------------------------------------
#  ?insilicoindividualargs

## ------------------------------------------------------------------------
#  class(mypop)
#  names(mypop)

## ------------------------------------------------------------------------
#  mypop$GenesVariants

## ------------------------------------------------------------------------
#  mypop2 = createInSilicoPopulation(3, myinsilicosystem, ngenevariants = 2)
#  mypop2$GenesVariants
#  
#  ## Creating a smaller system with only 3 genes
#  mysystem = createInSilicoSystem(G = 3, PC.p = 1)
#  
#  ## We will create only 1 variant of gene 1, 3 variants of gene 2 and
#  ## 2 variants of gene 3
#  nbvariants = c(1, 3, 2)
#  
#  qtlnames = c("qtlTCrate", "qtlRDrate",
#               "qtlTCregbind", "qtlRDregrate",
#               "qtlactivity", "qtlTLrate",
#               "qtlPDrate", "qtlTLregbind",
#               "qtlPDregrate", "qtlPTMregrate")
#  
#  genvariants = lapply(nbvariants, function(x){
#    matrix(1, nrow = length(qtlnames), ncol = x,
#           dimnames = list(qtlnames, 1:x))
#  })
#  names(genvariants) = mysystem$genes$id
#  
#  ## the 2nd variant of gene 2 has a mutation reducing its transcription rate by 3
#  genvariants$`2`["qtlTCrate", 2] = 0.33
#  ## and the 3rd variant has an increased translation rate
#  genvariants$`2`["qtlTLrate", 2] = 1.5
#  
#  ## The 2nd variant of gene 3 has a mutation decreasing the activity of
#  ## its active product
#  genvariants$`3`["qtlactivity", 2] = 0.7
#  
#  ## Allelic frequency of each variant
#  genvariants.freq = list('1' = c(1),
#                          '2' = c(0.6, 0.3, 0.1),
#                          '3' = c(0.9, 0.1))
#  
#  mypop3 = createInSilicoPopulation(10, mysystem,
#                                   genvariants = genvariants,
#                                   genvariants.freq = genvariants.freq)

## ------------------------------------------------------------------------
#  names(mypop$individualsList)

## ------------------------------------------------------------------------
#  mypop$individualsList$Ind1$haplotype
#  mypop$individualsList$Ind2$haplotype
#  mypop$individualsList$Ind3$haplotype

## ---- fig.width = 7, fig.height = 6--------------------------------------
#  plotMutations(mypop, myinsilicosystem, nGenesPerRow = 5)

## ---- fig.width = 6, fig.height = 6--------------------------------------
#  plotMutations(mypop, myinsilicosystem,
#                qtlEffectCoeffs = c("qtlTCrate", "qtlTLrate", "qtlRDrate", "qtlPDrate"),
#                inds = c("Ind1", "Ind2"),
#                alleles = "GCN2",
#                genes = 1:3)

## ------------------------------------------------------------------------
#  sim = simulateInSilicoSystem(myinsilicosystem, mypop, simtime = 1000, ntrials = 5)

## ------------------------------------------------------------------------
#  sim = simulateParallelInSilicoSystem(myinsilicosystem, mypop, simtime = 1000, ntrials = 5)

## ------------------------------------------------------------------------
#  sim$runningtime

## ------------------------------------------------------------------------
#  head(sim$Simulation)

## ---- fig.width = 7, fig.height = 6--------------------------------------
#  plotSimulation(sim$Simulation)

## ---- eval = F-----------------------------------------------------------
#  plotSimulation(sim$Simulation, yLogScale = F)

## ---- eval = F-----------------------------------------------------------
#  plotSimulation(sim$Simulation, inds = c("Ind1"), timeMin = 200, timeMax = 300)

## ---- fig.width = 6.5, fig.height = 6------------------------------------
#  plotHeatMap(sim$Simulation)

## ------------------------------------------------------------------------
#  simNoAllele = mergeAlleleAbundance(sim$Simulation)
#  head(simNoAllele)

## ------------------------------------------------------------------------
#  simNoPTM = mergePTMAbundance(simNoAllele)
#  head(simNoPTM)

## ------------------------------------------------------------------------
#  simNoComplex = mergeComplexesAbundance(simNoAllele)
#  head(simNoComplex)

