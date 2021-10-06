## ----label= 'Load a sample dataset', eval=TRUE-----------------------------
library(RTN)
data(dt4rtn)

## ----label='Create a new TNI object', eval=TRUE, results='hide'------------
#Input 1: 'expData', a named gene expression matrix (samples on cols)
#Input 2: 'regulatoryElements', a named vector with TF ids
#Input 3: 'rowAnnotation', an optional data frame with gene annotation
tfs <- dt4rtn$tfs[c("PTTG1","E2F2","FOXM1","E2F3","RUNX2")]
rtni <- tni.constructor(expData=dt4rtn$gexp, regulatoryElements=tfs, rowAnnotation=dt4rtn$gexpIDs)

## ----label='Permutation', eval=TRUE, results='hide'------------------------
rtni <- tni.permutation(rtni, verbose = FALSE)

## ----label='Bootstrap', eval=TRUE, results='hide'--------------------------
rtni <- tni.bootstrap(rtni) 

## ----label='Run DPI filter', eval=TRUE, results='hide'---------------------
rtni <- tni.dpi.filter(rtni) 

## ----label='Check summary', eval=TRUE, results='hide'----------------------
tni.get(rtni, what="summary")
refnet <- tni.get(rtni, what="refnet")
tnet <- tni.get(rtni, what="tnet")

## ----label='Get graph', eval=TRUE------------------------------------------
g <- tni.graph(rtni)

## ----label='Create a new TNA object (preprocess TNI-to-TNA)', eval=TRUE, results='hide'----
#Input 1: 'object', a TNI object with a pre-processed transcripional network
#Input 2: 'phenotype', a named numeric vector, usually with log2 differential expression values
#Input 3: 'hits', a character vector of gene ids considered as hits
#Input 4: 'phenoIDs', an optional data frame with anottation used to aggregate genes in the phenotype
rtna <- tni2tna.preprocess(object=rtni, 
                         phenotype=dt4rtn$pheno, 
                         hits=dt4rtn$hits, 
                         phenoIDs=dt4rtn$phenoIDs
                         )

## ----label='Run MRA analysis pipeline', eval=TRUE, results='hide'----------
rtna <- tna.mra(rtna)

## ----label='Run overlap analysis pipeline', eval=TRUE, results='hide'------
rtna <- tna.overlap(rtna)

## ----label='Run GSEA analysis pipeline pt. 1', eval=TRUE, results='hide'----
rtna <- tna.gsea1(rtna, stepFilter=FALSE, nPermutations=100)
# ps. default 'nPermutations' is 1000.

## ----label='Run GSEA analysis pipeline pt. 2', eval=TRUE, results='hide'----
rtna <- tna.gsea2(rtna, tfs="PTTG1", nPermutations=100)
# ps. default 'nPermutations' is 1000.

## ----label='Get results', eval=TRUE, results='hide'------------------------
tna.get(rtna, what="summary")
tna.get(rtna, what="mra")
tna.get(rtna, what="overlap")
tna.get(rtna, what="gsea1")
tna.get(rtna, what="gsea2")

## ----label='Plot GSEA', eval=FALSE, results='hide'-------------------------
#  tna.plot.gsea1(rtna, file="tna_gsea1", labPheno="abs(log2) diff. expression")
#  tna.plot.gsea2(rtna, file="tna_gsea2", labPheno="log2 diff. expression")

## ----label='Session information', eval=TRUE, echo=FALSE--------------------
sessionInfo()

