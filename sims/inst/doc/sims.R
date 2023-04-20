## ----options, echo=FALSE, eval=TRUE-----------------------
options(width = 60)

## ----loadPackage, eval=TRUE,  message=FALSE, warning=FALSE, results="hide"----
library("sims")
help("sims")

## ----listFunctions, eval=TRUE-----------------------------
ls("package:sims")

## ----loadJoslynDataset, eval=FALSE------------------------
#  data(joslyn)
#  help("joslyn")

## ----basicObjects, eval=TRUE------------------------------

## 1. Vocabulary of the ontology
vocabulary <- c("R", "B", "C", "K", "F", "G", "I", "E", "J", "H", 
"A", "D")

## 2. Links between terms (structure of the ontology)
origin <- c("B", "C", "K", "F", "G", "I", "I", "E", "J", "E", "J", 
"A", "A", "E", "H", "D", "D", "A")
terminus <- c("R", "R", "R", "B", "B", "B", "C", "C", "C", "K", 
"K", "F", "G", "I", "I", "E", "J", "H")
links <- data.frame(origin, terminus)
mat.g <- toMat(df = links, rnames = vocabulary, 
cnames = vocabulary)

print(mat.g)

## 3. Objects identifiers that are annotates in the ontology
object.ids <- letters[1:10]

## 4. Mapping from objects to terms (annotation of objects)
object <- c("b", "d", "f", "b", "g", "h", "i", "e", "a", "b", "c", 
"j")
term <- c("F", "F", "I", "E", "J", "J", "J", "H", "A", "A", "A", 
"D")
map <- data.frame(object, term)
mat.m <- toMat(df = map, rnames = object.ids, cnames = vocabulary)

print(mat.m)

## ----OOCclass, eval=FALSE---------------------------------
#  help("OOC")

## ----makeOOC, eval=TRUE-----------------------------------
joslyn.OOC <- toOOC(T = vocabulary, G = mat.g, O = object.ids, 
                    M = mat.m)
print(joslyn.OOC)

## ----helpsims, eval=FALSE---------------------------------
#  help("sims")

## ----ancestors, eval=TRUE---------------------------------
## Accessibility matrix
inv.IminusG <- inverseIminusG(joslyn.OOC)
A.mat <- getA(inv.IminusG)
print(A.mat)

## Ancestors
at <- ancestors(A.mat)
print(at)

## ----ICs, eval=TRUE---------------------------------------
resnik.sum <- resnikSummary(x = joslyn.OOC)
print(resnik.sum)

ic <- resnik.sum[, "ic"]

## ----sims, eval=TRUE--------------------------------------
## Computation of semantic similarity of Resnik
## (node-based approach)
sims.Res <- sims.nb(at, ic, method = "Res")
head(sims.Res)

## Computation of all semantic similarities of 
## edge-based approach
sims.all <- sims.nb(at, ic, method = "all")
head(sims.all)    

## ----summarySims, eval=TRUE-------------------------------
summarySims(sims.all)

## ----helpsimseb, eval=FALSE-------------------------------
#  help("sims.eb")

## ----sims.eb, eval=TRUE-----------------------------------
## Computation of semantic similarity of Resnik 
## (edge-based approach)
Resnik.eb <- sims.eb(x = joslyn.OOC, root = "R", at, 
                     method = "Rada")
head(Resnik.eb)

## Computation of all semantic similarities of 
## edge-based approach
sims.eb.all <- sims.eb(x = joslyn.OOC, root = "R", at, 
                       method = "all")
head(sims.eb.all)

## Summary of semantic similarities
summarySims(sims.eb.all)

## ----helpdistRada, eval=FALSE-----------------------------
#  help("distRada")

## ----summaryPaths, eval=TRUE------------------------------
sum.paths <- summaryPaths(x = joslyn.OOC, root = "R", len = TRUE)
head(sum.paths, 10)

## ----distRada, eval=TRUE----------------------------------
Rada <- distRada(sum.paths, at)
head(Rada)
summarySims(Rada)

## ----helppseudoDists, eval=FALSE--------------------------
#  help("pseudoDists")

## ----pseudoDists, eval=TRUE-------------------------------
## Computation of the pseudo-distance of 
## the minimum chain length
pd.hm <- pseudoDists(x = joslyn.OOC, root = "R", method = "hm")
head(pd.hm)

## Computation of all pseudo-distance
pd.all <- pseudoDists(x = joslyn.OOC, root = "R", method = "all")
head(pd.all)

## Summary of pseudo-distances
summarySims(pd.all)

## ----loadpratateIds, eval=TRUE----------------------------
data(prostateIds)
help("prostateIds")

## ----susbsetOfGenes, eval=TRUE----------------------------
## Entrez Gene ID's from Welsh et al. study
eg.we <- welsh01EntrezIDs[1:10]   

## Entrez Gene ID's from Singh et al. study
eg.sg <- singh01EntrezIDs[1:10]   

## ----organismPackage, eval=TRUE---------------------------
pckg <- "org.Hs.eg.db" 

## ----helpgosims, eval=FALSE-------------------------------
#  help("gosims")

## ----gosims, eval=TRUE------------------------------------
## All semantic similarities of node-based approach
all.nb <- gosims(eg = eg.we, ontology = "MF", pckg = pckg, 
                 type = "nb", method = "all")
summarySims(all.nb)                 

## ----gosimsAvsB, eval=TRUE--------------------------------
## Semantic similarity profiles computed with Resnik's measure
## from node-based approach
WEvsSG.nb <- gosimsAvsB(eg1 = eg.we, eg2 = eg.sg, ontology = "MF",
                        pckg = pckg, type = "nb", method = "Res")

## ----summarySimsAvsB, eval=TRUE---------------------------
summarySimsAvsB(WEvsSG.nb)

## ----plotHistSims, eval=TRUE, message=FALSE, warning=FALSE----
plotHistSims(x = WEvsSG.nb, freq = TRUE, 
             main = "Histogram of Semantic Similarities",
             xlab = "Semantic Similarity")

## ----sosimsProfiles, eval=TRUE----------------------------
gosimsProfiles(x = WEvsSG.nb, 
               col = c("tomato", "blue"), cex = 0.4, 
               top.labels = c("Welsh", "PAirs of GO IDs", "Singh"),
               main = "Semantic Similarity Profiles 
                       Between Welsh and Singh Studies",
               xlab = "Resnik")

## ----plotGODAG, eval=TRUE---------------------------------
## Induced subgraph associated with the list of genes from
## Welsh study
plotGODAG(eg1 = eg.we, eg2 = NULL, pckg = pckg, ontology = "MF")

## Induced subgraph associated with the list of genes from
## Singh study
plotGODAG(eg1 = eg.sg, eg2 = NULL, pckg = pckg, ontology = "MF")

## Induced subgraph associated with both lists of genes
plotGODAG(eg1 = eg.we, eg2 = eg.sg, pckg = pckg, ontology = "MF")

