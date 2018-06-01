# instructions to use the GoSimSem package in R
# marylaurenbenton | 2018


# you'll want to install the GoSemSim package
source("https://bioconductor.org/biocLite.R")
biocLite("GOSemSim")

# Load it into the workspace.
library(GoSemSim)


# read in or specify your data as a list of lists in R
# you'll want this to be GO term IDs to compare
# i.e. list(c(GO:###, GO:###, GO:###),
#           c(GO:###, GO:###, GO:###),
#           c(GO:###, GO:###, GO:###))
listOfLists <- list(list())  #todo YOUR DATA HERE


# save the identifiers/name for each list in a separate one.
# i.e. c('method1', 'method2', 'method3')
listOfNames <- c()  #todo YOUR DATA HERE


# create similarity matrix between GO terms
# @param df        list of lists of GO term IDs
# @param names     names of methods included; used to name final matrix
# @param ontology  GO ontology abbreviation. Defaults to "BP" (biological process).
# @param organism  organism for ontology. Defaults to "Human".
# @param measure   semantic similarity metrics. Defaults to "Wang".
# @param combine   method of combinining similarity. Defaults to "BMA".
makeGOSimMatrix <- function(df, names, ontology = "BP", organism = "Human", measure = "Wang", combine = "BMA") {
  simMatrix <- matrix(data=NA, nrow=length(df), ncol=length(df))
  for (i in 1:length(df)) {
    for (j in 1:i) {
      simMatrix[i,j] <- mgoSim(df[[i]], df[[j]], ont = ontology, organism = organism, measure = measure, combine = combine)
    }
  }
  rownames(simMatrix) <- names
  colnames(simMatrix) <- names
  return(simMatrix)
}


# mgoSim does the comparison between two lists of GO terms
# can add in the other parameters or just use the defaults specified by the function
simMatrix <- makeGOSimMatrix(listOfLists, listOfNames)


# now you can print, plot, or do whatever with your new matrix!! 
# if you want to do the same thing with gene lists instead, try the mgeneSim function
# further documentation found here: https://guangchuangyu.github.io/software/GOSemSim/
#                                   https://bioconductor.org/packages/release/bioc/html/GOSemSim.html
#                                   https://www.ncbi.nlm.nih.gov/pubmed/20179076

