## Overall plan: 
Determine hierarchy of the use case first + provide detail (input-output type, how many, etc)
add the related implicit use cases underneath the direct main cases

# USE CASE - direct 

* Determine biomanufacturing relevance for microbe-of-interest previously flown in space
* Determine biomanufacturing relevance for microbe-of-interest never flown in space

# USE CASE - implicit 

* Search NASA database for microbe
Function - Check if microbial data + diff expression analysis exist already in NASA database
Input - two set of strings (the Latin name Genus + species) of a microbe
Output - Boolean True/False (“Yes it exists in the NASA database” = flown to space or “Not in NASA database” = not yet flown to space)

* Fetch differential expression analysis of transcriptomics (for training the ML)
Function - Perform differential expression analysis using a package (DESEQ2)
Input - FASTA files containing transcripts
Output -  a .csv file of tabulated strings (gene IDs, gene functions) and floats (gene expression level and statistics e.g. p-value, log fold-change, etc)

* Search KEGG database for microbe (taxonomy)
Function - Check if microbial data exists in KEGG database using the “Search: GENOME” feature
Input - two set of strings (the Latin name Genus + species) of a microbe
Output - Boolean True/False (“Yes it exists, here is the microbe’s three-letter KEGG code” or “No it doesn’t not exist in KEGG”) 

* Comparison of taxonomy
Function - Comparison of microbes if NASA data were not available
Input - FASTA files containing transcripts
Output - a .csv file of tabulated strings (gene IDs, gene functions) and floats (gene expression level and statistics e.g. p-value, log fold-change, etc)	

* Assign input microbe a degree of relatedness to microbes present in NASA database
Function - Use Convolutional neural network (CNN) ML algorithm to predict relatedness
Input - float of point coordinates of microbe of interest’s gene expression level
Output - float percentage of predicted microbial relatedness

* Determine if transcriptional expression of microbe has been observed after replicating in space, under extreme conditions
Function - 
Input - 
Output -

* Generate visualization of predicted viability scores for input microbe in-space environment.
Function - Visualize predicted survivability scores as set number of sliding scale bars or spider/radar plots using matplotlib.
Input - floats of predicted scores generated based on degree of relatedness score(s)
Output - Set number of sliding scale bars or a spider/radar plot of survivability scores.

