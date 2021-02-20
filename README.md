# InSpace

## Use cases (direct):
* Determine biomanufacturing relevance for microbe-of-interest previously flown in space
* Determine biomanufacturing relevance for microbe-of-interest never flown in space

## Use cases (implicit): 
* Search KEGG database for microbe
** Function - check if microbial data exists in KEGG
** Input - two set of strings (the Latin name Genus + species) of a microbe
** Output - three/four letter code from KEGG (search our CSV) 

* Assign input microbe a degree of relatedness to microbes present in NASA database
** Function - use Convolutional neural network (CNN) ML algorithm to predict relatedness
** Input - float of point coordinates of microbe of interest’s gene expression level
** Output - float percentage of predicted microbial relatedness

* Generate visualization of predicted viability scores for input microbe in-space environment
** Function - visualize predicted survivability scores as set number of sliding scale bars or spider/radar plots using matplotlib
** Input - floats of predicted scores generated based on degree of relatedness score(s)
** Output - set number of sliding scale bars or a spider/radar plot of survivability scores
