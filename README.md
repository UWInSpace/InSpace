# InSpace

## Use cases (direct):
* Predicting relative viability of bacteria-of-interest (BOI) to spaceflight conditions
* Predicting expression profile of bacteria-of-interest (BOI) gene clusters in spaceflight conditions 

## Use cases (implicit): 
* Find matching genes/gene homologs between BOI and Test/train datasets 
  * Function - Use “package A” to connect FASTA from dataset ENTREZ IDs to NCBI BLAST of BOI
  * Input - *strings* of BOI’s name (Genus species strain)
  * Output - *boolean* T/F, *floats* of percent similarity, and *strings* of explanation to user 
   (e.g. T (40%): at least one (or more) match exists - calculating prediction, F (0%): “No match found, please specify another bacteria”) 

* Pull expression and significance data of matching genes from each NASA dataset 
  * User input: tagged ENTREZ IDs 
  * Function: match ENTREZ ID to desired expression level (log2FC) and significance (p values) 
  * Output: compiled data frame of expression level/significance data for matching genes 

* Determine if transcriptional data already exists for a bacteria-of-interest in the NASA GeneLab data repository.
  * Function - query the database for user-input bacteria along with transcriptional data
  * Input - two set of strings (the Latin name Genus + species) of bacteria
  * Output - message in the form of a sting indicating the transcriptional data does or does not already exist in the database, along with a visualization of viability scores if the bacterial data does exist

* Rank genes/gene homologs based on weighted user-specified conditions 
  * User input: compiled DF of expression/significance data, matching “strings” of user-specified conditions (gram +/-, microgravity, temperature, etc)
  * Function: determine weights of each gene/gene homologs for each dataset in respect to pulled database conditions (organized into new dataframe), ML clustering (python package) to both determine relevance of key pathways and visually represent overall changes
  * Output: clusters of predicted intact genes (neural network?), floats of percent overall predicted viability of organism 

* Assign input microbe a degree of relatedness to microbes present in NASA database
  * Function - use Convolutional neural network (CNN) ML algorithm to predict relatedness
  * Input - float of point coordinates of microbe of interest’s gene expression level
  * Output - float percentage of predicted microbial relatedness

* Generate visualization of predicted viability scores for input microbe in-space environment
  * Function - visualize predicted survivability scores as set number of sliding scale bars or spider/radar plots using matplotlib
  * Input - floats of predicted scores generated based on degree of relatedness score(s)
  * Output - set number of sliding scale bars or a spider/radar plot of survivability scores  

![image](https://user-images.githubusercontent.com/41084770/108634409-4df2af80-742e-11eb-85bc-301c1cf7c210.png)
