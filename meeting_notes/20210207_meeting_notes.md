## CURRENT QUESTIONS FOR DAVE/STEPHANIE:
* Julia: Lots of data, what do we do with it? Implicit use case: finding out if a microbe has been to space or not
* Rory: How can we navigate these large repositories of data?  Apache SPARK?  Hadoop?  Amazon EMR?
* Widi: diff expression analysis → how to learn it in a short amount of time?
* Ava: How can we limit our scope? Rn we are working with lots of data, might not be plausible with our current constraints for CPU and memory

## NEW SCOPE:
Answers question: To what extent will the gene expression of microbe X be modified under extreme conditions (near-zero gravity and increased radiation)?
Use the extent of altered gene expression for sequenced microbes from NASA GeneLab and relatedness of those microbes to microbe X using taxonomy to predict extent of microbe X’s altered gene expression.

## Summary of what we have: 
* Databases: 
	* NASA gene lab (transcriptomics data)
	* Taxonomy database
	* KEGG
* Input: Microbial latin name
* ML-Prediction: relatedness based on taxonomy →  gene cluster (weighing more on genes for survival-replication, stress, growth - from taxonomy)
* Output: likelihood of altered gene expression (predicted by ML)

