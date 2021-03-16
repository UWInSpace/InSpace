## Use cases (direct):
* Predict the change in trascriptional rate that would be observed for a user input bacterial protein if the encoding gene was exposed to spaceflight conditions.
  * Function - Generate predicted Log2 FC (Flight-to-Ground) value for a user input protein, with the log2 FC representing the change in trascription rate of a gene in spaceflight conditions relative to ground conditions.
  * Input - Protein accession number as an *int* or the amino acid sequence of the protein as a *string*, along with the user's email as a *string*.
  * Output - Predicted Log2 FC of the encoding gene for the user input protein as a *numpy array*.

* Predict the change in trascriptional rate that would be observed for a series of user input bacterial proteins if the encoding genes were exposed to spaceflight conditions.
  * Function - Generate predicted Log2 FC (Flight-to-Ground) value for a series of user input proteins, with the log2 FC representing the change in trascription rate of a gene in spaceflight conditions relative to ground conditions.
  * Input - A series of protein accession number as a *list* of *int* values or a series of amino acid sequences as a *list* of *string* values, along with the user's email as a *string*.
  * Output - Predicted Log2 FC values in the form of a *list* of *float* values.

## Use cases (implicit): 
* Retrieve amino acid sequence for protein of interest. 
  * Function - Use the Bio.Entrez package to search NCBI Protein database and retrieve protein information, then isolate and return the amino acid FASTA sequence.
  * Input - Protein accession number as an *int*, and the user's email address in the form of a *string*.
  * Output - The amino acid sequence for the given protein in the form of a *string*.

* Characterize various features of a given protein based on it's amino acid sequence.
  * Function - Calculate the mass and the number of nonpolar, polar, negative, and positive amino acids in the sequence, and use Bio.ProteinAnalysis functionality to characterize the isoelectric point and aromaticity of protein.
  * Input -  Amino acid sequence of protein as a *string*.
  * Output - Table containing all of the calculated values for the seven chosen features as a *Pandas dataframe*. 

* Train a bagging regression model using Scikit-learn.
  * Function - Train a bagging regressor model on differential trascription data obtained from the NASA GeneLab database.
  * Input - A *.csv* file with differential trascriptomics data for a selection of bacterial protein encoding genes and the amino acid sequences of those genes. 
  * Output - A model instance with trained parameters.

* Generate reference visualizations from training data.
  * Function - Provide seven scatter plots, each showing log2FC for one of the seven defined features.  Points are represent protein encoding genes and are differentiated by species of origin.  Only points from a user selected gene cluster.
  * Input - A *Pandas dataframe* with differential trascriptomics data for a selection of bacterial protein encoding genes and one of eight possible user defined gene clusters as a *string*.
  * Output - Seven scatter plots with points corresponding to genes within selected gene cluster.