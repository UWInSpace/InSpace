# Component Specifications

## `usr_seq`
- *Description:* Returns a sequence for either type of input
- *Inputs:* A user-defined amino acid sequence or protein accession number for NCBI protein database (int), an email address for NCBI queries (strings)
- *Output:* Amino acid sequence for user-defined protein (strings)



## `get_seq`
- *Description:* Retrieve the protein sequence (string of amino acids)
- *Inputs:* An email address for NCBI queries (strings), protein accession number for NCBI protein database (int)
- *Output:* Amino acid sequence for user-defined protein (strings)



## `get_prot_mass`
- *Description:* Tabulate approximate total mass of input amino acid string based on the sum of
    the average masses for every amino acid in the string
- *Inputs:* Amino acid sequence for user-defined protein (strings)
- *Output*: Protein mass value for user-defined protein (floats)



## `get_biopy_feat`
- *Description:* Taken-in dataframe with a column specified as 'SEQUENCE' and uses Biopython to generate additional features (aromaticity, isoelectric point) from each protein sequence to help better train a machine learning model to predict log2FC
- *Inputs:* Amino acid sequence for user-defined protein  (strings)
- *Output*: BioPython feature values (aromaticity, isoelectric point) for user-defined protein (int,floats)



## `count_aa_types`
- *Description:* Function takes the sequence, counts how many times amino acids occur per each of 4 groups (Polar/nonpolar, negative/positive) and returns values as percentages of the total amino acids in the sequence
- *Inputs:* Amino acid sequence for user-defined protein (strings)
- *Output*: Count AA type values (Polar/nonpolar, negative/positive) for user-defined protein (floats)


## `get_feat`
- *Description:* takes the sequence that the user is querying and returns a dataframe appended 
    with all of the features used in the predictive model
- *Inputs:*  sequence that the user is querying (strings)
- *Output*: a dataframe appended with all of the features used in the predictive model (df)


## `scale_input_featr`
- *Description:* X (our input features), and rescale based on min-max normalization
    it returns the normalized X
- *Inputs:* X input features pulled from the master csv (df)
- *Output*: normalized input features X (df)


## `bagging_regr`
- *Description:* Takes in the master csv and user input csv_file that should be inputted as strings. The test ratio, random state and n_estimator are set (from previous ML optimization) for the non-parameterized model to make predictions.
- *Inputs:* master csv and user input csv_file name (strings)
- *Output*: predicted output (log2FC) based on user input (float)


## `predict_log2fc`
- *Description:* Takes the user input of either a protein accession number or 
    an amino acid sequence and predicts how much the sequence will change after 
    being in space (near-zero gravity and increased radiation exposure).
- *Inputs:* a protein accession number (int) or an amino acid sequence (strings) and an email address (strings)
- *Output*: predicted output (log2FC) (float)


## `multi_pred_log2fc`
- *Description:* Takes the user input of either a protein accession number or 
    an amino acid sequence and predicts how much the sequence will change after 
    being in space (near-zero gravity and increased radiation exposure), function will iterate through and return a list of predicted log2fc's
- *Inputs:* a list of protein accession number (list of int) or an amino acid sequence (strings) and an email address (strings)
- *Output*: a list of predicted output (log2FC) (list of floats)
