Meeting Note Feb 17,2021

Use cases should be
* 1(ish) main use case for ‘whole thing’
* Should be asked from perspective of user: “what can I use the overall tool for?”
* Use cases for each function that support solving the ‘whole thing’
* Components ≠ Usecase
* Should be asked from perspective of user: “what can I use parts of the code for?”
* Make sure to push Use Cases to GitHub page!

Taxonomy - this is of high importance
* 16S RNA Sequences
* K-mer (4-mer) counts
* K-mer counts would be measuring distance of particular k-mers, so a CNN might not be the best bet.  Better to use standard feedforward (see  also CNN and Feedforward Neural Networks, and Intro to Feedforward Neural Networks)
* MLST

Characterize bacteria using seven house-keeping genes
* Categorical measurement using Jaccard Distance
* See also An introduction to k-mers for genome comparison and analysis

Sequencing things
* Use DESeq2 to estimate variance-mean dependence in count data from high-throughput sequencing assays and test for differential expression based on a model using the negative binomial distribution.
* Might not be that important for our project 

16S-sequence
* The prediction model
* Probably don’t want to use CNN
* Probably want regression, or deep neural network
