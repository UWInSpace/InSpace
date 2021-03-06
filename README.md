[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
![versions](https://img.shields.io/pypi/pyversions/pybadges.svg)

# InSpace
Interstellar iN-Situ Predictive Algorithm for Cellular Engineering (InSpace), made by UW students, aims to predict relative viability of bacteria-of-interest (BOI) to spaceflight conditions (data sourced from [NASA GeneLab](https://genelab-data.ndc.nasa.gov/genelab/projects)).

__Current Functionality__:
-  Retrieve protein sequence from NCBI (GI number)
-  Obtain properties of a protein sequence 
	- MW, polarity, isoelectricity, aromaticity 
-  Predict microbial log2FC (viability metric) for space environment
	-  Non-linear regression (Bagging Regressor)

## Justification
As spaceflight technologies advance and mission lengths increase, there is an interest to develop technologies that decrease cost and reduce risks associated with providing biologically relevant resources to longer-term space missions [1]. A potential solution is in-situ resource utilization or ISRU, coupled with the process of in-situ manufacturing, where materials located on-site can be utilized as precursors to create other space-based commodities [2]. 

In regards to the design of these complex processes, it is crucial to understand the effects on gene expression and protein dynamics in response to two of the largest environmental stressors: space radiation and gravity [3]. Our tool aims to provide researchers with a  machine learning model to predict gene expression of important gene clusters based on previously flown model organisms (Escherichia coli, Pseudomonas aeruginosa, and Bacillus subtilis). 


[1] Menezes A. A.; Cumbers J.; Hogan J. A.; & Arkin A. P. (2015). Towards synthetic biological approaches to resource utilization on space missions. ​Journal of The Royal Society Interface​​12​(102), 20140715. https://doi.org/10.1098/rsif.2014.0715

[2] Mahoney, E. In-Situ Resource Utilization http://www.nasa.gov/isru (accessed Jan 23, 2021).

[3] Simpson, J. A. (1983). Elemental and isotopic composition of the galactic cosmic rays. Annual Review of Nuclear and Particle Science​, ​33​(1), 323-382.

## Dependencies
1. [BioPython](https://anaconda.org/bioconda/biopython)
	- Bio.Entrez, Bio.SeqUtils.ProtParam packages
	- Retrieving amino acid sequences and feature generation for predictive model
2. [Numpy](https://anaconda.org/anaconda/numpy)
	- Performing mathematical operations
3. [Pandas](https://anaconda.org/anaconda/pandas)
	- Manipulating and handling dataframes
4. [Scikit-Learn](https://anaconda.org/anaconda/scikit-learn)
	- Building machine learning models (Bagging regressor, ensemble)
5. [Seaborn](https://anaconda.org/anaconda/seaborn)
	- Visualizing selected gene clusters in respect to predicted log2FC





![DIRECT_projectcuration](https://user-images.githubusercontent.com/66701908/111225945-af292100-859d-11eb-86b8-337ea766999d.jpg)

