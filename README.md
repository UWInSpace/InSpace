[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
![versions](https://img.shields.io/pypi/pyversions/pybadges.svg)

# InSpace
Interstellar iN-Situ Predictive Algorithm for Cellular Engineering (InSpace), made by UW students, aims to predict relative viability of bacteria-of-interest (BOI) to spaceflight conditions (data sourced from NASA GeneLab).

__Current Functionality__:
- [ ] A
- [ ] B
- [ ] C

## Justification
As spaceflight technologies advance and mission lengths increase, there is an interest to develop technologies that decrease cost and reduce risks associated with providing biologically relevant resources to longer-term space missions [1]. A potential solution is in-situ resource utilization or ISRU, coupled with the process of in-situ manufacturing, where materials located on-site can be utilized as precursors to create other space-based commodities [2]. 

In regards to the design of these complex processes, it is crucial to understand the effects on gene expression and protein dynamics in response to two of the largest environmental stressors: space radiation and gravity [3]. Our tool aims to provide researchers with a  machine learning model to predict gene expression of important gene clusters based on previously flown model organisms (Escherichia coli, Pseudomonas aeruginosa, and Bacillus subtilis). 


[1] Menezes A. A.; Cumbers J.; Hogan J. A.; & Arkin A. P. (2015). Towards synthetic biological approaches to resource utilization on space missions. ​Journal of The Royal Society Interface​​12​(102), 20140715. https://doi.org/10.1098/rsif.2014.0715

[2] Mahoney, E. In-Situ Resource Utilization http://www.nasa.gov/isru (accessed Jan 23, 2021).

[3] Simpson, J. A. (1983). Elemental and isotopic composition of the galactic cosmic rays. Annual Review of Nuclear and Particle Science​, ​33​(1), 323-382.

## Dependencies
1. [BioPython](https://anaconda.org/bioconda/biopython)
	- Bio.Entrez package
2. [DAVID](https://david.ncifcrf.gov/content.jsp?file=release.html)
	- Bioinformatics Resources
4. [Numpy](https://anaconda.org/anaconda/numpy)
	- Mathematical operations
5. [Pandas](https://anaconda.org/anaconda/pandas)
	- Handling dataframes
6. [Scikit-Learn](https://anaconda.org/anaconda/scikit-learn)
	- Machine learning models (ensemble)






![image](https://www.flickr.com/photos/192507756@N04/51041049311/in/dateposted-public/)

