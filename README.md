# tcrgrapher

R package for identifying condition associated T cell clones from a single 
repertoire.

## Installation

```R
install.packages('devtools')
library(devtools)
devtools::install_github("KseniaMIPT/tcrgrapher")
```

To use tcrgrapher OLGA is needed. For detailed information please visit
https://github.com/statbiophys/OLGA.

OLGA can be installed using pip or pip3.

```python
pip install olga
```

## Quick start

```R
library(tcrgrapher)
library(data.table)
sample <- fread('sample.txt')
df <- pipeline_OLGA(sample)
```

## Details

```pipeline_OLGA``` is the main function that takes table with CDR3 sequences as
an input. Table should have the following columns (names of the columns are not 
important but the following order is necessary)
* Read.count - Number of unique reads per cdr3 sequence
* freq - Clonotype frequency in the clonoset
* cdr3nt - CDR3 nucleotide sequence
* bestVGene - TRBV segment
* bestVGene - TRBD segment
* bestJGene - TRBJ segment
* VEnd - Position of the end of V segment in CDR3 sequence
* DStart - Position of the start of D segment in CDR3 sequence
* DEnd - Position of the end of D segment in CDR3 sequence
*JStart - Position of the start of J segment in CDR3 sequence

You can find default parameters and possible options bellow.

pipeline_OLGA(df, Q = 6.27, cores = 1, thres_counts = 1, N_neighbors_thres = 1,
p_adjust_method = "BH", chain = "mouseTRB")
* df - data.table
* Q - selection factor. 1/Q sequences pass selection in the thymus. The 
default value for mouses 6.27. If a human model is taken and Q is not changed 
manually Q = 27 is used
* cores - number of used cores, 1 by default
* thres_counts - Only sequences with number of counts above this threshold
are taken into account
* N_neighbors_thres - Only sequences with number of neighbors above
threshold are used to calculate generation probability
* p_adjust_method - One of the method from p.adjust from stats package
possible options: "bonferroni", "holm", "hochberg", "hommel", "BH" or "fdr",
"BY", "none". "BH" is a default method.
* chain - Statistical model selection. Possible options: "mouseTRB", "humanTRB",
"humanTRA".

Function returns the same table that was in input filtered by number
of counts and number of neighbors with additional columns. Additional columns
are the following
* D - Number of neighbors in clonoset. Neighbor is a similar sequence
with one mismatch
* VJ_n_total - Number of unique sequences with given VJ combination
* Pgen - Probability to be generated computed by OLGA
* Pgen_sum_corr - Sum of Pgen of all sequences similar to the given
with one mismatch
* Pgen_by_VJ - Conditional probability. Sum of probabilities to be
generated with given VJ combination
* p_val - p value under null hypothesis that number of sequence's
neighbors is not more than in the random model
* p_adjust - p value with multiple testing correction
