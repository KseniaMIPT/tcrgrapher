# tcrgrapher

R package for identifying condition associated T cell clones from a single 
repertoire.

## Installation

```R
install.packages('devtools')
library(devtools)
# to install the develop version
devtools::install_github("KseniaMIPT/tcrgrapher@develop")
```

To calculate generation probability, TCRgrapher can use OLGA or SONIA. 

For detailed information about OLGA, please visit
https://github.com/statbiophys/OLGA.

OLGA can be installed using pip or pip3.

```python
pip install olga
```

For detailed information about SONIA, please visit
https://github.com/statbiophys/SONIA

SONIA is a Python 2.7/3.6 software. It is available on PyPI and can be 
downloaded and installed through pip:

```python
pip install sonia
```

## Data loading

Data are stored in a TCRgrapher object. TCRgrapher is a S4 class that can be 
constructed by calling the TCRgrapher( ) function. TCRgrapher contains clonoset 
and metadata as data.tables.

There are three ways to initialize an object.

(1) Specify the path to the one file without metadata

If the path for only one file is specified, metadata will be produced automatically
and will contain one row and two columns: "file" and "sample_id".
Clonoset will have and additional column "sample_id" with one unique value and a
column "clone_id" with an unique id for every row.
Positions of "count", "cdr3nt", "cdr3aa", "V gene", "J gene" clonoset's columns must be specified.

```R
# See ?TCRgrapher
TCRgrObject <- TCRgrapher(file_path, 1, 3, 4, 5, 7) # positions of clonoset's columns
```

(2) Specify the path to the directory with files without metadata

If the path for the directory with files is specified and metadata is not specified,
metadata will be produced automatically and will have a number of rows equal
to the number of samples in the directory. All files will be merged into one
data.table with additional columns "sample_id"  and "clone_id" with
unique ids  for every row. Positions of "count", "cdr3nt", "cdr3aa",
"V gene", "J gene" clonoset's columns must be specified.

```R
TCRgrObject <- TCRgrapher(dir_path, 1, 3, 4, 5, 7)
```

(3) Specify the path to the directory with files and metadata

If the path for the directory with files and the path to the metadata are specified,
only files from the metadata column "file" will be taken into consideration.
Clonosets will be merged into one data.table with additional columns "sample_id" and
"clone_id" with unique ids for every row.
Values in the "sample_id" column will be the same as in the "sample_id"
metadata column. Positions of "count", "cdr3nt", "cdr3aa", "V gene", "J gene" clonoset's
columns and "file", "sampl_id" metadata's columns must be specified.

```R
TCRgrObject <- TCRgrapher(dir_path, 1, 3, 4, 5, 7, # positions of clonoset's columns
                         metadata_path, 1, 2)      # positions of metadtata's columns
```
```R
# get metadata
metadata(your_TCRgrapher_object)
# get clonoset
clonoset(your_TCRgrapher_object)
# add column
clonoset(your_TCRgrapher_object)$col <- your_vector
# subsetting (correct subsetting is important for edgeR analysis)
subset(your_TCRgrapher_object, vector_with_sample_ids_to_keep)
```

## ALICE pipeline

```ALICE_pipeline``` The function takes a TCRgrapher object as an input and performs
neighborhood enrichment analysis using the ALICE algorithm You can find default 
parameters and possible options below.

```R
tcrgrapher(TCRgrObject, Q_val = 6.27, cores = 1, thres_counts = 1, N_neighbors_thres = 1, 
          p_adjust_method = "BH", chain = 'mouseTRB', stats = 'OLGA', model= '-')
```
* TCRgrObject - TCRgrapher object that contains clonotype table
* Q - Selection factor. 1/Q sequences pass selection in a thymus. The 
default value for mice is 6.27. If a human model is taken and Q is not changed 
manually Q = 27 is used
* cores - number of used cores, 1 by default
* thres_counts - Only sequences with a number of counts above this threshold
are taken into account
* N_neighbors_thres - Only sequences with a number of neighbors above the 
threshold are used to calculate generation probability
* p_adjust_method - One of the methods from p.adjust from the stats package. 
Possible options: "bonferroni", "holm", "hochberg", "hommel", "BH" or "fdr",
"BY", "none". "BH" is a default method.
* chain - Statistical model selection. Possible options: "mouseTRB", "humanTRB",
"humanTRA".
* stats - Tool that will be used for generation probability calculation.
Possible options: "OLGA", "SONIA". "SONIA" also calculates Q for every sequence.
* model - Standard OLGA generation probability model is used by default. To set 
your one generation probability model, write "set_custom_model_VDJ 
<path_to_folder_with_model>". A generation probability model is usually IGOR output.
A folder should contain the following files: V_gene_CDR3_anchors.csv,
J_gene_CDR3_anchors.csv, model_marginals.txt, model_params.txt. Some models 
can be found in the folder "model".

Function returns TCRgrapher object with the same clonotype table that was in input
filtered by number of counts and number of neighbors with additional columns. 
Additional columns are the following
* D - Number of neighbors in clonoset. Neighbor is a similar sequence
with one mismatch
* VJ_n_total - Number of unique sequences with given VJ combination
* Pgen - Probability to be generated computed by OLGA
* Pgen_sum_corr - The sum of Pgen of all sequences similar to the given with one 
mismatch
* Pgen_by_VJ - Conditional probability. Sum of probabilities to be generated 
with given VJ combination
* p_val - p value under null hypothesis that number of sequence's
neighbors is not more than in the random model
* p_adjust - p value with multiple testing correction

## Additional functions for ALICE analysis

```pval_with_abundance(df)``` Function takes output of ALICE_pipeline function and 
adjusts p-value taking into account abundance of every clonotype. There are two 
additional columns in the output depending on count normalization: 
"pval_with_abundance_log2_counts" - Log2 is used for count normalization; 
"pval_with_abundance_counts" - no count normalization.

# edgeR analysis

```edgeR_pipeline``` function performs statistical analysis by edgeR to identify significantly
expanded clonotypes. First, it filters the data using relatively mild conditions
that better suit TCR repertoire data. Second, it normalize counts and estimate
dispersion by standard edgeR methods. Finally, it performs all pairwise comparisons
between groups and compares each group vs all others using glm and QL F-test.

To use edgeR_pipeline function you should create TCRgrapherCounts object. TCRgrapherCounts
is a subclass of TCRgrapher with additional slots: count_table and feature_info.
Count data could be aggregated in the following ways. First, clonotypes with the same
amino acid sequences and V gene would be merged together using default parameters 
```TCRgrapherCounts(TCRgrObject)```. Second, clonotypes with the same
amino acid sequences and VJ combination would be merged together using j_jene parameter
```TCRgrapherCounts(TCRgrObject, j_gene = TRUE)```. Third, clonotypes with the same
amino acid sequences regardless VJ would be merged. 
```TCRgrapherCounts(TCRgrObject, v_gene = FALSE, j_gene = FALSE)```.
Finally, clonotypes from the same connectivity component found by make_TCR_graph and
find_TCR_components functions would be groupped together by running 
```TCRgrapherCounts(TCRgrObject, cluster_id = TRUE)```

Metadata must contain a column that specifies levels of comparison. The name of the
column is the second parameter of the edgeR_pipeline function.

```R
library(edgeR)
# data loading
TCRgrObject <- TCRgrapher(file_path, 1, 3, 4, 5, 7, metadata_path, 1, 2)
# create a count table with aggregation by V segments
TCRgrCounts <- TCRgrapherCounts(TCRgrObject)
# run EdgeR pipeline 
edgeR_res <- edgeR_pipeline(TCRgrCounts, your_comparison)
```

## Additional functions

```make_TCR_graph(clonoset, v_gene = TRUE, j_gene = FALSE)``` 
Function makes graph from clonoset table with igraph package.
Every node of the graph is an unique clonotype from the table (one line). Edges
connects clonotypes with one amino acid mismatch or identical clonotypes if they
were in separate lines. If 'v_gene' and 'j_gene' are TRUE, edges connect only 
clonotypes with the same V gene or VJ combination.

```find_TCR_components(clonoset, graph)``` Function takes clonoset table and make_TCR_graph output
and returns the same clonoset table with additional column "cluster_id". All clusters of neighbors
with one mismatch have unique id. Function uses "components" function from igraph package.
