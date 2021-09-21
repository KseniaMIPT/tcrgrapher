# tcrgrapher

## Installation

```R
install.packages('devtools')
library(devtools)
devtools::install_github("KseniaMIPT/tcrgrapher")
```

## Quick start

```R
library(tcrgrapher)
library(data.table)
sample <- fread('sample.txt')
df <- pipeline_OLGA(sample)
```
