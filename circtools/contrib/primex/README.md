
<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Build Status](https://travis-ci.org/dieterich-lab/primex.svg?branch=master)](https://travis-ci.org/dieterich-lab/primex) [![codecov](https://codecov.io/gh/dieterich-lab/primex/branch/master/graph/badge.svg)](https://codecov.io/gh/dieterich-lab/primex/)

primex
======

An R package to ease primer design using Primer3.

Installation
------------

You can install `primex` from github with:

``` r
# install.packages("devtools")
devtools::install_github("dieterich-lab/primex")
```

Quick start
-----------

``` r
library(primex)

exonSeqs <- c(
  exon1 = paste0(
    "CTCACCATGGATGATGATATCGCCGCGCTCGTCGTCGACAACGGCTCCGGCATGTGCAAG",
    "GCCGGCTTCGCGGGCGACGATGCCCCCCGGGCCGTCTTCCCCTCCATCGTGGC",
    "AGGCACCAG"
    ), 
  exon2 = paste0(
    "GGCGTGATGGTGGGCATGGGTCAGAAGGATTCCTATGTGGGCGACGAGGCCCAGAGCAAG",
    "AGAGGCATCCTCACCCTGAAGTACCCCATCGAGCACGGCATCGTCACCAACTGGGACGAC",
    "ATGGAGAAAATCTGGCACCACACCTTCTACAATGAGCTGCGTGTGGCTCCCGAGGAGCAC",
    "CCCGTGCTGCTGACCGAGGCCCCCCTGAACCCCAAGGCCAACCGCGAGAAGATGACCCAG"
  )
)

seqOpts <- seqSettings(seqId = "transcript2", seq = exonSeqs) 

p3Opts  <- p3Settings() %>%  
  primerTm(min = 58, optimal = 63, max = 67)

primers <- design(seqOpts, p3Opts)
```

First two candidate pairs:

|                               | first                | second               |
|-------------------------------|:---------------------|:---------------------|
| PRIMER\_PAIR\_PENALTY         | 4.488355             | 4.494866             |
| PRIMER\_LEFT\_PENALTY         | 4.462257             | 0.032609             |
| PRIMER\_RIGHT\_PENALTY        | 0.026098             | 4.462257             |
| PRIMER\_LEFT\_SEQUENCE        | ACCAGGGCGTGATGGTGG   | ATGGATGATGATATCGCCGC |
| PRIMER\_RIGHT\_SEQUENCE       | CATGTCGTCCCAGTTGGTGA | CCACCATCACGCCCTGGT   |
| PRIMER\_LEFT                  | 117,18               | 6,20                 |
| PRIMER\_RIGHT                 | 244,20               | 134,18               |
| PRIMER\_LEFT\_TM              | 65.462               | 63.033               |
| PRIMER\_RIGHT\_TM             | 63.026               | 65.462               |
| PRIMER\_LEFT\_GC\_PERCENT     | 66.667               | 50.000               |
| PRIMER\_RIGHT\_GC\_PERCENT    | 55.000               | 66.667               |
| PRIMER\_LEFT\_SELF\_ANY\_TH   | 6.11                 | 0.00                 |
| PRIMER\_RIGHT\_SELF\_ANY\_TH  | 0.00                 | 0.00                 |
| PRIMER\_LEFT\_SELF\_END\_TH   | 0.00                 | 0.00                 |
| PRIMER\_RIGHT\_SELF\_END\_TH  | 0.00                 | 0.00                 |
| PRIMER\_LEFT\_HAIRPIN\_TH     | 43.01                | 27.69                |
| PRIMER\_RIGHT\_HAIRPIN\_TH    | 32.56                | 43.06                |
| PRIMER\_LEFT\_END\_STABILITY  | 9.4000               | 12.9000              |
| PRIMER\_RIGHT\_END\_STABILITY | 7.9000               | 7.9000               |
| PRIMER\_PAIR\_COMPL\_ANY\_TH  | 0.00                 | 0.62                 |
| PRIMER\_PAIR\_COMPL\_END\_TH  | 5.31                 | 0.00                 |
| PRIMER\_PAIR\_PRODUCT\_SIZE   | 128                  | 129                  |
