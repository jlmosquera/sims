# `sims` R-package

`sims` is an R-package for computing semantic similarities measures between terms of any ontology, and it is also for comparing lists of objects based on the terms of an ontology where these objects were annotated. The package was developed to deal with arbitrary ontologies, but it is particularly focused on the _Gene Ontology_ (GO).

__Contents__

1. Licensing
2. Overview
3. Installation
4. Documentation
5. Prerequisites / Dependencies
6. NOTE
7. Authors

## 1. Licensing

sims - An R-package for Computing Semantic Similarities
v 1.0.1 - Copyright (c) 2014 Jose Luis Mosquera - jlmosquera@gmail.com

Use and distribution subject to the terms of the GPL-2 license. See LICENSE for details.

## 2. Overview

`sims` package is devoted to compute semantic similarities measures between terms of an arbitrary ontology. Fourteen measures from two different approaches (i.e. node-based approach and edge-based appoach) have been implemented. Specifically:

 * __Node-based approach__

 It consists of seven semantic similarity measures based on the Information Content concept proposed by Resnik, Lin, Schlicker _et al._, Jiang and Conrath, Mazandu and Mulder, Pirro and Seco, and Pirro and Euzenat.
   
 * __Edge-based approach__
 
 There are implemented:

    * two semantic similarity measures proposed by Resnik, and Rada _et al._,
    * one distance measure proposed by Rada, and
    * four pseudo-distances proposed by Joslyn et al.

`sims` package was developed in order to deal with arbitrary ontologies, but it is particularly focused on the GO. For instance, there are some functions that allow building the refinements matrix (i.e. the accessibility matrix in terms of the Graph Theory), the mapping matrix (i.e a matrix that maps from Entrez Gene IDs to GO IDs), or plotting the DAG structure in order to compare two different list of GO terms.

The package can manage Entrez Gene IDs and GO IDs from any organism R-package. Moreover, given two lists of Entrez Gene IDs, `sims` package allows comparing statistically the semantic similarity profiles associated with such lists of genes and then plotting these profiles, as well as generating the induced _GO Directed Acyclic Graphs_ (GO DAG's).


## 3. Installation

In R console:

- Install 'devtools' package:
```R
    > install.packages("devtools")
```    
- Load the package 
```R
    > library(devtools)
``` 
- Install 'sims' package:
```R
    > install_github('sims')
``` 

## 4. Documentation

There are two manuals for this version of the package:

| sims-vignette.pdf | sins Vignette    |
| sims-manual.pdf   | Reference Manual |

## 5. Prerequisites / Dependencies

Depends : AnnotationDbi, expm, GOstats, plyr
Imports : Matrix, knitr, igraph, methods, plotrix, Rgraphviz, vegan
Suggests: org.Hs.eg.db

## 6. NOTE

`sims` is an R-package developed as a part of the PhD thesis called *_Methods and Models for the Analysis of Biological Significance Based on High-Throughput Data_* at the [University of Barcelona](http://www.ub.edu/web/ub/en/index.html?).

## 7. Authors

* Author and Mantainer: Jose Luis Mosquera, PhD (jlmosquera@gmail.com), PhD - Department of Statistics. Universtiy of Barcelona
* Advisor: Alex Sánchez, PhD - Department of Statistics. University of Barcelona
