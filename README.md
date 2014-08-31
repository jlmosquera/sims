# `sims` R-package

It is an R package is devoted to deal with arbitrary ontologies, compute semantic similarities between terms of the ontology and compare lists of objects annotated in these ontologies. The package is particularly focused on the _Gene Ontology_ (GO).

There are implemented fourteen measures from different approaches. Specifically:

 * __Node-based approach__

      There are implemented seven semantic similarity measures based on the Information Content concept proposed by 
      Resnik, Lin, Schlicker _et al._, Jiang and Conrath, Mazandu and Mulder, Pirro and Seco, and Pirro and Euzenat.
   
 * __Edge-based approach__
 
    There are implemented semantic similarities, a distance and pseudo-distances. More specifically,

    * two semantic similarity measures proposed by Resnik, and Rada _et al._,
    * one distance measure proposed by Rada, and
    * four pseudo-distances proposed by Joslyn et al.

The package can manage any ontology, but it is particularly focused on the Gene Ontology. In this regard, there are some functions that allow comparing statistically semantic cimilarity profiles and plotting these profiles as well as the induced GO DAG's.

## Installation

In R console:

- Install 'devtools' package:

    > install.packages("devtools")
    
- Load the package 

    > library(devtools)

- Install 'sims' package:

    > install_github('sims')

## NOTE

`sims` is an R-package developed as a part of a Thesis.

## Authors

 - Author and Mantainer: Jose Luis Mosquera (PhD candidate) - Department of Statistics. Universtiy of Barcelona
 - Thesis Advisor: Alex Sánchez, PhD - Department of Statistics. University of Barcelona
