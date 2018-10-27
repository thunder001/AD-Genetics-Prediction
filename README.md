# AD-Genetics-Prediction
A large-scale phenotype-based AD disease gene prediction

## Introduction
Alzheimer’s disease (AD) is a severe neurodegenerative disorder and has become a global public health problem. Intensive research has been conducted for AD. But the pathophysiology of AD is still not elucidated. Disease comorbidity often associates diseases with overlapping patterns of genetic markers. This may inform a common etiology and suggest essential protein targets. US Food and Drug Administration (FDA) Adverse Event Reporting System (FAERS) collects large-scale post-marketing surveillance data that provide a unique opportunity to investigate disease co-occurrence pattern. We aim to construct a heterogeneous network that integrates disease comorbidity network from FAERS with protein-protein interaction to prioritize the AD risk genes using network-based ranking algorithm.

## Methods

<p align="center">
  <img src="./figures/methods.jpg" width="400">
</p>

## Documentation for modules
1. __*assoc_rules:*__ Provides classes for generate and processing association rules from FAERS using FP-Growth algorithm
2. __*network:*__ Provides classes for DCN and DCN0PPI construction
3. __*graph_algorithm:*__ Provides classes for graph algorithms used in this project, including random walk with restart (rwr),random graph generation and *do novo* predition of AD risk genes
4. __*util:*__ Provides utility classes used in this project

## Results
This folder contains three files listed below.

1. DCN_PPI_net.txt: The bipartite network file that contains 18,919 nodes (1,059 disease nodes and 17,860 gene nodes) and 712,803 edges.
Format: UMLS_ID or gene symbol|UMLS_ID or gene symbol|weight
Note: The edge weigh is set to 1.0 since the network is undirected and unweighted.

2. disUMLS_name.txt: A disease node mapping file from UMLS ID name to disease concept name
Format: UMLS_ID|disease_name

3. AD_novel_genes.csv: A file that contains novel AD risk genes we predicted from DCN_PPI network
Format: Rank|Gene



