# AD-Genetics-Prediction
A large-scale phenotype-based AD disease gene prediction

## Introduction
Alzheimer’s disease (AD) is a severe neurodegenerative disorder and has become a global public health problem. Intensive research has been conducted for AD. But the pathophysiology of AD is still not elucidated. Disease comorbidity often associates diseases with overlapping patterns of genetic markers. This may inform a common etiology and suggest essential protein targets. US Food and Drug Administration (FDA) Adverse Event Reporting System (FAERS) collects large-scale post-marketing surveillance data that provide a unique opportunity to investigate disease co-occurrence pattern. We aim to construct a heterogeneous network that integrates disease comorbidity network from FAERS with protein-protein interaction to prioritize the AD risk genes using network-based ranking algorithm.

## Methods
![methods](./figures/methods.jpg =50x60)

## Documentation for modules
1. __*assoc_rules:*__ Provides classes for generate and processing association rules from FAERS using FP-Growth algorithm
2. __*network:*__ Provides classes for DCN and DCN0PPI construction
3. __*graph_algorithm:*__ Provides classes for graph algorithms used in this project, including random walk with restart (rwr),random graph generation and *do novo* predition of AD risk genes
4. __*util:*__ Provides utility classes used in this project
