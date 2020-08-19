# Convex Hull Method (CHM) 
Using concepts from computational geometry, these CHM implementations allow for the calculation of multi-dimensional production envelopes in large metabolic models. 

## Exact Implementation
This is a python implementation tested on Python 3.6.9 using `Qsopt` (version 2.5.10) for exact linear optimization. 

## Floating-point arithmetic Implementation
This is a MATLAB implementation tested on MATLAB R2020a using `gurobi` (version 9.0.2) for double-precision linear optimization. 

This implementation finds less extreme points overall, due to rounding, but scales easily to 6 dimensions (reactions of interest) on genome-scale metabolic models.

## Biological Examples

To compare precision and run time, we use both implementations to calculate production envelopes of the genome-scale metabolic model of *E. coli* [iJO1366](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3261703/) and its core model [EColiCore2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5206746/). 

To calculate multi-dimensional production envelopes of the syntropic exchange reactions of a microbial communities, we apply the floating point arithmetic implementation to a [biogas-producing community model](https://biotechnologyforbiofuels.biomedcentral.com/articles/10.1186/s13068-016-0429-x).










