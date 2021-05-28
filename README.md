# Convex Hull Method (CHM) 
Using concepts from computational geometry, these CHM implementations allow for the calculation of multi-dimensional production envelopes in large metabolic models. 

## Exact Implementation
This is a python implementation tested on Python 3.6.9 using `Qsopt` (version 2.5.10) for exact linear optimization. 

## Floating-point Arithmetic Implementation (Double-Precision) 

Here we provide a Python and a MATLAB implementation. 

The MATLAB implementation was tested on MATLAB R2020a using `gurobi` (version 9.0.2) for double-precision linear optimization.

The Python implementation was tested on Python 3.7.6 using `gurobi` (version 9.0.2) for double-precision linear optimization. 

The floating point arithmetic implementations find less extreme points overall, due to rounding, but scale easily to higher dimensions on genome-scale metabolic models (tested for up to 6 reactions on interest).

Preliminary run-time differences between the MATLAB and the Python implementation suggest that the MATLAB implementation runs slightly faster, however, both floating-point arithmetic implementations outperform the run-time of the exact implementation by orders of magnitude. 

The Python implementation contains an additional feature, compared to the MATLAB implementation, which stores the stoichiometric constraints of the gurobi model such that there is no need to initialize the gurobi model from scratch each time that a linear program is computed. 

## Biological Examples

To compare precision and run time, we use both implementations to calculate production envelopes of the genome-scale metabolic model of *E. coli* [iJO1366](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3261703/) and its core model [EColiCore2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5206746/). 

To calculate multi-dimensional production envelopes of the syntropic exchange reactions of a microbial communities, we apply the floating point arithmetic implementation to a [biogas-producing community model](https://biotechnologyforbiofuels.biomedcentral.com/articles/10.1186/s13068-016-0429-x).










