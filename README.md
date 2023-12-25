# clusterHIRE

The aim of this project is designed for detecting cell-type-specific risk-CpG sites in epigenome-wide association studies with unknown cancer subtypes. Previous name for this project is called MPCS.

The project contains the following directories:
1. clusterHIRE_package: including codes for HIRE's and clusterHIRE's algorithms written in R and C. There are different versions for clusterHIRE with different suffices (_v0, _v1, ...). Each version has some changes in C code.
2. GSE42861 and GSE77716: real data provided by NCBI.
3. Luo_GSE42861: results from GSE42861 using HIRE's algorithm.
4. cHIREewas: a R file which implements the clusterHIRE's algorithm.
