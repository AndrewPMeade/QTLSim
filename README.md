# QTLSim

QTLSim models the evolution of a quantitative trait locus (QTL) in a population under stabilising or directional selection. The QTL comprises genetic and environmental effects. Population size, the strength of selection, heritability, the size of the input (target) phenotypic variance and the size of the per generation mutational variance can all be set by the user. Haploid and a limited diploid modes are available. 

QTLSim does not model the genome or multiple genetic loci directly. For a more comprehensive, general purpose and flexible simulation package, see SLiM (https://messerlab.org/slim/). For a detailed description of the simulations, please see the accompanying paper. 

See ./Documents/QTLSimManual.pdf for more information. 

# System requirements
QTLSim is written in C and requires the GSL – GNU Scientific Library (https://www.gnu.org/software/gsl/). It has been tested on Windows 10, OS X (Intel and ARM) and Linux (Rocks 7.0). The program should run on any system that has a C99 compiler. 

# Installation guide
The program is a standalone binary and does not require installation. 

# Demo and Instructions for use
An example is provided “Example1.txt” which takes roughly 30 seconds to run. The program is run from the command line and takes a file which specifies the simulation options as a command line parameter, see "QTLSimManual.pdf” for a detailed description of simulation options and output. 

Previous simulation results can be reproduced by setting identical options and the same random number seed.  
