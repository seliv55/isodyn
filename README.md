# Isodyn
Version: 1.0

![Logo](text3923.png)

## Short Description

“C++”-program simulating the dynamics of metabolites and their isotopic isomers in central metabolic network using a corresponding kinetic model

## Description

“Isodyn” is a C++-program that supports an analysis of stable isotope tracer data to assess metabolic flux profiles in living cells. Isodyn simulates propagation of 13C label from artificially 13C enriched substrates into metabolites of central metabolic pathways forming various isotopic isomers (isotopomers). It fits the simulated dynamics of mass isotopomers to that observed experimentally, thus finding the parameters, which reflect the characteristics of corresponding biochemical reactions, and metabolic fluxes that correspond to the real fluxes for the analyzed biological object and conditions. Isodyn contains tools that check the goodness of fit and perform a statistical analysis of obtained metabolic fluxes. It using the change of metabolite concentrations in the analysis in addition to labeling data, as a restrictive factor. Optionally it uses manually corrected files provided by Midcor. It automatically adjusting the analysis to the provided experimental data, and automatically constructs the equations of model according to a list of reactions.


## Key features

- simulation of concentrations of 13C isotopomers originated from artificially 13C enriched substrates, obtained using mass spectrometry and corrected for natural occurring isotopes and peaks overlapping

## Functionality

- Post-processing
- Statistical Analysis
- Workflows

## Approaches

- Isotopic Labelling Analysis / 13C
    
## Data Analysis

- simulation of the mass isotopomer data, fitting them using a kinetic nodel, evaluation of metabolic fluxes corresponding to the best fit

## Instrument Data Types

- MS

## Tool Authors

- Vitaly Selivanov (Universitat de Barcelona)

## Container Contributors

- [Pablo Moreno](EBI)

## Website

- N/A

## Git Repository

-https://github.com/seliv55/wf/tree/master/isodyn

## Installation

- "IsoDyn" does not require installation. To run it in local computer it is sufficient to copy (clone) the repository. 
- If some of .cpp or .h files are modified, run ''' make clean && make '''

## Usage Instructions

A linux shell script "iso.sh" helps to run Isodyn. Here is its text with comments explaining the meaning of its parameters.

#!/bin/sh

fiso="SW620-Glucose"    # input data: file indicating 13C labeling of metabolites

fcon=xglc               # input data: file with measured concentrations

inpar="glc/1"           # input data: initial set of parameters to start

oudir="glc/"            # output directory

fstat="glc/statfl"      # path to write the results of fitting: mean and confidence intervals

fcmpr="glc/statfl"      # results of fitting for the conditions used for comparison.

manfi=77                # number of files to be saved during fitting

FNCKAS="0"              # a number of options forsing Isodyn to run in various modes. Default: make one simulation and stop

tst=yes                 # run Isodyn?

while getopts ":\a:b:i:\o:s:c:\m:FNCKAS" opt; do

  case $opt in
  
    a) fiso=$OPTARG;;
    
    b) fcon=$OPTARG;;
    
    i) inpar=$OPTARG;;
    
    o) oudir=$OPTARG;;
    
    s) fstat=$OPTARG;;
    
    c) fcmpr=$OPTARG;;
    
    m) manfi=$OPTARG;;
    
    F) FNCKAS=F;;     # fit data using Simulated Annealing algorithm
    
    N) FNCKAS=N;;     # find the number of degrees of freedom for estimation of goodness of fit.
    
    C) FNCKAS=C;;     # attempts to increase confidence intervals for fluxes
    
    K) FNCKAS=K;;     # special algorithm for fitting transketolase parameters
    
    A) FNCKAS=A;;     # special algorithm for fitting transaldolase parameters
    
    S) FNCKAS=S;;     # statistics on the results of fitting: mean and confidence intervals for fluxes
    
    *)
    
      echo "Invalid option: -$OPTARG" 
      
      cat help
      
      tst=no
      
      ;;
      
  esac
  
done

if [ $tst = yes ]

then                  # run Isodyn

./isodyn.out $fiso $fcon $inpar $oudir $fstat $fcmpr $manfi $FNCKAS

fi


- The following command forces Isodyn to run with default parameters
 
```
 ./iso.sh
```

- to run it with different options one of the two options should be used, either edit the script, or use the listed above options, for instance the following command can be used to fit data:
 
```
 ./iso.sh -F
```
To take concentrations from different file, for instance "abcd.efg":
```
 ./iso.sh -b abcd.efg
```
 
## The provided examples:
 
The provided file "SW620-Glucose", is an output of R-program Midcor (https://github.com/seliv55/midcor) that corrects mass spectra for natural isotope abundance and peaks overlapping. This output is designed for convenience of visual checking and manual edition, it can be edited or taken directly. It contains the mass isotopomer distributions measured for various metabolites of central pathways. The other infrmation necessary for simulations with Isodyn (path to the file with model parameters; path to output directory; max number of parameter files saved during optimization; path to the file where statistical data are saved) is presented in the file "xglc". Below are the examples of commands to run Isodyn in several modes.

## Publications

- 1: Selivanov VA, Vizán P, Mollinedo F, Fan TW, Lee PW, Cascante M.
Edelfosine-induced metabolic changes in cancer cells that precede the
overproduction of reactive oxygen species and apoptosis. BMC Syst Biol. 2010, 4:135.

- 2: de Mas IM, Selivanov VA, Marin S, Roca J, Orešič M, Agius L, Cascante M.
Compartmentation of glycogen metabolism revealed from 13C isotopologue
distributions. BMC Syst Biol. 2011, 5:175.

- 3: Selivanov VA, Marin S, Lee PW, Cascante M. Software for dynamic analysis of
tracer-based metabolomic data: estimation of metabolic fluxes and their
statistical analysis. Bioinformatics. 2006, 22(22):2806-12.

- 4: Selivanov VA, Meshalkina LE, Solovjeva ON, Kuchel PW, Ramos-Montoya A,
Kochetov GA, Lee PW, Cascante M. Rapid simulation and analysis of isotopomer
distributions using constraints based on enzyme mechanisms: an example from HT29 
cancer cells. Bioinformatics. 2005, 21(17):3558-64.

- 5: Selivanov VA, Puigjaner J, Sillero A, Centelles JJ, Ramos-Montoya A, Lee PW,
Cascante M. An optimized algorithm for flux estimation from isotopomer
distribution in glucose metabolites. Bioinformatics. 2004, 20(18):3387-97. 

