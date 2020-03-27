# Isodyn
Version: 1.0

![Logo](text3923.png)

## Short Description

“C++”-program simulating the dynamics of metabolites and their isotopic isomers in central metabolic network using a corresponding kinetic model

## Description

“Isodyn” is a C++-program that supports an analysis of stable isotope tracer data to assess metabolic flux profiles in living cells. Isodyn simulates propagation of 13C label from artificially 13C enriched substrates into metabolites of central metabolic pathways forming various isotopic isomers (isotopomers). It fits the simulated dynamics of mass isotopomers to that observed experimentally, thus finding the parameters, which reflect the characteristics of corresponding biochemical reactions, and metabolic fluxes that correspond to the real fluxes for the analyzed biological object and conditions. Isodyn contains tools that check the goodness of fit and perform a statistical analysis of obtained metabolic fluxes. It uses the change of metabolite concentrations in the analysis in addition to labeling data, as a restrictive factor. Optionally it uses manually corrected files provided by Midcor (https://github.com/seliv55/midcor). It automatically adjusts the analysis to the provided experimental data, and automatically constructs the equations of model according to a list of reactions.


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

```
#!/bin/sh
fiso="hct-glc" #13C labeling 
fcon=hct-conc #concentrations
inpar="output/2" #set of parameters to start
oudir="output"  #output directory
fstat="glc/statfl"  #statistics on the all sets of parameters
fcmpr="glc/statfl1" #statistics for matching conditions
manfi=99 #number of files saved
FNCKAS="0" #-Fit with Simulated annealing, find -Number of independent parameters,
tst=yes

while getopts ":a:b:i:o:s:c:m:FNCKAS" opt; do
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
```


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
 
The provided file "hct-glc", is an output of R-program Midcor (https://github.com/seliv55/midcor) that corrects mass spectra for natural isotope abundance and peaks overlapping. This output is designed for the convenience of visual checking and manual edition, it can be edited or taken directly. It contains the mass isotopomer distributions measured for various metabolites of central pathways. The other infrmation necessary for simulations with Isodyn (path to the file with model parameters; path to output directory; max number of parameter files saved during optimization; path to the file where statistical data are saved) is indicated in the script above.

## change the scheme of simulated reactions:

Users need to provide the file containing a desription of the simulated reactions in a special text format. This file, named 'model', should be presented in the folder 'reactions' that contains also the Python program 'rmod.py'. The latter reads 'model' and converts it into the modules of Isodyn containing ODEs for time course of total concentrations and isotopomers.

In the first part the file 'model' presents metabolites, which total and isotopomer concentrations are variables of the model:

Metab 3  npyr "Pyr"

Metab 2  ncoa  "CoA"

Metabolites are presented as a set of isotopomers. They are objects of a class, whos name is shown in the first column. Next in each row is the number of atoms in the carbon skeleton, then a word indicating the index of a given metabolite in the corresponding array, and finally the name of the given metabolite.
Internal and external metabolites are separated.

Then the simulated reactions are presented:

 lacin    y[nlac] ->  nlacc $  0 input lac lacc 0     (123->123)

 laccpyr  nlacc ->    npyr $   0 input lacc pyr 0     (123->123)

 cs0      nmal ncoa ->  ncit $  1 r3met condence cit mal coa 0 (12+3456->123456)
 
Each line starts from the name of reaction, then a list of substrates terminated by '->', then the list of products terminated by '$', flag (0 or 1) indicating reactions that the program can change by the same factor (if the flag is 1) for convenience of fitting. Then either the name of C++ function that simulates carbon trabsition in the reaction, or the name of a group of functions. In the latter case this name follows by a name of particular function in the group. Then the list of substrates and products, which are the parameters of the indicated function. Then the name of the reverse flux or 0 if the transformation is considered as irreversible. Finally the scheme of atom transition is chown.

Based on this description, the program 'rmod.py' constructs modules of Isodyn substituting the existing ones in the files 'nv.cpp', 'distr.cpp' (in the folder 'con512tpl') and 'nums.hh' (in the folder 'include'). Also it saves a file 'parameters' using the same format as the files with model parameters im the folder 'output', but corresponding to the modified scheme. The used old parameters file should be compared with the new generated file and the lines in the old file corresponding to the removed reactions should be removed and the lines corresponding to the added reactions should be provided. The same should be done for the list of metabolites with their initial values that are also presented in the parameters file. Then after recompilation ( make clean && make ) Isodyn uses the modifyed scheme of reactions. 


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

