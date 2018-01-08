# Isodyn
Version: 1.0

![Logo](text3923.png)

## Short Description

“C++”-program simulating the dynamics of metabolites and their isotopic isomers in central metabolic network using a corresponding kinetic model

## Description

“Isodyn” is a C++-program that performs an analysis of stable isotope tracer data to assess metabolic flux profiles in living cells. Isodyn simulates propagation of 13C label from artificially 13C enriched substrates into isotopic isomers (isotopomers) of metabolites of central metabolic pathways. It performs fitting the simulated dynamics of mass isotopomers to that observed experimentally, thus finding the parameters, which reflect the characteristics of corresponding biochemical reactions, and metabolic fluxes that correspond to the real fluxes under the analyzed biological object and conditions. Isodyn contains tools that check the goodness of fit and perform a statistical analysis of obtained metabolic fluxes.

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

Isodyn can run in several modes:

- The following command forces Isodyn to perform just one simulation of the metabolite labeling data presented in "labeling_data_file" and external concentrations and other input information (path to the file with model parameters; path to output directory; max number of parameter files saved during optimization; path to the file where statistical data are saved) presented in "additional_info_file" and stop:
 
```
  ./isodyn.out [labeling_data_file] [additional_info_file] 
```

- The addition of parameter "s" forces to calculate the confidence intervals for all fluxes, based on the previously saved in the output directory files with model parameters and metabolic fluxes, and save the results in file "statfl":
 
```
 ./isodyn.out [labeling_data_file] [additional_info_file] s 
```

- The addition of parameter "x" forces to recalculate the χ2 for the parameters sets previously saved in the output directory as files "1", "2", etc:
 
```
 ./isodyn.out [labeling_data_file] [additional_info_file] x 
```

- The addition of integer forces to perform optimization using Simulated Annealing algorithm minimizing χ2 and stop after saving "int_number" of files with optimized parameters in the output directory:
 
```
 ./isodyn.out [labeling_data_file] [additional_info_file] [int] 
```

- find a set of parameters, which change produces linearly independent changes iin the output corresponding to the experimental data:
 
```
 ./isodyn.out [labeling_data_file] [additional_info_file] g 
```

- optimise parameters for transketolase reaction:
 
```
 ./isodyn.out [labeling_data_file] [additional_info_file] tk 
```

- optimise parameters for transaldoase reaction:
 
```
 ./isodyn.out [labeling_data_file] [additional_info_file] ta 
```


 
## The provided examples:
 
The provided file "A549", is obtained applying tools supporting a workflow of primary analysis of raw data machine written in multipeak CDF files (R-programs RaMID or CDF2MID, and MIDcor). It contains the mass isotopomer distributions measured for various metabolites of central pathways. The other infrmation necessary for simulations with Isodyn (path to the file with model parameters; path to output directory; max number of parameter files saved during optimization; path to the file where statistical data are saved) is presented in the file "xglc". Below are the examples of commands to run Isodyn in several modes.

Single simulation:

```
 ./isodyn.out A549 xglc 
```

screenshot of a simulation of input data, shown only for unlabeled fraction (m0)

![screenshot](Screenshot.png)

 
- the other modes of Isodyn functioning are achieved by addition of one more parameter to the above command as described in "Usage instruction"

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

