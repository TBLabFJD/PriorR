# PriorR

Priorr is a prioritization program of disease-linked genetic variants devoloped within the Genetics&Genomics Department of La Fundacion Jimenez Diaz University Hospital. Priorr is conceived to analyse the output of the FJD-pipeline of SNVs or CNVs. This software program offers a number of useful functionalities for variant analysis such as: filtering by a virtual panel of genes. manual control of different population frequencies or pathogenicity predictors or filtering out variants that have been already found by another protocol.  

# License

VariantCallingFJD source code is provided under the [**Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)**](https://creativecommons.org/licenses/by-nc-sa/4.0/). VariantCallingFJD includes several third party packages provided under other open source licenses, please check them for additional details.

[![Licencia de Creative Commons](https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png)](http://creativecommons.org/licenses/by-nc-sa/4.0/)


# Dependencies

**Programming languages**

- **R v.3.6.1**

**R libraries:**

- shiny_1.3.2 
- DT_0.9
- data.table_1.14.0
- dplyr_1.0.7
- shinyWidgets_0.4.9
- lazyeval_0.2.2
- shinydashboard_0.7.1
- shinydashboardPlus_0.7.0
- shinyalert_1.0
- filesstrings_3.1.5
- shinyBS_0.61
- shinyjs_1.0
- optparse_1.6.4 


# Installation

Install dependencies

Clone code of Priorr using 
```sh
git clone https://github.com/TBLabFJD/PriorR.git
```

# Usage

## Preparation

Prepare the following folders

1. Program_files

This directory should be strucutured as follows:

<pre>
Program_files
    |
    |-- Users
    |-- Database
    |-- Dependencies
    |-- Paneles_Orphanet
</pre>

Subfolders explanation and content:

- Users: this folder contains the user information about the session. Should be created but empty when the program is installed.
- Database: this folder contains the internal database of variants of interest as they are cataloghed by the analysts. Should contain an empty table with column names (See example). 
- Dependencies: this folder contains the file that links genes to HPO files needed to find genes relevant for a disease (see uploaded files).
- Orphanet panels: this folder contains all the Orphanet panels (see uploaded files). 


2. Sample_files: should contain the variant files as they come out from the FJD-pipeline.

It should have the following structure:
<pre>
Samples_files
   |
   |-- Samples
   |-- Diff
</pre>
- Samples: contains samples to analyse in pvm.txt format (see example).
- Diff: contains diff files that indicate which variants are in the vcf commercial and which are not (see example).

## Running

**Windows system**:

Create a bat file to launch Priorr from windows (PRIORR.bat) and save it in your desktop:

Example:
```sh
Rscript="C:\Program Files\R\R-3.6.1\bin\Rscript.exe" # Route to Rscript.exe
%Rscript% \Route_to_Priorr\PRIORR2.3.R -p Route_to_program_files -s Route_to_Samples
```
**Linux system:**

Launch the following line in linux command line:
```sh
Rscript  /Route_to_Priorr/PRIORR2.3.R -p Route_to_program_files -s Route_to_Samples
```
## Intructions

Execute PriorR. 

Upload a variant file in format pvm.txt (see example) and modify the filters.



