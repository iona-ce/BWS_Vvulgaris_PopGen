# BWS_Vvulgaris_PopGen
Code used for the BWS Population Genetics study (currently under review) and other information (Data, Figures, Results). 

Data is also stored here: 

# In this GitHub repo:

## [Scripts](https://github.com/iona-ce/BWS_Vvulgaris_PopGen/tree/main/Scripts) 

[01_TidyingData_50Complete.R](https://github.com/iona-ce/BWS_Vvulgaris_PopGen/blob/main/Scripts/01_TidyingData_50Complete.R) : R code to clean raw data for subsequent analyses.

[02_ColonyDataPrep.R](https://github.com/iona-ce/BWS_Vvulgaris_PopGen/blob/main/Scripts/02_ColonyDataPrep.R): R code to prepare data for use with [COLONY](https://www.zsl.org/about-zsl/resources/software/colony) (Jones and Wang, 2019)

[03_RunningColony_Terminal.txt](https://github.com/iona-ce/BWS_Vvulgaris_PopGen/blob/main/Scripts/03_RunningColony_Terminal.txt) : Code to run COLONY in the terminal. 

[04_ColonyOutputs.R](https://github.com/iona-ce/BWS_Vvulgaris_PopGen/blob/main/Scripts/04_ColonyOutputs.R): R code to clean COLONY outputs; sibling analyses; sibling removals. 

[05_Using_GenAlEx.txt](https://github.com/iona-ce/BWS_Vvulgaris_PopGen/blob/main/Scripts/05_Using_GenAlEx.txt): Instructions on using [GenAlEx](https://biology-assets.anu.edu.au/GenAlEx/Welcome.html). 

[06_GenePop.R](https://github.com/iona-ce/BWS_Vvulgaris_PopGen/blob/main/Scripts/06_GenePop.R): R code to analyse data with `GenePop`: Hardy Weinberg, Linkage Disequilibrium, Fis. 

[07_StructurePrep.R](https://github.com/iona-ce/BWS_Vvulgaris_PopGen/blob/main/Scripts/07_StructurePrep.R): R code to prepare data for use with [STRUCTURE](https://web.stanford.edu/group/pritchardlab/structure.html) 

[08a_Structure_Terminal_2017.txt](https://github.com/iona-ce/BWS_Vvulgaris_PopGen/blob/main/Scripts/08a_Structure_Terminal_2017.txt): Code to run STRUCTURE in the terminal (including mainparams file edits) for 2017 data (national scale).

[08b_Structure_Terminal_2018.txt](https://github.com/iona-ce/BWS_Vvulgaris_PopGen/blob/main/Scripts/08b_Structure_Terminal_2018.txt): Same but for 2018 (regional scale) data. 

[09_PopHelper_StructureVisualisation.R](https://github.com/iona-ce/BWS_Vvulgaris_PopGen/blob/main/Scripts/09_PopHelper_StructureVisualisation.R): R code to visualise STRUCTURE outputs using `pophelper` and implement Evanno's method. 

[10_GenindObjects.R](https://github.com/iona-ce/BWS_Vvulgaris_PopGen/blob/main/Scripts/10_GenindObjects.R): R code to transform genetic data into genind objects.

[11_Mantel_Fst.R](https://github.com/iona-ce/BWS_Vvulgaris_PopGen/blob/main/Scripts/11_Mantel_Fst.R): R code to perform calculations of Fst and isolation by distance using Mantel tests.

[12_MapReprojection.R](https://github.com/iona-ce/BWS_Vvulgaris_PopGen/blob/main/Scripts/12_MapReprojection.R): R code to change map (NUTS areas) CRS.

[13_Mapping.R](https://github.com/iona-ce/BWS_Vvulgaris_PopGen/blob/main/Scripts/13_Mapping.R): R code to map data points. 
