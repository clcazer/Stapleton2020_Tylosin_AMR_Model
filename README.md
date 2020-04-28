#Modeling the effect of tylosin on macrolide-resistant enterococci in feedlots and reducing resistance transmission

The R Analysis folder provides input for the Matlab model and analyzes the output from the Matlab model. 
Within the scripts folder, "TYL model parameter distributions.R" is used before the Matlab model to determine the model parameter distributions. It uses two excel files in the data folder: "All Parameters + Citations_cc_2_13_20.xlsx" and "NARMS-Data-Cecal-2013-2017.xlsx".

The model is run in MatLab (Version R2019a), first using "TYL_Model_burn_in_calc.m" to determine the burn-in period for reaching equilibrium in the proportion of resistant enterococci in cattle. Then "master.m" is used to run all model scenarios in "TYL_Model_cc.m". Within "master.m", the filename variable may need editing to determine the saved location of results.

Results from MatLab (txt files) should be saved or copied to the R Analysis/data folder. Then analysis and figures can be created using the remaining R scripts. Additional data files used for validation are "Extracted Data_long form_meta.xlsx" (data from https://doi.org/10.1016/j.prevetmed.2020.104934) and "Schmidt validation data.xlsx".

---------------
General information

Title of the paper: Modeling the Effect of Tylosin on Macrolide-Resistant Enterococci in Feedlots and Reducing Resistance Transmission
Journal:
Authors and ORCIDs:
-G. Sean Stapleton:
-Casey L Cazer: 0000-0002-2290-1868
-Yrjo T. Grohn: 0000-0002-2721-9833

Corresponding author: G. Sean Stapleton, gss94@cornell.edu

---------------
Code Citation:
Cazer, Casey L., and Stapleton G. Sean. (2020). "Replication Package for Modeling the Effect of Tylosin on Macrolide-Resistant Enterococci in Feedlots and Reducing Resistance Transmission" (version 1) [R files and MatLab files]. DOI: 10.5281/zenodo.3724925
Code last updated: April 28, 2020

Code corresponding author: Casey L. Cazer, clc248@cornell.edu

Software and Version:
MatLab Version R2019a, run on Windows Server 2012 R2 Datacenter, 64-bit
R Version 3.6.1, run on Windows 10, 64-bit
