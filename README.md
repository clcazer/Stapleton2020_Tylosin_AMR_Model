#Modeling the effect of tylosin on macrolide-resistant enterococci in feedlots and reducing resistance transmission

The R Analysis folder provides input for the Matlab model and analyzes the output from the Matlab model. 
Within the scripts folder, "TYL model parameter distributions.Rmd" is used before the Matlab model to determine the model parameter distributions. It uses two excel files in the data folder: "All Parameters + Citations_cc_2_13_20.xlsx" and "NARMS-Data-Cecal-2013-2017.xlsx".

The model is run in MatLab, first using "TYL_Model_burn_in_calc.m" to determine the burn-in period for reaching equilibrium in the proportion of resistant enterococci in cattle. Then "master.m" is used to run all model scenarios in "TYL_Model_cc.m". Within "master.m", the filename variable may need editing to determine the saved location of results.

Results (txt files) should be saved or copied to the R Analysis/data folder. Then analysis and figures can be created using the remaining R and Rmd scripts. Additional data files used for validation are "Extracted Data_long form_meta.xlsx" (data from https://doi.org/10.1016/j.prevetmed.2020.104934) and "Schmidt validation data.xlsx".
