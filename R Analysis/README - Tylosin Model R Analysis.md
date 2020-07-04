# Tylosin-Model

This R project supports TYL_Model (matlab script). It analyzes esimates distributions of parameters used in the matlab model, then takes output from the model (.txt files in "data/") and performs a sensitivity analysis, creates figures (saved in "figures/", and returns descriptive analyses (saved in "results/"). Each script is intended to be run independently.

scripts:
"TYL model parameter distributions.Rmd": calculates best-fit distributions for mathematical model parameters from literature. Uses 'All Parameters + Citations_cc_2-13-20.xlsx" and "NARMS-Data-Cecal-2013-2017.csv" as data sources. Saves model parameter distributions (hyperparameters) to "results/parameter distributions.txt".

"TYL func boxplots.R": function to create a "functional boxplot", showing specficic percentiles of TYL model results over time

"TYL intervention vs CON intervention prop res cow.R": analyzing the proportion of resistant enterococci in cattle, comparing TYL intervention groups to CON intervention groups. Also analyzes the concentration of TYL in the large intestine and validation against feeding trials. Descriptive results are saved to "results/tyl intestine concentrations.txt", "results/validation.txt", "results/TYL effect.txt". Figures 2, 3, 4 are created and saved to "figures/".

"TYL intervention vs CON intervention prop res env.R": analyzing the proportion of resistant enterococci in the environment (feed, water, pen), comparing TYL intervention groups to CON intervention groups. Supplementary Figure 1 saved to "figures/".

"TYL intervention vs CON intervention ent conc.R": analyzing the concentration of resistant enterococci in cattle, pen, feed and water, comparing TYL intervention groups to CON intervention groups. Supplementary Figure 2 saved to "figures/". Descriptive results saved to "results/TYL effect on enterococci concentrations.txt".

"Intervention vs NI prop res cow.R": analyzing the proprtion of resistant enterococci in cattle, comparing intervention to no intervention groups (e.g. TYL intervention vs TYL no intervention; CON intervention vs CON no intervention). Supplementary Figure 3 saved to "figures/". Descriptive results saved to "results/intervention effect.txt".

"TYL sensitivity analysis.R": a sensitivity analysis of the TYL model outputs (from TYL no intervention scenario) to model parameters. Figure 5 saved to "figures/". Returns realized range of model parameters to "results/param_range.xlsx".

"TYL intervention vs CON intervention prop res cow_LG Pen.R": analyzing the proportion of resistant enterococci in cattle, comparing TYL intervention groups to CON intervention groups, in the model scenario with a large pen size (N=150 cattle). Also analyzes the concentration of TYL in the large intestine. Does not assess validation because the validation studies use a small pen size. Descriptive results are saved to "results/tyl intestine concentrations LGpen.txt", "results/TYL effect LGpen.txt". LG Pen versions of Figures 2 and 4 are created and saved to "figures/".

--------------------------------
$platform
[1] "i386-w64-mingw32"

$arch
[1] "i386"

$os
[1] "mingw32"

$system
[1] "i386, mingw32"

$status
[1] ""

$major
[1] "3"

$minor
[1] "6.1"

$year
[1] "2019"

$month
[1] "07"

$day
[1] "05"

$`svn rev`
[1] "76782"

$language
[1] "R"

$version.string
[1] "R version 3.6.1 (2019-07-05)"

$nickname
[1] "Action of the Toes"

--------------------------------