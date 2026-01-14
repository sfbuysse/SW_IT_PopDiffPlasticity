# SW_IT_PopDiffPlasticity

This contains the code used to analyze the 2021 and 2022 Experiments to identify population differentiation for plasticity between Sweden and Italy _Arabidopsis thaliana_ native genotypes. For each script, the .Rmd is where the code was written. It is then exported as a .md for ease of viewing on github and a .html for ease of viewing as a downloaded file. Note some code chunks are hidden in the .html for readability.

00_TrayLayouts shows how pots were randomized for each experiment. It was edited for readability and knitted on 10/17/2025.

01_CleanData fixes any data entry errors and looks at general distributions of each trait in each experiment to plan for traits that might need to be transformed to increase normality of residuals in linear models. It was edited for readability and knitted on 10/17/2025. The raw data was added to the github repo and this code was knitted on 1/14/2026.

02_Analysis includes all linear models and anovas. It was cleaned for readability and some unused analyses excluded on 10/21/2025. Was last knitted 1/14/2026.

03_Figures makes all figures at 600 dpi. It was cleaned for readability and some unused figures removed on 10/21/2025. Was last knitted 1/14/2026.

You can find this manuscript on biorxiv: https://doi.org/10.64898/2026.01.07.698234