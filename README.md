Respiratory Virus Mortality in Ontario, 1993-2024
Overview
This repository contains data and code for ecological regression analysis of population attributable mortality from respiratory viruses in Ontario, Canada over a 31-year period (1993-2024).
Citation
Manuscript: [Citation will be added upon publication]
Data and Code Repository: [DOI to be added]
If you use this data or code, please cite both the manuscript and this repository.
Study Summary
We estimated the population attributable fraction (PAF) of deaths due to influenza A, influenza B, respiratory syncytial virus (RSV), and SARS-CoV-2 in Ontario using ecological regression with negative binomial models. The analysis spans two distinct periods:
Pre-pandemic period (1993-2020): Influenza A accounted for 1.8% of all deaths (95% CI: 1.4-2.2%)
Pandemic period (2020-2024): SARS-CoV-2 accounted for 6.5% of all deaths (95% CI: 4.5-8.4%), representing a 3.6-fold higher burden than pre-pandemic influenza A
A key methodological finding was that Fourier seasonal adjustment substantially affected mortality estimates for RSV, revealing either paradoxical negative associations (PAR -0.45%) with seasonal adjustment or consistent positive associations (PAR 1.9%) without adjustment across both periods.
Test-adjusted SARS-CoV-2 incidence was estimated as described in Fisman DN, Greer AL, Brankston G, Hillmer M, O'Brien SF, Drews SJ, Tuite AR. COVID-19 Case Age Distribution: Correction for Differential Testing by Age. Ann Intern Med. 2021 Oct;174(10):1430-1438. https://doi.org/10.7326/M20-7003, and Bosco, S., Peng, A., Tuite, A.R. et al. Impact of adjustment for differential testing by age and sex on apparent epidemiology of SARS-CoV-2 infection in Ontario, Canada. BMC Infect Dis 25, 589 (2025). https://doi.org/10.1186/s12879-025-10968-6.
Repository Contents
respiratory-virus-mortality-ontario/
├── README.md
├── data/
│   ├── monthly_mortality_virology_data_ontario_1991_2024.csv
│   ├── meta_analysis_par_estimates.csv
│   └── DATA_DICTIONARY.md
├── code/
│   ├── 01_analysis.do
│   ├── 02_meta_analysis.do
│   └── 03_figures.R
├── figures/
└── manuscript/
Reproducing the Analysis
Requirements
Software:
Stata 17 or higher (for regression analysis)
Required packages: punaf, metan, suest
R 4.0 or higher (for figure generation)
Required packages: haven, dplyr, tidyr, ggplot2, patchwork
Step 1: Run Stata Analysis
Navigate to the repository directory and run the analysis scripts:
do code/01_analysis.do
do code/02_meta_analysis.do
Step 2: Generate Figures in R
source("code/03_figures.R")
Expected Runtime
Stata analysis: ~5-10 minutes
R figure generation: ~1-2 minutes
Data Sources
Mortality Data
Ontario Deaths Registry (1994-2024): Monthly death counts, Ontario Ministry of Health
Statistics Canada (1991-1993): Table 13-10-0392-01, Deaths and mortality rates by month
Virological Surveillance
FluWatch (1993-2024): Influenza A, B, and RSV percent positivity from sentinel laboratories. Source: Public Health Agency of Canada
SARS-CoV-2 Data
Adjusted case counts (March 2020 - August 2022): Test-adjusted case counts from Ontario COVID-19 Science Advisory Table. Methodology described in Fisman et al. (2021) Ann Intern Med 174(10):1430-1438 and Bosco et al. (2025) BMC Infect Dis 25:589
Percent positivity (September 2022 - March 2024): SARS-CoV-2 percent positivity from FluWatch sentinel surveillance
Methods Summary
Statistical Approach:
Negative binomial regression (to account for overdispersion)
Population offset (log-transformed)
Fourier seasonal adjustment (sine and cosine terms)
Stratification by pre-pandemic vs. pandemic period
Sensitivity analysis without Fourier adjustment
Wald tests to assess impact of seasonal adjustment
Key Findings
1. Pre-pandemic baseline: Influenza A accounted for 1.8% of deaths (95% CI 1.4-2.2%)
2. SARS-CoV-2 dominance: During the pandemic, SARS-CoV-2 accounted for 6.5% of deaths (95% CI 4.5-8.4%), 3.6-fold higher than pre-pandemic influenza A despite widespread vaccination and antiviral availability
3. Validation: Model-estimated COVID-19 deaths (18,052) closely matched reported COVID-19 deaths (18,603), representing 97% concordance
4. RSV paradox: With FFT: PAR -0.45% pre-pandemic, 0.78% pandemic. Without FFT: PAR 1.9% both periods. Interpretation: Seasonal adjustment may remove true RSV signal
5. Heterogeneity: Influenza A and RSV: High heterogeneity (I² >88%) driven by both time period and FFT method. SARS-CoV-2: Low heterogeneity (I² 2.5%) driven only by FFT method
Authors
David N. Fisman, MD, MPH, FRCPC Dalla Lana School of Public Health, University of Toronto
[Additional authors to be added]
Contact
For questions about the analysis or code:
David Fisman: david.fisman@utoronto.ca
License
This project is licensed under the MIT License - see the LICENSE file for details.
Acknowledgments
We acknowledge the Ontario Ministry of Health, Statistics Canada, and the Public Health Agency of Canada for providing the data that made this analysis possible. We thank the Ontario COVID-19 Science Advisory Table for development of testing-adjusted case count methodology.
Version History
v1.0 (December 2025): Initial release with manuscript submission
Last updated: December 2025
