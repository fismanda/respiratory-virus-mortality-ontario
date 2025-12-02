********************************************************************************
* META-ANALYSIS OF POPULATION ATTRIBUTABLE RISK ESTIMATES
* Assessment of heterogeneity by time period and seasonal adjustment method
********************************************************************************
* This code performs random-effects meta-analysis to assess the contribution
* of pandemic period and Fourier seasonal adjustment (FFT) to heterogeneity
* in population attributable risk (PAR) estimates
********************************************************************************

clear all
set more off

********************************************************************************
* 1. LOAD META-ANALYSIS DATA
********************************************************************************

* Import PAR estimates from primary analysis
import delimited "data/AJR_and_data_for_meta-regression.csv", clear

* Verify data structure
describe
list in 1/5

********************************************************************************
* 2. CREATE LABELS FOR FOREST PLOTS
********************************************************************************

* Already present in dataset:
* - fft_label: "yes" or "no" 
* - pandemic_label: "yes" or "no"
* - label: virus | FFT: [yes/no] | Pandemic: [yes/no]
* - subgroup: Fast Fourier Transform: [No/Yes] | Pandemic: [No/Yes]

********************************************************************************
* 3. OVERALL META-ANALYSIS BY VIRUS TYPE
********************************************************************************
* Forest plot showing all estimates stratified by virus type

metan par lcl ucl, effect(par) random by(virus) nooverall ///
    label(namevar=subgroup) xtitle("PAR (%)") ///
    xlabel(0,0.05,0.1,0.15,0.20) textsize(80) scheme(s1color)

graph export "figures/figure3_forest_plot.png", replace width(3000)

********************************************************************************
* 4. HETEROGENEITY ANALYSIS - ALL RESPIRATORY VIRUSES
********************************************************************************

* Overall heterogeneity
metan par lcl ucl if virus == "All Respiratory Viruses"

* By FFT status
metan par lcl ucl if virus == "All Respiratory Viruses" & fft == 1 
metan par lcl ucl if virus == "All Respiratory Viruses" & fft == 0

* By pandemic period
metan par lcl ucl if virus == "All Respiratory Viruses" & pandemic == 0
metan par lcl ucl if virus == "All Respiratory Viruses" & pandemic == 1

********************************************************************************
* 5. HETEROGENEITY ANALYSIS - INFLUENZA A
********************************************************************************

* Overall heterogeneity
metan par lcl ucl if virus == "Influenza A"

* By FFT status
metan par lcl ucl if virus == "Influenza A" & fft == 1
metan par lcl ucl if virus == "Influenza A" & fft == 0

* By pandemic period
metan par lcl ucl if virus == "Influenza A" & pandemic == 0
metan par lcl ucl if virus == "Influenza A" & pandemic == 1

********************************************************************************
* 6. HETEROGENEITY ANALYSIS - RESPIRATORY SYNCYTIAL VIRUS
********************************************************************************

* Overall heterogeneity
metan par lcl ucl if virus == "Respiratory Syncitial Virus"

* By FFT status
metan par lcl ucl if virus == "Respiratory Syncitial Virus" & fft == 1 
metan par lcl ucl if virus == "Respiratory Syncitial Virus" & fft == 0

* By pandemic period
metan par lcl ucl if virus == "Respiratory Syncitial Virus" & pandemic == 0
metan par lcl ucl if virus == "Respiratory Syncitial Virus" & pandemic == 1

********************************************************************************
* 7. HETEROGENEITY ANALYSIS - SARS-CoV-2 (pandemic only)
********************************************************************************

* Overall heterogeneity (by FFT only, as SARS-CoV-2 only exists in pandemic)
metan par lcl ucl if virus == "SARS-CoV-2"

* By FFT status
metan par lcl ucl if virus == "SARS-CoV-2" & fft == 1
metan par lcl ucl if virus == "SARS-CoV-2" & fft == 0

********************************************************************************
* 8. HETEROGENEITY ANALYSIS - INFLUENZA B
********************************************************************************

* Overall heterogeneity
metan par lcl ucl if virus == "Influenza B"

* By FFT status
metan par lcl ucl if virus == "Influenza B" & fft == 1
metan par lcl ucl if virus == "Influenza B" & fft == 0

* By pandemic period
metan par lcl ucl if virus == "Influenza B" & pandemic == 0
metan par lcl ucl if virus == "Influenza B" & pandemic == 1

********************************************************************************
* END OF META-ANALYSIS
********************************************************************************
* Key findings:
* - Influenza A: Heterogeneity driven by BOTH FFT and pandemic period (I²=93.7%)
* - RSV: Heterogeneity driven by BOTH FFT and pandemic period (I²=88.5%)
* - SARS-CoV-2: Low heterogeneity, driven only by FFT (I²=2.5%)
*   (Pandemic period not applicable - SARS-CoV-2 only present in pandemic)
********************************************************************************
