********************************************************************************
* RESPIRATORY VIRUS MORTALITY ANALYSIS - ONTARIO 1991-2026
* Population Attributable Risk Estimation Using Ecological Regression
* EXTENDED ANALYSIS WITH THREE TIME PERIODS - VERSION 4
********************************************************************************
* This code estimates the fraction of deaths in Ontario attributable to 
* influenza A, influenza B, RSV, and SARS-CoV-2 using negative binomial 
* regression with Fourier seasonal adjustment across three periods:
* 1. Pre-pandemic (Jan 1993 - Feb 2020)
* 2. PHEIC (March 2020 - April 2023) 
* 3. Post-PHEIC (May 2023 - Feb 2025)
*
* VERSION 4 CHANGES FROM V3:
* - Added Section 8: Sensitivity analysis excluding Mar-Jun 2020
*   (reviewer pre-emption strategy; tests robustness of combined PAF
*   to early chaotic period with healthcare disruption and coding issues)
*
* VERSION 3 APPROACH (MATCHES ORIGINAL PAPER):
* - Use RAW percent positivity for Flu A, Flu B, RSV, SARS % pos
* - Rescale SARS-CoV-2 adjusted cases by dividing by 10,000
* - In PHEIC model: Include BOTH adj_cases_10k AND monthly_sars_pct_pos
* - For SARS PAF: Set BOTH to zero in counterfactual
********************************************************************************

clear all
set more off

********************************************************************************
* 1. DATA PREPARATION
********************************************************************************

* Load monthly aggregated data
import delimited "/Users/davidfisman_new_mac/Dropbox/Family Room/Attributable Mortality/Updated/ontario_mortality_virology_EXTENDED_1991_2025_CORRECTED.csv", clear

* Destring dates
gen date2 = date(date, "20YMD")
format date2 %d
drop date
ren date2 date

* Create analysis variables
gen series_month = month + 12*(year-1991)

* Fourier transform terms for seasonal adjustment
gen sinmo = sin(6.28*month/12)
gen cosmo = cos(6.28*month/12)

* Define time periods
gen prepandemic = (date < date("March 1, 2020", "MDY"))
gen pheic = (date >= date("March 1, 2020", "MDY") & date <= date("April 30, 2023", "MDY"))
gen postpheic = (date >= date("May 1, 2023", "MDY"))

* Labels for clarity
label define period_lab 0 "Pre-pandemic" 1 "PHEIC" 2 "Post-PHEIC"
gen period = 0 if prepandemic == 1
replace period = 1 if pheic == 1  
replace period = 2 if postpheic == 1
label values period period_lab

********************************************************************************
* 2. VARIABLE PREPARATION - NO NORMALIZATION, JUST RAW VALUES
********************************************************************************

* Rescale SARS-CoV-2 adjusted cases for interpretable coefficients
gen adj_cases_10k = monthly_adj_cases/10000

* Replace missing values with zeros
replace monthly_sars_pct_pos = 0 if monthly_sars_pct_pos == .
replace adj_cases_10k = 0 if adj_cases_10k == .

* Restrict to confirmed complete mortality data (through February 2025)
keep if date <= date("February 28, 2025", "MDY")

display _newline
display "USING RAW PERCENT POSITIVITY VALUES (NO NORMALIZATION)"
display "SARS-CoV-2 adjusted cases rescaled by dividing by 10,000"
display _newline

********************************************************************************
* 3. PRE-PANDEMIC ANALYSIS (Jan 1993 - Feb 2020)
********************************************************************************

display _newline(2)
display "********************************************************************************"
display "PRE-PANDEMIC PERIOD: January 1993 - February 2020"
display "********************************************************************************"

* Long-term trend (centered on pre-pandemic midpoint)
gen linear_pre = series_month - 168 if prepandemic == 1
gen mosq_pre = linear_pre^2 if prepandemic == 1

*--- Model WITH Fourier seasonal adjustment ---
display _newline
display "MODEL WITH FOURIER SEASONAL ADJUSTMENT"
display "----------------------------------------"

nbreg deaths sinmo cosmo linear_pre mosq_pre ///
    monthly_flua_pct_pos monthly_flub_pct_pos monthly_rsv_pct_pos ///
    if prepandemic == 1, ///
    exposure(population) irr

estimates store prepan_fft

* Population Attributable Fractions
display _newline "Influenza A PAF:"
punaf if prepandemic == 1, atspec(monthly_flua_pct_pos=0) eform

display _newline "Influenza B PAF:"
punaf if prepandemic == 1, atspec(monthly_flub_pct_pos=0) eform

display _newline "RSV PAF:"
punaf if prepandemic == 1, atspec(monthly_rsv_pct_pos=0) eform

*--- Model WITHOUT Fourier (sensitivity) ---
display _newline
display "MODEL WITHOUT FOURIER (SENSITIVITY ANALYSIS)"
display "---------------------------------------------"

nbreg deaths linear_pre mosq_pre ///
    monthly_flua_pct_pos monthly_flub_pct_pos monthly_rsv_pct_pos ///
    if prepandemic == 1, ///
    exposure(population) irr

estimates store prepan_nofft

display _newline "Influenza A PAF (no Fourier):"
punaf if prepandemic == 1, atspec(monthly_flua_pct_pos=0) eform

display _newline "Influenza B PAF (no Fourier):"
punaf if prepandemic == 1, atspec(monthly_flub_pct_pos=0) eform

display _newline "RSV PAF (no Fourier):"
punaf if prepandemic == 1, atspec(monthly_rsv_pct_pos=0) eform

********************************************************************************
* 4. PHEIC PERIOD ANALYSIS (March 2020 - April 2023)
********************************************************************************

display _newline(2)
display "********************************************************************************"
display "PHEIC PERIOD: March 2020 - April 2023"
display "********************************************************************************"

* Long-term trend for PHEIC period
gen linear_pheic = series_month - 377 if pheic == 1
gen mosq_pheic = linear_pheic^2 if pheic == 1

*--- Model WITH Fourier seasonal adjustment ---
display _newline
display "MODEL WITH FOURIER SEASONAL ADJUSTMENT"
display "----------------------------------------"
display "NOTE: Model includes BOTH adj_cases_10k AND monthly_sars_pct_pos"
display "      adj_cases_10k has values March 2020-Aug 2022"
display "      monthly_sars_pct_pos has values Sept 2022-April 2023"
display _newline

nbreg deaths sinmo cosmo linear_pheic mosq_pheic ///
    monthly_flua_pct_pos monthly_flub_pct_pos monthly_rsv_pct_pos ///
    adj_cases_10k monthly_sars_pct_pos ///
    if pheic == 1, ///
    exposure(population) irr

estimates store pheic_fft

* Population Attributable Fractions
display _newline "Influenza A PAF:"
punaf if pheic == 1, atspec(monthly_flua_pct_pos=0) eform

display _newline "Influenza B PAF:"
punaf if pheic == 1, atspec(monthly_flub_pct_pos=0) eform

display _newline "RSV PAF:"
punaf if pheic == 1, atspec(monthly_rsv_pct_pos=0) eform

display _newline "SARS-CoV-2 PAF (setting BOTH measures to zero):"
punaf if pheic == 1, atspec(adj_cases_10k=0 monthly_sars_pct_pos=0) eform

*--- Model WITHOUT Fourier (sensitivity) ---
display _newline
display "MODEL WITHOUT FOURIER (SENSITIVITY ANALYSIS)"
display "---------------------------------------------"

nbreg deaths linear_pheic mosq_pheic ///
    monthly_flua_pct_pos monthly_flub_pct_pos monthly_rsv_pct_pos ///
    adj_cases_10k monthly_sars_pct_pos ///
    if pheic == 1, ///
    exposure(population) irr

estimates store pheic_nofft

display _newline "Influenza A PAF (no Fourier):"
punaf if pheic == 1, atspec(monthly_flua_pct_pos=0) eform

display _newline "Influenza B PAF (no Fourier):"
punaf if pheic == 1, atspec(monthly_flub_pct_pos=0) eform

display _newline "RSV PAF (no Fourier):"
punaf if pheic == 1, atspec(monthly_rsv_pct_pos=0) eform

display _newline "SARS-CoV-2 PAF (no Fourier, setting BOTH measures to zero):"
punaf if pheic == 1, atspec(adj_cases_10k=0 monthly_sars_pct_pos=0) eform

********************************************************************************
* 5. POST-PHEIC PERIOD ANALYSIS (May 2023 - Feb 2025)
********************************************************************************

display _newline(2)
display "********************************************************************************"
display "POST-PHEIC PERIOD: May 2023 - February 2025"
display "********************************************************************************"

* Long-term trend for post-PHEIC period
gen linear_post = series_month - 406 if postpheic == 1
gen mosq_post = linear_post^2 if postpheic == 1

*--- Model WITH Fourier seasonal adjustment ---
display _newline
display "MODEL WITH FOURIER SEASONAL ADJUSTMENT"
display "----------------------------------------"

nbreg deaths sinmo cosmo linear_post mosq_post ///
    monthly_flua_pct_pos monthly_flub_pct_pos monthly_rsv_pct_pos ///
    monthly_sars_pct_pos ///
    if postpheic == 1, ///
    exposure(population) irr

estimates store postpheic_fft

* Population Attributable Fractions
display _newline "Influenza A PAF:"
punaf if postpheic == 1, atspec(monthly_flua_pct_pos=0) eform

display _newline "Influenza B PAF:"
punaf if postpheic == 1, atspec(monthly_flub_pct_pos=0) eform

display _newline "RSV PAF:"
punaf if postpheic == 1, atspec(monthly_rsv_pct_pos=0) eform

display _newline "SARS-CoV-2 PAF:"
punaf if postpheic == 1, atspec(monthly_sars_pct_pos=0) eform

*--- Model WITHOUT Fourier (sensitivity) ---
display _newline
display "MODEL WITHOUT FOURIER (SENSITIVITY ANALYSIS)"
display "---------------------------------------------"

nbreg deaths linear_post mosq_post ///
    monthly_flua_pct_pos monthly_flub_pct_pos monthly_rsv_pct_pos ///
    monthly_sars_pct_pos ///
    if postpheic == 1, ///
    exposure(population) irr

estimates store postpheic_nofft

display _newline "Influenza A PAF (no Fourier):"
punaf if postpheic == 1, atspec(monthly_flua_pct_pos=0) eform

display _newline "Influenza B PAF (no Fourier):"
punaf if postpheic == 1, atspec(monthly_flub_pct_pos=0) eform

display _newline "RSV PAF (no Fourier):"
punaf if postpheic == 1, atspec(monthly_rsv_pct_pos=0) eform

display _newline "SARS-CoV-2 PAF (no Fourier):"
punaf if postpheic == 1, atspec(monthly_sars_pct_pos=0) eform

********************************************************************************
* 6. SUMMARY COMPARISON ACROSS PERIODS
********************************************************************************

display _newline(2)
display "********************************************************************************"
display "SUMMARY: COMPARISON ACROSS TIME PERIODS (WITH FOURIER ADJUSTMENT)"
display "********************************************************************************"

estimates table prepan_fft pheic_fft postpheic_fft, ///
    keep(monthly_flua_pct_pos monthly_flub_pct_pos monthly_rsv_pct_pos ///
         adj_cases_10k monthly_sars_pct_pos) ///
    b(%9.4f) se stats(N ll chi2)

display _newline
display "********************************************************************************"
display "ANALYSIS COMPLETE - VERSION 3 (RAW VALUES, NO NORMALIZATION)"
display "********************************************************************************"
display "This version matches the original paper methodology:"
display "  - Raw percent positivity for all viruses"
display "  - Rescaled adjusted cases (divided by 10,000)"
display "  - PHEIC model includes BOTH SARS measures"
display "  - SARS PAF sets BOTH measures to zero"
display _newline
display "Expected results:"
display "  - Pre-pandemic Flu A PAF should be ~1.8% (like original paper)"
display "  - PHEIC SARS-CoV-2 PAF should be substantial"
display "  - Post-PHEIC SARS-CoV-2 PAF = ??? (the key question!)"
display "********************************************************************************"

********************************************************************************
* 7. COMBINED PANDEMIC PERIOD ANALYSIS (March 2020 - Feb 2025)
********************************************************************************
* Non-stratified analysis across entire SARS-CoV-2 era for comparison

display _newline(2)
display "********************************************************************************"
display "COMBINED PANDEMIC PERIOD: March 2020 - February 2025"
display "********************************************************************************"
display "Non-stratified analysis combining PHEIC and post-PHEIC periods"
display _newline

* Define combined pandemic period
gen pandemic = (date >= date("March 1, 2020", "MDY"))

* Long-term trend for combined pandemic period
gen linear_pan = series_month - 390 if pandemic == 1
gen mosq_pan = linear_pan^2 if pandemic == 1

*--- Model WITH Fourier seasonal adjustment ---
display _newline
display "MODEL WITH FOURIER SEASONAL ADJUSTMENT"
display "----------------------------------------"

nbreg deaths sinmo cosmo linear_pan mosq_pan ///
    monthly_flua_pct_pos monthly_flub_pct_pos monthly_rsv_pct_pos ///
    adj_cases_10k monthly_sars_pct_pos ///
    if pandemic == 1, ///
    exposure(population) irr

estimates store pandemic_fft

* Population Attributable Fractions
display _newline "Influenza A PAF:"
punaf if pandemic == 1, atspec(monthly_flua_pct_pos=0) eform

display _newline "Influenza B PAF:"
punaf if pandemic == 1, atspec(monthly_flub_pct_pos=0) eform

display _newline "RSV PAF:"
punaf if pandemic == 1, atspec(monthly_rsv_pct_pos=0) eform

display _newline "SARS-CoV-2 PAF (setting BOTH measures to zero):"
punaf if pandemic == 1, atspec(adj_cases_10k=0 monthly_sars_pct_pos=0) eform

*--- Model WITHOUT Fourier (sensitivity) ---
display _newline
display "MODEL WITHOUT FOURIER (SENSITIVITY ANALYSIS)"
display "---------------------------------------------"

nbreg deaths linear_pan mosq_pan ///
    monthly_flua_pct_pos monthly_flub_pct_pos monthly_rsv_pct_pos ///
    adj_cases_10k monthly_sars_pct_pos ///
    if pandemic == 1, ///
    exposure(population) irr

estimates store pandemic_nofft

display _newline "Influenza A PAF (no Fourier):"
punaf if pandemic == 1, atspec(monthly_flua_pct_pos=0) eform

display _newline "Influenza B PAF (no Fourier):"
punaf if pandemic == 1, atspec(monthly_flub_pct_pos=0) eform

display _newline "RSV PAF (no Fourier):"
punaf if pandemic == 1, atspec(monthly_rsv_pct_pos=0) eform

display _newline "SARS-CoV-2 PAF (no Fourier, setting BOTH measures to zero):"
punaf if pandemic == 1, atspec(adj_cases_10k=0 monthly_sars_pct_pos=0) eform

display _newline(2)
display "********************************************************************************"
display "COMBINED PANDEMIC PERIOD SUMMARY"
display "********************************************************************************"
display "This analysis pools PHEIC and post-PHEIC for overall COVID PAF estimate"
display "N = 60 months (March 2020 - Feb 2025)"
display "Compare to stratified results to assess temporal changes"
display "********************************************************************************"

********************************************************************************
* 8. SENSITIVITY ANALYSIS: COMBINED PANDEMIC EXCLUDING MAR-JUN 2020
********************************************************************************
* Excludes first 4 months of pandemic (March-June 2020) to assess robustness
* of combined PAF estimate to early chaotic period (healthcare disruption,
* coding inconsistencies, testing capacity constraints)
* Suggested by reviewer pre-emption strategy (Alex, Feb 2026)

display _newline(2)
display "********************************************************************************"
display "SENSITIVITY: COMBINED PANDEMIC EXCLUDING MARCH-JUNE 2020"
display "********************************************************************************"
display "Rationale: Early pandemic months characterized by:"
display "  - Disrupted healthcare utilization and deferred care"
display "  - Limited testing capacity and inconsistent cause-of-death coding"
display "  - Extreme mortality volatility unrepresentative of endemic dynamics"
display "Excluding Mar-Jun 2020 (n=4 months); analysis period: Jul 2020 - Feb 2025"
display _newline

* Define restricted pandemic period (July 2020 onward)
gen pandemic_restricted = (date >= date("July 1, 2020", "MDY"))

* Long-term trend for restricted period
* Center on midpoint of Jul 2020 - Feb 2025 (approx series_month 354 + 30 = ~month 396)
gen linear_panr = series_month - 396 if pandemic_restricted == 1
gen mosq_panr = linear_panr^2 if pandemic_restricted == 1

*--- Model WITH Fourier seasonal adjustment ---
display _newline
display "MODEL WITH FOURIER SEASONAL ADJUSTMENT (Jul 2020 - Feb 2025, n=56)"
display "------------------------------------------------------------------------"

nbreg deaths sinmo cosmo linear_panr mosq_panr ///
    monthly_flua_pct_pos monthly_flub_pct_pos monthly_rsv_pct_pos ///
    adj_cases_10k monthly_sars_pct_pos ///
    if pandemic_restricted == 1, ///
    exposure(population) irr

estimates store panr_fft

display _newline "Influenza A PAF:"
punaf if pandemic_restricted == 1, atspec(monthly_flua_pct_pos=0) eform

display _newline "Influenza B PAF:"
punaf if pandemic_restricted == 1, atspec(monthly_flub_pct_pos=0) eform

display _newline "RSV PAF:"
punaf if pandemic_restricted == 1, atspec(monthly_rsv_pct_pos=0) eform

display _newline "SARS-CoV-2 PAF (setting BOTH measures to zero):"
punaf if pandemic_restricted == 1, atspec(adj_cases_10k=0 monthly_sars_pct_pos=0) eform

*--- Model WITHOUT Fourier (sensitivity of sensitivity) ---
display _newline
display "MODEL WITHOUT FOURIER (Jul 2020 - Feb 2025)"
display "---------------------------------------------"

nbreg deaths linear_panr mosq_panr ///
    monthly_flua_pct_pos monthly_flub_pct_pos monthly_rsv_pct_pos ///
    adj_cases_10k monthly_sars_pct_pos ///
    if pandemic_restricted == 1, ///
    exposure(population) irr

estimates store panr_nofft

display _newline "Influenza A PAF (no Fourier):"
punaf if pandemic_restricted == 1, atspec(monthly_flua_pct_pos=0) eform

display _newline "Influenza B PAF (no Fourier):"
punaf if pandemic_restricted == 1, atspec(monthly_flub_pct_pos=0) eform

display _newline "RSV PAF (no Fourier):"
punaf if pandemic_restricted == 1, atspec(monthly_rsv_pct_pos=0) eform

display _newline "SARS-CoV-2 PAF (no Fourier, setting BOTH measures to zero):"
punaf if pandemic_restricted == 1, atspec(adj_cases_10k=0 monthly_sars_pct_pos=0) eform

*--- Direct comparison: full vs restricted pandemic period ---
display _newline(2)
display "********************************************************************************"
display "COMPARISON: FULL vs RESTRICTED COMBINED PANDEMIC PERIOD"
display "********************************************************************************"
display "Full pandemic (Mar 2020-Feb 2025, n=60) vs"
display "Restricted pandemic (Jul 2020-Feb 2025, n=56)"
display _newline
display "Coefficient comparison (Fourier models):"

estimates table pandemic_fft panr_fft, ///
    keep(monthly_flua_pct_pos monthly_flub_pct_pos monthly_rsv_pct_pos ///
         adj_cases_10k monthly_sars_pct_pos) ///
    b(%9.4f) se stats(N ll chi2) ///
    title("Full pandemic vs Restricted pandemic (excl. Mar-Jun 2020)")

display _newline
display "Interpretation guide:"
display "  - If SARS-CoV-2 PAF similar in both: Result robust to early chaos"
display "  - If PAF drops but remains significant: Conservative estimate available"
display "  - If PAF drops substantially: Early period was driving estimate (flag)"
display "********************************************************************************"

********************************************************************************
* SECTION 9: GENERATE COUNTERFACTUAL PREDICTIONS FOR FIGURE 2
* Uses manual exp(_b[...]) approach to compute attributable death rates
* per 100,000 per month for each virus and period
* Output saved as figure2_data.dta and figure2_data.csv
********************************************************************************

display _newline(2)
display "********************************************************************************"
display "SECTION 9: GENERATING COUNTERFACTUAL PREDICTIONS FOR FIGURE 2"
display "********************************************************************************"
capture mkdir data

* Population denominator: person-months per 100,000
gen pop100k = population / 12 / 100000

* Initialize attributable rate variables
gen flua_attr_fft     = .
gen flub_attr_fft     = .
gen rsv_attr_fft      = .
gen rsv_attr_nofft    = .
gen sars_attr_fft     = .
gen combined_attr_fft = .

********************************************************************************
* PRE-PANDEMIC (Jan 1993 - Feb 2020)
********************************************************************************

estimates restore prepan_fft

* Fitted values (observed viral levels)
gen pred_pre = exp( ///
    _b[_cons] + ///
    _b[sinmo]*sinmo + _b[cosmo]*cosmo + ///
    _b[linear_pre]*linear_pre + _b[mosq_pre]*mosq_pre + ///
    _b[monthly_flua_pct_pos]*monthly_flua_pct_pos + ///
    _b[monthly_flub_pct_pos]*monthly_flub_pct_pos + ///
    _b[monthly_rsv_pct_pos]*monthly_rsv_pct_pos + ///
    ln(population)) if prepandemic == 1

* Flu A zeroed
gen cf_pre_fluA = exp( ///
    _b[_cons] + ///
    _b[sinmo]*sinmo + _b[cosmo]*cosmo + ///
    _b[linear_pre]*linear_pre + _b[mosq_pre]*mosq_pre + ///
    _b[monthly_flua_pct_pos]*0 + ///
    _b[monthly_flub_pct_pos]*monthly_flub_pct_pos + ///
    _b[monthly_rsv_pct_pos]*monthly_rsv_pct_pos + ///
    ln(population)) if prepandemic == 1

* Flu B zeroed
gen cf_pre_fluB = exp( ///
    _b[_cons] + ///
    _b[sinmo]*sinmo + _b[cosmo]*cosmo + ///
    _b[linear_pre]*linear_pre + _b[mosq_pre]*mosq_pre + ///
    _b[monthly_flua_pct_pos]*monthly_flua_pct_pos + ///
    _b[monthly_flub_pct_pos]*0 + ///
    _b[monthly_rsv_pct_pos]*monthly_rsv_pct_pos + ///
    ln(population)) if prepandemic == 1

* RSV zeroed (with FFT)
gen cf_pre_rsv = exp( ///
    _b[_cons] + ///
    _b[sinmo]*sinmo + _b[cosmo]*cosmo + ///
    _b[linear_pre]*linear_pre + _b[mosq_pre]*mosq_pre + ///
    _b[monthly_flua_pct_pos]*monthly_flua_pct_pos + ///
    _b[monthly_flub_pct_pos]*monthly_flub_pct_pos + ///
    _b[monthly_rsv_pct_pos]*0 + ///
    ln(population)) if prepandemic == 1

* All viruses zeroed
gen cf_pre_all = exp( ///
    _b[_cons] + ///
    _b[sinmo]*sinmo + _b[cosmo]*cosmo + ///
    _b[linear_pre]*linear_pre + _b[mosq_pre]*mosq_pre + ///
    _b[monthly_flua_pct_pos]*0 + ///
    _b[monthly_flub_pct_pos]*0 + ///
    _b[monthly_rsv_pct_pos]*0 + ///
    ln(population)) if prepandemic == 1

* Attributable rates (with FFT)
replace flua_attr_fft     = (pred_pre - cf_pre_fluA) / pop100k if prepandemic == 1
replace flub_attr_fft     = (pred_pre - cf_pre_fluB) / pop100k if prepandemic == 1
replace rsv_attr_fft      = (pred_pre - cf_pre_rsv)  / pop100k if prepandemic == 1
replace combined_attr_fft = (pred_pre - cf_pre_all)  / pop100k if prepandemic == 1

* RSV without FFT
estimates restore prepan_nofft

gen pred_pre_nofft = exp( ///
    _b[_cons] + ///
    _b[linear_pre]*linear_pre + _b[mosq_pre]*mosq_pre + ///
    _b[monthly_flua_pct_pos]*monthly_flua_pct_pos + ///
    _b[monthly_flub_pct_pos]*monthly_flub_pct_pos + ///
    _b[monthly_rsv_pct_pos]*monthly_rsv_pct_pos + ///
    ln(population)) if prepandemic == 1

gen cf_pre_rsv_nofft = exp( ///
    _b[_cons] + ///
    _b[linear_pre]*linear_pre + _b[mosq_pre]*mosq_pre + ///
    _b[monthly_flua_pct_pos]*monthly_flua_pct_pos + ///
    _b[monthly_flub_pct_pos]*monthly_flub_pct_pos + ///
    _b[monthly_rsv_pct_pos]*0 + ///
    ln(population)) if prepandemic == 1

replace rsv_attr_nofft = (pred_pre_nofft - cf_pre_rsv_nofft) / pop100k if prepandemic == 1

drop pred_pre cf_pre_fluA cf_pre_fluB cf_pre_rsv cf_pre_all pred_pre_nofft cf_pre_rsv_nofft

********************************************************************************
* PHEIC PERIOD (March 2020 - April 2023)
********************************************************************************

estimates restore pheic_fft

gen pred_pheic = exp( ///
    _b[_cons] + ///
    _b[sinmo]*sinmo + _b[cosmo]*cosmo + ///
    _b[linear_pheic]*linear_pheic + _b[mosq_pheic]*mosq_pheic + ///
    _b[monthly_flua_pct_pos]*monthly_flua_pct_pos + ///
    _b[monthly_flub_pct_pos]*monthly_flub_pct_pos + ///
    _b[monthly_rsv_pct_pos]*monthly_rsv_pct_pos + ///
    _b[adj_cases_10k]*adj_cases_10k + ///
    _b[monthly_sars_pct_pos]*monthly_sars_pct_pos + ///
    ln(population)) if pheic == 1

* Flu A zeroed
gen cf_pheic_fluA = exp( ///
    _b[_cons] + ///
    _b[sinmo]*sinmo + _b[cosmo]*cosmo + ///
    _b[linear_pheic]*linear_pheic + _b[mosq_pheic]*mosq_pheic + ///
    _b[monthly_flua_pct_pos]*0 + ///
    _b[monthly_flub_pct_pos]*monthly_flub_pct_pos + ///
    _b[monthly_rsv_pct_pos]*monthly_rsv_pct_pos + ///
    _b[adj_cases_10k]*adj_cases_10k + ///
    _b[monthly_sars_pct_pos]*monthly_sars_pct_pos + ///
    ln(population)) if pheic == 1

* Flu B zeroed
gen cf_pheic_fluB = exp( ///
    _b[_cons] + ///
    _b[sinmo]*sinmo + _b[cosmo]*cosmo + ///
    _b[linear_pheic]*linear_pheic + _b[mosq_pheic]*mosq_pheic + ///
    _b[monthly_flua_pct_pos]*monthly_flua_pct_pos + ///
    _b[monthly_flub_pct_pos]*0 + ///
    _b[monthly_rsv_pct_pos]*monthly_rsv_pct_pos + ///
    _b[adj_cases_10k]*adj_cases_10k + ///
    _b[monthly_sars_pct_pos]*monthly_sars_pct_pos + ///
    ln(population)) if pheic == 1

* RSV zeroed (with FFT)
gen cf_pheic_rsv = exp( ///
    _b[_cons] + ///
    _b[sinmo]*sinmo + _b[cosmo]*cosmo + ///
    _b[linear_pheic]*linear_pheic + _b[mosq_pheic]*mosq_pheic + ///
    _b[monthly_flua_pct_pos]*monthly_flua_pct_pos + ///
    _b[monthly_flub_pct_pos]*monthly_flub_pct_pos + ///
    _b[monthly_rsv_pct_pos]*0 + ///
    _b[adj_cases_10k]*adj_cases_10k + ///
    _b[monthly_sars_pct_pos]*monthly_sars_pct_pos + ///
    ln(population)) if pheic == 1

* SARS-CoV-2 zeroed (both measures)
gen cf_pheic_sars = exp( ///
    _b[_cons] + ///
    _b[sinmo]*sinmo + _b[cosmo]*cosmo + ///
    _b[linear_pheic]*linear_pheic + _b[mosq_pheic]*mosq_pheic + ///
    _b[monthly_flua_pct_pos]*monthly_flua_pct_pos + ///
    _b[monthly_flub_pct_pos]*monthly_flub_pct_pos + ///
    _b[monthly_rsv_pct_pos]*monthly_rsv_pct_pos + ///
    _b[adj_cases_10k]*0 + ///
    _b[monthly_sars_pct_pos]*0 + ///
    ln(population)) if pheic == 1

* All viruses zeroed
gen cf_pheic_all = exp( ///
    _b[_cons] + ///
    _b[sinmo]*sinmo + _b[cosmo]*cosmo + ///
    _b[linear_pheic]*linear_pheic + _b[mosq_pheic]*mosq_pheic + ///
    _b[monthly_flua_pct_pos]*0 + ///
    _b[monthly_flub_pct_pos]*0 + ///
    _b[monthly_rsv_pct_pos]*0 + ///
    _b[adj_cases_10k]*0 + ///
    _b[monthly_sars_pct_pos]*0 + ///
    ln(population)) if pheic == 1

replace flua_attr_fft     = (pred_pheic - cf_pheic_fluA) / pop100k if pheic == 1
replace flub_attr_fft     = (pred_pheic - cf_pheic_fluB) / pop100k if pheic == 1
replace rsv_attr_fft      = (pred_pheic - cf_pheic_rsv)  / pop100k if pheic == 1
replace sars_attr_fft     = (pred_pheic - cf_pheic_sars) / pop100k if pheic == 1
replace combined_attr_fft = (pred_pheic - cf_pheic_all)  / pop100k if pheic == 1

* RSV without FFT
estimates restore pheic_nofft

gen pred_pheic_nofft = exp( ///
    _b[_cons] + ///
    _b[linear_pheic]*linear_pheic + _b[mosq_pheic]*mosq_pheic + ///
    _b[monthly_flua_pct_pos]*monthly_flua_pct_pos + ///
    _b[monthly_flub_pct_pos]*monthly_flub_pct_pos + ///
    _b[monthly_rsv_pct_pos]*monthly_rsv_pct_pos + ///
    _b[adj_cases_10k]*adj_cases_10k + ///
    _b[monthly_sars_pct_pos]*monthly_sars_pct_pos + ///
    ln(population)) if pheic == 1

gen cf_pheic_rsv_nofft = exp( ///
    _b[_cons] + ///
    _b[linear_pheic]*linear_pheic + _b[mosq_pheic]*mosq_pheic + ///
    _b[monthly_flua_pct_pos]*monthly_flua_pct_pos + ///
    _b[monthly_flub_pct_pos]*monthly_flub_pct_pos + ///
    _b[monthly_rsv_pct_pos]*0 + ///
    _b[adj_cases_10k]*adj_cases_10k + ///
    _b[monthly_sars_pct_pos]*monthly_sars_pct_pos + ///
    ln(population)) if pheic == 1

replace rsv_attr_nofft = (pred_pheic_nofft - cf_pheic_rsv_nofft) / pop100k if pheic == 1

drop pred_pheic cf_pheic_fluA cf_pheic_fluB cf_pheic_rsv cf_pheic_sars cf_pheic_all
drop pred_pheic_nofft cf_pheic_rsv_nofft

********************************************************************************
* POST-PHEIC PERIOD (May 2023 - Feb 2025)
********************************************************************************

estimates restore postpheic_fft

gen pred_post = exp( ///
    _b[_cons] + ///
    _b[sinmo]*sinmo + _b[cosmo]*cosmo + ///
    _b[linear_post]*linear_post + _b[mosq_post]*mosq_post + ///
    _b[monthly_flua_pct_pos]*monthly_flua_pct_pos + ///
    _b[monthly_flub_pct_pos]*monthly_flub_pct_pos + ///
    _b[monthly_rsv_pct_pos]*monthly_rsv_pct_pos + ///
    _b[monthly_sars_pct_pos]*monthly_sars_pct_pos + ///
    ln(population)) if postpheic == 1

* Flu A zeroed
gen cf_post_fluA = exp( ///
    _b[_cons] + ///
    _b[sinmo]*sinmo + _b[cosmo]*cosmo + ///
    _b[linear_post]*linear_post + _b[mosq_post]*mosq_post + ///
    _b[monthly_flua_pct_pos]*0 + ///
    _b[monthly_flub_pct_pos]*monthly_flub_pct_pos + ///
    _b[monthly_rsv_pct_pos]*monthly_rsv_pct_pos + ///
    _b[monthly_sars_pct_pos]*monthly_sars_pct_pos + ///
    ln(population)) if postpheic == 1

* Flu B zeroed
gen cf_post_fluB = exp( ///
    _b[_cons] + ///
    _b[sinmo]*sinmo + _b[cosmo]*cosmo + ///
    _b[linear_post]*linear_post + _b[mosq_post]*mosq_post + ///
    _b[monthly_flua_pct_pos]*monthly_flua_pct_pos + ///
    _b[monthly_flub_pct_pos]*0 + ///
    _b[monthly_rsv_pct_pos]*monthly_rsv_pct_pos + ///
    _b[monthly_sars_pct_pos]*monthly_sars_pct_pos + ///
    ln(population)) if postpheic == 1

* RSV zeroed (with FFT)
gen cf_post_rsv = exp( ///
    _b[_cons] + ///
    _b[sinmo]*sinmo + _b[cosmo]*cosmo + ///
    _b[linear_post]*linear_post + _b[mosq_post]*mosq_post + ///
    _b[monthly_flua_pct_pos]*monthly_flua_pct_pos + ///
    _b[monthly_flub_pct_pos]*monthly_flub_pct_pos + ///
    _b[monthly_rsv_pct_pos]*0 + ///
    _b[monthly_sars_pct_pos]*monthly_sars_pct_pos + ///
    ln(population)) if postpheic == 1

* SARS-CoV-2 zeroed (% pos only in post-PHEIC)
gen cf_post_sars = exp( ///
    _b[_cons] + ///
    _b[sinmo]*sinmo + _b[cosmo]*cosmo + ///
    _b[linear_post]*linear_post + _b[mosq_post]*mosq_post + ///
    _b[monthly_flua_pct_pos]*monthly_flua_pct_pos + ///
    _b[monthly_flub_pct_pos]*monthly_flub_pct_pos + ///
    _b[monthly_rsv_pct_pos]*monthly_rsv_pct_pos + ///
    _b[monthly_sars_pct_pos]*0 + ///
    ln(population)) if postpheic == 1

* All viruses zeroed
gen cf_post_all = exp( ///
    _b[_cons] + ///
    _b[sinmo]*sinmo + _b[cosmo]*cosmo + ///
    _b[linear_post]*linear_post + _b[mosq_post]*mosq_post + ///
    _b[monthly_flua_pct_pos]*0 + ///
    _b[monthly_flub_pct_pos]*0 + ///
    _b[monthly_rsv_pct_pos]*0 + ///
    _b[monthly_sars_pct_pos]*0 + ///
    ln(population)) if postpheic == 1

replace flua_attr_fft     = (pred_post - cf_post_fluA) / pop100k if postpheic == 1
replace flub_attr_fft     = (pred_post - cf_post_fluB) / pop100k if postpheic == 1
replace rsv_attr_fft      = (pred_post - cf_post_rsv)  / pop100k if postpheic == 1
replace sars_attr_fft     = (pred_post - cf_post_sars) / pop100k if postpheic == 1
replace combined_attr_fft = (pred_post - cf_post_all)  / pop100k if postpheic == 1

* RSV without FFT
estimates restore postpheic_nofft

gen pred_post_nofft = exp( ///
    _b[_cons] + ///
    _b[linear_post]*linear_post + _b[mosq_post]*mosq_post + ///
    _b[monthly_flua_pct_pos]*monthly_flua_pct_pos + ///
    _b[monthly_flub_pct_pos]*monthly_flub_pct_pos + ///
    _b[monthly_rsv_pct_pos]*monthly_rsv_pct_pos + ///
    _b[monthly_sars_pct_pos]*monthly_sars_pct_pos + ///
    ln(population)) if postpheic == 1

gen cf_post_rsv_nofft = exp( ///
    _b[_cons] + ///
    _b[linear_post]*linear_post + _b[mosq_post]*mosq_post + ///
    _b[monthly_flua_pct_pos]*monthly_flua_pct_pos + ///
    _b[monthly_flub_pct_pos]*monthly_flub_pct_pos + ///
    _b[monthly_rsv_pct_pos]*0 + ///
    _b[monthly_sars_pct_pos]*monthly_sars_pct_pos + ///
    ln(population)) if postpheic == 1

replace rsv_attr_nofft = (pred_post_nofft - cf_post_rsv_nofft) / pop100k if postpheic == 1

drop pred_post cf_post_fluA cf_post_fluB cf_post_rsv cf_post_sars cf_post_all
drop pred_post_nofft cf_post_rsv_nofft

********************************************************************************
* SAVE OUTPUT FOR R FIGURES
********************************************************************************

display _newline "Saving figure2_data..."

* Keep only variables needed for Figure 2
keep date year month series_month population ///
     prepandemic pheic postpheic pandemic ///
     flua_attr_fft flub_attr_fft rsv_attr_fft rsv_attr_nofft ///
     sars_attr_fft combined_attr_fft

* Format date as string for CSV export
gen date_str = string(year) + "-" + string(month, "%02.0f") + "-15"

save "data/figure2_data.dta", replace
export delimited "data/figure2_data.csv", replace

display "Figure 2 data saved to data/figure2_data.dta and data/figure2_data.csv"
display "Variables: flua_attr_fft, rsv_attr_fft, rsv_attr_nofft, sars_attr_fft, combined_attr_fft"
display "All rates are deaths per 100,000 per month"
