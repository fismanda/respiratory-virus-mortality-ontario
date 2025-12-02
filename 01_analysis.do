********************************************************************************
* RESPIRATORY VIRUS MORTALITY ANALYSIS - ONTARIO 1991-2024
* Population Attributable Risk Estimation Using Ecological Regression
********************************************************************************
* This code estimates the fraction of deaths in Ontario attributable to 
* influenza A, influenza B, RSV, and SARS-CoV-2 using negative binomial 
* regression with Fourier seasonal adjustment
********************************************************************************

clear all
set more off

********************************************************************************
* 1. DATA PREPARATION
********************************************************************************

* Load monthly aggregated data
import delimited "data/monthly_mortality_virology_data_ontario_1991_2024.csv", clear

* Create analysis variables
gen series_month = month + 12*(year-1993)

* Fourier transform terms for seasonal adjustment
gen sinmo = sin(6.28*month/12)
gen cosmo = cos(6.28*month/12)

* Long-term trend terms
gen linear = series_month - (378/2)
gen mosq = linear^2

* Pandemic indicator (March 2020 onward)
gen pandemic = (date >= date("March 1, 2020", "MDY"))

* Rescale adjusted cases for interpretable coefficients
gen adj_cases_10k = monthly_adj_cases/10000

* Replace missing SARS-CoV-2 data with zeros for pre-pandemic period
replace monthly_sars_pct_pos = 0 if monthly_sars_pct_pos == .
replace adj_cases_10k = 0 if adj_cases_10k == .

* Restrict to series_month < 376 due to reporting delays
* (data through March 2024)
keep if series_month < 376

********************************************************************************
* 2. PRE-PANDEMIC ANALYSIS (before March 2020)
********************************************************************************

*-------------------------------------------------------------------------------
* 2a. Primary model WITH Fourier seasonal adjustment
*-------------------------------------------------------------------------------

nbreg deaths sinmo cosmo linear mosq ///
    monthly_flua_pct_pos monthly_flub_pct_pos monthly_rsv_pct_pos ///
    if pandemic == 0, ///
    exposure(population) irr

estimates store prepan_fft

* Model diagnostics
estat ic

* Predicted deaths from full model
predict pre_pandemic_deaths if pandemic == 0, n
gen pre_pandemic_death_rate = 100000*pre_pandemic_deaths/(population/12) if pandemic == 0

*--- Influenza A ---
* Population attributable fraction using PUNAF
punaf if pandemic == 0, atspec(monthly_flua_pct_pos=0) eform

* Counterfactual deaths without influenza A
gen pred_death_no_flua_pre_pan = exp(_b[_cons] + _b[linear]*linear + ///
    _b[mosq]*mosq + _b[sinmo]*sinmo + _b[cosmo]*cosmo + ///
    _b[monthly_flub_pct_pos]*monthly_flub_pct_pos + ///
    _b[monthly_rsv_pct_pos]*monthly_rsv_pct_pos + ///
    ln(population)) if pandemic == 0

gen pre_pandemic_death_rate_no_flua = 100000*pred_death_no_flua_pre_pan/(population/12) if pandemic == 0

* Attributable deaths and PAR
gen flua_attr_death_pre = pre_pandemic_deaths - pred_death_no_flua_pre_pan
gen flua_attr_death_rate_pre = 100000*flua_attr_death_pre/(population/12)
egen total_flua_deaths_pre = sum(flua_attr_death_pre)
egen total_deaths_pre = sum(pre_pandemic_deaths)
gen par_flua_pre = 1 - (total_deaths_pre - total_flua_deaths_pre)/total_deaths_pre

*--- Influenza B ---
punaf if pandemic == 0, atspec(monthly_flub_pct_pos=0) eform

gen pred_death_no_flub_pre_pan = exp(_b[_cons] + _b[linear]*linear + ///
    _b[mosq]*mosq + _b[sinmo]*sinmo + _b[cosmo]*cosmo + ///
    _b[monthly_flua_pct_pos]*monthly_flua_pct_pos + ///
    _b[monthly_rsv_pct_pos]*monthly_rsv_pct_pos + ///
    ln(population)) if pandemic == 0

gen pre_pandemic_death_rate_no_flub = 100000*pred_death_no_flub_pre_pan/(population/12) if pandemic == 0
gen flub_attr_death_pre = pre_pandemic_deaths - pred_death_no_flub_pre_pan
gen flub_attr_death_rate_pre = 100000*flub_attr_death_pre/(population/12)
egen total_flub_deaths_pre = sum(flub_attr_death_pre)
gen par_flub_pre = 1 - (total_deaths_pre - total_flub_deaths_pre)/total_deaths_pre

*--- RSV ---
punaf if pandemic == 0, atspec(monthly_rsv_pct_pos=0) eform

gen pred_death_no_rsv_pre_pan = exp(_b[_cons] + _b[linear]*linear + ///
    _b[mosq]*mosq + _b[sinmo]*sinmo + _b[cosmo]*cosmo + ///
    _b[monthly_flua_pct_pos]*monthly_flua_pct_pos + ///
    _b[monthly_flub_pct_pos]*monthly_flub_pct_pos + ///
    ln(population)) if pandemic == 0

gen pre_pandemic_death_rate_no_rsv = 100000*pred_death_no_rsv_pre_pan/(population/12) if pandemic == 0
gen rsv_attr_death_pre = pre_pandemic_deaths - pred_death_no_rsv_pre_pan
gen rsv_attr_death_rate_pre = 100000*rsv_attr_death_pre/(population/12)
egen total_rsv_deaths_pre = sum(rsv_attr_death_pre)
gen par_rsv_pre = 1 - (total_deaths_pre - total_rsv_deaths_pre)/total_deaths_pre

*--- All respiratory viruses ---
punaf if pandemic == 0, atspec(monthly_rsv_pct_pos=0 ///
    monthly_flub_pct_pos=0 monthly_flua_pct_pos=0) eform

gen pred_death_no_virus_pre_pan = exp(_b[_cons] + _b[linear]*linear + ///
    _b[mosq]*mosq + _b[sinmo]*sinmo + _b[cosmo]*cosmo + ///
    ln(population)) if pandemic == 0

gen pre_pandemic_death_rate_no_virus = 100000*pred_death_no_virus_pre_pan/(population/12) if pandemic == 0
gen virus_attr_death_pre = pre_pandemic_deaths - pred_death_no_virus_pre_pan
gen virus_attr_death_rate_pre = 100000*virus_attr_death_pre/(population/12)
egen total_virus_deaths_pre = sum(virus_attr_death_pre)
gen par_virus_pre = 1 - (total_deaths_pre - total_virus_deaths_pre)/total_deaths_pre

*-------------------------------------------------------------------------------
* 2b. Sensitivity analysis WITHOUT Fourier seasonal adjustment
*-------------------------------------------------------------------------------

nbreg deaths linear mosq ///
    monthly_flua_pct_pos monthly_flub_pct_pos monthly_rsv_pct_pos ///
    if pandemic == 0, ///
    exposure(population) irr

estimates store prepan_nofft

estat ic

* Predicted deaths from model without FFT
predict pre_pandemic_deaths_nofft if pandemic == 0, n
gen pre_pandemic_death_rate_nofft = 100000*pre_pandemic_deaths_nofft/(population/12) if pandemic == 0

*--- Influenza A (without FFT) ---
punaf if pandemic == 0, atspec(monthly_flua_pct_pos=0) eform

gen pred_death_noflua_prepan_nofft = exp(_b[_cons] + _b[linear]*linear + ///
    _b[mosq]*mosq + ///
    _b[monthly_flub_pct_pos]*monthly_flub_pct_pos + ///
    _b[monthly_rsv_pct_pos]*monthly_rsv_pct_pos + ///
    ln(population)) if pandemic == 0

gen prepan_deathrate_noflua_nofft = 100000*pred_death_noflua_prepan_nofft/(population/12) if pandemic == 0
gen flua_attr_pre_nofft = pre_pandemic_deaths_nofft - pred_death_noflua_prepan_nofft
gen flua_attr_rate_pre_nofft = 100000*flua_attr_pre_nofft/(population/12)
egen total_flua_deaths_pre_nofft = sum(flua_attr_pre_nofft)
egen total_deaths_pre_nofft = sum(pre_pandemic_deaths_nofft)
gen par_flua_pre_nofft = 1 - (total_deaths_pre_nofft - total_flua_deaths_pre_nofft)/total_deaths_pre_nofft

*--- RSV (without FFT) ---
punaf if pandemic == 0, atspec(monthly_rsv_pct_pos=0) eform

gen pred_death_norsv_prepan_nofft = exp(_b[_cons] + _b[linear]*linear + ///
    _b[mosq]*mosq + ///
    _b[monthly_flub_pct_pos]*monthly_flub_pct_pos + ///
    _b[monthly_flua_pct_pos]*monthly_flua_pct_pos + ///
    ln(population)) if pandemic == 0

gen prepan_deathrate_norsv_nofft = 100000*pred_death_norsv_prepan_nofft/(population/12) if pandemic == 0
gen rsv_attr_pre_nofft = pre_pandemic_deaths_nofft - pred_death_norsv_prepan_nofft
gen rsv_attr_rate_pre_nofft = 100000*rsv_attr_pre_nofft/(population/12)
egen total_rsv_deaths_pre_nofft = sum(rsv_attr_pre_nofft)
gen par_rsv_pre_nofft = 1 - (total_deaths_pre_nofft - total_rsv_deaths_pre_nofft)/total_deaths_pre_nofft

*--- Influenza B (without FFT) ---
punaf if pandemic == 0, atspec(monthly_flub_pct_pos=0) eform

gen pred_death_noflub_prepan_nofft = exp(_b[_cons] + _b[linear]*linear + ///
    _b[mosq]*mosq + ///
    _b[monthly_rsv_pct_pos]*monthly_rsv_pct_pos + ///
    _b[monthly_flua_pct_pos]*monthly_flua_pct_pos + ///
    ln(population)) if pandemic == 0

gen pre_pan_deathrate_noflub_nofft = 100000*pred_death_noflub_prepan_nofft/(population/12) if pandemic == 0
gen flub_attr_pre_nofft = pre_pandemic_deaths_nofft - pred_death_noflub_prepan_nofft
gen flub_attr_rate_pre_nofft = 100000*flub_attr_pre_nofft/(population/12)
egen total_flub_deaths_pre_nofft = sum(flub_attr_pre_nofft)
gen par_flub_pre_nofft = 1 - (total_deaths_pre_nofft - total_flub_deaths_pre_nofft)/total_deaths_pre_nofft

*--- All respiratory viruses (without FFT) ---
punaf if pandemic == 0, atspec(monthly_rsv_pct_pos=0 ///
    monthly_flub_pct_pos=0 monthly_flua_pct_pos=0) eform

gen pred_death_no_virus_pre_nofft = exp(_b[_cons] + _b[linear]*linear + ///
    _b[mosq]*mosq + ln(population)) if pandemic == 0

gen pre_death_rate_no_virus_nofft = 100000*pred_death_no_virus_pre_nofft/(population/12) if pandemic == 0
gen virus_attr_death_pre_nofft = pre_pandemic_deaths_nofft - pred_death_no_virus_pre_nofft
gen virus_attr_death_rate_pre_nofft = 100000*virus_attr_death_pre_nofft/(population/12)
egen total_virus_deaths_pre_nofft = sum(virus_attr_death_pre_nofft)
gen par_virus_pre_nofft = 1 - (total_deaths_pre_nofft - total_virus_deaths_pre_nofft)/total_deaths_pre_nofft

********************************************************************************
* 3. PANDEMIC ANALYSIS (March 2020 onward)
********************************************************************************

* First confirm pandemic increased mortality after adjusting for trends/seasonality
nbreg deaths sinmo cosmo linear mosq pandemic, exposure(population) irr
* Expected result: IRR ~1.050 (95% CI 1.022-1.079), P<0.001

*-------------------------------------------------------------------------------
* 3a. Primary model WITH Fourier seasonal adjustment
*-------------------------------------------------------------------------------

nbreg deaths sinmo cosmo linear mosq ///
    monthly_flua_pct_pos monthly_flub_pct_pos monthly_rsv_pct_pos ///
    adj_cases_10k monthly_sars_pct_pos ///
    if pandemic == 1, ///
    exposure(population) irr

estimates store pan_fft

estat ic

* Predicted deaths from full model
predict pandemic_deaths if pandemic == 1, n
gen pandemic_death_rate = 100000*pandemic_deaths/(population/12) if pandemic == 1

*--- SARS-CoV-2 ---
* Note: adj_cases_10k covers March 2020-August 2022 (Science Table)
*       monthly_sars_pct_pos covers September 2022-March 2024 (FluWatch)

punaf if pandemic == 1, atspec(adj_cases_10k=0 monthly_sars_pct_pos=0) eform

gen pred_death_no_sars_pan = exp(_b[_cons] + _b[linear]*linear + ///
    _b[mosq]*mosq + _b[sinmo]*sinmo + _b[cosmo]*cosmo + ///
    _b[monthly_flua_pct_pos]*monthly_flua_pct_pos + ///
    _b[monthly_flub_pct_pos]*monthly_flub_pct_pos + ///
    _b[monthly_rsv_pct_pos]*monthly_rsv_pct_pos + ///
    ln(population)) if pandemic == 1

gen pandemic_death_rate_no_sars = 100000*pred_death_no_sars_pan/(population/12) if pandemic == 1
gen sars2_attr_death = pandemic_deaths - pred_death_no_sars_pan if pandemic == 1
gen sars2_attr_death_rate = 100000*sars2_attr_death/(population/12) if pandemic == 1
egen total_sars2_deaths = sum(sars2_attr_death) if pandemic == 1
egen total_deaths_pan = sum(pandemic_deaths) if pandemic == 1
gen par_sars2_pan = 1 - (total_deaths_pan - total_sars2_deaths)/total_deaths_pan if pandemic == 1

*--- Influenza B (pandemic) ---
punaf if pandemic == 1, atspec(monthly_flub_pct_pos=0) eform

gen pred_death_no_flub_pan = exp(_b[_cons] + _b[linear]*linear + ///
    _b[mosq]*mosq + _b[sinmo]*sinmo + _b[cosmo]*cosmo + ///
    _b[monthly_flua_pct_pos]*monthly_flua_pct_pos + ///
    _b[monthly_sars_pct_pos]*monthly_sars_pct_pos + ///
    _b[adj_cases_10k]*adj_cases_10k + ///
    _b[monthly_rsv_pct_pos]*monthly_rsv_pct_pos + ///
    ln(population)) if pandemic == 1

gen pandemic_death_rate_no_flub = 100000*pred_death_no_flub_pan/(population/12) if pandemic == 1
gen flub_attr_death_pan = pandemic_deaths - pred_death_no_flub_pan
gen flub_attr_death_rate_pan = 100000*flub_attr_death_pan/(population/12) if pandemic == 1
egen total_flub_deaths_pan = sum(flub_attr_death_pan) if pandemic == 1
gen par_flub_pan = 1 - (total_deaths_pan - total_flub_deaths_pan)/total_deaths_pan

*--- RSV (pandemic) ---
punaf if pandemic == 1, atspec(monthly_rsv_pct_pos=0) eform

gen pred_death_no_rsv_pan = exp(_b[_cons] + _b[linear]*linear + ///
    _b[mosq]*mosq + _b[sinmo]*sinmo + _b[cosmo]*cosmo + ///
    _b[monthly_flua_pct_pos]*monthly_flua_pct_pos + ///
    _b[monthly_sars_pct_pos]*monthly_sars_pct_pos + ///
    _b[adj_cases_10k]*adj_cases_10k + ///
    _b[monthly_flub_pct_pos]*monthly_flub_pct_pos + ///
    ln(population)) if pandemic == 1

gen pandemic_death_rate_no_rsv = 100000*pred_death_no_rsv_pan/(population/12) if pandemic == 1
gen rsv_attr_death_pan = pandemic_deaths - pred_death_no_rsv_pan
gen rsv_attr_death_rate_pan = 100000*rsv_attr_death_pan/(population/12)
egen total_rsv_deaths_pan = sum(rsv_attr_death_pan)
gen par_rsv_pan = 1 - (total_deaths_pan - total_rsv_deaths_pan)/total_deaths_pan

*--- Influenza A (pandemic) ---
punaf if pandemic == 1, atspec(monthly_flua_pct_pos=0) eform

gen pred_death_no_flua_pan = exp(_b[_cons] + _b[linear]*linear + ///
    _b[mosq]*mosq + _b[sinmo]*sinmo + _b[cosmo]*cosmo + ///
    _b[monthly_flub_pct_pos]*monthly_flub_pct_pos + ///
    _b[monthly_sars_pct_pos]*monthly_sars_pct_pos + ///
    _b[adj_cases_10k]*adj_cases_10k + ///
    _b[monthly_rsv_pct_pos]*monthly_rsv_pct_pos + ///
    ln(population)) if pandemic == 1

gen pandemic_death_rate_no_flua = 100000*pred_death_no_flua_pan/(population/12) if pandemic == 1
gen flua_attr_death_pan = pandemic_deaths - pred_death_no_flua_pan
gen flua_attr_death_rate_pan = 100000*flua_attr_death_pan/(population/12)
egen total_flua_deaths_pan = sum(flua_attr_death_pan)
gen par_flua_pan = 1 - (total_deaths_pan - total_flua_deaths_pan)/total_deaths_pan

*--- All respiratory viruses (pandemic) ---
punaf if pandemic == 1, atspec(monthly_flua_pct_pos=0 monthly_flub_pct_pos=0 ///
    adj_cases_10k=0 monthly_sars_pct_pos=0 monthly_rsv_pct_pos=0) eform

gen pred_death_no_virus_pan = exp(_b[_cons] + _b[linear]*linear + ///
    _b[mosq]*mosq + _b[sinmo]*sinmo + _b[cosmo]*cosmo + ///
    ln(population)) if pandemic == 1

gen pandemic_death_rate_no_virus = 100000*pred_death_no_virus_pan/(population/12) if pandemic == 1
gen virus_attr_death_pan = pandemic_deaths - pred_death_no_virus_pan
gen virus_attr_death_rate_pan = 100000*virus_attr_death_pan/(population/12)
egen total_virus_deaths_pan = sum(virus_attr_death_pan)
gen par_virus_pan = 1 - (total_deaths_pan - total_virus_deaths_pan)/total_deaths_pan

*-------------------------------------------------------------------------------
* 3b. Sensitivity analysis WITHOUT Fourier seasonal adjustment
*-------------------------------------------------------------------------------

nbreg deaths linear mosq ///
    monthly_flua_pct_pos monthly_flub_pct_pos monthly_rsv_pct_pos ///
    monthly_sars_pct_pos adj_cases_10k ///
    if pandemic == 1, ///
    exposure(population) irr

estimates store pan_nofft

estat ic

* Predicted deaths from model without FFT
predict pandemic_deaths_nofft if pandemic == 1, n
gen pandemic_death_rate_nofft = 100000*pandemic_deaths_nofft/(population/12) if pandemic == 1

*--- SARS-CoV-2 (without FFT) ---
punaf if pandemic == 1, atspec(monthly_sars_pct_pos=0 adj_cases_10k=0) eform

gen pred_death_nosars2_pan_nofft = exp(_b[_cons] + _b[linear]*linear + ///
    _b[mosq]*mosq + ///
    _b[monthly_flub_pct_pos]*monthly_flub_pct_pos + ///
    _b[monthly_flua_pct_pos]*monthly_flua_pct_pos + ///
    _b[monthly_rsv_pct_pos]*monthly_rsv_pct_pos + ///
    ln(population)) if pandemic == 1

gen pan_deathrate_nosars2_nofft = 100000*pred_death_nosars2_pan_nofft/(population/12) if pandemic == 1
gen sars2_attr_pan_nofft = pandemic_deaths_nofft - pred_death_nosars2_pan_nofft
gen sars2_attr_rate_pan_nofft = 100000*sars2_attr_pan_nofft/(population/12)
egen total_sars2_deaths_pan_nofft = sum(sars2_attr_pan_nofft)
egen total_deaths_pan_nofft = sum(pandemic_deaths_nofft) if pandemic == 1
gen par_sars2_pan_nofft = 1 - (total_deaths_pan_nofft - total_sars2_deaths_pan_nofft)/total_deaths_pan_nofft

*--- RSV (without FFT, pandemic) ---
punaf if pandemic == 1, atspec(monthly_rsv_pct_pos=0) eform

gen pred_death_norsv_pan_nofft = exp(_b[_cons] + _b[linear]*linear + ///
    _b[mosq]*mosq + ///
    _b[monthly_flub_pct_pos]*monthly_flub_pct_pos + ///
    _b[monthly_flua_pct_pos]*monthly_flua_pct_pos + ///
    _b[monthly_sars_pct_pos]*monthly_sars_pct_pos + ///
    _b[adj_cases_10k]*adj_cases_10k + ///
    ln(population)) if pandemic == 1

gen pan_deathrate_norsv_nofft = 100000*pred_death_norsv_pan_nofft/(population/12) if pandemic == 1
gen rsv_attr_pan_nofft = pandemic_deaths_nofft - pred_death_norsv_pan_nofft
gen rsv_attr_rate_pan_nofft = 100000*rsv_attr_pan_nofft/(population/12)
egen total_rsv_deaths_pan_nofft = sum(rsv_attr_pan_nofft)
gen par_rsv_pan_nofft = 1 - (total_deaths_pan_nofft - total_rsv_deaths_pan_nofft)/total_deaths_pan_nofft

*--- Influenza A (without FFT, pandemic) ---
punaf if pandemic == 1, atspec(monthly_flua_pct_pos=0) eform

gen p_death_noflua_nofft = exp(_b[_cons] + _b[linear]*linear + ///
    _b[mosq]*mosq + ///
    _b[monthly_flub_pct_pos]*monthly_flub_pct_pos + ///
    _b[monthly_sars_pct_pos]*monthly_sars_pct_pos + ///
    _b[adj_cases_10k]*adj_cases_10k + ///
    _b[monthly_rsv_pct_pos]*monthly_rsv_pct_pos + ///
    ln(population)) if pandemic == 1

gen pan_deathrate_noflua_nofft = 100000*p_death_noflua_nofft/(population/12) if pandemic == 1
gen flua_attr_pan_nofft = pandemic_deaths_nofft - p_death_noflua_nofft
gen flua_attr_rate_pan_nofft = 100000*flua_attr_pan_nofft/(population/12)
egen total_flua_deaths_pan_nofft = sum(flua_attr_pan_nofft)
gen par_flua_pan_nofft = 1 - (total_deaths_pan_nofft - total_flua_deaths_pan_nofft)/total_deaths_pan_nofft

*--- Influenza B (without FFT, pandemic) ---
punaf if pandemic == 1, atspec(monthly_flub_pct_pos=0) eform

gen p_death_noflub_nofft = exp(_b[_cons] + _b[linear]*linear + ///
    _b[mosq]*mosq + ///
    _b[monthly_flua_pct_pos]*monthly_flua_pct_pos + ///
    _b[monthly_sars_pct_pos]*monthly_sars_pct_pos + ///
    _b[adj_cases_10k]*adj_cases_10k + ///
    _b[monthly_rsv_pct_pos]*monthly_rsv_pct_pos + ///
    ln(population)) if pandemic == 1

gen pan_deathrate_noflub_nofft = 100000*p_death_noflub_nofft/(population/12) if pandemic == 1
gen flub_attr_pan_nofft = pandemic_deaths_nofft - p_death_noflub_nofft
gen flub_attr_rate_pan_nofft = 100000*flub_attr_pan_nofft/(population/12)
egen total_flub_deaths_pan_nofft = sum(flub_attr_pan_nofft)
gen par_flub_pan_nofft = 1 - (total_deaths_pan_nofft - total_flub_deaths_pan_nofft)/total_deaths_pan_nofft

*--- All respiratory viruses (without FFT, pandemic) ---
punaf if pandemic == 1, atspec(monthly_rsv_pct_pos=0 monthly_flub_pct_pos=0 ///
    monthly_flua_pct_pos=0 monthly_sars_pct_pos=0 adj_cases_10k=0) eform

gen pred_death_no_virus_pan_nofft = exp(_b[_cons] + _b[linear]*linear + ///
    _b[mosq]*mosq + ln(population)) if pandemic == 1

gen pan_death_rate_no_virus_nofft = 100000*pred_death_no_virus_pan_nofft/(population/12) if pandemic == 1
gen virus_attr_death_pan_nofft = pandemic_deaths_nofft - pred_death_no_virus_pan_nofft
gen virus_attr_death_rate_pan_nofft = 100000*virus_attr_death_pan_nofft/(population/12)
egen total_virus_deaths_pan_nofft = sum(virus_attr_death_pan_nofft)
gen par_virus_pan_nofft = 1 - (total_deaths_pan_nofft - total_virus_deaths_pan_nofft)/total_deaths_pan_nofft

********************************************************************************
* 4. WALD TESTS FOR IMPACT OF FOURIER SEASONAL ADJUSTMENT
********************************************************************************

*-------------------------------------------------------------------------------
* 4a. Pre-pandemic: Test if coefficients differ between FFT and no-FFT models
*-------------------------------------------------------------------------------

* Combine estimates using seemingly unrelated estimation
suest prepan_fft prepan_nofft

* Test if flu A coefficient differs between models
test [prepan_fft_mean]monthly_flua_pct_pos = [prepan_nofft_mean]monthly_flua_pct_pos

* Test if flu B coefficient differs between models
test [prepan_fft_mean]monthly_flub_pct_pos = [prepan_nofft_mean]monthly_flub_pct_pos

* Test if RSV coefficient differs between models
test [prepan_fft_mean]monthly_rsv_pct_pos = [prepan_nofft_mean]monthly_rsv_pct_pos

*-------------------------------------------------------------------------------
* 4b. Pandemic: Test if coefficients differ between FFT and no-FFT models
*-------------------------------------------------------------------------------

* Combine estimates
suest pan_fft pan_nofft

* Test if flu A coefficient differs
test [pan_fft_mean]monthly_flua_pct_pos = [pan_nofft_mean]monthly_flua_pct_pos

* Test if flu B coefficient differs
test [pan_fft_mean]monthly_flub_pct_pos = [pan_nofft_mean]monthly_flub_pct_pos

* Test if RSV coefficient differs
test [pan_fft_mean]monthly_rsv_pct_pos = [pan_nofft_mean]monthly_rsv_pct_pos

* Test if SARS-CoV-2 (adjusted cases) coefficient differs
test [pan_fft_mean]adj_cases_10k = [pan_nofft_mean]adj_cases_10k

* Test if SARS-CoV-2 (% positivity) coefficient differs
test [pan_fft_mean]monthly_sars_pct_pos = [pan_nofft_mean]monthly_sars_pct_pos

********************************************************************************
* END OF ANALYSIS
********************************************************************************
* Results can be exported for meta-analysis using 02_meta_analysis.do
********************************************************************************
