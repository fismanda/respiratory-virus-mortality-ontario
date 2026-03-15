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
