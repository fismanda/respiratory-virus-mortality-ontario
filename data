# Data Dictionary

## Monthly Mortality and Virology Data (1991-2024)

**File:** `monthly_mortality_virology_data_ontario_1991_2024.csv`

This file contains monthly aggregated data for Ontario from January 1991 through June 2024, comprising the input variables needed to reproduce all analyses in the manuscript.

### Time Variables

| Variable | Type | Description |
|----------|------|-------------|
| `year` | Numeric | Calendar year (1991-2024) |
| `month` | Numeric | Calendar month (1-12) |
| `date` | Date | Mid-point of month (15th day), format: YYYY-MM-DD |

### Mortality Data

| Variable | Type | Description |
|----------|------|-------------|
| `deaths` | Numeric | Total all-cause deaths in Ontario for the month |
| `population` | Numeric | Ontario population estimate for the year (annual estimates) |

**Data Sources:**
- Deaths (1991-1993): Statistics Canada, Table 13-10-0392-01
- Deaths (1994-2024): Ontario Deaths Registry, Ontario Ministry of Health
- Population: Statistics Canada annual estimates

**Note:** Ontario Deaths Registry and Statistics Canada death counts show mean concordance of 0.3% where both are available (1994-2024).

### Virological Surveillance

| Variable | Type | Description | Time Coverage |
|----------|------|-------------|---------------|
| `monthly_flua_pct_pos` | Numeric | Influenza A percent positivity from FluWatch sentinel surveillance | 1993-2024 |
| `monthly_flub_pct_pos` | Numeric | Influenza B percent positivity from FluWatch sentinel surveillance | 1993-2024 |
| `monthly_rsv_pct_pos` | Numeric | RSV percent positivity from FluWatch sentinel surveillance | 1993-2024 |
| `monthly_sars_pct_pos` | Numeric | SARS-CoV-2 percent positivity from FluWatch sentinel surveillance | Sep 2022-2024 |

**Data Source:** Public Health Agency of Canada, FluWatch sentinel laboratory surveillance program

**Missing Values:**
- Influenza and RSV data: Missing for January 1991 - December 1992 (FluWatch surveillance began in 1993)
- **H1N1 pandemic exclusion:** Data for November 2009 - August 2010 excluded from analysis due to disruption of routine surveillance during the 2009 H1N1 pandemic
- SARS-CoV-2 positivity: Zero/missing before September 2022 (added to FluWatch surveillance in fall 2022)

### SARS-CoV-2 Case Data

| Variable | Type | Description | Time Coverage |
|----------|------|-------------|---------------|
| `monthly_adj_cases` | Numeric | Test-adjusted SARS-CoV-2 case counts | March 2020 - August 2022 |

**Data Source:** Ontario COVID-19 Science Advisory Table surveillance data

**Methodology:** Test-adjusted case counts correct for differential testing rates by age and sex. Methods described in:
- Fisman DN, et al. COVID-19 Case Age Distribution: Correction for Differential Testing by Age. Ann Intern Med. 2021;174(10):1430-1438. https://doi.org/10.7326/M20-7003
- Bosco S, et al. Impact of adjustment for differential testing by age and sex on apparent epidemiology of SARS-CoV-2 infection in Ontario, Canada. BMC Infect Dis. 2025;25:589. https://doi.org/10.1186/s12879-025-10968-6

**Note on SARS-CoV-2 exposure measures:**
- March 2020 - August 2022: Models use `monthly_adj_cases` (Science Table data)
- September 2022 - March 2024: Models use `monthly_sars_pct_pos` (FluWatch data)
- This transition reflects the end of the Science Table's surveillance activities in August 2022 and the addition of SARS-CoV-2 to FluWatch sentinel surveillance

### Data Completeness

| Time Period | Deaths | Population | Flu A/B | RSV | SARS-CoV-2 Cases | SARS-CoV-2 % Pos |
|-------------|--------|------------|---------|-----|------------------|------------------|
| 1991-1992 | ✓ | ✓ | ✗ | ✗ | ✗ | ✗ |
| 1993-Oct 2009 | ✓ | ✓ | ✓ | ✓ | ✗ | ✗ |
| Nov 2009-Aug 2010 | ✓ | ✓ | ✗* | ✗* | ✗ | ✗ |
| Sep 2010-Feb 2020 | ✓ | ✓ | ✓ | ✓ | ✗ | ✗ |
| Mar 2020-Aug 2022 | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ |
| Sep 2022-Mar 2024 | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ |

*H1N1 pandemic period excluded from analysis due to surveillance disruption

---

## Meta-Analysis Data

**File:** `meta_analysis_par_estimates.csv`

This file contains population attributable risk (PAR) estimates stratified by virus type, time period (pre-pandemic vs. pandemic), and seasonal adjustment method (with vs. without Fourier transforms).

### Variables

| Variable | Type | Description |
|----------|------|-------------|
| `pandemic` | Binary | Time period: 0 = pre-pandemic (1993-2020), 1 = pandemic (2020-2024) |
| `virus` | Character | Virus type: "Influenza A", "Influenza B", "Respiratory Syncitial Virus", "SARS-CoV-2", "All Respiratory Viruses" |
| `fft` | Binary | Seasonal adjustment: 0 = without Fourier terms, 1 = with Fourier terms |
| `par` | Numeric | Population attributable risk (proportion, 0-1 scale) |
| `lcl` | Numeric | Lower 95% confidence limit |
| `ucl` | Numeric | Upper 95% confidence limit |
| `lnpar` | Numeric | Natural log of PAR (for meta-analysis) |
| `selnpar` | Numeric | Standard error of ln(PAR) |
| `fft_label` | Character | "yes" or "no" (for plotting) |
| `pandemic_label` | Character | "yes" or "no" (for plotting) |
| `subgroup` | Character | Combined label for forest plot stratification |

**Notes:**
- SARS-CoV-2 estimates only available for pandemic period (by definition)
- Some PAR estimates are negative (particularly RSV with Fourier adjustment), indicating paradoxical associations when seasonal adjustment removes true viral signal
- Used in `02_meta_analysis.do` to generate Figure 3 and assess heterogeneity

---

## Notes on Reproducibility

### Data Truncation
The analysis uses data through March 2024 (`series_month < 376` in code) to account for reporting delays in death registration.

### Missing Data Handling
- Missing virological surveillance data (1991-1992): These months are excluded from regression models
- **H1N1 pandemic period (November 2009 - August 2010):** Excluded from analysis due to disruption of routine influenza and RSV surveillance during the 2009 H1N1 pandemic. This period is visible as a gap in Figure 2.
- SARS-CoV-2 exposure variables are zero/missing for periods when data are not available, appropriately handled by time-stratified models

### Pandemic Definition
The pandemic period begins March 1, 2020, corresponding to the emergence of SARS-CoV-2 in Ontario and implementation of public health measures that affected both disease transmission and healthcare-seeking behavior.

---

## Data Provenance

All data in this repository are derived from publicly available sources or have been aggregated and anonymized. Individual-level data are subject to data sharing agreements and cannot be publicly released. For questions about data sources or methodology, contact the corresponding author.

**Last Updated:** December 2025
