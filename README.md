# Population Attributable Mortality — Respiratory Viruses, Ontario

**Fisman DN, Grima AA, Wilson NJ, Mann S, Tuite AR, Lee CE.**  
Population Attributable Mortality Associated with Respiratory Viruses in Ontario.  
*Submitted to The Lancet, March 2026.*

---

## Repository contents

| File | Description |
|------|-------------|
| `01_analysis_EXTENDED_v4.do` | Main Stata analysis (negative binomial models, punaf PAF estimation, meta-analysis) |
| `ontario_mortality_virology_EXTENDED_1991_2025_CORRECTED.csv` | Monthly all-cause mortality and viral surveillance data, Ontario 1991–2025 |
| `AJR_and_data_for_meta-regression.csv` | PAF estimates and CIs for random-effects meta-analysis (all confirmed from punaf output) |
| `figure2.R` | R script: Figure 2 (model-estimated virus-attributable mortality over time) |
| `figure3_forest_plot.R` | R script: Figure 3 (forest plot, PAF meta-analysis by virus and modeling approach) |
| `figure3_forest_plot.png` | Figure 3 output (300 dpi) |

---

## Data sources

- **Mortality:** Ontario Deaths Registry (1994–February 2025); Statistics Canada (1991–1993)
- **Viral surveillance:** FluWatch (Public Health Agency of Canada), influenza A/B and RSV percent positivity, 1993–present
- **SARS-CoV-2:** Test-adjusted case counts from Ontario CCM system (March 2020–August 2022); FluWatch percent positivity (September 2022 onward)
- **Population denominators:** Statistics Canada

---

## Analysis

All primary analyses conducted in **Stata 18**. PAFs estimated using `punaf` (counterfactual approach). Random-effects meta-analyses using `metan` (DerSimonian-Laird). Figures produced in **R** (ggplot2).

### Key model periods
| Period | n (months) | Notes |
|--------|-----------|-------|
| Pre-pandemic | 308 | January 1993–February 2020 |
| Combined pandemic | 60 | March 2020–February 2025 — **primary** |
| PHEIC | 38 | March 2020–April 2023 — exploratory |
| Post-PHEIC | 22 | May 2023–February 2025 — exploratory |
| Restricted pandemic | 56 | July 2020–February 2025 — sensitivity |

---

## Key results (primary models, with Fourier seasonal adjustment)

| Virus | Period | PAF | 95% CI |
|-------|--------|-----|--------|
| Influenza A | Pre-pandemic | 1.8% | 1.4–2.3% |
| SARS-CoV-2 | Combined pandemic | 6.1% | 4.2–8.0% |
| SARS-CoV-2 | PHEIC | 7.3% | 5.0–9.6% |
| SARS-CoV-2 | Post-PHEIC | 9.8% | 1.1–17.7% |

---

## Ethics

University of Toronto Research Ethics Board, Protocol #41690.

---

## Archive

This repository is archived at Zenodo: [https://zenodo.org/records/17887720](https://zenodo.org/records/17887720)  
*(Update this URL after creating the new Zenodo version)*
