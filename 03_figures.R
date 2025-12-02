################################################################################
# RESPIRATORY VIRUS MORTALITY - FIGURES 1 AND 2
# Creates publication-quality figures for mortality and virology surveillance
################################################################################

library(haven)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

################################################################################
# LOAD DATA
################################################################################

# Load the analysis dataset (output from 01_analysis.do)
# This should be the .dta file AFTER running the Stata analysis
dat <- read_dta("data/MONTHLY_VIROLOGY_AND_DEATH_DATA_FOR_ANALYSIS.dta")

# Convert Stata date format to R date
dat <- dat %>%
  mutate(date = as.Date(date, origin = "1960-01-01"))

################################################################################
# FIGURE 1: MORTALITY AND VIROLOGICAL SURVEILLANCE
################################################################################

#-------------------------------------------------------------------------------
# Panel A: Deaths per 100,000
#-------------------------------------------------------------------------------

df_deaths <- dat %>%
  mutate(year = as.integer(format(date, "%Y")),
         source = case_when(
           series_month < 376 & year >= 1994 ~ "Ontario Deaths Registry",
           year < 1994 ~ "Statistics Canada"
         )) %>%
  filter(!is.na(source))

p1 <- ggplot(df_deaths, aes(x = date, y = deaths100k, color = source)) +
  geom_line() +
  geom_vline(xintercept = as.Date("2020-03-15"), linetype = "dashed") +
  labs(y = "Deaths per 100,000", x = "Date", color = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")

#-------------------------------------------------------------------------------
# Panel B: Influenza A and B percent positivity
#-------------------------------------------------------------------------------

df_flu <- dat %>%
  select(date, monthly_flua_pct_pos, monthly_flub_pct_pos) %>%
  pivot_longer(cols = c(monthly_flua_pct_pos, monthly_flub_pct_pos),
               names_to = "virus", values_to = "positivity") %>%
  mutate(virus = recode(virus,
                        monthly_flua_pct_pos = "Influenza A",
                        monthly_flub_pct_pos = "Influenza B"))

p2 <- ggplot(df_flu, aes(x = date, y = positivity, color = virus)) +
  geom_line() +
  geom_vline(xintercept = as.Date("2020-03-15"), linetype = "dashed") +
  labs(y = "Percent Positive", x = "Date", color = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")

#-------------------------------------------------------------------------------
# Panel C: RSV percent positivity
#-------------------------------------------------------------------------------

df_rsv <- dat %>%
  select(date, monthly_rsv_pct_pos)

p3 <- ggplot(df_rsv, aes(x = date, y = monthly_rsv_pct_pos)) +
  geom_line(color = "orange") +
  geom_vline(xintercept = as.Date("2020-03-15"), linetype = "dashed") +
  labs(y = "RSV Percent Positive", x = "Date") +
  theme_minimal() +
  theme(legend.position = "bottom")

#-------------------------------------------------------------------------------
# Panel D: SARS-CoV-2 adjusted cases and percent positivity
#-------------------------------------------------------------------------------

# Prepare SARS-CoV-2 data
df_sars <- dat %>%
  select(date, adj_cases_10k, monthly_sars_pct_pos)

# Adjusted cases: March 2020 - August 2022
df_cases <- df_sars %>%
  filter(date >= as.Date("2020-03-15") & date <= as.Date("2022-08-31"))

# Percent positivity: September 2022 - August 2024
df_pos <- df_sars %>%
  filter(date >= as.Date("2022-09-01") & date <= as.Date("2024-08-31"))

# Calculate scaling factor for dual y-axis
scale_factor <- max(df_cases$adj_cases_10k, na.rm = TRUE) /
  max(df_pos$monthly_sars_pct_pos, na.rm = TRUE)

p4 <- ggplot() +
  geom_line(data = df_cases,
            aes(x = date, y = adj_cases_10k, color = "Adjusted Cases per 10,000")) +
  geom_line(data = df_pos,
            aes(x = date, y = monthly_sars_pct_pos * scale_factor,
                color = "SARS-2 Positivity (%)")) +
  scale_y_continuous(
    name = "Adjusted Cases per 10,000",
    sec.axis = sec_axis(~./scale_factor, name = "SARS-2 Positivity (%)")
  ) +
  labs(x = "Date", color = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")

#-------------------------------------------------------------------------------
# Combine panels and save Figure 1
#-------------------------------------------------------------------------------

figure1 <- (p1 + p2) / (p3 + p4) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))

# View
print(figure1)

# Save to file
ggsave("figures/figure1_mortality_virology.png", figure1, 
       width = 10, height = 8, dpi = 300)
ggsave("figures/figure1_mortality_virology.tiff", figure1, 
       width = 10, height = 8, dpi = 300)

################################################################################
# FIGURE 2: MODEL-ESTIMATED VIRUS-ATTRIBUTABLE MORTALITY
################################################################################

# Note: The data gaps (2009-2010) are for H1N1 pandemic exclusion period

#-------------------------------------------------------------------------------
# Panel A: Combined viruses, Influenza A, and SARS-CoV-2
#-------------------------------------------------------------------------------

pA <- ggplot(dat) +
  # All viruses combined (blue)
  geom_line(data = dat %>% filter(pandemic == 0 & date < as.Date("2009-11-15")),
            aes(x = date, y = virus_attr_death_rate_pre, color = "All viruses")) +
  geom_line(data = dat %>% filter(pandemic == 0 & date > as.Date("2010-08-15")),
            aes(x = date, y = virus_attr_death_rate_pre, color = "All viruses")) +
  geom_line(data = dat %>% filter(pandemic == 1),
            aes(x = date, y = virus_attr_death_rate_pan, color = "All viruses")) +
  
  # Influenza A (brown)
  geom_line(data = dat %>% filter(pandemic == 0 & date < as.Date("2009-11-15")),
            aes(x = date, y = flua_attr_death_rate_pre, color = "Influenza A")) +
  geom_line(data = dat %>% filter(pandemic == 0 & date > as.Date("2010-08-15")),
            aes(x = date, y = flua_attr_death_rate_pre, color = "Influenza A")) +
  geom_line(data = dat %>% filter(pandemic == 1),
            aes(x = date, y = flua_attr_death_rate_pan, color = "Influenza A")) +
  
  # SARS-CoV-2 (green, pandemic only)
  geom_line(data = dat %>% filter(pandemic == 1),
            aes(x = date, y = sars2_attr_death_rate, color = "SARS-CoV-2")) +
  
  # Manual color scale
  scale_color_manual(values = c(
    "All viruses" = scales::alpha("blue", 0.6),
    "Influenza A" = scales::alpha("brown", 0.6),
    "SARS-CoV-2" = scales::alpha("darkgreen", 0.6)
  )) +
  
  # Formatting
  geom_vline(xintercept = as.Date("2020-03-15"), linetype = "dashed", color = "black") +
  scale_y_continuous(limits = c(0, 250), breaks = seq(0, 250, 50)) +
  scale_x_date(date_breaks = "10 years", date_labels = "%Y") +
  labs(y = "Death Rate per 100,000", x = "Date", color = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  annotate("text", x = as.Date("1993-01-01"), y = 250, label = "A",
           size = 8, hjust = 0, vjust = 1, fontface = "bold")

#-------------------------------------------------------------------------------
# Panel B: RSV with and without Fourier seasonal adjustment
#-------------------------------------------------------------------------------

pB <- ggplot(dat) +
  # RSV with FFT (blue)
  geom_line(data = dat %>% filter(pandemic == 0 & date < as.Date("2009-11-15")),
            aes(x = date, y = rsv_attr_death_rate_pre, color = "RSV (with FFT)")) +
  geom_line(data = dat %>% filter(pandemic == 0 & date > as.Date("2010-08-15")),
            aes(x = date, y = rsv_attr_death_rate_pre, color = "RSV (with FFT)")) +
  geom_line(data = dat %>% filter(pandemic == 1),
            aes(x = date, y = rsv_attr_death_rate_pan, color = "RSV (with FFT)")) +
  
  # RSV without FFT (black)
  geom_line(data = dat %>% filter(pandemic == 0 & date > as.Date("2010-08-15")),
            aes(x = date, y = rsv_attr_rate_pre_nofft, color = "RSV (no FFT)")) +
  geom_line(data = dat %>% filter(pandemic == 0 & date < as.Date("2009-11-15")),
            aes(x = date, y = rsv_attr_rate_pre_nofft, color = "RSV (no FFT)")) +
  geom_line(data = dat %>% filter(pandemic == 1),
            aes(x = date, y = rsv_attr_rate_pan_nofft, color = "RSV (no FFT)")) +
  
  # Manual color scale
  scale_color_manual(values = c(
    "RSV (with FFT)" = scales::alpha("blue", 0.6),
    "RSV (no FFT)" = scales::alpha("black", 0.6)
  )) +
  
  geom_vline(xintercept = as.Date("2020-03-15"), linetype = "dashed", color = "black") +
  scale_y_continuous(limits = c(-20, 260), breaks = seq(-20, 250, 50)) +
  scale_x_date(date_breaks = "10 years", date_labels = "%Y") +
  labs(y = "RSV Death Rate per 100,000", x = "Date", color = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  annotate("text", x = as.Date("1993-01-01"), y = 260, label = "B",
           size = 8, hjust = 0, vjust = 1, fontface = "bold")

#-------------------------------------------------------------------------------
# Combine panels and save Figure 2
#-------------------------------------------------------------------------------

figure2 <- pA + pB + 
  plot_layout(ncol = 2) &
  theme(plot.margin = margin(0, 0, 0, 0, "cm"))

# View
print(figure2)

# Save to file
ggsave("figures/figure2_attributable_mortality.png", figure2, 
       width = 12, height = 6, dpi = 300)
ggsave("figures/figure2_attributable_mortality.tiff", figure2, 
       width = 12, height = 6, dpi = 300)

################################################################################
# END OF FIGURE GENERATION
################################################################################
# Figure 3 (forest plot) is generated in Stata using 02_meta_analysis.do
################################################################################
