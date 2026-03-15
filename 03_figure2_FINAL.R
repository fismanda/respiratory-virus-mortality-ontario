################################################################################
# FIGURE 2: MODEL-ESTIMATED VIRUS-ATTRIBUTABLE MORTALITY
# Reads counterfactual prediction data from Stata (figure2_data.csv)
# Three periods: pre-pandemic, PHEIC, post-PHEIC
# H1N1 gap (Nov 2009 - Aug 2010) shown as visual break
################################################################################

library(dplyr)
library(ggplot2)
library(patchwork)

################################################################################
# LOAD DATA
################################################################################

dat <- read.csv("data/figure2_data.csv", stringsAsFactors = FALSE)
dat$date <- as.Date(dat$date_str)

# H1N1 exclusion gap
h1n1_start <- as.Date("2009-11-01")
h1n1_end   <- as.Date("2010-08-31")

# Set attributable rates to NA in H1N1 gap so lines break cleanly
dat <- dat %>%
  mutate(across(c(flua_attr_fft, flub_attr_fft, rsv_attr_fft,
                  rsv_attr_nofft, combined_attr_fft),
                ~ ifelse(date >= h1n1_start & date <= h1n1_end, NA, .)))

# SARS only exists pandemic onward — already NA pre-2020 by construction
# Combined for pre-pandemic has no SARS component, which is correct

# Shading rectangles
h1n1_rect <- data.frame(xmin = h1n1_start, xmax = h1n1_end,
                         ymin = -Inf, ymax = Inf)

pandemic_start <- as.Date("2020-03-15")
pheic_end      <- as.Date("2023-04-30")

################################################################################
# PANEL A: Combined viruses, Influenza A, SARS-CoV-2
################################################################################

# Reasonable y ceiling — clip extreme outliers for display
# (single-month spikes in post-PHEIC from high % positivity months)
y_ceiling_A <- 200

pA <- ggplot(dat, aes(x = date)) +

  # H1N1 gap shading
  geom_rect(data = h1n1_rect,
            aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = y_ceiling_A),
            inherit.aes = FALSE, fill = "gray80", alpha = 0.4) +
  annotate("text", x = as.Date("2010-03-15"), y = y_ceiling_A * 0.92,
           label = "H1N1\ngap", size = 2.8, color = "gray50", vjust = 1) +

  # Combined all viruses (blue)
  geom_line(aes(y = pmin(combined_attr_fft, y_ceiling_A),
                color = "All viruses"),
            na.rm = FALSE, linewidth = 0.7) +

  # Influenza A (brown)
  geom_line(aes(y = pmin(flua_attr_fft, y_ceiling_A),
                color = "Influenza A"),
            na.rm = FALSE, linewidth = 0.7) +

  # SARS-CoV-2 (green, pandemic only — NAs elsewhere break the line)
  geom_line(aes(y = pmin(sars_attr_fft, y_ceiling_A),
                color = "SARS-CoV-2"),
            na.rm = FALSE, linewidth = 0.7) +

  # Period markers
  geom_vline(xintercept = pandemic_start,
             linetype = "dashed", linewidth = 0.6, color = "black") +
  geom_vline(xintercept = pheic_end,
             linetype = "dotted", linewidth = 0.5, color = "gray40") +

  scale_color_manual(
    values = c("All viruses" = "#4472C4",
               "Influenza A" = "#843C0C",
               "SARS-CoV-2"  = "#2E8B57"),
    breaks = c("All viruses", "Influenza A", "SARS-CoV-2")
  ) +

  scale_y_continuous(limits = c(0, y_ceiling_A),
                     breaks = seq(0, y_ceiling_A, 50)) +
  scale_x_date(date_breaks = "5 years", date_labels = "%Y",
               limits = c(as.Date("1993-01-01"), as.Date("2025-03-01"))) +

  labs(y = "Death Rate per 100,000", x = NULL, color = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position  = "bottom",
        panel.grid.minor = element_blank(),
        plot.tag = element_text(face = "bold", size = 13))

################################################################################
# PANEL B: RSV with and without Fourier adjustment
################################################################################

y_min_B   <- -25
y_ceil_B  <- 120

pB <- ggplot(dat, aes(x = date)) +

  # H1N1 gap shading
  geom_rect(data = h1n1_rect,
            aes(xmin = xmin, xmax = xmax, ymin = y_min_B, ymax = y_ceil_B),
            inherit.aes = FALSE, fill = "gray80", alpha = 0.4) +
  annotate("text", x = as.Date("2010-03-15"), y = y_ceil_B * 0.88,
           label = "H1N1\ngap", size = 2.8, color = "gray50", vjust = 1) +

  # Zero reference line
  geom_hline(yintercept = 0, linetype = "solid",
             color = "gray70", linewidth = 0.4) +

  # RSV with FFT (blue)
  geom_line(aes(y = pmax(pmin(rsv_attr_fft,  y_ceil_B), y_min_B),
                color = "RSV (with FFT)"),
            na.rm = FALSE, linewidth = 0.7) +

  # RSV without FFT (black)
  geom_line(aes(y = pmax(pmin(rsv_attr_nofft, y_ceil_B), y_min_B),
                color = "RSV (no FFT)"),
            na.rm = FALSE, linewidth = 0.7) +

  # Period markers
  geom_vline(xintercept = pandemic_start,
             linetype = "dashed", linewidth = 0.6, color = "black") +
  geom_vline(xintercept = pheic_end,
             linetype = "dotted", linewidth = 0.5, color = "gray40") +

  scale_color_manual(
    values = c("RSV (with FFT)" = "#4472C4",
               "RSV (no FFT)"   = "#222222")
  ) +

  scale_y_continuous(limits = c(y_min_B, y_ceil_B),
                     breaks = seq(-20, 120, 20)) +
  scale_x_date(date_breaks = "5 years", date_labels = "%Y",
               limits = c(as.Date("1993-01-01"), as.Date("2025-03-01"))) +

  labs(y = "RSV Death Rate per 100,000", x = NULL, color = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position  = "bottom",
        panel.grid.minor = element_blank(),
        plot.tag = element_text(face = "bold", size = 13))

################################################################################
# COMBINE AND SAVE
################################################################################

figure2 <- (pA + pB) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 13))

# Add shared x-axis label
figure2 <- figure2 & labs(x = "Date")

print(figure2)

ggsave("figures/figure2_attributable_mortality.png",  figure2,
       width = 12, height = 5.5, dpi = 300)
ggsave("figures/figure2_attributable_mortality.tiff", figure2,
       width = 12, height = 5.5, dpi = 300)

cat("Figure 2 saved.\n")

################################################################################
# QUICK SANITY CHECK - print range of attributable rates by period
################################################################################

cat("\n--- Sanity check: attributable rate ranges ---\n")
dat %>%
  group_by(period = case_when(
    prepandemic == 1 ~ "Pre-pandemic",
    pheic       == 1 ~ "PHEIC",
    postpheic   == 1 ~ "Post-PHEIC"
  )) %>%
  summarise(
    flua_min  = round(min(flua_attr_fft,     na.rm = TRUE), 1),
    flua_max  = round(max(flua_attr_fft,     na.rm = TRUE), 1),
    sars_min  = round(min(sars_attr_fft,     na.rm = TRUE), 1),
    sars_max  = round(max(sars_attr_fft,     na.rm = TRUE), 1),
    rsv_fft_min  = round(min(rsv_attr_fft,   na.rm = TRUE), 1),
    rsv_fft_max  = round(max(rsv_attr_fft,   na.rm = TRUE), 1),
    rsv_nofft_max = round(max(rsv_attr_nofft, na.rm = TRUE), 1)
  ) %>%
  print()
