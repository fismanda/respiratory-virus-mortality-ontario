################################################################################
# Figure 3 — Forest plot of population attributable fractions (PAF)
# Descriptive forest plot (NO meta-analysis): stratum-specific PAF + 95% CI.
# Values are taken from the definitive Tables 1 and 2 (with/without FFT) plus
# the four exploratory SARS-CoV-2 sub-period estimates (PHEIC / post-PHEIC).
# Reads data/AJR_forest_final.csv
################################################################################

library(readr)
library(dplyr)
library(ggplot2)

# ---- 1. Load (values stored as proportions; convert to %) ----
d <- read_csv("data/AJR_forest_final.csv", show_col_types = FALSE) %>%
  mutate(par = par * 100, lcl = lcl * 100, ucl = ucl * 100)

# ---- 2. Tidy names + stratum labels ----
d <- d %>%
  mutate(fft_lab = ifelse(fft == 1, "with FFT", "no FFT"))

d <- d %>%
  mutate(stratum = ifelse(
    virus == "SARS-CoV-2",
    paste0(subperiod, ", ", fft_lab),
    paste0(ifelse(pandemic == 1, "Pandemic", "Pre-pandemic"), ", ", fft_lab)
  ))

# ---- 3. Ordering ----
virus_order <- c("All Respiratory Viruses", "Influenza A", "Influenza B",
                 "Respiratory Syncytial Virus", "SARS-CoV-2")
nonsars_order <- c("Pre-pandemic, no FFT", "Pre-pandemic, with FFT",
                   "Pandemic, no FFT", "Pandemic, with FFT")
sars_order <- c("Combined, no FFT", "Combined, with FFT",
                "PHEIC, no FFT", "PHEIC, with FFT",
                "Post-PHEIC, no FFT", "Post-PHEIC, with FFT")

ord_index <- function(virus, stratum) {
  if (virus == "SARS-CoV-2") match(stratum, sars_order)
  else match(stratum, nonsars_order)
}
d$ord <- mapply(ord_index, d$virus, d$stratum)
d$virus_grp <- factor(d$virus, levels = virus_order)
d <- d %>% arrange(virus_grp, ord)

# ---- 4. y positions with headers + gaps ----
rows <- list(); header_rows <- list(); ypos <- 0
for (vg in virus_order) {
  sub <- d %>% filter(virus_grp == vg)
  if (nrow(sub) == 0) next
  ypos <- ypos - 1
  header_rows[[length(header_rows)+1]] <- data.frame(virus_grp = vg, y = ypos)
  for (i in seq_len(nrow(sub))) {
    ypos <- ypos - 1
    rows[[length(rows)+1]] <- cbind(sub[i, ], y = ypos)
  }
  ypos <- ypos - 0.5
}
plotd   <- bind_rows(rows)
headers <- bind_rows(header_rows)

plotd <- plotd %>%
  mutate(period = ifelse(pandemic == 1, "Pandemic era", "Pre-pandemic"))

# ---- 5. Plot ----
p <- ggplot(plotd, aes(x = par, y = y, colour = period)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_errorbarh(aes(xmin = lcl, xmax = ucl), height = 0.25, linewidth = 0.6) +
  geom_point(size = 2.6) +
  geom_text(data = headers, aes(x = -7, y = y, label = virus_grp),
            inherit.aes = FALSE, fontface = "bold", hjust = 0, size = 3.6) +
  geom_text(data = plotd, aes(x = -7, y = y, label = stratum),
            inherit.aes = FALSE, hjust = 0, size = 3.0, colour = "grey25") +
  geom_text(data = plotd,
            aes(x = 19, y = y,
                label = sprintf("%.1f%% (%.1f, %.1f)", par, lcl, ucl)),
            inherit.aes = FALSE, hjust = 0, size = 3.0) +
  scale_colour_manual(values = c("Pre-pandemic" = "#4C72B0",
                                 "Pandemic era" = "#C44E52"),
                      name = NULL) +
  scale_x_continuous("Population attributable fraction (%)",
                     limits = c(-7.5, 26),
                     breaks = seq(0, 18, 2)) +
  theme_minimal(base_size = 11) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(colour = "grey92"),
    axis.line.x = element_line(colour = "grey70"),
    legend.position = "bottom",
    plot.margin = margin(10, 14, 10, 10)
  )

ggsave("figure3_forest_plot.png", p, width = 9.5, height = 9, dpi = 300)
message("Wrote figure3_forest_plot.png")
