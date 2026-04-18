# ================================================================
# Script 8: Epidemiological Analysis & Early Warning Evaluation
# ================================================================

library(sf); library(dplyr); library(ggplot2); library(zoo)
library(lubridate); library(gridExtra); library(scales); library(tidyr)
library(viridis); library(conflicted); library(grid)


conflict_prefer("lag", "dplyr")
conflicted::conflicts_prefer(dplyr::filter)
output_dir <- "C:/Users/odeyr/OneDrive/Documents/R/R projects and exercises/2025 Mathematical Statistics Research  Project/Wastewater modelling/EDA_Output"
fig_dir <- file.path(output_dir, "figures")
res_dir <- file.path(output_dir, "results")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

# Clear environment
rm(list = ls()[!ls() %in% c("output_dir", "fig_dir", "res_dir")])

cat("=== SCRIPT 8: EPIDEMIOLOGICAL ANALYSIS & EARLY WARNING EVALUATION ===\n\n")

# --- Load all required data ---

load(file.path(output_dir, "master_analysis.RData"))

rk_results <- read.csv(file.path(res_dir, "rk_weekly_cv_results.csv"))
rk_overall <- read.csv(file.path(res_dir, "rk_weekly_cv_overall.csv"))

comp_results <- read.csv(file.path(output_dir, "competing_models_results.csv"))

load(file.path(output_dir, "variable_selection_results.RData"))

cat("Data loaded successfully\n")
cat("  - Analysis data:", nrow(analysis_data), "rows\n")
cat("  - RK CV folds:", nrow(rk_results), "\n")
cat("  - Competing models:", nrow(comp_results), "\n")
cat("  - RK columns:", paste(names(rk_overall), collapse = ", "), "\n")
cat("  - Comp columns:", paste(names(comp_results), collapse = ", "), "\n")

# --- Aggregate to municipality level ---
daily <- analysis_data %>%
  group_by(report_date) %>%
  summarise(total_cases = sum(daily_cases_sp, na.rm = TRUE),
            mean_gene = mean(gene_copies_ml, na.rm = TRUE), .groups = "drop") %>%
  arrange(report_date) %>%
  mutate(cases_7ma = rollapplyr(total_cases, 7, mean, fill = NA))

weekly <- analysis_data %>%
  group_by(weekly_number) %>%
  summarise(week_start = min(report_date),
            total_cases = sum(daily_cases_sp, na.rm = TRUE),
            mean_gene = mean(gene_copies_ml, na.rm = TRUE), .groups = "drop") %>%
  arrange(week_start) %>% 
  filter(!is.na(weekly_number))

n_w <- nrow(weekly)
cat("  - Weekly aggregated:", n_w, "weeks\n")

# --- Wave definitions ---
waves <- tibble(
  wave = paste0("Wave ", 1:4),
  start = as.Date(c("2020-07-01", "2020-11-01", "2021-05-01", "2021-11-01")),
  end   = as.Date(c("2020-10-31", "2021-03-31", "2021-09-30", "2022-01-31"))
)
waves$duration <- as.numeric(waves$end - waves$start + 1)
wcols <- c("#E8D44D", "#E07B39", "#D1495B", "#8B5CF6")
wcols_n <- setNames(wcols, waves$wave)

assign_wave <- function(d) {
  for (i in seq_len(nrow(waves))) {
    if (!is.na(d) && d >= waves$start[i] && d <= waves$end[i]) return(waves$wave[i])
  }
  return(NA_character_)
}
daily$wave <- sapply(daily$report_date, assign_wave)
weekly$wave <- sapply(weekly$week_start, assign_wave)

# ================================================================
# PART A: EPIDEMIOLOGICAL METRICS (Figures 19-22)
# ================================================================
cat("\n=== PART A: EPIDEMIOLOGICAL METRICS (Figures 19-22) ===\n")

# --- Figure 19: Effective Reproduction Number (Rt) ---
cat("  Computing Rt...\n")
daily$sum7     <- rollapplyr(daily$total_cases, 7, sum, fill = NA)
daily$sum7_lag <- dplyr::lag(daily$sum7, 7)
daily$rt       <- ifelse(daily$sum7_lag > 0, daily$sum7 / daily$sum7_lag, NA)
daily$rt       <- pmax(pmin(daily$rt, 4), 0)

p19 <- ggplot(daily, aes(x = report_date)) +
  geom_rect(data = waves, aes(xmin = start, xmax = end, fill = wave,
                              ymin = -Inf, ymax = Inf),
            alpha = 0.15, inherit.aes = FALSE) +
  geom_ribbon(aes(ymin = 0, ymax = pmin(rt, 1)), fill = "#2ECC71", alpha = 0.3, na.rm = TRUE) +
  geom_ribbon(aes(ymin = 1, ymax = pmax(rt, 1)), fill = "#E74C3C", alpha = 0.3, na.rm = TRUE) +
  geom_line(aes(y = rt), color = "#2C3E50", linewidth = 0.6, na.rm = TRUE) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey30") +
  scale_fill_manual(values = wcols_n, guide = guide_legend(override.aes = list(alpha = 0.3))) +
  labs(
       x = NULL, y = "Rt", fill = NULL) +
  theme_minimal(base_size = 11) +
  coord_cartesian(ylim = c(0, 4)) +
  theme(legend.position = "top")
ggsave(file.path(fig_dir, "fig19_reproduction_number.png"), p19, width = 10, height = 5, dpi = 600)

rm(p19); cat("  [OK] Figure 19 saved\n")

# --- Figure 20: Growth Rate & Doubling Time ---
cat("  Computing growth rates...\n")
daily$growth_rate <- (daily$cases_7ma - dplyr::lag(daily$cases_7ma, 7)) /
  (dplyr::lag(daily$cases_7ma, 7) + 1e-10) * 100
daily$growth_ma14 <- rollapplyr(daily$growth_rate, 14, mean, fill = NA)
daily$doubling_time <- ifelse(daily$growth_rate > 0,
                              7 * log(2) / log(1 + daily$growth_rate / 100), NA)
daily$dt_ma14 <- rollapplyr(daily$doubling_time, 14, mean, fill = NA)

p20a <- ggplot(daily, aes(x = report_date)) +
  geom_rect(data = waves, aes(xmin = start, xmax = end, fill = wave,
                              ymin = -Inf, ymax = Inf),
            alpha = 0.12, inherit.aes = FALSE) +
  geom_col(aes(y = growth_rate), fill = "#3498DB", alpha = 0.4, width = 1) +
  geom_line(aes(y = growth_ma14), color = "#E74C3C", linewidth = 0.7, na.rm = TRUE) +
  labs(x = NULL, y = "Growth Rate (%)") +
  theme_minimal(base_size = 10) +
  scale_fill_manual(values = wcols_n, guide = "none")

p20b <- ggplot(daily %>% filter(growth_rate > 0 & doubling_time < 200),
               aes(x = report_date)) +
  geom_rect(data = waves, aes(xmin = start, xmax = end, fill = wave,
                              ymin = -Inf, ymax = Inf),
            alpha = 0.12, inherit.aes = FALSE) +
  geom_point(aes(y = doubling_time), color = "#8E44AD", alpha = 0.3, size = 0.8) +
  geom_line(aes(y = dt_ma14), color = "#E74C3C", linewidth = 0.7, na.rm = TRUE) +
  geom_hline(yintercept = 7, linetype = "dotted", color = "grey50") +
  geom_hline(yintercept = 14, linetype = "dashed", color = "grey50") +
  labs(x = NULL, y = "Days") +
  theme_minimal(base_size = 10) +
  scale_fill_manual(values = wcols_n, guide = "none")

p20 <- grid.arrange(p20a, p20b, nrow = 2)
ggsave(file.path(fig_dir, "fig20_growth_doubling.png"), p20, width = 10, height = 7, dpi = 600)

rm(p20, p20a, p20b); cat("  [OK] Figure 20 saved\n")

# --- Figure 21: Case Fatality Rate ---
cat("  Computing CFR...\n")
has_deaths <- "daily_deaths_sp" %in% names(analysis_data) &&
  sum(!is.na(analysis_data$daily_deaths_sp)) > 0

if (!has_deaths) {
  deaths_path <- file.path(dirname(output_dir),
                           "WASTEWATER DATA/Cleaned datasets/ekurhuleni_covid_cleaned.csv")
  deaths_df <- tryCatch(read.csv(deaths_path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (!is.null(deaths_df) && "eku_total_deaths" %in% names(deaths_df)) {
    deaths_df$report_date <- as.Date(deaths_df$report_date)
    
    deaths_df <- deaths_df %>%
      arrange(report_date) %>%
      mutate(daily_deaths = c(NA, diff(eku_total_deaths))) %>%
      filter(!is.na(daily_deaths)) %>%
      mutate(daily_deaths = pmax(daily_deaths, 0))
    dd <- deaths_df %>% group_by(report_date) %>%
      summarise(total_deaths = sum(daily_deaths, na.rm = TRUE), .groups = "drop")
    daily <- daily %>% left_join(dd, by = "report_date")
    daily$total_deaths[is.na(daily$total_deaths)] <- 0
    has_deaths <- TRUE
    cat("    Deaths loaded from CSV (cumulative -> daily via diff)\n")
  }
}

if (has_deaths && "total_deaths" %in% names(daily)) {
  daily$cfr <- rollapplyr(daily$total_deaths, 30, sum, fill = NA) /
    rollapplyr(daily$total_cases, 30, sum, fill = NA) * 100
  # Cap CFR at reasonable range for display
  daily$cfr <- pmin(daily$cfr, 15)
  
  p21 <- ggplot(daily, aes(x = report_date, y = cfr)) +
    geom_rect(data = waves, aes(xmin = start, xmax = end, fill = wave,
                                ymin = -Inf, ymax = Inf),
              alpha = 0.15, inherit.aes = FALSE) +
    geom_line(color = "#C0392B", linewidth = 0.7, na.rm = TRUE) +
    geom_smooth(color = "darkred", linewidth = 0.5, se = TRUE, alpha = 0.2, na.rm = TRUE) +
    scale_fill_manual(values = wcols_n, guide = guide_legend(override.aes = list(alpha = 0.3))) +
    labs(
         x = NULL, y = "CFR (%)") +
    theme_minimal(base_size = 11) + theme(legend.position = "top")
} else {
  cat("  [!] Deaths data not available; producing placeholder\n")
  p21 <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, size = 5, hjust = 0.5,
             label = "Deaths data not available at subplace level.\nMunicipality-level CFR requires external death data.") +
    theme_void()
}
ggsave(file.path(fig_dir, "fig21_cfr.png"), p21, width = 10, height = 5, dpi = 600)
rm(p21); cat("  [OK] Figure 21 saved\n")

# --- Figure 22: Wave Summary Table ---

wave_list <- list()
for (i in 1:nrow(waves)) {
  w <- waves[i, ]
  mask <- daily$report_date >= w$start & daily$report_date <= w$end
  daily_in_wave <- daily[mask, ]
  
  peak_idx <- which.max(daily_in_wave$total_cases)
  
  wave_list[[i]] <- data.frame(
    wave = w$wave,
    start = as.character(w$start),
    end = as.character(w$end),
    duration_days = w$duration,
    total_cases = sum(daily_in_wave$total_cases, na.rm = TRUE),
    peak_daily_cases = if (sum(!is.na(daily_in_wave$total_cases)) > 0) 
      max(daily_in_wave$total_cases, na.rm = TRUE) else NA_real_,
    peak_date = if (length(peak_idx) > 0 && nrow(daily_in_wave) > 0) 
      as.character(daily_in_wave$report_date[peak_idx]) else NA_character_,
    mean_rt = mean(daily_in_wave$rt, na.rm = TRUE),
    mean_growth = mean(daily_in_wave$growth_rate, na.rm = TRUE),
    wave_cfr = if ("cfr" %in% names(daily_in_wave)) 
      mean(daily_in_wave$cfr, na.rm = TRUE) else NA_real_,
    stringsAsFactors = FALSE
  )
}
wave_summary <- bind_rows(wave_list)

cat("\n  --- Wave Summary ---\n")
for (i in seq_len(nrow(wave_summary))) {
  w <- wave_summary[i, ]
  cat(sprintf("  %s (%s to %s, %d days):\n", w$wave, w$start, w$end, w$duration_days))
  cat(sprintf("    Total cases: %s  |  Peak: %s on %s  |  Mean Rt: %.2f  |  CFR: %.2f%%\n",
              format(w$total_cases, big.mark = ","),
              format(w$peak_daily_cases, big.mark = ","),
              w$peak_date,
              w$mean_rt,
              ifelse(is.na(w$wave_cfr), NA, w$wave_cfr)))
}

ws_print <- wave_summary %>% 
  mutate(across(where(is.numeric), ~round(., 2)))

p22 <- tableGrob(ws_print, rows = NULL, theme = ttheme_minimal(base_size = 9))
p22 <- arrangeGrob(p22)
ggsave(file.path(fig_dir, "fig22_wave_summary.png"), p22, width = 16, height = 4, dpi = 600)
rm(p22); cat("  [OK] Figure 22 saved\n")


cat("\n--- WAVE-BY-WAVE EPIDEMIOLOGICAL SUMMARY (LaTeX) ---\n")
cat("\\begin{tabular}{lrrrrrrr}\n\\toprule\n")
cat("Wave & Start & End & Duration & Total Cases & Peak & Mean Rt & CFR (\\%) \\\\\n\\midrule\n")
for (i in seq_len(nrow(wave_summary))) {
  w <- wave_summary[i, ]
  cat(sprintf("%s & %s & %s & %d & %s & %s & %.2f & %s \\\\\n",
              w$wave, w$start, w$end, w$duration_days,
              format(w$total_cases, big.mark = ","),
              format(w$peak_daily_cases, big.mark = ","),
              w$mean_rt,
              ifelse(is.na(w$wave_cfr), "---", sprintf("%.2f", w$wave_cfr))))
}
cat("\\bottomrule\n\\end{tabular}\n\n")

# ================================================================
# PART B: EARLY WARNING EVALUATION (Figures 23-26)
# ================================================================
cat("\n=== PART B: EARLY WARNING EVALUATION (Figures 23-26) ===\n")

# --- Define outbreak threshold ---
case_p75 <- quantile(weekly$total_cases, 0.75, na.rm = TRUE)
weekly$outbreak <- weekly$total_cases > case_p75
cat(sprintf("Outbreak threshold (75th pct): %.0f weekly cases\n", case_p75))
cat(sprintf("Outbreak weeks: %d of %d (%.1f%%)\n", 
            sum(weekly$outbreak), n_w, 100 * mean(weekly$outbreak)))

# --- Evaluate Gene Threshold Method ---
cat("  Evaluating gene threshold method...\n")
thresholds <- quantile(weekly$mean_gene, c(.50, .60, .70, .75, .80, .90), na.rm = TRUE)

eval_df <- data.frame(
  method = "Gene Threshold",
  threshold = thresholds,
  sens = NA, spec = NA, ppv = NA, npv = NA, tp = NA, fp = NA, tn = NA, fn = NA
)

for (t in seq_along(thresholds)) {
  gh <- weekly$mean_gene >= thresholds[t]
  warned <- sapply(seq_len(n_w), function(i) {
    if (i <= 4) return(FALSE)
    any(gh[(i - 4):(i - 1)], na.rm = TRUE)
  })
  tp <- sum(weekly$outbreak & warned);  fn <- sum(weekly$outbreak & !warned)
  fp <- sum(!weekly$outbreak & warned); tn <- sum(!weekly$outbreak & !warned)
  eval_df$sens[t] <- tp / (tp + fn)
  eval_df$spec[t] <- tn / (tn + fp)
  eval_df$ppv[t]  <- tp / (tp + fp + 1e-10)
  eval_df$npv[t]  <- tn / (tn + fn + 1e-10)
  eval_df$tp[t] <- tp; eval_df$fp[t] <- fp; eval_df$tn[t] <- tn; eval_df$fn[t] <- fn
}
eval_df$youden <- eval_df$sens + eval_df$spec - 1

cat("\n  Threshold evaluation results:\n")
print(eval_df %>% select(threshold, tp, fp, fn, tn, sens, spec, youden) %>%
        mutate(across(where(is.numeric), ~round(., 3))))

# Handle case where Youden optimisation fails
if (all(is.na(eval_df$youden))) {
  cat("\n  [!] No threshold achieved both sensitivity and specificity.\n")
  valid_sens <- which(!is.na(eval_df$sens) & eval_df$sens > 0)
  if (length(valid_sens) > 0) {
    opt_i <- valid_sens[which.max(eval_df$sens[valid_sens])]
  } else {
    opt_i <- 3
  }
  opt_thr <- thresholds[opt_i]
} else {
  opt_i   <- which.max(eval_df$youden)
  opt_thr <- eval_df$threshold[opt_i]
}

# Force tp/fp/fn/tn to 0 if NA (happens when no gene data exists for lookback)
eval_df$tp[opt_i] <- ifelse(is.na(eval_df$tp[opt_i]), 0, eval_df$tp[opt_i])
eval_df$fp[opt_i] <- ifelse(is.na(eval_df$fp[opt_i]), 0, eval_df$fp[opt_i])
eval_df$fn[opt_i] <- ifelse(is.na(eval_df$fn[opt_i]), 0, eval_df$fn[opt_i])
eval_df$tn[opt_i] <- ifelse(is.na(eval_df$tn[opt_i]), 0, eval_df$tn[opt_i])

cat(sprintf("Using gene threshold: %.2f copies/mL\n", opt_thr))
cat(sprintf("  -> TP=%d, FP=%d, FN=%d, TN=%d\n",
            eval_df$tp[opt_i], eval_df$fp[opt_i], eval_df$fn[opt_i], eval_df$tn[opt_i]))

# --- Evaluate Lagged Cases Baseline ---
cat("  Evaluating lagged cases baseline...\n")
lagged_warn <- c(FALSE, weekly$outbreak[-n_w])
bl_tp <- sum(weekly$outbreak & lagged_warn)
bl_fn <- sum(weekly$outbreak & !lagged_warn)
bl_fp <- sum(!weekly$outbreak & lagged_warn)
bl_tn <- sum(!weekly$outbreak & !lagged_warn)
bl_sens <- bl_tp / (bl_tp + bl_fn + 1e-10)
bl_spec <- bl_tn / (bl_tn + bl_fp + 1e-10)
bl_ppv <- bl_tp / (bl_tp + bl_fp + 1e-10)
bl_npv <- bl_tn / (bl_tn + bl_fn + 1e-10)
bl_youden <- bl_sens + bl_spec - 1

# --- Figure 23: Outbreak Detection Performance ---
cat("  Generating Figure 23...\n")
evl <- eval_df %>% 
  pivot_longer(cols = c(sens, spec, ppv, npv),
               names_to = "metric", values_to = "value")

p23a <- ggplot(evl %>% filter(metric %in% c("sens", "spec")),
               aes(x = threshold, y = value, color = metric)) +
  geom_line(linewidth = .9, na.rm = TRUE) + geom_point(size = 2, na.rm = TRUE) +
  scale_color_manual(values = c(sens = "#E74C3C", spec = "#2ECC71")) +
  labs(x = "Gene threshold", y = "Proportion") +
  theme_minimal(base_size = 10) + coord_cartesian(ylim = c(0, 1))

p23b <- ggplot(evl %>% filter(metric %in% c("ppv", "npv")),
               aes(x = threshold, y = value, color = metric)) +
  geom_line(linewidth = .9, na.rm = TRUE) + geom_point(size = 2, na.rm = TRUE) +
  scale_color_manual(values = c(ppv = "#E67E22", npv = "#3498DB")) +
  labs(x = "Gene threshold", y = "Proportion") +
  theme_minimal(base_size = 10) + coord_cartesian(ylim = c(0, 1))

p23c <- ggplot(eval_df %>% filter(!is.na(sens), !is.na(spec)),
               aes(x = 1 - spec, y = sens)) +
  geom_line(color = "#8E44AD", linewidth = .9) +
  geom_point(size = 2.5, color = "#8E44AD") +
  geom_point(data = eval_df[opt_i, ], size = 4, shape = 21,
             fill = "gold", stroke = 1.5, color = "#8E44AD") +
  geom_abline(slope = 1, linetype = "dashed", color = "grey50") +
  labs(
       x = "1 - Specificity", y = "Sensitivity") +
  theme_minimal(base_size = 10) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  annotate("text", x = 0.5, y = 0.5, size = 4, color = "grey40", hjust = 0.5,
           label = sprintf("Sens = %.0f%% | Spec = %.0f%%",
                           eval_df$sens[opt_i] * 100, eval_df$spec[opt_i] * 100))

# use geom_tile with visible borders and minimum alpha so 0-cells are visible
cm <- matrix(c(eval_df$tp[opt_i], eval_df$fn[opt_i],
               eval_df$fp[opt_i], eval_df$tn[opt_i]), 2, 2,
             dimnames = list(Actual = c("Outbreak", "No Outbreak"),
                             Predicted = c("Warned", "Not Warned")))
cm_df <- as.data.frame(as.table(cm))

p23d <- ggplot(cm_df, aes(x = Predicted, y = Actual, fill = Freq)) +
  geom_tile(color = "black", linewidth = 0.5) +
  geom_text(aes(label = Freq), size = 6, fontface = "bold",
            color = ifelse(cm_df$Freq == 0, "grey50", "white")) +
  scale_fill_gradient(low = "#FDDDAC", high = "#E74C3C",
                      limits = c(0, max(cm_df$Freq, na.rm = TRUE))) +
  
  theme_minimal(base_size = 10) + theme(legend.position = "none")

p23 <- grid.arrange(p23a, p23b, p23c, p23d, nrow = 2
                    )
ggsave(file.path(fig_dir, "fig23_early_warning.png"), p23, width = 12, height = 8, dpi = 600)
rm(p23, p23a, p23b, p23c, p23d, cm_df, evl); cat("  [OK] Figure 23 saved\n")

# --- Figure 24: Lead Time Distribution ---
cat("  Computing lead times...\n")
gene_high_opt <- weekly$mean_gene >= opt_thr
lead_data <- data.frame(wave = character(), lead = integer(), stringsAsFactors = FALSE)

for (i in seq_len(nrow(waves))) {
  wks <- which(weekly$wave == waves$wave[i] & weekly$outbreak)
  if (length(wks) == 0) {
    lead_data <- rbind(lead_data, data.frame(wave = waves$wave[i], lead = NA_integer_))
    next
  }
  onset <- wks[1]
  prior <- max(1, onset - 8)
  warns <- which(gene_high_opt[prior:(onset - 1)]) + prior - 1
  lt <- if (length(warns) > 0) onset - warns[1] else NA_integer_
  lead_data <- rbind(lead_data, data.frame(wave = waves$wave[i], lead = lt))
}

cat("  Lead time results:\n")
print(lead_data)

# Handle all-NA case properly
lead_data$lead_display <- ifelse(is.na(lead_data$lead), 0, lead_data$lead)
lead_data$lead_label  <- ifelse(is.na(lead_data$lead), "No signal", as.character(lead_data$lead))

# Ensure factor levels
lead_data$wave <- factor(lead_data$wave, levels = waves$wave)

# Set a minimum y-max even when all zeros
max_lead <- max(c(lead_data$lead_display, 4), na.rm = TRUE)  # Minimum 4 for visibility

p24 <- ggplot(lead_data, aes(x = wave, y = lead_display, fill = wave)) +
  geom_col(width = 0.6, show.legend = FALSE, alpha = 0.85) +
  geom_text(aes(label = lead_label),
            vjust = -0.2,
            size = 4, 
            color = "grey40") +
  scale_fill_manual(values = wcols_n) +
  labs(x = NULL, y = "Lead Time (weeks)") +
  theme_minimal(base_size = 11) +
  coord_cartesian(ylim = c(0, max_lead * 1.2)) +
  annotate("text", x = 2.5, y = max_lead * 0.5, size = 3.5, color = "grey30",
           label = "Note: No early warning signal detected at municipality level.\nGene threshold never exceeded before outbreak onsets.\nSpatial RK model (Script 5) provides subplace-level predictions.")

ggsave(file.path(fig_dir, "fig24_lead_time.png"), p24, width = 8, height = 5.5, dpi = 600)
rm(p24); cat("  [OK] Figure 24 saved\n")

# --- Figure 25: Dual-Axis Surveillance Timeline ---
cat("  Generating surveillance timeline...\n")
timeline <- weekly %>%
  mutate(cases_smooth = rollapplyr(total_cases, 3, mean, fill = NA),
         gene_warn = mean_gene >= opt_thr)

max_cases <- max(timeline$total_cases, na.rm = TRUE)
max_gene  <- max(timeline$mean_gene, na.rm = TRUE)
scale_f <- if (max_gene > 0 && !is.na(max_gene) && !is.na(max_cases)) {
  max_cases / max_gene
} else {
  1
}

warn_rects <- timeline %>% filter(gene_warn) %>%
  transmute(xmin = week_start, xmax = week_start + 7)
has_warnings <- nrow(warn_rects) > 0

wave_rects <- waves %>% mutate(ymin = -Inf, ymax = Inf)

p25 <- ggplot(timeline, aes(x = week_start)) +
  geom_rect(data = wave_rects, aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax),
            alpha = 0.08, fill = "grey60", inherit.aes = FALSE) +
  geom_col(aes(y = total_cases), fill = "#3498DB", alpha = 0.35, width = 7) +
  geom_line(aes(y = cases_smooth), color = "#2C3E50", linewidth = 0.7, na.rm = TRUE) +
  geom_line(aes(y = mean_gene * scale_f), color = "#E74C3C", linewidth = 1, na.rm = TRUE) +
  geom_hline(yintercept = case_p75, linetype = "dashed", color = "darkblue", linewidth = 0.5) +
  geom_vline(xintercept = waves$start, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  scale_y_continuous(name = "Weekly Cases") +
  scale_y_continuous(name = "Mean Gene Copies/mL", 
                     sec.axis = sec_axis(trans = ~ . / scale_f)) +
  labs(
       x = NULL) +
  theme_minimal(base_size = 11)

if (has_warnings) {
  p25 <- p25 + geom_rect(data = warn_rects, aes(xmin = xmin, xmax = xmax,
                                                ymin = -Inf, ymax = Inf), alpha = 0.07, fill = "red",
                         inherit.aes = FALSE)
}

for (i in seq_len(nrow(lead_data))) {
  if (!is.na(lead_data$lead[i])) {
    p25 <- p25 + annotate("text", x = waves$start[i] + 14, y = Inf,
                          label = sprintf("%d-wk lead", lead_data$lead[i]),
                          vjust = 1.5, hjust = 0, size = 3, color = "#C0392B", fontface = "italic")
  }
}
ggsave(file.path(fig_dir, "fig25_surveillance_timeline.png"), p25, width = 12, height = 6, dpi = 600)
rm(p25, warn_rects, wave_rects); cat("  [OK] Figure 25 saved\n")

# --- Figure 26: Comprehensive Early Warning Method Comparison ---
cat("  Generating Figure 26...\n")

ew_comparison <- data.frame(
  Method = c("Gene Threshold (Optimal)", "Lagged Cases (1-week Baseline)"),
  Sensitivity = c(eval_df$sens[opt_i], bl_sens),
  Specificity = c(eval_df$spec[opt_i], bl_spec),
  PPV = c(eval_df$ppv[opt_i], bl_ppv),
  NPV = c(eval_df$npv[opt_i], bl_npv),
  Youden_J = c(eval_df$youden[opt_i], bl_youden),
  TP = c(eval_df$tp[opt_i], bl_tp),
  FP = c(eval_df$fp[opt_i], bl_fp),
  FN = c(eval_df$fn[opt_i], bl_fn),
  TN = c(eval_df$tn[opt_i], bl_tn)
)

cat("\n--- EARLY WARNING METHOD COMPARISON ---\n")
print(as.data.frame(ew_comparison %>% mutate(across(where(is.numeric), ~round(., 3)))))

ew_long <- ew_comparison %>%
  select(Method, Sensitivity, Specificity, PPV, NPV, Youden_J) %>%
  pivot_longer(cols = -Method, names_to = "Metric", values_to = "Value")

p26 <- ggplot(ew_long, aes(x = Metric, y = Value, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), alpha = 0.85, width = 0.6) +
  geom_text(aes(label = sprintf("%.2f", Value)), 
            position = position_dodge(width = 0.7), vjust = -0.3, size = 3) +
  scale_fill_manual(values = c("Gene Threshold (Optimal)" = "#E74C3C", 
                               "Lagged Cases (1-week Baseline)" = "#3498DB")) +
  labs(
       x = NULL, y = "Proportion") +
  coord_cartesian(ylim = c(0, 1.15)) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.text.x = element_text(angle = 0))

ggsave(file.path(fig_dir, "fig26_early_warning_comparison.png"), p26, width = 9, height = 6, dpi = 600)
rm(p26); cat("  [OK] Figure 26 saved\n")

# ================================================================
# PART C: MODEL PREDICTION PERFORMANCE (Figure 27)
# ================================================================
cat("\n=== PART C: MODEL PREDICTION PERFORMANCE (Figure 27) ===\n")

all_regression <- bind_rows(
  rk_overall %>%
    select(Model, Mean_MAE, Mean_RMSE, Mean_sMAPE, Mean_R2) %>%
    rename(MAE = Mean_MAE, RMSE = Mean_RMSE, sMAPE = Mean_sMAPE, R2 = Mean_R2) %>%
    mutate(Validation = "Spatial 10-fold CV"),
  
  comp_results %>%
    select(Model, MAE, RMSE, sMAPE, R2) %>%
    mutate(Validation = "Temporal 80/20 Test")
)

p27a <- ggplot(all_regression, aes(x = reorder(Model, R2), y = R2, fill = Validation)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), alpha = 0.85) +
  geom_text(aes(label = sprintf("%.3f", R2)),
            position = position_dodge(width = 0.7), vjust = -0.3, size = 3) +
  scale_fill_manual(values = c("Spatial 10-fold CV" = "#2ECC71", 
                               "Temporal 80/20 Test" = "#3498DB")) +
  labs(
       x = NULL, y = expression(R^2)) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8),
        legend.position = "bottom")

p27b <- ggplot(all_regression, aes(x = reorder(Model, MAE), y = MAE, fill = Validation)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), alpha = 0.85) +
  geom_text(aes(label = sprintf("%.4f", MAE)),
            position = position_dodge(width = 0.7), vjust = -0.3, size = 3) +
  scale_fill_manual(values = c("Spatial 10-fold CV" = "#2ECC71", 
                               "Temporal 80/20 Test" = "#3498DB")) +
  labs(
       x = NULL, y = "Mean Absolute Error") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8),
        legend.position = "bottom")

p27 <- grid.arrange(p27a, p27b, nrow = 1
                )
ggsave(file.path(fig_dir, "fig27_model_prediction_comparison.png"), p27, 
       width = 14, height = 6, dpi = 600)
rm(p27, p27a, p27b); cat("  [OK] Figure 27 saved\n")

# ================================================================
# CONSOLE SUMMARY & SAVE RESULTS
# ================================================================
cat("\n============================================================\n")
cat("  COMPLETE EPIDEMIOLOGICAL & EARLY WARNING SUMMARY\n")
cat("============================================================\n")

cat("\n--- EPIDEMIOLOGICAL METRICS ---\n")
cat(sprintf("  Total study period: %s to %s (%d days)\n",
            min(daily$report_date), max(daily$report_date),
            as.numeric(max(daily$report_date) - min(daily$report_date))))
cat(sprintf("  Total cases: %s\n", format(sum(daily$total_cases, na.rm = TRUE), big.mark = ",")))
cat(sprintf("  Mean daily cases: %.1f\n", mean(daily$total_cases, na.rm = TRUE)))
cat(sprintf("  Peak daily cases: %s on %s\n", 
            format(max(daily$total_cases, na.rm = TRUE), big.mark = ","),
            daily$report_date[which.max(daily$total_cases)]))
cat(sprintf("  Mean Rt (overall): %.2f\n", mean(daily$rt, na.rm = TRUE)))
if ("cfr" %in% names(daily)) {
  cat(sprintf("  Mean CFR (30-day rolling): %.2f%%\n", mean(daily$cfr, na.rm = TRUE)))
}

cat("\n--- EARLY WARNING PERFORMANCE ---\n")
cat(sprintf("  Outbreak threshold: > %.0f weekly cases (75th percentile)\n", case_p75))
cat(sprintf("  Optimal gene threshold: %.2f copies/mL\n", opt_thr))
cat(sprintf("  Gene Method — Sens: %.3f | Spec: %.3f | PPV: %.3f | NPV: %.3f | Youden: %.3f\n",
            eval_df$sens[opt_i], eval_df$spec[opt_i], eval_df$ppv[opt_i], 
            eval_df$npv[opt_i], eval_df$youden[opt_i]))
cat(sprintf("  Lagged Cases  — Sens: %.3f | Spec: %.3f | PPV: %.3f | NPV: %.3f | Youden: %.3f\n",
            bl_sens, bl_spec, bl_ppv, bl_npv, bl_youden))
cat(sprintf("  Sensitivity improvement: %+.1f percentage points\n",
            (eval_df$sens[opt_i] - bl_sens) * 100))

cat("\n--- LEAD TIME BY WAVE ---\n")
for (i in seq_len(nrow(lead_data))) {
  cat(sprintf("  %s: %s\n", lead_data$wave[i],
              ifelse(is.na(lead_data$lead[i]), "No lead time detected",
                     sprintf("%d weeks", lead_data$lead[i]))))
}
cat(sprintf("  Median lead time: %s\n",
            ifelse(all(is.na(lead_data$lead)), "N/A",
                   sprintf("%.1f weeks", median(lead_data$lead, na.rm = TRUE)))))

cat("\n--- REGRESSION MODEL PERFORMANCE ---\n")
print(as.data.frame(all_regression %>% 
                      select(Validation, Model, R2, MAE, RMSE) %>%
                      mutate(across(where(is.numeric), ~round(., 4))) %>%
                      arrange(Validation, desc(R2))))

# --- Save all results ---
save(eval_df, wave_summary, lead_data, opt_thr, daily, weekly,
     ew_comparison, all_regression,
     file = file.path(output_dir, "epi_metrics_results.RData"))

write.csv(eval_df, file.path(res_dir, "early_warning_evaluation.csv"), row.names = FALSE)
write.csv(lead_data, file.path(res_dir, "lead_times.csv"), row.names = FALSE)
write.csv(wave_summary, file.path(res_dir, "wave_summary.csv"), row.names = FALSE)
write.csv(ew_comparison, file.path(res_dir, "early_warning_comparison.csv"), row.names = FALSE)
write.csv(all_regression, file.path(res_dir, "all_models_performance.csv"), row.names = FALSE)

cat("  SCRIPT 8 COMPLETE: ALL FIGURES AND RESULTS SAVED\n")
cat("  Figures 19-27 ->", fig_dir, "\n")
cat("  CSV tables     ->", res_dir, "\n")
cat("  RData          -> epi_metrics_results.RData\n")