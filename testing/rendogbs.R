library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
pt_files <- list.files("input/Populus_tremula/", pattern = "*\\.csv")
pt_files <- paste("input/Populus_tremula/", pt_files, sep = "")
pt_datasets <- lapply(X = pt_files, FUN = read.table, skip = 6, header = T, sep = ",")
# pt_df <- do.call(what = rbind, args = pt_datasets)
pt_df <- bind_rows(pt_datasets, .id = "source")
ggplot(data = pt_df, mapping = aes(x = factor(gc_bin_pct), y = matched_n)) +
  geom_col() +
  facet_wrap(.~source) +
  theme_classic()
pt_df_sum <- pt_df |> 
              group_by(gc_bin_pct) |> 
              summarise(min_matched_n = min(matched_n),
                        max_matched_n = max(matched_n),
                        mean_matched_n = mean(matched_n),
                        median_matched_n = median(matched_n),
                        sd_matched_n = sd(matched_n),
                        min_matched_pct_of_all_matched = min(matched_pct_of_all_matched),
                        max_matched_pct_of_all_matched = max(matched_pct_of_all_matched),
                        mean_matched_pct_of_all_matched = mean(matched_pct_of_all_matched),
                        median_matched_pct_of_all_matched = median(matched_pct_of_all_matched),
                        sd_matched_pct_of_all_matched = sd(matched_pct_of_all_matched),
                        min_unmatched_n = min(unmatched_n),
                        max_unmatched_n = max(unmatched_n),
                        mean_unmatched_n = mean(unmatched_n),
                        median_unmatched_n = median(unmatched_n),
                        sd_unmatched_n = sd(unmatched_n),
                        min_unmatched_pct_of_all_unmatched = min(unmatched_pct_of_all_unmatched),
                        max_unmatched_pct_of_all_unmatched = max(unmatched_pct_of_all_unmatched),
                        mean_unmatched_pct_of_all_unmatched = mean(unmatched_pct_of_all_unmatched),
                        median_unmatched_pct_of_all_unmatched = median(unmatched_pct_of_all_unmatched),
                        sd_unmatched_pct_of_all_unmatched = sd(unmatched_pct_of_all_unmatched),
                        min_total_n = min(total_n),
                        max_total_n = max(total_n),
                        mean_total_n = mean(total_n),
                        median_total_n = median(total_n),
                        sd_total_n = sd(total_n),
                        min_recovery_rate_in_bin = min(recovery_rate_in_bin),
                        max_recovery_rate_in_bin = max(recovery_rate_in_bin),
                        mean_recovery_rate_in_bin = mean(recovery_rate_in_bin),
                        median_recovery_rate_in_bin = median(recovery_rate_in_bin),
                        sd_recovery_rate_in_bin = sd(recovery_rate_in_bin)
              )
pt_df_sum_long <- pt_df_sum |> 
                      pivot_longer(cols = 2:31, names_to = "variable") |> 
                      arrange(variable, gc_bin_pct)

parameter <- c("matched_n", "matched_pct_of_all_matched", "recovery_rate_in_bin",
               "total_n", "unmatched_n", "unmatched_pct_of_all_unmatched")                     
parameters <- rep(parameter, each = 20)
parameters <- rep(parameters, times = 5)
pt_df_sum_long$parameters <- parameters
pt_df_sum_long$test <- str_extract(string = pt_df_sum_long$variable, pattern = "^[[:alpha:]]{2,6}")
pt_df_sum_long$variable <- NULL
pt_df_sum_wide <- pt_df_sum_long |> 
                      mutate(gc_bin_pct = factor(gc_bin_pct, levels = c("0-5", "5-10", "10-15", "15-20", "20-25", "25-30", "30-35",
                                                                        "35-40", "40-45", "45-50", "50-55", "55-60","60-65", "65-70",
                                                                        "70-75","75-80", "80-85", "85-90", "90-95", "95-100"),
                                                 labels = c("0-5", "5-10", "10-15", "15-20", "20-25", "25-30", "30-35",
                                                            "35-40", "40-45", "45-50", "50-55", "55-60","60-65", "65-70", "70-75",
                                                            "75-80", "80-85", "85-90", "90-95", "95-100")
                                                 ),
                             parameter = factor(parameters, levels = c("matched_n", "matched_pct_of_all_matched",
                                                                       "unmatched_n", "unmatched_pct_of_all_unmatched",
                                                                       "total_n", "recovery_rate_in_bin"))
                             ) |> 
                      pivot_wider(names_from = "test")
ggplot(pt_df_sum_wide, aes(x = gc_bin_pct, y = mean)) +
  geom_point(shape = 21, fill = "gold2", colour = "black", size = 2, alpha = 0.4) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0, colour = "purple") +
  facet_wrap(~parameter, scales = "free_y", nrow = 3) +
  labs(y = "mean +/- sd") + 
  theme_classic() + 
  theme(axis.text.x.bottom = element_text(angle = 65,hjust = 1))
ggsave(filename = "graphic/preliminary_test.jpg", height = 180, width = 140, units = "mm")  
