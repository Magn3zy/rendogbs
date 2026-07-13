library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(patchwork)
library(ggtext)

# extraction nr. 1
path_to_all_samples <- list.files("input", recursive = TRUE, pattern = "*\\.csv", full.names = TRUE)
list_names <- str_remove(string = path_to_all_samples, pattern = "input\\/") |> 
  str_remove(pattern = "\\.csv") |> 
  str_replace(pattern = "\\/", replacement = "\\-")
# test <- readLines(con = "input/Anas_platyrhynchos/SRR24837080.csv")
# grep(pattern = "^gc_bin_pct", x = test)
# stringr::str_which(test, pattern = "^gc_bin_pct") # grep ekvivalent
# file_upload(x = "input/Anas_platyrhynchos/SRR24837080.csv")
metadata_extraction<- function(x){
  upload <- readLines(x, n = 7)
  metadata_names <- stringr::str_extract(string = upload, pattern = "^[[:alpha:]_]+")
  metadata_values <- stringr::str_extract(string = upload, pattern = "[[:digit:]]+")
  meta_dataset <- data.frame(meta_names = metadata_names, meta_values = metadata_values)
  return(meta_dataset)
}
# test <- readLines("input/Anas_platyrhynchos/SRR24836806.csv", n = 7)
# vysledovka <- str_split_fixed(test, pattern = ",", n = 2)
meta_list <- lapply(X = path_to_all_samples, FUN = metadata_extraction)
names(meta_list) <- list_names
meta_df <- bind_rows(meta_list, .id = "source")
meta_df_sum <- meta_df |> 
  separate(col = source, into = c("species", "project_id"),sep = "-") |> 
  filter(meta_names %in% c("Predicted_fragments", "Matched", "Total_reads")) |> 
  mutate(meta_values = as.numeric(meta_values)) |> 
  pivot_wider(names_from = meta_names, values_from = meta_values) |> 
  mutate(tr_vs_predicted = round(Total_reads/Predicted_fragments, digits = 1), 
         pct_matched = round((Matched / Predicted_fragments)*100, digits = 2)) |> 
  select(species, tr_vs_predicted, pct_matched) |> 
  group_by(species) |> 
  summarise(mean_tr_vs_predicted = mean(tr_vs_predicted),
            sd_tr_vs_predicted = sd(tr_vs_predicted),
            min_tr_vs_predicted = min(tr_vs_predicted),
            max_tr_vs_predicted = max(tr_vs_predicted),
            mean_pct_matched = mean(pct_matched),
            sd_pct_matched = sd(pct_matched),
            min_pct_matched = min(pct_matched),
            max_pct_matched = max(pct_matched)
  ) |> 
  mutate(methylation_sensitive = c(F, T, F, F, F, T, T, T, F, F, F, T, T, T, F, T, F, T, T, F)) |> 
  arrange(methylation_sensitive, species) |> 
  mutate(species_num = rep(1:10, times = 2),
         species_legend = paste(species_num, ": ", species, sep = "")) 
sl_nonsensitive <- tibble(lab = c("1: *Anas platyrhynchos*<br>
                                      2: *Camellia sinensis assamica*<br>
                                      3: *Coffea arabica*<br>
                                      4: *Crassostrea virginica*<br>
                                      5: *Oncorhynchus mykiss*<br>
                                      6: *Oryza sativa*<br>
                                      7: *Populus tremula*<br>
                                      8: *Quercus rubra*<br>
                                      9: *Scatophagus argus*<br>
                                      10: *Sparus aurata*"
                                  ) 
                                    )
sl_sensitive <- tibble(lab = c("1: *Brassica napus*<br>
                                2: *Fragaria \u00D7 ananassa*<br>
                                3: *Halyomorpha halys*<br>
                                4: *Labeo rohita*<br>
                                5: *Prunus persica*<br>
                                6: *Quercus cerris*<br>
                                7: *Quercus ilex*<br>
                                8: *Salmo trutta*<br>
                                9: *Sesamum indicum*<br>
                               10: *Solanum lycopersicum*")
                       )


p1 <-  ggplot(data = meta_df_sum [meta_df_sum$methylation_sensitive == 0,], 
              mapping = aes(x = mean_pct_matched,
                            y = mean_tr_vs_predicted)) +
  geom_errorbar(aes(xmin = min_pct_matched, xmax = max_pct_matched,colour = "Range"),
                alpha = 0.2,
                linewidth = 0.5,
                linetype = 1,
                width = 30
  ) +
  geom_errorbar(aes(ymin = min_tr_vs_predicted, ymax = max_tr_vs_predicted, colour = "Range"),
                alpha = 0.2,
                linewidth = 0.5,
                linetype = 1,
                width = 0.5
  ) +
  geom_errorbar(aes(xmin = mean_pct_matched - sd_pct_matched, xmax = mean_pct_matched + sd_pct_matched, colour = "SD"), 
                linewidth = 0.5,
                width = 30
  ) +
  geom_errorbar(aes(ymin = mean_tr_vs_predicted - sd_tr_vs_predicted, ymax = mean_tr_vs_predicted + sd_tr_vs_predicted, colour = "SD"), 
                linewidth = 0.5,
                width = 0.5) +
  geom_point(shape = 21, colour = "black", fill = "#00796B", size = 3) +
  scale_colour_manual(name = "", values = c("Range" = "purple", "SD" = "#4DB6AC"))+
  # geom_point(aes(colour = reorder(species_legend, species_num))) +
  geom_text(aes(label = species_num), size = 2, colour = "ivory") +
  labs(x = "Matched", y = "Total reads / Predicted fragments ratio", title = "Methylation non-sensitive") +
  # annotate("text", label = meta_df_sum$species_legend [meta_df_sum$methylation_sensitive == FALSE], 
  #          x = 50, y = seq(from = 1000, to = 1900, by = 100))+
  geom_richtext(data = sl_nonsensitive, aes(x = 26, y = 1800, label = lab), 
                hjust = 0, 
                size = 3,
                label.colour = "gray")+
  scale_x_continuous(breaks = seq(from = 30, to = 100, by = 10), 
                     labels = paste(seq(from = 30, to = 100, by = 10), "%", sep = ""))+
  scale_y_continuous(breaks = seq(from = 0, to = 2500, by = 500), 
                     labels = paste(seq(from = 0, to = 2500, by = 500), "\u00D7", sep = ""))+
  theme_bw() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.1,0.4),
        legend.background = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10),
        panel.grid.major = element_line(colour = "gray", linetype = 3, linewidth = 0.3),
        # panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 6, face = "bold"),
        axis.title = element_text(size = 8)
  )
# ggsave(filename = "graphic/test_non-sensitive.jpg", units = "mm", height = 120, width = 150)

# methylation sensitive
# palette_yes = {
#   1: "#FDD9D2",
#   2: "#F28E7F",
#   3: "#D95F4A"
# }
# 
p2 <- ggplot(data = meta_df_sum [meta_df_sum$methylation_sensitive == 1,], 
             mapping = aes(x = mean_pct_matched,
                           y = mean_tr_vs_predicted)) +
  geom_errorbar(aes(xmin = min_pct_matched, xmax = max_pct_matched,colour = "Range"),
                alpha = 0.2,
                linewidth = 0.5,
                linetype = 1,
                width = 30
  ) +
  geom_errorbar(aes(ymin = min_tr_vs_predicted, ymax = max_tr_vs_predicted, colour = "Range"),
                alpha = 0.2,
                linewidth = 0.5,
                linetype = 1,
                width = 0.5
  ) +
  geom_errorbar(aes(xmin = mean_pct_matched - sd_pct_matched, xmax = mean_pct_matched + sd_pct_matched, colour = "SD"), 
                linewidth = 0.5,
                width = 30
  ) +
  geom_errorbar(aes(ymin = mean_tr_vs_predicted - sd_tr_vs_predicted, ymax = mean_tr_vs_predicted + sd_tr_vs_predicted, colour = "SD"), 
                linewidth = 0.5,
                width = 0.5) +
  geom_point(shape = 21, colour = "black", fill = "#D95F4A", size = 3) +
  scale_colour_manual(name = "", values = c("Range" = "purple", "SD" = "#F28E7F"))+
  # geom_point(aes(colour = reorder(species_legend, species_num))) +
  geom_text(aes(label = species_num), size = 2, colour = "ivory") +
  labs(x = "Matched", y = "Total reads / Predicted fragments ratio", title = "Methylation sensitive") +
  geom_richtext(data = sl_sensitive, aes(x = 16, y = 1125, label = lab),
                hjust = 0,
                size = 3,
                label.colour = "gray")+
  scale_x_continuous(breaks = seq(from = 20, to = 100, by = 10),
                     labels = paste(seq(from = 20, to = 100, by = 10), "%", sep = ""))+
  scale_y_continuous(breaks = seq(from = 0, to = 1400, by = 200),
                     labels = paste(seq(from = 0, to = 1400, by = 200), "\u00D7", sep = ""))+
  theme_bw() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.1,0.39),
        legend.background = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10),
        panel.grid.major = element_line(colour = "gray", linetype = 3, linewidth = 0.3),
        # panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 6, face = "bold"),
        axis.title = element_text(size = 8)
  )

# ggsave(filename = "graphic/test_sensitive.jpg", units = "mm", height = 120, width = 150)
p1 + p2 + plot_annotation(tag_levels = list(c("(A)", "(B)"))) & theme(plot.tag = element_text(size = 8))
ggsave(filename = "graphic/multi_predicted_vs_matched.jpg", units = "mm", height = 120, width = 300)


# extraction nr. 2
path_to_all_samples <- list.files("input", recursive = TRUE, pattern = "*\\.csv", full.names = TRUE)
list_names <- str_remove(string = path_to_all_samples, pattern = "input\\/") |> 
                str_remove(pattern = "\\.csv") |> 
                str_replace(pattern = "\\/", replacement = "\\-")
# test <- readLines(con = "input/Anas_platyrhynchos/SRR24837080.csv")
# grep(pattern = "^gc_bin_pct", x = test)
# stringr::str_which(test, pattern = "^gc_bin_pct") # grep ekvivalent
# file_upload(x = "input/Anas_platyrhynchos/SRR24837080.csv")
file_upload <- function(x){
                  upload <- readLines(x)
                  metadata_skip <- grep(pattern = "^gc_bin_pct", x = upload)
                  dataset <- read.table(file = x, skip = metadata_skip - 1, sep = ",", header = T, fill = TRUE)
                  return(dataset)
                  }
list_samples <- lapply(X = path_to_all_samples, FUN = file_upload)
names(list_samples) <- list_names
df_samples <- bind_rows(list_samples, .id = "source")
df_samples$matched_pct_of_all_unmatched <- NULL
df_samples <- df_samples |> 
                  separate(col = source, into = c("species", "project_id"),sep = "-")
                  
df_samples_sum <- df_samples |> 
  mutate(gc_bin_pct = factor(gc_bin_pct, levels = c("0-5", "5-10", "10-15", "15-20", "20-25", "25-30", "30-35",
                                                     "35-40", "40-45", "45-50", "50-55", "55-60","60-65", "65-70",
                                                     "70-75","75-80", "80-85", "85-90", "90-95", "95-100"))) |>
  group_by(species, gc_bin_pct) |>
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
df_samples_sum_long <- df_samples_sum |> 
  pivot_longer(cols = 3:32, names_to = "variable") |> 
  arrange(variable, gc_bin_pct)

parameter_name <- c("matched_n", "matched_pct_of_all_matched", "recovery_rate_in_bin",
               "total_n", "unmatched_n", "unmatched_pct_of_all_unmatched")                     
parameter <- rep(parameter_name, each = 400) # nr samples per variable
parameter <- rep(parameter, times = 5) # deskriptivnich parametru - min, max, mean, median, sd
df_samples_sum_long$parameter <- parameter
df_samples_sum_long$descriptor <- str_extract(string = df_samples_sum_long$variable, pattern = "^[[:alpha:]]{2,6}")
df_samples_sum_long$variable <- NULL
df_samples_sum_wide <- df_samples_sum_long |> 
                        mutate(parameter = factor(parameter, levels = c("matched_n", "matched_pct_of_all_matched",
                                            "unmatched_n", "unmatched_pct_of_all_unmatched",
                                            "total_n", "recovery_rate_in_bin")
                                            )
                               ) |> 
  pivot_wider(names_from = "descriptor") 
df_samples_sum_wide <- df_samples_sum_wide |> 
                        arrange(parameter, gc_bin_pct, species)

df_samples_sum_wide <- df_samples_sum_wide |> 
  mutate(gc_bin_pct_num = as.numeric(gc_bin_pct)) 

ggplot(data = df_samples_sum_wide, mapping = aes(x = gc_bin_pct, y = mean)) + 
geom_point(shape = 21, fill = "gold2", colour = "black", size = 2, alpha = 0.4) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0, colour = "purple") +
  facet_grid(parameter~species, scales = "free_y") +
  labs(y = "mean +/- sd") + 
  theme_classic() + 
  theme(axis.text.x.bottom = element_text(angle = 65,hjust = 1))
ggsave(filename = "graphic/multi_preliminary_test.jpg", height = 360, width = 1080, units = "mm")  



# meta_df_sum |> select(species, methylation_sensitive) |> arrange(methylation_sensitive, species)
# pt_files <- list.files("input/Populus_tremula/", pattern = "*\\.csv")
# pt_files <- paste("input/Populus_tremula/", pt_files, sep = "")
# pt_datasets <- lapply(X = pt_files, FUN = read.table, skip = 6, header = T, sep = ",")
# pt_df <- do.call(what = rbind, args = pt_datasets)
pt_df <- bind_rows(pt_datasets, .id = "source")
# ggplot(data = pt_df, mapping = aes(x = factor(gc_bin_pct), y = matched_n)) +
#   geom_col() +
#   facet_wrap(.~source) +  
#   theme_classic()
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
 pt_df_sum_wide$parameters <- NULL                     
pt_df_sum_wide <- pt_df_sum_wide |> 
                    mutate(gc_bin_pct_num = as.numeric(gc_bin_pct)) 
  ggplot
  geom_point(shape = 21, fill = "gold2", colour = "black", size = 2, alpha = 0.4) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0, colour = "purple") +
  facet_wrap(~parameter, scales = "free_y", nrow = 3) +
  labs(y = "mean +/- sd") + 
  theme_classic() + 
  theme(axis.text.x.bottom = element_text(angle = 65,hjust = 1))
ggsave(filename = "graphic/preliminary_test.jpg", height = 180, width = 140, units = "mm")  

# 
getwd()
BiocManager::install("Biostrings")
library(Biostrings)
phylo_path <- list.files(path = "input/18S_rRNA_fasta", full.names = TRUE)
phylo_list <- lapply(X = phylo_path, FUN = readDNAStringSet)
phylo_merge <- do.call(c, args = phylo_list)
writeXStringSet(phylo_merge, "output/phylo_merge.fasta")
