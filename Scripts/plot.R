test_df_all = fread("~/tre1_eqtl/alleqtls_analysis/finalsets/testset_07.30/sv_comparisonresults_all_test.csv")
control_df_all = fread("~/tre1_eqtl/alleqtls_analysis/analysis_08.25/control_selection/control_snps_quantiles/sv_comparisonresults_all_control.csv")
test_df_all$group = "test"
control_df_all$group ="control"


combined_df = rbind(test_df_all, control_df_all)

combined_df <- combined_df %>%
  mutate(sv_category = ifelse(sv_diff_type == "same_svs", "same_svs", "partiallydifferent_svs"))

counts <- combined_df %>%
  count(group, sv_category) %>%
  tidyr::pivot_wider(names_from = group, values_from = n, values_fill = 0)

contingency_table
contingency_table2

# Make contingency table
contingency_table <- table(combined_df$group, combined_df$sv_category)

# Fisher's Exact Test
fisher_test_result <- fisher.test(contingency_table)

print(fisher_test_result)


test_df_del = fread("~/tre1_eqtl/alleqtls_analysis/testset_07.30/SV_overlap_test_del/sv_comparisonresults_del_test.csv")
control_df_del = fread("~/tre1_eqtl/alleqtls_analysis/noneandtrans_control/SV_overlap_control_del/sv_comparisonresults_delonly_control.csv")

test_df_del$group = "test"
control_df_del$group = "control"

combined_df2 = rbind(test_df_del, control_df_del)

combined_df2 <- combined_df2 %>%
  mutate(sv_category = ifelse(sv_diff_type == "same_svs", "same_svs", "partiallydifferent_svs"))

counts2 <- combined_df2 %>%
  count(group, sv_category) %>%
  tidyr::pivot_wider(names_from = group, values_from = n, values_fill = 0)

contingency_table2 <- table(combined_df2$group, combined_df2$sv_category)

fisher_test_result2 <- fisher.test(contingency_table2)

print(fisher_test_result2)


library(data.table)
library(dplyr)
library(ggplot2)
library(ggtext)


library(RColorBrewer)


display.brewer.all()
spectral_colors <- brewer.pal(n = 11, name = "Spectral")
custom_colors <- c("control" = spectral_colors[4],  # reddish
                   "test" = spectral_colors[9])    

combined_df <- as.data.frame(combined_df)
combined_df$sv_category <- factor(combined_df$sv_category, levels = c("same_svs", "partiallydifferent_svs"))
total_counts1 <- combined_df %>% count(group) %>% rename(total = n)
total_counts1 <- combined_df %>% count(group, name = "total")


pval1 <- signif(fisher.test(contingency_table)$p.value, 3)

label_df1 <- combined_df %>%
  count(group, sv_category) %>%
  left_join(total_counts1, by = "group") %>%
  mutate(label = paste0(n, "/", total))

#plot for counts

ggplot(label_df1, aes(x = sv_category, y = n, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = label), 
            position = position_dodge(width = 0.9), 
            vjust = -0.3, size = 5, fontface = "bold") +
  scale_fill_manual(values = custom_colors) +
  labs(x = "SV Difference Category", y = "Count", fill = "Group",
       title = "SV Overlap: Test vs Control",
       subtitle = paste0("Fisher's Exact Test p = ", pval1)) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_markdown(),
    axis.text = element_text(color = "black"),
    legend.position = "top"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))


library(dplyr)
library(ggplot2)
library(scales)   # for percent_format()

##plot for proportion
theme_set(theme_classic(base_size = 19))
theme_update(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"),
             plot.title = element_text(hjust = 0.5), plot.subtitle=element_text(hjust=0.5))

# Build labels with counts and percent, and add proportion
label_df1 <- combined_df %>%
  count(group, sv_category) %>%
  left_join(total_counts1, by = "group") %>%
  mutate(prop  = n / total,
         label = paste0(n, "/", total, " (", percent(prop, accuracy = 0.1), ")"))

# Plot proportions, keep count labels
ggplot(label_df1, aes(x = sv_category, y = prop, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
  geom_text(aes(label = label),
            position = position_dodge(width = 0.9),
            vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = custom_colors) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.15))) +
  labs(x = "SV Difference Category",
       y = "Proportion within group",
       fill = "Group",
       title = "SV Overlap: Test vs Control",
       subtitle = paste0("Fisher's Exact Test p = ", pval1))



contingency_table


spectral_colors

combined_df2$sv_category <- factor(combined_df2$sv_category, levels = c("same_svs", "partiallydifferent_svs"))
total_counts2 <-  combined_df2 %>% count(group, name = "total")


pval2 <- signif(fisher.test(contingency_table2)$p.value, 3)

label_df2 <- combined_df2 %>%
  count(group, sv_category) %>%
  left_join(total_counts2, by = "group") %>%
  mutate(label = paste0(n, "/", total))

ggplot(label_df2, aes(x = sv_category, y = n, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = label), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 4) +
  scale_fill_manual(values = custom_colors) +
  labs(x = "SV Difference Category", y = "Count", fill = "Group",
       title = "SV Overlap (DEL): Test vs Control",
       subtitle = paste0("Fisher's Exact Test p = ", pval2)) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_markdown(),
    axis.text = element_text(color = "black"),
    legend.position = "top"
  )














