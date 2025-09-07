library(ggplot2)
library(scales)
library(dplyr)

# Counts
total_test <- 7859
perfect_test <- 535
nonperfect_test <- total_test - perfect_test

total_control <- 2942
perfect_control <- 56
nonperfect_control <- total_control - perfect_control

# Fisher’s test
cont_table <- matrix(c(perfect_test, nonperfect_test,
                       perfect_control, nonperfect_control),
                     nrow = 2, byrow = TRUE)
cont_table

rownames(cont_table) <- c("Test", "Control")
colnames(cont_table) <- c("Perfect", "NonPerfect")

fisher_result <- fisher.test(cont_table)
pval1 <- signif(fisher_result$p.value, 3)
pval1


label_df1 <- data.frame(
  group = c("Test", "Control"),
  perfect = c(perfect_test, perfect_control),
  total   = c(total_test, total_control)
) %>%
  mutate(prop = perfect / total,
         label = paste0(perfect, "/", total, " (", percent(prop, accuracy = 0.1), ")"))

# Plot
theme_set(theme_classic(base_size = 19))
theme_update(axis.text.x = element_text(colour = "black"),
             axis.text.y = element_text(colour = "black"),
             plot.title = element_text(hjust = 0.5),
             plot.subtitle = element_text(hjust = 0.5))

ggplot(label_df1, aes(x = group, y = prop, fill = group)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = label), vjust = -0.3, size = 4) +
  scale_fill_manual(values = c("Test" = "#1f78b4", "Control" = "#33a02c")) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Group",
       y = "Proportion of Perfect Matches",
       fill = "Group",
       title = "SV Overlap: Test vs Control",
       subtitle = paste0("Fisher's Exact Test p = ", pval1))
