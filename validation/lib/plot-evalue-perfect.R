library(dplyr)
library(ggplot2)
library(readr)

args <- commandArgs(trailingOnly=TRUE)

csv_file <- args[1]
out_file <- args[2]

csv <- read_csv(csv_file, col_types="ccdc") %>%
       filter(Dataset == "perfect") %>%
       mutate(Within=(Gene == Read))

print(csv %>%
      filter(Within == FALSE) %>%
      group_by(Gene) %>%
      summarise(Evalue=quantile(Evalue, 0.001)))

pdf(out_file, width=10, height=6)
ggplot(csv, aes(x=Within, y=log(Evalue))) +
  geom_boxplot() +
  facet_grid(. ~ Gene) +
  labs(x="Reads within gene", y="HMMER log(E-value)") +
  theme_light()
dev.off()

