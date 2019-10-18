library(dplyr)
library(ggplot2)
library(readr)

args <- commandArgs(trailingOnly=TRUE)

csv_file <- args[1]
out_file <- args[2]

csv <- read_csv(csv_file, col_types="ccdc") %>%
       filter(Dataset == "reference") %>%
       mutate(Within=(Gene == Read))

print(csv %>%
      filter(Within == FALSE) %>%
      group_by(Gene) %>%
      summarise(Evalue=quantile(Evalue, 0.001)))

png(out_file, width=10, height=6, units="in", res=150)
ggplot(csv, aes(x=Within, y=log(Evalue))) +
  geom_boxplot() +
  facet_grid(. ~ Gene) +
  labs(title="E-value distributions for hmmscan alignments of gene matching vs. non-matching reference sequences",
       x="Gene matching sequences", y="log(E-value)") +
  theme_light()
dev.off()

