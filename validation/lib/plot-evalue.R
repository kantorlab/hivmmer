library(dplyr)
library(ggplot2)
library(readr)

args <- commandArgs(trailingOnly=TRUE)

csv_file <- args[1]
title    <- args[2]
out_file <- args[3]

csv <- read_csv(csv_file, col_types="ccddddc") %>%
       mutate(Evalue=log(DBLength*Length/(2^Score)),
              Dataset=case_when(Dataset == "reference" & Gene == Read ~ "matching reference",
                                Dataset == "reference"                ~ "mis-matching reference",
                                TRUE                                  ~ Dataset))

thresholds <- csv %>%
              filter(Dataset == "UniVec" | Dataset == "random" | Dataset == "mis-matching reference") %>%
              group_by(Gene, Dataset) %>%
              summarise(Evalue=quantile(Evalue, 0.001))
print(thresholds %>% group_by(Gene) %>% summarise(Evalue=min(Evalue)))

png(out_file, width=10, height=10, units="in", res=150)
ggplot(csv, aes(x=Evalue, color=Dataset)) +
  geom_step(stat="ecdf", size=0.25) +
  geom_vline(data=thresholds, aes(xintercept=Evalue, color=Dataset), size=0.25) +
  facet_grid(Gene ~ .) +
  labs(title=paste("E-value distributions for hmmscan alignments to", title),
       x="log(E-value)", y="Cumulative Density") +
  ylim(0.001, 1) +
  theme_light() +
  theme(legend.position="top")
dev.off()

