library(dplyr)
library(ggplot2)
library(readr)

args <- commandArgs(trailingOnly=TRUE)

csv_file <- args[1]
out_file <- args[2]

csv <- read_csv(csv_file, col_types="ccdc") %>%
       mutate(Evalue=pmax(-20, log(Evalue)),
	      Within=(Gene == Read))

thresholds <- csv %>%
              filter(Dataset == "UniVec" | Dataset == "random" | (Dataset == "perfect" & Within == FALSE)) %>%
              group_by(Gene, Dataset) %>%
              summarise(Evalue=quantile(Evalue, 0.001))
print(thresholds %>% group_by(Gene) %>% summarise(Evalue=min(exp(Evalue))))

pdf(out_file, width=6, height=10)
ggplot(csv, aes(x=Evalue, color=Dataset)) +
  geom_step(stat="ecdf", size=0.25) +
  geom_vline(data=thresholds, aes(xintercept=Evalue, color=Dataset), size=0.25) +
  facet_grid(Gene ~ .) +
  labs(x="HMMER log(E-value)", y="Cumulative Density") +
  xlim(-20, 0) +
  ylim(0.001, 0.999) +
  theme_light() +
  theme(legend.position="top")
dev.off()

