#Calculate Fst, Dxy, and nucleotide diveristy between D2 and D3 populations with PopGenome
library(ape)
library(tidyverse)
library(PopGenome)
library(dplyr)
library(ggplot2)

#Load the whole-genome alignment
laspinema <- readData("./alignment.fa", format = "fasta", include.unknown = TRUE, FAST = TRUE)

#Set populations (make text file with two columns, one with sample names and second one with population)
laspinema_info <- read_delim("./laspinema_info.txt", delim = "\t")
populations <- split(laspinema_info$ind, laspinema_info$pop)
laspinema <- set.populations(laspinema, populations, diploid = F)

#Preparing the dataset to calculate measurements in 50kb window size with a step of 12.5kb
laspinema_sw_50 <- sliding.window.transform(laspinema, width = 50000, jump = 12500, type = 2)
laspinema_sw_50 <- diversity.stats(laspinema_sw_50, pi = TRUE)
laspinema_nuc_div_sw_50 <- laspinema_sw_50@nuc.diversity.within
laspinema_nuc_div_sw_50 <- laspinema_nuc_div_sw_50/50000
position <- seq(from = 1, to = 7804270, by = 12500) + 50000
position <- seq(from = 1, to = 7804270, by = 12500)
window_stop <- position + 50000
position <- position[which(window_stop < 7804270)]
window_stop <- window_stop[which(window_stop < 7804270)]
windows <- data.frame(start=position, stop=window_stop, mid= position + (window_stop - position)/2)
laspinema_nd_sw_50 <- data.frame(windows, laspinema_nuc_div_sw_50, row.names = NULL)
laspinema_nd_sw_50 <- as_tibble(laspinema_nd_sw_50)

#Calculating Fst
laspinema_sw_50_fst <- F_ST.stats(laspinema_sw_50, mode = "nucleotide")
nd_50 <- laspinema_sw_50_fst@nuc.diversity.within/50000
pops <- c("D2", "D3")
colnames(nd_50) <- paste0(pops, "_pi")
Fst_value_50 <- t(laspinema_sw_50_fst@nuc.F_ST.pairwise)
dxy_50 <- get.diversity(laspinema_sw_50_fst, between = T) [[2]]/50000
xx <- colnames(Fst_value_50)
xx <- sub("pop1", "D2", xx)
xx <- sub("pop1", pops[1], xx)
xx <- sub("pop2", pops[2], xx)
xx <- sub("/", "_", xx)
colnames(Fst_value_50) <- paste0(xx, "_fst")
colnames(dxy_50) <- paste0(xx, "_dxy")

#Data frame with all the values; Fst, Dxy, piD2, piD1
laspinema_data_50 <- as_tibble(data.frame(windows, nd_50, Fst_value_50, dxy_50))
laspinema_data_50 %>% select(contains("pi")) %>% summarise_all(mean)
write.csv(laspinema_data_50, "./new.csv")

#remove the row with the region of inflated Fst (2862501 mid position)
#load new file without the region
new <- read.csv("./new.csv", row.names=1)

#plot all the values
combined <- new %>% select(mid, D2_D3_dxy, D2_D3_fst)
combined_g <- gather(combined, -mid, key="stat", value="value")
combined <- ggplot(combined_g, aes(mid/10^6, value, color=stat)) + geom_line(size=1)
combined <- combined + facet_grid(stat~., scales = "free_y")
combined <- combined + xlab("Position (Mb)")
combined + theme_bw() + theme(legend.position = "none", text = element_text(size=20)) + geom_point(size=1.5)

#Box plots pi
pi_g <- new %>% select(contains("pi")) %>% gather(key = "species", value = "pi")
ggplot(pi_g, aes(species, pi)) + geom_boxplot() + theme_light() + xlab(NULL)

#Wilcox test
D2pi <- pull(new, D2_pi)
D3pi <- pull(new, D3_pi)
wilcox.test(D3pi, D2pi, exact = FALSE, alternative = "less", conf.int = 0.95)

#########Wilcoxon rank sum test with continuity correction
######data:  D3pi and D2pi
####W = 10648, p-value < 2.2e-16
