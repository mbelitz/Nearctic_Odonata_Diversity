library(dplyr)
library(ggplot2)

# read in rangeSize data
gb <- read.csv("data/rangeSize_byTraits.csv")

ll <- ggplot(gb) +
  geom_violin(mapping = aes(x = LL, y = rangeSize, fill = LL)) +
  geom_boxplot(mapping = aes(x = LL, y = rangeSize), fill = NA) +
  scale_fill_manual(values = c("#BE64AC", "#5AC8C8")) +
  scale_y_log10() +
  labs(y = "Range size", x = "",
       fill = "") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.position = "none")

# do this for forest now 
fnf <- ggplot(gb) +
  geom_violin(mapping = aes(x = Fst, y = rangeSize, fill = Fst)) +
  geom_boxplot(mapping = aes(x = Fst, y = rangeSize), fill = NA) +
  scale_fill_manual(values = c("#73AE80", "#6C83B5"))+
  scale_y_log10() +
  labs(y = "Range size", x = "",
       fill = "") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.position = "none")

cp <- cowplot::plot_grid(ll, fnf, labels = c("A", "B"))
ggsave("figures/nearctic/Figure5_RangeSize.png", dpi = 450, width = 8, height = 4.25)

ggplot(gb) +
  geom_violin(mapping = aes(x = Fst, y = rangeSize, fill = Fst)) +
  geom_boxplot(mapping = aes(x = Fst, y = rangeSize), fill = NA) +
  geom_violin(mapping = aes(x = LL, y = rangeSize, fill = LL)) +
  geom_boxplot(mapping = aes(x = LL, y = rangeSize), fill = NA) +
  scale_fill_manual(values = c("#73AE80", "#6C83B5",
                               "#BE64AC", "#5AC8C8"))+
  scale_y_continuous(expand = c(0,0)) +
  labs(y = "Range size", x = "",
       fill = "") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.position = "none")

ggsave("figures/nearctic/Figure5_RangeSize_V2.png", dpi = 450, width = 10, height = 4)

