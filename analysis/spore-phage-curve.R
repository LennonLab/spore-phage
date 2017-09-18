library(tidyverse)
library(cowplot)
spores <- read.csv("data/2017-07-27_sporulation-curve-phage.csv")

#spores[which(spores$cfus.heat > spores$cfus.control),] <- NA

spores.summ <- spores %>% 
  group_by(treatment, time) %>%
  summarize(cont.mean = mean(cfus.control), cont.sd = sd(cfus.control),
            heat.mean = mean(cfus.heat), heat.sd = sd(cfus.heat)) %>%
  group_by(treatment, time) %>%
  mutate(proportion.spores = heat.mean / cont.mean, 
         proportion.spores.sd = sum(cont.sd, heat.sd))

spores.compare <- na.omit(spores.summ[which(spores$time %in% c(6, 54)),])
spores.compare$treatment <- c("Control", "Control", "+ Inactivated Phage", "+ Inactivated Phage")
spores.compare.plot <- ggplot(spores.compare) +
  geom_point(aes(x = time, y = proportion.spores, color = treatment), size = 3) +
  scale_color_brewer(type = "qual", name = "Treatment", palette = 3, direction = -1) +
  theme_cowplot() +
  ylab("Proportion Spores") + xlab("Time (hrs)")

ggsave("figures/spore-proportion-phage.pdf")
