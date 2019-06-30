library(readr)
library(psych)
library(tidyverse)
library(car)

c_auto_results <- read_csv("data/results/c_auto_results.csv")
c_auto_results <- as.data.frame(c_auto_results)

#remove outliers
c_auto_results <- c_auto_results %>% filter(non_adpt_cv < 5)

#distribution
TEMP <- c_auto_results$non_adpt_cv - c_auto_results$adpt_cv
TEMP
hist(TEMP, main= "Model Comparison of Variation", xlab="Differences in CV")

#t-test
t.test(c_auto_results$non_adpt_se, c_auto_results$apt_se, paired = TRUE, alternative = "two.sided")
wilcox.test(c_auto_results$non_adpt, c_auto_results$adpt, exact = F, paired=TRUE, conf.int = T, correct = T, conf.level = .95)

#boxplot
d_viz <- gather(c_auto_results, 'non_adpt_se', 'apt_se', key = "model", value = "predictions")
ggplot(d_viz, aes(x=model, y=predictions)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90))

#uniformity
hist(c_auto_results$non_adpt)
hist(c_auto_results$adpt)

##qq
qqPlot(c_auto_results$non_adpt, distribution = "unif")
qqPlot(c_auto_results$adpt, distribution = "unif")

#kstest
ks.test(c_auto_results$non_adpt, "punif")
ks.test(c_auto_results$adpt, "punif")

