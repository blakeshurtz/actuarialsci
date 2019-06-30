library(readxl)
library(psych)

ttest <- read_excel("data/misc/ttest.xlsx", col_names = FALSE)
colnames(ttest) <- c("non_adpt", "apt")

###test
t.test(ttest$non_adpt, ttest$apt, paired = TRUE, alternative = "two.sided")

###wilcoxon
wilcox.test(ttest$non_adpt, ttest$apt, exact = F, paired=TRUE, conf.int = T, correct = T, conf.level = .95)
