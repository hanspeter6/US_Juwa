### Principal Components Analysis

### read in the data
j_dat <- read.table("leukemia_array.txt")

# transpose the set
j_dat_tr <- as.data.frame(t(j_dat))

#consider some samples
j_dat_tr[1,1:3]
j_dat_tr[9,1:3]

# consider some stats of the samples
sample_means <- apply(j_dat_tr, 1, mean)
plot(sample_means, col = c(rep(3,8), rep(2,8)), pch = 19)
sample_sds <- apply(j_dat_tr, 1, sd)
plot(sample_sds, col = c(rep(3,8), rep(2,8)), pch = 19)
sample_medians <- apply(j_dat_tr,1, median)
plot(sample_medians, col = c(rep(3,8), rep(2,8)), pch = 19)

hist(j_dat$p1_good)
hist(log(scale(j_dat$p1_good)))
hist(scale(j_dat$p1_good))
hist(scale(log(j_dat$p1_good)))

# consider some distributions of variables (gene expressions)
var_means <- apply(j_dat_tr, 2, mean)
plot(var_means)
var_sds <- apply(j_dat_tr,2, sd)
plot(var_sds)
var_medians <- apply(j_dat_tr,2, median)
plot(var_medians)

hist(j_dat_tr$`1405_i_at`, breaks = 10)
hist(j_dat_tr$`200019_s_at`, breaks = 10)
hist(j_dat_tr$`202148_s_at`, breaks = 10)

hist(var_means)
hist(var_sds)
hist(var_medians)

range(var_sds)
boxplot(var_sds)

# First log transform the data and then scale (think about value of log-transform...makes for more difficult interpretation, but does it add value..)
nu_j <- as.data.frame(scale(log(j_dat_tr)))

hist(nu_j$`1405_i_at`, breaks = 10)
hist(nu_j$`200019_s_at`, breaks = 10)
hist(nu_j$`202148_s_at`, breaks = 10)

#considering correlation of genes
#j_cor <- cor(j_dat_tr)
#library(caret)
#cor_08 <- findCorrelation(j_cor, cutoff = 0.8) # 20 000+
#cor_95 <- findCorrelation(j_cor, cutoff = 0.95) # 3802
# appears to be many strong correlations

# doing pca on the unscaled data
pca_unscaled <- prcomp(j_dat_tr, center = FALSE, scale. = FALSE, retx = TRUE)
screeplot(pca_unscaled, type = "lines")
sum_unscaled <- summary(pca_unscaled)
plot(sum_unscaled$importance[3,])

# doing pca on scaled data
pca_scaled <- prcomp(j_dat_tr, center = TRUE, scale. = TRUE, retx = TRUE)
screeplot(pca_scaled, type = "lines")
sum_scaled <- summary(pca_scaled)
plot(sum_scaled$importance[3,])

# log-transforming the original dataset
j_dat_log <- as.data.frame(log(t(j_dat)))

# doing pca on log transformed data
pca_log <- prcomp(j_dat_log, center = FALSE, scale. = FALSE, retx = TRUE)
screeplot(pca_log, type = "lines")
sum_log <- summary(pca_log)
plot(sum_log$importance[3,])

# doing same on log-transformed but standardised data
pca_log_tr <- prcomp(j_dat_log, center = TRUE, scale. = TRUE, retx = TRUE)
screeplot(pca_log_tr, type = "lines")
sum_log_tr <- summary(pca_log_tr)
plot(sum_log_tr$importance[3,])

# SECOND SECTION

# consider variance across samples for each gene expression
vars <- apply(j_dat_tr, 2, FUN = var)
top_100 <- sort(vars, decreasing = TRUE)[1:100]
ind_top_100 <- which(vars %in% top_100)
j_100 <- j_dat_tr[,ind_top_100]

# try pca on unscaled untransformed dataset
pca_100_unscaled <- prcomp(j_100, center = FALSE, scale. = FALSE, retx = TRUE)
screeplot(pca_100_unscaled, type = "lines")
sum_100_unscaled <- summary(pca_100_unscaled)
plot(sum_100_unscaled$importance[3,])

# do same for scaled and log-transformed
pca_100_logscaled <- prcomp(log(j_100), center = TRUE, scale. = TRUE, retx = TRUE)
screeplot(pca_100_logscaled, type = "lines")
sum_100_logscaled <- summary(pca_100_logscaled)
plot(sum_100_logscaled$importance[3,])

# some thoughts
# consider variance of these 100 variables as percentage of total variance:
sum(vars[ind_top_100])/sum(vars) 49.2 % of total

# if get rid of highly correlated variables...use caret and ...@0.9
# takes time...do when can...

