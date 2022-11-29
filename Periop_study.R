# Load required packages
library(PPRaman)
library(baseline)
library(prospectr)
library(stringr)
library(dplyr)
library(purrr)
library(scales)
library(EMSC)
library(pracma)
library(MASS)

# load in periop data
fn <- list.files("Peri_op", full.names = TRUE)

# snipping ids from filename
patient_num <- substr(basename(fn), 1, 4)
library("stringr")
patient_num <- str_remove(patient_num, "_")

cycle <- substr(basename(fn), nchar(basename(fn))-6, nchar(basename(fn))-4)
cycle <- str_remove(cycle, "_")

# load emsc references
#em_ref <- list.files("/Users/Sue/OneDrive - Swansea University/Chemotherapy_Study_SC/Water Mean/", full.names = TRUE)
#water_ref <- lapply(em_ref, read.table)[[1]]
#water <- interp1(water_ref$V1, water_ref$V2, seq(609, 1720, by=1.085), 
#                "linear")
#water_smooth <- sav_gol(water,  fl=9, forder=5)

# load in propofol
propofol_ref <- list.files("/Users/suechandler/Desktop/Propofol", full.names=TRUE)
propofol<- lapply(propofol_ref, read.table)[[1]]
propofol <- lapply(data, read_in_shift_general, abscissa = seq(609, 1717, by=1.092))

png("Propofol.png", units="cm", width=16, height=10, res=450)
plot_ly(propofol, type="l", xlab = expression("Wavenumber / cm"^-1), ylab = expression("Intensity"))
title(font.main=4, main="Propofol Spectra")
abline(h=0, col = 3, lty = 2)
dev.off()

plot_ly(x=propofol[,1], y=propofol[,2],type="scatter", mode="lines")


# read all the files in and stores as a new dataframe with the patient name

data <- lapply(fn, read.table)
shifted_spec <- lapply(data, read_in_shift_general, abscissa = seq(609, 1717, by=1.092))

all_data <- do.call("rbind", shifted_spec)

# plotter func
matplot.spec <- function(x, y, col, ...){
  matplot(x, y, ylab = "Intensity / a.u.", xlab = expression("Wavenumber / cm"^-1), type = "l", col=col, ...)
}

wavenumber <- as.data.frame(seq(609, 1717, by=1.092))

#process spectra
processed_spec <-  process_spectra(all_data, rm_bl = "pol", poly_bl_order = 4, smoothing_opt = 1,
                                   poly_smooth_order = 4, filter_length = 9, norm_meth = "norm_snv")

#proc_modpoly <- baseline.modpolyfit(all_data, degree = 4)$corrected #polynomail order = 4
#proc_norm <- norm_snv(proc_modpoly)
#proc_norm <- norm_p(proc_norm)
proc_norm <- processed_spec
library(data.table)
wavenumber <- data.frame(seq(609, 1717, by=1.092)[-c(1:((9 + 1)/2), (ncol(all_data) - ((9 + 1)/2)-1):ncol(all_data))])
matplot.spec(wavenumber[,1], t(proc_norm), col = ifelse(cycle %like% "i" & patient_num == "PT1", "red", "black"))

#patient_1 <- proc_norm[patient_num == "PT1",]
#id_pt1 <- cycle[patient_num == "PT1"]
 
# computing standard deviate based on area under the spectra
# using trapezium rule from Pracma package

auc <- trapz(x=wavenumber[,1], as.numeric(proc_norm[1,]))
auc_all <- apply(proc_norm, 1, function(y){trapz(wavenumber[,1], y)})
# boxplot - see a few outliers we need to remove
boxplot(auc_all)
summary(auc_all) 
# mean/ IQR, range
#compute standard deviation of all the auc values
stand_dev <- sd(auc_all)
# upper and lower bounds - 2 SD encompasses around 95% of the data
ub <- mean(auc_all) + (2*stand_dev)
lb <- mean(auc_all) - (2*stand_dev)

# colours for labels
palette("polychrome 36")
cols <- palette.colors(n=36, palette = "polychrome 36", recycle = FALSE)

# colour based on whether above or below the upper/lower bounds
png("outliers.png", unit="cm", width=16, height=10, res=450)
matplot.spec(wavenumber[,1], t(proc_norm[,]), lty=1, lwd=0.5, col= ifelse(auc_all > ub, cols[[1]], 
                                                          ifelse(auc_all < lb, cols[[1]],
                                                               cols[[3]])))
legend("topleft", legend = c("Spectra within bounds", "Outlying spectra"), 
       col = c(cols[[3]], cols[[1]]),  lty = 1)
title(main = "All patient spectra.")
dev.off()

# Plot patient spectra prior to removing outliers
png("patient_16.png", unit="cm", width=16, height=10, res=450)
matplot.spec(wavenumber[,1], t(proc_norm[patient_num == "PT16",]), lty=1, lwd=0.5, col= ifelse(cycle[patient_num == "PT16"] %like% "i", cols[[1]], 
                                                                               ifelse(cycle[patient_num == "PT16"] %like% "pre", cols[[2]],
                                                                                      cols[[3]])))
legend("topleft", legend = c("Pre operative"), 
       col = c(cols[[2]]),  lty = 1)
title(main = "Patient 16 spectra by timepoint")
dev.off()


ID <- paste0(patient_num, "-", cycle)

# which are being removed
ID[auc_all > ub | auc_all < lb]


# remove outliers
outliers_rm <- proc_norm[!(auc_all > ub | auc_all < lb),]
patient_num <- patient_num[!(auc_all > ub | auc_all < lb)]
cycle <- cycle[!(auc_all > ub | auc_all < lb)]

ID <- paste0(patient_num, "-", cycle)

# comment this in to remove patients 1 and 9
#outliers_rm <- outliers_rm[!(patient_num == "PT1" | patient_num == "PT9"| patient_num == "PT14"| patient_num == "PT8"), ]
#patient_num <- patient_num[!(patient_num == "PT1" | patient_num == "PT9"| patient_num == "PT14"| patient_num == "PT8")]
#cycle <- cycle[!(patient_num == "PT1" | patient_num == "PT9"| patient_num == "PT14"| patient_num == "PT8")]
# new id's after removing outliers

#EMSC included in preprocessing
#emsc_data <- EMSC(outliers_rm, reference = colMeans(outliers_rm))$corrected
emsc_data <- outliers_rm #EMSC currently removed (to reaadd add line above and remove this line)
emsc_data <- data.frame(emsc_data)
emsc_data$patient <- patient_num
emsc_data$cycle <- cycle
mean_patient <- aggregate(.~patient + cycle, emsc_data, mean)
ids <- mean_patient[,c(1:2)]
mean_patient <- mean_patient[,-c(1:2)]

# creating new data frame for each patient
pts <- unique(ids$patient)
# section creates individual patient dataframes with data and labels
for (i in 1:length(pts)){
  sub2 <- mean_patient[ids$patient == pts[i],]
  assign(pts[i], sub2)
  env = .GlobalEnv
  
  print(i)
}


for (i in 1:length(pts)){
  sub2 <- ids$cycle[ids$patient == pts[i]]
  assign(paste0("label", pts[i]), sub2)
  env = .GlobalEnv
  
  print(i)
}


#matplot.spec(wavenumber[,1], t(PT1[labelPT1 == "pre",]), col="red")
tiff("patient_8.tiff", units="in", width=7, height=4, res=600)
matplot.spec(wavenumber, t(PT8[,]), col= ifelse(labelPT8 %like% "i", cols[[1]], 
                                                 ifelse(labelPT8 == "pre", cols[[2]],
                                                        cols[[3]])), lty=1)
legend("topleft", legend = c("Pre-op", "Intra-op", "Post-op"), col=c(cols[[2]], cols[[1]], cols[[3]]),
       lty = 1)
dev.off()

# Principle component analysis
pca_data <- prcomp(mean_patient)
score <-pca_data$x
loadings <- pca_data$rotation
# subset to get rid of the intra op
scores_pre_pos <- score[!(ids$cycle %like% "i"),]
scores_pro_pos_labels <- ids$cycle[!(ids$cycle %like% "i")]
scores_pat <- ids$patient[!(ids$cycle %like% "i")]
loadings_pre_pos <- loadings[,!(ids$cycle %like% "i")]
loadings_pro_pos_labels <- ids$cycle[!(ids$cycle %like% "i")]

# PC1  vs. PC2 Pre-op/Post-op
plot(scores_pre_pos[,3], scores_pre_pos[,4], col = ifelse(scores_pro_pos_labels == "pre", cols[[2]], cols[[3]]))
png("pc_pre_post.png", units="cm", width=16, height=10, res=450)
pairs(scores_pre_pos[,1:6], col = ifelse(scores_pro_pos_labels %like% "i", cols[[1]], 
                                ifelse(scores_pro_pos_labels == "pre", cols[[2]],
                                       cols[[3]])), pch = 19)
title(main = "PC scores 1-6 pre versus post operative",line = +3)
dev.off()

# pc1-pc6coloured by each stage pre/i/post
pairs(score[,1:6], col = ifelse(cycle %like% "i", cols[[1]], 
                          ifelse(cycle == "pre", cols[[2]],
                                 cols[[3]])), pch = 19)

#
mean_cycle <- aggregate(.~ID, data.frame(mean_patient, "ID" = ids$cycle), mean)
mean_cycle_id <- mean_cycle$ID
mean_cycle <- mean_cycle[,-1]]
matplot.spec(wavenumber, t(mean_cycle), col= ifelse(mean_cycle_id %like% "i", cols[[1]], 
                                                 ifelse(mean_cycle_id == "pre", cols[[2]],
                                                        cols[[3]])))
legend("topleft", legend = c("Pre-op", "Intra-op", "Post-op"), col=c(cols[[2]], cols[[1]], cols[[3]]),
       lty = 1)

# loading on PC1 and PC2, and PC3
plot(x=wavenumber[,1], y =loadings_pre_pos[,1], type = "l")
plot(x=wavenumber[,1], y =loadings_pre_pos[,2], type = "l")
plot(x=wavenumber[,1], y =loadings_pre_pos[,3], type = "l")
plot(x=wavenumber[,1], y =loadings_pre_pos[,4], type = "l")
plot(x=wavenumber[,1], y =loadings_pre_pos[,5], type = "l")
plot(x=wavenumber[,1], y =loadings_pre_pos[,6], type = "l")
plot(x=wavenumber[,1], y =loadings_pre_pos[,7], type = "l")
plot(x=wavenumber[,1], y =loadings_pre_pos[,8], type = "l")
plot(x=wavenumber[,1], y =loadings_pre_pos[,9], type = "l")
plot(x=wavenumber[,1], y =loadings_pre_pos[,9], type = "l")
plot(x=wavenumber[,1], y =loadings_pre_pos[,10], type = "l")

png("PC1loading plus mean spec.png", units="cm", width=16, height=10, res=450)
plot(wavenumber[,1],colMeans(mean_patient), col  = "red", type = "l",xlab = expression("Wavenumber / cm"^-1), ylab = expression("PC-1 Loading "), ylim=c(-0.1, 1.2))
lines(x=wavenumber[,1], y =((loadings_pre_pos[,1]-min(loadings_pre_pos[,1]))/(max(loadings_pre_pos[,1]) - min(loadings_pre_pos[,1]))), col = "blue")
dev.off()

png("PC2loading plus mean spec.png", units="cm", width=16, height=10, res=450)
plot(wavenumber[,1],colMeans(mean_patient), col  = "red", type = "l", xlab = expression("Wavenumber / cm"^-1), ylab = expression("PC-2 Loading "), ylim=c(-0.1, 1.2))
lines(x=wavenumber[,1], y =((loadings_pre_pos[,2]-min(loadings_pre_pos[,2]))/(max(loadings_pre_pos[,2]) - min(loadings_pre_pos[,2]))), col = "blue")
dev.off()

png("PC3loading plus mean spec.png", units="cm", width=16, height=10, res=450)
plot(wavenumber[,1],colMeans(mean_patient), col  = "red", type = "l", xlab = expression("Wavenumber / cm"^-1), ylab = expression("PC-3 Loading "), ylim=c(-0.1, 1.2))
lines(x=wavenumber[,1], y =((loadings_pre_pos[,3]-min(loadings_pre_pos[,3]))/(max(loadings_pre_pos[,3]) - min(loadings_pre_pos[,3]))), col = "blue")
dev.off()

# explained variance
summ_pc <- summary(pca_data)
plot(summ_pc$importance[3,])
min(which(summ_pc$importance[3,]>=0.95)) # number of PCs which contain 95% of the variance
# LDA in R - pre/post
pc <- as.matrix(scores_pre_pos[,c(1:52)])
row.names(pc) <- scores_pat
pc_top <- cbind(pc, ifelse(scores_pro_pos_labels == "pre", 1, 0))
colnames(pc_top)[53] <- "diagnosis"
pc_top <- data.frame(pc_top)
pc_snip <- pc_top[,c("PC1", "PC2", "diagnosis")]
pc_snip$diagnosis <- factor(pc_snip$diagnosis)

library(MASS)
# Fit the model
model <- lda(diagnosis~., data = pc_snip)
lda.data <- cbind(pc_snip, predict(model)$x)
library(ggplot2)
ggplot(lda.data, aes(x=PC1, y=PC2)) +
  geom_point(aes(color = diagnosis))

library(klaR)
partimat(diagnosis ~ PC1 + PC2, data = pc_snip, method = "lda")


# Make predictions
predictions <- model %>% predict(test.transformed)
# Model accuracy
mean(predictions$class==test.transformed$Species)

# t test between prepost
t_score <- list()
p_val <- list()
normality <-list()
for (i in 1:ncol(mean_patient)){
  test <- wilcox.test(mean_patient[ids$cycle == "pre",i], mean_patient[ids$cycle %in% c("p1", "p2", "p3"),i])
  normality[[i]] <- c(shapiro.test(mean_patient[ids$cycle == "pre",i])$p.val, shapiro.test(mean_patient[ids$cycle %in% c("p1", "p2", "p3"),i])$p.val)
  t_score[[i]] <- test$statistic
  p_val[[i]] <- test$p.value
  print(i)
}
#system("say Done")
normality[[929]]
output_t <- data.frame("Wavenumber"=wavenumber[,1], "T" = unlist(t_score), "P"=unlist(p_val),"Significance"= ifelse(unlist(p_val)<0.0001, "****",
                                                                                                                    ifelse(unlist(p_val)<0.001, "***",
                                                                                                                         ifelse(unlist(p_val)<0.01, "**",
                                                                                                                                  ifelse(unlist(p_val)<0.05, "*", "ns")))))



# t test between pre and intra op 1

t_score_cy <- list()
p_val_cy <- list()
normality_cy <-list()
for (i in 1:ncol(mean_patient)){
  test <- wilcox.test(mean_patient[ids$cycle == "pre",i], mean_patient[ids$cycle %in% c("i1"),i])
  normality_cy[[i]] <- c(shapiro.test(mean_patient[ids$cycle == "pre",i])$p.val, shapiro.test(mean_patient[ids$cycle %in% c("i1"),i])$p.val)
  t_score_cy[[i]] <- test$statistic
  p_val_cy[[i]] <- test$p.value
  print(i)
}
#system("say Done")
output_t_cy <- data.frame("Wavenumber"=wavenumber[,1], "T" = unlist(t_score_cy), "P"=unlist(p_val_cy),"Significance"= ifelse(unlist(p_val_cy)<0.0001, "****",
                                                                                                                    ifelse(unlist(p_val_cy)<0.001, "***",
                                                                                                                           ifelse(unlist(p_val_cy)<0.01, "**",
                                                                                                                                  ifelse(unlist(p_val_cy)<0.05, "*", "ns")))))

# So for the comparison with different intra ops we might want to remove the 
# two controls PT1 and PT9
rm_non_cancers <- mean_patient[!(ids$patient %in% c("PT1", "PT9")),]
rm_non_cancers_ids <- ids[!(ids$patient %in% c("PT1", "PT9")),]




plot(wavenumber[,1], mean_cycle[mean_cycle_id == "pre",]-mean_cycle[mean_cycle_id == "i1",], type="l")
plot(wavenumber[,1], mean_cycle[mean_cycle_id == "pre",], type="l", col = cols[[1]])
lines(wavenumber[,1], mean_cycle[mean_cycle_id == "i1",], type="l", col = cols[[2]])


#f <- paste(names(train_raw.df)[31], "~", paste(names(train_raw.df)[-31], collapse=" + "))
# k fold CV
k <- 5
n <-  floor(nrow(pc_top)/k)
#errors <- rep(NA, k)
#mean_pat$ID =1:nrow(mean_pat)
# shuffle the order randomly - usually due to reading in you end up with cancers and controls
# bunched together at the top or bottom
random_order <- sample(unique(mean_pat$ID)) # shuffl
#rand_5 <- rep(random_order, each = 5)
mean_pat$ID <- factor(mean_pat$ID, levels = random_order)
mean_pat$flag <- factor(ifelse(mean_pat$flag == 1, 0, 1))

# random subset to ensure n=equal

shuffled_train <- mean_pat[order(mean_pat$ID),]
k_fold_res <- list()
# kfold cross val
for (i in 1:k){
  # create the fold index
  s1 = ((i-1)*n+1)
  s2 = (i*n)
  subset = s1:s2
  # subset test and train set
  train_cv = shuffled_train[-subset, 2:(ncol(shuffled_train))]
  test_cv = shuffled_train[subset,2:(ncol(shuffled_train))]
  
  # fit model of choice
  fit = randomForest(flag~., train_cv)
  
  # predict left out subset
  pred = predict(fit, test_cv[,])
  # store confusion matrix for each fold
  k_fold_res[[i]] <- table("Observed" = test_cv$flag, "Predicted" = pred)
  print(paste("Computing fold", i))
}
# compute metrics for each fold to output
all_sen <- unlist(lapply(k_fold_res, sensitivity))
all_spec <- unlist(lapply(k_fold_res, specificity))
all_ppv <- unlist(lapply(k_fold_res, posPredValue))
all_npv <- unlist(lapply(k_fold_res, negPredValue))
all_acc <- unlist(lapply(k_fold_res, function(x){(x[1,1] + x[2,2])/ sum(x)}))
k_folds <- data.frame("Sensitivity" = all_sen, "Specificity" = all_spec, "PPV" = all_ppv,
                      "NPV" = all_npv, "Accuracy" = all_acc)



# checking if difference between lower grade/early state/control
# 9, 14, 8 exclude, rest are high grade, test against patient 1 (control)
p_val_pt1 <- list()
for (i in 1:ncol(mean_patient)){
  p_val_pt1[[i]] <- wilcox.test(x = mean_patient[ids$patient != c("PT1", "PT9", "PT8", "PT14") & ids$cycle %in% c("pre"),
                                                 i],mu=PT1[labelPT1 == "pre",i])$p.val
  shapiro.test(mean_patient[ids$patient != "PT1" & ids$cycle %in% c("pre"),
                            i])
}

pt1_pval <- data.frame("Wavenumber"=wavenumber[,1], "P_val"=unlist(p_val_pt1),"Significance"= ifelse(unlist(p_val_pt1)<0.0001, "****",
                                                                                                   ifelse(unlist(p_val_pt1)<0.001, "***",
                                                                                                          ifelse(unlist(p_val_pt1)<0.01, "**",
                                                                                                                 ifelse(unlist(p_val_pt1)<0.05, "*", "ns")))))
write.csv(pt1_pval,file = "pvaluespt1post.csv",fileEncoding = "macroman")                                                           
