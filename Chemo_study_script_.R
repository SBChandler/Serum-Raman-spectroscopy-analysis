# Load required packages
library(PPRaman)
library(baseline)
library(prospectr)
library(stringr)
library(dplyr)
library(purrr)
library(plotly)
#remotes::install_local("PPRaman_2.0.0.tar.gz")
#install.packages("PPRaman_2.0.0.tar.gz", repos=NULL, type="source")

# load in chemo data
fn <- list.files("All Patients/", full.names = TRUE)

# need to apply RC to some of the data - read in ehre
rc_dat <- as.matrix(read.table(file = "RC_mean_trimmed_oct2021-2.txt")[,1])

wavenumber<-data.frame(seq(609, 1713, by=1.089))
#wavenumber <- data.frame(seq(609, 1717, by=1.092))

# read all the files in and stores as a new dataframe with the patient name
for (i in 1:length(fn)){
  fni <- list.files(paste0(fn[i], "/785 txt/"), full.names = TRUE)
  all_spec <- lapply(fni, read.table)
  spec <- lapply(all_spec, read_in_shift_general, abscissa = seq(609, 1713, by=1.089))
  mat_spec <- do.call("rbind", spec)
  if(any(mat_spec > 1e09)){
     mat_spec_2 <- mat_spec
  } else{
     mat_spec_2 <- sweep(mat_spec, MARGIN=2, rc_dat, "*")
  }
  assign(str_replace_all(basename(fn[i]), fixed(" "), ""), as.matrix(mat_spec_2))
  env = .GlobalEnv
  
  print(i)
}


# process spectra - loops through each of the patient dataframes
for (i in 1:length(fn)){
  spectra <-eval(parse(text = str_replace_all(basename(fn[i]),
                                                       fixed(" "), "")))
  spectra <- process_spectra(spectra, rm_bl = "pol", poly_bl_order = 5, smoothing_opt = 1,
                                 poly_smooth_order = 4, filter_length = 9, norm_meth = "norm_snv")
  
  assign(paste0("Processed_", str_replace_all(basename(fn[i]), fixed(" "), "")), spectra)
  env = .GlobalEnv
  print(i)
}

# plotter func
matplot.spec <- function(x, y, col, ...){
  matplot(x, y, ylab = "Intensity / a.u.", xlab = expression("Wavenumber / cm"^-1), type = "l", col=col, ...)
}

# defining some labels 
for (i in 1:length(fn)){
  fni <- list.files(paste0(fn[i], "/785 txt/"))
  assign(paste0("label_", str_replace_all(basename(fn[i]), fixed(" "), "")),
         substr(fni, start=1, stop=8))
  env = .GlobalEnv
  print(i)
}

# colours for labels
palette("polychrome 36")
cols <- palette.colors(n=36, palette = "polychrome 36", recycle = FALSE)

cols_cycle <- function(label){
  ifelse(label == "Cycle_1_", cols[[16]],
         ifelse(label == "Cycle_2_", cols[[3]],
                ifelse(label == "Cycle_3_", cols[[4]],
                       ifelse(label == "Cycle_4_", cols [[5]],
                              ifelse(label == "Cycle_5_", cols[[7]],
                                     ifelse(label =="Cycle_6_", cols[[11]],
                                            ifelse(label == "Cycle_7_", cols[[14]],
                                                   ifelse(label == "Cycle_8_", cols[[15]],
                                                          ifelse(label == "Cycle_9_", cols[[17]],
                                                                 ifelse(label == "Cycle_10", cols[[18]],
                                                                        ifelse(label == "Cycle_11", cols[[32]],
                                                                               ifelse(label == "Cycle_0_", cols[[34]], "NA"))))))))))))
}

lab_cycle <- function(label){
  ifelse(label == "Cycle_1_", "1",
         ifelse(label == "Cycle_2_", "2",
                ifelse(label == "Cycle_3_", "3",
                       ifelse(label == "Cycle_4_", "C4",
                              ifelse(label == "Cycle_5_", "C5",
                                     ifelse(label =="Cycle_6_", "C6",
                                            ifelse(label == "Cycle_7_", "C7",
                                                   ifelse(label == "Cycle_8_", "C8",
                                                          ifelse(label == "Cycle_9_", "C9",
                                                                 ifelse(label == "Cycle_10", "C10",
                                                                        ifelse(label == "Cycle_11", "C11",
                                                                               ifelse(label == "C0", "red", "NA"))))))))))))
}


num_cycle <- function(label){
  ifelse(label == "Cycle_1_", 1,
         ifelse(label == "Cycle_2_", 2,
                ifelse(label == "Cycle_3_", 3,
                       ifelse(label == "Cycle_4_", 4,
                              ifelse(label == "Cycle_5_", 5,
                                     ifelse(label =="Cycle_6_", 6,
                                            ifelse(label == "Cycle_7_", 7,
                                                   ifelse(label == "Cycle_8_", 8,
                                                          ifelse(label == "Cycle_9_", 9,
                                                                 ifelse(label == "Cycle_10", 10,
                                                                        ifelse(label == "Cycle_11", 11,
                                                                               ifelse(label == "Cycle_0_", 0,
                                                                                      ifelse(label=="Cycle_12", 12,
                                                                                             ifelse(label == "Cycle_13", 13, "NA"))))))))))))))
}



wavenumber<- seq(609, 1713, by=1.089)[-c(1:((9 + 1)/2), (ncol(Processed_Patient7) - ((9 + 3)/2)):ncol(Processed_Patient7))]
png("patient_7.png", units="in", width=7, height=4, res=450)
matplot.spec(wavenumber, t(Processed_Patient6), col = cols_cycle(label_Patient6))
dev.off()

# Group A combining
GA_reg <- rbind(Processed_Patient23)
GA_prog <- rbind(Processed_patient10, Processed_Patient13, Processed_Patient14, Processed_Patient17, Processed_Patient24, Processed_patient38)
GA_stab <- rbind(Processed_Patient18, Processed_Patient19, Processed_Patient26)

# GA labs
GA_reg_col <- cols_cycle(c(label_Patient23))
GA_prog_col <- cols_cycle(c(label_patient10, label_Patient13, label_Patient14, label_Patient17, label_Patient24, label_patient38))
GA_stab_col <- cols_cycle(c(label_Patient18, label_Patient19, label_Patient26))

# Group B combining
GB_no <- rbind(Processed_patient1, Processed_patient5, Processed_Patient7, Processed_Patient8, Processed_patient9,
            Processed_Patient15, Processed_Patient16, Processed_Patient20, Processed_Patient25, Processed_Patient28, Processed_Patient29, Processed_patient2, Processed_Patient12)
GB_rec <- rbind(Processed_Patient35, Processed_patient39)
                
GB_no_col <- cols_cycle(c(label_patient1, label_patient5, label_Patient7, label_Patient8, label_patient9,
                       label_Patient15, label_Patient16, label_Patient20, label_Patient25, label_Patient28, label_Patient29,label_patient2, label_Patient12))

GB_rec_col <- cols_cycle(c(label_Patient35, label_patient39))

# PCA
png("pt_39_PCA.png", units="cm", width=18, height=12, res=400)
pca_1 <- prcomp(Processed_patient39, center = TRUE)
pairs(pca_1$x[,1:6], col = cols_cycle(label_patient39), pch = 20, cex = 0.9, oma=c(3,3,3,15))
par(xpd = TRUE)

legend("right", legend = c("Pre chemotherapy", "Pre chemotherapy cycle 2", "Pre chemotherapy cycle 3", "Pre chemotherapy cycle 4"), 
       col= unique(cols_cycle(label_patient39)), pch = 20, cex = 0.55)
title(font.sub=4, sub ="Patient 39: PCA of spectra by chemotherapy cycle", adj=0.35,  outer=FALSE)
dev.off()

# pre vs post PCA

ga_all <-rbind(GA_reg, GA_prog, GA_stab)
ga_all_pre <- ga_all[lab_cycle(c(GA_reg_col, GA_prog_col, GA_stab_col)) == "C0",]


# All data 
all_dat <- rbind(Processed_patient1, Processed_patient2, Processed_patient3, Processed_patient4, Processed_patient5,
                 Processed_Patient6, Processed_Patient7, Processed_Patient8, Processed_patient9, Processed_patient10,
                 Processed_Patient11, Processed_Patient12, Processed_Patient13, Processed_Patient14,
                 Processed_Patient15, Processed_Patient16, Processed_Patient17, Processed_Patient18, Processed_Patient19, Processed_Patient20,
                 Processed_Patient23,Processed_Patient24, Processed_Patient25, Processed_Patient26, Processed_Patient27,
                  Processed_Patient28, Processed_Patient29,Processed_Patient30, Processed_patient31,Processed_Patient33, Processed_Patient34, Processed_Patient35,
                 Processed_patient37, Processed_patient38, Processed_patient39)
all_col <- cols_cycle(c(label_patient1, label_patient2, label_patient3, label_patient4, label_patient5,
                        label_Patient6, label_Patient7, label_Patient8, label_patient9, label_patient10,
                        label_Patient11, label_Patient12, label_Patient13, label_Patient14,
                        label_Patient15, label_Patient16, label_Patient17, label_Patient18, label_Patient19, label_Patient20,
                        label_Patient23, label_Patient24, label_Patient25, label_Patient26, label_Patient27,
                        label_Patient28, label_Patient29, label_Patient30, label_patient31, label_Patient33,label_Patient34,label_Patient35,
                        label_patient37, label_patient38, label_patient39))
all_lab <- lab_cycle(c(label_patient1, label_patient2, label_patient3, label_patient4, label_patient5,
                        label_Patient6, label_Patient7, label_Patient8, label_patient9, label_patient10,
                        label_Patient11, label_Patient12, label_Patient13, label_Patient14,
                        label_Patient15, label_Patient16, label_Patient17, label_Patient18, label_Patient19, label_Patient20,
                        label_Patient23,  label_Patient24, label_Patient25, label_Patient26,label_Patient27,
                        label_Patient28, label_Patient29, label_Patient30, label_patient31, label_Patient33,label_Patient34,label_Patient35,
                       label_patient37, label_patient38, label_patient39))

# patient_cols 

all_dat$cycle <- num_cycle(c(label_patient1, label_patient2, label_patient3, label_patient4, label_patient5,
                             label_Patient6, label_Patient7, label_Patient8, label_patient9, label_patient10,
                             label_Patient11, label_Patient12, label_Patient13, label_Patient14,
                             label_Patient15, label_Patient16, label_Patient17, label_Patient18, label_Patient19, label_Patient20,
                             label_Patient23,  label_Patient24,  label_Patient25, label_Patient26,label_Patient27,
                             label_Patient28, label_Patient29, label_Patient30,  label_patient31, label_Patient33, label_Patient34, label_Patient35,
                             label_patient37, label_patient38, label_patient39))

patient_lab <- c(rep("p1", times =nrow(patient1)), rep("p2", times =nrow(patient2)), rep("p3", times =nrow(patient3)), rep("p4", times =nrow(patient4)), rep("p5", times =nrow(patient5)),
                 rep("p6", times =nrow(Patient6)), rep("p7", times =nrow(Patient7)), rep("p8", times =nrow(Patient8)), rep("p9", times =nrow(patient9)), rep("p10", times =nrow(patient10)),
                 rep("p11", times =nrow(Patient11)), rep("p12", times =nrow(Patient12)), rep("p13", times =nrow(Patient13)), rep("p14", times =nrow(Patient14)), rep("p15", times =nrow(Patient15)),
                 rep("p16", times =nrow(Patient16)), rep("p17", times =nrow(Patient17)), rep("p18", times =nrow(Patient18)), rep("p19", times =nrow(Patient19)), rep("p20", times =nrow(Patient20)), rep("p23", times =nrow(Patient23)), rep("p24", times =nrow(Patient24)), 
                 rep("p25", times =nrow(Patient25)), rep("p26", times =nrow(Patient26)), rep("p27", times =nrow(Patient27)), rep("p28", times =nrow(Patient28)), rep("p29", times =nrow(Patient29)), rep("p30", times =nrow(Patient30)),
                 rep("p31", times =nrow(patient31)), rep("p33", times =nrow(Patient33)), rep("p34", times =nrow(Patient34)), rep("p35", times =nrow(Patient35)), 
                 rep("p37", times =nrow(patient37)), rep("p38", times =nrow(patient38)),
                 rep("p39", times =nrow(patient39)))

all_dat$plab <- patient_lab
# grab the minum and maximum for each i.e. 1st and last reading

mean_vals <- aggregate(.~plab+cycle, all_dat, mean)
max_cycle <- mean_vals %>% group_by(plab) %>% top_n(1, cycle) # looks for max cycle by patient ID using dplyr
min_cycle <- mean_vals[mean_vals$cycle == 1,] # subset for Cycle 1

# mean based on cycle only (across all patients) either this line or 214 not both
# mean_cycle <- aggregate(.~cycle, all_dat[,-1004], mean)

# example - how to select just patient 5 and cycle 1
matplot(t(mean_vals[mean_vals$plab == "p5" & mean_vals$cycle == "1", -c(1:2)]),type = "l")
# add patient 5 cycle end too
matplot(t(mean_vals[mean_vals$plab == "p5" & mean_vals$cycle == "4", -c(1:2)]),type = "l", add = TRUE, col ="red")

# lets do the same plot but with the max_cycle and min_cycle
matplot(t(min_cycle[min_cycle$plab == "p39", -c(1:2)]),type = "l")
# add patient 5 cycle end too
matplot(t(max_cycle[max_cycle$plab == "p39", -c(1:2)]),type = "l", add = TRUE, col ="red")

# difference spectra
matplot(wavenumber,t(max_cycle[max_cycle$plab == "p39", -c(1:2)] - min_cycle[min_cycle$plab == "p39", -c(1:2)]),type = "l", col ="red")
lines(x=wavenumber, y=rep(0, times=length(wavenumber)), lty =2)



# pca C0, Cmax
dat_pc <- mean_vals[mean_vals$cycle == "1",]
ga <- c("p10", "p13", "p14", "p17","p24", "p38","p18", "p19","p26", "p23")
dat_pc_ga <- dat_pc[dat_pc$plab %in% ga,]
gb <- c("p1","p2", "p3", "p4", "p5", "p6", "p7", "p8", "p9", "p11", "p12",
        "p15", "p16", "p20", "p25", "p27", "p28", "p29", "p30", "p31","p33", "p34", "p37","p35", "p39")

dat_pc_gb <- dat_pc[dat_pc$plab %in% gb,]

#gc <- c("p35", "p39")
#dat_pc_gc <- dat_pc[dat_pc$plab %in% gc,]

# to look at sep ie RR group or MD group groups only use this line
mean_vals <- mean_vals[mean_vals$plab %in% ga,]
max_cycle <- max_cycle[max_cycle$plab %in% ga,]

# example - plot by cycle (will plot whatever I have in the MEAN VALS CYCLE!!)
cycle1_only <- mean_vals[mean_vals$cycle == "1",]
cycle0_only <- mean_vals[mean_vals$cycle == "1",]
max_min_cycles <- rbind(cycle0_only, max_cycle)

#matplot(t(cycle0_only[,-c(1:2)]), type="l", col = ifelse(cycle0_only$plab %in% ga, cols[[14]], cols[[10]]), lty = 1)
matplot(t(cycle0_only[,-c(1:2)]), type="l", col = cols[1:length(gb)], lty = 1)

# when plotting  to look at groups separately group RR (with mean vals changed use these to plot i.e line 241 run) 
matplot(t(max_cycle[,-c(1:2)]), type="l", col = ifelse(max_cycle$plab %in% c("p35", "p39"), cols[[14]], cols[[10]]), lty = 1)
legend("topleft", cex= 0.8, legend= c("Recurrent Disease", "No Recurrence") , col = c(cols[[14]], cols[[10]]), lty = 1)

matplot(t(cycle1_only[,-c(1:2)]), type="l", col = ifelse(cycle1_only$plab %in% c("p35", "p39"), cols[[14]], cols[[10]]), lty = 1)
legend("topleft", cex = 0.8, legend= c("Recurrent Disease", "No Recurrence") , col = c(cols[[14]], cols[[10]]), lty = 1)

#taking the means of the RR no recurrence and RR recurrence
cycle1_only$secondary_group <- ifelse(cycle1_only$plab %in% c("p35", "p39"), "Recurrent disease", "No recurrence")
max_cycle$secondary_group <- ifelse(max_cycle$plab %in% c("p35", "p39"), "Recurrent disease", "No recurrence")
meanGA_2groups_C1 <- aggregate(.~secondary_group, cycle1_only[,-c(1:2)], mean)
meanGA_2groups_Cmax <- aggregate(.~secondary_group, max_cycle[,-c(1:2)], mean)

png("meanspec_RRgrp_recvno_C1.png", units="cm", width=16, height=10, res=450)
matplot(wavenumber,t(meanGA_2groups_C1[,-1]), type="l", col= c(cols[[14]], cols[[10]]), lty = 1, xlab= expression("Wavenumber / cm"^-1), ylab = "Intensity/a.u")
legend("topleft", cex = 0.6, legend= c("No Recurrence", "Recurrent disease") , col = c(cols[[14]], cols[[10]]), lty = 1)
title(font.main=4, cex.main = 0.8, main ="Mean spectra pre adjuvant chemotherapy between recurrence and no recurrence patients")
dev.off()

png("meanspec_RRgrp_recvno_Cmax.png", units="cm", width=16, height=10, res=450)
matplot(t(meanGA_2groups_Cmax[,-1]), type="l", col = c(cols[[14]], cols[[10]]), lty = 1, xlab= expression("Wavenumber / cm"^-1), ylab = "Intensity/a.u")
legend("topleft", cex =0.655 , legend= c("No Recurrence", "Recurrent Disease") , col = c(cols[[14]], cols[[10]]), lty = 1)
title(font.main=4, cex.main = 0.6, main ="Mean spectra before final chemotherapy cycle between recurrence and no recurrence patients")
dev.off()

#difference spectrum between recurrent and no recurrent
png("diffspec_RRgrp_recvsno_C1.png", units="cm", width=16, height=10, res=450)
matplot(wavenumber, t(meanGA_2groups_C1[1,-1] - meanGA_2groups_C1[2,-1]), type="l", xlab = expression("Wavenumber / cm"^-1), ylab = "Difference Spectrum")
abline(h=0, col = "red", lty = 2)
title(font.main=4, cex.main = 0.6, main ="Difference spectrum pre adjuvant chemotherapy between recurrence and no recurrence patients")
dev.off()

png("diffspec_RRgrp_recvsno_Cmax.png", units="cm", width=16, height=10, res=450)
matplot(wavenumber, t(meanGA_2groups_Cmax[1,-1] - meanGA_2groups_Cmax[2,-1]), type="l", xlab = expression("Wavenumber / cm"^-1), ylab = "Difference Spectrum")
abline(h=0, col = "red", lty = 2)
title(font.main=4, cex.main = 0.6, main ="Difference spectrum before final cycle of chemotherapy between recurrence and no recurrence patients")
dev.off()


fig <- plot_ly(x=wavenumber, y=as.numeric(t(meanGA_2groups_C1[1,-1]))-as.numeric(t(meanGA_2groups_C1[2,-1])), type="scatter", mode="lines")
fig
# layout(title = 'Spectra',  xaxis = list(title = 'Wavenumber / cm^-1'), 
#  yaxis = list(title = 'Intensity/a.u.'))

fig <- plot_ly(x=wavenumber, y=as.numeric(t(meanGA_2groups_Cmax[1,-1]))-as.numeric(t(meanGA_2groups_Cmax[2,-1])), type="scatter", mode="lines")
fig

# for group A (Metastatic disease group) MAKE SURE STUFF ABOVE IS RUNNING RIGHT GROUP ETC
# when plotting  to look at groups separately group MD (with mean vals changed use these to plot i.e line 241 run) 
matplot(wavenumber, t(max_cycle[,-c(1:2)]), type="l", col = ifelse(max_cycle$plab %in% c("p18", "p19", "p26"), cols[[14]], cols[[10]]), lty = 1)
legend("topleft", cex= 0.8, legend= c("Stable disease", "Progressive disease") , col = c(cols[[14]], cols[[10]]), lty = 1)

matplot(wavenumber, t(cycle1_only[,-c(1:2)]), type="l", col = ifelse(cycle1_only$plab %in% c("p18", "p19", "p26"), cols[[14]], cols[[10]]), lty = 1)
legend("topleft", cex = 0.8, legend= c("Stable disease", "Progressive disease") , col = c(cols[[14]], cols[[10]]), lty = 1)

#means of the group
cycle1_only$secondary_group <- ifelse(cycle1_only$plab %in% c("p18", "p19", "p26"), "Stable disease", "Progressive disease")
max_cycle$secondary_group <- ifelse(max_cycle$plab %in% c("p18", "p19", "p26"), "Stable disease", "Progressive disease")
meanGA_2groups_C1 <- aggregate(.~secondary_group, cycle1_only[,-c(1:2)], mean)
meanGA_2groups_Cmax <- aggregate(.~secondary_group, max_cycle[,-c(1:2)], mean)

png("meanspec_MDgrp_stabvsprog_C1.png", units="cm", width=16, height=10, res=450)
matplot(wavenumber, t(meanGA_2groups_C1[,-1]), type="l", col= c(cols[[14]], cols[[10]]), lty = 1, xlab= expression("Wavenumber / cm"^-1), ylab = "Intensity/a.u")
legend("topleft", cex = 0.6, legend= c("Stable disease", "Progressive disease") , col = c(cols[[14]], cols[[10]]), lty = 1)
title(font.main=4, cex.main = 0.8, main ="Mean spectra pre adjuvant chemotherapy for stable and progressive disease patients")
dev.off()

png("meanspec_MDgrp_stabvsprog_Cmax.png", units="cm", width=16, height=10, res=450)
matplot(wavenumber, t(meanGA_2groups_Cmax[,-1]), type="l", col = c(cols[[14]], cols[[10]]), lty = 1, xlab= expression("Wavenumber / cm"^-1), ylab = "Intensity/a.u")
legend("topleft", cex =0.6 , legend= c("Stable disease", "Progressive disease") , col = c(cols[[14]], cols[[10]]), lty = 1)
title(font.main=4, cex.main = 0.8, main ="Mean spectra at final cycle for stable and progressive disease patients")
dev.off()

#difference spectrum between recurrent and no recurrent
png("diffspec_MDgrp_stabvsprog_C1.png", units="cm", width=16, height=10, res=450)
matplot(wavenumber, t(meanGA_2groups_C1[1,-1] - meanGA_2groups_C1[2,-1]), type="l", xlab = expression("Wavenumber / cm"^-1), ylab = "Difference Spectrum")
abline(h=0, col = "red", lty = 2)
title(font.main=4, cex.main = 0.6, main ="Difference spectrum pre adjuvant chemotherapy between stable and progressive disease patients")
dev.off()

png("diffspec_MDgrp_stavvsprog_Cmax.png", units="cm", width=16, height=10, res=450)
matplot(wavenumber, t(meanGA_2groups_Cmax[1,-1] - meanGA_2groups_Cmax[2,-1]), type="l", xlab = expression("Wavenumber / cm"^-1), ylab = "Difference Spectrum")
abline(h=0, col = "red", lty = 2)
title(font.main=4, cex.main = 0.6, main ="Difference spectrum at final cycle of chemotherapy between stable disease and progressive disease patients")
dev.off()


fig <- plot_ly(x=wavenumber, y=as.numeric(t(meanGA_2groups_C1[1,-1]))-as.numeric(t(meanGA_2groups_C1[2,-1])), type="scatter", mode="lines")
fig
# layout(title = 'Spectra',  xaxis = list(title = 'Wavenumber / cm^-1'), 
#  yaxis = list(title = 'Intensity/a.u.'))

fig <- plot_ly(x=wavenumber, y=as.numeric(t(meanGA_2groups_Cmax[1,-1]))-as.numeric(t(meanGA_2groups_Cmax[2,-1])), type="scatter", mode="lines")
fig

# MAKE SURE LINE 204 is correct not by group from below here lets plot the mean of group a and b at cycle 0 (or what cycle its set to above) and the standard deviations
cycle_0_plotting <- cycle0_only
cycle_0_plotting$group <- ifelse(cycle0_only$plab %in% ga, "Group A", "Group B")
# now compute means for eacy group and standard deviation
mean_c0 <- aggregate(.~group, data =cycle_0_plotting[,-c(1,2)], mean )
sd_c0 <- aggregate(.~group, data =cycle_0_plotting[,-c(1,2)], sd )

# now lets plot
library(scales)
plot.new
png("cycle_1grpAvgrpB.png", units="in", width=7, height=4, res=450)
matplot(wavenumber,t(mean_c0[,-c(1)]), type="l", col = ifelse(mean_c0$group == "Group A", cols[[14]], cols[[10]]), lty = 1, ylab= "Intensity a.u", xlab = expression("Wavenumber / cm"^-1))
polygon(c(wavenumber, rev(wavenumber)), c(mean_c0[1,-1]+ sd_c0[1,-1], rev(mean_c0[1,-1]- sd_c0[1,-1])), border=NA, col = alpha(cols[[14]], 0.4))
polygon(c(wavenumber, rev(wavenumber)), c(mean_c0[2,-1]+ sd_c0[2,-1], rev(mean_c0[2,-1]- sd_c0[2,-1])), border=NA, col = alpha(cols[[10]], 0.4))
legend("topleft", col = c(cols[[14]], cols[[10]]), legend=c("MD group", "RR group"), lty=1)
title(font.main=4, main ="MD and RR groups Mean and SD spectra pre adjuvant chemotherapy")
dev.off()

#prog_group <- pca_minmax <- prcomp(rbind(dat_pc_ga[,-c(1:2, 34:182)], dat_pc_gb[,-c(1:2, 34:182)]), scale=FALSE) # spectra were separating on the region ~650-810cm^-1 

pca_max <- prcomp(rbind(dat_pc_ga[,-c(1:2)], dat_pc_gb[,-c(1:2)]), scale=FALSE) 
summary(pca_max)

png("RRRversusMD_C1_PCA.png", units="cm", width=16, height=10, res=450)
pairs(pca_max$x[,1:6], col = ifelse( rbind(dat_pc_ga, dat_pc_gb)$plab %in% ga, cols[[14]],
      ifelse(rbind(dat_pc_ga, dat_pc_gb)$plab %in% gb, cols[[10]], cols[[7]])), pch = 20, cex = 0.8, oma=c(3,3,3,15))
par(xpd = TRUE)
legend("right", fill = c(cols[[14]], cols[[10]]), legend = c("MD group", "RR group"), cex=0.6)
dev.off()

summary(pca_max)

png("RRRvsMDPC4loading.png", units="cm", width=16, height=10, res=450)
plot(wavenumber, pca_min$rotation[,4], col = "blue", type = "l", xlab = expression("Wavenumber / cm"^-1), ylab = expression(""))
abline(h=0, col = 3, lty = 2)
title(font.main=4, main="PC 4 loading")
dev.off()

pairs(pca_minmax$x[,1:6], col = ifelse(rbind(dat_pc_ga, dat_pc_gb, dat_pc_gc)$plab %in% ga, cols[[3]],
                                      cols[[4]]), pch = 19, oma=c(3,3,3,15))

pairs(pca_minmax$x[,1:6], col = ifelse(rbind(dat_pc_ga, dat_pc_gb)$cycle == "1", cols[[1]],
                                       cols[[2]]), pch = 19, oma=c(3,3,3,15))
par(xpd = TRUE)
legend("bottomright", fill = c(cols[[3]], cols[[7]], cols[[11]]), legend = c("Metastatic disease", "RRR no recurrence", "RRR recurrence"))


plot(wavenumber, pca_max$rotation[,1], type = "l")
matplot.spec(wavenumber, t(rbind(dat_pc_ga, dat_pc_gb)), col = c(rep(cols[[1]], times = 11), rep(cols[[2]], times =16)))

matplot.spec(wavenumber, t(rbind(colMeans((dat_pc_ga[,-c(1:2)])), colMeans((dat_pc_gb[-c(1:2)])))), col = c(rep(cols[[3]], times = 1), rep(cols[[7]], times =1)), lty = 1, lwd = 1.2)
legend("topleft", col = c(cols[[3]], cols[[7]]), legend = c("Group A", "Group B"), lty = 1)


# pca C0, Cmax
dat_pc_minmax <- rbind(max_cycle, min_cycle)
ga <- c("p23", "p10", "p13", "p14", "p17","p24", "p38","p18", "p19", "p26" )
dat_pc_ga_minmax <- dat_pc_minmax[dat_pc_minmax$plab %in% ga,]
gb <- c("p1","p2", "p3", "p4", "p5", "p6", "p7", "p8", "p9", "p11", "p12",
        "p15", "p16", "p20", "p25", "p27", "p28", "p29", "p30", "p31","p33", "p34", "p37", "p35", "p39")

dat_pc_gb_minmax <- dat_pc_minmax[dat_pc_minmax$plab %in% gb,]

#gc <- c("p35", "p39")
#dat_pc_gc <- dat_pc[dat_pc$plab %in% gc,]

#prog_group <- pca_minmax <- prcomp(rbind(dat_pc_ga[,-c(1:2, 34:182)], dat_pc_gb[,-c(1:2, 34:182)]), scale=FALSE) # spectra were separating on the region ~650-810cm^-1 

pca_max_min <- prcomp(rbind(dat_pc_ga_minmax[,-c(1:2)], dat_pc_gb_minmax[,-c(1:2)]), scale=FALSE) 
summary(pca_max_min)

png("RRRversusMD_Cmax_PCA.png", units="cm", width=16, height=10, res=450)
pairs(pca_max_min$x[,1:6], col = ifelse( rbind(dat_pc_ga_minmax, dat_pc_gb_minmax)$plab %in% ga & rbind(dat_pc_ga_minmax, dat_pc_gb_minmax)$cycle == "1" , cols[[3]],
                                        ifelse(rbind(dat_pc_ga_minmax, dat_pc_gb_minmax)$plab %in% gb & rbind(dat_pc_ga_minmax, dat_pc_gb_minmax)$cycle == "1", cols[[7]],
                                               ifelse( rbind(dat_pc_ga_minmax, dat_pc_gb_minmax)$plab %in% ga & rbind(dat_pc_ga_minmax, dat_pc_gb_minmax)$cycle != "1", cols[[4]], cols[[10]]))), pch = 19, oma=c(3,3,3,15))
par(xpd = TRUE)
legend("right", fill = c(cols[[3]], cols[[7]]), legend = c("Metastatic disease", "Risk reduction group"), cex=0.6)
dev.off()

# Plots of groups at cycle 1 +/- SD
f1_sd <- as.numeric(apply(dat_pc_ga[,-c(1:2)], 2, sd))
mean_ga <- colMeans(dat_pc_ga[,-c(1:2)])
matplot.spec(wavenumber, mean_ga, col =cols[[3]], lty = 1, lwd = 2)
library(scales)
polygon(c(wavenumber,rev(wavenumber)),c(mean_ga+f1_sd,rev(mean_ga-f1_sd)),col=alpha("black", 0.4), border = NA)
legend("topleft", col = c(cols[[3]], "grey"), legend = c("Group A - Metastatic Disease", "± Standard Deviation"), lty = 1, lwd = c(1, 5))

f2_sd <- as.numeric(apply(dat_pc_gb[,-c(1:2)], 2, sd))
mean_gb <- colMeans(dat_pc_gb[-c(1:2)])
matplot.spec(wavenumber, mean_gb, col =cols[[2]], lty = 1, lwd = 2)
polygon(c(wavenumber,rev(wavenumber)),c(mean_gb+f2_sd,rev(mean_gb-f2_sd)),col=alpha("black", 0.4), border = NA)
legend("topleft", col = c(cols[[2]], "grey"), legend = c("Group B - No Recurrence", "± Standard Deviation"), lty = 1, lwd = c(1, 5))

matplot.spec(wavenumber, mean_ga - mean_gb, col ="black", lty = 1, lwd = 1)
lines(c(600, 1750), c(0,0), lty = 2)


# pca C0, Cmax
dat_pc <- rbind(max_cycle, min_cycle)
dat_pc_ga <- dat_pc[dat_pc$plab %in% c("p10", "p13", "p14", "p17", "p38"),-c(1:2)]
dat_pc_gb <- dat_pc[dat_pc$plab %in% c("p1", "p3", "p5", "p7", "p8", "p9", "p11",
                                       "p15", "p16", "p20", "p25", "p28", "p29", "p31", "p37", "p39",
                                       "p2", "p12", "p23", "p18", "p19", "p26"),-c(1:2)]
#prog_group <- 
pca_minmax <- data.frame(prcomp(rbind(dat_pc_ga[-c(1:2)], dat_pc_gb[,-c(1:2)]), scale=FALSE)$x)


pairs(pca_minmax[,1:6], col = c(rep(cols[[1]], times = nrow(dat_pc_ga)), rep(cols[[2]], times =nrow(dat_pc_gb))), pch = 19, oma=c(3,3,3,15))
par(xpd = TRUE)
legend("bottomright", fill = c(cols[[1]], cols[[2]]), legend = c("Progressive", "Other"))
  
pca_minmax <- pca_minmax[,-ncol(pca_minmax)]
pca_minmax$flag <- factor(c(rep(1, nrow(dat_pc_ga)), rep(0, nrow(dat_pc_gb))))

f <- paste(names(pca_minmax)[27], "~", paste(names(pca_minmax)[1:2], collapse=" + "))
library(MA)
wdbc_raw.lda <- lda(as.formula(paste(f)), data = pca_minmax[,c(1:6, 27)])
predict(wdbc_raw.lda, pca_minmax)

plot(wdbc_raw.lda)

plot(wdbc_raw.lda$xlevels)


#HCA
dist_mat <- dist(pca_minmax, method = 'euclidean')
hclust_avg <- hclust(dist_mat, method = 'average')
plot(hclust_avg)
cut_avg <- cutree(hclust_avg, k = 5)

plot(hclust_avg)
rect.hclust(hclust_avg , k = 3, border = 2:6)
abline(h = 3, col = 'red')

library(dendextend)
avg_dend_obj <- as.dendrogram(hclust_avg)
avg_col_dend <- color_branches(avg_dend_obj, h = 1.4)
plot(avg_col_dend)


seeds_df_cl <- mutate(pca_minmax, cluster = cut_avg)
count(seeds_df_cl,cluster)



# mann whitney

# perform mann whitney u-test

res_tmp <- list()
t_val <- list()
p_val <- list()
for (i in 1:1002){
  res_tmp[[i]] <- wilcox.test(x=data.frame(dat_pc_ga[,-c(1:2)])[,i], y = data.frame(dat_pc_gb[,-c(1:2)])[,i])
  t_val[[i]] <- res_tmp[[i]]$statistic
  p_val[[i]] <- res_tmp[[i]]$p.value
  print(i)
}  


plot(wavenumber, unlist(p_val), type = "l")
lines(c(600, 1800), c(0.05, 0.05), col = "red", lty = 2)
wavenumber[which(unlist(p_val) < 0.05)]


png("MD_RRR_mean_ttest_C1.png", units="cm", width=16, height=10, res=450)
matplot.spec(wavenumber, t(rbind(colMeans(dat_pc_ga[,-c(1:2)]), colMeans(dat_pc_gb[,-c(1:2)]))), col = c(rep(cols[[14]], times = 1), rep(cols[[10]], times =1)), lty = 1, lwd = 1.2)
legend("topleft", col = c(cols[[14]], cols[[10]]), legend = c("MD group", "RR group "), lty = 1)
title(font.main=4, main ="MD and RR mean spectra pre adjuvant chemotherapy")

for (i in 1:length(p_val)){
  if(p_val[[i]] <0.05){
    lines(c(wavenumber[i], wavenumber[i]), c(0, 1), lty=1, col = alpha(cols[[1]], 0.2))
  }
}
dev.off()

# test again the CRC model 150 vs. 150
# read in the model
model <- readRDS("SNV_5th_ntree=1000_mtry=5.rds")
feats <- rownames(model[[1]]$importance) # extract top 100 features the model was built on

# select the columns from feature selection, and keep the IDs and diagnosis
yourdata <- mean_vals[c(feats, "cycle", "plab"),]

# test model on your data, the model is 100 RFs so we test again each individually and average
# make sure to remove the ID column before running through the model
library(randomForest)
res <- list()
for (i in 1:100){
  res[[i]] <- predict(model[[i]], yourdata[,-c(ncol(yourdata)-1, ncol(yourdata))], type = "prob")[,2]
}
# combine the model results and compute mean
all_res <- do.call("rbind", res)
mean_res <- colMeans(all_res)
rounded_res <- round(mean_res)


table(rounded_res, yourdata$cycle)
table(rounded_res, yourdata$plab)

library(reshape2)
output <- data.frame("Patient" = yourdata$plab, "Cycle"=yourdata$cycle, "Model_output"=rounded_res)

test <- dcast(output, Patient~Cycle)
View(test)

write.csv(test, "Model_Chemo_Results_by_Cycle_all.csv")
