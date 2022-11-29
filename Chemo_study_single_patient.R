pca_1 <- prcomp(Processed_patient39, center = TRUE)
pairs(pca_1$x[,1:6], col = cols_cycle(label_patient39), pch = 19, oma=c(3,3,3,15))
par(xpd = TRUE)

pca_2 <- prcomp(Processed_Patient35, center = TRUE)
pairs(pca_1$x[,1:6], col = cols_cycle(label_Patient35), pch = 19, oma=c(3,3,3,15))
par(xpd = TRUE)

pc_sum <-summary(pca_1)
pc_sum2 <-summary(pca_2)
png("pca_pt39.png", units="cm", width=16, height=10, res=450)
plot(pca_1$x[,1], pca_1$x[,6], pch = 19, col = (cols_cycle(label_patient39)), xlab = paste0("PC1 - " ,round(pc_sum$importance[2,1]*100, 2),
                                                                                         "% Explained Variance"), ylab = paste0("PC6 - " , round(pc_sum$importance[2,6]*100, 2),
                                                                                                                             "% Explained Variance"), oma=c(3,3,3,15), xlim= c (-0.55,1.20), ylim=c (-0.20,0.20))

legend("topleft", cex=0.6, fill =unique(cols_cycle(label_patient39)), legend = c("Cycle 1", "Cycle 2", "Cycle 3", "Cycle 4"))

title (main= "PCA analysis by cycle of chemotherapy")
dev.off()

png("pca_pt39and39_.png", units="cm", width=16, height=10, res=450)
plot(pca_1$x[,1], pca_1$x[,6], pch = 20, col = (cols_cycle(label_patient39)), xlab = "PC1 ",  ylab = "PC6", oma=c(3,3,3,15), xlim= c (-0.55,1.85), ylim=c (-0.30,0.30))
points(pca_2$x[,1], pca_2$x[,6], pch = 20, col = (cols_cycle(label_Patient35)),xlab = "PC1 " , ylab = "PC6", oma=c(3,3,3,15), xlim=c (-0.55,1.85), ylim=c (-0.30,0.30))

legend("topleft", cex=0.5, fill =unique(cols_cycle(label_Patient35)), legend = c("Cycle 1", "Cycle 2", "Cycle 3", "Cycle 4", "Cycle 5", "Cycle 6", "Cycle 7"))

title (main= "PCA analysis by cycle of patients with recurrent disease")
dev.off()


png("patient_39.png", units="in", width=7, height=4, res=450)
matplot.spec(wavenumber, t(Processed_patient39[label_patient39 %in% c("Cycle_2_", "Cycle_3_"),]), col = c(cols[[7]], cols[[3]]))
legend("topright", c("Cycle 2", "Cycle 3"), col=c(cols[[7]], cols[[3]]))
title(main ="Patient 39 cycle 2 and cycle 3",line = +3)
dev.off()

# example - how to select just patient 5 and cycle 1
png("pt5_cyc1.png", units="in", width=7, height=4, res=450)
matplot(t(mean_vals[mean_vals$plab == "p5" & mean_vals$cycle == "1", -c(1:2)]),type = "l")
title(main ="Patient 5 cycle 1", xlab = "Wavenumber", ylab = "Intensity")
dev.off()

# Interactive plots to help identify peaks
library(plotly)

#multiple from same patient
fig <--plot_ly(x=wavenumber, y=as.numeric(t(mean_vals[mean_vals$plab == "p39" & mean_vals$cycle == "1", -c(1:2)])), type="scatter", mode="lines",line=list(width=0.8))
fig %>% add_trace(x=wavenumber, y=as.numeric(t(mean_vals[mean_vals$plab == "p39" & mean_vals$cycle == "2", -c(1:2)])), type="scatter", mode= "lines", line=list(width=0.8))
fig %>% add_trace(x=wavenumber, y=as.numeric(t(mean_vals[mean_vals$plab == "p39" & mean_vals$cycle == "3", -c(1:2)])), type="scatter", mode= "lines", line=list(width=0.8))
fig %>% add_trace(x=wavenumber, y=as.numeric(t(mean_vals[mean_vals$plab == "p39" & mean_vals$cycle == "4", -c(1:2)])), type="scatter", mode= "lines", line=list(width=0.8))
fig


# example - how to select just patient 5 and cycle 1
matplot(t(mean_vals[mean_vals$plab == "p" & mean_vals$cycle == "1", -c(1:2)]),type = "l", col= "red")
# add patient 5 cycle end too
matplot(t(mean_vals[mean_vals$plab == "p9" & mean_vals$cycle == "2", -c(1:2)]),type = "l", add = TRUE, col= "green")
matplot(t(mean_vals[mean_vals$plab == "p9" & mean_vals$cycle == "3", -c(1:2)]),type = "l", add = TRUE, col= "black")
matplot(t(mean_vals[mean_vals$plab == "p9" & mean_vals$cycle == "4", -c(1:2)]),type = "l", add = TRUE, col= "blue")
matplot(t(mean_vals[mean_vals$plab == "p8" & mean_vals$cycle == "4", -c(1:2)]),type = "l", add = TRUE, col ="yellow")
matplot(t(mean_vals[mean_vals$plab == "p7" & mean_vals$cycle == "5", -c(1:2)]),type = "l", add = TRUE, col= "orange" )
matplot(t(mean_vals[mean_vals$plab == "p3" & mean_vals$cycle == "7", -c(1:2)]),type = "l", add = TRUE, col= "cyan" )
matplot(t(mean_vals[mean_vals$plab == "p2" & mean_vals$cycle == "8", -c(1:2)]),type = "l", add = TRUE, col= "magenta")
matplot(t(mean_vals[mean_vals$plab == "p2" & mean_vals$cycle == "9", -c(1:2)]),type = "l", add = TRUE, col= "gray" )
matplot(t(mean_vals[mean_vals$plab == "p2" & mean_vals$cycle == "10", -c(1:2)]),type = "l", add = TRUE, col= "red" )
matplot(t(mean_vals[mean_vals$plab == "p2" & mean_vals$cycle == "10", -c(1:2)]),type = "l", add = TRUE, col= "pink" )

png("allmeanspec.png", units="in", width=7, height=4, res=450)
matplot(t(mean_vals[-c(1:2)]), type = "l", col= "red", ylab = "Intensity / a.u.", xlab = expression("Wavenumber / cm"^-1))
title(main ="All Patient Spectra",line = +3)
dev.off()
