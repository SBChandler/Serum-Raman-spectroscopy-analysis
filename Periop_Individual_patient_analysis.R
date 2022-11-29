# INdividual patient analysis - colours are based on the labels 
png("pt_2_indiv.png", units="cm", width=32, height=20, res=450)
matplot.spec(wavenumber, t(PT2[,]), lty=1, lwd=0.7, col = ifelse(labelPT2 == "pre", cols[[6]],
                                                ifelse(labelPT2 == "i1", cols[[13]],
                                                  ifelse(labelPT2 == "i2", cols[[3]],
                                                      ifelse(labelPT2 == "i3", cols[[4]],
                                                             ifelse(labelPT2 == "p1", cols[[1]], "black"))))))
                  				      
legend("topleft", legend = c("Pre operative", "Vessel clamped", "1 hour post clamp", "2 hours post clamp", "Post operative day 1"), 
                 col=c(cols[[6]], cols[[13]], cols[[3]], cols[[4]], cols[[1]]), lty=1, lwd=1, cex=0.8)
title(font.main=4, main = "Patient 2 spectra by timepoint")
dev.off()

#plotting the difference spectra of patient pre and post
png("pt_4_diff_D1_D2.png", units="cm", width=16, height=10, res=450)
plot(wavenumber[,1], PT4[labelPT4=="p1",] -PT4[labelPT4=="p2",], lty = 1 , type ="l", 
     xlab = expression("Wavenumber / cm"^-1), ylab = "Difference Spectrum")
abline(h=0, col = "red", lty = 2)
title(font.main=4, main ="Patient 4: Post operative day 1 and post operative day 2")
dev.off()

# lines function adds to the current plot, whereas plot produces a new plot
# function to normalise the spectrum
# normal <- function(x){(x - min(x))/(max(x)-min(x))}
# lines(wavenumber[,1], normal(PT1[labelPT1=="pre",] -PT1[labelPT1=="11",]), lty = 1 , type ="l")

#pca patient by patient
pc_pat <- prcomp(PT17)


# pair plot, we're selecting the $x from pc_pat which is the scores, and 1:6 is the first 6 PCs
# might need to change 1:6, 
png("pt_17_PCA.png", units="cm", width=16, height=10, res=450)
pairs(pc_pat$x[,1:3], pch=20, cex=1, col= ifelse(labelPT17 == "pre", cols[[6]],
                                          ifelse(labelPT17 == "i1", cols[7],
                                                    ifelse(labelPT17 == "p1", cols[10],"black"))), oma=c(3,3,3,15))
                                                  
                                                        
# just for pretty plotting
par(xpd = TRUE)

legend("right", legend = c("Pre operative", "Pre vessel clamp", "Post operative day 3"), 
       col=c(cols[[6]], cols[[7]], cols[[10]]), pch = 19, cex=0.6)
title(font.sub=4, sub ="Patient 17: PCA analysis", adj=0.35,  outer=FALSE)
dev.off()

# can do the absolute value too to check the actual amount comparisons
plot(wavenumber[,1], pc_pat$rotation[,3], type = "l")
plot(wavenumber[,1], abs(pc_pat$rotation[,4]), type = "l")

# % variance explained on the different principle components
summary(pc_pat)
cumpro <- cumsum(pc_pat$sdev^2 / sum(pc_pat$sdev^2))
plot(cumpro[], xlab = "PC #", ylab = "Amount of explained variance", type = "b", main = "Cumulative variance plot")

# loading on PC 1, change the ,1 to 2, 3, 4, for PC2, 3, 4 etc
png("pt17_PC1loading.png", units="cm", width=16, height=10, res=450)
plot(wavenumber[,1], pc_pat$rotation[,1], col = "blue", type = "l", xlab = expression("Wavenumber / cm"^-1), ylab = expression(""))
abline(h=0, col = 3, lty = 2)
title(font.main=4, main="PC 1 loading: 94.35% of variance")
dev.off()

png("pt17_PC2loading.png", units="cm", width=16, height=10, res=450)
plot(wavenumber[,1], pc_pat$rotation[,2], col = "blue", type = "l", xlab = expression("Wavenumber / cm"^-1), ylab = expression(""))
abline(h=0, col = 3, lty = 2)
title(font.main=4, main="PC 2 loading: 5.65% of variance")
dev.off()

png("pt17_PC3loading.png", units="cm", width=16, height=10, res=450)
plot(wavenumber[,1], pc_pat$rotation[,3], col = "blue", type = "l", xlab = expression("Wavenumber / cm"^-1), ylab = expression(""))
abline(h=0, col = 3, lty = 2)
title(font.main=4, main="PC 3 loading: 0.71% of variance")
dev.off()

png("pt17_PC4loading.png", units="cm", width=16, height=10, res=450)
plot(wavenumber[,1], pc_pat$rotation[,4], col = "blue", type = "l", xlab = expression("Wavenumber / cm"^-1), ylab = expression(")"))
abline(h=0, col = 3, lty = 2)
title(font.main=4, main="PC 4 loading: 0.85% of variance")
dev.off()



# Interactive plots to help identify peaks
library(plotly)
#basic x/y plot for a patient with only 1 observation
plot_ly(x=wavenumber[,1], y=as.numeric(t(PT3[,])), type="scatter", mode="lines")

# Now lets try with a patient with multiple observations
plot_ly(x=wavenumber[,1], y=as.numeric(propofol), type="scatter", mode="lines")

# difference spectrum
fig <- plot_ly(x=wavenumber[,1], y=as.numeric(t(PT2[labelPT2 == "pre",]))-as.numeric(t(PT2[labelPT2 == "i1",])), type="scatter", mode="lines")%>%
layout(title = 'Spectra',  xaxis = list(title = 'Wavenumber / cm^-1'), 
       yaxis = list(title = 'Intensity/a.u.'))
#fig <- fig %>% add_trace(x=wavenumber[,1], y=as.numeric(t(PT9[labelPT9 == "p1",])), type="scatter", mode= "lines")
#fig <- fig %>% add_trace(x=wavenumber[,1], y=as.numeric(t(PT3[labelPT3 == "p2",])), type="scatter", mode= "lines")
fig

#multiple from same patient
fig <- plot_ly(x=wavenumber[,1], y=as.numeric(t(PT12[labelPT12 == "pre",])), name = "Pre operative", type="scatter", mode="lines", line=list(width=1))%>%
  layout(title = 'Patient 12 Spectra',  xaxis = list(title = 'Wavenumber / cm^-1'), 
         yaxis = list(title = 'Intensity/a.u.'))
fig <- fig %>% add_trace(x=wavenumber[,1], y=as.numeric(t(PT12[labelPT12 == "i1",])), name = "Vessel clamped", type="scatter", mode= "lines", line=list(width=0.8))
fig <- fig %>% add_trace(x=wavenumber[,1], y=as.numeric(t(PT12[labelPT12 == "i2",])), name = "1 hour post clamp", type="scatter", mode= "lines", line=list(width=0.8))
fig <- fig %>% add_trace(x=wavenumber[,1], y=as.numeric(t(PT12[labelPT12 == "i3",])), name = "2 hours post clamp", type="scatter", mode= "lines", line=list(width=0.8))
#fig <- fig %>% add_trace(x=wavenumber[,1], y=as.numeric(t(PT10[labelPT10 == "i4",])), name = "", type= "scatter", mode= "lines", line=list(width=0.8))
fig <- fig %>% add_trace(x=wavenumber[,1], y=as.numeric(t(PT12[labelPT12 == "p1",])), name = "Day 3 post operative", type="scatter", mode= "lines", line=list(width=0.8))
fig

# now how to plot a loadinginteractively
fig <- plot_ly(x=wavenumber[,1], y=pc_pat$rotation[,1], name ="PC 1 loading", type="scatter", mode="lines")
fig <- fig %>% add_trace(x=wavenumber[,1], y=pc_pat$rotation[,2], name= "PC 2 loading", type="scatter", mode="lines")
fig <- fig %>% add_trace(x=wavenumber[,1], y=pc_pat$rotation[,3], name= "PC 3 loading", type="scatter", mode="lines")
fig
