load("./rdatas/corridor.Rdata")

devtools::install_github('Evomod-R', 'diogro')
library(EvomodR)
library(ggplot2)
library(reshape2)

modules.corridor = AVGRatioPlot(main.data.corridor, TRUE)
modules.corridor = modules.corridor + theme_bw() +
  theme(legend.position = c(0, 1),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill="transparent"))
ggsave("~/lg_avg_corridor.tiff", width= 8.7, height = 13.5, units =  "cm", dpi = 600)

modules.div = AVGRatioPlot(main.data.div.sel, TRUE)
modules.div = modules.div + theme_bw() +
  theme(legend.position = c(0, 0.25),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill="transparent"))
ggsave("~/lg_avg_div.tiff", width= 15, height = 10, units =  "cm", dpi = 600)


load("./rdatas/div.sel.Rdata")


auto = LastGenStatMultiPlotWithMean(main.data.div.sel, Autonomy, "Autonomy")
auto = auto + theme_bw() +
theme(axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = c(1, 0),
      legend.justification = c(1, 0),
      legend.background = element_rect(fill="transparent"),
      legend.title = element_text("")) +
scale_colour_discrete(name = "")
ggsave("~/lg_auto.tiff", width= 8.7, height = 10, units =  "cm", dpi = 600)

avg.ratio = AVGRatioPlot(main.data.div.sel)
avg.ratio = avg.ratio  + theme_bw() +
theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = c(1, 0),
        legend.justification = c(1, 0),
        legend.background = element_rect(colour = "black"))
ggsave("~/lg_avgratio.tiff", width= 8.7, height = 10, units =  "cm", dpi = 600)

eigen.var = LastGenMultiStatMultiPlot(main.data.div.sel, function(x) EigenVar(x, 10), "Eigenvalues (% variation)")
eigen.var= eigen.var + theme_bw() +
theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
  #theme(legend.position = c(1, 0),
        #legend.justification = c(1, 0),
        #legend.background = element_rect(colour = "black"))
ggsave("~/lg_eigen_var.tiff", width= 8.7, height = 10, units =  "cm", dpi = 600)


tiff("~/divergent_plot.tiff", height = 7, width = 17.8, units="cm", res = 600)
divergent_plot = grid.arrange(avg.ratio, eigen.var, auto, ncol = 3)
dev.off()

corr.omega = LastGenStatMultiPlot(main.data.div.sel, CalcCorrOmega, "Fitness Surface Correlation (Mantel)")
corr.omega = corr.omega + theme_bw() +
theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = c(1, 0),
        legend.justification = c(1, 0),
        legend.background = element_rect(colour = "black"))
ggsave("~/lg_corr_omega.tiff", width= 8.7, height = 10, units =  "cm", dpi = 600)

corr.omega.RS = LastGenStatMultiPlot(main.data.div.sel, CalcCorrOmegaRS, "Fitness Surface Correlation (RS)")
corr.omega.RS= corr.omega.RS + theme_bw() +
theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = c(1, 0),
        legend.justification = c(1, 0),
        legend.background = element_rect(colour = "black"))
ggsave("~/lg_corr_omega_RS.tiff", width= 8.7, height = 10, units =  "cm", dpi = 600)

corr.omega.Krz = LastGenStatMultiPlot(main.data.div.sel, CalcCorrOmegaKrz, "Krznowski Correlation for 2 PC")
corr.omega.Krz= corr.omega.Krz + theme_bw() +
theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = c(1, 0),
        legend.justification = c(1, 0),
        legend.background = element_rect(colour = "black"))
ggsave("~/lg_corr_omega_Krz.tiff", width= 8.7, height = 10, units =  "cm", dpi = 600)

corr.omega.eigenvector = LastGenMultiStatMultiPlot(main.data.div.sel, CalcCorrOmegaEigenVector, "Eigenvector Correlation")
corr.omega.eigenvector= corr.omega.eigenvector + theme_bw() +
theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1)) +
scale_colour_discrete(name = "Eigenvector") +
geom_hline(aes(yintercept = 0.7), color = "black") +
  theme(legend.position = c(1, 0),
        legend.justification = c(1, 0),
        legend.background = element_rect(colour = "transparent"))
ggsave("~/lg_corr_omega_evecs.tiff", width= 8.7, height = 10, units =  "cm", dpi = 600)

tiff("~/comparison_plot.tiff", height = 17, width = 14, units="cm", res = 600)
comparison_plot = grid.arrange(corr.omega, corr.omega.RS, corr.omega.Krz, corr.omega.eigenvector, ncol = 2)
dev.off()

tiff("~/small_comparison_plot.tiff", height = 10, width = 14, units="cm", res = 600)
comparison_plot = grid.arrange(corr.omega.Krz, corr.omega.eigenvector, ncol = 2)
dev.off()

tiff("~/subspace_plot.tiff", height = 7, width = 10, units="cm", res = 600)
subspace_plot = grid.arrange(corr.omega.Krz, corr.omega.eigenvector, ncol = 2)
dev.off()


#r2 = StatMultiPlot(main.data.div.sel, MapCalcR2, "Mean Squared Correlations") + theme_bw()
#ggsave("~/ts.r2.tiff")
#flex = StatMultiPlot(main.data.div.sel, CalcIsoFlex, "Directional Flexibility") + theme_bw()
#ggsave("~/ts.flex.tiff")
#evol = StatMultiPlot(main.data.div.sel, CalcIsoEvol, "Directional Evolvability") + theme_bw()
#ggsave("~/ts.evol.tiff")
#auto = StatMultiPlot(main.data.div.sel, CalcIsoAuto, "Directional Autonomy") + theme_bw()
#ggsave("~/ts.auto.tiff")
#avg.ratio = StatMultiPlot(main.data.div.sel, CalcAVGRatio, "AVGRatio") + theme_bw()
#ggsave("~/ts.avgratio.tiff")
#corr.omega = StatMultiPlot(main.data.div.sel, CalcCorrOmega, "Fitness Surface Correlation") + theme_bw()
#ggsave("~/ts.corr.omega.tiff")

load("./rdatas/stabilizing.Rdata")
#stab.corr.omega = NoSelStatMultiPlot(main.data.stabilizing, CalcCorrOmega, "Fitness Surface Correlation") + theme_bw()
#ggsave("~/tiffs/ts.stab.corr.omega.tiff")
#stab.r2 = NoSelStatMultiPlot(main.data.stabilizing, MapCalcR2, "Mean Squared Correlations") + theme_bw()
#ggsave("~/tiffs/ts.stab.r2.tiff")
#stab.flex = NoSelStatMultiPlot(main.data.stabilizing, function(mat.list) CalcIsoStat(mat.list, Flexibility), "Directional Flexibility") + theme_bw()
#ggsave("~/tiffs/ts.stab.flex.tiff")
#stab.evol = NoSelStatMultiPlot(main.data.stabilizing, function(mat.list) CalcIsoStat(mat.list, Evolvability), "Directional Evolvability") + theme_bw()
#ggsave("~/tiffs/ts.stab.evol.tiff")
#stab.auto = NoSelStatMultiPlot(main.data.stabilizing, function(mat.list) CalcIsoStat(mat.list, Autonomy), "Directional Autonomy") + theme_bw()
#ggsave("~/tiffs/ts.stab.auto.tiff")
#stab.avg.ratio = NoSelStatMultiPlot(main.data.stabilizing, CalcAVGRatio, "AVGRatio") + theme_bw()
#ggsave("~/tiffs/ts.stab.avgratio.tiff")

load("./rdatas/drift.Rdata")
#drift.corr.omega = NoSelStatMultiPlot(main.data.drift, CalcCorrOmega, "Fitness Surface Correlation") + theme_bw()
#ggsave("~/tiffs/ts.drift.corr.omega.tiff")
#drift.r2 = NoSelStatMultiPlot(main.data.drift, MapCalcR2, "Mean Squared Correlations") + theme_bw()
#ggsave("~/tiffs/ts.drift.r2.tiff")
#drift.flex = NoSelStatMultiPlot(main.data.drift, function(mat.list) CalcIsoStat(mat.list, Flexibility), "Directional Flexibility") + theme_bw()
#ggsave("~/tiffs/ts.drift.flex.tiff")
#drift.evol = NoSelStatMultiPlot(main.data.drift, function(mat.list) CalcIsoStat(mat.list, Evolvability), "Directional Evolvability") + theme_bw()
#ggsave("~/tiffs/ts.drift.evol.tiff")
#drift.auto = NoSelStatMultiPlot(main.data.drift, function(mat.list) CalcIsoStat(mat.list, Autonomy), "Directional Autonomy") + theme_bw()
#ggsave("~/tiffs/ts.drift.auto.tiff")
#drift.avg.ratio = NoSelStatMultiPlot(main.data.drift, CalcAVGRatio, "AVGRatio") + theme_bw()
#ggsave("~/tiffs/ts.drift.avgratio.tiff")

#drift.stab.corr.omega = NoSelStatMultiPlotMultiPop(main.data.drift, main.data.stabilizing , CalcCorrOmega, "Fitness Surface Correlation") + theme_bw()
#ggsave("~/tiffs/ts.drift.stab.corr.omega.tiff")
#drift.stab.r2 = NoSelStatMultiPlotMultiPop(main.data.drift, main.data.stabilizing, MapCalcR2, "Mean Squared Correlations") + theme_bw()
#ggsave("~/tiffs/ts.drift.stab.r2.tiff")
#drift.stab.flex = NoSelStatMultiPlotMultiPop(main.data.drift, main.data.stabilizing, function(mat.list) CalcIsoStat(mat.list, Flexibility), "Directional Flexibility") + theme_bw()
#ggsave("~/tiffs/ts.drift.stab.flex.tiff")
#drift.stab.evol = NoSelStatMultiPlotMultiPop(main.data.drift, main.data.stabilizing, function(mat.list) CalcIsoStat(mat.list, Evolvability), "Directional Evolvability") + theme_bw()
#ggsave("~/tiffs/ts.drift.stab.evol.tiff")
#drift.stab.auto = NoSelStatMultiPlotMultiPop(main.data.drift, main.data.stabilizing, function(mat.list) CalcIsoStat(mat.list, Autonomy), "Directional Autonomy") + theme_bw()
#ggsave("~/tiffs/ts.drift.stab.auto.tiff")
drift.stab.avg.ratio = NoSelStatMultiPlotMultiPop(main.data.drift, main.data.stabilizing, AVGRatioSimple, "AVGRatio")
save(drift.stab.avg.ratio, file= "./rdatas/drift.stab.avg.ratio.Rdata")
drift.stab.avg.ratio = drift.stab.avg.ratio + theme_bw() +
  theme(legend.position = c(1, 1),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill="transparent"))
ggsave("~/ts_drift_stab_avgratio.tiff", width= 12, height = 10, units =  "cm", dpi = 600)

burn.in.pop = ReadFolder("burn_in", sel.type = "burn.in", direct.sel=F)
burn.in.avg = PlotCorrs(burn.in.pop$p.cor)
burn.in.avg = burn.in.avg + theme_bw() +
  theme(legend.position = c(0, 1),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill="transparent"))
ggsave("~/burnin.p.avg.corr.tiff")

library(mvtnorm)
library(cpcbp)
omega = as.matrix(read.table("~/projects/evomod/c-gsl/input/omega.csv"))
P = main.data.div.sel[[198]]$p.cov[[10000]]
write.csv2(rbind(omega, P),"~/Desktop/cpc-mats.csv"  )
Ppop = rmvnorm(20, sigma = P)
Omegapop = rmvnorm(20, sigma = omega)
pop = rbind(Ppop, Omegapop)
f = rep(c("p", "o"), e = 20)
out = phillips.cpc(pop, f)

if(!require(gridExtra)){install.packages('gridExtra'); library(gridExtra)}
if(!require(ellipse)){install.packages('ellipse'); library(ellipse)}

library(xtable)
xtable(cbind(eigen(omega)$vectors[,1:2], eigen(P)$vectors[,1:2]))

omega[omega == 0] <- 0.000001
eigen_omega = eigen(omega)

cov_mat_1 = main.data.div.sel[[198]]$p.cov[[1]]
pop_mat_omega_1 = t(eigen_omega$vectors) %*% cov_mat_1 %*% eigen_omega$vectors

cov_mat_2 = main.data.div.sel[[198]]$p.cov[[10000]]
pop_mat_omega_2= t(eigen_omega$vectors) %*% cov_mat_2 %*% eigen_omega$vectors

plot(ellipse((pop_mat_omega_2[1:2, 1:2])), col = 'red')
points(ellipse((pop_mat_omega_1[1:2, 1:2])))
