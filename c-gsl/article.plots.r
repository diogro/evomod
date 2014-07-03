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
ggsave("~/lg.avg.corridor.tiff", width= 8.7, height = 10, units =  "cm", dpi = 600)


load("./rdatas/div.sel.Rdata")

#nd = LastGenStatMultiPlot(main.data.div.sel, MapEffectiveDimension, "Effective Dimensionality") + theme_bw()
#ggsave("~/lg.nd.tiff")
#r2 = LastGenStatMultiPlot(main.data.div.sel, MapCalcR2, "Mean Squared Correlations") + theme_bw()
#ggsave("~/lg.r2.tiff")
#flex = LastGenStatMultiPlotWithMean(main.data.div.sel, Flexibility, "Flexibility") + theme_bw()
#ggsave("~/lg.flex.tiff")
#evol = LastGenStatMultiPlot(main.data.div.sel, function(mat.list) CalcIsoStat(mat.list, Evolvability), "Evolvability") + theme_bw()
#ggsave("~/lg.evol.tiff")

install.packages('gridExtra')
library(gridExtra)
auto = LastGenStatMultiPlotWithMean(main.data.div.sel, Autonomy, "Autonomy")
auto = auto + theme_bw() +
theme(axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = c(1, 0),
      legend.justification = c(1, 0),
      legend.background = element_rect(fill="transparent"),
      legend.title = element_text("")) +
scale_colour_discrete(name = "")
ggsave("~/lg.auto.tiff", width= 8.7, height = 10, units =  "cm", dpi = 600)

avg.ratio = AVGRatioPlot(main.data.div.sel)
avg.ratio = avg.ratio  + theme_bw() +
theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = c(1, 0),
        legend.justification = c(1, 0),
        legend.background = element_rect(colour = "black"))
ggsave("~/lg.avgratio.tiff", width= 8.7, height = 10, units =  "cm", dpi = 600)

eigen.var = LastGenMultiStatMultiPlot(main.data.div.sel, function(x) EigenVar(x, 3), "Eigenvalues")
eigen.var= eigen.var + theme_bw() +
theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
  #theme(legend.position = c(1, 0),
        #legend.justification = c(1, 0),
        #legend.background = element_rect(colour = "black"))
ggsave("~/lg.eigen.var.tiff", width= 8.7, height = 10, units =  "cm", dpi = 600)


tiff("~/divergent_plot.tiff", height = 7, width = 17.8, units="cm", res = 600)
divergent_plot = grid.arrange(avg.ratio, eigen.var, auto, ncol = 3)
dev.off()

corr.omega = LastGenStatMultiPlot(main.data.div.sel, CalcCorrOmega, "Fitness Surface Correlation")
corr.omega = corr.omega + theme_bw() +
theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = c(1, 0),
        legend.justification = c(1, 0),
        legend.background = element_rect(colour = "black"))
ggsave("~/lg.corr.omega.tiff", width= 8.7, height = 10, units =  "cm", dpi = 600)

corr.omega.RS = LastGenStatMultiPlot(main.data.div.sel, CalcCorrOmegaRS, "Fitness Surface Correlation")
corr.omega.RS= corr.omega.RS + theme_bw() +
theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = c(1, 0),
        legend.justification = c(1, 0),
        legend.background = element_rect(colour = "black"))
ggsave("~/lg.corr.omega.RS.tiff", width= 8.7, height = 10, units =  "cm", dpi = 600)

corr.omega.Krz = LastGenStatMultiPlot(main.data.div.sel, CalcCorrOmegaKrz, "Krznowski Correlation for 2 PC")
corr.omega.Krz= corr.omega.Krz + theme_bw() +
theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = c(1, 0),
        legend.justification = c(1, 0),
        legend.background = element_rect(colour = "black"))
ggsave("~/lg.corr.omega.Krz.tiff", width= 8.7, height = 10, units =  "cm", dpi = 600)

corr.omega.eigenvector = LastGenMultiStatMultiPlot(main.data.div.sel, CalcCorrOmegaEigenVector, "Eigenvector Correlation")
corr.omega.eigenvector= corr.omega.eigenvector + theme_bw() +
theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1)) +
scale_colour_discrete(name = "Eigenvector") +
  theme(legend.position = c(1, 0),
        legend.justification = c(1, 0),
        legend.background = element_rect(colour = "black"))
ggsave("~/lg.corr.omega.evecs.tiff", width= 8.7, height = 10, units =  "cm", dpi = 600)

tiff("~/comparison_plot.tiff", height = 17, width = 14, units="cm", res = 600)
comparison_plot = grid.arrange(corr.omega, corr.omega.RS, corr.omega.Krz, corr.omega.eigenvector, ncol = 2)
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
##Creating relevant plot objects

#source('mod.indicator.r')
### Actual Figures

library(ggplot2)
library(scales)
load("~/Dropbox/labbio/articles/EvoMod\ -\ article/images/article.plots.rdata")
PlotFormat.Small = function(x) ggsave(x, current.plot, width=3.27, height=5)
PlotFormat.Med = function(x) ggsave(x, current.plot, width=4.86, height=5)
PlotFormat.Large = function(x) ggsave(x, current.plot, width=6.83, height=4)
PlotFormat.Large.Short = function(x) ggsave(x, current.plot, width=6.83, height=3)

## figure 1

burn.in.avg = PlotCorrs(burn.in.pop$p.cor)
#PlotFormat.Small("~/Dropbox/labbio/articles/EvoMod\ -\ article/images/burnin.p.avg.corr.tiff")

## Figure 2
# Stabilizing Drift
#drift.stab.avg.ratio = NoSelStatMultiPlotMultiPop(main.data.drift, main.data.stabilizing, CalcAVGRatio, "AVGRatio") + theme_bw()
PlotFormat.Large.Short("~/Dropbox/labbio/articles/EvoMod\ -\ article/images/ts.drift.stab.avgratio.tiff")

## Figure 3
# AVG directional
#avg.ratio = AVGRatioPlot(main.data.div.sel) + theme_bw()
PlotFormat.Large.Short("~/Dropbox/labbio/articles/EvoMod\ -\ article/images/lg.avgratio.tiff")

## Figure 4
#corr.omega = LastGenStatMultiPlot(main.data.div.sel, CalcCorrOmega, "Fitness Surface Correlation") + theme_bw()
#PlotFormat.Large("~/lg.corr.omega.tiff")

## Figure 5
#auto = LastGenStatMultiPlotWithMean(main.data.div.sel, Autonomy, "Autonomy") + theme_bw()
#PlotFormat.Large.Short("~/lg.auto.tiff")

## Figure 6

