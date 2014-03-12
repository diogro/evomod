#load("./corridor.Rdata")

#load("./div.sel.Rdata")

#nd = LastGenStatMultiPlot(main.data.div.sel, MapEffectiveDimension, "Effective Dimensionality") + theme_bw()
#ggsave("~/lg.nd.tiff")
#r2 = LastGenStatMultiPlot(main.data.div.sel, MapCalcR2, "Mean Squared Correlations") + theme_bw()
#ggsave("~/lg.r2.tiff")
#flex = LastGenStatMultiPlotWithMean(main.data.div.sel, Flexibility, "Flexibility") + theme_bw()
#ggsave("~/lg.flex.tiff")
#evol = LastGenStatMultiPlot(main.data.div.sel, function(mat.list) CalcIsoStat(mat.list, Evolvability), "Evolvability") + theme_bw()
#ggsave("~/lg.evol.tiff")
#auto = LastGenStatMultiPlotWithMean(main.data.div.sel, Autonomy, "Autonomy") + theme_bw()
#ggsave("~/lg.auto.tiff")
#avg.ratio = AVGRatioPlot(main.data.div.sel) + theme_bw()
#ggsave("~/lg.avgratio.tiff")
#corr.omega = LastGenStatMultiPlot(main.data.div.sel, CalcCorrOmega, "Fitness Surface Correlation") + theme_bw()
#ggsave("~/lg.corr.omega.tiff")

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

#load("./stabilizing.Rdata")
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

#load("./drift.Rdata")
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
#drift.stab.avg.ratio = NoSelStatMultiPlotMultiPop(main.data.drift, main.data.stabilizing, CalcAVGRatio, "AVGRatio") + theme_bw()
#ggsave("~/tiffs/ts.drift.stab.avgratio.tiff")

#burn.in.pop = ReadFolder("burn_in", sel.type = "burn.in", direct.sel=F)
#burn.in.avg = AVGRatio(burn.in.pop$p.cor, 0, F, num.cores = 10, generations = 1:20000)
#burn.in.avg['generation'] = burn.in.avg['generation']
#burn.in.avg['Probability'] = NULL
#burn.in.avg['Selection_Strength'] = NULL
#m.avg = melt(burn.in.avg[,-c(2,5)], id.vars = c('.id', 'generation'))
#m.avg = m.avg[!((m.avg['.id'] != "Full Integration") & (m.avg['variable'] == "AVG-")),]
#m.avg = m.avg[!((m.avg['.id'] == "Full Integration") & (m.avg['variable'] == "AVG+")),]
#avg.plot = ggplot(m.avg, aes(generation,
                             #value,
                             #group=interaction(variable, generation, .id),
                             #colour=interaction(.id, variable))) +
#layer(geom="point") +
#labs(x="Generation",
     #y="Average Correlation",
     #color = "Module") +
#scale_colour_discrete(labels=c("Within Module 1",
                               #"Within Module 2",
                               #"Between Modules")) + theme_bw()
#ggsave("~/Dropbox/labbio/articles/EvoMod\ -\ article/images/burnin.p.avg.corr.tiff")
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

#burn.in.avg = AVGRatio(burn.in.pop$p.cor, 0, F, num.cores = 10, generations = 1:20000)
#burn.in.avg['generation'] = burn.in.avg['generation']
#burn.in.avg['Probability'] = NULL
#burn.in.avg['Selection_Strength'] = NULL
#m.avg = melt(burn.in.avg[,-c(2,5)], id.vars = c('.id', 'generation'))
#m.avg = m.avg[!((m.avg['.id'] != "Full Integration") & (m.avg['variable'] == "AVG-")),]
#m.avg = m.avg[!((m.avg['.id'] == "Full Integration") & (m.avg['variable'] == "AVG+")),]
#avg.plot = ggplot(m.avg, aes(generation,
                             #value,
                             #group=interaction(variable, generation, .id),
                             #colour=interaction(.id, variable))) +
#layer(geom="point") +
#labs(x="Generation",
     #y="Average Correlation",
     #color = "Module") +
#scale_colour_discrete(labels=c("Within Module 1",
                               #"Within Module 2",
                               #"Between Modules")) + theme_bw()
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

