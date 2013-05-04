##Creating relevant plot objects

#source('pop.functions.r')
#n.traits <- 10
#pop.path <- "output/burn_in"
#burnin.plots <- PlotPop(pop.path, n.traits)
#file.name = "p.corr.dat"
#DriftAVGPlot = AVGRatioAVGPlot(file.name, pattern = "Drift*", n.traits)
#file.name = "p.corr.dat"
#StabilizingAVGPlot = AVGRatioAVGPlot(file.name, pattern = "Stabilizing*", n.traits)
#phen.multi.plot.10000 = PhenotipeMultiPlot("DivSel-Rep-", n.traits)
#pop.path <- "./output/DivSel-Rep-0-0.0038"
#DivSelExample <- PlotPop(pop.path, n.traits)
#p.cor.w.multi.plot.10000 = CorrOmegaMultiPlot (file.name, "DivSel-Rep-", n.traits, Label=F)
#pop.path <- "./output/CoridorSel-0.005/"
#CorridorExample <- PlotPop(pop.path, n.traits)
#save.image("~/Dropbox/labbio/articles/EvoMod\ -\ article/images/article.plots.rdata")

### Actual Figures

library(ggplot2)
library(scales)
load("~/Dropbox/labbio/articles/EvoMod\ -\ article/images/article.plots.rdata")
PlotFormat.Small = function(x) ggsave(x, current.plot, width=3.27, height=5)
PlotFormat.Med = function(x) ggsave(x, current.plot, width=4.86, height=5)
PlotFormat.Large = function(x) ggsave(x, current.plot, width=6.83, height=4)

## figure 1

current.plot = burnin.plots$p.avg.corr + scale_fill_continuous(guide = guide_legend()) + theme_bw() + theme_bw() +
                          theme_bw() + theme(legend.position="bottom") +
                          guides(colour = guide_legend(nrow = 3))
PlotFormat.Small("~/Dropbox/labbio/articles/EvoMod\ -\ article/images/burnin.p.avg.corr.tiff")

## Figure 2
current.plot =DriftAVGPlot + theme_bw()
PlotFormat.Large("~/Dropbox/labbio/articles/EvoMod\ -\ article/images/drift.AVG.plot.tiff")

## Figure 3
# Stabilizing

current.plot =StabilizingAVGPlot + theme_bw()
PlotFormat.Large("~/Dropbox/labbio/articles/EvoMod\ -\ article/images/stabilizing.AVG.plot.tiff")

## Figure 4
# Stabilizing

## Figure 5

current.plot =phen.multi.plot.10000 + scale_fill_continuous(guide = guide_legend()) +
                          theme_bw() + theme(legend.position="bottom", legend.text = element_text(angle = 45))
PlotFormat.Small("~/Dropbox/labbio/articles/EvoMod\ -\ article/images/meanphenotype-div.tiff")

## Figure 6
current.plot =DivSelExample$p.avg.corr + scale_fill_continuous(guide = guide_legend()) +
                          theme_bw() + theme(legend.position="bottom") +
                          guides(colour = guide_legend(nrow = 3))
PlotFormat.Med("~/Dropbox/labbio/articles/EvoMod\ -\ article/images/divsel.p.avg.corr.tiff")

## Figure 7

current.plot =p.cor.w.multi.plot.10000 + scale_fill_continuous(guide = guide_legend()) +
                          theme_bw() + theme(legend.position="bottom", legend.text = element_text(angle = 45))
PlotFormat.Large("~/Dropbox/labbio/articles/EvoMod\ -\ article/images/divsel.p.omega.corr.tiff")

## Figure 8

current.plot =CorridorExample$p.avg.corr + scale_fill_continuous(guide = guide_legend()) +
                      theme_bw() + theme(legend.position="bottom") +
                      guides(colour = guide_legend(nrow = 3))
PlotFormat.Med("~/Dropbox/labbio/articles/EvoMod\ -\ article/images/corridor.p.avg.corr.tiff")
