### Actual Figures

source('pop.functions.r')
library(scales)
n.traits <- 10

##Creating relevant plot objects
burnin.plots <- PlotPop(pop.path, n.traits)
file.name = "p.corr.dat"
DriftAVGPlot = AVGRatioAVGPlot(file.name, pattern = "Drift*", n.traits)
file.name = "p.corr.dat"
StabilizingAVGPlot = AVGRatioAVGPlot(file.name, pattern = "Stabilizing*", n.traits)
phen.multi.plot.10000 = PhenotipeMultiPlot("DivSel-Rep-", n.traits)
pop.path <- "./output/DivSel-Rep-0-0.0038"
DivSelExample <- PlotPop(pop.path, n.traits)
p.cor.w.multi.plot.10000 = CorrOmegaMultiPlot (file.name, "DivSel-Rep-", n.traits, Label=F)
pop.path <- "./output/CoridorSel-0.005/"
CorridorExample <- PlotPop(pop.path, n.traits)
save.image("~/Dropbox/labbio/articles/EvoMod\ -\ article/images/article.plots.rdata")


load("~/Dropbox/labbio/articles/EvoMod\ -\ article/images/article.plots.rdata")
PlotFormat = function(x) tiff(x)
## figure 1

pop.path <- "output/burn_in"
PlotFormat("~/Dropbox/labbio/articles/EvoMod\ -\ article/images/burnin.p.avg.corr.tiff")
burnin.plots$p.avg.corr + scale_fill_continuous(guide = guide_legend()) + theme_bw() + theme_bw() +
                          theme_bw() + theme(legend.position="bottom") +
                          guides(colour = guide_legend(nrow = 2))
dev.off(dev.cur())

## Figure 2
PlotFormat("~/Dropbox/labbio/articles/EvoMod\ -\ article/images/drift.AVG.plot.tiff")
DriftAVGPlot + theme_bw()
dev.off(dev.cur())

## Figure 3
# Stabilizing

PlotFormat("~/Dropbox/labbio/articles/EvoMod\ -\ article/images/stabilizing.AVG.plot.tiff")
StabilizingAVGPlot + theme_bw()
dev.off(dev.cur())

## Figure 4
# Stabilizing

## Figure 5

PlotFormat("~/Dropbox/labbio/articles/EvoMod\ -\ article/images/meanphenotype-div.tiff")
phen.multi.plot.10000 + scale_fill_continuous(guide = guide_legend()) +
                          theme_bw() + theme(legend.position="bottom", legend.text = element_text(angle = 90))
dev.off(dev.cur())

## Figure 6
PlotFormat("~/Dropbox/labbio/articles/EvoMod\ -\ article/images/divsel.p.avg.corr.tiff")
DivSelExample$p.avg.corr + scale_fill_continuous(guide = guide_legend()) +
                          theme_bw() + theme(legend.position="bottom") +
                          guides(colour = guide_legend(nrow = 2))
dev.off(dev.cur())

## Figure 7

PlotFormat("~/Dropbox/labbio/articles/EvoMod\ -\ article/images/divsel.p.omega.corr.tiff")
p.cor.w.multi.plot.10000 + scale_fill_continuous(guide = guide_legend()) +
                          theme_bw() + theme(legend.position="bottom", legend.text = element_text(angle = 90))
dev.off(dev.cur())

## Figure 8

PlotFormat("~/Dropbox/labbio/articles/EvoMod\ -\ article/images/corridor.p.avg.corr.tiff")
CorridorExample$p.avg.corr + scale_fill_continuous(guide = guide_legend()) +
                      theme_bw() + theme(legend.position="bottom") +
                      guides(colour = guide_legend(nrow = 2))
dev.off(dev.cur())
