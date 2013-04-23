source('pop.functions.r')
n.traits <- 10
file.name = "p.corr.dat"
pattern = "DivSel-Short-20"
pop.path <- "output/burn_in"
burnin.plots <- PlotPop(pop.path, n.traits)

drift.pop.number  <- 0:299
drift.folders = paste("Drift-", drift.pop.number, sep='')
drift.plots = vector('list', length(drift.pop.number))
for (i in 1:length(drift.pop.number)){
    pop.folder = paste("output/", drift.folders[i], sep = '')
    print(pop.folder)
    drift.plots[[i]]  <- PlotPop(pop.folder, n.traits)
}
names(drift.plots) = drift.folders

file.name = "p.corr.dat"
p.cor.w.multi.plot.10000 = CorrOmegaMultiPlot (file.name, "DivSel-Rep-", n.traits, Label=F)
p.cor.w.multi.plot.100 = CorrOmegaMultiPlot (file.name, "DivSel-Short-100-", n.traits, Label=F)
p.cor.w.multi.plot.1000 = CorrOmegaMultiPlot (file.name, "DivSel-Short-1000-", n.traits, Label=F)
p.cor.w.multi.plot.20 = CorrOmegaMultiPlot (file.name, "DivSel-Short-20", n.traits, Label=F)
p.within.multi.plot.10000 = WithInMultiPlot (file.name, "DivSel-Rep", n.traits)
p.within.multi.plot.100 = WithInMultiPlot (file.name, "DivSel-Short-100-", n.traits)
p.within.multi.plot.1000 = WithInMultiPlot (file.name, "DivSel-Short-1000-", n.traits)
p.within.multi.plot.20 = WithInMultiPlot (file.name, "DivSel-Short-20", n.traits)
file.name = "g.corr.dat"
g.cor.w.multi.plot.10000 = CorrOmegaMultiPlot (file.name, "DivSel-Rep", n.traits, Label=F)
g.cor.w.multi.plot.100 = CorrOmegaMultiPlot (file.name, "DivSel-Short-100-", n.traits, Label=F)
g.cor.w.multi.plot.1000 = CorrOmegaMultiPlot (file.name, "DivSel-Short-1000-", n.traits, Label=F)
g.cor.w.multi.plot.20 = CorrOmegaMultiPlot (file.name, "DivSel-Short-20", n.traits, Label=F)
g.within.multi.plot.10000 = WithInMultiPlot (file.name, "DivSel-Rep", n.traits)
g.within.multi.plot.100 = WithInMultiPlot (file.name, "DivSel-Short-100-", n.traits)
g.within.multi.plot.1000 = WithInMultiPlot (file.name, "DivSel-Short-1000-", n.traits)
g.within.multi.plot.20 = WithInMultiPlot (file.name, "DivSel-Short-20-", n.traits)
phen.multi.plot.10000 = PhenotipeMultiPlot("DivSel-Rep-", n.traits)
phen.multi.plot.100 = PhenotipeMultiPlot("DivSel-Short-100-", n.traits)
phen.multi.plot.1000 = PhenotipeMultiPlot("DivSel-Short-1000-", n.traits)
phen.multi.plot.20 = PhenotipeMultiPlot("DivSel-Short-20-", n.traits)
save.image("multiplots.rdata")

### Actual Figures
source('pop.functions.r')
n.traits <- 10
PlotFormat = function(x) tiff(x)

## figure 1

pop.path <- "output/burn_in"
burnin.plots <- PlotPop(pop.path, n.traits)
PlotFormat("~/Dropbox/labbio/articles/EvoMod\ -\ article/images/burnin.p.avg.corr.tiff")
burnin.plots$p.avg.corr
dev.off(dev.cur())

## Figure 2
file.name = "p.corr.dat"
DriftAVGPlot = AVGRatioAVGPlot(file.name, pattern = "Drift*", n.traits)
PlotFormat("~/Dropbox/labbio/articles/EvoMod\ -\ article/images/drift.AVG.plot.tiff")
DriftAVGPlot
dev.off(dev.cur())

## Figure 3
# Stabilizing

## Figure 4
# Stabilizing

## Figure 5

phen.multi.plot.10000 = PhenotipeMultiPlot("DivSel-Rep-", n.traits)
PlotFormat("~/Dropbox/labbio/articles/EvoMod\ -\ article/images/meanphenotype-div.tiff")
phen.multi.plot.10000
dev.off(dev.cur())

## Figure 6
pop.path <- "./output/DivSel-Rep-0-0.0038"
DivSelExample <- PlotPop(pop.path, n.traits)
PlotFormat("~/Dropbox/labbio/articles/EvoMod\ -\ article/images/divsel.p.avg.corr.tiff")
DivSelExample$p.avg.corr
dev.off(dev.cur())

## Figure 7

p.cor.w.multi.plot.10000 = CorrOmegaMultiPlot (file.name, "DivSel-Rep-", n.traits, Label=F)
PlotFormat("~/Dropbox/labbio/articles/EvoMod\ -\ article/images/divsel.p.omega.corr.tiff")
p.cor.w.multi.plot.10000
dev.off(dev.cur())

## Figure 8

pop.path <- "./output/CoridorSel-0.005/"
CorridorExample <- PlotPop(pop.path, n.traits)
PlotFormat("~/Dropbox/labbio/articles/EvoMod\ -\ article/images/corridor.p.avg.corr.tiff")
Corridor$p.avg.corr
dev.off(dev.cur())
