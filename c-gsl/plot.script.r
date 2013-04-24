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
