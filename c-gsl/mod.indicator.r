library(Morphometrics)
library(ggplot2)
library(reshape2)
library(EvomodR)

non.cor.div.sel = ReadPattern(pattern = "DivSel-Rep", sel.type = "Divergent")
save(non.cor.div.sel, file="./rdatas/non.cor.div.sel.Rdata")
#main.data.div.sel = ReadPattern()
#save(main.data.div.sel, file="./rdatas/div.sel.Rdata")
#main.data.corridor = ReadPattern("Corridor", sel.type = "corridor")
#save(main.data.corridor, file='corridor.Rdata')
#main.data.stabilizing = ReadPattern("Stabilizing", sel.type = "Stabilizing", direct.sel = F)
#save(main.data.stabilizing, file='stabilizing.Rdata')
#main.data.drift = ReadPattern("Drift", sel.type = "drift", direct.sel = F)
#save(main.data.drift, file='./rdatas/drift.Rdata')

#load("./rdatas/drift.Rdata")
#load("./rdatas/corridor.Rdata")
#load("./rdatas/div.sel.Rdata")
#load("./rdatas/stabilizing.Rdata")

#time.series.drift.mantel = ldply(main.data.drift, function(x) TimeSeriesMantel(x$p.cor))
#save(time.series.drift.mantel, file = "./rdatas/ts.mantel.Rdata")
#time.series.stab.mantel = ldply(main.data.stabilizing, function(x) TimeSeriesMantel(x$p.cor))
#save(time.series.drift.mantel, time.series.stab.mantel, file = "./rdatas/ts.mantel.Rdata")
load("./rdatas/ts.mantel.Rdata")
stab = melt(time.series.stab.mantel)[,-1]
drift = melt(time.series.drift.mantel)[,-1]
stab.mantel = ddply(stab, 'variable', function(x) c(colMeans(x[2]), quantile(x[,2], 0.025), quantile(x[,2], 0.975), 'Stabilizing'))
drift.mantel = ddply(drift, 'variable', function(x) c(colMeans(x[2]), quantile(x[,2], 0.025), quantile(x[,2], 0.975), 'Drift'))
data.avg = rbind(drift.mantel, stab.mantel)
names(data.avg) = c("generation", "stat_mean", "stat_lower", "stat_upper", "Selection_scheme")
data.avg$generation = as.numeric(as.character(data.avg$generation))
data.avg[2:4] = llply(data.avg[2:4], as.numeric)
time.series  <- ggplot(data.avg, aes(generation, stat_mean, color=Selection_scheme)) +
                    geom_smooth(aes(ymin = stat_lower, ymax = stat_upper, color=Selection_scheme), data=data.avg, stat="identity") +
                    scale_y_continuous("sequential mantel")  +
                    scale_x_continuous("Generation")
