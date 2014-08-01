library(Morphometrics)
library(ggplot2)
library(reshape2)
library(EvomodR)

#main.data.div.sel = ReadPattern()
#save(main.data.div.sel, file="./rdatas/div.sel.Rdata")
#non.cor.corridor = ReadPattern("Corridor", sel.type = "corridor")
#save(non.cor.corridor, file='./rdatas/non.cor.corridor.Rdata')
non.cor.stabilizing = ReadPattern("Stabilizing", sel.type = "Stabilizing", direct.sel = F)
save(non.cor.stabilizing, file='./rdatas/non.cor.stabilizing.Rdata')
non.cor.drift = ReadPattern("Drift", sel.type = "drift", direct.sel = F)
save(non.cor.drift, file='./rdatas/non.cor.drift.Rdata')

