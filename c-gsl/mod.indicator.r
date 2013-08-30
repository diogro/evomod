source("~/projects/lem-usp-R/matrix.func.r")

ReadMatrices  <- function(input.file, n.traits){
    data.init = read.table(input.file)
    n.corrs = (n.traits*n.traits-n.traits)/2
    gen.number = data.init[seq(1,length(data.init[,1]),n.corrs+1),]
    raw.trait.means = data.init[-seq(1,length(data.init[,1]),n.corrs+1),]
    generations = length(gen.number)
    cor.matrices = list()
    for(gen in 1:generations){
        lower = 1+((gen-1)*n.corrs)
        upper = gen*n.corrs
        current.mat = matrix(1, n.traits, n.traits)
        current.mat[upper.tri(current.mat)] = raw.trait.means[lower:upper]
        current.mat[lower.tri(current.mat)] = t(current.mat)[lower.tri(current.mat)]
        cor.matrices[[gen]] = current.mat
    }
    names(cor.matrices) = gen.number
    return(cor.matrices)
}

ReadVariances  <- function(input.file, n.traits){
    data.init = read.table(input.file)
    gen.number = data.init[seq(1,length(data.init[,1]),n.traits+1),]
    raw.trait.means = data.init[-seq(1,length(data.init[,1]),n.traits+1),]
    generations = length(gen.number)
    var.vectors = list()
    for(gen in 1:generations){
        lower = 1+((gen-1)*n.traits)
        upper = gen*n.traits
         raw.trait.means[lower:upper]
        var.vectors[[gen]] = raw.trait.means[lower:upper]
    }
    names(var.vectors) = gen.number
    return(var.vectors)
}

CalcCovar  <- function(corr, vars){
    num.gens = length(vars)
    covs = list()
    for(i in 1:num.gens)
        covs[[i]] = corr[[i]]*outer(vars[[i]], vars[[i]])
    names(covs) = names(vars)
    return(covs)
}

ReadFolder  <- function(input.folder, n.traits = 10, sel.type, direct.sel = T){
    input.folder = paste("./output", input.folder, sep="/")
    input.file = paste(input.folder, "p.corr.dat", sep = '/')
    p.cor = ReadMatrices(input.file, n.traits)
    input.file = paste(input.folder, "g.corr.dat", sep = '/')
    g.cor = ReadMatrices(input.file, n.traits)
    input.file = paste(input.folder, "p.var.dat", sep = '/')
    p.var = ReadVariances(input.file, n.traits)
    input.file = paste(input.folder, "g.var.dat", sep = '/')
    g.var = ReadVariances(input.file, n.traits)
    input.file = paste(input.folder, "h.var.dat", sep = '/')
    h.var = ReadVariances(input.file, n.traits)
    input.file = paste(input.folder, "phenotype.dat", sep = '/')
    phenotype = ReadVariances(input.file, n.traits)

    p.cov = CalcCovar(p.cor, p.var)
    g.cov = CalcCovar(g.cor, g.var)

    if(direct.sel){
        aux.file = paste(input.folder, "pop.parameters.txt", sep="/")
        parameters = scan(aux.file, character())
        index = which("theta"==parameters)+2
        selection.strength = as.numeric(parameters[index])
    }
    else{
        selection.strength = 0.
    }
    out.list = list(p.cor = p.cor,
                    g.cor = g.cor,
                    p.var = p.var,
                    g.var = g.var,
                    h.var = h.var,
                    p.cov = p.cov,
                    g.cov = g.cov,
                    selection.type = sel.type,
                    selection.strength = selection.strength,
                    generation = as.numeric(names(p.var)))
    return(out.list)
}

CalcIsoStat  <- function(mat.list, Stat){
    betas = rep(c(1, 0), each=dim(mat.list[[1]])[1]/2)
    betas = Normalize(betas)
    out = lapply(mat.list, function(mat) Stat(betas, mat))
    return(unlist(out))
}

CalcMeanStat  <- function(mat.list, Stat, nsk = 1000){
    n.traits = dim(mat.list[[1]])[1]
    beta.mat <- array (rnorm (n.traits * nsk), c(n.traits, nsk))
    beta.mat <- apply (beta.mat, 2, Normalize)
    out = lapply(mat.list, function(mat) { return(mean (apply (beta.mat, 2, Stat, cov.matrix = mat)))})
    return(unlist(out))
}

Evolvability  <- function(betas, cov.matrix) betas%*%(cov.matrix%*%betas)/sum(diag(cov.matrix))
Autonomy  <- function(betas, cov.matrix) 1/(betas%*%solve(cov.matrix, betas))

MapCalcR2  <- function(mat.list){
    r2.list = lapply(mat.list, CalcR2)
    return(unlist(r2.list))
}

CalcAVGRatio <- function(mat.list){
    n.traits = dim(mat.list[[1]])[1]
    modularity.hipot = t(rbind(rep(c(1,0), each = n.traits/2), rep(c(0,1), each = n.traits/2)))
    no.hip <- dim (modularity.hipot) [2]
    m.hip.array <- array (0, c(n.traits, n.traits, no.hip + 1))
    for (N in 1:no.hip){
        for (L in 1:n.traits){
            for (M in 1:n.traits){
                m.hip.array[L,M,N] <- ifelse (modularity.hipot[L,N] & modularity.hipot[M,N], 1, 0)
            }
        }
    }
    m.hip.array[,,no.hip+1] <- as.integer (as.logical (apply (m.hip.array, c(1,2), sum)))
    m.hip.array = m.hip.array[,,3]
    mask = as.logical(m.hip.array[lower.tri(m.hip.array)])
    AVGRatio <- unlist(lapply(mat.list, function(x) mean(x[lower.tri(x)][mask])/mean(x[lower.tri(x)][!mask])))
    return(AVGRatio)
}

CalcCorrOmega <- function(mat.list){
    n.traits = dim(mat.list[[1]])[1]
    omega = as.matrix(read.table ("input/omega.csv", header=F, sep=' '))[1:n.traits, 1:n.traits]
    omega = omega[upper.tri(omega)]
    corr.omega <- lapply(mat.list, function(x) cor(x[upper.tri(x)], omega))
    return(unlist(corr.omega))
}

ReadPattern <- function(pattern = "DivSel-Rep-*",
                        n.traits = 10,
                        sel.type = 'divergent',
                        direct.sel = T){
    folders  <- dir("output/", pattern)
    main.data = list()
    for (pop in 1:(length(folders))){
        print(pop)
        main.data[[pop]] = ReadFolder(folders[pop], n.traits, sel.type, direct.sel)
    }
    names(main.data) = folders
    return(main.data)
}

StatMultiPlot <- function(pop.list, MapStatFunction, y.axis, n.traits = 10){
    require(ggplot2)
    generation.vector = pop.list[[1]]$generation
    n.gen = length(generation.vector)
    n.pop = length(pop.list)
    data.avg = array(dim=c(n.gen*n.pop, 3))
    for (pop in 1:n.pop){
        stat <- MapStatFunction((pop.list[[pop]]$p.cov))
        print(pop)
        lower = 1+((pop-1)*n.gen)
        upper = pop*n.gen
        label.vector = rep(as.numeric(pop.list[[pop]]$selection.strength), n.gen)
        data.avg[lower:upper,1] = generation.vector
        data.avg[lower:upper,2] = stat
        data.avg[lower:upper,3] = label.vector
    }
    data.avg = data.frame(as.numeric(data.avg[,1]), as.numeric(data.avg[,2]), data.avg[,3])
    names(data.avg) = c("generation", "stat", "Selection_Strength")
    time.series  <- ggplot(data.avg, aes(generation, stat, group = Selection_Strength, color=Selection_Strength)) +
                    layer(geom = "smooth") + scale_y_continuous(y.axis)
    return(time.series)
}

LastGenStatMultiPlot  <- function(pop.list, MapStatFunction, y.axis, n.traits = 10){
    require(ggplot2)
    generation.vector = pop.list[[1]]$generation
    n.gen = length(generation.vector)
    n.pop = length(pop.list)
    data.avg = array(dim=c(n.pop, 2))
    for (pop in 1:n.pop){
        stat <- MapStatFunction(list(pop.list[[pop]]$p.cov[[n.gen]]))
        print(pop)
        lower = pop
        label.vector = as.numeric(pop.list[[pop]]$selection.strength)
        data.avg[pop,1] = stat
        data.avg[pop,2] = label.vector
    }
    data.avg = data.frame(as.numeric(data.avg[,1]), as.numeric(data.avg[,2]))
    names(data.avg) = c("stat", "Selection_Strength")
    time.series  <- ggplot(data.avg, aes(Selection_Strength, stat, group = Selection_Strength)) +
                    layer(geom = "boxplot") + scale_y_continuous(y.axis) + scale_x_continuous("Selection Strength")
    return(time.series)
}

LastGenStatMultiPlotWithMean  <- function(pop.list, Stat, y.axis, n.traits = 10){
    require(ggplot2)
    require(reshape2)
    generation.vector = pop.list[[1]]$generation
    n.gen = length(generation.vector)
    n.pop = length(pop.list)
    data.avg = array(dim=c(n.pop, 3))
    for (pop in 1:n.pop){
        direct.stat <- CalcIsoStat(list(pop.list[[pop]]$p.cov[[n.gen]]), Stat)
        mean.stat <- CalcMeanStat(list(pop.list[[pop]]$p.cov[[n.gen]]), Stat)
        print(pop)
        lower = pop
        label.vector = as.numeric(pop.list[[pop]]$selection.strength)
        data.avg[pop,1] = direct.stat
        data.avg[pop,2] = mean.stat
        data.avg[pop,3] = label.vector
    }
    data.avg = data.frame(as.numeric(data.avg[,1]), as.numeric(data.avg[,2]), as.numeric(data.avg[,3]))
    names(data.avg) = c("Directional", "Mean", "Selection_Strength")
    data.avg = melt(data.avg, c("Selection_Strength"))
    time.series  <- ggplot(data.avg, aes(Selection_Strength, value, group=interaction(Selection_Strength, variable),color=variable)) +
                    layer(geom = "boxplot") + scale_y_continuous(y.axis) + scale_x_continuous("Selection Strength")
    return(time.series)
}


NoSelStatMultiPlot <- function(pop.list, Stat, y.axis, n.traits = 10){
    require(ggplot2)
    generation.vector = pop.list[[1]]$generation
    n.gen = length(generation.vector)
    n.pop = length(pop.list)
    data.avg = array(dim=c(n.gen*n.pop, 2))
    for (pop in 1:n.pop){
        stat <- CalcIsoStat((pop.list[[pop]]$p.cov), Stat)
        print(pop)
        lower = 1+((pop-1)*n.gen)
        upper = pop*n.gen
        data.avg[lower:upper,1] = generation.vector
        data.avg[lower:upper,2] = stat
    }
    data.avg = data.frame(as.numeric(data.avg[,1]), as.numeric(data.avg[,2]))
    names(data.avg) = c("generation", "stat")
    time.series  <- ggplot(data.avg, aes(generation, stat)) +
                    geom_point(alpha=1/500) +
                    geom_smooth() +
                    stat_smooth(geom="ribbon")+
                    scale_y_continuous(y.axis)  +
                    scale_x_continuous("Generation")
    return(time.series)
}

NoSelStatMultiPlotMultiPop <- function(drift.list, stab.list, StatMap, y.axis, n.traits = 10){
    require(ggplot2)
    require("plyr")
    generation.vector = drift.list[[1]]$generation
    n.gen = length(generation.vector)
    n.pop.drift = length(drift.list)
    n.pop.stab = length(stab.list)
    data.avg = array(dim=c(n.gen*n.pop.drift*n.pop.stab, 3))
    for (pop in 1:n.pop.drift){
        stat <- StatMap((drift.list[[pop]]$p.cov))
        print(pop)
        lower = 1+((pop-1)*n.gen)
        upper = pop*n.gen
        label.vector = rep("Drift", n.gen)
        data.avg[lower:upper,1] = generation.vector
        data.avg[lower:upper,2] = stat
        data.avg[lower:upper,3] = label.vector
    }
    for (pop in (n.pop.drift+1):(n.pop.drift + n.pop.stab)){
        stat <- StatMap((stab.list[[pop-n.pop.drift]]$p.cov))
        print(pop)
        lower = 1+((pop-1)*n.gen)
        upper = pop*n.gen
        label.vector = rep("Stabilizing", n.gen)
        data.avg[lower:upper,1] = generation.vector
        data.avg[lower:upper,2] = stat
        data.avg[lower:upper,3] = label.vector
    }
    data.avg = data.frame(as.numeric(data.avg[,1]),
                          as.numeric(data.avg[,2]),
                          as.character(data.avg[,3]))
    names(data.avg) = c("generation", "stat", "Selection_scheme")
    time.series  <- ggplot(data.avg, aes(generation, stat, group=Selection_scheme, color=Selection_scheme)) +
                    #geom_point(aes(colour=Selection_scheme), size=0.2, alpha=1/5000) +
                    geom_smooth() +
                    stat_smooth(span=0.99)+
                    scale_y_continuous(y.axis)  +
                    scale_x_continuous("Generation")
    return(time.series)
}



#main.data.div.sel = ReadPattern()
#save(main.data.div.sel, file="./div.sel.Rdata")
#main.data.corridor = ReadPattern("Coridor")
#save(main.data.corridor, file='corridor.Rdata')
#main.data.stabilizing = ReadPattern("Stabilizing", sel.type = "Stabilizing", direct.sel = F)
#save(main.data.stabilizing, file='stabilizing.Rdata')
#main.data.drift = ReadPattern("Drift", sel.type = "drift", direct.sel = F)
#save(main.data.drift, file='drift.Rdata')

#load("./div.sel.Rdata")

#avg.ratio = LastGenStatMultiPlot(main.data.div.sel, CalcAVGRatio, "AVGRatio") + theme_bw()
#avg.ratio = StatMultiPlot(main.data.div.sel, CalcAVGRatio, "AVGRatio") + theme_bw()

#r2 = LastGenStatMultiPlot(main.data.div.sel, MapCalcR2, "Mean Squared Correlations") + theme_bw()
#ggsave("~/lg.r2.tiff")
#flex = LastGenStatMultiPlotWithMean(main.data.div.sel, Flexibility, "Flexibility") + theme_bw()
#ggsave("~/lg.flex.tiff")
#evol = LastGenStatMultiPlot(main.data.div.sel, function(mat.list) CalcIsoStat(mat.list, Evolvability), "Evolvability") + theme_bw()
#ggsave("~/lg.evol.tiff")
#auto = LastGenStatMultiPlotWithMean(main.data.div.sel, Autonomy, "Autonomy") + theme_bw()
#ggsave("~/lg.auto.tiff")
#avg.ratio = LastGenStatMultiPlot(main.data.div.sel, CalcAVGRatio, "AVGRatio") + theme_bw()
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

load("./stabilizing.Rdata")
#stab.corr.omega = NoSelStatMultiPlot(main.data.stabilizing, CalcCorrOmega, "Fitness Surface Correlation") + theme_bw()
#ggsave("~/tiffs/ts.stab.corr.omega.tiff")
#stab.r2 = NoSelStatMultiPlot(main.data.stabilizing, MapCalcR2, "Mean Squared Correlations") + theme_bw()
#ggsave("~/tiffs/ts.stab.r2.tiff")
#stab.flex = NoSelStatMultiPlot(main.data.stabilizing, CalcIsoFlex, "Directional Flexibility") + theme_bw()
#ggsave("~/tiffs/ts.stab.flex.tiff")
#stab.evol = NoSelStatMultiPlot(main.data.stabilizing, CalcIsoEvol, "Directional Evolvability") + theme_bw()
#ggsave("~/tiffs/ts.stab.evol.tiff")
#stab.auto = NoSelStatMultiPlot(main.data.stabilizing, CalcIsoAuto, "Directional Autonomy") + theme_bw()
#ggsave("~/tiffs/ts.stab.auto.tiff")
#stab.avg.ratio = NoSelStatMultiPlot(main.data.stabilizing, CalcAVGRatio, "AVGRatio") + theme_bw()
#ggsave("~/tiffs/ts.stab.avgratio.tiff")

load("./drift.Rdata")
#drift.corr.omega = NoSelStatMultiPlot(main.data.drift, CalcCorrOmega, "Fitness Surface Correlation") + theme_bw()
#ggsave("~/tiffs/ts.drift.corr.omega.tiff")
#drift.r2 = NoSelStatMultiPlot(main.data.drift, MapCalcR2, "Mean Squared Correlations") + theme_bw()
#ggsave("~/tiffs/ts.drift.r2.tiff")
#drift.flex = NoSelStatMultiPlot(main.data.drift, CalcIsoFlex, "Directional Flexibility") + theme_bw()
#ggsave("~/tiffs/ts.drift.flex.tiff")
#drift.evol = NoSelStatMultiPlot(main.data.drift, CalcIsoEvol, "Directional Evolvability") + theme_bw()
#ggsave("~/tiffs/ts.drift.evol.tiff")
#drift.auto = NoSelStatMultiPlot(main.data.drift, CalcIsoAuto, "Directional Autonomy") + theme_bw()
#ggsave("~/tiffs/ts.drift.auto.tiff")
#drift.avg.ratio = NoSelStatMultiPlot(main.data.drift, CalcAVGRatio, "AVGRatio") + theme_bw()
#ggsave("~/tiffs/ts.drift.avgratio.tiff")

#drift.stab.corr.omega = NoSelStatMultiPlotMultiPop(main.data.drift, main.data.stabilizing , CalcCorrOmega, "Fitness Surface Correlation") + theme_bw()
drift.list = main.data.drift[1:20]
stab.list = main.data.stabilizing[1:20]
StatMap = CalcCorrOmega
y.axis = "Fitness Surface Correlation"
#ggsave("~/tiffs/ts.drift.stab.corr.omega.tiff")
ldply(
