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

ReadFolder  <- function(input.folder, n.traits, sel.type){
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

    aux.file = paste(input.folder, "pop.parameters.txt", sep="/")
    parameters = scan(aux.file, character())
    index = which("theta"==parameters)+2
    selection.strength = as.numeric(parameters[index])
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

CalcIsoFlex  <- function(mat.list){
    betas = rep(c(1, 0), each=dim(mat.list[[1]])[1]/2)
    betas = Normalize(betas)
    out = lapply(mat.list, function(mat) Flexibility(betas, mat))
    return(unlist(out))
}

CalcIsoEvol  <- function(mat.list){
    betas = rep(c(1, 0), each=dim(mat.list[[1]])[1]/2)
    betas = Normalize(betas)
    out = lapply(mat.list, function(mat) betas%*%(mat%*%betas)/sum(diag(mat)))
    return(unlist(out))
}

CalcIsoAuto  <- function(mat.list){
    betas = rep(c(1, 0), each=dim(mat.list[[1]])[1]/2)
    betas = Normalize(betas)
    out = lapply(mat.list, function(mat) 1/(betas%*%solve(mat, betas)))
    return(unlist(out))
}

MapCalcR2  <- function(mat.list){
    r2.list = lapply(mat.list, CalcR2)
    return(unlist(r2.list))
}

AVGRatioCalc <- function(input.file, n.traits){
    data.corr = SetDataFrame(input.file, n.traits, T)
    data.avg = data.frame()
    data.corr$module = as.character(data.corr$module)
    data.corr$module[data.corr$module != "between module"] = "within module"
    module.name = unique(data.corr$module)
    for (i in 1:length(module.name)){
        aux.data = data.corr[data.corr$module==module.name[i],]
        aux.data = as.vector(tapply(aux.data$main, aux.data$generation, mean, module = module.name[i]))
        data.avg = rbind(data.avg, data.frame(generation = seq(data.corr$generation[1], data.corr$generation[length(data.corr$generation)]), main = aux.data, module = module.name[i]))
    }
    AVGRatio <- abs(data.avg$main[data.avg$module == "within module"])/abs(data.avg$main[data.avg$module == "between module"])
    return(AVGRatio)
}

CalcCorrOmega <- function(mat.list){
    omega = as.matrix(read.table ("input/omega.csv", header=F, sep=' '))[1:n.traits, 1:n.traits]
    omega = omega[upper.tri(omega)]
    corr.omega <- mapply(mat.list, function(x) cor(x, omega))
    return(unlist(corr.omega))
}

ReadPattern <- function(pattern = "DivSel-Rep-*", n.traits = 10, sel.type = 'divergent'){
    folders  <- dir("output/", pattern)
    main.data = list()
    for (pop in 1:(length(folders))){
        main.data[[pop]] = ReadFolder(folders[pop], n.traits, sel.type)
    }
    names(main.data) = folders
    return(main.data)
}

#main.data.div.sel = ReadPattern()
#save(main.data.div.sel, file="./div.sel.Rdata")
#main.data.corridor = ReadPattern("Coridor*")
#save(main.data.corridor, file='corridor.Rdata')
#main.data.stabilizing = ReadPattern("Stabilizing*")
#save(main.data.stabilizing, file='stabilizing.Rdata')

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
        label.vector = as.numeric(pop.list[[pop]]$selection.strengh)
        data.avg[pop,1] = stat
        data.avg[pop,2] = label.vector
    }
    data.avg = data.frame(as.numeric(data.avg[,1]), as.numeric(data.avg[,2]))
    names(data.avg) = c("stat", "Selection_Strength")
    time.series  <- ggplot(data.avg, aes(Selection_Strength, stat, group = Selection_Strength)) +
                    layer(geom = "boxplot") + scale_y_continuous(y.axis) + scale_x_continuous("Selection Strength")
    return(time.series)
}

f
load("./div.sel.Rdata")
corr.omega = LastGenStatMultiPlot(main.data.div.sel, CalcCorrOmega, "Fitness Surface Correlation") + theme_bw()
#load("./div.sel.Rdata")
#r2 = LastGenStatMultiPlot(main.data.div.sel, MapCalcR2, "Mean Squared Correlations") + theme_bw()
#ggsave("~/lg.r2.tiff")
#flex = LastGenStatMultiPlot(main.data.div.sel, CalcIsoFlex, "Directional Flexibility") + theme_bw()
#ggsave("~/lg.flex.tiff")
#evol = LastGenStatMultiPlot(main.data.div.sel, CalcIsoEvol, "Directional Evolvability") + theme_bw()
#ggsave("~/lg.evol.tiff")
#auto = LastGenStatMultiPlot(main.data.div.sel, CalcIsoAuto, "Directional Autonomy") + theme_bw()
#ggsave("~/lg.auto.tiff")
#evol = StatMultiPlot(main.data.div.sel, CalcIsoEvol,"Directional Evolvability") + theme_bw()
#ggsave("~/evol.tiff")
