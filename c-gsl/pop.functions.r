SetDataFrame  <- function(input.file, n.traits, corr.plot){
    data.init = read.table(input.file)
    if(corr.plot){
        aux.trait = rep(c(1, -2), each = n.traits/2)
        aux.trait = aux.trait%*%t(aux.trait)
        aux.trait = aux.trait[upper.tri(aux.trait)]
        aux.selection = aux.trait
        aux.trait[aux.trait==1] = "within module 1"
        aux.trait[aux.trait==4] = "within module 2"
        aux.trait[aux.trait==-2] = "between module"
        aux.selection[aux.selection==1] = "Positive"
        aux.selection[aux.selection==4] = "Negative"
        aux.selection[aux.selection==-2] = "Divergent"
        n.traits = (n.traits*n.traits-n.traits)/2
    }
    else{
        aux.trait = rep(c(1, -2), each = n.traits/2)
        aux.trait[aux.trait==1] = "module 1"
        aux.trait[aux.trait==-2] = "module 2"
    }
    gen.number = data.init[seq(1,length(data.init[,1]),n.traits+1),]
    raw.trait.means = data.init[-seq(1,length(data.init[,1]),n.traits+1),]
    generations = length(gen.number)
    gen.number = rep(gen.number, each = n.traits)
    if(corr.plot){
        sel.scheme = rep(aux.selection, generations)
    }
    else{
        sel.scheme = rep(rep(c("positive", "negative"), each=n.traits/2), generations)
    }
    module = as.character(rep(aux.trait, generations))
    trait.name = as.factor(rep(1:n.traits, generations))
    data.clean = data.frame(raw.trait.means, trait.name, gen.number, sel.scheme, module)
    names(data.clean) = c("main", "trait", "generation", "selection", "module")
    return(data.clean)
}

time.series.plot  <-  function(input.file, y.axis, n.traits, selection = T, corr.plot = F){
    require(ggplot2)
    data.clean = SetDataFrame(input.file, n.traits, corr.plot)
    if (selection)
        time.series <- ggplot(data.clean, aes(generation, main, group = trait, color = selection)) + layer (geom = "point") + scale_y_continuous(y.axis)
    else if(corr.plot)
        time.series <- ggplot(data.clean, aes(generation, main, group = trait, color = module)) + layer (geom = c("point", "line")) + scale_y_continuous(y.axis)
    else
        time.series <- ggplot(data.clean, aes(generation, main, group = trait, color = trait)) + layer (geom = "point") + scale_y_continuous(y.axis)
    return(time.series)
}

PlotPop  <- function (pop.path, n.traits){
    mean.phenotype.plot  <- time.series.plot (paste(pop.path, "phenotype.dat", sep = '/'), "mean phenotype", n.traits, T)
    g.var.plot <- time.series.plot (paste(pop.path, "g.var.dat", sep = '/'), "genetic variance", n.traits, F)
    p.var.plot <- time.series.plot (paste(pop.path, "p.var.dat", sep = '/'), "phenotypic variance", n.traits, F)
    h.var.plot <- time.series.plot (paste(pop.path, "h.var.dat", sep = '/'), "heritability", n.traits, F)
    g.corr.plot <- time.series.plot (paste(pop.path, "g.corr.dat", sep = '/'), "genetic correlations", n.traits, T, T)
    g.avg.corr.plot <- AVGCorrPlot (paste(pop.path, "g.corr.dat", sep = '/'), "mean genetic correlations", n.traits)
    p.corr.plot <- time.series.plot (paste(pop.path, "p.corr.dat", sep = '/'), "phenotypic correlations", n.traits, T, T)
    p.avg.corr.plot <- AVGCorrPlot (paste(pop.path, "p.corr.dat", sep = '/'), "mean phenotypic correlations", n.traits)
    plots = list (
          g.var = g.var.plot,
          p.var = p.var.plot,
          h.var = h.var.plot,
          mean.phenotype = mean.phenotype.plot,
          g.corr = g.corr.plot,
          p.corr = p.corr.plot,
          g.avg.corr = g.avg.corr.plot,
          p.avg.corr = p.avg.corr.plot
          )
    return (plots)
}

PlotTiffManyPops  <-  function(list.plots, file.name){
    require(ggplot2)
    require(gridExtra)
    dir.create("output/images")
    file.names = names(list.plots[[1]])
    for (i in 1:length(file.names)){
        tiff(paste("output/images/", file.name, ".", file.names[i], ".tiff", sep = ''), width = 1080, height = 1980)
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(3, 2)))
        vplayout <- function(x, y)
            viewport(layout.pos.row = x, layout.pos.col = y)
        print(list.plots[[1]][[i]], vp = vplayout(1, 1))
        print(list.plots[[2]][[i]], vp = vplayout(1, 2))
        print(list.plots[[3]][[i]], vp = vplayout(2, 1))
        print(list.plots[[4]][[i]], vp = vplayout(2, 2))
        print(list.plots[[5]][[i]], vp = vplayout(3, 1))
        print(list.plots[[6]][[i]], vp = vplayout(3, 2))
        dev.off(dev.cur())
    }
}

PlotPngSinglePop <-  function(list.plots, file.name){
    require(ggplot2)
    require(gridExtra)
    dir.create("output/images")
    file.names = names(list.plots)
    for (i in 1:length(file.names)){
        png(paste("output/images/", file.name, ".", file.names[i], ".png", sep = ''), width = 700 , height =  300)
        #postscript(paste("output/images/", file.name, ".", file.names[i], ".eps", sep = ''), width = 7 , height =  3, horizontal=F )
        print(list.plots[[i]])
        dev.off(dev.cur())
    }
}

PlotPngManyPops  <-  function(list.plots, file.name){
    require(ggplot2)
    require(gridExtra)
    dir.create("output/images")
    file.names = names(list.plots[[1]])
    for (i in 1:length(file.names)){
        png(paste("output/images/", file.name, ".", file.names[i], ".png", sep = ''), width = 1080, height = 1980)
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(3, 2)))
        vplayout <- function(x, y)
            viewport(layout.pos.row = x, layout.pos.col = y)
        print(list.plots[[1]][[i]], vp = vplayout(1, 1))
        print(list.plots[[2]][[i]], vp = vplayout(1, 2))
        print(list.plots[[3]][[i]], vp = vplayout(2, 1))
        print(list.plots[[4]][[i]], vp = vplayout(2, 2))
        print(list.plots[[5]][[i]], vp = vplayout(3, 1))
        print(list.plots[[6]][[i]], vp = vplayout(3, 2))
        dev.off(dev.cur())
    }
}

AVGCorrPlot <- function(input.file, y.axis, n.traits){
    require(ggplot2)
    data.corr = SetDataFrame(input.file, n.traits, T)
    data.avg = data.frame()
    module.name = unique(data.corr$module)
    for (i in 1:length(module.name)){
        aux.data = data.corr[data.corr$module==module.name[i],]
        aux.data = as.vector(tapply(aux.data$main, aux.data$generation, mean, module = module.name[i]))
        data.avg = rbind(data.avg, data.frame(generation = seq(data.corr$generation[1], data.corr$generation[length(data.corr$generation)]), main = aux.data, module = module.name[i]))
    }
    time.series <- ggplot(data.avg, aes(generation, main, group = module, color = module)) + layer (geom = "line") + scale_y_continuous(y.axis)
    return(time.series)
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

AVGRatioSinglePlot  <- function(input.file, n.traits){
    require(ggplot2)
    y.axis = "AVGRatio"
    data.corr = SetDataFrame(input.file, n.traits, T)
    AVGRatio <- AVGRatioCalc(input.file, n.traits)
    data.avg = data.frame(generation = seq(data.corr$generation[1], data.corr$generation[1]+length(AVGRatio)-1), main = AVGRatio)
    time.series <- ggplot(data.avg, aes(generation, main)) + layer (geom = "line") + scale_y_continuous(y.axis)
    return(time.series)
}

AVGRatioAVGPlot <- function(file.name, pattern = "Drift*", n.traits){
    require(ggplot2)
    y.axis = "AVG Ratio"
    folders  <- dir("output/", pattern)
    aux.file = paste("output", folders[1], file.name, sep="/")
    data.corr = SetDataFrame(aux.file, n.traits, T)
    generation.vector = seq(data.corr$generation[1], data.corr$generation[length(data.corr$generation)])
    n.gen = length(generation.vector)
    n.pop = length(folders)
    data.avg = array(dim=c(n.gen*n.pop, 2))
    for (pop in 1:(length(folders))){
        aux.file = paste("output", folders[pop], file.name, sep="/")
        AVGRatio <- AVGRatioCalc(aux.file, n.traits)
        lower = 1+((pop-1)*n.gen)
        upper = pop*n.gen
        print(lower)
        print(upper)
        data.avg[lower:upper,1] = generation.vector
        data.avg[lower:upper,2] = AVGRatio
    }
    data.avg = as.data.frame(data.avg)
    names(data.avg) = c("generation", "AVGRatio")
    time.series  <- ggplot(data.avg, aes(generation, AVGRatio)) + scale_y_continuous(y.axis) +
                    geom_point(alpha=1/500)+geom_smooth() + stat_smooth(geom="ribbon")
    return(time.series)
}

CorrOmegaCalc <- function(input.file, n.traits){
    data.corr = SetDataFrame(input.file, n.traits, T)
    omega = as.matrix(read.table ("input/omega.csv", header=F, sep=' '))[1:n.traits, 1:n.traits]
    omega = omega[upper.tri(omega)]
    corr.omega <- tapply(data.corr$main, data.corr$generation, function(x) cor(x, omega))
    return(corr.omega)
}

CorrOmegaSinglePlot  <- function(input.file, n.traits){
    require(ggplot2)
    y.axis = "Matrix Correlation with Omega"
    data.corr = SetDataFrame(input.file, n.traits, T)
    corr.omega <- CorrOmegaCalc(input.file, n.traits)
    data.avg = data.frame(generation = seq(data.corr$generation[1], data.corr$generation[1]+length(corr.omega)-1), main = corr.omega)
    time.series <- ggplot(data.avg, aes(generation, main)) + layer (geom = "line") + scale_y_continuous(y.axis)
    return(time.series)
}

CorrOmegaMultiPlot <- function(file.name, pattern = "DivSel*", n.traits, Label = F){
    require(ggplot2)
    y.axis = "Matrix Correlation with Omega"
    folders  <- dir("output/", pattern)
    aux.file = paste("output", folders[1], file.name, sep="/")
    data.corr = SetDataFrame(aux.file, n.traits, T)
    generation.vector = seq(data.corr$generation[1], data.corr$generation[length(data.corr$generation)])
    n.gen = length(generation.vector)
    n.pop = length(folders)
    data.avg = array(dim=c(n.gen*n.pop, 3))
    for (pop in 1:(length(folders))){
        aux.file = paste("output", folders[pop], file.name, sep="/")
        corr.omega <- CorrOmegaCalc(aux.file, n.traits)
        lower = 1+((pop-1)*n.gen)
        upper = pop*n.gen
        print(lower)
        print(upper)
        if (Label){label.vector = rep(folders[pop], n.gen)}
        else {
            aux.file.name = "pop.parameters.txt"
            aux.file = paste("output", folders[pop], aux.file.name, sep="/")
            parameters = scan(aux.file, character())
            index = which("theta"==parameters)+2
            label.vector = rep(as.numeric(parameters[index]), n.gen)
        }
        data.avg[lower:upper,1] = generation.vector
        data.avg[lower:upper,2] = corr.omega
        data.avg[lower:upper,3] = label.vector
    }
    data.avg = data.frame(as.numeric(data.avg[,1]), as.numeric(data.avg[,2]), data.avg[,3])
    names(data.avg) = c("generation", "corr.omega", "Selection_Strengh")
    time.series  <- ggplot(data.avg, aes(generation, corr.omega, group = Selection_Strengh, color=Selection_Strengh)) +
                    layer(geom = "smooth") + scale_y_continuous(y.axis)
    return(time.series)
}

PhenotipeMultiPlot <- function(pattern = "DivSel*", n.traits){
    require(ggplot2)
    file.name = "phenotype.dat"
    y.axis = "Mean phenotype norm"
    folders  <- dir("output/", pattern)
    aux.file = paste("output", folders[1], file.name, sep="/")
    data.phenotype = SetDataFrame(aux.file, n.traits, F)
    generation.vector = rep(seq(data.phenotype$generation[1], data.phenotype$generation[length(data.phenotype$generation)]), each=n.traits)
    n.gen = length(generation.vector)
    n.pop = length(folders)
    data.avg = array(dim=c(n.gen*n.pop, 3))
    for (pop in 1:(length(folders))){
        aux.file = paste("output", folders[pop], file.name, sep="/")
        ploted.data = abs(as.numeric(SetDataFrame(aux.file, n.traits, F)[,1]))
        lower = 1+((pop-1)*n.gen)
        upper = pop*n.gen
        print(lower)
        print(upper)
        aux.file.name = "pop.parameters.txt"
        aux.file = paste("output", folders[pop], aux.file.name, sep="/")
        parameters = scan(aux.file, character())
        index = which("theta"==parameters)+2
        label.vector = rep(as.numeric(parameters[index]), n.gen)
        data.avg[lower:upper,1] = generation.vector
        data.avg[lower:upper,2] = ploted.data
        data.avg[lower:upper,3] = label.vector
    }
    data.avg = data.frame(as.numeric(data.avg[,1]), as.numeric(data.avg[,2]), data.avg[,3])
    names(data.avg) = c("generation", "ploted.data", "Selection_Strengh")
    time.series  <- ggplot(data.avg, aes(generation, ploted.data, group = Selection_Strengh, color=Selection_Strengh)) +
                    layer(geom = "smooth") + scale_y_continuous(y.axis)
    return(time.series)
}
