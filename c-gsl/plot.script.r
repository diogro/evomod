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
    module = rep(aux.trait, generations)
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
        time.series <- ggplot(data.clean, aes(generation, main, group = trait, color = module)) + layer (geom = "point") + scale_y_continuous(y.axis)
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
    g.avg.corr.plot <- AVGRatioPlot (paste(pop.path, "g.corr.dat", sep = '/'), "mean genetic correlations", n.traits)
    p.corr.plot <- time.series.plot (paste(pop.path, "p.corr.dat", sep = '/'), "phenotypic correlations", n.traits, T, T)
    p.avg.corr.plot <- AVGRatioPlot (paste(pop.path, "p.corr.dat", sep = '/'), "mean phenotypic correlations", n.traits)
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

PlotPng  <-  function(list.plots, file.name){
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

AVGRatioPlot <- function(input.file, y.axis, n.traits){
    require(ggplot2)
    data.corr = SetDataFrame(input.file, n.traits, T)
    data.avg = data.frame()
    module.name = unique(data.corr$module)
    for (i in 1:length(module.name)){
        aux.data = data.corr[data.corr$module==module.name[i],]
        aux.data = as.vector(tapply(aux.data$main, aux.data$generation, mean, module = module.name[i]))
        data.avg = rbind(data.avg, data.frame(generation = seq(data.corr$generation[1], data.corr$generation[length(data.corr$generation)]), main = aux.data, module = module.name[i]))
    }
    time.series <- ggplot(data.avg, aes(generation, main, group = module, color = module)) + layer (geom = "point") + scale_y_continuous(y.axis)
    return(time.series)
}

n.traits <- 10
pop.path <- "output/burn_in"
burnin.plots <- PlotPop(pop.path, n.traits)

sel.strengths  <- seq(20, 200, 30)/10000
div.folders = paste("DivSel-", sel.strengths, sep='')
div.plots = vector('list', length(sel.strengths))
for (i in 1:length(sel.strengths)){
    pop.folder = paste("output/", div.folders[i], sep = '')
    print(pop.folder)
    div.plots[[i]]  <- PlotPop(pop.folder, n.traits)
}
names(div.plots) = div.folders
corridor.folders = paste("CoridorSel-", sel.strengths, sep='')
corridor.plots = vector('list', length(sel.strengths))
for (i in 1:length(sel.strengths)){
    pop.folder = paste("output/", corridor.folders[i], sep = '')
    print(pop.folder)
    corridor.plots[[i]]  <- PlotPop(pop.folder, n.traits)
}
names(corridor.plots) = corridor.folders

PlotPng  (corridor.plots, "corridor")
PlotPng  (div.plots, "divergent")
