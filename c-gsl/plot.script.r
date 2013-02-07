time.series.plot  <-  function(input.file, y.axis, n.traits, selection = T, corr.plot = F){
    require(ggplot2)
    data.init = read.table(input.file)
    if(corr.plot){
        aux.trait = rep(c(1, -1), each = n.traits/2)
        aux.trait = aux.trait%*%t(aux.trait)
        aux.trait = aux.trait[lower.tri(aux.trait)]
        aux.trait[aux.trait==1] = "within module"
        aux.trait[aux.trait==-1] = "between module"
        n.traits = (n.traits*n.traits-n.traits)/2
    }
    gen.number = data.init[seq(1,length(data.init[,1]),n.traits+1),]
    raw.trait.means = data.init[-seq(1,length(data.init[,1]),n.traits+1),]
    generations = length(gen.number)
    if(corr.plot){
        sel.scheme = rep(aux.trait, generations)
    }
    else{
        sel.scheme = rep(rep(c("positive", "negative"), each=n.traits/2), generations)
    }
    trait.name = as.factor(rep(1:n.traits))
    data.clean = data.frame(raw.trait.means, trait.name, gen.number, sel.scheme)
    names(data.clean) = c("main", "trait", "generation", "selection")
    if (selection)
        time.series <- ggplot(data.clean, aes(generation, main, group = trait, color = selection)) + layer (geom = "point") + scale_y_continuous(y.axis)
    else
        time.series <- ggplot(data.clean, aes(generation, main, group = trait, color = trait)) + layer (geom = "point") + scale_y_continuous(y.axis)
    return(time.series)
}
plot.pop  <- function (pop.path, n.traits){
    mean.phenotype.plot  <- time.series.plot (paste(pop.path, "phenotype.dat", sep = '/'), "mean phenotype", 10, F)
    g.var.plot <- time.series.plot (paste(pop.path, "g.var.dat", sep = '/'), "genetic variance", 10, T)
    p.var.plot <- time.series.plot (paste(pop.path, "p.var.dat", sep = '/'), "phenotypic variance", 10, F)
    h.var.plot <- time.series.plot (paste(pop.path, "h.var.dat", sep = '/'), "heritability", 10, F)
    g.corr.plot <- time.series.plot (paste(pop.path, "g.corr.dat", sep = '/'), "genetic correlations", 10, T, T)
    p.corr.plot <- time.series.plot (paste(pop.path, "p.corr.dat", sep = '/'), "phenotypic correlations", 10, T, T)
    plots = list (
          g.var = g.var.plot,
          p.var = p.var.plot,
          h.var = h.var.plot,
          mean.phenotype = mean.phenotype.plot,
          g.corr = g.corr.plot,
          p.corr = p.corr.plot
          )
    return (plots)
}
n.traits <- 10
pop.path <- "output/burn_in"
burnin.plots <- plot.pop(pop.path, n.traits)
