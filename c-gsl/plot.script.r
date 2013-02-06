time.series.plot  <-  function(input.file, y.axis, n.traits, selection = T){
    require(ggplot2)
    data.init = read.table(input.file)
    raw.trait.means = data.init[-seq(1,length(data.init[,1]),n.traits+1),]
    generations = length(raw.trait.means)/n.traits
    sel.scheme = rep(rep(c("positive", "negative"), each=n.traits/2), generations)
    gen.number = rep(1:generations, each=n.traits)
    trait.name = as.factor(rep(1:n.traits))
    data.clean = data.frame(raw.trait.means, trait.name, gen.number, sel.scheme)
    names(data.clean) = c("main", "trait", "generation", "selection")
    if (selection)
        time.series <- ggplot(data.clean, aes(generation, main, group = trait, color = selection)) + layer (geom = "line") + scale_y_continuous(y.axis)
    else
        time.series <- ggplot(data.clean, aes(generation, main, group = trait, color = trait)) + layer (geom = "line") + scale_y_continuous(y.axis)
    return(time.series)
}
n.traits <- 10
pop.path <- "output"
mean.phenotype.plot  <- time.series.plot (paste(pop.path, "phenotype.dat", sep = '/'), "mean phenotype", 10, T)
g.var.plot <- time.series.plot (paste(pop.path, "g.var.dat", sep = '/'), "genetic variance", 10, F)
p.var.plot <- time.series.plot (paste(pop.path, "p.var.dat", sep = '/'), "phenotipic variance", 10, F)
h.var.plot <- time.series.plot (paste(pop.path, "h.var.dat", sep = '/'), "heritability", 10, F)
