library(reshape2)
library(ggplot2)
library(EvomodR)

AVGWrap = function(p.cor, n.e, generations = 1:length(p.cor)){
  burn.in.avg = AVGRatio(p.cor, 40/10000, F, num.cores = 4, generations = generations)
  burn.in.avg['Probability'] = NULL
  burn.in.avg['Selection_Strength'] = NULL
  m.avg = melt(burn.in.avg[,-c(2,5)], id.vars = c('.id', 'generation'))
  m.avg = m.avg[!((m.avg['.id'] != "Full Integration") & (m.avg['variable'] == "AVG-")),]
  m.avg = m.avg[!((m.avg['.id'] == "Full Integration") & (m.avg['variable'] == "AVG+")),]
  m.avg['Population_Size'] = n.e
  return(m.avg)
}

NeDataFrame = function(pop.list, generations = 1:length(p.cor)){
  m.avg = ldply(pop.list, 
                function(x) AVGWrap(x$p.cor[generations], 
                                    x$n.e, generations), 
                .progress = 'text'
                )
  return(m.avg)
}

NeFacetPlot = function(m.avg){
  avg.plot = ggplot(m.avg, aes(generation,
                               value,
                               group=interaction(variable, generation, .id),
                               colour=interaction(.id, variable))) +
    layer(geom="point") +
    labs(x="Generation",
         y="Average Correlation",
         color = "Module") +
    scale_colour_discrete(labels=c("Within Module 1",
                                   "Within Module 2",
                                   "Between Modules")) + 
    facet_wrap( ~ Population_Size, ncol = 4) + theme_bw()
  return(avg.plot)
}

