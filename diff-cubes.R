# need these libraries for plotting, but not generation
library(ggraph)
library(igraph)
library(ggpubr)

# given a data frame with From, To, P[robability] columns, work out the fluxes through each edge
# i.e. P*source node probability
# very inefficient and won't scale!
fluxes = function(df) {
  # initial flux column
  df$F = 0
  # probability of occupancy of each node -- initially unit mass at 0
  stateprobs = rep(0, 8)
  stateprobs[0+1] = 1
  # take three steps across the cube
  for(step in 1:3) {
    for(i in 1:nrow(df)) {
      # find rows where flux is zero and source node probability is nonzero, and populate, updating node probabilities
      if(df$F[i] == 0 & stateprobs[df$From[i]+1] != 0) {
        flux = df$P[i]*stateprobs[df$From[i]+1]
        df$F[i] = flux
        stateprobs[df$To[i]+1] = stateprobs[df$To[i]+1]+flux
      }
    }
  }
  return(df)
}

# given a data frame with From, To, P[robability] columns, normalise P so that outgoing edges from each node sum to 1
normalise = function(df) {
  for(i in 0:4) {
    norm = sum(df$P[which(df$From==i)])
    df$P[which(df$From==i)] = df$P[which(df$From==i)]/norm
  }
  return(df)
}

# ggraph visualisation for a given cube

cubic.graph = function(df) {
  g.df = graph_from_data_frame(df)
  return(ggraph(g.df, layout="sugiyama") + geom_edge_link(aes(edge_width=F, edge_alpha=F)) + 
    geom_node_point() + geom_node_label(aes(label=name),size=3) +
    scale_edge_width(limits=c(0,1)) + scale_edge_alpha(limits=c(0,1)) +
    theme_graph()+ggtitle(df$label) )
}
#change F to P above to get figures of the probability weightings instead of the flux weightings

# manually set up the L=3 hypercube and initialise data structures
froms = c(0, 0, 0, 1, 1, 2, 2, 4, 4, 3, 5, 6)
tos =   c(1, 2, 4, 3, 5, 3, 6, 5, 6, 7, 7, 7)
df.template = data.frame(From=froms, To=tos, P=NA)
df.collection = data.frame()
g.collection = list()

## the following examples pull the template dataframe, assign transition probabilities reflecting a given model, compute fluxes, add to a growing collection, and plot to a growing list
## uniform null case
df = df.template
df$P = 1
df = fluxes(normalise(df))
df$label = "Null"
g.collection[[length(g.collection)+1]] = cubic.graph(df)
df.collection = rbind(df.collection, df)

## single pathway
df = df.template
df$P = c(1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1)
df = fluxes(normalise(df))
df$label = "One path"
g.collection[[length(g.collection)+1]] = cubic.graph(df)
df.collection = rbind(df.collection, df)

## single pathway, some noise
df = df.template
df$P = c(1, 0.1, 0.1, 1, 0.1, 1, 1, 1, 1, 1, 1, 1)
df = fluxes(normalise(df))
df$label = "One path-small-noise"
g.collection[[length(g.collection)+1]] = cubic.graph(df)
df.collection = rbind(df.collection, df)

## single pathway, more noise
df = df.template
df$P = c(1, 0.5, 0.5, 1, 0.5, 1, 1, 1, 1, 1, 1, 1)
df = fluxes(normalise(df))
df$label = "One path-more-noise"
g.collection[[length(g.collection)+1]] = cubic.graph(df)
df.collection = rbind(df.collection, df)

## double pathway
df = df.template
df$P = c(1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1)
df = fluxes(normalise(df))
df$label = "Two paths"
g.collection[[length(g.collection)+1]] = cubic.graph(df)
df.collection = rbind(df.collection, df)

## double pathway, some noise
df = df.template
df$P = c(1, 0.1, 1, 1, 0.1, 1, 1, 0.1, 1, 1, 1, 1)
df = fluxes(normalise(df))
df$label = "Two paths-small-noise"
g.collection[[length(g.collection)+1]] = cubic.graph(df)
df.collection = rbind(df.collection, df)

## double pathway, more noise
df = df.template
df$P = c(1, 0.5, 1, 1, 0.5, 1, 1, 0.5, 1, 1, 1, 1)
df = fluxes(normalise(df))
df$label = "Two paths-more-noise"
g.collection[[length(g.collection)+1]] = cubic.graph(df)
df.collection = rbind(df.collection, df)

## single then branch
df = df.template
df$P = c(1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1)
df = fluxes(normalise(df))
df$label = "Branch"
g.collection[[length(g.collection)+1]] = cubic.graph(df)
df.collection = rbind(df.collection, df)

## single then branch, some noise
df = df.template
df$P = c(1, 0.1, 0.1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
df = fluxes(normalise(df))
df$label = "Branch-small-noise"
g.collection[[length(g.collection)+1]] = cubic.graph(df)
df.collection = rbind(df.collection, df)

## single then branch, more noise
df = df.template
df$P = c(1, 0.5, 0.5, 1, 1, 1, 1, 1, 1, 1, 1, 1)
df = fluxes(normalise(df))
df$label = "Branch-more-noise"
g.collection[[length(g.collection)+1]] = cubic.graph(df)
df.collection = rbind(df.collection, df)

## other single
df = df.template
df$P = c(0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1)
df = fluxes(normalise(df))
df$label = "One path2"
g.collection[[length(g.collection)+1]] = cubic.graph(df)
df.collection = rbind(df.collection, df)

## other single, some noise
df = df.template
df$P = c(0.1, 0.1, 1, 1, 1, 1, 1, 0.1, 1, 1, 1, 1)
df = fluxes(normalise(df))
df$label = "One path2-small-noise"
g.collection[[length(g.collection)+1]] = cubic.graph(df)
df.collection = rbind(df.collection, df)

## other single, more noise
df = df.template
df$P = c(0.5, 0.5, 1, 1, 1, 1, 1, 0.5, 1, 1, 1, 1)
df = fluxes(normalise(df))
df$label = "One path2-more-noise"
g.collection[[length(g.collection)+1]] = cubic.graph(df)
df.collection = rbind(df.collection, df)

## other single then branch
df = df.template
df$P = c(0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
df = fluxes(normalise(df))
df$label = "Branch2"
g.collection[[length(g.collection)+1]] = cubic.graph(df)
df.collection = rbind(df.collection, df)

## other single then branch, some noise
df = df.template
df$P = c(0.1, 0.1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
df = fluxes(normalise(df))
df$label = "Branch2-small-noise"
g.collection[[length(g.collection)+1]] = cubic.graph(df)
df.collection = rbind(df.collection, df)

## other single then branch, more noise
df = df.template
df$P = c(0.5, 0.5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
df = fluxes(normalise(df))
df$label = "Branch2-more-noise"
g.collection[[length(g.collection)+1]] = cubic.graph(df)
df.collection = rbind(df.collection, df)

# output plot
png("./Data_files_hypercubes/Outputs_github/diff-cubes.png", width=800*2, height=800*2, res=72*2)
print(ggarrange(plotlist=g.collection))
dev.off()

# output cubes
write.csv(df.collection, "./Data_files_hypercubes/Outputs_github/diff-cubes.csv", row.names=FALSE, quote=FALSE)
