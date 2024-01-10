# analyse anti-microbial resistance data from BV-BRC

# before running:
# - visit https://www.bv-brc.org/view/Taxonomy/1773#view_tab=genomes
# - find the bacterium of interest (e.g. Mycobacterium tuberculosis)
# - download the "Genomes" dataset in CSV format and save it to BVBRC_genome_[label].csv
# - download the "AMR Phenotypes" dataset in CSV format and save it to BVBRC_genome_amr_[label].csv
# - choose thresholds (see below)
# - the default setup assumes that you have downloaded tuberculosis data

library(countrycode)
library(pheatmap)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(Rcpp)
library(RcppArmadillo)
sourceCpp("hyperhmm-r.cpp")

source("hypercube-plots.R")

##### system-specific details
# records.threshold is the minimum value of records for a drug to be included in the analysis
# plot.thresh and plot.thresh.cs are the probability thresholds for plotting an edge in the graph visualisations for continents and countries respectively
species.label = "TB"; records.threshold = 100000; plot.thresh = 0.001; plot.thresh.cs = 0.001
#species.label = "Kp"; records.threshold = 400000; plot.thresh = 0.001; plot.thresh.cs = 0.0001

#################
## curation of data

# read dataframe with country of origin records for different genomes
df.1 = read.csv(paste0("BVBRC_genome_", species.label, ".csv"))
# read dataframe with (long-form) drug phenotypes for different genomes
df.2 = read.csv(paste0("BVBRC_genome_amr_", species.label, ".csv"))

# get unique IDs in both sets
uniq.1 = unique(df.1$Genome.ID)
uniq.2 = unique(df.2$Genome.ID)

# get unique IDs and drugs
uids = intersect(uniq.1, uniq.2)
drugs = unique(df.2$Antibiotic)
table(df.2$Antibiotic)

# slow and inefficient -- loop through unique IDs, pull the country from dataframe 1, and append to dataframe 2
df.2$Country = NA
for(i in 1:length(uids)) {
  refs.2 = which(df.2$Genome.ID == uids[i])
  ref.1 = which(df.1$Genome.ID == uids[i])
  df.2$Country[refs.2] = df.1$Isolation.Country[ref.1]
}

# also inefficient -- convert long-form to wide-form data
# initialise dataframe with just ID and country
df = data.frame(ID = df.2$Genome.ID, Country = df.2$Country)
drug.sums = data.frame(drug = drugs, count = 0)
# loop through drugs
for(drug in drugs) {
  # initialise set of phenotypes
  ic = rep(NA, nrow(df))
  # get rows displaying resistance and set our phenotype values accordingly
  res.refs = df.2$Genome.ID[which(df.2$Antibiotic == drug & (df.2$Resistant.Phenotype == "Resistant" | df.2$Resistant.Phenotype == "Intermediate"))]
  new.refs = which(df$ID %in% res.refs)
  ic[new.refs] = 1
  # get rows displaying suscpetibility and set our phenotype values accordingly
  sus.refs = df.2$Genome.ID[which(df.2$Antibiotic == drug & (df.2$Resistant.Phenotype == "Susceptible"))]
  new.refs = which(df$ID %in% sus.refs)
  ic[new.refs] = 0
  # add to dataframe (and count non-missing data)
  df[[drug]] = ic
  drug.sums$count[drug.sums$drug == drug] = length(which(!is.na(ic)))
}

# here are the set of drugs we looked at with TB in HyperTraPS/HyperHMM papers
# they are all in this set
drugs.interest = c("isoniazid", "rifampin",
                   "streptomycin", "ethambutol",
                   "pyrazinamide", "prothionamide",
                   "ofloxacin", "moxifloxacin",
                   "capreomycin", "amikacin")

# keep only those drugs with >400k non-missing entries
keep.cols = c(1, 2, 2+which(drug.sums$count > records.threshold))
final.df = df[,keep.cols]

# prune any rows with absent drug records
nas_in_each_row <- apply(final.df[,3:ncol(final.df)], 1, function(row) any(is.na(row)))
final.df = final.df[!nas_in_each_row,]

# get continent name from country name
final.df$Continent = countrycode(final.df$Country, "country.name", "continent")
clean.df = final.df[!is.na(final.df$Continent),]

table(clean.df$Continent)
write.csv(clean.df, paste0("bvbrc-clean-", species.label, ".csv"))

#################
## inference of pathways for continents

# go through set of continents
uniq.conts = unique(clean.df$Continent)
# initialise lists for data matrices and inference output
data.mat = list()
output.inf = list()
# loop through continents
for(i in 1:length(uniq.conts)) {
  print(uniq.conts[i])
  # get set of unique drug barcodes for this continent
  data.mat[[i]] = as.matrix(clean.df[clean.df$Continent==uniq.conts[i],3:(ncol(clean.df)-1)])
  rownames(data.mat[[i]]) = NULL
  data.mat[[i]] = unique(data.mat[[i]])
  # store a single, max lik HyperHMM run 
  output.inf[[i]] = HyperHMM(obs = data.mat[[i]], nboot = 0)
}

# make plots
plots = list()
for(i in 1:length(uniq.conts)) {
  plots[[i]] = plot.hypercube.flux(output.inf[[i]], thresh=plot.thresh) + ggtitle(uniq.conts[i])
}

# plot hypercubes for each continent
g.continents = ggarrange(plotlist = plots)
sf = 3
png(paste0("all-plots-continent-", species.label, ".png"), width=2000*sf, height=2000*sf, res=72*sf)
print(g.continents)
dev.off()

# output full sets of transitions 
pruned.df = data.frame()
for(i in 1:length(output.inf)) {
  tmp = output.inf[[i]]$transitions
  tmp$label = uniq.conts[i]
  pruned.df = rbind(pruned.df, tmp)
}

write.csv(pruned.df, paste0("transition-graphs-pruned-continent-", species.label, ".csv"))

################
## inference of pathways for countries

# go through set of countries
uniq.countries = unique(clean.df$Country)
# initialise lists for data matrices and inference output
data.mat.cs = list()
output.inf.cs = list()
# loop through continents
for(i in 1:length(uniq.countries)) {
  print(uniq.countries[i])
  # get set of unique drug barcodes for this continent
  data.mat.cs[[i]] = as.matrix(clean.df[clean.df$Country==uniq.countries[i],3:(ncol(clean.df)-1)])
  rownames(data.mat.cs[[i]]) = NULL
  data.mat.cs[[i]] = unique(data.mat.cs[[i]])
  # store a single, max lik HyperHMM run 
  output.inf.cs[[i]] = HyperHMM(obs = data.mat.cs[[i]], nboot = 0)
}

# make plots (for alphabetically-ordered countries)
ordered.cs = output.inf.cs[order(uniq.countries)]
ordered.labels = uniq.countries[order(uniq.countries)]
plots.cs = list()
for(i in 1:length(uniq.countries)) {
  plots.cs[[i]] = plot.hypercube.flux(ordered.cs[[i]], thresh=plot.thresh.cs) + ggtitle(ordered.labels[i])
}

# plot all inferred hypercubes by country
g.countries = ggarrange(plotlist = plots.cs)
sf = 3
png(paste0("all-plots-country-", species.label, ".png"), width=5000*sf, height=5000*sf, res=72*sf)
print(g.countries)
dev.off()

sf = 3
png(paste0("all-plots-country-smaller-", species.label, ".png"), width=2000*sf, height=2000*sf, res=72*sf)
print(g.countries)
dev.off()

plots.cs.nolegend = list()
for(i in 1:length(uniq.countries)) {
  plots.cs.nolegend[[i]] = plot.hypercube.flux(ordered.cs[[i]], thresh=plot.thresh.cs) + 
    ggtitle(ordered.labels[i]) + theme(legend.position = "none")
}

# output full sets of transitions 
pruned.df = data.frame()
for(i in 1:length(output.inf.cs)) {
  tmp = output.inf.cs[[i]]$transitions
  tmp$label = uniq.countries[i]
  pruned.df = rbind(pruned.df, tmp)
}

write.csv(pruned.df, paste0("transition-graphs-pruned-country-", species.label, ".csv"))

#################
### PCA for country-level parameterisations

# pull together parameterisations into a matrix (each row a parameterisation for a country)
all.sets = matrix(0, nrow=length(output.inf.cs), ncol=length(output.inf.cs[[1]]$transitions$Flux))
for(i in 1:length(uniq.countries)) {
  all.sets[i,] = output.inf.cs[[i]]$transitions$Flux
}

# do PCA and extract first two PCs (for some reason we must convert to numeric)
pca_result <- prcomp(all.sets, scale. = TRUE)
pca_components <- pca_result$x[, 1:2]
pca_data <- as.data.frame(cbind(PC1 = pca_components[, 1], 
                                PC2 = pca_components[, 2], 
                                labels = uniq.countries))
pca_data$PC1 = as.numeric(pca_data$PC1)
pca_data$PC2 = as.numeric(pca_data$PC2)

# plot these PCs
g.pca = ggplot(pca_data, aes(x=PC1, y=PC2, label=labels)) + geom_point() + 
  geom_text_repel(size=2, max.overlaps=30) + theme_light()
sf = 3
png(paste0("PCA-", species.label, ".png"), width=600*sf, height=400*sf, res=72*sf)
print(g.pca)
dev.off()

#####################
### MDS of l1 distances for country-level parameterisations

# create distance matrix of l1 norms
l1.mat = matrix(0, nrow=length(uniq.countries), ncol=length(uniq.countries))
for(i in 1:length(uniq.countries)) {
  for(j in 1:length(uniq.countries)) {
    l1.mat[i,j] = sum(abs(output.inf.cs[[i]]$transitions$Flux - output.inf.cs[[j]]$transitions$Flux))
  }
}

# do MDS embedding
mds_result <- cmdscale(l1.mat, k = 2)
mds.df = data.frame(label=uniq.countries,
                    x = mds_result[,1],
                    y = mds_result[,2])
# plot MDS embedding
g.mds = ggplot(mds.df, aes(x=x,y=y,label=label)) + geom_point() + 
  geom_text_repel(size=2, max.overlaps=20) + 
  xlab("MDS coord 1") + ylab("MDS coord 2") + theme_light()
sf = 3
png(paste0("MDS-", species.label, ".png"), width=600*sf, height=400*sf, res=72*sf)
print(g.mds)
dev.off()

png(paste0("MDS-PCA-", species.label, ".png"), width=1000*sf, height=400*sf, res=72*sf)
print(ggarrange(g.mds, g.pca))
dev.off()

# heatmap visualisation of l1 distances
to.hm = l1.mat
row.names(to.hm) = uniq.countries

png(paste0("l1-dists-", species.label, ".png"), width=800*sf, height=800*sf, res=72*sf)
pheatmap(to.hm)
dev.off()

########################
### compute full WFCC curves

# a function to return a dataframe containing the WFCC curve
WFCC = function(g.1, g.2) {
  # set of t values to consider
  fVals = sort(unique(c(0,g.1,g.2,1)))
  wfcc = data.frame(fVals = fVals, WFCC = 0)
  # loop through t values
  for(i in 1:nrow(wfcc)) {
    print(i)
    # loop through edges
    for(j in 1:length(g.1)) {
      # increment WFCC score in the interval of disjointedness for this edge
      w.lo = min(g.1[j], g.2[j])
      w.hi = max(g.1[j], g.2[j])
      if(wfcc$fVals[i] >= w.lo & wfcc$fVals[i] < w.hi) {
        wfcc$WFCC[i] = wfcc$WFCC[i]+1
      }
    } 
  }
  return(wfcc)
}

# particular case study: compare Hong Kong (basically uniform/null) with a set of other countries
# with different evolutionary pathways (from visual inspection)
res.df = data.frame()
ref.name = "Hong Kong"
test.names = c("Kenya", "Sweden", "Norway", "Bulgaria", "Italy", "Spain", "Belgium")
# loop over comparisons
for(i in 1:length(test.names)) {
  this.test = test.names[i]
  print(c(this.test, ref.name))
  g.1 = output.inf.cs[[which(uniq.countries == ref.name)]]$transitions$Flux
  g.2 = output.inf.cs[[which(uniq.countries == this.test)]]$transitions$Flux
  this.wfcc = WFCC(g.1, g.2)
  this.wfcc$ref = ref.name
  this.wfcc$test = this.test
  res.df = rbind(res.df, this.wfcc)
}

# plot WFCC in log space
g.wfcc.log = ggplot(res.df[res.df$fVals > 0,], aes(x=fVals, y=WFCC, color=factor(test, levels=test.names))) + 
  geom_step(linewidth=3, alpha=0.5, direction="hv") + 
  scale_y_log10() + theme_light() 

# plot zoomed-in WFCC in lin space
g.wfcc.lin = ggplot(res.df[res.df$fVals >= 0,], aes(x=fVals, y=WFCC, color=factor(test, levels=test.names))) + 
  geom_step(linewidth=1, alpha=1, direction="hv") + coord_cartesian(xlim = c(-0.001,0.025)) +
  theme_light()

sf = 3
png(paste0("WFCC-", species.label, ".png"), width=1000*sf, height=500*sf, res=72*sf)
print(ggarrange(g.wfcc.log, g.wfcc.lin))
dev.off()

# interrogate some particular instances.
# most single-trait cases have streptomycin as the first drug, but not all
one.hits = c("Brazil", "Colombia", "Ecuador", "Finland", 
             "Ghana", "Oman", "South Korea", "Sweden", 
             "Turkmenistan", "Viet Nam")
data.mat.cs[which(uniq.countries %in% one.hits)]




