library(ggplot2)
#Synthetic hypercubes for fluxes
df = read.csv("./Data_files_hypercubes/Distance_matrix_flux.csv")
dist_matrix = as.matrix(df[1:nrow(df), 2:ncol(df)])
mds_result <- cmdscale(dist_matrix, k = 2)
mds.df = data.frame(label=colnames(df)[2:ncol(df)],
                    x = mds_result[,1],
                    y = mds_result[,2])
svg("./Data_files_hypercubes/Outputs_github/MDS_flux.svg")
ggplot(mds.df, aes(x=x,y=y,label=label)) + geom_point() + geom_text(nudge_y=-0.05) + 
  xlab("MDS coord 1") + ylab("MDS coord 2") + theme_light()
dev.off()

#Synthetic hypercubes for probabilities
df = read.csv("./Data_files_hypercubes/Distance_matrix_probability.csv")
dist_matrix = as.matrix(df[1:nrow(df), 2:ncol(df)])
mds_result <- cmdscale(dist_matrix, k = 2)
mds.df = data.frame(label=colnames(df)[2:ncol(df)],
                    x = mds_result[,1],
                    y = mds_result[,2])
svg("./Data_files_hypercubes/Outputs_github/MDS_probability.svg")
ggplot(mds.df, aes(x=x,y=y,label=label)) + geom_point() + geom_text(nudge_y=-0.05) + 
  xlab("MDS coord 1") + ylab("MDS coord 2") + theme_light()
dev.off()

#Klebsiella for fluxes
df = read.csv("./Data_files_hypercubes/distance_matrix_kelbsiella_f.csv")
dist_matrix = as.matrix(df[1:nrow(df), 2:ncol(df)])
mds_result <- cmdscale(dist_matrix, k = 2)
mds.df = data.frame(label=colnames(df)[2:ncol(df)],
                    x = mds_result[,1],
                    y = mds_result[,2])
svg("./Data_files_hypercubes/Outputs_github/klebsiella_f_MDS.svg")
ggplot(mds.df, aes(x=x,y=y,label=label)) + geom_point() + geom_text(nudge_y=-0.05) + 
  xlab("MDS coord 1") + ylab("MDS coord 2") + theme_light()
dev.off()

#Klebsiella for probabilities
df = read.csv("./Data_files_hypercubes/distance_matrix_kelbsiella_p.csv")
dist_matrix = as.matrix(df[1:nrow(df), 2:ncol(df)])
mds_result <- cmdscale(dist_matrix, k = 2)
mds.df = data.frame(label=colnames(df)[2:ncol(df)],
                    x = mds_result[,1],
                    y = mds_result[,2])
svg("./Data_files_hypercubes/Outputs_github/klebsiella_p_MDS.svg")
ggplot(mds.df, aes(x=x,y=y,label=label)) + geom_point() + geom_text(nudge_y=-0.05) + 
  xlab("MDS coord 1") + ylab("MDS coord 2") + theme_light()
dev.off()
