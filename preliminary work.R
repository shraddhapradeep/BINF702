library(BiocManager)
library(HMP16SData)
library(phyloseq)
library(magrittr)
library(ggplot2)
library(tibble)
library(dplyr)
library(dendextend)
library(circlize)
library(ExperimentHub)
library(gridExtra)
library(cowplot)
library(readr)
library(haven)
library(SummarizedExperiment)
library(vegan)
library(cluster)



V13<-V13()



v13metadata(V13()) # phyla tree, throwing an error

head(rownames(V13))# name of the otu 
head(colData(V13)) # more info for each participant 
head(rowData(V13),n=10) #otu and the different taxonomic annotations 



#Exploratory data analysis
#graph showing all the body sites collection in the package 


list(V13 = V13(), V35 = V35()) %>%
  table_one() %>%
  kable_one()


# pie chat of hmp body site of collection
v13_body_sites <- colData(v13)$HMP_BODY_SITE
v35_body_sites <- colData(v35)$HMP_BODY_SITE


#adding body sites together 

v13_v15_bodysites<-c(v13_body_sites,v35_body_sites)
v13_v15_bodysites<-na.omit(v13_v15_bodysites)

# frequency of bodysites 
body.freq<- table(v13_v15_bodysites)
view(body.freq)
#pie chat  of body sites
pie(body.freq,
    main = "Distribution of Body Sites",
    labels = paste(names(body.freq), ": ", body.freq),
    col = rainbow(length(body.freq)))
legend("topright", legend = names(body.freq), fill = rainbow(length(body.freq)))

##Statistical analysis Permonova or Anosim


V13_oral <-
  V13() %>%
  subset(select = HMP_BODY_SITE == "Oral")
V13_skin <-
  V13() %>%
  subset(select = HMP_BODY_SITE == "Skin")
V13_Urogenital<-
  V13() %>%
  subset(select = HMP_BODY_SITE == "Urogenital Tract")
V13_Airways <-
  V13() %>%
  subset(select = HMP_BODY_SITE == "Airways")
V13_Gastrointestinal<-
  V13() %>%
  subset(select = HMP_BODY_SITE == "Gastrointestinal Tract")


#Creating phyloseq objects
V13_oral_phyloseq <-
  as_phyloseq(V13_oral)

V13_skin_phyloseq <-
  as_phyloseq(V13_skin)

V13_Urogenital_phyloseq <-
  as_phyloseq(V13_Urogenital)

V13_Airways_phyloseq <-
  as_phyloseq(V13_Airways)

V13_Gastrointestinal_phyloseq <-
  as_phyloseq(V13_Gastrointestinal)

#used to subset because it will take a lot of computational power to analyze all the data
sample_samples <- function(x, size) {  
  sampled_names <-
    sample_names(x) %>%
    sample(size)
  
  prune_samples(sampled_names, x)
}

#sub-setting each body site phyloseq object is then sampled to contain only 50 samples
V13_oral_phyloseq %<>%
  sample_samples(25)

V13_skin_phyloseq %<>%
  sample_samples(25)

V13_Urogenital_phyloseq %<>%
  sample_samples(25)

V13_Airways_phyloseq %<>%
  sample_samples(25)

V13_Gastrointestinal_phyloseq %<>%
  sample_samples(25)




#Merging into one object
V13_phyloseq <-
  merge_phyloseq(V13_oral_phyloseq, V13_skin_phyloseq, V13_Urogenital_phyloseq, V13_Airways_phyloseq, V13_Gastrointestinal_phyloseq)

#filtering out Low abundant species data
V13_phyloseq %<>%
  taxa_sums() %>%
  is_greater_than(0) %>%
  prune_taxa(V13_phyloseq) # prune_taxa is a phyloseq object 

par(mfrow = c(1, 3))
plot_richness(V13_phyloseq, x='HMP_BODY_SITE', color = 'HMP_BODY_SITE', measures = 'Shannon')
plot_richness(V13_phyloseq, x='HMP_BODY_SITE', color = 'HMP_BODY_SITE', measures = 'Simpson')
plot_richness(V13_phyloseq, x='HMP_BODY_SITE', color = 'HMP_BODY_SITE', measures = 'Observed')

#Alpha diversity using the richness_measures in the phyloseq package 
richness_measures <-
  c("Observed", "Shannon", "Simpson")

#Creates boxplots of the alpha diverity 
V13_phyloseq %>%
  plot_richness(x = "HMP_BODY_SITE", color = "HMP_BODY_SITE", measures = richness_measures) +
  stat_boxplot(geom = "errorbar") +
  geom_boxplot() +
  scale_fill_manual(values = c("red", "blue", "green", "yellow", "purple"),  # Example colors
                    labels = c("Body Site 1", "Body Site 2", "Body Site 3", "Body Site 4", "Body Site 5")) +  # Example labels
  theme_bw() +
  theme(axis.title.x = element_blank())




class(V13_phyloseq)
dist_matrix<-distance(V13_phyloseq, method='bray')# using the Phyloseq distance matrix 
cluster<- hclust(dist_matrix) %>%
  as.dendrogram()
cluster_dendrogram <- as.dendrogram(cluster)

cluster_results <- hclust(dist_matrix)

#or  if phyloseq distance method doesn't work, 

# Calculate the Bray-Curtis dissimilarity matrix
#diss_matrix2 <- vegdist(otu_table(V13_phyloseq), method = "bray") #calculating distance using the distance function in the vegan package
#cluster2<-hclust(diss_matrix2)

#cluster_dendrogram2 <- as.dendrogram(cluster2)

V13_sample_data<- as.data.frame(sample_data(V13_phyloseq)) #sample data of v13 object as a data frame 

#denogram
#adding the col labels to the sample data frame for the dendrogram 
## Adding a color match to each body site in the sample data. color numbers can be gotten from https://www.color-hex.com/popular-colors.php
V13_sample_data$labels_col <- ifelse(V13_sample_data$HMP_BODY_SITE == "Oral", "#F8766D", 
                                     ifelse(V13_sample_data$HMP_BODY_SITE == "Gastrointestinal Tract", "#ff80ed",
                                            ifelse(V13_sample_data$HMP_BODY_SITE == "Airways", "#065535",
                                                   ifelse(V13_sample_data$HMP_BODY_SITE == "Urogenital Tract", "#ffd700",
                                                          ifelse(V13_sample_data$HMP_BODY_SITE == "Skin", "#00ffff", "#00BFC4")))))

plot(cluster_dendrogram, label = V13_sample_data$HMP_BODY_SITE, col = V13_sample_data$labels_col) # The Bray Curtis is much neater

# Get the OTU table or abundance data
otu_table <- as.matrix(otu_table(V13_phyloseq))
 

# Get the sample metadata
sample_data <- sample_data(V13_phyloseq)


# Convert dissimilarity matrix to a numeric matrix
numeric_dist_matrix <- as.matrix(dist_matrix)

# Plot heatmap of the numeric matrix
heatmap(numeric_dist_matrix)

# Convert into data frame of objects
V13_sample_data_df <- as.data.frame(V13_sample_data)

# Extract the data frame from the sample_data object
V13_sample_data_df <- V13_sample_data_df@.Data


V13_sample_data_df <- data.frame(V13_sample_data_df)


# Extract the list components
col1 <- V13_sample_data_df[[1]]
col2 <- V13_sample_data_df[[2]]
col3 <- V13_sample_data_df[[3]]
col4 <- V13_sample_data_df[[4]]
col5 <- V13_sample_data_df[[5]]
col6 <- V13_sample_data_df[[6]]
col7 <- V13_sample_data_df[[7]]
col8 <- V13_sample_data_df[[8]]

# Create a data frame with appropriate column names
V13_sample_data_df <- data.frame(
  RSID = col1,
  VISITNO = col2,
  SEX = col3,
  RUN_CENTER = col4,
  HMP_BODY_SITE = col5,
  HMP_BODY_SUBSITE = col6,
  SRS_SAMPLE_ID = col7,
  labels_col = col8
)




# Perform PERMANOVA analysis to examine whether the groupings were statistically significant
permanova_result <- adonis2(numeric_dist_matrix ~ as.factor(V13_sample_data_df$HMP_BODY_SITE), data = V13_sample_data_df, strata = NULL, permutations = 999) # this is if distance method from phyloseq works, else use the one below


summary(permanova_result)




#Principle Coordinates Analysis
##ordinate data is Phyloseq's method of scaling  data 
###PCoA is a method used to visualize and explore the similarity or dissimilarity of samples based on multivariate data, such as microbiome composition.
v13_ordinate<- ordinate(V13_phyloseq,'PCoA', distance="bray")
V13_phyloseq%>%
  plot_ordination(v13_ordinate,color = 'HMP_BODY_SITE', shape= 'HMP_BODY_SITE')+
  theme_bw()+
  theme(legend.position = 'bottom')

#silhouette score 
cluster_assignments <- cutree(cluster_results, 5) #cut into 5 for 5 groups 
silhouette_scores <- silhouette(cluster_assignments, numeric_dist_matrix)
summary(silhouette_scores)
#1: well matched, 0 is close to decision boundary, -1 is could be in wrong cluster
