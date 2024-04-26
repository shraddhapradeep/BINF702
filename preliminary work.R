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




v13<-V13()







#Exploratory data analysis
#graph showing all the body sites collection in the package 

list(v13 = V13()) %>%
  table_one() %>%
  kable_one()


# pie chat of hmp body site of collection
v13_body_sites <- colData(v13)$HMP_BODY_SITE

#adding body sites together 

v13_v15_bodysites<-(v13_body_sites)
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
  prune_taxa(V13_phyloseq)# prune_taxa is a phyloseq object 

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


#BEta diversity analysis 
dist_matrix<- phyloseq::distance(V13_phyloseq, method='bray')# using the Phyloseq distance matrix 

#Hierarchical Clustering
cluster<- hclust(dist_matrix) 




V13_sample_data<- as.data.frame(sample_data(V13_phyloseq)) #sample data of v13 object as a data frame 

#denogram
#adding the col labels to the sample data frame for the dendrogram 
## Adding a color match to each body site in the sample data. color numbers can be gotten from https://www.color-hex.com/popular-colors.php
V13_sample_data$labels_col <- ifelse(V13_sample_data$HMP_BODY_SITE == "Oral", "#F8766D", 
                                     ifelse(V13_sample_data$HMP_BODY_SITE == "Gastrointestinal Tract", "#ff80ed",
                                            ifelse(V13_sample_data$HMP_BODY_SITE == "Airways", "#065535",
                                                   ifelse(V13_sample_data$HMP_BODY_SITE == "Urogenital Tract", "#ffd700",
                                                          ifelse(V13_sample_data$HMP_BODY_SITE == "Skin", "#00ffff", "#00BFC4")))))



# Create a subset of the dendrogram for better visualization
subset_dendrogram <- cutree(cluster, k = 5)
# Plot dendrogram with colored labels
plot(cluster, main = "Hierarchical Clustering Dendrogram",
     labels = V13_sample_data$HMP_BODY_SITE, col = V13_sample_data$labels_col)




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


dev.off()

#Principle Coordinates Analysis
##ordinate data is Phyloseq's method of scaling  data 
###PCoA is a method used to visualize and explore the similarity or dissimilarity of samples based on multivariate data, such as microbiome composition.
v13_ordinate<- ordinate(V13_phyloseq,'PCoA', distance="bray", scale=TRUE)
V13_phyloseq%>%
  plot_ordination(v13_ordinate,color = 'HMP_BODY_SITE', shape= 'HMP_BODY_SITE')+
  theme_bw()+
  theme(legend.position = 'bottom')



summary(v13_ordinate)




# Extract eigenvalues from the v13_ordinate object
eigenvalues <- v13_ordinate$values

# Plot the eigenvalues
barplot(eigenvalues[, 1], main = "Eigenvalues", xlab = "Axis", ylab = "Eigenvalue", names.arg = 1:nrow(eigenvalues))

#silhouette score 
cluster_assignments <- cutree(cluster, 5) #cut into 5 for 5 groups 
silhouette_scores <- silhouette(cluster_assignments, numeric_dist_matrix)
summary(silhouette_scores)
#1: well matched, 0 is close to decision boundary, -1 is could be in wrong cluster


#Prediction using random forests model 

df <- psmelt(V13_phyloseq)%>%
  na.omit(df)


#Normalied abundance 
df_normalized <- df %>%
  group_by(SRS_SAMPLE_ID) %>%
  mutate(`Relative Abundance` = Abundance / sum(Abundance)) %>%
  ungroup()
    

df_grouped <- df %>%
  group_by(HMP_BODY_SITE, PHYLUM) %>%
  summarize(total_abundance = sum(Abundance)) %>%
  ungroup()

# Rank the genera within each body site based on total abundance
df_ranked <- df_grouped %>%
  group_by(HMP_BODY_SITE) %>%
  mutate(rank = dense_rank(desc(total_abundance))) %>%
  ungroup()

# Select the top 10 genera for each body site Select the top 10 genera for each body site and rename the rank column
df_top10 <- df_ranked %>%
  filter(rank <= 10) %>%
  select(HMP_BODY_SITE, PHYLUM, total_abundance,rank)

# Optionally, you can arrange the data frame by body site and rank for better visualization
df_top10 <- df_top10 %>%
  arrange(HMP_BODY_SITE, rank)

#color code 
bangcolors <- c("#CC79A7", "#D55E00", "#0072B2", "#F0E442", "#009E73", "#56B4E9",
                "#E69F00", "#000000", "#CC002E", "#6A3D9A", "#FFD966", "#B2DF8A",
                "#66C2A5", "#FF7F00", "#8C564B", "#1F78B4", "#FDD0A2", "#E41A1C",
                "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628", "#999999",
                "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
                "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#A1DAB4")


# Convert to data frame
df_top10 <- as.data.frame(df_top10)
ggplot(data=df_top10) +
  geom_bar(aes(x=HMP_BODY_SITE, y=total_abundance, fill=PHYLUM), stat="identity")+
  scale_fill_manual(values=bangcolors) +
  theme_minimal() +
  labs(title = "Top 10 Genera by Abundance in Each Body Site")
  


         

#factor categories 
df_normalized <- df_normalized %>% 
  mutate(HMP_BODY_SITE = factor(HMP_BODY_SITE))%>% 
  mutate(SEX= factor(SEX)) %>% 
  mutate(CLASS= factor(CLASS))%>%
  mutate(PHYLUM= factor(PHYLUM))%>%
  mutate(ORDER= factor(ORDER))%>%
  filter(Abundance > 0) #filter out low abundant species 


df_normalized<- as.data.frame(df_normalized) # convert into a data frame 

library(randomForest)


#selecting columns for analysis 
v13_2 <- dplyr::select(df_normalized, HMP_BODY_SITE,OTU,Abundance,CLASS,
                       PHYLUM, 
                       ORDER, 
                       FAMILY, 
                       GENUS)
set.seed(0)
train<- sample(1:nrow(v13_2),0.8*nrow(v13_2)) # selecting 80% random rows for training 
V13_train<-(v13_2[train,])
V13_test<-(v13_2[-train,])

rfmodel <- randomForest(HMP_BODY_SITE ~ ., data = V13_train,ntree = 500, importance = TRUE, proximity = TRUE )
print(rfmodel)

# Making predictions on the training set
v13_train_pred <- predict(rfmodel, V13_train)

# Checking the accuracy on the training set
confusion_matrix<-table(v13_train_pred, V13_train$HMP_BODY_SITE)
accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
print(accuracy)


# Making predictions on the test set
v13_test_pred <- predict(rfmodel, V13_test)

# Checking the accuracy
confusion_matrix2<-table(v13_test_pred, V13_test$HMP_BODY_SITE)
accuracy_test <- sum(diag(confusion_matrix2)) / sum(confusion_matrix2)

# Print and plot the variable-importance measures
importance(rfmodel)

# Open a PDF device
pdf("Variable_Importance_Plottt.pdf", width = 10, height = 8)

# Plot variable importance
print(varImpPlot(rfmodel, n.var = 5, main = "5 Variable Importance"))

# Close the PDF device
dev.off()



