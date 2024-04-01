library(BiocManager)
library(phyloseq)


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



v13<-V13()
v35<-V35()
v13

metadata(V13()) # phylo tree

head(rownames(V13())) # name of the otu 
head(colData(v13)) # more info for each participant 
head(rowData(V13()),n=10) #otu and the different taxonomic anotations 
V13[1,1]
 

#Exploratory data analysis
#graph showing all the body sites collection in the package 


list(V13 = v13(), v35 = V35()) %>%
  table_one() %>%
  kable_one()


# pie chat of hmp body site of collection
v13_body_sites <- colData(v13)$HMP_BODY_SITE
v35_body_sites <- colData(v35)$HMP_BODY_SITE


#adding bodysites together 
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

#Statistical analysis Permonova or Anosim

view(v13_body_sites)



#Given the heterogeneity in patient age and region, we performed an analysis of similarities (ANOSIM) to examine whether the groupings were statistically significant
V35_stool <-
  V35() %>%
  subset(select = HMP_BODY_SUBSITE == "Stool")
colData(V35_stool)
rowData(V35_stool)




