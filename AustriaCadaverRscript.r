#AMS Cadaver Script Jan 2020
#J Receveur

###########
#Import
############
#Function from https://gist.github.com/grabear/018e86413b19b62a6bb8e72a9adba349
#Parse Silva function
parse_taxonomy_silva_128 <- function(char.vec){
  # Use default to assign names to elements in case problem with greengenes prefix
  char.vec = parse_taxonomy_default(char.vec)
  # Check for unassigned taxa
  if (char.vec["Rank1"] == "Unassigned") {
    char.vec <- c(Rank1="D_0__Unassigned", Rank2="D_1__Unassigned", Rank3="D_2__Unassigned", Rank4="D_3__Unassigned",
                  Rank5="D_4__Unassigned", Rank6="D_5__Unassigned", Rank7="D_6__Unassigned")
  }
  # Define the meaning of each prefix according to SILVA taxonomy
  Tranks = c(D_0="Kingdom", D_1="Phylum", D_2="Class", D_3="Order", D_4="Family", D_5="Genus", D_6="Species")
  # Check for prefix using regexp, warn if there were none. trim indices, ti
  ti = grep("[[:alpha:]]\\_[[:digit:]]{1}\\_\\_", char.vec)
  if( length(ti) == 0L ){
    warning(
      "No silva prefixes were found. \n",
      "Consider using parse_taxonomy_delfault() instead if true for all OTUs. \n",
      "Dummy ranks may be included among taxonomic ranks now."
    )
    # Will want to return without further modifying char.vec
    taxvec = char.vec
    # Replace names of taxvec according to prefix, if any present...
  } else {
    # Format character vectors for Ambiguous taxa
    if( length(ti) < 7 ){
      for (key in names(char.vec)) {
        if ( char.vec[key] == "Ambiguous_taxa" ) {
          tax_no <- (as.numeric(substr(key, 5, 5)) - 1)
          char.vec[key] = sprintf("D_%s__Ambiguous_taxa", tax_no)
        }
      }
      # Reset the trimmed indicies if Ambiguous taxa
      ti = grep("[[:alpha:]]\\_[[:digit:]]{1}\\_\\_", char.vec)
    }
    # Remove prefix using sub-"" regexp, call result taxvec
    taxvec = gsub("[[:alpha:]]\\_[[:digit:]]{1}\\_\\_", "", char.vec)
    # Define the ranks that will be replaced
    repranks = Tranks[substr(char.vec[ti], 1, 3)]
    # Replace, being sure to avoid prefixes notK present in Tranks
    names(taxvec)[ti[!is.na(repranks)]] = repranks[!is.na(repranks)]
  }
  return(taxvec)
}



library(vegan)
library(MASS)
library(ggplot2)
library(plyr)
library(dplyr)
library(magrittr)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)
library(randomForest)
library(knitr)
library(ape)
library(ggpubr)
library(Rmisc)
library(multcompView)
set.seed(10)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#000000","#CC79A7")
theme_set(theme_bw(base_size = 12)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
biom=import_biom("AustriaCadaverJan2020.biom",parseFunction=parse_taxonomy_silva_128)
biom
tax_table(biom)
#tax_table(biom) <- tax_table(biom)[,-c(5:10,14)]#remove dummy ranks

metadata=read.csv("AustriaBodyFarmMetadataWDivJan2020.csv",header = TRUE)
metadata$Timepoint<-factor(metadata$Timepoint, levels = c("4-Aug","1-Sep","17-Nov"))
metadata$Timepoint2<-factor(metadata$Timepoint2, levels = c("Aug","Sep","Nov"))

metadata$Timepoint2
metadata$Combined<-paste0(metadata$Donor," ",metadata$SampleType)
ddply(metadata, c("Donor","SampleType"), summarise, N    = length(id))




head(metadata)
tree=read_tree("AustriaCadaverTree.nwk")

sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$id
physeq=merge_phyloseq(biom,sampdat,tree)
physeq


head(tax_table(physeq))
#plot_heatmap(physeq)

head(metadata)
###############
#Taxonomy
################
GPr  = transform_sample_counts(physeq, function(x) x / sum(x) ) #transform samples based on relative abundance
PhylumAll=tax_glom(GPr, "Phylum")



PhylumLevel = filter_taxa(PhylumAll, function(x) mean(x) > 3e-3, TRUE) #filter out any taxa lower tha 0.3%
FamilyAll=tax_glom(GPr,"Family")
FamilyLevel = filter_taxa(FamilyAll, function(x) mean(x) > 3e-2, TRUE) #filter out any taxa lower tha 3%
GenusAll=tax_glom(GPr,"Genus")
GenusLevel = filter_taxa(GenusAll, function(x) mean(x) > 1e-2, TRUE) #filter out any taxa lower tha 1%

df <- psmelt(PhylumAll)
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Phylum"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
Trtdata

df <- psmelt(FamilyLevel)
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Family"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
Trtdata

df <- psmelt(GenusLevel)
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Genus"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
Trtdata



#############3
#Heatmaps
############

GPr2  = transform_sample_counts(physeq, function(x) x / sum(x)*100 ) #transform samples based on relative abundance
PhylumAll2=tax_glom(GPr2, "Phylum")


Heatmap<-plot_heatmap(PhylumAll2,"MDS","unifrac","SampleType","Phylum",sample.order= "SampleType")+facet_wrap(Donor~Timepoint2,scale="free_x")

theme_set(theme_bw(base_size = 10)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

#heatmap(otu_table(PhylumAll))
Heatmap
dev.off()
tiff("Figures/HeatmapPhylumPlot0.3.tiff", width = 174, height = 220, units = 'mm', res = 1000)
Heatmap
dev.off()

Phylum2Subset<-subset_samples(PhylumAll2, SampleType!="soil")

Heatmap2<-plot_heatmap(Phylum2Subset,"MDS","unifrac","SampleType","Phylum",sample.order= "SampleType")+facet_wrap(Donor~Timepoint2,scale="free_x")

Heatmap2
dev.off()
tiff("Figures/HeatmapPhylumPlot2.tiff", width = 174, height = 220, units = 'mm', res = 1000)
Heatmap2
dev.off()


GPr2  = transform_sample_counts(physeq, function(x) x / sum(x)*100 ) #transform samples based on relative abundance
FamilyAll2=tax_glom(GPr2, "Family")
HeatmapFamily<-plot_heatmap(FamilyAll2,"MDS","unifrac","SampleType","Family",sample.order= "SampleType")+facet_wrap(Donor~Timepoint2,scale="free_x")
HeatmapFamily

dev.off()
tiff("Figures/HeatmapFamilyPlot.tiff", width = 174, height = 220, units = 'mm', res = 1000)
HeatmapFamily
dev.off()

FamilyAllSubset<-subset_samples(FamilyAll2,SampleType!="soil")
HeatmapFamily2<-plot_heatmap(FamilyAllSubset,"MDS","unifrac","SampleType","Family",sample.order= "SampleType")+facet_wrap(Donor~Timepoint2,scale="free_x")
HeatmapFamily2

dev.off()
tiff("Figures/HeatmapFamilyNoSoil.tiff", width = 174, height = 220, units = 'mm', res = 1000)
HeatmapFamily2
dev.off()


############################
##################3Heatmaps subset taxa
##############################
GPr  = transform_sample_counts(physeq, function(x) x / sum(x) *100) #transform samples based on relative abundance
PhylumAll=tax_glom(GPr, "Phylum")



PhylumLevel2 = filter_taxa(PhylumAll, function(x) mean(x) > 3e-1, TRUE) #filter out any taxa lower tha 1%


Heatmap<-plot_heatmap(PhylumLevel2,"MDS","unifrac","SampleType","Phylum",sample.order= "SampleType")+facet_wrap(Donor~Timepoint2,scale="free_x")

theme_set(theme_bw(base_size = 10)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

#heatmap(otu_table(PhylumAll))
Heatmap
dev.off()
tiff("Figures/HeatmapPhylumPlot0.3.tiff", width = 174, height = 220, units = 'mm', res = 1000)
Heatmap
dev.off()

Phylum2Subset<-subset_samples(PhylumLevel2, SampleType!="soil")

Heatmap2<-plot_heatmap(Phylum2Subset,"MDS","unifrac","SampleType","Phylum",sample.order= "SampleType")+facet_wrap(Donor~Timepoint2,scale="free_x")

Heatmap2
dev.off()
tiff("Figures/HeatmapPhylumNoSoil0.3.tiff", width = 174, height = 220, units = 'mm', res = 1000)
Heatmap2
dev.off()


GPr2  = transform_sample_counts(physeq, function(x) x / sum(x)*100 ) #transform samples based on relative abundance
FamilyAll2=tax_glom(GPr2, "Family")
FamilyLevel2 = filter_taxa(FamilyAll2, function(x) mean(x) > 3, TRUE) #filter out any taxa lower tha 3%

HeatmapFamily<-plot_heatmap(FamilyLevel2,"MDS","unifrac","SampleType","Family",sample.order= "SampleType")+facet_wrap(Donor~Timepoint2,scale="free_x")
HeatmapFamily

dev.off()
tiff("Figures/HeatmapFamilyPlot3Per.tiff", width = 174, height = 220, units = 'mm', res = 1000)
HeatmapFamily
dev.off()








###########
#

GPr2  = transform_sample_counts(physeq, function(x) x / sum(x)*100 ) #transform samples based on relative abundance
FamilyAll2=tax_glom(GPr2, "Family")
FamilyLevel2 = filter_taxa(FamilyAll2, function(x) mean(x) > 0.3, TRUE) #filter out any taxa lower tha 3%

HeatmapFamily<-plot_heatmap(FamilyLevel2,"MDS","unifrac","SampleType","Family",sample.order= "SampleType")+facet_wrap(Donor~Timepoint2,scale="free_x")
HeatmapFamily

dev.off()
tiff("Figures/HeatmapFamilyPlot0.3Per.tiff", width = 250, height = 400, units = 'mm', res = 1000)
HeatmapFamily
dev.off()

GPr2  = transform_sample_counts(physeq, function(x) x / sum(x)*100 ) #transform samples based on relative abundance
FamilyAll2=tax_glom(GPr2, "Family")
FamilyLevel3 = filter_taxa(FamilyAll2, function(x) mean(x) > 0.03, TRUE) #filter out any taxa lower tha 3%

FamilyAllSubset<-subset_samples(FamilyLevel3,SampleType!="soil")
HeatmapFamily2<-plot_heatmap(FamilyAllSubset,"MDS","unifrac","SampleType","Family",sample.order= "SampleType")+facet_wrap(Donor~Timepoint2,scale="free_x")
#HeatmapFamily2$scales$scales[[2]]$name <- "Family"
#print(HeatmapFamily2)
HeatmapFamily2

dev.off()
tiff("Figures/HeatmapFamily0.03NoSoil.tiff", width = 174, height = 400, units = 'mm', res = 600)
HeatmapFamily2
dev.off()

#################
GPr2  = transform_sample_counts(physeq, function(x) x / sum(x)*100 ) #transform samples based on relative abundance
FamilyAll2=tax_glom(GPr2, "Family")
FamilyLevel2 = filter_taxa(FamilyAll2, function(x) mean(x) > 0.03, TRUE) #filter out any taxa lower tha 3%

HeatmapFamily<-plot_heatmap(FamilyLevel2,"MDS","unifrac","SampleType","Family",sample.order= "SampleType")+facet_wrap(Donor~Timepoint2,scale="free_x")
HeatmapFamily

dev.off()
tiff("Figures/HeatmapFamilyPlot.03Per.tiff", width = 250, height = 400, units = 'mm', res = 1000)
HeatmapFamily
dev.off()

GPr2  = transform_sample_counts(physeq, function(x) x / sum(x)*100 ) #transform samples based on relative abundance
FamilyAll2=tax_glom(GPr2, "Family")
FamilyLevel3 = filter_taxa(FamilyAll2, function(x) mean(x) > 0.03, TRUE) #filter out any taxa lower tha 3%

FamilyAllSubset<-subset_samples(FamilyLevel3,SampleType!="soil")
HeatmapFamily2<-plot_heatmap(FamilyAllSubset,"MDS","unifrac","SampleType","Family",sample.order= "SampleType")+facet_wrap(Donor~Timepoint2,scale="free_x")
#HeatmapFamily2$scales$scales[[2]]$name <- "Family"
#print(HeatmapFamily2)
HeatmapFamily2

dev.off()
tiff("Figures/HeatmapFamily0.03NoSoil.tiff", width = 174, height = 400, units = 'mm', res = 600)
HeatmapFamily2
dev.off()



################################3

GPr2  = transform_sample_counts(physeq, function(x) x / sum(x)*100 ) #transform samples based on relative abundance
GenusAll2=tax_glom(GPr2, "Genus")
GenusLevel2 = filter_taxa(GenusAll2, function(x) mean(x) > 1, TRUE) #filter out any taxa lower tha 3%

HeatmapGenus<-plot_heatmap(GenusLevel2,"MDS","unifrac","SampleType","Genus",sample.order= "SampleType")+facet_wrap(Donor~Timepoint2,scale="free_x")
HeatmapGenus

dev.off()
tiff("Figures/HeatmapGenusPlot1Per.tiff", width = 174, height = 220, units = 'mm', res = 1000)
HeatmapGenus
dev.off()

GenusAllSubset<-subset_samples(GenusLevel2,SampleType!="soil")
HeatmapGenus2<-plot_heatmap(GenusAllSubset,"MDS","unifrac","SampleType","Genus",sample.order= "SampleType")+facet_wrap(Donor~Timepoint2,scale="free_x")
HeatmapGenus2

dev.off()
tiff("Figures/HeatmapGenusNoSoil1Per.tiff", width = 174, height = 220, units = 'mm', res = 1000)
HeatmapGenus2
dev.off()
################
#Phylum level plot
df <- psmelt(PhylumLevel)
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Phylum"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
Trtdata
#compare_means(Abundance~Timepoint, data=df,group.by="Phylum",method = "kruskal.test",p.adjust.method="fdr")
#Phylum level
df <- psmelt(PhylumLevel)
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Phylum","SampleType","Donor","Timepoint2"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance)
)
head(Trtdata)
write.csv(Trtdata,"PhylumLevelWOOther.csv")
Trtdata2<-read.csv("PhylumLevelWithOther.csv",header=T)
Trtdata2$Timepoint2<-factor(Trtdata2$Timepoint2, levels = c("Aug","Sep","Nov"))

PhylumPlot=ggplot(Trtdata2, aes(x=SampleType,y=mean))+geom_bar(aes(fill = Phylum),colour="black", stat="identity")+
  facet_grid(~Phylum)+xlab("Sample Location")+ylab("Relative Abundance (%, > 0.3%)") + theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.title.x=element_blank())+facet_grid(Donor~Timepoint2)  +scale_fill_viridis_d()
#scale_x_discrete(labels=c("0 hrs", "24 hrs", "48 hrs","72 hrs"))+ theme(axis.text.x = element_text(angle = 45, hjust = 1))#+ scale_fill_manual(values=cbPalette)
PhylumPlot

theme_set(theme_bw(base_size = 13)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

dev.off()
tiff("Figures/PhylumPlot.tiff", width = 174, height = 174, units = 'mm', res = 1200)
PhylumPlot
dev.off()


#PhylumPlot$
#Phylum level by sample type


PhylumPlot2=ggplot(Trtdata2, aes(x=Timepoint2,y=mean))+geom_bar(aes(fill = Phylum),colour="black", stat="identity")+
  facet_grid(~Phylum)+xlab("Timepoint")+ylab("Relative Abundance (%, > 0.3%)") + theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.title.x=element_blank())+facet_grid(Case~SampleType)+scale_fill_viridis_d()
#scale_x_discrete(labels=c("0 hrs", "24 hrs", "48 hrs","72 hrs"))+ theme(axis.text.x = element_text(angle = 45, hjust = 1))#+ scale_fill_manual(values=cbPalette)
PhylumPlot2
theme_set(theme_bw(base_size = 10)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

dev.off()
tiff("Figures/PhylumPlot2.tiff", width = 174, height = 174, units = 'mm', res = 1200)
PhylumPlot2
dev.off()



theme_set(theme_bw(base_size = 13)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

#PhylumBySampleTypeSubset
TrtdataSubset<-subset(Trtdata,SampleType=="mouth"|SampleType=="nose"|SampleType=="rectum"|SampleType=="ears"|SampleType=="eyes")
PhylumPlot3=ggplot(TrtdataSubset, aes(x=Timepoint2,y=mean))+geom_bar(aes(fill = Phylum),colour="black", stat="identity")+
  facet_grid(~Phylum)+xlab("Timepoint")+ylab("Relative Abundance (%, > 0.3%)") + theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.title.x=element_blank())+facet_grid(Case~SampleType)+  scale_fill_manual(values=cbPalette)
#scale_x_discrete(labels=c("0 hrs", "24 hrs", "48 hrs","72 hrs"))+ theme(axis.text.x = element_text(angle = 45, hjust = 1))#+ scale_fill_manual(values=cbPalette)
PhylumPlot3


dev.off()
tiff("Figures/PhylumPlot3.tiff", width = 174, height = 174, units = 'mm', res = 1200)
PhylumPlot3
dev.off()


df <- psmelt(PhylumLevel)
df$Abundance=df$Abundance*100
head(df)
Trtdata <- ddply(df, c("Phylum","SampleType","Timepoint2"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
head(Trtdata)
write.csv(Trtdata,"PhylumLevelWOOtherDonorMerged.csv")
#################################3
df <- psmelt(FamilyLevel)
df$Abundance=df$Abundance*100
head(df)
Trtdata <- ddply(df, c("Family","SampleType","Timepoint2"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
head(Trtdata)
write.csv(Trtdata,"FamilyLevelWOOtherDonorMerged.csv")
############################
df <- psmelt(GenusLevel)
df$Abundance=df$Abundance*100
head(df)
Trtdata <- ddply(df, c("Genus","SampleType","Timepoint2"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
head(Trtdata)
write.csv(Trtdata,"GenusLevelWOOtherDonorMerged.csv")

###################
#Family level
################
dfAll<-psmelt(FamilyAll)
levels(unique(dfAll$Family))
df <- psmelt(FamilyLevel)
df$Abundance=df$Abundance*100


Trtdata <- ddply(df, c("Family","SampleType","Case","Timepoint2"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance)
)
Trtdata
write.csv(Trtdata,"FamilyLevelWOOther.csv")
Trtdata2<-read.csv("FamilyLevelWithOther.csv",header=T)
Trtdata2$Timepoint2<-factor(Trtdata2$Timepoint2, levels = c("Aug","Sep","Nov"))

FamilyPlot=ggplot(Trtdata2, aes(x=Timepoint2,y=mean))+geom_bar(aes(fill = Family),colour="black", stat="identity")+
  facet_grid(~Family)+xlab("Sample Type")+ylab("Relative Abundance (%, > 3%)") + theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.title.x=element_blank())+facet_grid(Case~SampleType)+scale_fill_viridis_d()#+  scale_fill_manual(values=cbPalette)
#scale_x_discrete(labels=c("0 hrs", "24 hrs", "48 hrs","72 hrs"))+ theme(axis.text.x = element_text(angle = 45, hjust = 1))#+ scale_fill_manual(values=cbPalette)
FamilyPlot

theme_set(theme_bw(base_size = 9)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

dev.off()
tiff("Figures/FamilyPlot.tiff", width = 174, height = 174, units = 'mm', res = 1200)
FamilyPlot
dev.off()

Trtdata <- ddply(df, c("Family","SampleType"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
head(Trtdata)
write.csv(Trtdata,"FamilyLevelWOOtherTimepointAndDonorMerged.csv")



df <- psmelt(GenusLevel)
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Genus","SampleType"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
head(Trtdata)
write.csv(Trtdata,"GenusLevelWOOtherTimepointsAndDonorMerged.csv")

#PCoA
StatSubset<-subset_samples(physeq,SampleType=="mouth"|SampleType=="nose"|SampleType=="ears")

GPdist=phyloseq::distance(StatSubset, "wunifrac")
beta=betadisper(GPdist, sample_data(StatSubset)$SampleType)
permutest(beta)
boxplot(beta)


adonis(GPdist ~ SampleType, as(sample_data(StatSubset), "data.frame"))
StatSubset<-subset_samples(physeq,SampleType!="soil")

ord=ordinate(StatSubset,"PCoA", "wunifrac")
ordplot=plot_ordination(StatSubset, ord,"samples", color="SampleType",shape="Timepoint2")+geom_point(size=4)+  scale_color_manual(values=cbPalette)
ordplot
theme_set(theme_bw(base_size = 12)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

dev.off()
tiff("Figures/PCoANoSoil.tiff", width = 174, height = 174, units = 'mm', res = 1200)
ordplot
dev.off()





###############
#Random Forest
##################
#All samples
GenusAll
GenusAllSubset2<-subset_taxa(GenusAll, Genus != "NA")
GenusAllSubset2
ForestData=GenusAll#Change this one so you dont have to rewrite all variables
predictors=t(otu_table(ForestData))
response <- as.factor(sample_data(ForestData)$Timepoint)
rf.data <- data.frame(response, predictors)
MozzieForest <- randomForest(response~., data = rf.data, ntree = 1000)
print(MozzieForest)#returns overall Random Forest results
imp <- importance(MozzieForest)#all the steps that are imp or imp. are building a dataframe that contains info about the taxa used by the Random Forest testto classify treatment 
imp <- data.frame(predictors = rownames(imp), imp)
imp.sort <- arrange(imp, desc(MeanDecreaseGini))
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
imp.20 <- imp.sort[1:20, ]#22

#write.csv(imp.20,"Top20GenusLevelPredictorsAllTimepoint.csv")

ImpWGenus<-read.csv("Top20GenusLevelPredictorsAllTimepoint.csv",header=T)
ImpWGenus$FamilyGenus<-paste0(ImpWGenus$Family,": ",ImpWGenus$Genus)
ImpWGenus<-arrange(ImpWGenus, desc(MeanDecreaseGini))
ImpWGenus$Genus <- factor(ImpWGenus$Genus, levels = ImpWGenus$Genus)
ggplot(ImpWGenus, aes(x = Genus, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip()


#Stat subset

GenusSubset<-subset_samples(GenusAll,SampleType=="mouth"|SampleType=="nose"|SampleType=="ears")
GenusAllSubset2<-subset_taxa(GenusSubset, Genus != "NA")

ForestData=GenusSubset#Change this one so you dont have to rewrite all variables
predictors=t(otu_table(ForestData))
response <- as.factor(sample_data(ForestData)$SampleType)
rf.data <- data.frame(response, predictors)
MozzieForest <- randomForest(response~., data = rf.data, ntree = 1000)
print(MozzieForest)#returns overall Random Forest results
imp <- importance(MozzieForest)#all the steps that are imp or imp. are building a dataframe that contains info about the taxa used by the Random Forest testto classify treatment 
imp <- data.frame(predictors = rownames(imp), imp)
imp.sort <- arrange(imp, desc(MeanDecreaseGini))
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)

imp.20 <- imp.sort[1:20, ]
imp.20
#write.csv(imp.20,"Top20GenusLevelPredictorsSubsetTimepoint.csv")
ImpWGenus<-read.csv("Top20GenusLevelPredictorsSubsetTimepoint.csv",header=T)
ImpWGenus$FamilyGenus<-paste0(ImpWGenus$Family,": ",ImpWGenus$Genus)
ImpWGenus<-arrange(ImpWGenus, desc(MeanDecreaseGini))
ImpWGenus$Genus <- factor(ImpWGenus$Genus, levels = ImpWGenus$Genus)

ggplot(ImpWGenus, aes(x = Genus, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip()


df <- psmelt(GenusSubset)
df<- subset(df,Genus %in% ImpWGenus$Genus[1:10])
df

df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Genus","SampleType","Case","Timepoint2"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance)
)
cdataplot=ggplot(Trtdata, aes(x=SampleType,y=mean))+geom_bar(aes(fill = Genus),colour="black", stat="identity")+ facet_grid(Case~Timepoint2)+xlab("Sample Type")+ylab("Relative Abundance (%)") + theme(axis.text.x = element_text(angle = 45, hjust = 1))+theme(axis.title.x=element_blank())+scale_fill_viridis_d()
cdataplot

############
#Richness
###########
FaithAll<-ggplot(metadata, aes(x=Timepoint2,y=faith_pd,shape=Donor,group=Donor,color=Donor))+geom_point(size=2)+geom_line()+
  facet_grid(~Phylum)+xlab("Sample Location")+ylab("Faith's Phylogenetic Diversity") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  theme(axis.title.x=element_blank(),axis.text.x = element_text(size = 7))+facet_wrap(~ SampleType)+  scale_fill_manual(values=cbPalette)#+scale_x_discrete(labels=c("8/4" = "Aug 4","9/1"= "Sep 1","11/17"="Nov 17"))
FaithAll
theme_set(theme_bw(base_size = 11)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

dev.off()
tiff("Figures/Faith'sDivAll.tiff", width = 84, height = 84, units = 'mm', res = 1200)
FaithAll
dev.off()


Faiths<-ggplot(Subset, aes(x=Timepoint2,y=faith_pd,shape=Donor,group=Donor,color=Donor))+geom_point(size=2)+geom_line()+
  facet_grid(~Phylum)+xlab("Sample Location")+ylab("Faith's Phylogenetic Diversity") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  theme(axis.title.x=element_blank(),axis.text.x = element_text(size = 7))+facet_wrap(~ SampleType)+  scale_fill_manual(values=cbPalette)#+scale_x_discrete(labels=c("8/4" = "Aug 4","9/1"= "Sep 1","11/17"="Nov 17"))
Faiths
dev.off()
tiff("Figures/Faith'sDivSubset.tiff", width = 84, height = 84, units = 'mm', res = 1200)
Faiths
dev.off()



Subset<-subset(metadata,SampleType=="mouth"|SampleType=="nose"|SampleType=="rectum"|SampleType=="ears")
Faiths<-ggplot(Subset, aes(x=Timepoint2,y=faith_pd,shape=Donor,group=Donor,color=Donor))+geom_point(size=2)+geom_line()+
  facet_grid(~Phylum)+xlab("Sample Location")+ylab("Faith's Phylogenetic Diversity") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  theme(axis.title.x=element_blank(),axis.text.x = element_text(size = 7))+facet_wrap(~ SampleType)+  scale_fill_manual(values=cbPalette)#+scale_x_discrete(labels=c("8/4" = "Aug 4","9/1"= "Sep 1","11/17"="Nov 17"))
Faiths
dev.off()
tiff("Figures/Faith'sDivSubset.tiff", width = 84, height = 84, units = 'mm', res = 1200)
Faiths
dev.off()



Subset<-subset(metadata,SampleType=="mouth"|SampleType=="nose"|SampleType=="rectum"|SampleType=="ears")
Shannon<-ggplot(Subset, aes(x=Timepoint,y=shannon,shape=Donor,group=Donor,color=Donor))+geom_point(size=2)+geom_line()+
  facet_grid(~Phylum)+xlab("Sample Location")+ylab("Shannon Diversity") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  theme(axis.title.x=element_blank(),axis.text.x = element_text(size = 7))+facet_wrap(~SampleType)+  scale_fill_manual(values=cbPalette)+scale_x_discrete(labels=c("8/4" = "Aug 4","9/1"= "Sep 1","11/17"="Nov 17"))

dev.off()
tiff("Figures/ShannonDivSubset.tiff", width = 84, height = 84, units = 'mm', res = 1200)
Shannon
dev.off()

dev.off()
tiff("Figures/DiversityCombinedSubset.tiff", width = 174, height = 84, units = 'mm', res = 1200)
ggarrange(Faiths,Shannon,
          labels = c("A", "B"), nrow = 1)
dev.off()


InvSimpson<-plot_richness(physeq, x="Timepoint2",color="Donor",shape="Donor", measures=c("InvSimpson"))+ylab("Inverse Simpson")+facet_wrap(~SampleType)+geom_point()+scale_fill_manual(cbPalette)+geom_point(size=5)
metadata$InvSimpson<-InvSimpson$data$value
Observed<-plot_richness(physeq, x="Timepoint2",color="Donor",shape="Donor", measures=c("Observed"))+facet_wrap(~SampleType)+geom_point()+scale_fill_manual(cbPalette)+geom_point(size=5)
metadata$ObservedSpecies<-Observed$data$value
head(metadata)
write.csv(metadata,"AustriaCadaverMetadataWDiversityMetrics.csv")
head(metadata)
