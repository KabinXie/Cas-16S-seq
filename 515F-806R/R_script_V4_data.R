# Cas-16S-seq
# Compare the data of 515F-806R (V4) amplicon and three gRNA were used.
#
# Summarizing the OTU table 
# Kabin Xie
# 2020.1.30

###################################
## The script include following 
## A. Load the data 
## B. Process the soil samples
## C. Process the leaf samples
###################################

library(ggplot2)
library(cowplot)
library(gridExtra)
library(phyloseq)
library(vegan)

########################################################################################
# A. Load the data    
#   1. load the data
#   2. loade the sequences and build the phyloseq object
########################################################################################

## A1. load the otu table, sample information, raw-OTU processing summaries  
# remove the the # symbol at the header line before load the otu table
# Pay attention to col names if it start with digit. An X is added if the colname/rowname 
# started with digit.

setwd("~/mystation/microbiome/Cas-16S-seq/SLY20191224_leaf_soil_sample/V4_leaf_soil_sample/analyze")
otutab<-read.table(file="../SLY_results/all.otutab100.V4.txt", head=T, check.names=FALSE)
head(otutab)
rownames(otutab)<-otutab[,1]

# In our analysis pipeline, OTUs (clustered by 100%) which are aligned to rice mito 16S rRNA
# with no more than 2 mismatches are considered as rice rDNA (all.otu100.mt2.txt)
OTU.mt<-read.table(file="../SLY_results/all.otu100.v4.mtOTU.txt", head=T)
OTU.cp<-read.table(file="../SLY_results/all.otu100.v4.cpOTU.txt", head=T)
table(OTU.mt$OTU_ID %in% OTU.cp$OTU_ID)
table(!otutab$OTU_ID %in% OTU.mt$OTU_ID)
table(!otutab$OTU_ID %in% OTU.cp$OTU_ID)

otu.matrix<-as.matrix(otutab[!otutab$OTU_ID %in% OTU.mt$OTU_ID,-1])  #-1, remove first name column

head(otu.matrix)
dim(otu.matrix)
dim(otutab)

# get the percentage of rice mt-rDNA for each sample
mito.no<-colSums(otutab[otutab$OTU_ID %in% OTU.mt$OTU_ID, -1])
mito.ratio<-mito.no/apply(otutab[,-1],2, sum)
mito.ratio

chlo.no<-colSums(otutab[otutab$OTU_ID %in% OTU.cp$OTU_ID, -1])
chlo.ratio<-chlo.no/apply(otutab[,-1],2, sum)
chlo.ratio

bact.no<-colSums(otutab[!otutab$OTU_ID %in% 
                          c(as.character(OTU.cp$OTU_ID), as.character(OTU.mt$OTU_ID)),-1])
bact.ratio<-bact.no/apply(otutab[,-1],2, sum)

chlo_mito_bact_ratio<-rbind(mito.ratio,chlo.ratio,bact.ratio)
colSums(chlo_mito_bact_ratio)
barplot(chlo_mito_bact_ratio[,7:12])
chlo_mito_bact_no<-rbind(mito.no,chlo.no,bact.no)

## load the sample meta data table
sample.info<-read.table("../SLY_results/vsearch_sample_info_v4.csv", head=T, sep=",")
rownames(sample.info)<-sample.info$Seq_ID
chlo_mito_bact_ratio<-t(chlo_mito_bact_ratio)
as.character(sample.info$Seq_ID)==as.character(rownames(chlo_mito_bact_ratio)) # check the order
sample.info<-cbind(sample.info, chlo_mito_bact_ratio[rownames(sample.info),])

## read NGS-OTU processing data (vsearch-pipline summary) ##
vs.summary<-read.table(file="../SLY_results/V4_vsearch.summary.csv", head=TRUE, sep=",")
head(vs.summary)
rownames(vs.summary)<-vs.summary$Seq_ID
sample.info<-cbind(sample.info, vs.summary[rownames(sample.info), c(3,4,5)])
OTU_No<-data.frame(OTU_No=colSums(otutab[,-1]>0))

sample.info<-cbind(sample.info, OTU_No=OTU_No[rownames(sample.info),])
sample.info<-sample.info[with(sample.info, order(cp.gRNA, Compartment,Cas9, Cycle_No)),]
sample.info$Total<-colSums(otutab[,-1])[sample.info$Plant_sample]
write.table(sample.info, file="../SLY_results/sample.info.v4.txt", sep="\t")
rm(vs.summary, OTU_No)

sample.info.soil<-sample.info[sample.info$Compartment == "Soil",]
sample.info.leaf<-sample.info[sample.info$Compartment == "leaf",]
## Here is an example to examine the data
## number of dereplication sequences per samples.
qplot(interaction(Cas9, Compartment, Cycle_No), derep/merged, data=sample.info, geom=c("boxplot", "jitter"))
p.reads<-qplot(interaction(Cas9, Compartment, Cycle_No), 100*derep/filtered, 
               fill = Cas9, data=sample.info)
p.reads<-p.reads+geom_boxplot();
#p.reads<-p.reads+geom_jitter();
p.reads<-p.reads+theme(text=element_text(size=16))
p.reads<-p.reads+theme_grey();
p.reads
ggsave(filename = "../SLY_results/derep_ratio.pdf", p.reads,
       width = 10, height = 10, 
       units = "cm")


##########################################################
#  A2. Build the phyloseq object                         #
#    assign taxanomy and assemble the phyloseq object    #
##########################################################
library(phyloseq)
library(gridExtra)
library(DECIPHER)
library(dada2)
library(Biostrings)
library(ShortRead)
#get OTU sequences and assign taxanomy
otu.seq<-readDNAStringSet("../SLY_results/all.otus100.v4.fasta", format="fasta")
table((width(otu.seq)))

#taxonomy assign taken more than 4 hours with 7 cores.
set.seed(197)
#taxa <- assignTaxonomy(otu.seq, "~/mystation/microbiome/silva_nr_v132_train_set.fa.gz", multithread=TRUE, minBoot=60)
#save the workspace here, RStudio may crash here. (set n to read n seq each time, save memory)
#saveRDS(taxa, file="../SLY_results/v4_taxa_assign.RData")
taxa<-readRDS("../SLY_results/v4_taxa_assign.RData")  # Load the taxa information
#taxa.specis<-addSpecies(taxa, "~/mystation/microbiome/silva_species_assignment_v132.fa.gz", n=1000, verbose=F)
#saveRDS(taxa.specis, file="species_assign.RData")

## rowname of taxa is the OTU sequence which is named by fasta header ##
taxa.rownames.matrix<-matrix(unlist(strsplit(names(rownames(taxa)), ";")), ncol=2, byrow=T)
rownames(taxa)<-taxa.rownames.matrix[,1]

otu.seq.bact<-otu.seq[names(otu.seq)[taxa[,5]!="Mitochondria"]]
writeXStringSet(otu.seq.bact, filepath = "../SLY_results/otu100.bact.fasta")

rm(otu.seq.bact)

ps <- phyloseq(otu_table(otutab[,-1], taxa_are_rows=T), 
               sample_data(sample.info), 
               tax_table(taxa)) 

# subset root and soil samples #
ps.leaf<-phyloseq(otu_table(otutab[, as.character(sample.info.leaf$Seq_ID)], taxa_are_rows = T),
                  sample_data(sample.info.leaf),
                  tax_table(taxa))
ps.leaf<-prune_taxa(rowSums(otu_table(ps.leaf))>0, ps.leaf )
ps.soil<-phyloseq(otu_table(otutab[,as.character(sample.info.soil$Seq_ID)],taxa_are_rows = T), 
                  sample_data(sample.info.soil),
                  tax_table(taxa))
ps.soil<-prune_taxa(rowSums(otu_table(ps.soil))>0, ps.soil)

sample_sums(prune_taxa(rowSums(otu_table(ps.soil))>0,ps.soil))
ps_richness<-estimate_richness(ps, measures =c("Shannon","Simpson"))
sample.info<-cbind(sample.info, estimate_richness(ps, measures =c("Shannon","Simpson"))[as.character(sample.info$Seq_ID),])

#remove the rice Mt sequences
ps.noMt<-prune_taxa(!rownames(otu_table(ps)) %in% OTU.mt$OTU_ID, ps)
# Examine other mitochondria OTUs ()
un_id_Mt<-prune_taxa(data.frame(tax_table(ps.noMt))$Family == "Mitochondria",ps.noMt)
sample_sums(un_id_Mt);
taxa_sums(un_id_Mt) # They are low abundant (<10).
rm(un_id_Mt)

ps.noM.noCh<-prune_taxa(!rownames(otu_table(ps.noMt)) %in% OTU.cp$OTU_ID, ps.noMt)
ps_richness<-estimate_richness(ps.noMt, measures =c("Shannon","Simpson"))
colnames(ps_richness)<-c("Shannon_NoMtCp","Simpson_NoMtCp")
sample.info<-cbind(sample.info, ps_richness[as.character(sample.info$Seq_ID),])
write.table(sample.info, file="../SLY_results/sample.info.with.alpha.txt", sep="\t")

#########################################################
##  B. Process the soil samples                        ##
##  v4_OTUs in soil samples 
##  1. Rarefaction curve to examine OTU diversity      ##
##  2. Normalize the soil data by rarefaction
##  3.OTUs detected in Cas+ and/or Cas9- data (Barplot) 
##  4. Bacterial Diversities comparison 
##  5. statistical testing of diversity metrics
##  6. Scatterplot of OTU abundance
##  7. Identify Differential OTUs (Deseq2)
#########################################################
ps.soil

##### B1. Rarefication curve for soil samples #########
rarefy.observed<-data.frame()
rarefy.otu<-data.frame()
sampleSums(ps.soil) # min: 32288
set.seed(99)
sample_size=50000
for( i in seq(1, sample_size, by = 500)) {
  rarefy.ps.soil<-rarefy_even_depth(ps.soil,i, replace=F)
  # abundance threshold 10, TUs less than this threshold are NOT counted
  if(i>50000) {
    rarefy.ps.soil<-prune_taxa(taxa_sums(rarefy.ps.soil)>10,rarefy.ps.soil)
  }
  
  rarefy.richness<-estimate_richness(rarefy.ps.soil)
  t<-data.frame(Rarefy=i,OTU.No=rarefy.richness[,"Observed"])
  rownames(t)<-rownames(rarefy.richness)
  t<-cbind(sample.info.soil[rownames(t), c("Plant_sample","mt.gRNA", "Cas9")], t)
  rarefy.observed<-rbind(rarefy.observed, t)
  rm(t)
}

tail(rarefy.observed)
rarefy.observed$Cas9<-relevel(rarefy.observed$Cas9, ref="No")
p.rarefy<-ggplot(aes(x=Rarefy/1000, y=OTU.No, colour=Cas9),data=rarefy.observed)
p.rarefy<-p.rarefy+geom_smooth(fill="white",weight=0.5)+facet_wrap(~Plant_sample)
p.rarefy<-p.rarefy+labs(x="Sequences (x 1000)", y="OTU Number")
p.rarefy
ggsave(filename = "../SLY_results/Rare_Curve_soil_v4_OTU100_threshold10.pdf", p.rarefy,
       width = 15, height = 5, 
       units = "in")
rm(t,rarefy.ps.soil, p.rarefy, rarefy.observed, rarefy.otu,rarefy.richness)


####################################################
##  B2 Normalize the soil data by rarefaction  ###
##  Normalize to 16000 / sample, same to 799F-1193R data

# rarefaction sample_size 16000
min(sampleSums(ps.soil))
rarefy.ps.soil<-rarefy_even_depth(ps.soil,16000, replace=F, rngseed = 17)
#save the rarefied dataset
#saveRDS(rarefy.ps.soil, file="../SLY_results/rarefyl.ps.soil_v4_20200303.RData")  #rarefy 32000
#load the previously rarefied soil dataset
rarefy.ps.soil<-readRDS("../SLY_results/rarefyl.ps.soil_v4_20200303.RData") 
sample_sums(rarefy.ps.soil)

# just check the OTUs of rice mt-rDNA and cp-rDNA
tax_table(rarefy.ps.soil)[c("OTU_1","OTU_3"),] 

# Check all mitochrondria OTUs (18 OTUs)
prune_taxa(data.frame(tax_table(rarefy.ps.soil))$Family == "Mitochondria",rarefy.ps.soil)
# check the rice mitochrondria OTUs (13 OTUs)
prune_taxa(rownames(otu_table(rarefy.ps.soil)) %in% OTU.mt$OTU_ID, rarefy.ps.soil)

## NOTE: in following analysis, the rice OTUs are removed from the rarefied soil dataset

## remove the rice mt-rRNA OTU, resultes 17866 OTUs
rarefy.ps.soil<-prune_taxa(!rownames(otu_table(rarefy.ps.soil)) %in% OTU.mt$OTU_ID, rarefy.ps.soil)
#remove rice cp-rRNA OTUs
rarefy.ps.soil<-prune_taxa(!rownames(otu_table(rarefy.ps.soil)) %in% OTU.cp$OTU_ID, rarefy.ps.soil)# remove other mitochondria OTU, resulting 14089 OTUs
               
#rarefy.ps.soil<-prune_taxa(data.frame(tax_table(rarefy.ps.soil))$Family != "Mitochondria",rarefy.ps.soil)
#prune_taxa(data.frame(tax_table(rarefy.ps.soil))$Order == "Chloroplast",rarefy.ps.soil) # 43


# If want to remove all chloroplast OTUs
# rarefy.ps.soil<-prune_taxa(data.frame(tax_table(rarefy.ps.soil))$Phylum != "Chloroplast",rarefy.ps.soil)

## go through the data
sample_sums(rarefy.ps.soil)

table(rowMedians(otu_table(rarefy.ps.soil)))
table(otu_table(rarefy.ps.soil))
apply(otu_table(rarefy.ps.soil),2,max) 
##  output of number of most abundant OTU in each sample
# V4BS4 V4BS5 V4BS6 V4BS1 V4BS2 V4BS3 
# 469   410   525   436   465   569 

###########################################################
##  B3. Compare the OTUs that are presented in Cas9+/- samples 
##  Examine the coverage of Cas9+ and Cas9- soil samples
##  

## split sample info according gRNA
lst.sample<-split(sample.info.soil, as.character(sample.info.soil$mt.gRNA))
gRNA_id<-names(lst.sample)

## OTUs less the threshold labeled as no detected.
count.overlap <- function(t, threshold=0) { #make sure column 1 is Cas9 plus
  a<-length(t[t[,1]>threshold,1])
  b<-length(t[t[,2]>threshold,1])
  a_b<-length(t[t[,2]>threshold&t[,1]>threshold,1]);
  total<-a+b-a_b
  count.ov<-c(total, a-a_b, b-a_b, a_b);
  names(count.ov)<-c("Total", colnames(t),"OV")
  count.ov
}

overlap.info.all<-data.frame(
  gRNA  = as.character(),
  Sample= as.character(), 
  Cat   = as.character(),
  Num   = as.numeric(),
  RA_thresh = as.numeric()
  #Ov    = as.numeric()
)
p.ov<-ggplot()
p.all<-ggplot()
## set up the rathreshold range
## a.k.a, thred of OTU count = c(6.4, 3.2, 2, 0)

thresh=c(0.0005, 0.0002, 0.0001, 0)
for(n in 1:length(thresh)){
  ra.threshold=thresh[n];
  
  for(j in 1:length(gRNA_id)) {
    subset.sample<-lst.sample[[gRNA_id[j]]]
    subset.otutab<-otu_table(subset_samples(rarefy.ps.soil, mt.gRNA==gRNA_id[j]))
    subset.otutab<-subset.otutab[apply(subset.otutab, 1, sum) > 0, ]
    
    # convert to relative abundance (RA, %)
    subset.otutab<-transform_sample_counts(subset.otutab, function(x) {x/sum(x)})
    sample_sums(subset.otutab)
    
    #split according to samples
    subset.sample$Plant_sample<-as.character(subset.sample$Plant_sample)
    split_sample_info<-split(subset.sample, subset.sample$Plant_sample)
    overlap.info<-data.frame(
      gRNA  = as.character(),
      Sample= as.character(), 
      Cat   = as.character(),
      Num   = as.numeric(),
      RA_thresh = as.numeric()
      #Ov    = as.numeric()
    )
    for (i in 1:length(split_sample_info)) {
      t<-split_sample_info[[i]]
      t.otu<-subset.otutab[,colnames(subset.otutab) %in% t$Seq_ID]
      
      ## for each soil sample, keep OTUs whose mean(RA)> ra.threshold
      abund.filter<-apply(t.otu,1, mean)>ra.threshold
      t.otu<-t.otu[abund.filter,]
      #t.otu<-data.frame(t.otu)
      colnames(t.otu)<-t[colnames(t.otu),"Cas9"]  # sample_info rownames = t.otu colnames
      overlap<-count.overlap(t.otu, threshold = 0)
      overlap.df<-data.frame(gRNA=gRNA_id[j],Sample= names(split_sample_info)[i], Cat=names(overlap)[c(2,4,3)], Num=overlap[c(2,4,3)])
      overlap.df$RA_thresh<-ra.threshold
      #Replace Yea and No with Cas9+ and Cas9-
      overlap.df$Cat<-gsub("Yes", "Cas9+", overlap.df$Cat)
      overlap.df$Cat<-gsub("No", "Cas9-", overlap.df$Cat)
      overlap.info<-rbind(overlap.info, overlap.df,make.row.names=F)
      overlap.info.all<-rbind(overlap.info.all, overlap.df,make.row.names=F)
    }
    overlap.info$Cat<-factor(overlap.info$Cat, levels=c("Cas9-","OV","Cas9+"))
    x<-(n-1)*length(gRNA_id)+j
    p.ov[[x]]<-ggplot(data=overlap.info, aes(x=Sample, y=Num, fill=Cat))+geom_bar(stat="identity", width = 0.6)
    p.ov[[x]]<-p.ov[[x]]+labs(title =(paste(gRNA_id[j],"cp-gRNA578")), x=paste(split_sample_info[[1]][1,1],"samples"), y= paste("OTU(RA):",ra.threshold))
    p.ov[[x]]<-p.ov[[x]]+geom_text(aes(label=Num),vjust=1.0, color="black")
    p.ov[[x]]<-p.ov[[x]]+theme(text = element_text(size=6), plot.title = element_text(hjust = 0.5))
    p.ov[[x]]<-p.ov[[x]]+theme_gray()
  }  
}
p.ov.all<-plot_grid(p.ov[[1]],p.ov[[2]],p.ov[[3]], p.ov[[4]],
                    nrow = 2, align = 'v')
p.ov.all
#plot with bar chart
ggsave(filename = "../SLY_results/Figure. Venn of observed v4_OTU at RA_threshold_0202.pdf", p.ov.all,
       width = 15, height = 15, 
       units = "in")
write.table(overlap.info.all, file="../SLY_results/Table overlap v4 OTU at RA_threshold_0202.csv")
effective_OTU_No<-overlap.info.all[overlap.info.all$RA_thresh==2e-4,]

# be careful to sample order
sample.info.soil.rarefy$Effective_OTU<-c(with(effective_OTU_No, Num[Cat=="OV"]+Num[Cat=="Cas9-"]),
    with(effective_OTU_No, Num[Cat=="OV"]+Num[Cat=="Cas9+"]))

rm(overlap.info.all, t, abund.filter, t.otu, p.ov,p.ov.all, x, thresh,
   overlap, overlap.df, overlap.info,  count.overlap, effective_OTU_No)



#########################################################################
## B4 compare the microbial diversities between Cas9+ and Cas9- samples

## if want to remove low abundant OTU with average frequency <=1
## rarefy.ps.soil<-prune_taxa(rowMeans(otu_table(rarefy.ps.soil))>2,rarefy.ps.soil) 

# all OTUs are used to estimate the bacterial diversities
sample_sums(rarefy.ps.soil)
summary(rowSums(otu_table(rarefy.ps.soil)))
soil_richness<-estimate_richness(rarefy.ps.soil)
sample.info.soil.rarefy<-cbind(sample.info.soil, soil_richness[rownames(sample.info.soil),])

p.richness<-qplot(Cas9, Observed, fill = Cas9, data=sample.info.soil.rarefy)
p.richness<-p.richness+geom_boxplot()
p.richness

pslog.soil<-transform_sample_counts(prune_taxa(rowMeans2(otu_table(rarefy.ps.soil))>2, rarefy.ps.soil), function(x) log(1+x))
out.log.soil<-ordinate(pslog.soil, method="PCoA", distance="bray")
p.ord.soil<-plot_ordination(pslog.soil, out.log.soil, color="Cas9", shape="gRNA")+geom_point(size=I(2))
p.ord.soil

## summarize soil samples at genus/family levels
rarefy.ps.soil
rarefy.ps.soil.genus<-tax_glom(rarefy.ps.soil, "Genus", NArm=TRUE)
rarefy.ps.soil.family<-tax_glom(rarefy.ps.soil, "Family", NArm=TRUE)
sample_sums(rarefy.ps.soil.genus)
sample_sums(rarefy.ps.soil.family)

rarefy.ps.soil.family<-prune_taxa(taxa_sums(rarefy.ps.soil.family)>6,rarefy.ps.soil.family)
## Show the heatmap of abundance at family level for rarefied soil samples
heatmap(otu_table(rarefy.ps.soil.family))
heatmap(otu_table(rarefy.ps.soil.family), Rowv=NA)
p.heat<-plot_heatmap(rarefy.ps.soil.family)
#p.heat<-plot_heatmap(rarefy.ps.soil.family, 
#                     sample.order="Plant_sample",
#                     taxa.label="Family (Whole Rare Dataset)", taxa.order="Phylum",
#                     low="#000033", high="#CCFF66",
#                     na.value = "white")
p.heat
ggsave(filename = "../SLY_results/Fig.v4_soil-rarefy-heatmap_family.pdf", p.heat,
       width = 5, height = 10, 
       units = "in")
rm(p.heat)

#######################################################
## B5. statistical testing of soil data
sample.info.soil.rarefy<-sample.info.soil.rarefy[order(sample.info.soil.rarefy[,"Cas9"]),]

## statistical testing differencec between Cas+ and Cas-
## observed OTU
t.test(sample.info.soil.rarefy$Observed[sample.info.soil.rarefy$Cas9=="No"],
       sample.info.soil.rarefy$Observed[sample.info.soil.rarefy$Cas9=="Yes"],
       paired = T, alternative = "two.sided", conf.level = 0.95)

## effective OTU
t.test(sample.info.soil.rarefy$Effective_OTU[sample.info.soil.rarefy$Cas9=="No"],
       sample.info.soil.rarefy$Effective_OTU[sample.info.soil.rarefy$Cas9=="Yes"],
       paired = T, alternative = "two.sided", conf.level = 0.95)

## test Shannon diversities
t.test(sample.info.soil.rarefy$Shannon[sample.info.soil.rarefy$Cas9=="No"],
       sample.info.soil.rarefy$Shannon[sample.info.soil.rarefy$Cas9=="Yes"],
       paired = T, alternative = "two.sided", conf.level = 0.95)
## test Simpson diversities
t.test(sample.info.soil.rarefy$Simpson[sample.info.soil.rarefy$Cas9=="No"],
       sample.info.soil.rarefy$Simpson[sample.info.soil.rarefy$Cas9=="Yes"],
       paired = T, alternative = "two.sided", conf.level = 0.95)

## test OTU composition of Cas9+ and Cas9- soil data
otu.matrix<-data.frame(otu_table(rarefy.ps.soil))
head(otu.matrix)
otu.df<-t(otu.matrix)
colnames(otu.df)<-rownames(otu.matrix)

## make sure otu.df data have same order to sample info dateframe
otu.df<-otu.df[as.character(sample.info.soil.rarefy$Seq_ID),]
adonis2(otu.df ~ Plant_sample+Cas9, data=sample.info.soil.rarefy, method="bray",permutations = 999)
rm(otu.df, otu.matrix)

## save soil sample information
write.table(sample.info.soil.rarefy, file="../SLY_results/Table. summary of v4 rarefied soil data.csv", sep = ",")

###############################################################
##  B6 scatter plot of OTU abundance in Cas9+/- data

# here, the mean threshold was not used by set to 0
mean.threshold=0  

#split sample info according gRNA #####
lst.sample<-split(sample.info.soil, as.character(sample.info.soil$cp.gRNA))
gRNA_id<-names(lst.sample)

# calculate the corealtion coefficient value
#otu.matrix<-as.matrix(otutab[,-1])  #Mt OTUs are included here
otu.matrix<-as.matrix(otu_table(rarefy.ps.soil))
dim(otu.matrix)  #check the data
cor.df<-cor(otu.matrix[,colnames(otu.matrix) %in% sample.info.soil$Seq_ID])
write.table(cor.df, file="../SLY_results/corr-R.txt")  
rm(cor.df)

p<-ggplot()
for(j in 1:length(gRNA_id)) { 
  #subset.sample<-lst.sample[[gRNA_id[j]]]
  #subset.seq<-colnames(otu.matrix) %in% subset.sample$Seq_ID  #Get sequence_ID
  subset.otutab<-otu.matrix
  
  ## filter low adundant OTUs ##
  subset.otutab<-subset.otutab[apply(subset.otutab, 1, sum) > 0, ]
  filter.abund<-apply(subset.otutab, 1,mean)>mean.threshold
  filter.freq<-apply(subset.otutab>0, 1, min)>0  # at least 1 read in 6 samples
  subset.otutab.filtered<-subset.otutab[filter.abund & filter.freq, ]
  dim(subset.otutab.filtered)
  dim(subset.otutab)
  cor(subset.otutab.filtered)
  cor(subset.otutab)
  head(subset.otutab)
  ## compare OTU abundances between Cas9 +/- samples ###
  
  subset.sample$Plant_sample<-as.character(subset.sample$Plant_sample)
  
  # get each pair of Cas9+/- samples
  split_sample_info<-split(subset.sample, subset.sample$Plant_sample)
  t.otu.all<-data.frame(plus=as.numeric(), minus=as.numeric(),rep=as.character())
  for (i in 1:length(split_sample_info)) {
    t<-split_sample_info[[i]]
    t.otu<-subset.otutab.filtered[,colnames(subset.otutab.filtered) %in% t$Seq_ID]
    t.otu<-data.frame(t.otu)
    colnames(t.otu)<-t[colnames(t.otu),"Cas9"]
    t.otu$rep<-names(split_sample_info)[i]
    t.otu.all<-rbind(t.otu.all, t.otu)
    #ggplot(x=log(plus), y=log(minus), data=data.frame(t.otu), geom=c("point"), size=I(0.5),asp = 1)
  }
  p<-qplot(x=log2(Yes), y=log2(No), data=t.otu.all, 
                geom=c("point"), 
                size=I(2),
                alpha=0.5,
                asp = 1, 
                shape= factor(rep), 
                col=factor(rep)) + geom_smooth(method="lm", se=F)
  p<-p+labs(title =gRNA_id[j], x="log2(Abund) Cas9+", y= "log2(Abund) Cas9-")
  p<-p+theme(text = element_text(size=16), plot.title = element_text(hjust = 0.5))
  
  p<-p+theme_gray()
  
}

#plot with figrure 4d
p
ggsave(filename = "../SLY_results/Fig. v4_rarefy_soil_otu-abund.pdf", p,
       width = 5, height = 5, 
       units = "in")
rm(p, t.otu.all, t.otu, split_sample_info, subset.otutab, subset.otutab.filtered,
   filter.abund,filter.freq)



#############################################################
## B7 Identify Differential OTUs between Cas+ and Cas- data

library(DESeq2)
otu_table(ps.soil)[c("OTU_1", "OTU_3")]  #examine two rice rDNA OTU

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

## Analyzing the differential OTUs in soil samples
## each gRNA have three biological repeats.
sig.alpha=0.01  # set significant p value
ps.soil
prune_taxa(rowSums(otu_table(ps.soil))>0,ps.soil)
sample_data(ps.soil)

Dif.OTU<-data.frame();
dds<-prune_taxa(rowSums(otu_table(ps.soil))>0,ps.soil)
keep <- rowSums(otu_table(dds)) >= 0  #set the threshold to filter OTUs
table(keep)
dds <- prune_taxa(keep,dds)
sample_data(dds)$Cas9<-relevel(sample_data(dds)$Cas9,ref="No")
ps.dds<-phyloseq_to_deseq2(dds, ~Cas9)
geoMeans<- apply(counts(ps.dds),1, gm_mean)
dds<- estimateSizeFactors(ps.dds, geoMeans = geoMeans)
dds<- DESeq(dds, test="Wald",fitType="parametric")
res<-results(dds)
res[c("OTU_1", "OTU_3"),]
summary(res)  # zero 
plotMA(res)
## Get 0 significant changed OTUs here.
sigtab = res[(res$padj < sig.alpha & !is.na(res$padj)), ]
rm(dds, res, sigtab,keep)


######################################################################
##   C. Processing the Leaf samples
##     1. Summarize the data
##     2. Rarefaction curve 
##     3. Summarize the normalized leaf data (Rarefy to 50k sequences)
##     4. Statistical test the difference of microbial diversities
##     5. Compare bacterial diversities at different taxa level
##     6. Analyze differential OTU in root samples (DEseq2)
######################################################################

## Note. After normalization, Cas9- data obtain 135,282,159 bacterial 
## sequences, while Cas9+ data obtain >30k bacterial sequences.

############################
## C1. summarize the data
 
## Check the OTU Number and mt-rDNA ratio with raw data
sample.info.leaf # this data is not rarefied
ps.leaf

p.OTU.No<-qplot(interaction(Cas9), OTU_No, fill = Cas9, data=sample.info.leaf)
p.OTU.No<-p.OTU.No+geom_boxplot()
p.OTU.No

mt.cp.bact.ratio<-with(sample.info.leaf, data.frame(Plant_sample = rep(Plant_sample,3),
                                  Seq_ID = rep(Seq_ID,3), 
                                  Cas9 = rep(Cas9,3),
                                  Cat = rep(c("mito","chlo","bact"),1, each=6),
                                  Percentage = c(-mito.ratio, -chlo.ratio,bact.ratio)))
mt.cp.bact.ratio<-mt.cp.bact.ratio[order(mt.cp.bact.ratio$Seq_ID),]

p<-ggplot(data=mt.cp.bact.ratio,
        aes(x=Seq_ID, y=Percentage*100, colour=Cat, fill=Cat)) + 
        geom_bar(stat = "identity", width = 0.75)+
        theme(axis.text.x = element_text(angle = 90, hjust = 1))+
        labs(x="samples",y="Percentage(%)")+
        geom_text(aes(label=sprintf("%.1f %%",Percentage*100)),vjust=1.0, color="black")+
        guides(size=guide_legend())
p
ggsave(filename="../SLY_results/Fig 6A cp_mt_Bar_plot.pdf", p,
       width = 8, height = 15, 
       units = "cm")

rm(p,mt.cp.bact.ratio)

####################################################
## C2 Rarefaction curve 
##   OTU count vs sequencing depth: Rarefaction 
## 

# makesure the mito OTU is included in analyze.
sample_size<-min(sample_sums(ps.leaf))  
rarefy.observed<-data.frame()
rarefy.otu<-data.frame()
for( i in seq(1, sample_size, by = 500)) {
  # setting the thredshold here, keep OTUs whose abundance more than xx.
  rarefy.otu<-otu_table(rarefy_even_depth(prune_taxa(rowSums(otu_table(ps.leaf))>0,ps.leaf), sample.size=i, replace=T))  
  # change to ps.root, prune_taxa(rowSums(otu_table(ps.root))>10,ps.root), ps.root.noMt
  t<-data.frame(Rarefy=i,OTU.No=colSums(rarefy.otu>0))
  t<-cbind(sample.info.leaf[rownames(t), c("Plant_sample","mt.gRNA", "Cas9", "Cycle_No")], t)
  rarefy.observed<-rbind(rarefy.observed, t)
  rm(t)
}
tail(rarefy.observed)
rarefy.observed$mt.gRNA<-relevel(rarefy.observed$mt.gRNA, ref="No")
t<-rarefy.observed[rarefy.observed$Cycle_No=="8",]
p.rarefy<-ggplot(aes(x=Rarefy/1000, y=OTU.No, colour=mt.gRNA, Cycle_No),data=t)
p.rarefy<-p.rarefy+geom_smooth(fill="white",weight=0.5)+facet_wrap(~Plant_sample)
p.rarefy<-p.rarefy+labs(x="Sequences (x 1000)", y="OTU Number")
p.rarefy
ggsave(filename = "../SLY_results/Fig.Rare_Curve_leaf_v4_OTU100_threshold_0.pdf", p.rarefy,
       width = 15, height = 5, 
       units = "in")
rm(t, p.rarefy, rarefy.observed, rarefy.otu)

####################################################
## C3 Summarize the normalized leaf data
##    50k sequences

## Rarefication of leaf samples
sampleSums(ps.leaf)  #
set.seed(157) 
#rarefy to min sampleSums, 50000
#rarefy.ps.leaf<-rarefy_even_depth(ps.leaf, sample.size=50000,replace=F,rngsee=157)
#saveRDS(rarefy.ps.leaf, file="../SLY_results/Rarefied_leaf_v4.RData")
rarefy.ps.leaf<-readRDS("../SLY_results/Rarefied_leaf_v4.RData")

colSums(otu_table(rarefy.ps.leaf)>0)

sample.info.leaf.rarefy<-sample.info.leaf[,c(1:7,11:15)]

leaf.mt<-prune_taxa(taxa_names(rarefy.ps.leaf) %in% OTU.mt$OTU_ID,rarefy.ps.leaf)
mito.no<-sample_sums(leaf.mt)
mito.ratio<-mito.no/sample_sums(rarefy.ps.leaf)
mito.ratio

leaf.cp<-prune_taxa(taxa_names(rarefy.ps.leaf) %in% OTU.cp$OTU_ID,rarefy.ps.leaf)
chlo.no<-sample_sums(leaf.cp)
chlo.ratio<-chlo.no/sample_sums(rarefy.ps.leaf)
chlo.ratio

bact.no<-sample_sums(rarefy.ps.leaf)-sample_sums(leaf.cp)-sample_sums(leaf.mt)
bact.ratio<-bact.no/sample_sums(rarefy.ps.leaf)

chlo_mito_bact_ratio<-rbind(mito.ratio,chlo.ratio,bact.ratio)
colSums(chlo_mito_bact_ratio)
chlo_mito_bact_ratio<-t(chlo_mito_bact_ratio)
as.character(sample.info.leaf.rarefy$Seq_ID)==as.character(rownames(chlo_mito_bact_ratio)) # check the order
sample.info.leaf.rarefy<-cbind(sample.info.leaf.rarefy, chlo_mito_bact_ratio[rownames(sample.info.leaf.rarefy),])
rm(chlo_mito_bact_no,chlo_mito_bact_ratio,bact.no,bact.ratio,
   chlo.no,chlo.ratio,mito.no,mito.ratio)

mt.cp.bact.ratio<-with(sample.info.leaf.rarefy, data.frame(Plant_sample = rep(Plant_sample,3),
                                                    Seq_ID = rep(Seq_ID,3), 
                                                    Cas9 = rep(Cas9,3),
                                                    Cat = rep(c("mito","chlo","bact"),1, each=6),
                                                    Percentage = c(-mito.ratio, -chlo.ratio,bact.ratio)))
mt.cp.bact.ratio<-mt.cp.bact.ratio[order(mt.cp.bact.ratio$Seq_ID),]

p<-ggplot(data=mt.cp.bact.ratio,
          aes(x=Seq_ID, y=Percentage*100, colour=Cat, fill=Cat)) + 
  geom_bar(stat = "identity", width = 0.75)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(x="samples(rarefaction 5E+4)",y="Percentage(%)")+
  geom_text(aes(label=sprintf("%.1f %%",abs(Percentage)*100)),vjust=1.0, color="black")+
  guides(size=guide_legend())
p
ggsave(filename="../SLY_results/Fig 6A cp_mt_Bar_plot_rarefy.pdf", p,
       width = 8, height = 15, 
       units = "cm")

rm(p,mt.cp.bact.ratio)


rarefy.ps.leaf
## remove the cp and mt to estimate bacterial diversities
## Total: 
## Os rDNA: cp: 4701, mt: 2963.
prune_taxa(taxa_names(rarefy.ps.leaf) %in% OTU.cp$OTU_ID,rarefy.ps.leaf)
prune_taxa(taxa_names(rarefy.ps.leaf) %in% OTU.mt$OTU_ID,rarefy.ps.leaf)
## Total: 16143  Chloroplast: 4707 taxa, mito: 3039
prune_taxa(data.frame(tax_table(rarefy.ps.leaf))$Order == "Chloroplast",rarefy.ps.leaf)
prune_taxa(data.frame(tax_table(rarefy.ps.leaf))$Family == "Mitochondria",rarefy.ps.leaf)

rarefy.ps.leaf.nohost<-prune_taxa(!taxa_names(rarefy.ps.leaf) %in% OTU.cp$OTU_ID,rarefy.ps.leaf)
rarefy.ps.leaf.nohost<-prune_taxa(!taxa_names(rarefy.ps.leaf.nohost) %in% OTU.mt$OTU_ID,rarefy.ps.leaf.nohost)
rarefy.ps.leaf.bacteria<-prune_taxa(data.frame(tax_table(rarefy.ps.leaf.nohost))$Family != "Mitochondria",rarefy.ps.leaf.nohost) 
#rarefy.ps.leaf.bacteria<-prune_taxa(data.frame(tax_table(rarefy.ps.leaf.bacteria))$Order != "Chloroplast",rarefy.ps.leaf.bacteria) 
rarefy.ps.leaf.bacteria  #7915 taxa

# filter low abundant OTU
abund.filter<-apply(otu_table(rarefy.ps.leaf.bacteria),1, median)
rarefy.ps.leaf.bacteria.filter<-prune_taxa(abund.filter>=1,rarefy.ps.leaf.bacteria)
# alternative filter, temove OTU with sum(abundance) <=20
#rarefy.ps.leaf.bacteria<-prune_taxa(taxa_sums(rarefy.ps.leaf.bacteria)>20,rarefy.ps.leaf.bacteria)

sample.info.leaf.rarefy<-cbind(sample.info.leaf.rarefy,
                               estimate_richness(rarefy.ps.leaf.bacteria)[as.character(sample.info.leaf.rarefy$Seq_ID),])
sample.info.leaf.rarefy$Effective_OTU<-estimate_richness(rarefy.ps.leaf.bacteria.filter)[as.character(sample.info.leaf.rarefy$Seq_ID),"Observed"]

p<-ggplot(aes(x=Seq_ID, y=Effective_OTU, fill=Cas9),
          data=sample.info.leaf.rarefy)+
          geom_bar(position="dodge", stat="identity", width = 0.75)+theme_grey()
p<-p+geom_text(aes(label=Effective_OTU),vjust=1.0, color="black")
p<-p+theme(axis.text.x = element_text(angle = 90, hjust = 0.5))
p               
ggsave(filename="../SLY_results/Fig 6B observed_bact_OTU_leaf_v4_rarefy_median filter2.pdf", p,
       width = 8, height = 8, 
       units = "cm")
rm(p)

############################################################################
##  C4 test the difference of microbial diversities          ###
##

sample.info.leaf.rarefy$Cas9<-relevel(sample.info.leaf.rarefy$Cas9,ref="No")

# bacterial sequence fraction
t.test(sample.info.leaf.rarefy$bact.ratio[sample.info.leaf.rarefy$Cas9=="No"],
       sample.info.leaf.rarefy$bact.ratio[sample.info.leaf.rarefy$Cas9=="Yes"],
       paired = T, alternative = "two.sided", conf.level = 0.95)


#Observed OTU number
t.test(sample.info.leaf.rarefy$Observed[sample.info.leaf.rarefy$Cas9=="No"],
       sample.info.leaf.rarefy$Observed[sample.info.leaf.rarefy$Cas9=="Yes"],
       paired = T, alternative = "two.sided", conf.level = 0.95)

#effective OTU number
t.test(sample.info.leaf.rarefy$Effective_OTU[sample.info.leaf.rarefy$Cas9=="No"],
       sample.info.leaf.rarefy$Effective_OTU[sample.info.leaf.rarefy$Cas9=="Yes"],
       paired = T, alternative = "two.sided", conf.level = 0.95)

#shannon index
t.test(log2(sample.info.leaf.rarefy$Shannon[sample.info.leaf.rarefy$Cas9=="No"]),
       log2(sample.info.leaf.rarefy$Shannon[sample.info.leaf.rarefy$Cas9=="Yes"]),
       paired = T, alternative = "two.sided", conf.level = 0.95)
#Simpson index
t.test(sample.info.leaf.rarefy$Simpson[sample.info.leaf.rarefy$Cas9=="No"],
       sample.info.leaf.rarefy$Simpson[sample.info.leaf.rarefy$Cas9=="Yes"],
       paired = T, alternative = "two.sided", conf.level = 0.95)

## test OTU composition of Cas9+ and Cas9- soil data

otu.matrix<-data.frame(otu_table(rarefy.ps.leaf.bacteria))
head(otu.matrix)
otu.df<-t(otu.matrix)
colnames(otu.df)<-rownames(otu.matrix)

## make sure otu.df data have same order to sample info dateframe
otu.df<-otu.df[as.character(sample.info.leaf.rarefy$Seq_ID),]
adonis2(otu.df ~ Plant_sample+Cas9, data=sample.info.leaf.rarefy, method="bray")
rm(otu.df, otu.matrix)

## save soil sample information
write.table(sample.info.leaf.rarefy, file="../SLY_results/Table v4_OTU summary of rarefied leaf data_median filter.csv", sep = ",")

###############################################################
## C5 compare the bacterial abundance at different taxa level
##  

## Show the heatmap of abundance at family level for rarefied soil samples

rarefy.ps.leaf.bacteria.family<-tax_glom(rarefy.ps.leaf.bacteria.filter, "Family", NArm=F)
rarefy.ps.leaf.bacteria.genus<-tax_glom(rarefy.ps.leaf.bacteria.filter, "Genus", NArm=F)

sample_sums(rarefy.ps.leaf.bacteria)
taxa_sums(rarefy.ps.leaf.bacteria.family)
rarefy.ps.leaf.bacteria.family
table(data.frame(tax_table(rarefy.ps.leaf.bacteria))[, "Family"])

## Draw the heatmap of bacterial families
p.heat<-plot_heatmap(rarefy.ps.leaf.bacteria.family, "NMDS","bray",
                     sample.order="Cas9",
                     #taxa.order = "Family",
                     taxa.label="Family", low="#000033", high="#CCFF66", na.value="white")

p.heat<-p.heat+scale_fill_viridis_c(trans="log2",option = "D",na.value = grey(0.2))
p.heat<-p.heat+labs(x=paste("sample","(median>=1)"))
p.heat
ggsave(filename = "../SLY_results/Fig.v4_leaf-rarefy-heatmap_family_4.pdf", p.heat,
       width = 5, height = 10, 
       units = "in")
rm(p.heat)

## Instead, make the heatplot using ggplot2
psmelt(rarefy.ps.leaf.bacteria.family)
p.heat2<-ggplot(psmelt(rarefy.ps.leaf.bacteria.family), aes(x=Sample, y=Family,fill = Abundance)) + 
        geom_tile(linejoin = "bevel")
#try different color 
# low-high, ""#000033"-"#66CCFF"
#p.heat2<-p.heat2+scale_fill_gradient(low = "#000033", high = "#66CCFF", trans="log2",na.value = "white")
p.heat2+scale_fill_viridis_c(trans="log2",option = "D")
rm(p.heat2)

## export the bacterial OTU/genus/family
## modify seq names in otu.seq (OTU_1;size=xxx to OTU_1)
otu_names<-data.frame(do.call(rbind, strsplit(names(otu.seq),";")))
names(otu.seq)<-otu_names[,1]

leaf.otu<-data.frame(OTU_ID=NA,otu_table(rarefy.ps.leaf.bacteria.filter))
leaf.otu$OTU_ID<-rownames(leaf.otu)
leaf.taxa<-tax_table(rarefy.ps.leaf.bacteria.filter)
leaf.otu<-cbind(leaf.otu, leaf.taxa[rownames(leaf.otu),])
dim(leaf.otu)
leaf.otu$Seq<-as.character(otu.seq[rownames(leaf.otu)])
write.table(leaf.otu,file="../SLY_results/Table Leaf_v4_OTU_abund_median_filtered.txt", sep="\t", row.names = F)
rm(leaf.otu, otu_names, leaf.taxa)
# Exprot the OTU clustered at Family level
leaf.family<-data.frame(otu_table(rarefy.ps.leaf.bacteria.family))
leaf.taxa<-tax_table(rarefy.ps.leaf.bacteria.family)[,1:5]
leaf.family<-cbind(leaf.family, leaf.taxa[rownames(leaf.family),])
leaf.family$Seq<-as.character(otu.seq[rownames(leaf.family)])
write.table(leaf.family,file="../SLY_results/Table Leaf_v4_Family_abund_median_filter.txt", sep="\t")
rm(leaf.family, leaf.taxa)

## Phylum level
pslog<-transform_sample_counts(rarefy.ps.leaf.bacteria, function(x) log(1+x))
out.log<-ordinate(pslog, method="NMDS", distance="bray")
p.ord<-plot_ordination(pslog, out.log, color="Cas9", shape="Plant_sample")+geom_point(size=I(2))
p.ord
ggsave(filename = "../SLY_results/v4_leaf_rarefied_plot_Ordination_NMD.pdf", p.ord, width=3, height=3, units= "in")
rm(pslog, p.ord, out.log)

rarefy.ps.leaf.bacteria.phyla<-tax_glom(rarefy.ps.leaf.bacteria,taxrank = "Phylum")
rarefy.ps.leaf.bacteria.phyla<-transform_sample_counts(rarefy.ps.leaf.bacteria.phyla, function(x) {x/sum(x)} )
rarefy.ps.leaf.bacteria.phyla<-psmelt(rarefy.ps.leaf.bacteria.phyla)

#ps.leaf.phyla<-ps.leaf.phyla[ps.leaf.phyla$Abundance>0.01,]
table(as.character(rarefy.ps.leaf.bacteria.phyla$Phylum))
summary(ps.leaf.phyla$Abundance)

phylum_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861"
)

#
p.phylum<-ggplot(data=rarefy.ps.leaf.bacteria.phyla,aes(x=Seq_ID, y=Abundance, colour=Phylum, fill=Phylum)) + 
  geom_bar(stat = "identity")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(x="Cas9",y="Abundance")
p.phylum
#ggsave(filename = "Relative_abundance_Phyl_otu100_10_rarefy.pdf", p.phylum,
#       width = 5, height = 10, 
#       units = "in")
rm(ps.root.phyla, p.phylum,phylum_colors)


##########################################################
##  C6 Analyzing differential OTU in root samples
##  DESseq2

## This code is not valid since the sequencing depgh varied >10 times
## between Cas9+ and Cas9-

library(DESeq2)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

## Remove Mito and chlo sequences
#ps.leaf.noMt<-prune_taxa(!taxa_names(rarefy.ps.leaf) %in% OTU.mt$OTU_ID, rarefy.ps.leaf)
#ps.leaf.noMt.noCp<-prune_taxa(!taxa_names(ps.leaf.noMt) %in% OTU.cp$OTU_ID, ps.leaf.noMt)
#ps.leaf.noMt.noCp<-prune_taxa(rowSums(otu_table(ps.leaf.noMt.noCp))>10,ps.leaf.noMt.noCp)
#sample_data(ps.leaf.noMt.noCp)
#sample_data(ps.leaf.noMt.noCp)$Cas9<-relevel(sample_data(ps.leaf.noMt.noCp)$Cas9,ref="No")
#dds<-phyloseq_to_deseq2(ps.leaf.noMt.noCp, ~Cas9)

# load ps.leaf for DESeq 
sample_data(ps.leaf)$Cas9<-relevel(sample_data(ps.leaf)$Cas9,ref="No")
dds<-phyloseq_to_deseq2(ps.leaf, ~Cas9)
geoMeans<- apply(counts(dds),1, gm_mean)
dds<- estimateSizeFactors(dds, geoMeans = geoMeans)
dds<- DESeq(dds, fitType="parametric", test="Wald")
res=results(dds)
plotMA(res)
sigtab = res[(res$padj < 0.01 & !is.na(res$padj)), ]  #15 OTUs
#sigtab = res[(res$pvalue < 0.01 & !is.na(res$pvalue)), ]  #29 OTUs

sigtab<-data.frame(sigtab)
dim(sigtab)

sigtab<-cbind(sigtab,tax_table(ps.leaf)[rownames(sigtab),])
sigtab$Seq<-as.character(otu.seq[rownames(sigtab),])
write.table(sigtab,file = "../SLY_results/Table v4_leaf_diff_OTUs.csv",sep=",")
# all reduced OTU in Cas9+ are chloroplast or mitochondrion rRNA 
sigtab[sigtab$log2FoldChange<0,1:11]

#checking  following OTUs
sigtab["OTU_6",]    #Stenotrophomonas, high abundant
sigtab["OTU_4",]    #Achromobacter

## Done!
