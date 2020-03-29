# Cas-16S-seq
# This R script is used to compare the 799F-1193R 
# data of soil and root data. 
# The OTU is clustered using VSEARCH pipeline at 
# 100% identities.
#  
# Kabin Xie
# 2019.9.5

# revised
# 2020.1.30

##################################################
## Analysis procedures                          ##
## 1. Load the otu table and sample information ##  
## 2. Build the phyloseq object                 ##
## 3. Process the soil samples                  ##
## 4. Process the root samples                  ##
## 5. Identify differential OTUs (DESeq2)       ##
##################################################

library(ggplot2)
library(cowplot)
library(gridExtra)
library(phyloseq)
library(vegan)
library(agricolae)

########################################################################################
# 1. Load the otu table and sample information                                         #
########################################################################################

# remove the the # symbol at the header line before load the otu table
# Pay attention to col names if it start with digit. An X is added if the colname/rowname 
# started with digit.

setwd("~/mystation/microbiome/Cas-16S-seq/Analyze/vsearch100")
otutab<-read.table(file="../vsearch.out/all.otutab100.txt", head=T, check.names=FALSE)
head(otutab)
rownames(otutab)<-otutab[,1]

# In our analysis pipeline, OTUs (clustered by 100%) which are aligned to rice mito 16S rRNA
# with no more than 2 mismatches are considered as rice rDNA (all.otu100.mt2.txt)
OTU.mt<-read.table(file="all.otu100.mt2.txt", head=T)
table(!otutab$OTU_ID %in% OTU.mt$OTU_ID)
otu.matrix<-as.matrix(otutab[!otutab$OTU_ID %in% OTU.mt$OTU_ID,-1])  #-1, remove first name column

head(otu.matrix)
dim(otu.matrix)
dim(otutab)

# get the percentage of rice mt-rDNA for each sample
mito.no<-colSums(otutab[otutab$OTU_ID %in% OTU.mt$OTU_ID, -1])
mito.ratio<-mito.no/apply(otutab[,-1],2, sum)
mito.ratio
mito_bac<-rbind(mito.no, apply(otu.matrix,2, sum))
rownames(mito_bac)<-c("Mito", "Bact")
barplot(as.matrix(mito_bac[,1:10]))
mito_bac<-rbind(mito.ratio, 1-mito.ratio)
row.names(mito_bac)<-c("Mito", "bact")

## load the sample meta data table
sample.info<-read.table("../vsearch_sample_info.csv", head=T, sep=",")
rownames(sample.info)<-sample.info$Seq_ID
mito_bac<-t(mito_bac)
as.character(sample.info$Seq_ID)==as.character(rownames(mito_bac)) # check the order
sample.info<-cbind(sample.info, mito_bac[rownames(sample.info),])

## read NGS-OTU processing data (vsearch-pipline summary) ##
vs.summary<-read.table(file="../vs.pipline.summary100.csv", head=TRUE, sep=",")
head(vs.summary)
rownames(vs.summary)<-vs.summary$sample
sample.info<-cbind(sample.info, vs.summary[rownames(sample.info), c(2,3,4,5)])
OTU_No<-data.frame(OTU_No=colSums(otutab[,-1]>0))
 
# combined sample meta data, NGS-OTU processing data, and mt-rDNA data
sample.info<-cbind(sample.info,OTU_No=OTU_No[rownames(sample.info),], Total=colSums(otutab[,-1])[rownames(sample.info)])
sample.info<-sample.info[with(sample.info, order(gRNA, Compartment,Cas9, Cycle_No)),]
write.table(sample.info, file="sample.info.combined100.txt", sep="\t")
rm(vs.summary, OTU_No)

#remove the negative control
# bS1_4 is one of the negative controls that were used  to monitor the contamination
sample.info=sample.info[sample.info$Seq_ID!="bS1_4",]  #remove the empty control
sample.info.soil<-sample.info[sample.info$Compartment == "Soil",]
sample.info.root<-sample.info[sample.info$Compartment == "Root",]
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
ggsave(filename = "derep_ratio.pdf", p.reads,
       width = 10, height = 10, 
       units = "cm")
rm(p.reads)

##########################################################
#  2. Build the phyloseq object                          #
#    assign taxanomy and assemble the phyloseq object    #
##########################################################
library(phyloseq)
library(gridExtra)
library(DECIPHER)
library(dada2)
library(Biostrings)
library(ShortRead)
#get OTU sequences and assign taxanomy
otu.seq<-readDNAStringSet("../vsearch.out/all.otus100.fasta", format="fasta")
table((width(otu.seq)))

#taxonomy assign taken more than 4 hours with 7 cores.
set.seed(7)
#taxa <- assignTaxonomy(otu.seq, "~/mystation/microbiome/silva_nr_v132_train_set.fa.gz", multithread=TRUE, minBoot=60)
#save the workspace here, RStudio may crash here. (set n to read n seq each time, save memory)
#saveRDS(taxa, file="taxa_assign.RData")
taxa<-readRDS("taxa_assign.RData")  # Load the taxa information
#taxa.specis<-addSpecies(taxa, "~/mystation/microbiome/silva_species_assignment_v132.fa.gz", n=1000, verbose=F)
#saveRDS(taxa.specis, file="species_assign.RData")

## rowname of taxa is the OTU sequence which is named by fasta header ##
taxa.rownames.matrix<-matrix(unlist(strsplit(names(rownames(taxa)), ";")), ncol=2, byrow=T)
rownames(taxa)<-taxa.rownames.matrix[,1]

otu.seq.bact<-otu.seq[names(otu.seq)[taxa[,5]!="Mitochondria"]]
writeXStringSet(otu.seq.bact, filepath = "otu100.bact.fasta")
rm(otu.seq.bact)

ps <- phyloseq(otu_table(otutab[,-1], taxa_are_rows=T), 
               sample_data(sample.info), 
               tax_table(taxa)) 

# subset root and soil samples #
ps.root<-phyloseq(otu_table(otutab[, as.character(sample.info.root$Seq_ID)], taxa_are_rows = T),
                  sample_data(sample.info.root),
                  tax_table(taxa))
ps.root<-prune_taxa(rowSums(otu_table(ps.root))>0, ps.root )
ps.soil<-phyloseq(otu_table(otutab[,as.character(sample.info.soil$Seq_ID)],taxa_are_rows = T), 
                  sample_data(sample.info.soil),
                  tax_table(taxa))
ps.soil<-prune_taxa(rowSums(otu_table(ps.soil))>0, ps.soil )

sample_sums(prune_taxa(rowSums(otu_table(ps.soil))>10,ps.soil))
ps_richness<-estimate_richness(ps, measures =c("Shannon","Simpson"))
sample.info<-cbind(sample.info, estimate_richness(ps, measures =c("Shannon","Simpson"))[as.character(sample.info$Seq_ID),])

#remove the rice Mt sequences
ps.noMt<-prune_taxa(!rownames(otu_table(ps)) %in% OTU.mt$OTU_ID, ps)
ps_richness<-estimate_richness(ps.noMt, measures =c("Shannon","Simpson"))
colnames(ps_richness)<-c("Shannon_NoMt","Simpson_NoMt")
sample.info<-cbind(sample.info, ps_richness[as.character(sample.info$Seq_ID),])
write.table(sample.info, file="sample.info.with.alpha.txt", sep="\t")

#########################################################
##  3. Processe soil samples                         ##
##     analyze OTUs in soil sample (Fig. 4)            ##
#########################################################
ps.soil

##### 3.1. Rarefication curve for soil samples #########
rarefy.observed<-data.frame()
rarefy.otu<-data.frame()
set.seed(99)
sample_size=40000
for( i in seq(1, sample_size, by = 500)) {
  rarefy.ps.soil<-rarefy_even_depth(ps.soil,i, replace=F)
  # abundance threshold
  if(i>10) {
    rarefy.ps.soil<-prune_taxa(taxa_sums(rarefy.ps.soil)>10,rarefy.ps.soil)
  }
  
  rarefy.richness<-estimate_richness(rarefy.ps.soil)
  t<-data.frame(Rarefy=i,OTU.No=rarefy.richness[,"Observed"])
  rownames(t)<-rownames(rarefy.richness)
  t<-cbind(sample.info.soil[rownames(t), c("Plant_sample","gRNA", "Cas9")], t)
  rarefy.observed<-rbind(rarefy.observed, t)
  rm(t)
}

tail(rarefy.observed)
rarefy.observed$gRNA<-relevel(rarefy.observed$gRNA, ref="No")
p.rarefy<-ggplot(aes(x=Rarefy/1000, y=OTU.No, colour=Cas9),data=rarefy.observed)
p.rarefy<-p.rarefy+geom_smooth(fill="white",weight=0.5)+facet_wrap(~gRNA+Plant_sample)
p.rarefy<-p.rarefy+labs(x="Sequences (x 1000)", y="OTU Number")
p.rarefy
ggsave(filename = "Rare_Curve_soil_OTU100_threshold10.pdf", p.rarefy,
       width = 15, height = 5, 
       units = "in")
rm(t,rarefy.ps.soil, p.rarefy, rarefy.observed, rarefy.otu)


####################################################
##  3.2 Normalize the soil data by rarefaction  ###
##  Normalize to 16000 / sample

# rarefaction sample_size 16000
rarefy.ps.soil<-rarefy_even_depth(ps.soil,16000, replace=F, rngseed = 7)
#save the rarefied dataset
#saveRDS(rarefy.ps.soil, file="rarefyl.ps.soil_20200128.RData")
#load the previously rarefied soil dataset
rarefy.ps.soil<-readRDS("rarefyl.ps.soil_20200128.RData") 
sample_sums(rarefy.ps.soil)

#remove the mitochondrion
rarefy.ps.soil<-prune_taxa(data.frame(tax_table(rarefy.ps.soil))$Family != "Mitochondria",rarefy.ps.soil)
rarefy.ps.soil<-prune_taxa(!rownames(otu_table(rarefy.ps.soil)) %in% OTU.mt$OTU_ID, rarefy.ps.soil )

sample_sums(rarefy.ps.soil)
table(rowMedians(otu_table(rarefy.ps.soil)))

###########################################################
##  3.3 
##  Examine the detected of OTU between Cas9+ and Cas9- soil samples

## split sample info according gRNA
lst.sample<-split(sample.info.soil, as.character(sample.info.soil$gRNA))
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
## a.k.a, thred of OTU count = c(8.0, 3.2, 1.6, 0)
thresh=c(0.0005, 0.0002, 0.0001, 0)
for(n in 1:length(thresh)){
  ra.threshold=thresh[n];
  
  for(j in 1:length(gRNA_id)) {
    subset.sample<-lst.sample[[gRNA_id[j]]]
    subset.otutab<-otu_table(subset_samples(rarefy.ps.soil, gRNA==gRNA_id[j]))
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
    p.ov[[x]]<-p.ov[[x]]+labs(title =gRNA_id[j], x=paste(split_sample_info[[1]][1,1],"samples"), y= paste("OTU(RA):",ra.threshold))
    p.ov[[x]]<-p.ov[[x]]+geom_text(aes(label=Num),vjust=1.0, color="black")
    p.ov[[x]]<-p.ov[[x]]+theme(text = element_text(size=6), plot.title = element_text(hjust = 0.5))
    p.ov[[x]]<-p.ov[[x]]+theme_gray()
  }  
}
p.ov.all<-plot_grid(p.ov[[1]],p.ov[[2]],p.ov[[3]], 
                    p.ov[[4]],p.ov[[5]],p.ov[[6]],
                    p.ov[[7]],p.ov[[8]],p.ov[[9]], 
                    p.ov[[10]],p.ov[[11]],p.ov[[12]],
                    ncol = 3, align = 'v')

#plot with bar chart
p.ov.all
ggsave(filename = "Figure C2.overlapOTU vs RA_threshold_0128.pdf", p.ov.all,
       width = 15, height = 15, 
       units = "in")
write.table(overlap.info.all, file="Table C2 overlapOTU vs RA_threshold_0128.csv")

effective_OTU_No<-overlap.info.all[overlap.info.all$RA_thresh==2e-4,]
rm(overlap.info.all, t, abund.filter, t.otu, p.ov,p.ov.all, x, thresh,
   overlap, overlap.df, overlap.info,  count.overlap)



#########################################################################
## 3.4 compare the microbial diversities between Cas9+ and Cas9- samples

## Above plot indicates that RA > 1e-4 per sample is suitable to filter low-abund OTU
## it is close to mean(OTU) across all samples > 1
## We do not filter low-abundant OTU when calculate bacterial diversities
## In case to remove low abundant OTU, try following code 
rarefy.ps.soil.filter<-prune_taxa(taxa_sums(rarefy.ps.soil)>18,rarefy.ps.soil) 
#rarefy.ps.soil<-prune_taxa(rowMeans2(otu_table(rarefy.ps.soil))>1,rarefy.ps.soil) 

## estimate bacterial diversities
sample_sums(rarefy.ps.soil)
summary(rowSums(otu_table(rarefy.ps.soil)))
soil_richness<-estimate_richness(rarefy.ps.soil)
sample.info.soil.rarefy<-cbind(sample.info.soil[,1:14], soil_richness[rownames(sample.info.soil),])

sample.info.soil.rarefy$Effective_OTU<-NA
sample.info.soil.rarefy$Effective_OTU[sample.info.soil.rarefy$Cas9=="No"]<-with(effective_OTU_No, Num[Cat=="OV"]+Num[Cat=="Cas9-"])
                                         
sample.info.soil.rarefy$Effective_OTU[sample.info.soil.rarefy$Cas9=="Yes"]<-with(effective_OTU_No, Num[Cat=="OV"]+Num[Cat=="Cas9+"])
rm(effective_OTU_No)

p.richness<-qplot(interaction(Cas9, gRNA), Observed, fill = Cas9, data=sample.info.soil.rarefy)
p.richness<-p.richness+geom_boxplot()
p.richness

pslog.soil<-transform_sample_counts(prune_taxa(rowMeans2(otu_table(rarefy.ps.soil))>2, rarefy.ps.soil), function(x) log(1+x))
out.log.soil<-ordinate(pslog.soil, method="PCoA", distance="bray")
p.ord.soil<-plot_ordination(pslog.soil, out.log.soil, color="Cas9", shape="gRNA")+geom_point(size=I(2))
p.ord.soil

## Coverage of soil samples at genus/family levels
rarefy.ps.soil
rarefy.ps.soil.genus<-tax_glom(rarefy.ps.soil, "Genus", NArm=TRUE)
rarefy.ps.soil.family<-tax_glom(rarefy.ps.soil, "Family", NArm=TRUE)
sample_sums(rarefy.ps.soil.genus)
sample_sums(rarefy.ps.soil.family)
rarefy.ps.soil.family

## Show the heatmap of abundance at family level for rarefied soil samples
heatmap(otu_table(rarefy.ps.soil.family), Rowv=NA)
p.heat<-plot_heatmap(rarefy.ps.soil.family, "NMDS","bray",
             sample.order="Plant_sample",
             taxa.label="Phylum", taxa.order="Phylum", low="#000033", high="#CCFF66")
p.heat
ggsave(filename = "Fig 4C.soil-rarefy-heatmap0128_phylumname.pdf", p.heat,
       width = 5, height = 10, 
       units = "in")
rm(p.heat)

## statistical testing of soil data
sample.info.soil.rarefy<-sample.info.soil.rarefy[order(sample.info.soil.rarefy[,"Cas9"]),]

## test number of detected OTU
t.test(sample.info.soil.rarefy$Observed[sample.info.soil.rarefy$Cas9=="No"],
       sample.info.soil.rarefy$Observed[sample.info.soil.rarefy$Cas9=="Yes"],
       paired = T, alternative = "two.sided", conf.level = 0.95)

## test hypothesis greater
t.test(sample.info.soil.rarefy$Observed[sample.info.soil.rarefy$Cas9=="No"],
       sample.info.soil.rarefy$Observed[sample.info.soil.rarefy$Cas9=="Yes"],
       paired = T, alternative = "greater", conf.level = 0.95)

## Try wilcox.test, same conclusion.
## Cas9 treatment detected more OTUs (p = 0.037)
wilcox.test(sample.info.soil.rarefy$Observed[sample.info.soil.rarefy$Cas9=="No"],
            sample.info.soil.rarefy$Observed[sample.info.soil.rarefy$Cas9=="Yes"],
            paired = T, alternative = "greater", conf.level = 0.95)

## Effective OTUs
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

## permutational ANOVA to test OTU diversities 
otu.matrix<-data.frame(otu_table(rarefy.ps.soil))
head(otu.matrix)
otu.df<-t(otu.matrix)
colnames(otu.df)<-rownames(otu.matrix)

## make sure otu.df data have same order to sample info dateframe
otu.df<-otu.df[as.character(sample.info.soil.rarefy$Seq_ID),]
adonis2(otu.df ~ Plant_sample+Cas9, data=sample.info.soil.rarefy, method="bray")
rm(otu.df, otu.matrix)

## save soil sample information
write.table(sample.info.soil.rarefy[with(sample.info.soil.rarefy, order(gRNA, Cas9, Plant_sample)),], 
            file="Table S3. summary of rarefied soil data.csv", sep = ",")

###############################################################
##  3.6 scatter plot of OTU abundance in Cas9+/- data

# here, the mean threshold was not used by set to 0
mean.threshold=0  

#split sample info according gRNA #####
lst.sample<-split(sample.info.soil, as.character(sample.info.soil$gRNA))
gRNA_id<-names(lst.sample)

# calculate the corealtion coefficient value
#otu.matrix<-as.matrix(otutab[,-1])  #Mt OTUs are included here
otu.matrix<-as.matrix(otu_table(rarefy.ps.soil))
dim(otu.matrix)  #check the data
cor.df<-cor(otu.matrix[,colnames(otu.matrix) %in% sample.info.soil$Seq_ID])
write.table(cor.df, file="corr-R.txt")  
rm(cor.df)

p<-ggplot()
for(j in 1:length(gRNA_id)) { 
  subset.sample<-lst.sample[[gRNA_id[j]]]
  subset.seq<-colnames(otu.matrix) %in% subset.sample$Seq_ID  #Get sequence_ID
  subset.otutab<-subset(otu.matrix,select=subset.seq)
  
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
    colnames(t.otu)<-t[colnames(t.otu),5]
    t.otu$rep<-names(split_sample_info)[i]
    t.otu.all<-rbind(t.otu.all, t.otu)
    #ggplot(x=log(plus), y=log(minus), data=data.frame(t.otu), geom=c("point"), size=I(0.5),asp = 1)
  }
  p[[j]]<-qplot(x=log2(Yes), y=log2(No), data=t.otu.all, 
                geom=c("point"), 
                size=I(2),
                alpha=0.5,
                asp = 1, 
                shape= factor(rep), 
                col=factor(rep)) + geom_smooth(method="lm", se=F)
  p[[j]]<-p[[j]]+labs(title =gRNA_id[j], x="log2(Abund) Cas9+", y= "log2(Abund) Cas9-")
  p[[j]]<-p[[j]]+theme(text = element_text(size=16), plot.title = element_text(hjust = 0.5))
  
  p[[j]]<-p[[j]]+theme_gray()
  
}

#plot figrure 4B 
p.all<-plot_grid(p[[1]],p[[2]],p[[3]], nrow=1, labels=c("D", "E","F"),align = 'v', label_size = 20)
p.all
ggsave(filename = "Fig 4B. rarefy_soil_otu-paired.pdf", p.all,
       width = 15, height = 5, 
       units = "in")
rm(p.all, p)


######################################################################
##     4. Proccess the Root samples
##  
######################################################################
## 4.1 summarize the data
## Figure 5A-raw mito ratio
## Figure 5B-raw observed OTUs

## Check the OTU Number and mt-rDNA ratio with raw data
sample.info.root<-sample.info[sample.info$Compartment == "Root", ]
p.OTU.No<-qplot(interaction(Cas9, gRNA,Cycle_No), OTU_No, fill = Cas9, data=sample.info.root)
p.OTU.No<-p.OTU.No+geom_boxplot()
p.OTU.No

p.ratio<-ggplot(aes(y=Mito*100, x=interaction(gRNA, Cycle_No), fill=Cas9),ylim(0,100), data=sample.info.root)
p.ratio<-p.ratio+geom_boxplot()
p.ratio<-p.ratio+theme_grey()
p.ratio
ggsave(filename="Figure5A.pdf", plot_grid(p.ratio, p.OTU.No, rel_widths = c(1, 1)),
       width = 18, height = 6, 
       units = "cm")

rm(p.ratio, p.reads, p.OTU.No)


####################################################
##### 4.2 OTU count vs sequencing depth: Rarefaction 
## Fig 5b Rarefy curve

# makesure the mito OTU is included in analyze.
 sample_size<-min(sample_sums(ps.root))  
 rarefy.observed<-data.frame()
 
 for( i in seq(1, sample_size, by = 500)) {
  rarefy.otu<-otu_table(rarefy_even_depth(prune_taxa(rowSums(otu_table(ps.root))>0,ps.root), sample.size=i, replace=T))  
  # change to ps.root, prune_taxa(rowSums(otu_table(ps.root))>10,ps.root), ps.root.noMt
  t<-data.frame(Rarefy=i,OTU.No=colSums(rarefy.otu>0))
  t<-cbind(sample.info.root[rownames(t), c("Plant_sample","gRNA", "Cas9", "Cycle_No")], t)
  rarefy.observed<-rbind(rarefy.observed, t)
  rm(t)
 }
 tail(rarefy.observed)
 rarefy.observed$gRNA<-relevel(rarefy.observed$gRNA, ref="No")
 t<-rarefy.observed[rarefy.observed$Cycle_No=="8",]
 p.rarefy<-ggplot(aes(x=Rarefy/1000, y=OTU.No, colour=gRNA, Cycle_No),data=t)
 p.rarefy<-p.rarefy+geom_smooth(fill="white",weight=0.5)+facet_wrap(~Plant_sample)
 p.rarefy<-p.rarefy+labs(x="Sequences (x 1000)", y="OTU Number")
 p.rarefy
 ggsave(filename = "Figure5C Rare_Curve_8cycles_OTU100_threshold_0.pdf", p.rarefy,
        width = 15, height = 5, 
        units = "in")
rm(t, p.rarefy, rarefy.observed, rarefy.otu)


############################################################
## estimate bacterial diversities at same seq depth

## Rarefication of root samples
set.seed(157)
ps.root.rarefy<-rarefy_even_depth(ps.root, sample.size=41500,replace=F)
#saveRDS(ps.root.rarefy, file="rarefyl.ps.root_20200228.RData")
ps.root.rarefy<-readRDS(file = "rarefyl.ps.root_20200228.RData")
#Observed.OTU<-colSums(otu_table(ps.root.rarefy)>0)
prune_taxa(taxa_sums(ps.root.rarefy)>10,ps.root.rarefy) 
# >10, get 4150 taxa
# >24, get 2020 taxa
## Remove rice Mito sequences
## ps.root.noMt<-prune_taxa(!taxa_names(ps.root) %in% OTU.mt$OTU_ID, ps.root) # raw data
ps.root.noMt<-prune_taxa(!taxa_names(ps.root.rarefy) %in% OTU.mt$OTU_ID, ps.root.rarefy)

ps.root.noMt.filter<-prune_taxa(rowSums(otu_table(ps.root.noMt))>24,ps.root.noMt)
sample_sums(ps.root.noMt)
min(sample_sums(ps.root.noMt))

Observed.OTU<-colSums(otu_table(ps.root.noMt)>0)
sample.info.root.rarefy<-sample.info.root[,1:14]
sample.info.root.rarefy$Observed<-Observed.OTU[as.character(sample.info.root.rarefy$Seq_ID)]
sample.info.root.rarefy<-cbind(sample.info.root.rarefy, 
                               estimate_richness(ps.root.noMt, measures =c("Shannon","Simpson"))[as.character(sample.info.root.rarefy$Seq_ID),])
effective_OTU_No<-colSums(otu_table(ps.root.noMt.filter)>0)
sample.info.root.rarefy$Effective_OTU<-effective_OTU_No[rownames(sample.info.root.rarefy)]

## Get the mito ration in rarefied dataset
## Os mt-rDNA: 1411 OTUs. couple of OTUs match rice rDNA with 97% identities not include.
prune_taxa(rownames(otu_table(ps.root.rarefy)) %in% OTU.mt$OTU_ID, ps.root.rarefy )
## Mitochondriia: 1700 OTUs
ps.root.rarefy.mito<-prune_taxa(data.frame(tax_table(ps.root.rarefy))$Family == "Mitochondria",ps.root.rarefy)
mito.ratio<-sample_sums(ps.root.rarefy.mito)/sample_sums(ps.root.rarefy)
sample.info.root.rarefy$Mito<-mito.ratio[as.character(sample.info.root.rarefy$Seq_ID)]
sample.info.root.rarefy$bact<-1-sample.info.root.rarefy$Mito
rm(ps.root.rarefy.mito,mito.ratio)

write.table(sample.info.root.rarefy[with(sample.info.root.rarefy, order(Cycle_No,Cas9,gRNA,Plant_sample)),],
            file="./Table S5 sample_info_root_rarefy_bact.csv", sep = ",")
write.table(sample.info.root[with(sample.info.root, order(Cycle_No, Cas9,gRNA, Plant_sample)),],
            file="./Table S5-1 sample_info_root_raw.csv", sep=",")
## Make the plot
sample.info.root.rarefy
p.OTU.No<-ggplot(aes(y=Effective_OTU, x=interaction(Cas9, gRNA,Cycle_No), fill = Cas9), data=sample.info.root.rarefy)
p.OTU.No<-p.OTU.No+geom_boxplot()+theme_grey()
p.OTU.No<-p.OTU.No+labs(x="gRNA", y="Effective OTU", title="OTU number (Rarefied)")
p.OTU.No

p.ratio<-ggplot(aes(y=Mito*100, x=interaction(gRNA, Cycle_No), fill=Cas9),ylim(0,100), data=sample.info.root.rarefy)
p.ratio<-p.ratio+geom_boxplot()
p.ratio<-p.ratio+theme_grey()+labs(x="gRNA", y="Mitochondrion (%)", title="Mito ratio (rarefied)")
p.ratio
ggsave(filename="Figure5A rarefy_to_41k.pdf", plot_grid(p.ratio, p.OTU.No, rel_widths = c(1, 1)),
       width = 18, height = 6, 
       units = "cm")


############################################################################
#############  4.3 statistical test of microbiome diversity.################

#### 4.3.1 stitistical test of diversities from raw data     ################
library(agricolae)
sample.info.root$gRNA<-relevel(sample.info.root$gRNA,ref="No")
sample.info.root$Cycle_No<-as.factor(sample.info.root$Cycle_No)

## filtered sequences number
aov.test<-aov(filtered~gRNA+Cycle_No, data=sample.info.root)
summary(aov.test)
TukeyHSD(aov.test)
outHSD<-HSD.test(aov.test, "gRNA")
outHSD
outHSD<-HSD.test(aov.test, "Cycle_No")
outHSD

## total OTU frequencies
aov.test<-aov(Total~gRNA+Cycle_No, data=sample.info.root)
summary(aov.test)
TukeyHSD(aov.test)

# Mito relative abundance
aov.test<-aov(Mito~gRNA+Cycle_No, data=sample.info.root)
summary(aov.test)
TukeyHSD(aov.test)
outHSD<-HSD.test(aov.test,"gRNA")
outHSD
outHSD<-HSD.test(aov.test, "Cycle_No")
outHSD

# raw OTU number
aov.test<-aov(OTU_No~gRNA+Cycle_No, data=sample.info.root)
summary(aov.test)
TukeyHSD(aov.test)
outHSD<-HSD.test(aov.test,"gRNA")
outHSD

#############################################################################
## 4.3.2 stastical test of rarefied root dataset
sample.info.root.rarefy$gRNA<-relevel(sample.info.root.rarefy$gRNA,ref="No")
sample.info.root.rarefy$Cycle_No<-as.factor(sample.info.root.rarefy$Cycle_No)
# Mito relative abundance
aov.test<-aov(Mito~gRNA+Cycle_No, data=sample.info.root.rarefy)
summary(aov.test)
TukeyHSD(aov.test)
outHSD<-HSD.test(aov.test,"gRNA")
outHSD
outHSD<-HSD.test(aov.test, "Cycle_No")
outHSD
# Mito ratio between gRNAs
aov.test<-aov(Mito~gRNA+Cycle_No, data=sample.info.root.rarefy[sample.info.root.rarefy$Cas9=="Yes",])
summary(aov.test)
TukeyHSD(aov.test)

#OTU number
aov.test<-aov(Observed~gRNA+Cycle_No, data=sample.info.root.rarefy)
summary(aov.test)
TukeyHSD(aov.test)
outHSD<-HSD.test(aov.test,"gRNA")
outHSD

# effective OTU number
aov.test<-aov(Effective_OTU~gRNA+Cycle_No, data=sample.info.root.rarefy)
summary(aov.test)
TukeyHSD(aov.test)

# test variation of effective OTU number between mt-gRNA
aov.test<-aov(Effective_OTU~gRNA+Cycle_No, data=sample.info.root.rarefy[sample.info.root.rarefy$Cas9=="Yes",])
summary(aov.test)
TukeyHSD(aov.test)

# Shannon
aov.test<-aov(Shannon~gRNA+Cycle_No, data=sample.info.root.rarefy)
summary(aov.test)
TukeyHSD(aov.test)
outHSD<-HSD.test(aov.test,"gRNA")
outHSD

# Simpson
aov.test<-aov(Simpson~gRNA+Cycle_No, data=sample.info.root.rarefy)
summary(aov.test)
TukeyHSD(aov.test)
outHSD<-HSD.test(aov.test,"gRNA")
outHSD

## permutational ANOVA to test OTU diversities 
otu.matrix<-data.frame(otu_table(ps.root.noMt), check.names = FALSE)
head(otu.matrix)
otu.df<-t(otu.matrix)
colnames(otu.df)<-rownames(otu.matrix)

## make sure otu.df data have same order to sample info dateframe
otu.df<-otu.df[as.character(sample.info.root.rarefy$Seq_ID),]
adonis2(otu.df ~ Plant_sample+gRNA+Cycle_No, data=sample.info.root.rarefy, method="bray")
rm(otu.df, otu.matrix)




#####################################################
## 4.4 compare the abundance at phylum level
## Fig 5e. 
## 


pslog<-transform_sample_counts(ps.root.noMt, function(x) log(1+x))
out.log<-ordinate(pslog, method="NMDS", distance="bray")
p.ord<-plot_ordination(pslog, out.log, color="gRNA", shape="Plant_sample")+geom_point(size=I(2))

ps.root.phyla<-tax_glom(ps.root.noMt, taxrank = "Phylum")
ps.root.phyla<-transform_sample_counts(ps.root.phyla, function(x) {x/sum(x)} )
ps.root.phyla<-psmelt(ps.root.phyla)

ps.root.phyla<-ps.root.phyla[ps.root.phyla$Abundance>0.01,]
table(as.character(ps.root.phyla$Phylum))
summary(ps.root.phyla$Abundance)

## compare bacterial abundance at phylum level, Cas9+ vs Cas9-
aov.test<-aov(Abundance~gRNA+Cycle_No, data=ps.root.phyla[ps.root.phyla$Phylum=="Proteobacteria",])
summary(aov.test)
TukeyHSD(aov.test)
aov.test<-aov(Abundance~gRNA+Cycle_No, data=ps.root.phyla[ps.root.phyla$Phylum=="Bacteroidetes",])
summary(aov.test)
TukeyHSD(aov.test)
aov.test<-aov(Abundance~gRNA+Cycle_No, data=ps.root.phyla[ps.root.phyla$Phylum=="Firmicutes",])
summary(aov.test)
TukeyHSD(aov.test)
aov.test<-aov(Abundance~gRNA+Cycle_No, data=ps.root.phyla[ps.root.phyla$Phylum=="Actinobacteria",])
summary(aov.test)
TukeyHSD(aov.test)


write.table(file="Root.Phylum.comparison.txt",ps.root.phyla, sep="\t")

phylum_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861"
)

#Figure 5d
p.phylum<-ggplot(data=ps.root.phyla,aes(x=gRNA, y=Abundance, colour=Phylum, fill=Phylum)) + 
  geom_bar(stat = "identity")+facet_wrap(~Cycle_No + Plant_sample)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(x="Cas9",y="Rarefied Ralative Abundance (>1%)")
p.phylum
ggsave(filename = "Relative_abundance_Phyl_otu100_10_rarefy_41k.pdf", p.phylum,
       width = 5, height = 10, 
       units = "in")
rm(ps.root.phyla, p.phylum,phylum_colors)

###########################################################################
#  
#    5. Identify differential OTUs (DESeq2) 2019.8.13
###########################################################################
library(DESeq2)

## DESeq generates errors like ".. colData.. is not data.matrix, it is 
## cause by conflict with a package "agricolae" loaded in this environment.

#################################3
## 5.1 analyzing soil dataset

ps.soil.deseq2<-phyloseq_to_deseq2(ps.soil, ~ Cas9+gRNA)
ps.soil.deseq2

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

Deseq.analyze= function(dds, sig.alpha=0.05) {
  geoMeans<- apply(counts(dds),1, gm_mean)
  dds<- estimateSizeFactors(dds, geoMeans = geoMeans)
  dds<- DESeq(dds, fitType="local")
  res<-results(dds)
  #summary(res)
  sigtab = res[(res$padj < sig.alpha & !is.na(res$padj)), ]
  sigtab = cbind(gRNA=i, OTU=rownames(sigtab), as(sigtab, "data.frame"), as(tax_table(ps.soil)[rownames(sigtab), ], "matrix"))
  sigtab
}

## Analyzing the differential OTUs in soil samples
## each gRNA have three biological repeats.
   sig.alpha=0.01  # set significant p value
   ps.soil
   prune_taxa(rowSums(otu_table(ps.soil))>0,ps.soil)
   sample_data(ps.soil)
   gRNA_id=levels(sample_data(ps.soil)$gRNA)
   Dif.OTU<-data.frame();
for(i in gRNA_id) {
   dds<-prune_samples(sample_data(ps.soil)$gRNA == i,ps.soil)
   keep <- rowSums(otu_table(dds)) >= 10  #set the threshold to filter OTUs
   table(keep)
   dds <- prune_taxa(keep,dds)
   sample_data(dds)$Cas9<-relevel(sample_data(dds)$Cas9,ref="No")
   
   ps.dds<-phyloseq_to_deseq2(dds, ~Cas9)
   
   geoMeans<- apply(counts(ps.dds),1, gm_mean)
   dds<- DESeq2::estimateSizeFactors(ps.dds, geoMeans = geoMeans)
   dds<- DESeq(dds, test="Wald",fitType="parametric")
   res<-results(dds)
   #summary(res)
   #plotMA(res)
   #sigtab = res[(res$padj < sig.alpha & !is.na(res$padj)), ]
   sigtab = res[(res$pvalue < sig.alpha & !is.na(res$pvalue)), ]
   if(dim(DataFrame(sigtab))[1]>0) {
      sigtab = cbind(gRNA=i, OTU=rownames(sigtab), as(sigtab, "data.frame"), as(tax_table(ps.soil)[rownames(sigtab), ], "matrix"))
      sigtab
      #sigtab = Deseq.analyze(dds, sig.alpha)
      Dif.OTU<-rbind(Dif.OTU, sigtab, make.row.names=F)
   }
}
Dif.OTU
OTU.id<-matrix(unlist(strsplit(names(otu.seq),";size=",fixed=TRUE)), ncol=2, byrow=T)
names(otu.seq)<-OTU.id[,1]
Dif.OTU$Sequence<-as.character(otu.seq[Dif.OTU$OTU])

table(Dif.OTU$OTU)
write.table(file="Dif.OTU100.soil.threshold10.txt", Dif.OTU, sep="\t")
p<-ggplot(aes(y=Mito*100, x=interaction(gRNA, Cas9), fill=Cas9),ylim(0,100), data=sample.info.soil)
p<-p+geom_boxplot()
p<-p+theme_grey()
p
ggsave(filename = "Mito_soil_frequency.pdf", p,
       width = 15, height = 10, 
       units = "in")

rm(p,Dif.OTU, OTU.id)

##########################################################
##  5.2 Analyzing differential OTU in root samples
##
#  Analyze differential OTUs between Cas9+ vs - per gRNA for each sample
#  results. 8 and 15 cycles were used as technicalrepeat.
# Fig. 5d

sample_data(ps.root)
# recontruct ps.root.noMt using original ps.rool object
ps.root.noMt<-prune_taxa(!taxa_names(ps.root) %in% OTU.mt$OTU_ID, ps.root)
ps.root.noMt<-prune_taxa(data.frame(tax_table(ps.root))$Family != "Mitochondria", ps.root)

# remove low abundant reads
ps.root.noMt<-prune_taxa(taxa_sums(ps.root.noMt)>10, ps.root.noMt)
sample_id<-levels(sample_data(ps.root)$Plant_sample)
p_FCvsP<-ggplot();
p_FCvsPadj<-ggplot();
Dif.OTU<-data.frame() # output the result in this dataframe.
n=0; sig.alpha=0.01 # set to 0.05 if use padj
for(plant_sample in sample_id){  
   ps.sel<-prune_samples(sample_data(ps.root.noMt)$Plant_sample == plant_sample,ps.root.noMt)
   for (i in 1:(length(gRNA_id))) {
     ps.sel.gRNA<-prune_samples(sample_data(ps.sel)$gRNA %in% c("No", gRNA_id[i]), ps.sel)
     ## choose one prune method
     #ps.sel.gRNA<-prune_taxa(data.frame(tax_table(ps.sel.gRNA))$Family != "Mitochondria", ps.sel.gRNA)
     keep <- rowSums(otu_table(ps.sel.gRNA)) > 0
     table(keep)
     ps.sel.gRNA <- prune_taxa(keep,ps.sel.gRNA)
     sample_data(ps.sel.gRNA)$gRNA<-relevel(sample_data(ps.sel.gRNA)$gRNA,ref="No")
     dds<-phyloseq_to_deseq2(ps.sel.gRNA, ~gRNA)
     geoMeans<- apply(counts(dds),1, gm_mean)
     dds<- estimateSizeFactors(dds, geoMeans = geoMeans)
     dds<- DESeq(dds, fitType="parametric", test="Wald")
     res<-results(dds)
     #res["OTU_446",]
     #sigtab = res[(res$padj < sig.alpha & !is.na(res$padj)), ]
     sigtab = res[(res$pvalue<sig.alpha & ! is.na(res$pvalue)),]
     sigtab = cbind(gRNA=gRNA_id[i], OTU=rownames(sigtab), as(sigtab, "data.frame"), 
                    as(tax_table(ps.root.noMt)[rownames(sigtab), ], "matrix"))
     sigtab$Plant_sample=plant_sample     
     Dif.OTU<-rbind(Dif.OTU,sigtab)
   #summary(res)
   #plotMA(res)
   #col=res$padj>0
   p_FCvsP[[i+n]]<-ggplot(data=data.frame(res), aes(x=log2FoldChange, y=-log10(pvalue)))  
   p_FCvsP[[i+n]]<-p_FCvsP[[i+n]]+geom_point(aes(colour=pvalue>0.01),size=I(2),alpha=0.7)
   p_FCvsP[[i+n]]<-p_FCvsP[[i+n]]+geom_vline(xintercept=0, color= "grey")
   p_FCvsP[[i+n]]<-p_FCvsP[[i+n]]+labs(title = paste(plant_sample,"\n",  gRNA_id[i], " vs CK"))
   p_FCvsPadj[[i+n]]<-ggplot(data=data.frame(res), aes(x=log2FoldChange,y=-log10(padj)))
   p_FCvsPadj[[i+n]]<- p_FCvsPadj[[i+n]]+geom_point(aes(colour=padj>0.01),size=I(2),alpha=0.7)
   p_FCvsPadj[[i+n]]<- p_FCvsPadj[[i+n]]+geom_vline(xintercept=0, color= "grey")
   p_FCvsPadj[[i+n]]<- p_FCvsPadj[[i+n]]+labs(title = paste(plant_sample,"\n",gRNA_id[i], " vs CK "),hjust=0.5)
   }
   n=n+3
}    
dim(Dif.OTU)
Dif.OTU$Sequence<-as.character(otu.seq[Dif.OTU$OTU])
Dif.OTU[Dif.OTU$log2FoldChange< 0,] # no OTU have LogFC < 0

# Get the abundance of differential OTUs LFC in case use ps.root as input
unique(Dif.OTU$OTU)  # 29 OTUs in total.
write.table(otu_table(ps.root)[unique(as.character(Dif.OTU$OTU))], file="Dif_OTU_NoMt_rawReads.txt", sep="\t")
otu_table(ps.root.noMt)[unique(as.character(Dif.OTU$OTU))]/sample_sums(ps.root.noMt)
summary(otu_table(ps.root)[,1]/sample_sums(ps.root.noMt)[1])
write.table(otu_table(ps.root)[unique(as.character(Dif.OTU$OTU))]/sample_sums(ps.root.noMt), file="Dif_OTU_NoMt_RA.txt", sep="\t")

t<-as.character(unique(Dif.OTU$OTU[Dif.OTU$log2FoldChange>0]))
t<-otu.matrix[t,]
colSums(t)/colSums(otu.matrix)  #relative abundance of all differential OTU logFC>1.
cbind(sample.info[rownames(t),1:6],t)

write.table(Dif.OTU, file="Table S6.Dif_OTU_Root_pval0.01.csv", sep = ",", row.names = F)    #-w-ChrMr to store resultes including Mt OTU in analysis


p.all<-plot_grid(p_FCvsP[[1]],p_FCvsP[[4]],p_FCvsP[[7]],
                 p_FCvsP[[2]],p_FCvsP[[5]],p_FCvsP[[8]],
                 p_FCvsP[[3]],p_FCvsP[[6]],p_FCvsP[[9]],nrow=3)

ggsave(filename = "Fig. 5DDif_OTU100_root_LFCvsPvalue_NoMt.pdf", p.all,
       width = 15, height = 10, 
       units = "in")


## Done!!!