## This short R script is used to analyze chimera frequencies in root samples
## The chimeric sequences were detected using vsearch --uchime_denovo function
## 
## Kabin Xie, 2020.5.1


library(ggplot2)
library(cowplot)
library(agricolae)
# load the sample information data
vs.summary<-read.table(file="../vs.pipline.summary100.csv", head=TRUE, sep=",")
sample.info<-read.table("../vsearch_sample_info.csv", head=T, sep=",")
sample.info<-cbind(sample.info, vs.summary[rownames(sample.info), c(2,3,4,5)])
# Load the chimera table (same to OTU tables)
# remove the "#" symbol and replace "OTU ID" by OTU_ID in the first line before loading
chimeras.tab<-read.table(file="../vsearch.out/all.chimeras100.txt", head=T, check.names=FALSE)

# Check the Table
head(chimeras.tab)
chimera.abund<-(rowSums(chimeras.tab[,-1]))
table(chimera.abund)

# Remove low abundant noise, if neccessary.
# 
chi_No<-colSums(chimeras.tab[chimera.abund>0,-1])
sum(chimera.abund)
sum(chi_No) 

# combine the sample.info data and chimeras table
sample.info<-cbind(sample.info, chi_No=chi_No[sample.info$Seq_ID])
sample.info$chi_ratio<-sample.info$chi_No/sample.info$filtered
sample.info=sample.info[sample.info$Seq_ID!="bS1_4",]  #remove the empty control

# Comparing data in the root dataset
sample.info.root<-sample.info[sample.info$Compartment=="Root",]

# ANOVA testing
# Higher quality sequences ratio
sample.info.root$Cycle_No<-as.factor(sample.info.root$Cycle_No)

aov.test<-aov(filtered/merged~gRNA, data=sample.info.root[sample.info.root$Cycle_No == 15,])
summary(aov.test)
aov.test<-aov(filtered/merged~gRNA, data=sample.info.root[sample.info.root$Cycle_No == 8,])
summary(aov.test)

aov.test<-aov(filtered/merged~Cycle_No, data=sample.info.root)
summary(aov.test)
aov.test<-aov(filtered/merged~gRNA+Cycle_No, data=sample.info.root)
summary(aov.test)
TukeyHSD(aov.test)

# chimera rate
aov.test<-aov(chi_ratio~gRNA, data=sample.info.root[sample.info.root$Cycle_No == 15,])
summary(aov.test)
aov.test<-aov(chi_ratio~gRNA, data=sample.info.root[sample.info.root$Cycle_No == 8,])
summary(aov.test)
TukeyHSD(aov.test)

outHSD<-HSD.test(aov.test,"gRNA")
outHSD
aov.test<-aov(chi_ratio~gRNA+Cycle_No, data=sample.info.root)
summary(aov.test)
TukeyHSD(aov.test)
outHSD<-HSD.test(aov.test,"gRNA")
outHSD

# calculate mean and sd for each treatment
chimera.table<-aggregate(chi_ratio~gRNA+Cycle_No, data=sample.info.root, mean)
chimera.table$chi_sd<-aggregate(chi_ratio~gRNA+Cycle_No, data=sample.info.root, sd)$chi_ratio
chimera.table$filtered<-aggregate(filtered/merged~gRNA+Cycle_No, data=sample.info.root, mean)[,3]
chimera.table$filtered.sd<-aggregate(filtered/merged~gRNA+Cycle_No, data=sample.info.root, sd)[,3]

# output the data
write.table(file= "chimera_summmary.csv", chimera.table, sep=",")
write.csv(file="chimera_data.csv",sample.info, row.names = F)

# Draw the figures
p.filter<-ggplot(aes(y= 100*filtered/merged, x=interaction(gRNA, Cycle_No), fill = Cas9),
                 data=sample.info.root,ylim=c(0,100))
p.filter<-p.filter+geom_boxplot()+theme_grey()+ylim(0,100)
p.filter<-p.filter+labs(title ="Qual filtered sequences", y= "Filtered sequences (%)")

p.chimera<-ggplot(aes(y=100*chi_No/merged, x=interaction(gRNA, Cycle_No), fill=Cas9), 
                ylim(0,100),data=sample.info.root)
p.chimera<-p.chimera+geom_boxplot()+theme_grey()+ylim(0,100);
p.chimera<-p.chimera+labs(title = "Chimeric sequences", y="Chimera sequences (%)")

plot_grid(p.filter, p.chimera, rel_widths = c(1, 1), ncol=1)
ggsave(filename="Fig.Chimera.pdf", plot_grid(p.filter, p.chimera, rel_widths = c(1, 1), ncol=1),
       width = 8.5, height = 11, 
       units = "in")
