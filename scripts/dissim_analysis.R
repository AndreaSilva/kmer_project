# NOTE! All paths are relative to the $PROJECT path.
# Script for processing of real metagenomes;
# include: normalisation of k-mer frequencies, counting Bray-Curtis Dissimilarity, and save matrices

rm(list= ls(all.names = TRUE))
# Includes file with methods
source("functions.R")
library(ecodist) 
library(scales)
library(vegan)
library(data.table)
library(reshape)

setwd("/data2/bio/uploads/Dubinkina_2015_suppl_data/scripts")

###########################    calculate and save k-mer dissimilarities     ##########################
method<-"BC"   # used metric BC or COR
data_dir<-"/data2/bio/uploads/Dubinkina_2015_suppl_data/data/Kmers/"    # path to files with k-mers 
subset_names<-c("synt_diff", "synt_similar")    # sim experiment
time<-data.frame()
for (i in 5:12){        # in range of K
  time<-rbind(time,calc_dissimilarity(data_dir, subset_names[1], i, method))
}
subset_names<-c("China", "China_tax_mapped", "China_fun_mapped", "USA", "USA_tax_mapped", "USA_fun_mapped")
calc_dissimilarity(data_dir, subset_names, 11, method)    # for real datasets 

#####################         Additional Figure 2 (comparison of kmer and TAX for simulations)         ##############################
################ redraw par wth for
load("../output/dist_BC_synt_diff_K11.RData")
load("../data/synt_diff_TAX_org.RData")
BC_TAX_diff<-bcdist(t(Bac_data))
plotPath<-paste("../output/synt_diff_BC_TAX_vs_BC_kmer_K11.pdf", sep="")
pdf(plotPath, height= 5, width = 5)
plot(BC_TAX_diff, Dist, col=alpha("red",0.3), pch=16, xlim=c(0, 1), ylim=c(0, 0.4), xlab="BC TAX", ylab="BC kmer" )
cor_t_k<-paste("r = ",round(cor(BC_TAX_diff, Dist, method="spearman"),digits=2), sep="")
mantel(as.matrix(BC_TAX_diff), as.matrix(Dist), permutations=999, method="spearman")
legend("bottomright",  cor_t_k, pch = 16, pt.bg = "white", lty = 0, col = "white", bty="n")
dev.off()

load("../output/dist_BC_synt_similar_K11.RData")
load("../data/synt_similar_TAX_org.RData")
BC_TAX_similar<-bcdist(t(Bac_data))
plotPath<-paste("../output/synt_similar_BC_TAX_vs_BC_kmer_K11.pdf", sep="")
pdf(plotPath, height= 5, width = 5)
plot(BC_TAX_similar, Dist, col=alpha("blue",0.3), pch=16, xlim=c(0, 1), ylim=c(0, 0.4), xlab="BC TAX", ylab="BC kmer" )
cor_t_k<-paste("r = ",round(cor(BC_TAX_similar, Dist, method="spearman"),digits=2), sep="")
mantel(as.matrix(BC_TAX_similar), as.matrix(Dist), permutations=999, method="spearman")
legend("bottomright",  cor_t_k, pch = 16, pt.bg = "white", lty = 0, col = "white", bty="n")
dev.off()

#####################    Additional Figure 3 graph correlation and calculation time for all k     #############################
A<-vector()
k=0
for (i in 5:12){
  k=k+1
  file_name<-paste("../output/dist_BC_synt_diff_K", i,".RData", sep="")
  load(file=file_name)
  A[k]<-cor(BC_TAX_diff, Dist, method="spearman")
}

B<-vector()
k=0
for (i in 5:12){
  k=k+1
  file_name<-paste("../output/dist_BC_synt_similar_K", i,".RData", sep="")
  load(file=file_name)
  B[k]<-cor(BC_TAX_similar, Dist, method="spearman")
}
C<-c(5,6,7,8,9,10,11,12)
plotPath<-paste("../output/AdditionalFigure3.pdf", sep="")
pdf(plotPath, height= 5, width = 5)
par(mar=c(5,5,4,5))
plot(C, B, col=alpha("blue"), pch=15, ylim=c(0.4, 1), xlab="Value of k", ylab="Value of Spearman correlation coefficient r" )
lines(C,B, col='blue', lwd=2)
points(C,A, col=alpha("red"), pch=16,)
lines(C,A, col='red', lwd=2)
norm_time<-0.4+rowMeans(time)/rowMeans(time)[8]*0.6
points(C,norm_time, col=alpha("green"), pch=17,)
lines(C, norm_time, col='green', lwd=2)
axis(4, at=seq(0.4,1,by=0.1), labels=seq(0,770/60, by=2))
mtext("Calculation time in minutes", side=4, line=2,)
legend("topleft", c("High-diversity", "Low-diversity", "Time") , pch = c(16,15,17), pt.bg = "white", lty = 0, col = c("red", "blue", "green"), bty="n")
dev.off()

#############################           Additional Figure 1 snp comparison             #######################
load("../data/synt_diff_SNP_simulation_K10_dist.RData")
load("../data/synt_diff_TAX_org_SNP_simulation.RData")

# calculate taxonomic bc distance matrix
tax.d <- bcdist(tax.m)

# select projects
kmer.prjs <- c("nmut", "u001", "u01", "p001", "p01")

# figure params
x.lb <- "BC Tax" # x axis lab
y.lb <- "BC kmer" # y axis lab
x.lim <- c(0, 1) # x axis limits
y.lim <- c(0, 0.3) # y axis limits
round.cor <- 3 # round correlation digits

## Uniform figure

pdf("../output/UniformFigure1.pdf", width = 9, height = 3.3)
par(ps = 15)

# select projects to plot and set titles
plot.nms <- c("nmut" = "No SNP",
              "u001" = "Mean SNP rate: 1%", 
              "u01" = "Mean SNP rate: 10%")
layout(matrix(1:length(plot.nms),
              ncol = length(plot.nms),
              byrow = T))

for (cprj in names(plot.nms)) {
  plot(x = tax.d, y = kmer.dist.list[[cprj]],col=alpha("red",0.1), pch=16,
       xlab = x.lb, ylab = y.lb,
       xlim = x.lim, ylim = y.lim,
       main = plot.nms[cprj])
  
  cor.sp <- round(cor(x = tax.d, y = kmer.dist.list[[cprj]],
                      method = "spearman"), d = round.cor) 
  legend("topleft",
         legend = substitute(expression(r == cval), list(cval = cor.sp))[[2]],
         bty = 'n')
}
dev.off()

## Phylogenetic figure

cairo_pdf("../output/AdditionalFigure1.pdf", width = 9, height = 3.3)
par(ps = 15)

# select projects to plot and set titles
plot.nms <- c("nmut" = "No SNP",
              "p001" = "Mean SNP rate: 1%", 
              "p01" = "Mean SNP rate: 10%")
layout(matrix(1:length(plot.nms),
              ncol = length(plot.nms),
              byrow = T))

for (cprj in names(plot.nms)) {
  plot(x = tax.d, y = kmer.dist.list[[cprj]], col=alpha("red",0.1), pch=16,
       xlab = x.lb, ylab = y.lb,
       xlim = x.lim, ylim = y.lim,
       main = plot.nms[cprj])
  
  cor.sp <- round(cor(tax.d, kmer.dist.list[[cprj]],
                      method = "spearman"), d = round.cor) 
  legend("topleft",
         legend = substitute(expression(r == cval), list(cval = cor.sp))[[2]],
         bty = 'n')
}
dev.off()

# plotPath<-paste("../output/AdditionalFigure1.pdf", sep="")
# pdf(plotPath, height= 5, width = 9)
# par(mfrow=c(1,2))
# plot(tax.d, kmer.d, col=alpha("red",0.1), pch=16, xlim=c(0, 1), ylim=c(0, 0.3), xlab="BC TAX", ylab="BC kmer", main="Not mutated genomes" )
# cor_t_k<-paste("r = ",round(cor(tax.d,kmer.d, method="spearman"),digits=2), sep="")
# mantel(as.matrix(tax.d), as.matrix(kmer.d), permutations=999, method="spearman")
# legend("bottomright",  cor_t_k, pch = 16, pt.bg = "white", lty = 0, col = "white", bty="n")
# 
# load("../data/Kmers/synt_diff_mutated_K10.RData")
# kmer.d<-bcdist(kmer.m.m)
# plot(tax.d, kmer.d, col=alpha("red",0.1), pch=16, xlim=c(0, 1), ylim=c(0, 0.3), xlab="BC TAX", ylab="BC kmer", main="Mutated genomes"  )
# cor_t_k<-paste("r = ",round(cor(tax.d,kmer.d, method="spearman"),digits=2), sep="")
# mantel(as.matrix(tax.d), as.matrix(kmer.d), permutations=999, method="spearman")
# legend("bottomright",  cor_t_k, pch = 16, pt.bg = "white", lty = 0, col = "white", bty="n")
# 
# dev.off()

###################   load and calculate reference-based dissimilarities and draw pcoa for all methods  ###############
load("../output/dist_BC_all_K11.RData")
p1<-draw_pcoa(Dist, "kmer")
p1
load("../data/TAX_org_composition.RData")
d.BC_TAX_org<-bcdist(TAX_org)
p2<-draw_pcoa(d.BC_TAX_org, "TAX_org")
p2
load("../data/TAX_genus_composition.RData")
# d.BC_TAX_genus<-bcdist(TAX_genus)
# draw_pcoa(d.BC_TAX_genus, "TAX_genus")

load("../data/COG_composition.RData")
d.BC_COG<-bcdist(COG)
p3<-draw_pcoa(d.BC_COG, "COG")
p3

load("../data/MetaPhlAn_composition.RData")
d.BC_MetaPhlAn<-bcdist(Metaphlan)
p4<-draw_pcoa(d.BC_MetaPhlAn, "MetaPhlAn")
p4

load("../data/WG_UniFrac.RData")
p5<-draw_pcoa(as.dist(WG_Unifrac), "WG_UniFrac")
p5

###############     compare k-mer with different reference-based methods      ################
### k-mer
load("../output/dist_BC_USA_K11.RData")
df.all_dist_USA <- as.data.frame(as.vector(Dist))
colnames(df.all_dist_USA) <- "d_BC_kmer"
USA<-colnames(as.matrix(Dist))
load("../output/dist_BC_China_K11.RData")
df.all_dist_China <- as.data.frame(as.vector(Dist))
colnames(df.all_dist_China) <- "d_BC_kmer"
China<-colnames(as.matrix(Dist))
### COG
df.all_dist_USA$d_BC_COG <- bcdist(COG[which(rownames(COG)%in%USA),])
df.all_dist_China$d_BC_COG <- bcdist(COG[which(rownames(COG)%in%China),])
### TAX org
df.all_dist_USA$d_BC_TAX_org <- bcdist(TAX_org[which(rownames(TAX_org)%in%USA),])
df.all_dist_China$d_BC_TAX_org <- bcdist(TAX_org[which(rownames(TAX_org)%in%China),])
### TAX genus
df.all_dist_USA$d_BC_TAX_genus <- bcdist(TAX_genus[which(rownames(TAX_genus)%in%USA),])
df.all_dist_China$d_BC_TAX_genus <- bcdist(TAX_genus[which(rownames(TAX_genus)%in%China),])
### MetaPhlAn
df.all_dist_USA$d_BC_MetaPhlAn <- bcdist(Metaphlan[which(rownames(Metaphlan)%in%USA),])
df.all_dist_China$d_BC_MetaPhlAn <- bcdist(Metaphlan[which(rownames(Metaphlan)%in%China),])
### WG UniFrac
df.all_dist_USA$d_WG_UniFrac <- as.dist(WG_Unifrac[which(rownames(WG_Unifrac)%in%USA),which(colnames(WG_Unifrac)%in%USA)])
df.all_dist_China$d_WG_UniFrac <- as.dist(WG_Unifrac[which(rownames(WG_Unifrac)%in%China),which(colnames(WG_Unifrac)%in%China)])

# correlation
cor_all_vs_all <- matrix(0,6,6) # matrix(0,7,7)
colnames(cor_all_vs_all) <- c("BC k-mer", "BC COG", "BC TAX (org)", "BC TAX (genus)",
                            "BC MetaPhlAn (org)", "WG UniFrac")
rownames(cor_all_vs_all) <- c("BC k-mer", "BC COG", "BC TAX (org)", "BC TAX (genus)",
                            "BC MetaPhlAn (org)", "WG UniFrac")
for (i in 1:ncol(cor_all_vs_all)){
  #i=1
  for (j in 1:i){
    cor_all_vs_all[i,j]<-cor(df.all_dist_USA[,i], df.all_dist_USA[,j], method="spearman")
  }
}   # SRS
zzz
for (i in 1:ncol(cor_all_vs_all)){
  #i=1
  for (j in i:ncol(cor_all_vs_all)){
    cor_all_vs_all[i,j]<-cor(df.all_dist_China[,i], df.all_dist_China[,j], method="spearman")
  }
}   # SRR
diag(cor_all_vs_all)<-0
data_table_labels<-round(cor_all_vs_all, digits=2)
data_table<-as.data.frame(cor_all_vs_all)
library(gplots)
library(spatstat)
#### plot result as table
plotPath<-paste("../output/table_correlation_metrics.pdf", sep="")
pdf(plotPath, height= 5, width = 6)
heatmap.2(as.matrix(data_table),dendrogram="none", Rowv=F, Colv=F,cellnote=as.matrix(data_table_labels),
          notecol="black",col=redblue(100),scale="none",key=TRUE, keysize=1.5,
          density.info="none", trace="none", cexRow=0.7,cexCol=0.7, srtRow=30,srtCol=30, mar=c(5,7))  # redraw with ggplot heatmap
dev.off()

###################################  Additional Figure 4  graphs metric vs metric    #########################################
load("../data/Map_stats.RData")
plotPath<-paste("../output/boxplot_mapping_stats.pdf", sep="")
pdf(plotPath, height= 5, width = 6)
par(mfrow=c(1,2))
boxplot(as.vector(Map_percent[which(rownames(Map_percent)%in%USA),1]), 
        as.vector(Map_percent[which(rownames(Map_percent)%in%China),1]), col=c("#CC79A7","#009E73"), 
        names=c("USA", "China"), ylim=c(1,100), main="TAX")     # taxonomic
boxplot(as.vector(Map_percent[which(rownames(Map_percent)%in%USA),2]), 
        as.vector(Map_percent[which(rownames(Map_percent)%in%China),2]), col=c("#CC79A7","#009E73"),
        names=c("USA", "China"), ylim=c(1,100),main="COG")      # functional
dev.off()

#######################   Figure 2 graphs for functional and taxonomic k-mers    #########################

plotPath<-paste("../output/Figure2.pdf", sep="")
pdf(plotPath, height= 8, width = 8)
par(mfrow=c(2,2))
### TAX
m1_vs_m2(df.all_dist_USA$d_BC_TAX_org, df.all_dist_USA$d_BC_kmer,
         df.all_dist_China$d_BC_TAX_org, df.all_dist_China$d_BC_kmer, "BC TAX org", "BC kmer", 1)
load("../output/dist_BC_USA_K11.RData")
A<-as.matrix(bcdist(TAX_org[USA,]))
B<-as.matrix(Dist)
AB.m<-cbind(melt(A),melt(B))
AB.mean<-AB.m[,c(3,6)]
AB.mean<-AB.mean[-which(rowSums(AB.mean)==0),]
fit<-lm(AB.mean[,2]~AB.mean[,1])
abline(fit, col="black")
### COG
m1_vs_m2(df.all_dist_USA$d_BC_COG, df.all_dist_USA$d_BC_kmer,
         df.all_dist_China$d_BC_COG, df.all_dist_China$d_BC_kmer, "BC COG", "kmer", 1)
load("../output/dist_BC_USA_K11.RData")
A<-as.matrix(bcdist(COG[USA,]))
B<-as.matrix(Dist)
AB.m<-cbind(melt(A),melt(B))
AB.mean<-AB.m[,c(3,6)]
AB.mean<-AB.mean[-which(rowSums(AB.mean)==0),]
fit<-lm(AB.mean[,2]~AB.mean[,1])
abline(fit, col="black")

### TAX
load("../output/dist_BC_USA_tax_mapped_K11.RData")
dist_USA_tax_mapped <- Dist
load("../output/dist_BC_China_tax_mapped_K11.RData")
dist_China_tax_mapped <- Dist
m1_vs_m2(df.all_dist_USA$d_BC_TAX_org, dist_USA_tax_mapped,
         df.all_dist_China$d_BC_TAX_org, dist_China_tax_mapped, "BC TAX org", "kmer", 1)
### COG
load("../output/dist_BC_USA_fun_mapped_K11.RData")
dist_USA_fun_mapped <- Dist
load("../output/dist_BC_China_fun_mapped_K11.RData")
dist_China_fun_mapped <- Dist
m1_vs_m2(df.all_dist_USA$d_BC_COG, dist_USA_fun_mapped,
         df.all_dist_China$d_BC_COG, dist_China_fun_mapped, "BC COG", "kmer", 1)
dev.off()

####################          Additional Figure 4 metric with metric 3 remaining     ######################
plotPath<-paste("../output/AdditionalFigure4.pdf", sep="")
pdf(plotPath, height= 4, width = 12)
par(mfrow=c(1,3))
m1_vs_m2(df.all_dist_USA$d_BC_TAX_genus, df.all_dist_USA$d_BC_kmer,
         df.all_dist_China$d_BC_TAX_genus, df.all_dist_China$d_BC_kmer, "BC TAX genus", "BC kmer", 1)
m1_vs_m2(df.all_dist_USA$d_BC_MetaPhlAn, df.all_dist_USA$d_BC_kmer,
         df.all_dist_China$d_BC_MetaPhlAn, df.all_dist_China$d_BC_kmer, "BC MetaPhlAn", "BC kmer", 1)
m1_vs_m2(df.all_dist_USA$d_WG_UniFrac, df.all_dist_USA$d_BC_kmer,
         df.all_dist_China$d_WG_UniFrac, df.all_dist_China$d_BC_kmer, "WG UniFrac", "BC kmer", 0.2)
dev.off()

#######################          graphs for phages experiment            #########################
load("../data/phage_USA.RData")
n.phage<-rownames(phage)[1:6]
n.not_phage<-rownames(phage)[7:nrow(phage)]

load("../output/dist_BC_USA_K11.RData")
m.dist_USA<-as.matrix(Dist)
m.BC_TAX_org_USA<-as.matrix(bcdist(TAX_org[which(rownames(TAX_org)%in%USA),]))

m.mycol <- matrix(rep("grey", nrow(phage)*nrow(phage)), nrow(phage), nrow(phage))
rownames(m.mycol)  <- rownames(m.BC_TAX_org_USA)
colnames(m.mycol)  <- colnames(m.BC_TAX_org_USA)
m.mycol[which(rownames(m.mycol)%in%n.phage),  which(colnames(m.mycol)%in%n.not_phage)]<-"red"
m.mycol[which(rownames(m.mycol)%in%n.not_phage), which(colnames(m.mycol)%in%n.phage)]<-"red"
ind <- which(m.mycol == "red")
par(mfrow=c(1,1))
plotPath<-paste("../output/Phage.pdf", sep="")
pdf(plotPath, height= 5, width = 5)
plot(m.BC_TAX_org_USA, m.dist_USA, col=alpha(m.mycol, 0.1), pch=17,
     xlab="BC TAX org", ylab="BC kmer", xlim=c(0,1), ylim=c(0,0.7))
points(m.BC_TAX_org_USA[ind], m.dist_USA[ind], col=alpha(m.mycol[ind],0.7), pch=17)
legend("bottomright", c("USA all","USA with/without crAss") , pch = 17, pt.bg = "white", lty = 0, col = c("grey","red"), bty = "n")
dev.off()


########## add pie chart

######## chi-square test
M<-cbind(c(738,8),c(6864,8142))
chisq.test(M)

####################        check percentage of mapping         ######################
#################### bacterial composition
load("../data/Map_stats.RData")
par(mfrow=c(1,2))
#### TAX 
boxplot(as.vector(Map_percent[USA,1]), as.vector(Map_percent[China,1]), col=c("blue", "red"), names=c("USA", "China"), ylim=c(1,100), main="Taxonomic")
#### COG
boxplot(Map_percent[USA,2], Map_percent[China,2], col=c("blue", "red"), names=c("USA", "China"), ylim=c(1,100),main="Functional")
wilcox.test(Map_percent[,1], Map_percent[,2], alternative = "less", paired = T) 





########### tmp lm analysis for USA k-mers from mapped reads
load("../output/dist_BC_USA_K11.RData")
A<-as.matrix(bcdist(TAX_org[USA,]))
B<-as.matrix(Dist)
library(reshape)
AB.m<-cbind(melt(A),melt(B))

AB.mean<-AB.m[,c(3,6)]
AB.mean<-AB.mean[-which(rowSums(AB.mean)==0),]
fit<-lm(AB.mean[,2]~AB.mean[,1])
abline(fit, col="black")
#y=0.1492x+0.1381

subset<-AB.m[which(AB.m[,6]>0.4),]

load("../output/dist_BC_USA_tax_mapped_K11.RData")
tmp<-as.matrix(Dist)
rownames(tmp)<-gsub("America_", "", rownames(tmp))
rownames(tmp)<-gsub(".reference-mapped.sorted.bam.reference-mapped.fasta.11", "", rownames(tmp))
colnames(tmp)<-rownames(tmp)
AB.map.m<-melt(tmp)

y.assesed<-0.1492*AB.m[subset[,1]:subset[,2],3]+0.1381
deltas1<-(AB.m[subset[,1]:subset[,2],6]-y.assesed)
d1<-mean(deltas1)
d1
sd1<-sd(deltas1)
deltas2<-(AB.map.m[subset[,1]:subset[,2],3]-y.assesed)
d2<-mean(deltas2)
d2
sd2<-sd(deltas2)

d2/d1
wilcox.test(deltas1, deltas2, alternative = "two.sided", paired = T) 
mean((deltas1-deltas2)/deltas1)
sd((deltas1-deltas2)/deltas1)

############################## similar analysis for functional mapping    #######
load("../output/dist_BC_USA_K11.RData")
A<-as.matrix(bcdist(COG[USA,]))
B<-as.matrix(Dist)
library(reshape)
AB.m<-cbind(melt(A),melt(B))

AB.mean<-AB.m[,c(3,6)]
AB.mean<-AB.mean[-which(rowSums(AB.mean)==0),]
fit<-lm(AB.mean[,2]~AB.mean[,1])
abline(fit, col="black")
#y=0.71375*x+0.06469

subset<-AB.m[which(AB.m[,6]>0.4),]

load("../output/dist_BC_USA_fun_mapped_K11.RData")
tmp<-as.matrix(Dist)
rownames(tmp)<-gsub("America_", "", rownames(tmp))
rownames(tmp)<-gsub(".reference-mapped.sorted.bam.reference-mapped.fasta.11", "", rownames(tmp))
colnames(tmp)<-rownames(tmp)
AB.map.m<-melt(tmp)

y.assesed<-0.1492*AB.m[subset[,1]:subset[,2],3]+0.1381
deltas1<-(AB.m[subset[,1]:subset[,2],6]-y.assesed)
d1<-mean(deltas1)
d1
sd1<-sd(deltas1)
sd1
deltas2<-(AB.map.m[subset[,1]:subset[,2],3]-y.assesed)
d2<-mean(deltas2)
d2
sd2<-sd(deltas2)
sd2
wilcox.test(deltas1, deltas2, alternative = "two.sided", paired = T) 
mean((deltas2-deltas1)/deltas1)
sd((deltas2-deltas1)/deltas1)
