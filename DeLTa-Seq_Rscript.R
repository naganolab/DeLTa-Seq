#### Fig.2a
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbeeswarm)
data = read.table("qPCR/reducing_reagent_RTcheck.txt", fill = T, sep = "\t", header = T,stringsAsFactors = T)
data$condition = factor(data$condition, levels = c("RT(-)","DTT 5 mM","DTT 10 mM","DTT 100 mM", "2Me 1%", "2Me 2.5%", "2Me 5%", "2Me 10%"))
out = NULL
for (con in unique(data$condition)) {
  tmp = data[data$condition==con,]
  mean = mean(as.numeric(tmp$Cp))
  # tmp = tmp[diff!=max(diff),]
  tmp2 = cbind(data.frame(condition = tmp$condition), data.frame(Ct =  as.numeric(tmp$Cp)))
  out = rbind(out, tmp2)
}
mean = out %>%
  group_by(condition) %>%
  summarise(mean=mean(Ct))
out$Ct = mean[1,]$mean-out$Ct
mean$mean =mean[1,]$mean-mean$mean 


orange = colorRampPalette(c("white", "orange"))(100)
green = colorRampPalette(c("white", "green"))(100)
pdf("Fig2a.pdf", width = 17)
ggplot(out, aes(x = condition, y = Ct, group=condition, color = condition)) +
  geom_beeswarm(size=4, cex =4, color = "black",alpha = 0)+
  theme_bw() + 
  theme(legend.position = "bottom")+
  scale_color_manual(values = c("gray",green[c(30,60,100)],orange[c(30,40,70,100)]))+
  geom_spoke(data=mean,aes(x=as.numeric(condition)-0.5, y=mean),angle = 0, radius = 1, size = 2)+
  geom_beeswarm(size=4, cex =4, color = "black")
dev.off()

#### Fig. 2b-e,g
Ath.act2 = read.table("qPCR/Ath_act2_qPCR.txt", header = T, sep = "\t")
Ath.act2$Samples = factor(Ath.act2$Samples, levels = c("NC","DTT 10 mM", "DTT 50 mM", "DTT 100 mM","2Me 2.5%"))
tmp = data.frame(DiffCt = -(Ath.act2$MeanCp-mean(Ath.act2$MeanCp[Ath.act2$Samples=="NC"])))
Ath.act2= cbind(Ath.act2, tmp)
mean = Ath.act2 %>% group_by(Samples) %>% summarise(mean = mean(DiffCt))
pdf("Fig2b.pdf", width = 5)
ggplot(Ath.act2, aes(x = Samples, y = DiffCt, group=Samples, color = Samples)) +
  geom_beeswarm(size=4, cex =4, alpha =0)+
  scale_y_continuous(breaks = -1:11, limits = c(-0.5,10.5))+
  scale_color_manual(values = c("gray",green[c(30,60,100)],orange[c(40)]))+ 
  theme_bw() +
  ylab("Ct")+
  geom_spoke(data=mean,aes(x=as.numeric(Samples)-0.5, y=mean),angle = 0, radius = 1, size = 2)+
  theme(legend.position = "NA")+
  geom_beeswarm(size=4, cex =4, color = "black")
dev.off()


Osa.act = read.table("qPCR/Osa_act_qPCR.txt", header = T, sep = "\t")
Osa.act = Osa.act[Osa.act$Condition!="Std",]
Osa.act$Condition = factor(Osa.act$Condition, levels = c("NC","DTT 10 mM", "DTT 50 mM", "DTT 100 mM","2Me 2.5%","Std"))
mean.leaf = Osa.act[Osa.act$tissue=="leaf",] %>% group_by(Condition) %>% summarise(mean = mean(MeanCp))
mean.root = Osa.act[Osa.act$tissue=="root",] %>% group_by(Condition) %>% summarise(mean = mean(MeanCp))
Osa.act.leaf = Osa.act[Osa.act$tissue=="leaf",]
Osa.act.root = Osa.act[Osa.act$tissue=="root",]
Osa.act.leaf$MeanCp = mean.leaf$mean[1]-Osa.act.leaf$MeanCp
Osa.act.root$MeanCp = mean.root$mean[1]-Osa.act.root$MeanCp

mean.leaf[,2] = mean.leaf$mean[1] -mean.leaf[,2]
mean.root[,2] = mean.root$mean[1]-mean.root[,2] 
mean.root$Condition = factor(mean.root$Condition)

pdf("Fig2c.pdf", width = 5)
ggplot(Osa.act.leaf, aes(x = Condition, y = MeanCp, group=Condition, color = Condition)) +
  geom_beeswarm(size=4, cex =4, color = "black", alpha = 0)+
  scale_y_continuous(breaks = -1:10, limits = c(-0.5, 5))+
  scale_color_manual(values = c("gray",green[c(30,60,100)],orange[c(40)]))+ 
  theme_bw()+
  ylab("Ct")+
  geom_spoke(data=mean.leaf,aes(x=as.numeric(Condition)-0.5, y=mean),angle = 0, radius = 1, size = 2)+
  theme(legend.position = "NA")+
  geom_beeswarm(size=4, cex =4, color = "black")
dev.off()

pdf("Fig2d.pdf", width = 3)
ggplot(Osa.act.root, aes(x = Condition, y = MeanCp, group=Condition, color = Condition)) +
  geom_beeswarm(size=4, cex =4, color = "black", alpha = 0)+
  scale_color_manual(values = c("gray",green[100]))+ 
  scale_y_continuous(breaks = -1:10, limits = c(-1.5,6))+
  theme_bw()+
  ylab("Ct")+
  geom_spoke(data=mean.root,aes(x=as.numeric(Condition)-0.5, y=mean),angle = 0, radius = 1, size = 2)+
  theme(legend.position = "NA")+
  geom_beeswarm(size=4, cex =4, color = "black")
dev.off()

Tae.act = read.table("qPCR/Tae_ACT_Cp.txt", header = T, sep = "\t", skip = 1)
Tae.act$Name = factor(Tae.act$Name, levels = c("NC","DTT 100 mM"))
tmp = data.frame(DiffCt = -(Tae.act$Cp-mean(Tae.act$Cp[Tae.act$Name=="NC"])))
Tae.act= cbind(Tae.act, tmp)
mean = Tae.act %>% group_by(Name) %>% summarise(mean = mean(DiffCt))
pdf("Fig2e.pdf", width = 3)
ggplot(Tae.act, aes(x = Name, y = DiffCt, group=Name, color = Name)) +
  geom_beeswarm(size=4, cex =4, alpha =0)+
  scale_color_manual(values = c("gray",green[c(30,60,100)],orange[c(40)]))+ 
  theme_bw() +
  ylab("Ct")+
  scale_y_continuous(breaks = -2:10, limits = c(-1,3))+
  geom_spoke(data=mean,aes(x=as.numeric(Name)-0.5, y=mean),angle = 0, radius = 1, size = 2)+
  theme(legend.position = "NA")+
  geom_beeswarm(size=4, cex =4, color = "black")
dev.off()

Osa.act = read.table("qPCR/Osa_act_qPCR.txt", header = T, sep = "\t")
data = Osa.act[Osa.act$Condition=="Std",]
data = cbind(data, data.frame(Mesured = 1/2^data$MeanCp))
data = cbind(data,data.frame(rep = as.character(1:3)))
lm = lm(Mesured ~ X, data)
pred = data.frame( X = range(data$X),  Mesured= predict(lm,newdata = data.frame(X=range(data$X))))

pdf("Fig2g.pdf", width = 14, height = 5)
ggplot(data, aes(x = X, y = Mesured, color = rep)) +
  geom_point(dat = data,aes(shape=rep), size = 8)+
  geom_line(data = pred , aes(x = X, y = Mesured), color = "red")+
  scale_y_continuous(breaks = c(7.5*10^-8, 5*10^-8, 2.5*10^-8))+
  theme_bw() + 
  scale_shape_discrete(solid=F)+
  theme(legend.position = "bottom")
dev.off()
cor(data$Mesured, predict(lm))^2

#### Fig. 3
rawcnt = NULL
files = dir("RSEM_out_lasy-Dlasy/")
for (file in files[-(1:2)]) {
  tmp = read.table(sprintf("RSEM_out_lasy-Dlasy/%s",file), header = T, stringsAsFactors = F)
  rawcnt = cbind(rawcnt, tmp[,5])
}
colnames(rawcnt) = as.numeric(unlist(lapply(strsplit(files[-(1:2)],"_"),function(x){return(gsub("f0","",x[3]))})))
rownames(rawcnt) = tmp[,2]
rawcnt=rawcnt[,colSums(rawcnt)>50000]
at = read.table("lasy-Dlasy_sample_attribute.txt", header = F, stringsAsFactors = T, fill = T)
at = at[colnames(rawcnt),]
load("20180511_Araport11_genes.201606.transcript.rep_ERCC_Virus7457.desc")
data = rownames(des)[des$NormalizationGroup=="data"]
rRNA = rownames(des)[des$NormalizationGroup=="rRNA"]

rpm = t(t(rawcnt)/colSums(rawcnt[data,])*10^6)
log2rpm = log2(rpm+1)

des = des[is.element(rownames(des), rownames(log2rpm)),]
source("scripts/colname2rgb.R")
library(dplyr)
library(reshape2)
lasy.sample = as.character(at$V1[at$V2=="Lasy-Seq"])
Dlasy.sample = as.character(at$V1[at$V2=="Dlasy-Seq"])
purification = c(rep("Purified",length(lasy.sample)),rep("Lysate", length(Dlasy.sample)))
cols = c(rep("gray",length(lasy.sample)),rep("green", length(Dlasy.sample)))
total = colSums(rawcnt) 
rRNA = colSums(rawcnt[rownames(des)[des$NormalizationGroup=="rRNA"],])

rRNA.ratio = rRNA/100000*100
rRNA.ratio.df = data.frame(rRNAratio = rRNA.ratio[c(lasy.sample, Dlasy.sample)], sample = factor(purification, levels=c("Purified","Lysate")))
rRNA.ratio.df.mean = melt(rRNA.ratio.df %>%
  group_by(sample) %>%
  summarise(mean(rRNAratio)))
g1 = ggplot(data = rRNA.ratio.df, aes(y = rRNAratio, x = sample, fill = sample)) +
  geom_bar(data = rRNA.ratio.df.mean,stat = "identity", aes(x = sample, y = value)) +
  theme_bw()+
  geom_beeswarm()+
  theme(legend.position = "none")+
  ylim(c(0,5))
Num.of.detected.gene = colSums(rawcnt>0)
Num.of.detected.gene.df = data.frame(Num = Num.of.detected.gene[c(lasy.sample, Dlasy.sample)], sample = factor(purification, levels=c("Purified","Lysate")))
Num.of.detected.gene.df.mean = melt(Num.of.detected.gene.df %>%
  group_by(sample) %>%
  summarise(mean(Num)))
g2 = ggplot(data = Num.of.detected.gene.df, aes(y = Num, x = sample, fill = sample)) +
  geom_bar(data = Num.of.detected.gene.df.mean,stat = "identity", aes(x=sample, y = value)) +
  theme_bw()+
  geom_quasirandom()+
  theme(legend.position = "none")

mapping.ratio.df = data.frame("mapping.ratio" = total/100000*100, sample = factor(purification, levels=c("Purified","Lysate"))) 
mapping.ratio.df.mean = melt(mapping.ratio.df %>%
  group_by(sample) %>%
  summarise(mean(mapping.ratio )))
g3 = ggplot(data = mapping.ratio.df, aes(y = mapping.ratio, x = sample, fill = sample))+
  geom_bar(data = mapping.ratio.df.mean,stat = "identity", aes(x=sample, y = value)) +
  theme_bw()+
  theme(legend.position = "none")+
  ylim(c(0,100))+
  geom_quasirandom()


library(gridExtra)
pdf("Fig3abc.pdf", width = 14, height = 3.5)
grid.arrange(g3,g1,g2,ncol = 3)
dev.off()
t.test(total[at[,2]=="Dlasy-Seq"],total[at[,2]=="Lasy-Seq"])
t.test(Num.of.detected.gene[at[,2]=="Dlasy-Seq"],Num.of.detected.gene[at[,2]=="Lasy-Seq"])
t.test(rRNA.ratio[at[,2]=="Dlasy-Seq"],rRNA.ratio[at[,2]=="Lasy-Seq"])

###DEG analysis
library(TCC)                          
param_FDR <- 0.05                     
data <- cbind(rawcnt[,at$V2=="Lasy-Seq"],rawcnt[,at$V2=="Dlasy-Seq"])
data.cl <- c(rep(1, sum(at$V2=="Lasy-Seq")), rep(2, sum(at$V2=="Dlasy-Seq")))
tcc <- new("TCC", data, data.cl)       

tcc <- calcNormFactors(tcc, norm.method="tmm", test.method="edger",
                       iteration=3, FDR=0.1, floorPDEG=0.05)

tcc <- estimateDE(tcc, test.method="edger", FDR=param_FDR)
result <- getResult(tcc, sort=FALSE)
rownames(result) = result$gene_id
sum(tcc$stat$q.value < param_FDR)   

DEG  =  tcc$gene_id[tcc$stat$q.value < param_FDR]
Dlasy.enriched.gnl = result[DEG,]
Lasy.enriched.gnl = as.character(Dlasy.enriched.gnl$gene_id[Dlasy.enriched.gnl$m.value<0])
Dlasy.enriched.gnl = as.character(Dlasy.enriched.gnl$gene_id[Dlasy.enriched.gnl$m.value>0])
fai = read.table("Araport11_genes.201606.transcript.rep_DELTA.fa.fai", stringsAsFactors = F)
rownames(fai) = fai$V1
length.df = data.frame(length = fai[,2], flag = "others", stringsAsFactors = F)
rownames(length.df) = fai[,1]
length.df[Dlasy.enriched.gnl,2] = "High in Dlasy-Seq (FDR<0.05)"
length.df[Lasy.enriched.gnl,2] = "Low in Dlasy-Seq (FDR<0.05)"

#log2rpm plot
rawcnt.marge = NULL
for (m in unique(at[,2])) {
  cat(m)
  cat("\n")
  tmp = apply(rawcnt[,as.character(at[,1])[at[,2]==m]], 1, sum)
  rawcnt.marge = cbind(rawcnt.marge,tmp)
}
colnames(rawcnt.marge) = unique(at[,2])
rownames(rawcnt.marge) = rownames(rawcnt)
dim(rawcnt.marge)


DEG.nonhost = DEG[grep("G",DEG, invert = T)]
DEG.high = DEG[result[DEG,]$m.value>0]
DEG.low = DEG[result[DEG,]$m.value<0]
log2rpm = t(log2((t(rawcnt.marge)/colSums(rawcnt.marge)*10^6)+1))
pdf("Fig3d.pdf")
plot(log2rpm[,"Lasy-Seq"],log2rpm[,"Dlasy-Seq"], main = sprintf("Pearson's correlation: %s",round(cor(log2rpm[rownames(des)[des$NormalizationGroup=="data"],"Lasy-Seq"],log2rpm[rownames(des)[des$NormalizationGroup=="data"],"Dlasy-Seq"]),3)), pch = 16, col = "gray", xlab = "Relative expression level (log2rpm)\ncDNA from Purified RNA", ylab = "Relative expression level (log2rpm)\ncDNA from lysate")
points(log2rpm[DEG.high,"Lasy-Seq"],log2rpm[DEG.high,"Dlasy-Seq"], col = "red", pch = 16)
points(log2rpm[DEG.low,"Lasy-Seq"],log2rpm[DEG.low,"Dlasy-Seq"], col = "blue", pch = 16)
lines(c(-100,100),c(-100,100), col = "black")
dev.off()
length.df$flag = factor(length.df$flag, levels = c("others", "High in Dlasy-Seq (FDR<0.05)", "Low in Dlasy-Seq (FDR<0.05)")) 
pdf("Fig3e.pdf")
g1 = ggplot(length.df[length.df$flag=="High in Dlasy-Seq (FDR<0.05)",], aes(x=length)) +
  geom_histogram(aes(y = ..ncount../sum(..ncount..)), position = "identity",  bins = 100, fill = "red") +
  xlim(c(0,6000))+
  ylim(c(0,0.15))+
  theme_bw() + 
  theme(legend.position = "bottom")
g2 = ggplot(length.df[length.df$flag=="Low in Dlasy-Seq (FDR<0.05)",], aes(x=length)) +
  geom_histogram(aes(y = ..ncount../sum(..ncount..)), position = "identity",  bins = 100, fill = "blue") +
  xlim(c(0,6000))+
  ylim(c(0,0.15))+
  theme_bw() + 
  theme(legend.position = "bottom")
g3 = ggplot(length.df[length.df$flag=="others",], aes(x=length)) +
  geom_histogram(aes(y = ..ncount../sum(..ncount..)), position = "identity",  bins = 100, fill = "gray") +
  xlim(c(0,6000))+
  ylim(c(0,0.15))+
  theme_bw() + 
  theme(legend.position = "bottom")
grid.arrange(g1,g2,g3, ncol = 1)
dev.off()

mean(length.df$length[length.df$flag=="others"]) # 1352.166
mean(length.df$length[length.df$flag== "High in Dlasy-Seq (FDR<0.05)"]) # 798.7264
# Fisher's test
t =1000
mat = matrix(c(nrow(fai), sum(fai$V2<t), length(Dlasy.enriched.gnl), sum(fai[Dlasy.enriched.gnl,]$V2<t)), byrow = T, nco = 2)
fisher.test(mat)

# Check the relationship between length and expression level
pdf("Osa_length_ex_plot.pdf")
plot(length.df[rownames(log2rpm),]$length, log2rpm[,2], ylim = c(0,15), pch = 16, col = "lightgray")
points(length.df[DEG.high,]$length, log2rpm[DEG.high,2], col = "blue")
points(length.df[DEG.low,]$length, log2rpm[DEG.low,2], col = "red")
dev.off()

#### Fig.S2
# RNA-Seq of yeast from lysate and purified
files = list.files("yeast/")
rawcnt.yeast = NULL
for (i in 1:length(files)) {
  tmp = read.table(sprintf("yeast/%s", files[i]), header = T)
  rawcnt.yeast = cbind(rawcnt.yeast, tmp$expected_count)
}
rownames(rawcnt.yeast) = tmp$gene_id
plot(colSums(rawcnt.yeast))
rawcnt.yeast = rawcnt.yeast[,colSums(rawcnt.yeast)>10000]
colnames(rawcnt.yeast) = c(sprintf("Y.lysate-%s",1:6),sprintf("Y.purified-%s",1:6))
library(TCC)
param_FDR <- 0.05
data <- rawcnt.yeast[,1:12]
data.cl <- rep(c(1,2),each = 6 )
tcc <- new("TCC", data, data.cl)
tcc <- calcNormFactors(tcc, norm.method="tmm", test.method="edger", iteration=3, FDR=0.1, floorPDEG=0.05)
tcc <- estimateDE(tcc, test.method="edger", FDR=param_FDR)
result <- getResult(tcc, sort=FALSE)
rownames(result) = result$gene_id
result = cbind(TCC::getNormalizedData(tcc),result[,-1])
sum(tcc$stat$q.value < param_FDR)      #条件を満たす遺伝子数を表示

DEG  =  tcc$gene_id[tcc$stat$q.value < param_FDR]
DEG.high = DEG[result[DEG,]$m.value<0]
DEG.low = DEG[result[DEG,]$m.value>0]

fai = read.table("Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.fai", stringsAsFactors = F)
rownames(fai) = fai$V1
length.df = data.frame(length = fai[,2], flag = "others", stringsAsFactors = F)
rownames(length.df) = fai[,1]
length.df[DEG.high,2] = "High in Dlasy-Seq (FDR<0.05)"
length.df[DEG.low,2] = "Low in Dlasy-Seq (FDR<0.05)"
length.df$flag = factor(length.df$flag, levels = c("others", "High in Dlasy-Seq (FDR<0.05)", "Low in Dlasy-Seq (FDR<0.05)")) 

#log2rpm plot
rawcnt.marge = NULL
tmp = apply(rawcnt.yeast[,1:6], 1, sum)
rawcnt.marge = cbind(rawcnt.marge,tmp)
tmp = apply(rawcnt.yeast[,7:12], 1, sum)
rawcnt.marge = cbind(rawcnt.marge,tmp)


log2rpm = t(log2((t(rawcnt.marge)/colSums(rawcnt.marge)*10^6)+1))
pdf("FigS2a.pdf")
plot(log2rpm[,2],log2rpm[,1], main = sprintf("Pearson's correlation: %s",cor(log2rpm[,1], log2rpm[,2])), ylab = "Dlasy", xlab = "Lasy" , pch = 16, ylim = c(0,15), xlim = c(0,15), col = "gray")
points(log2rpm[DEG.high,2],log2rpm[DEG.high,1], col = "red", pch = 16)
points(log2rpm[DEG.low,2],log2rpm[DEG.low,1], col = "blue", pch = 16)
lines(c(-100,100),c(-100,100), col = "black")
dev.off()

pdf("FigS2c.pdf")
g1 = ggplot(length.df[length.df$flag=="High in Dlasy-Seq (FDR<0.05)",], aes(x=length)) +
  geom_histogram(aes(y = ..ncount../sum(..ncount..)), position = "identity",  bins = 100, fill = "red") +
  xlim(c(0,8000))+
  ylim(c(0,0.2))+
  theme_bw() + 
  theme(legend.position = "bottom")
g2 = ggplot(length.df[length.df$flag=="Low in Dlasy-Seq (FDR<0.05)",], aes(x=length)) +
  geom_histogram(aes(y = ..ncount../sum(..ncount..)), position = "identity",  bins = 100, fill = "blue") +
  xlim(c(0,8000))+
  ylim(c(0,0.2))+
  theme_bw() + 
  theme(legend.position = "bottom")
g3 = ggplot(length.df[length.df$flag=="others",], aes(x=length)) +
  geom_histogram(aes(y = ..ncount../sum(..ncount..)), position = "identity",  bins = 100, fill = "gray") +
  xlim(c(0,8000))+
  ylim(c(0,0.2))+
  theme_bw() + 
  theme(legend.position = "bottom")
grid.arrange(g1,g2,g3, ncol = 1)
dev.off()
mean(length.df$length[length.df$flag=="others"]) # 1537.608
mean(length.df$length[length.df$flag!="others"]) # 798.7264

# RNA-Seq of zebrafish from lysate and purified
files = list.files("zebrafish/")
rawcnt.zebra = NULL
for (i in 1:length(files)) {
  tmp = read.table(sprintf("zebrafish/%s", files[i]), header = T)
  rawcnt.zebra = cbind(rawcnt.zebra, tmp$expected_count)
}
rownames(rawcnt.zebra) = tmp$gene_id
colnames(rawcnt.zebra) = c(sprintf("Z.lysate.%s",1:6),sprintf("Z.purified.%s",1:6))
library(TCC)
param_FDR <- 0.05
data <- rawcnt.zebra
data.cl <- rep(c(1,2),each = 6 )
tcc <- new("TCC", data, data.cl)  
tcc <- calcNormFactors(tcc, norm.method="tmm", test.method="edger", iteration=3, FDR=0.1, floorPDEG=0.05)
tcc <- estimateDE(tcc, test.method="edger", FDR=param_FDR)
result <- getResult(tcc, sort=FALSE) 
rownames(result) = result$gene_id
sum(tcc$stat$q.value < param_FDR) 
DEG  =  tcc$gene_id[tcc$stat$q.value < param_FDR]
DEG.high = DEG[result[DEG,]$m.value<0]
DEG.low = DEG[result[DEG,]$m.value>0]

fai = read.table("Danio_rerio.GRCz11.cdna.all.fa.fai", stringsAsFactors = F)
rownames(fai) = fai$V1
length.df = data.frame(length = fai[,2], flag = "others", stringsAsFactors = F)
rownames(length.df) = fai[,1]
length.df[DEG.high,2] = "High in Dlasy-Seq (FDR<0.05)"
length.df[DEG.low,2] = "Low in Dlasy-Seq (FDR<0.05)"

#log2rpm plot
rawcnt.marge = NULL
tmp = apply(rawcnt.zebra[,1:6], 1, sum)
rawcnt.marge = cbind(rawcnt.marge,tmp)
tmp = apply(rawcnt.zebra[,7:12], 1, sum)
rawcnt.marge = cbind(rawcnt.marge,tmp)

log2rpm.zebra = t(log2((t(rawcnt.zebra)/colSums(rawcnt.zebra)*10^6)+1))

log2rpm = t(log2((t(rawcnt.marge)/colSums(rawcnt.marge)*10^6)+1))
pdf("FigS2b.pdf")
plot(log2rpm[,2],log2rpm[,1], main = sprintf("Pearson's correlation: %s",cor(log2rpm[,1], log2rpm[,2])), ylab = "Dlasy", xlab = "Lasy" , pch = 16, ylim = c(0,15), xlim = c(0,15), col = "gray")
points(log2rpm[DEG.high,2],log2rpm[DEG.high,1], col = "red", pch = 16)
points(log2rpm[DEG.low,2],log2rpm[DEG.low,1], col = "blue", pch = 16)
lines(c(-100,100),c(-100,100), col = "black")
dev.off()

length.df$flag = factor(length.df$flag, levels = c("others", "High in Dlasy-Seq (FDR<0.05)", "Low in Dlasy-Seq (FDR<0.05)")) 
pdf("FigS2d.pdf")
g1 = ggplot(length.df[length.df$flag=="High in Dlasy-Seq (FDR<0.05)",], aes(x=length)) +
  geom_histogram(aes(y = ..ncount../sum(..ncount..)), position = "identity",  bins = 100, fill = "red") +
  xlim(c(0,20000))+
  ylim(c(0,0.1))+
  theme_bw() + 
  theme(legend.position = "bottom")
g2 = ggplot(length.df[length.df$flag=="Low in Dlasy-Seq (FDR<0.05)",], aes(x=length)) +
  geom_histogram(aes(y = ..ncount../sum(..ncount..)), position = "identity",  bins = 100, fill = "blue") +
  xlim(c(0,20000))+
  ylim(c(0,0.1))+
  theme_bw() + 
  theme(legend.position = "bottom")
g3 = ggplot(length.df[length.df$flag=="others",], aes(x=length)) +
  geom_histogram(aes(y = ..ncount../sum(..ncount..)), position = "identity",  bins = 100, fill = "gray") +
  xlim(c(0,20000))+
  ylim(c(0,0.1))+
  theme_bw() + 
  theme(legend.position = "bottom")
grid.arrange(g1,g2,g3, ncol = 1)
dev.off()

mean(length.df$length[length.df$flag=="others"]) # 1537.608
mean(length.df$length[length.df$flag!="others"]) # 798.7264

load("ensemble_zebra_transcript2gene")
DEG.high.gene = zebra_transcript2gene[DEG.high]
maternal.genes = read.csv("DEV095091TableS4.csv", header = T, stringsAsFactors = F, skip = 1)
maternal.zygotic.genes = read.csv("DEV095091TableS5.csv", header = T, stringsAsFactors = F, skip = 1)
zygotic.genes = read.csv("DEV095091TableS6.csv", header = T, stringsAsFactors = F, skip = 1)
maternal.genes = unique(maternal.genes$Gene.ID)
maternal.zygotic.genes = unique(maternal.zygotic.genes$Gene.ID)
zygotic.genes = unique(zygotic.genes$Gene.ID)
DEG.high.gene = unique(DEG.high.gene)
sum(is.element(DEG.high.gene, maternal.genes))
sum(is.element(DEG.high.gene, maternal.zygotic.genes))
fai.gene = unique(zebra_transcript2gene[fai$V1])
protein.coding = c(maternal.genes, maternal.zygotic.genes, zygotic.genes)
fai.gene = fai.gene[is.element(fai.gene, protein.coding)]
DEG.high.gene = DEG.high.gene[is.element(DEG.high.gene, protein.coding)]

category = c(zygotic.genes)
category = maternal.zygotic.genes
category = c(maternal.zygotic.genes, maternal.genes)
category = c(maternal.genes)

DEG.high.maternal = sum(is.element(DEG.high.gene, category))
DEG.high.zygotic = sum(!is.element(DEG.high.gene, category ))
fai.maternal = sum(is.element(fai.gene, category))
fai.zygotic = sum(!is.element(fai.gene, category))
fisher.test(matrix(c(DEG.high.maternal, fai.maternal, DEG.high.zygotic, fai.zygotic), ncol = 2, byrow = T))

#### Fig.4a
# Load target
target = read.table("target_list.txt", header = F, stringsAsFactors = F)[,1]
target = matrix(unlist(strsplit(target,"\\.")), byrow = T, ncol = 2)[,1]
target = unique(target)
target.short = target

id = read.table("Araport11_genes.201606.transcript_ERCC_Virus3981.id", stringsAsFactors = F)[,1]
tmp = data.frame("gene.id" = id[grep("^AT",id)], "short.gene.id" = matrix(unlist(strsplit(id[grep("^AT",id)],"\\.")), byrow = T, ncol = 2)[,1], stringsAsFactors = F)
tmp = tmp[is.element(tmp$short.gene.id,target),]
target = tmp$gene.id

# load rawcnt
files = list.files("RSEM_out_UMI/")
files = files[grep("genes", files)]
files = files[grep("Index", files)]
DeLTa.files = files[grep("Index1", files)]
DeLTa.total = sort(as.integer(unique(lapply(strsplit(DeLTa.files,"_"), FUN = function(x){return(x[2])}))))
DeLTa.rawcnt = list()
for (t in DeLTa.total){
  rawcnt = NULL
  for (r in 1:96) {
    file = files[grep(sprintf("Index1%02d_%s_",r,t), files)]
    file = sprintf("RSEM_out_UMI/%s",file)
    tmp = read.table(file, header = T, stringsAsFactors = F, row.names = 1)
    tmp = tmp[target,]$expected_count/sum(tmp[id[grep("^AT",id)],]$expected_count)*10^6
    tmp2 = NULL
    for(gn in target.short){
      tmp2 = c(tmp2,log2(sum(tmp[grep(gn, target)])+1))
    }
    tmp = tmp2
    rawcnt = cbind(rawcnt, tmp)
  }
  rawcnt = list(rawcnt)
  names(rawcnt) = t
  DeLTa.rawcnt = c(DeLTa.rawcnt, rawcnt)
}

# Load full data
DeLTa.all = read.table("RSEM_out_UMI/DeLTa.genes.results", header = T, row.names = 1)
DeLTa.all$expected_count = DeLTa.all$expected_count/sum(DeLTa.all[id[grep("^AT",id)],]$expected_count)*10^6
tmp2 = NULL
tmp = DeLTa.all[target,]$expected_count
for(gn in target.short){
  tmp2 = c(tmp2,log2(sum(tmp[grep(gn, target)])+1))
}
names(tmp2) = target.short
DeLTa.all = tmp2

# Lasy.all = read.table("RSEM_out_UMI/Lasy.genes.results", header = T, row.names = 1)
# Lasy.all$expected_count = Lasy.all$expected_count/sum(Lasy.all[id[grep("^AT",id)],]$expected_count)*10^6
# tmp2 = NULL
# tmp = Lasy.all[target,]$expected_count
# for(gn in target.short){
#   tmp2 = c(tmp2,log2(sum(tmp[grep(gn, target)])+1))
# }
# names(tmp2) = target.short
# Lasy.all = tmp2

DeLTa.cors.all = NULL
for(i in 1:length(DeLTa.rawcnt)){
  DeLTa.cors.all = cbind(DeLTa.cors.all, apply(DeLTa.rawcnt[[i]], 2, FUN = function(x){
    return(cor(x,DeLTa.all))
  })
  )
}
colnames(DeLTa.cors.all) = DeLTa.total

# Lasy.cors.all = NULL
# for(i in 1:length(Lasy.rawcnt)){
#   Lasy.cors.all = cbind(Lasy.cors.all, apply(Lasy.rawcnt[[i]], 2, FUN = function(x){
#     return(cor(x,Lasy.all))
#   })
#   )
# }
# colnames(Lasy.cors.all) = Lasy.total


# plot
library(reshape2)
library(ggplot2)
DeLTa.cors.all = melt(DeLTa.cors.all)
DeLTa.cors.all = data.frame("Method" = "DeLTa-Seq", "total" = DeLTa.cors.all[,2], "Cor" = DeLTa.cors.all[,3])
data = DeLTa.cors.all
data.ave = data %>%
  group_by(total) %>%
  summarise_all(mean)
data.ave= cbind(data.ave, data.frame("x" = "1"))

pdf("Fig4a.pdf", width = 4)
ggplot(data, aes(x=as.factor(total),y=Cor, group = total))+
  geom_boxplot(fill = "red")+
  ylim(c(0.75,1))+
  theme_bw()
dev.off()

##### Fig. 4bc
# Load target
target = read.table("target_list.txt", header = F, stringsAsFactors = F)[,1]
target = matrix(unlist(strsplit(target,"\\.")), byrow = T, ncol = 2)[,1]
target = unique(target)
target.short = target

id = read.table("Araport11_genes.201606.transcript_ERCC_Virus3981.id", stringsAsFactors = F)[,1]
tmp = data.frame("gene.id" = id[grep("^AT",id)], "short.gene.id" = matrix(unlist(strsplit(id[grep("^AT",id)],"\\.")), byrow = T, ncol = 2)[,1], stringsAsFactors = F)
tmp = tmp[is.element(tmp$short.gene.id,target),]
target = tmp$gene.id

# load rawcnt
files = list.files("RSEM_out_UMI/")
files = files[grep("genes", files)]
DeLTa.files = files[grep("DeLTa", files)]
Lasy.files = files[grep("Lasy", files)]
DeLTa.total = sort(as.integer(unique(lapply(strsplit(DeLTa.files,"_"), FUN = function(x){return(x[3])}))))
Lasy.total = sort(as.integer(unique(lapply(strsplit(Lasy.files,"_"), FUN = function(x){return(x[3])}))))

Lasy.rawcnt = list()
for (t in Lasy.total){
  rawcnt = NULL
  for (r in 1:12) {
    file = files[grep(sprintf("Lasy_%02d_%s_",r,t), files)]
    file = sprintf("RSEM_out_UMI//%s",file)
    tmp = read.table(file, header = T, stringsAsFactors = F, row.names = 1)
    tmp = tmp[target,]$expected_count/sum(tmp[id[grep("^AT",id)],]$expected_count)*10^6
    tmp2 = NULL
    for(gn in target.short){
      tmp2 = c(tmp2,log2(sum(tmp[grep(gn, target)])+1))
    }
    tmp = tmp2
    names(tmp) = target.short
    rawcnt = cbind(rawcnt, tmp)
  }
  rawcnt = list(rawcnt)
  names(rawcnt) = t
  Lasy.rawcnt = c(Lasy.rawcnt, rawcnt)
}

DeLTa.rawcnt = list()
for (t in DeLTa.total){
  rawcnt = NULL
  for (r in 1:12) {
    file = files[grep(sprintf("DeLTa_%02d_%s_",r,t), files)]
    file = sprintf("RSEM_out_UMI/%s",file)
    tmp = read.table(file, header = T, stringsAsFactors = F, row.names = 1)
    tmp = tmp[target,]$expected_count/sum(tmp[id[grep("^AT",id)],]$expected_count)*10^6
    tmp2 = NULL
    for(gn in target.short){
      tmp2 = c(tmp2,log2(sum(tmp[grep(gn, target)])+1))
    }
    tmp = tmp2
    rawcnt = cbind(rawcnt, tmp)
  }
  rawcnt = list(rawcnt)
  names(rawcnt) = t
  DeLTa.rawcnt = c(DeLTa.rawcnt, rawcnt)
}


# Against marged rawcnt
DeLTa.all = read.table("RSEM_out_UMI/DeLTa.genes.results", header = T, row.names = 1)
DeLTa.all$expected_count = DeLTa.all$expected_count/sum(DeLTa.all[id[grep("^AT",id)],]$expected_count)*10^6
tmp2 = NULL
tmp = DeLTa.all[target,]$expected_count
for(gn in target.short){
  tmp2 = c(tmp2,log2(sum(tmp[grep(gn, target)])+1))
}
names(tmp2) = target.short
DeLTa.all = tmp2

Lasy.all = read.table("RSEM_out_UMI/Lasy.genes.results", header = T, row.names = 1)
Lasy.all$expected_count = Lasy.all$expected_count/sum(Lasy.all[id[grep("^AT",id)],]$expected_count)*10^6
tmp2 = NULL
tmp = Lasy.all[target,]$expected_count
for(gn in target.short){
  tmp2 = c(tmp2,log2(sum(tmp[grep(gn, target)])+1))
}
names(tmp2) = target.short
Lasy.all = tmp2

DeLTa.cors.all = NULL
for(i in 1:length(DeLTa.rawcnt)){
  DeLTa.cors.all = cbind(DeLTa.cors.all, apply(DeLTa.rawcnt[[i]], 2, FUN = function(x){
    return(cor(x,DeLTa.all))
  })
  )
}
colnames(DeLTa.cors.all) = DeLTa.total

Lasy.cors.all = NULL
for(i in 1:length(Lasy.rawcnt)){
  Lasy.cors.all = cbind(Lasy.cors.all, apply(Lasy.rawcnt[[i]], 2, FUN = function(x){
    return(cor(x,Lasy.all))
  })
  )
}
colnames(Lasy.cors.all) = Lasy.total


# plot
library(reshape2)
Lasy.cors.all = melt(Lasy.cors.all)
Lasy.cors.all = data.frame("Method" = "Lasy-Seq", "total" = Lasy.cors.all[,2], "Cor" = Lasy.cors.all[,3])

DeLTa.cors.all = melt(DeLTa.cors.all)
DeLTa.cors.all = data.frame("Method" = "DeLTa-Seq", "total" = DeLTa.cors.all[,2], "Cor" = DeLTa.cors.all[,3])

DeLTa.total.decoy = Lasy.total[!is.element(Lasy.total, DeLTa.total)]
DeLTa.cors.all.decoy = data.frame("Method" = "DeLTa-Seq", "total" = DeLTa.total.decoy, "Cor" = 0.5)
Lasy.total.decoy = DeLTa.total[!is.element(DeLTa.total, Lasy.total)]
Lasy.cors.all.decoy = data.frame("Method" = "Lasy-Seq", "total" = Lasy.total.decoy, "Cor" = 0.5)

data = rbind(Lasy.cors.all, DeLTa.cors.all, DeLTa.cors.all.decoy, Lasy.cors.all.decoy)
library(dplyr)
data.mean = data %>% 
  group_by(interaction(Method, total)) %>%
  summarise(mean =mean(Cor))


library(ggplot2)
pdf("Fig4c.pdf", width = 14)
ggplot(data, aes(x=as.factor(total),y=Cor, group = interaction(Method,total), fill = Method))+
  geom_boxplot()+
  # ylim(c(0.6,1))+
  # scale_color_manual(values = c("gray","red"))+
  scale_fill_manual(values = c("gray","red"))+
  theme_bw()
dev.off()

# plot each gene
library(dplyr)
library(ggbeeswarm)
pdf("Fig4b.pdf",width = 14)
for(g in 1:length(target.short)){
  data = NULL
  gn = target.short[g]
  cat(gn)
  cat("\n")
  for(i in 1:length(DeLTa.rawcnt)){
    DeLTa = DeLTa.rawcnt[[i]][g,]
    data = rbind(data,data.frame("Method" = "DeLTa", "log2rpm" = DeLTa ,"total" = names(DeLTa.rawcnt)[i]))
  }
  for(i in 1:length(Lasy.rawcnt)){
    Lasy = Lasy.rawcnt[[i]][g,]
    data = rbind(data,data.frame("Method" = "Lasy", "log2rpm"= Lasy, "total" = names(Lasy.rawcnt)[i]))
  }
  data$Method = factor(data$Method, levels = c( "Lasy","DeLTa"))
  data.each = data
  data = data %>%
    group_by(interaction(Method, total)) %>%
    summarize(mean = mean(log2rpm), sd = sd(log2rpm))
  tmp = matrix(unlist(strsplit(as.character(data$`interaction(Method, total)`),"\\.")), byrow = T, ncol = 2)
  data = data.frame("Method" = tmp[,1], "Total" = tmp[,2], "Mean" = data$mean, "SD" = data$sd, stringsAsFactors = F)
  level = sort(unique(as.integer(data$Total)))
  data$Total = factor(data$Total, levels = level)
  data.each$total = factor(as.character(data.each$total), levels = level)
  data$Method = factor(data$Method, levels = c("Lasy","DeLTa"))
  gg =ggplot(data, aes(x = Total, y = Mean, group = Method, fill = Method))+
    # geom_beeswarm(data = data.each, aes(x=total, y = log2rpm, color = Method), size = 0.5)+
    geom_line()+
    geom_violin(data = data.each[data.each$Method=="Lasy",], aes(x=total, y = log2rpm, group = total))+
    geom_violin(data = data.each[data.each$Method=="DeLTa",], aes(x=total, y = log2rpm, group = total))+
    geom_line()+
    ggtitle(gn)+
    scale_fill_manual(values = c("gray","red"))+
    ylim(c(0,20))+
    theme_bw()
  
  print(gg)
}
dev.off()

#### Fig. S3
#UMI conversion ratio
umi = read.table("UMI_conversion.txt", header = T)
library(ggplot2)
pdf("FigS3.pdf")
ggplot(umi, aes(x = 1:192, y = uniq/all))+
  geom_point(color = rep(c("gray","red"), each = 96))+
  ylim(c(0,1))+
  theme_bw()
dev.off()

mean((umi$uniq/umi$all)[1:96])
sd((umi$uniq/umi$all)[1:96])
mean((umi$uniq/umi$all)[97:192])
sd((umi$uniq/umi$all)[97:192])

#### Fig.5
count = read.table("input_volume_read.count", row.names = 2)
lib.ok = rownames(count)[count[,1] > 10^4]
files.ok = sprintf("%s.fq.genes.results",unlist(lapply(strsplit(lib.ok,"\\."),function(x){return(x[1])})))
files = dir("RSEM_out_input_volume/")
files = files[is.element(files,files.ok)]
at = matrix(unlist(lapply(strsplit(files,"_"),FUN = function(x){
  return(c(x[1],x[2]))
})),byrow = T,ncol=2)
colnames(at) = c("input","rep")

libs = unique(sprintf("%s-%s",at[,1],at[,2]))
file = files[1]
tmp = read.table(sprintf("RSEM_out_input_volume/%s",file), stringsAsFactors = F, header = T)
rawcnt = matrix(0,ncol = length(libs), nrow = nrow(tmp))
colnames(rawcnt) = libs
rownames(rawcnt) = rownames(tmp)

for (file in files) {
  cat(file)
  cat("\n")
  v = strsplit(file,"_")
  lib = sprintf("%s-%s",v[[1]][1],v[[1]][2])
  tmp = read.table(sprintf("RSEM_out_input_volume/%s",file), stringsAsFactors = F, header = T)
  rawcnt[,lib] = rawcnt[,lib] + tmp$expected_count
}
rownames(rawcnt) = tmp$gene_id
at = data.frame("lib" = colnames(rawcnt), input = unlist(lapply(strsplit(colnames(rawcnt),"-"), FUN = function(x){return(x[1])})))
total.read = colSums(rawcnt)
log2rpm = t(log2(t(rawcnt)/colSums(rawcnt*10^6)+1))
tmp = table(at$input)
cols = c(rep("red",tmp[1]),rep("blue",tmp[2]),rep("green",tmp[3]),rep("orange",tmp[4]))
summary.df = data.frame("input" = factor(c(rep(0.01,tmp[1]),rep(0.1,tmp[2]),rep(1,tmp[3]),rep(10,tmp[4]))), "mapped_read"=total.read,  "target_ratio" = total.read/10000*100, "Num.of.detected.genes" = colSums(rawcnt>0))
pallete = c( "green", "orange", "blue","red")

library(dplyr)
summary.df %>%
  group_by(input) %>%
  summarise(mean = mean(target_ratio))
summary.df %>%
  group_by(input) %>%
  summarise(mean = mean(Num.of.detected.genes))


pdf("Fig5ab.pdf", height = 3.5)
g_base = ggplot(summary.df, aes(x=input))+
  theme_bw()+
  scale_fill_manual(values = pallete)+
  scale_color_manual(values = pallete)+
  xlab("Input RNA amount (ng)")

g = g_base +   geom_boxplot(aes(y = target_ratio,color=input)) + 
  ylim(c(0,100))+
  geom_jitter(aes(y =  target_ratio,color=input))
print(g)

g = g_base +   geom_boxplot(aes(y = Num.of.detected.genes,color=input))+
  ylim(c(0,100)) +
  geom_jitter(aes(y = Num.of.detected.genes,color=input))
print(g)

dev.off()

cors = NULL
inputs =unique(at[,2])
for(input in inputs[1:3]){
  for(i in 145:149){
    tmp = log2rpm[,at[,2]==input]
    tmp = data.frame("cor"=cor(tmp,log2rpm[,i]), "input" = input)
    cors = rbind(cors,tmp)
  }
}
for(i in 145:149){
  for (m in 146:149) {
    if(i<m){
      tmp = data.frame("cor"=cor(log2rpm[,m],log2rpm[,i]), "input" = inputs[4])
      cors = rbind(cors,tmp)
    }
  }
}

pdf("Fig5c.pdf", height = 3.5)
ggplot(cors, aes(x=input, color = input,y = cor))+
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values = pallete)+
  scale_color_manual(values = pallete)+
  ylab("Pearson's correlation")+
  xlab("Input RNA amount (ng)")+
  geom_jitter()

dev.off()

############# Fig.6 DeLTa-Seq VS Lasy-Seq (individual replicate 96)
# Load target
target = read.table("target_list.txt", header = F, stringsAsFactors = F)[,1]
target = matrix(unlist(strsplit(target,"\\.")), byrow = T, ncol = 2)[,1]
target = unique(target)
target.short = target

id = read.table("Araport11_genes.201606.transcript_ERCC_Virus3981.id", stringsAsFactors = F)[,1]
tmp = data.frame("gene.id" = id[grep("^AT",id)], "short.gene.id" = matrix(unlist(strsplit(id[grep("^AT",id)],"\\.")), byrow = T, ncol = 2)[,1], stringsAsFactors = F)
tmp = tmp[is.element(tmp$short.gene.id,target),]
target = tmp$gene.id

# load rawcnt
files = list.files("RSEM_out_UMI/")
files = files[grep("genes", files)]
files = files[grep("Index", files)]
DeLTa.files = files[grep("Index1", files)]
Lasy.files = files[grep("Index0", files)]
DeLTa.total = sort(as.integer(unique(lapply(strsplit(DeLTa.files,"_"), FUN = function(x){return(x[2])}))))
Lasy.total = sort(as.integer(unique(lapply(strsplit(Lasy.files,"_"), FUN = function(x){return(x[2])}))))

Lasy.rawcnt = list()
for (t in Lasy.total){
  rawcnt = NULL
  for (r in 1:96) {
    file = files[grep(sprintf("Index0%02d_%s_",r,t), files)]
    file = sprintf("RSEM_out_UMI/%s",file)
    tmp = read.table(file, header = T, stringsAsFactors = F, row.names = 1)
    tmp = tmp[target,]$expected_count/sum(tmp[id[grep("^AT",id)],]$expected_count)*10^6
    tmp2 = NULL
    for(gn in target.short){
      tmp2 = c(tmp2,log2(sum(tmp[grep(gn, target)])+1))
    }
    tmp = tmp2
    names(tmp) = target.short
    rawcnt = cbind(rawcnt, tmp)
  }
  rawcnt = list(rawcnt)
  names(rawcnt) = t
  Lasy.rawcnt = c(Lasy.rawcnt, rawcnt)
}

DeLTa.rawcnt = list()
for (t in DeLTa.total){
  rawcnt = NULL
  for (r in 1:96) {
    file = files[grep(sprintf("Index1%02d_%s_",r,t), files)]
    file = sprintf("RSEM_out_UMI/%s",file)
    tmp = read.table(file, header = T, stringsAsFactors = F, row.names = 1)
    tmp = tmp[target,]$expected_count/sum(tmp[id[grep("^AT",id)],]$expected_count)*10^6
    tmp2 = NULL
    for(gn in target.short){
      tmp2 = c(tmp2,log2(sum(tmp[grep(gn, target)])+1))
    }
    tmp = tmp2
    rawcnt = cbind(rawcnt, tmp)
  }
  rawcnt = list(rawcnt)
  names(rawcnt) = t
  DeLTa.rawcnt = c(DeLTa.rawcnt, rawcnt)
}

# Against marged rawcnt
DeLTa.all = read.table("RSEM_out_UMI/DeLTa.genes.results", header = T, row.names = 1)
DeLTa.all$expected_count = DeLTa.all$expected_count/sum(DeLTa.all[id[grep("^AT",id)],]$expected_count)*10^6
tmp2 = NULL
tmp = DeLTa.all[target,]$expected_count
for(gn in target.short){
  tmp2 = c(tmp2,log2(sum(tmp[grep(gn, target)])+1))
}
names(tmp2) = target.short
DeLTa.all = tmp2

Lasy.all = read.table("RSEM_out_UMI/Lasy.genes.results", header = T, row.names = 1)
Lasy.all$expected_count = Lasy.all$expected_count/sum(Lasy.all[id[grep("^AT",id)],]$expected_count)*10^6
tmp2 = NULL
tmp = Lasy.all[target,]$expected_count
for(gn in target.short){
  tmp2 = c(tmp2,log2(sum(tmp[grep(gn, target)])+1))
}
names(tmp2) = target.short
Lasy.all = tmp2

DeLTa.cors.all = NULL
for(i in 1:length(DeLTa.rawcnt)){
  DeLTa.cors.all = cbind(DeLTa.cors.all, apply(DeLTa.rawcnt[[i]], 2, FUN = function(x){
    return(cor(x,DeLTa.all))
  })
  )
}
colnames(DeLTa.cors.all) = DeLTa.total

Lasy.cors.all = NULL
for(i in 1:length(Lasy.rawcnt)){
  Lasy.cors.all = cbind(Lasy.cors.all, apply(Lasy.rawcnt[[i]], 2, FUN = function(x){
    return(cor(x,Lasy.all))
  })
  )
}
colnames(Lasy.cors.all) = Lasy.total


# plot
library(reshape2)
Lasy.cors.all = melt(Lasy.cors.all)
Lasy.cors.all = data.frame("Method" = "Lasy-Seq", "total" = Lasy.cors.all[,2], "Cor" = Lasy.cors.all[,3])

DeLTa.cors.all = melt(DeLTa.cors.all)
DeLTa.cors.all = data.frame("Method" = "DeLTa-Seq", "total" = DeLTa.cors.all[,2], "Cor" = DeLTa.cors.all[,3])


DeLTa.total.decoy = Lasy.total[!is.element(Lasy.total, DeLTa.total)]
DeLTa.cors.all.decoy = data.frame("Method" = "DeLTa-Seq", "total" = DeLTa.total.decoy, "Cor" = 0.4)
Lasy.total.decoy = DeLTa.total[!is.element(DeLTa.total, Lasy.total)]
Lasy.cors.all.decoy = data.frame("Method" = "Lasy-Seq", "total" = Lasy.total.decoy, "Cor" = 0.4)

data = rbind(Lasy.cors.all, DeLTa.cors.all, DeLTa.cors.all.decoy, Lasy.cors.all.decoy)

library(ggplot2)
pdf("Fig6a.pdf")
ggplot(data, aes(x=as.factor(total),y=Cor, group = interaction(Method,total), fill = Method))+
  geom_boxplot()+
  # ylim(c(0.6,1))+
  # scale_color_manual(values = c("gray","red"))+
  scale_fill_manual(values = c("red","gray"))+
  theme_bw()
dev.off()

############# Fig. 6bc DeLTa-Seq VS Lasy-Seq (Marged fasta)
# Load target
target = read.table("target_list.txt", header = F, stringsAsFactors = F)[,1]
target = matrix(unlist(strsplit(target,"\\.")), byrow = T, ncol = 2)[,1]
target = unique(target)
target.short = target

id = read.table("Araport11_genes.201606.transcript_ERCC_Virus3981.id", stringsAsFactors = F)[,1]
tmp = data.frame("gene.id" = id[grep("^AT",id)], "short.gene.id" = matrix(unlist(strsplit(id[grep("^AT",id)],"\\.")), byrow = T, ncol = 2)[,1], stringsAsFactors = F)
tmp = tmp[is.element(tmp$short.gene.id,target),]
target = tmp$gene.id

# load rawcnt
files = list.files("RSEM_out_UMI/")
files = files[grep("genes", files)]
DeLTa.files = files[grep("DeLTa", files)]
Lasy.files = files[grep("Lasy", files)]
DeLTa.total = sort(as.integer(unique(lapply(strsplit(DeLTa.files,"_"), FUN = function(x){return(x[3])}))))
Lasy.total = sort(as.integer(unique(lapply(strsplit(Lasy.files,"_"), FUN = function(x){return(x[3])}))))

Lasy.rawcnt = list()
for (t in Lasy.total){
  rawcnt = NULL
  for (r in 1:12) {
    file = files[grep(sprintf("Lasy_%02d_%s_",r,t), files)]
    file = sprintf("data/uniq_subsampling/%s",file)
    tmp = read.table(file, header = T, stringsAsFactors = F, row.names = 1)
    tmp = tmp[target,]$expected_count/sum(tmp[id[grep("^AT",id)],]$expected_count)*10^6
    tmp2 = NULL
    for(gn in target.short){
      tmp2 = c(tmp2,log2(sum(tmp[grep(gn, target)])+1))
    }
    tmp = tmp2
    names(tmp) = target.short
    rawcnt = cbind(rawcnt, tmp)
  }
  rawcnt = list(rawcnt)
  names(rawcnt) = t
  Lasy.rawcnt = c(Lasy.rawcnt, rawcnt)
}

DeLTa.rawcnt = list()
for (t in DeLTa.total){
  rawcnt = NULL
  for (r in 1:12) {
    file = files[grep(sprintf("DeLTa_%02d_%s_",r,t), files)]
    file = sprintf("RSEM_out_UMI/%s",file)
    tmp = read.table(file, header = T, stringsAsFactors = F, row.names = 1)
    tmp = tmp[target,]$expected_count/sum(tmp[id[grep("^AT",id)],]$expected_count)*10^6
    tmp2 = NULL
    for(gn in target.short){
      tmp2 = c(tmp2,log2(sum(tmp[grep(gn, target)])+1))
    }
    tmp = tmp2
    rawcnt = cbind(rawcnt, tmp)
  }
  rawcnt = list(rawcnt)
  names(rawcnt) = t
  DeLTa.rawcnt = c(DeLTa.rawcnt, rawcnt)
}



# Against marged rawcnt
DeLTa.all = read.table("RSEM_out_UMI/DeLTa.genes.results", header = T, row.names = 1)
DeLTa.all$expected_count = DeLTa.all$expected_count/sum(DeLTa.all[id[grep("^AT",id)],]$expected_count)*10^6
tmp2 = NULL
tmp = DeLTa.all[target,]$expected_count
for(gn in target.short){
  tmp2 = c(tmp2,log2(sum(tmp[grep(gn, target)])+1))
}
names(tmp2) = target.short
DeLTa.all = tmp2

Lasy.all = read.table("RSEM_out_UMI/Lasy.genes.results", header = T, row.names = 1)
Lasy.all$expected_count = Lasy.all$expected_count/sum(Lasy.all[id[grep("^AT",id)],]$expected_count)*10^6
tmp2 = NULL
tmp = Lasy.all[target,]$expected_count
for(gn in target.short){
  tmp2 = c(tmp2,log2(sum(tmp[grep(gn, target)])+1))
}
names(tmp2) = target.short
Lasy.all = tmp2

DeLTa.cors.all = NULL
for(i in 1:length(DeLTa.rawcnt)){
  DeLTa.cors.all = cbind(DeLTa.cors.all, apply(DeLTa.rawcnt[[i]], 2, FUN = function(x){
    return(cor(x,DeLTa.all))
  })
  )
}
colnames(DeLTa.cors.all) = DeLTa.total

Lasy.cors.all = NULL
for(i in 1:length(Lasy.rawcnt)){
  Lasy.cors.all = cbind(Lasy.cors.all, apply(Lasy.rawcnt[[i]], 2, FUN = function(x){
    return(cor(x,Lasy.all))
  })
  )
}
colnames(Lasy.cors.all) = Lasy.total


# plot
library(reshape2)
Lasy.cors.all = melt(Lasy.cors.all)
Lasy.cors.all = data.frame("Method" = "Lasy-Seq", "total" = Lasy.cors.all[,2], "Cor" = Lasy.cors.all[,3])

DeLTa.cors.all = melt(DeLTa.cors.all)
DeLTa.cors.all = data.frame("Method" = "DeLTa-Seq", "total" = DeLTa.cors.all[,2], "Cor" = DeLTa.cors.all[,3])

DeLTa.total.decoy = Lasy.total[!is.element(Lasy.total, DeLTa.total)]
DeLTa.cors.all.decoy = data.frame("Method" = "DeLTa-Seq", "total" = DeLTa.total.decoy, "Cor" = 0.5)
Lasy.total.decoy = DeLTa.total[!is.element(DeLTa.total, Lasy.total)]
Lasy.cors.all.decoy = data.frame("Method" = "Lasy-Seq", "total" = Lasy.total.decoy, "Cor" = 0.5)

data = rbind(Lasy.cors.all, DeLTa.cors.all, DeLTa.cors.all.decoy, Lasy.cors.all.decoy)



library(ggplot2)
pdf("Fig.6c.pdf", width = 14)
ggplot(data, aes(x=as.factor(total),y=Cor, group = interaction(Method,total), fill = Method))+
  geom_boxplot()+
  # ylim(c(0.6,1))+
  # scale_color_manual(values = c("gray","red"))+
  scale_fill_manual(values = c("red","gray"))+
  theme_bw()
dev.off()

# plot each gene
library(dplyr)
library(ggbeeswarm)
pdf("Fig.6b.pdf",width = 14)
for(g in 1:length(target.short)){
  data = NULL
  gn = target.short[g]
  cat(gn)
  cat("\n")
  for(i in 1:length(DeLTa.rawcnt)){
    DeLTa = DeLTa.rawcnt[[i]][g,]
    data = rbind(data,data.frame("Method" = "DeLTa", "log2rpm" = DeLTa ,"total" = names(DeLTa.rawcnt)[i]))
  }
  for(i in 1:length(Lasy.rawcnt)){
    Lasy = Lasy.rawcnt[[i]][g,]
    data = rbind(data,data.frame("Method" = "Lasy", "log2rpm"= Lasy, "total" = names(Lasy.rawcnt)[i]))
  }
  data$Method = factor(data$Method, levels = c( "Lasy","DeLTa"))
  data.each = data
  data = data %>%
    group_by(interaction(Method, total)) %>%
    summarize(mean = mean(log2rpm), sd = sd(log2rpm))
  tmp = matrix(unlist(strsplit(as.character(data$`interaction(Method, total)`),"\\.")), byrow = T, ncol = 2)
  data = data.frame("Method" = tmp[,1], "Total" = tmp[,2], "Mean" = data$mean, "SD" = data$sd, stringsAsFactors = F)
  level = sort(unique(as.integer(data$Total)))
  data$Total = factor(data$Total, levels = level)
  data.each$total = factor(as.character(data.each$total), levels = level)
  data$Method = factor(data$Method, levels = c("Lasy","DeLTa"))
  gg =ggplot(data, aes(x = Total, y = Mean, group = Method, fill = Method))+
    # geom_beeswarm(data = data.each, aes(x=total, y = log2rpm, color = Method), size = 0.5)+
    geom_line()+
    geom_violin(data = data.each[data.each$Method=="Lasy",], aes(x=total, y = log2rpm, group = total))+
    geom_violin(data = data.each[data.each$Method=="DeLTa",], aes(x=total, y = log2rpm, group = total))+
    geom_line()+
    ggtitle(gn)+
    scale_fill_manual(values = c("gray","red"))+
    ylim(c(0,20))+
    theme_bw()
  
  print(gg)
}
dev.off()

# Fig S2ac RNA-Seq of zebrafish from lysate and purified
files = list.files("zebrafish/")
rawcnt.zebra = NULL
for (i in 1:length(files)) {
  tmp = read.table(sprintf("zebrafish/%s", files[i]), header = T)
  rawcnt.zebra = cbind(rawcnt.zebra, tmp$expected_count)
}
rownames(rawcnt.zebra) = tmp$gene_id
colnames(rawcnt.zebra) = c(sprintf("Z.lysate.%s",1:6),sprintf("Z.purified.%s",1:6))
library(TCC)
param_FDR <- 0.05 
data <- rawcnt.zebra
data.cl <- rep(c(1,2),each = 6 )
tcc <- new("TCC", data, data.cl)
tcc <- calcNormFactors(tcc, norm.method="tmm", test.method="edger", iteration=3, FDR=0.1, floorPDEG=0.05)
tcc <- estimateDE(tcc, test.method="edger", FDR=param_FDR)
result <- getResult(tcc, sort=FALSE)
rownames(result) = result$gene_id
sum(tcc$stat$q.value < param_FDR) 

DEG  =  tcc$gene_id[tcc$stat$q.value < param_FDR]
DEG.high = DEG[result[DEG,]$m.value<0]
DEG.low = DEG[result[DEG,]$m.value>0]
length(DEG.high)
length(DEG.low)

fai = read.table("Danio_rerio.GRCz11.cdna.all.fa.fai", stringsAsFactors = F)
rownames(fai) = fai$V1
length.df = data.frame(length = fai[,2], flag = "others", stringsAsFactors = F)
rownames(length.df) = fai[,1]
length.df[DEG.high,2] = "High in Dlasy-Seq (FDR<0.05)"
length.df[DEG.low,2] = "Low in Dlasy-Seq (FDR<0.05)"

#log2rpm plot
rawcnt.marge = NULL
tmp = apply(rawcnt.zebra[,1:6], 1, sum)
rawcnt.marge = cbind(rawcnt.marge,tmp)
tmp = apply(rawcnt.zebra[,7:12], 1, sum)
rawcnt.marge = cbind(rawcnt.marge,tmp)

log2rpm.zebra = t(log2((t(rawcnt.zebra)/colSums(rawcnt.zebra)*10^6)+1))

log2rpm = t(log2((t(rawcnt.marge)/colSums(rawcnt.marge)*10^6)+1))
pdf("FigS2a.pdf")
plot(log2rpm[,2],log2rpm[,1], main = sprintf("Pearson's correlation: %s",cor(log2rpm[,1], log2rpm[,2])), ylab = "Dlasy", xlab = "Lasy" , pch = 16, ylim = c(0,15), xlim = c(0,15), col = "gray")
points(log2rpm[DEG.high,2],log2rpm[DEG.high,1], col = "red", pch = 16)
points(log2rpm[DEG.low,2],log2rpm[DEG.low,1], col = "blue", pch = 16)
lines(c(-100,100),c(-100,100), col = "black")
dev.off()

length.df$flag = factor(length.df$flag, levels = c("others", "High in Dlasy-Seq (FDR<0.05)", "Low in Dlasy-Seq (FDR<0.05)")) 
pdf("FigS2c.pdf")
g1 = ggplot(length.df[length.df$flag=="High in Dlasy-Seq (FDR<0.05)",], aes(x=length)) +
  geom_histogram(aes(y = ..ncount../sum(..ncount..)), position = "identity",  bins = 100, fill = "red") +
  xlim(c(0,20000))+
  ylim(c(0,0.1))+
  theme_bw() + 
  theme(legend.position = "bottom")
g2 = ggplot(length.df[length.df$flag=="Low in Dlasy-Seq (FDR<0.05)",], aes(x=length)) +
  geom_histogram(aes(y = ..ncount../sum(..ncount..)), position = "identity",  bins = 100, fill = "blue") +
  xlim(c(0,20000))+
  ylim(c(0,0.1))+
  theme_bw() + 
  theme(legend.position = "bottom")
g3 = ggplot(length.df[length.df$flag=="others",], aes(x=length)) +
  geom_histogram(aes(y = ..ncount../sum(..ncount..)), position = "identity",  bins = 100, fill = "gray") +
  xlim(c(0,20000))+
  ylim(c(0,0.1))+
  theme_bw() + 
  theme(legend.position = "bottom")
grid.arrange(g1,g2,g3, ncol = 1)
dev.off()

mean(length.df$length[length.df$flag=="others"]) # 1537.608
mean(length.df$length[length.df$flag!="others"]) # 798.7264
# Fisher's test
t =1000
mat = matrix(c(nrow(fai), sum(fai$V2<t), length(DEG.high), sum(fai[DEG.high,]$V2<t)), byrow = T, nco = 2)
fisher.test(mat)


#Fig S2bd RNA-Seq of yeast from lysate and purified
files = list.files("yeast")
rawcnt.yeast = NULL
for (i in 1:length(files)) {
  tmp = read.table(sprintf("yeast/%s", files[i]), header = T)
  rawcnt.yeast = cbind(rawcnt.yeast, tmp$expected_count)
}
rownames(rawcnt.yeast) = tmp$gene_id
pplot(colSums(rawcnt.yeast))
rawcnt.yeast = rawcnt.yeast[,colSums(rawcnt.yeast)>10000]
colnames(rawcnt.yeast) = c(sprintf("Y.lysate-%s",1:6),sprintf("Y.purified-%s",1:6),sprintf("Y.hot-phenol-%s",1:6))
library(TCC)
param_FDR <- 0.05
data <- rawcnt.yeast[,1:12]
data.cl <- rep(c(1,2),each = 6 )
tcc <- new("TCC", data, data.cl)
tcc <- calcNormFactors(tcc, norm.method="tmm", test.method="edger", iteration=3, FDR=0.1, floorPDEG=0.05)
tcc <- estimateDE(tcc, test.method="edger", FDR=param_FDR)
result <- getResult(tcc, sort=FALSE)
rownames(result) = result$gene_id
result = cbind(TCC::getNormalizedData(tcc),result[,-1])
write.csv(result, file = "lysate-VS-purified.csv")
sum(tcc$stat$q.value < param_FDR)

DEG  =  tcc$gene_id[tcc$stat$q.value < param_FDR]
DEG.high = DEG[result[DEG,]$m.value<0]
DEG.low = DEG[result[DEG,]$m.value>0]
length(DEG.high)
length(DEG.low)

fai = read.table("Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.fai", stringsAsFactors = F)
rownames(fai) = fai$V1
length.df = data.frame(length = fai[,2], flag = "others", stringsAsFactors = F)
rownames(length.df) = fai[,1]
length.df[DEG.high,2] = "High in Dlasy-Seq (FDR<0.05)"
length.df[DEG.low,2] = "Low in Dlasy-Seq (FDR<0.05)"
length.df$flag = factor(length.df$flag, levels = c("others", "High in Dlasy-Seq (FDR<0.05)", "Low in Dlasy-Seq (FDR<0.05)")) 

#log2rpm plot
rawcnt.marge = NULL
tmp = apply(rawcnt.yeast[,1:6], 1, sum)
rawcnt.marge = cbind(rawcnt.marge,tmp)
tmp = apply(rawcnt.yeast[,7:12], 1, sum)
rawcnt.marge = cbind(rawcnt.marge,tmp)


log2rpm = t(log2((t(rawcnt.marge)/colSums(rawcnt.marge)*10^6)+1))
pdf("FigS2b.pdf")
plot(log2rpm[,2],log2rpm[,1], main = sprintf("Pearson's correlation: %s",cor(log2rpm[,1], log2rpm[,2])), ylab = "Dlasy", xlab = "Lasy" , pch = 16, ylim = c(0,15), xlim = c(0,15), col = "gray")
points(log2rpm[DEG.high,2],log2rpm[DEG.high,1], col = "red", pch = 16)
points(log2rpm[DEG.low,2],log2rpm[DEG.low,1], col = "blue", pch = 16)
lines(c(-100,100),c(-100,100), col = "black")
dev.off()

pdf("FigS2d.pdf")
g1 = ggplot(length.df[length.df$flag=="High in Dlasy-Seq (FDR<0.05)",], aes(x=length)) +
  geom_histogram(aes(y = ..ncount../sum(..ncount..)), position = "identity",  bins = 100, fill = "red") +
  xlim(c(0,8000))+
  ylim(c(0,0.2))+
  theme_bw() + 
  theme(legend.position = "bottom")
g2 = ggplot(length.df[length.df$flag=="Low in Dlasy-Seq (FDR<0.05)",], aes(x=length)) +
  geom_histogram(aes(y = ..ncount../sum(..ncount..)), position = "identity",  bins = 100, fill = "blue") +
  xlim(c(0,8000))+
  ylim(c(0,0.2))+
  theme_bw() + 
  theme(legend.position = "bottom")
g3 = ggplot(length.df[length.df$flag=="others",], aes(x=length)) +
  geom_histogram(aes(y = ..ncount../sum(..ncount..)), position = "identity",  bins = 100, fill = "gray") +
  xlim(c(0,8000))+
  ylim(c(0,0.2))+
  theme_bw() + 
  theme(legend.position = "bottom")
grid.arrange(g1,g2,g3, ncol = 1)

dev.off()
mean(length.df$length[length.df$flag=="others"]) # 1537.608
mean(length.df$length[length.df$flag!="others"]) # 798.7264
# Fisher's test
t =1000
mat = matrix(c(nrow(fai), sum(fai$V2<t), length(DEG.high), sum(fai[DEG.high,]$V2<t)), byrow = T, nco = 2)
fisher.test(mat)


