#code for Figure 1
chrlens <- c(598660471, 700547350, 498638509, 787782082,  812755788, 656544405, 754128162, 851934019, 619618552, 
              754227511, 673810255, 518332611, 713360525, 714697677, 569951140, 622669697, 731188232, 495380293, 744491536, 764072961, 642921167)
chrnames <- c("1A","1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D", "4A", "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D", "7A", "7B","7D")
dfcovg <- read.table("coverage.txt", header=F, row.names = 1, fill=TRUE)
rownames(dfcovg) <- substr(rownames(dfcovg),2,3)
dfcovg <- dfcovg[chrnames,]
sums <- rowSums(dfcovg)

#This code is for generating supplementary Figure 1

par(mar=c(0,2,0,0), mfrow = c(23,1))
plot(c(1:1000),dfcovg["1A",]/100000,type="l", xaxt="n", yaxt="n"); abline(v=213, col="green"); mtext(side = 2, text = "1A", line = 0.5)
plot(c(1:1000),dfcovg["1B",]/100000,type="l", xaxt="n", yaxt="n"); abline(v=240, col="green"); mtext(side = 2, text = "1B", line = 0.5)
plot(c(1:1000),dfcovg["1D",]/100000,type="l", xaxt="n", yaxt="n"); abline(v=175, col="green"); mtext(side = 2, text = "1D", line = 0.5)
plot(c(1:1000),dfcovg["2A",]/100000,type="l", xaxt="n", yaxt="n"); abline(v=335, col="green"); mtext(side = 2, text = "2A", line = 0.5)
plot(c(1:1000),dfcovg["2B",]/100000,type="l", xaxt="n", yaxt="n"); abline(v=358, col="green"); mtext(side = 2, text = "2B", line = 0.5)
plot(c(1:1000),dfcovg["2D",]/100000,type="l", xaxt="n", yaxt="n"); abline(v=273, col="green"); mtext(side = 2, text = "2D", line = 0.5)
plot(c(1:1000),dfcovg["3A",]/100000,type="l", xaxt="n", yaxt="n"); abline(v=330, col="green");  mtext(side = 2, text = "3A", line = 0.5)
plot(c(1:1000),dfcovg["3B",]/100000,type="l", xaxt="n", yaxt="n"); abline(v=355, col="green"); mtext(side = 2, text = "3B", line = 0.5)
plot(c(1:1000),dfcovg["3D",]/100000,type="l", xaxt="n", yaxt="n"); abline(v=245, col="green"); mtext(side = 2, text = "3D", line = 0.5)
plot(c(1:1000),dfcovg["4A",]/100000,type="l", xaxt="n", yaxt="n"); abline(v=313, col="green"); mtext(side = 2, text = "4A", line = 0.5)
plot(c(1:1000),dfcovg["4B",]/100000,type="l", xaxt="n", yaxt="n"); abline(v=300, col="green"); mtext(side = 2, text = "4B", line = 0.5)
plot(c(1:1000),dfcovg["4D",]/100000,type="l", xaxt="n", yaxt="n"); abline(v=205, col="green");mtext(side = 2, text = "4D", line = 0.5)
plot(c(1:1000),dfcovg["5A",]/100000,type="l", xaxt="n", yaxt="n"); abline(v=253, col="green"); mtext(side = 2, text = "5A", line = 0.5)
plot(c(1:1000),dfcovg["5B",]/100000,type="l", xaxt="n", yaxt="n"); abline(v=200, col="green"); mtext(side = 2, text = "5B", line = 0.5)
plot(c(1:1000),dfcovg["5D",]/100000,type="l", xaxt="n", yaxt="n"); abline(v=190, col="green"); mtext(side = 2, text = "5D", line = 0.5)
plot(c(1:1000),dfcovg["6A",]/100000,type="l", xaxt="n", yaxt="n"); abline(v=290, col="green"); mtext(side = 2, text = "6A", line = 0.5)
plot(c(1:1000),dfcovg["6B",]/100000,type="l", xaxt="n", yaxt="n"); abline(v=333, col="green"); mtext(side = 2, text = "6B", line = 0.5)
plot(c(1:1000),dfcovg["6D",]/100000,type="l", xaxt="n", yaxt="n"); abline(v=235, col="green"); mtext(side = 2, text = "6D", line = 0.5)
plot(c(1:1000),dfcovg["7A",]/100000,type="l", xaxt="n", yaxt="n"); abline(v=365, col="green"); mtext(side = 2, text = "7A", line = 0.5)
plot(c(1:1000),dfcovg["7B",]/100000,type="l", xaxt="n", yaxt="n"); abline(v=310, col="green"); mtext(side = 2, text = "7B", line = 0.5)
plot(c(1:1000),dfcovg["7D",]/100000,type="l", yaxt="n"); abline(v=340, col="green"); mtext(side = 2, text = "7D", line = 0.5)
mtext(side = 1, text = "Chromosome Position(1Mbp)", line = 2.5)
     

bed <- read.table("iwgsc_refseqv2.1_annotation_200916_HC.bed", sep="\t", header=F)
min(bed$V2)
max(bed$V2)
min(bed$V3)
max(bed$V3)

write.table(bed, file="iwgsc_refseqv2.1_annotation_200916_HC2.bed", sep="\t", quote = F,  row.names = F, col.names=F)
#####################################################################################
##############################
#this code is for visualizing gaps graphs
#java -jar MFbio.jar --task mafinfo --srcdir /CoreDetector_outdir/msa.maf.gz
chromspecies <- c("CS","spelt","mace","julius","norin61","lancer","mattis","jagger","landmark","arina","stanley", "paragon")
species <- c("CS","spelt","mace","julius","norin61","lancer","mattis","jagger","landmark","arina","stanley", "paragon","weebil","robigus","claire","cadenza")
gaps <- read.table("gaps.txt", sep="\t", header=F)

dim(gaps)
#13,833,261   6   gaps larger than 100bp 
lratio <- function(read){
  #returns ratio of lowercase characters 
  lread <- tolower(read)
  lst1<- unlist(strsplit(read,""))
  lst2<- unlist(strsplit(lread,""))
  return(length(which(lst1==lst2))/nchar(read))
  
}
colnames(gaps) <- c("hsp", "specie", "chrom", "start", "realstart", "len")
length(unique(gaps$hsp))
#661,022  hsp's that contain gaps larger than 100bp
#aggregate based on chromosome, start postion and number of cultivars related to that gap
gap_agg1 <- aggregate(specie~hsp+start, gaps , FUN=function(v){return(length(unique(v)))})
gap_agg1 <- gap_agg1[order(gap_agg1$hsp, gap_agg1$start),]
max(gap_agg1$specie)
#15
#unique insertion
uins <- gap_agg1[which(gap_agg1$specie==15)  ,]
dim(uins)
#602754
gapmerged <- merge(gaps, uins, by.x=c(1,4), by.y=c(1,2))
dim(gapmerged)
#9,041,310       7

merged10k <- gapmerged[gapmerged$len>9000,]
dim(merged10k)[1]/15
#7424 number of unique insertions larger than 10K
#now for each gap alignment we find the missing one that would be insertion specie
agg_func <- function(x){
  #cat(paste(length(x[,3]),","))
  if(dim(x)[1]==15){
    spc <- species[which(!species%in%x[,3])]
    #find chr and position of CS
    vv <- data.frame()
    if (spc!="CS"){
      vv <- x[x[,3]=="CS", c(1,2,4,6)]
    }else{
      #stanley is closest in terms of CNV to CS
      vv <- x[x[,3]=="stanley", c(1,2,4,6)]
    }
    return (list(vv[1,1], vv[1,2], spc,  vv[1,3], vv[1,4]))
  }else
    return (list(NA,NA,NA,NA,NA))
}
res <- by(merged10k, merged10k[,c(1,2)], FUN=agg_func, simplify=FALSE)
res <- unlist(as.list(res))
length(res)
#31,110
insdf <- data.frame(matrix(res, ncol=5, byrow =T))
dim(insdf)
#6,222  5
insdf <- insdf[-which(rowSums(is.na(insdf))==5),]
dim(insdf)
#6217  5
colnames(insdf) <- c("hsp","start","specie","chrom","len")
insdf$hsp <- as.numeric(insdf$hsp)
insdf$start <- as.numeric(insdf$start)
insdf$len <- as.numeric(insdf$len)
max(insdf$len)
#104,737
insdf <- insdf[order(insdf$hsp, insdf$start),]

write.table(insdf, file="uins.txt", quote=F, sep="\t", row.names = F, col.names=F)

insdf[insdf$chrom=="1A" & insdf$specie=="CS",]
#hsp start specie chrom   len
#27023 21872     CS    1A 11611
merged10k[merged10k$hsp==27023,]

insdf[which(insdf$len==104737),]
#stanley 1D 15,289,886 104,737
insdf <- insdf[order(-insdf$len),]
# first one in 1A is hsp 89955
merged10k[merged10k$hsp==89955,]
#......
#89955  1201  norin61    1A 510489366 86450       15
#......

#to see more details of that insertion
# we extract it using task extractmafsegment and the result is in the file stanley_uins.fa that later blasted to NCBI and the result is in manuscript
#hsp start specie.x chrom realstart    len specie.y
#290943 15800    arina    1D  14479177 104737       15
merged10k[merged10k$len==104737,]

gapmerged[gapmerged$len==9999,]

insdf[which(insdf$len==97957),]
#3159  arina    3B 247881852 97957
merged10k[merged10k$len==97957,]
#hsp :812033  start:5328 len:97957


#see if chinese spring has unique insertion larger than 10K
insdf[insdf$specie=="CS",]
#5023     CS    7B 181434262 17128
#5154     CS    1A  49406678 11611
#5720     CS    3B  44205678 14810
#5751     CS    4A  75774923 12440

#this is for manuscript Figure  4c
sum_ins <- aggregate(len~specie+chrom, insdf , FUN=length)
sum_ins <- sum_ins[sum_ins$specie%in%chromspecies,]
min(sum_ins$len)
#1
max(sum_ins$len)
#250
tab <-  xtabs(len ~ specie + chrom, data = sum_ins)

temat <- as.matrix(tab)
heatmap(temat, scale='none' )

#Figure 4a of manuscript
gaps_12 <- gaps[which(gaps$specie%in%chromspecies), ]
dim(gaps)
#13,833,261
sum_gaps <- aggregate(len~specie+chrom, gaps_12 , FUN=length)
min(sum_gaps$len)
#15,569
max(sum_gaps$len)
#61,801
tab <-  xtabs(len ~ specie + chrom, data = sum_gaps)

temat <- as.matrix(tab)
heatmap(temat, scale='none' )
#Figure 4b of manuscript
gaps_12_10k <- gaps_12[which(gaps_12$len>10000),]
dim(gaps_12_10k)
#134,644   6
sum_gaps <- aggregate(len~specie+chrom, gaps_12_10k , FUN=length)
min(sum_gaps$len)
#68
max(sum_gaps$len)
#1717
tab <-  xtabs(len ~ specie + chrom, data = sum_gaps)

temat <- as.matrix(tab)
heatmap(temat, scale='none' )


#gaps larger that 10K
length(which(gaps_12$len>10000))
#134,644
#gaps larger that 30K
length(which(gaps$len>30000))
#5548
################################################
#Find largest unique insertion in cactus result
gaps <- read.table("gaps_1A_cactus.txt", sep="\t", header=F)

dim(gaps)
#118, 383   6   gaps larger than 100bp 
lratio <- function(read){
  #returns ratio of lowercase characters 
  lread <- tolower(read)
  lst1<- unlist(strsplit(read,""))
  lst2<- unlist(strsplit(lread,""))
  return(length(which(lst1==lst2))/nchar(read))
  
}
colnames(gaps) <- c("hsp", "specie", "chrom", "start", "realstart", "len")
length(unique(gaps$hsp))
#40,053  hsp's that contain gaps larger than 100bp
#aggregate based on chromosome, start postion and number of cultivars related to that gap
gap_agg1 <- aggregate(specie~hsp+start, gaps , FUN=function(v){return(length(unique(v)))})
gap_agg1 <- gap_agg1[order(gap_agg1$hsp, gap_agg1$start),]
max(gap_agg1$specie)
#15
uins <- gap_agg1[which(gap_agg1$specie==15)  ,]
dim(uins)
#706
gapmerged <- merge(gaps, uins, by.x=c(1,4), by.y=c(1,2))
dim(gapmerged)
#10,624
max(gapmerged$len)
#9999
#now for each gap alignment we find the missing one that would be insertion specie
agg_func <- function(x){
  #cat(paste(length(x[,3]),","))
  if(dim(x)[1]==15){
    spc <- species[which(!species%in%x[,3])]
    #find chr and position of CS
    vv <- data.frame()
    if (spc!="CS"){
      vv <- x[x[,3]=="CS", c(1,2,4,6)]
    }else{
      #stanley is closest in terms of CNV to CS
      vv <- x[x[,3]=="stanley", c(1,2,4,6)]
    }
    return (list(vv[1,1], vv[1,2], spc,  vv[1,3], vv[1,4]))
  }else
    return (list(NA,NA,NA,NA,NA))
}
res <- by(gapmerged, gapmerged[,c(1,2)], FUN=agg_func, simplify=FALSE)
res <- unlist(as.list(res))
length(res)
#3530
insdf <- data.frame(matrix(res, ncol=5, byrow =T))
dim(insdf)
#706  5
insdf <- insdf[-which(rowSums(is.na(insdf))==5),]
dim(insdf)
#690  5
colnames(insdf) <- c("hsp","start","specie","chrom","len")
insdf$hsp <- as.numeric(insdf$hsp)
insdf$start <- as.numeric(insdf$start)
insdf$len <- as.numeric(insdf$len)
max(insdf$len)
#9999
insdf <- insdf[order(insdf$hsp, insdf$start),]

write.table(insdf, file="uins_cactus.txt", quote=F, sep="\t", row.names = F, col.names=F)

insdf[which(insdf$len==9999),]
# hsp      start  specie chrom  len
# 3787     1     CS     1A     9999
insdf <- insdf[order(-insdf$len),]

#to see more details of that insertion
gapmerged[gapmerged$len==9999,]

gaps <- gaps[gaps$hsp==3787,]
#####################################################################################################
#This code is for finding genes in core alignment
#java -jar  MFbio.jar  --task maf2bed --srcdir /outDir/core_msa.maf.gz  --destdir core_CS.bed --p1 CS
#Now intersection of core and genes
#bedtools intersect -wao -a /scratchdata1/users/a1195806/mario/references/CS2.1/whole/iwgsc_refseqv2.1_annotation_200916_HC.bed  -b /scratchdata1/users/a1195806/mario/projects/msa/wheat/minimap_out_sm_15/core_CS.bed > /scratchdata1/users/a1195806/mario/projects/msa/wheat/minimap_out_sm_15/intersect.bed
inter <- read.table("intersect.bed", sep="\t", header=F)
dim(inter)
#133,641      7
#aggregate based on gene region (v1,v2,v3) and return sum of core coverage for that region
agg <- aggregate(inter$V7, by=list(inter$V1, inter$V2, inter$V3), FUN=sum, na.rm=TRUE)
dim(agg)
# 105463      3
dim(agg[agg$x==0,])
#8929   3  8.4% of genes not covered at all
8929/105463
sum(agg$x)
# 308174096

sum(agg$Group.3-agg$Group.2)
#368986606
308174096/368986606
#83.5% genic region covered 

length(which((agg$Group.3-agg$Group.2 -2)<= agg$x))
#68665
68665/105463
#0.651   65% of gene are fully covered
#now draw the graph of gene area  ; Supplementary Figure 2
agg$million <- floor(agg$Group.2/1000000)
aggm <- aggregate(agg$x, by=list(agg$Group.1, agg$million), FUN=sum, na.rm=TRUE)
max(aggm$x)
#156273
aggm[aggm$x==max(aggm$x),]
#7B     592 156273
ymax <- 150000 
dev.off();
plot.new()
par(mar=c(0,2,0,0), mfrow = c(23,1))
plot(aggm[aggm$Group.1=="1A","Group.2"],aggm[aggm$Group.1=="1A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=215, col="green"); mtext(side = 2, text = "1A", line = 0.5)
plot(aggm[aggm$Group.1=="1B","Group.2"],aggm[aggm$Group.1=="1B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=240, col="green"); mtext(side = 2, text = "1B", line = 0.5)
plot(aggm[aggm$Group.1=="1D","Group.2"],aggm[aggm$Group.1=="1D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=175, col="green"); mtext(side = 2, text = "1D", line = 0.5)
plot(aggm[aggm$Group.1=="2A","Group.2"],aggm[aggm$Group.1=="2A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=335, col="green"); mtext(side = 2, text = "2A", line = 0.5)
plot(aggm[aggm$Group.1=="2B","Group.2"],aggm[aggm$Group.1=="2B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=360, col="green"); mtext(side = 2, text = "2B", line = 0.5)
plot(aggm[aggm$Group.1=="2D","Group.2"],aggm[aggm$Group.1=="2D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=275, col="green"); mtext(side = 2, text = "2D", line = 0.5)
plot(aggm[aggm$Group.1=="3A","Group.2"],aggm[aggm$Group.1=="3A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=330, col="green");  mtext(side = 2, text = "3A", line = 0.5)
plot(aggm[aggm$Group.1=="3B","Group.2"],aggm[aggm$Group.1=="3B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=355, col="green"); mtext(side = 2, text = "3B", line = 0.5)
plot(aggm[aggm$Group.1=="3D","Group.2"],aggm[aggm$Group.1=="3D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=245, col="green"); mtext(side = 2, text = "3D", line = 0.5)
plot(aggm[aggm$Group.1=="4A","Group.2"],aggm[aggm$Group.1=="4A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=250, col="green"); mtext(side = 2, text = "4A", line = 0.5)
plot(aggm[aggm$Group.1=="4B","Group.2"],aggm[aggm$Group.1=="4B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=300, col="green"); mtext(side = 2, text = "4B", line = 0.5)
plot(aggm[aggm$Group.1=="4D","Group.2"],aggm[aggm$Group.1=="4D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=205, col="green");mtext(side = 2, text = "4D", line = 0.5)
plot(aggm[aggm$Group.1=="5A","Group.2"],aggm[aggm$Group.1=="5A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=253, col="green"); mtext(side = 2, text = "5A", line = 0.5)
plot(aggm[aggm$Group.1=="5B","Group.2"],aggm[aggm$Group.1=="5B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=200, col="green"); mtext(side = 2, text = "5B", line = 0.5)
plot(aggm[aggm$Group.1=="5D","Group.2"],aggm[aggm$Group.1=="5D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=190, col="green"); mtext(side = 2, text = "5D", line = 0.5)
plot(aggm[aggm$Group.1=="6A","Group.2"],aggm[aggm$Group.1=="6A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=290, col="green"); mtext(side = 2, text = "6A", line = 0.5)
plot(aggm[aggm$Group.1=="6B","Group.2"],aggm[aggm$Group.1=="6B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=333, col="green"); mtext(side = 2, text = "6B", line = 0.5)
plot(aggm[aggm$Group.1=="6D","Group.2"],aggm[aggm$Group.1=="6D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=235, col="green"); mtext(side = 2, text = "6D", line = 0.5)
plot(aggm[aggm$Group.1=="7A","Group.2"],aggm[aggm$Group.1=="7A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=368, col="green"); mtext(side = 2, text = "7A", line = 0.5)
plot(aggm[aggm$Group.1=="7B","Group.2"],aggm[aggm$Group.1=="7B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=305, col="green"); mtext(side = 2, text = "7B", line = 0.5)
plot(aggm[aggm$Group.1=="7D","Group.2"],aggm[aggm$Group.1=="7D","x"],type="l",  yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=345, col="green"); mtext(side = 2, text = "7D", line = 0.5)
mtext(side = 1, text = "Chromosome Position(1Mbp)", line = 2.5)

notcovered <- agg[agg$x==0,]
notcovered$Group.2 <- notcovered$Group.2 + 1

#find name of genes not covered at all

gen <- read.table("iwgsc_refseqv2.1_annotation_200916_HC.gen", sep="\t", header=F)
gen$V2 <- substr(gen$V2,4,5)
nmerged <- merge(notcovered, gen , by.x=c(1,2), by.y=c(2,3), all.x = T)
dim(nmerged)
#8956   9
length(unique(nmerged$V1))
#8948
write.table(nmerged$V1, file="notfoundgenes.txt", quote=F)
#List of genes in notfoundgenes.txt  was used for gene ontoly analysis at www.comparativego.com


#############################################################################################
#This code is for Transposable Element analysis and graphs
te <- read.table("Tae.Chinese_Spring.refSeqv2.1.TE_annotation.txt", sep="\t", header=F)
dim(te)
#4731260       6
tab <- table(te$V1,te$V5)
addmargins(tab)
tab2 <- table(te$V5, te$V6)
tab2
te$V2 <- as.integer(te$V2) #this will avoid scientific notation
#te$V1 <- substr(te$V1,4,5)
write.table(te[,1:3], file="Tae.Chinese_Spring.refSeqv2.1.TE_annotation.bed", sep="\t", quote = F,  row.names = F, col.names=F)
#bedtools intersect -wao -a /scratchdata1/users/a1195806/mario/references/CS2.1/whole/Tae.Chinese_Spring.refSeqv2.1.TE_annotation.bed  -b /scratchdata1/users/a1195806/mario/projects/msa/wheat/minimap_out_sm_15_7/msa_cnv_sorted.bed > /scratchdata1/users/a1195806/mario/projects/msa/wheat/minimap_out_sm_15_7/intersect_te.bed
inter_te <- read.table("intersect_te.bed", sep="\t", header=F)
dim(inter_te)
#6,408,589       7

#aggregate based on TE region and return sum of core coverage of that region
aggte <- aggregate(inter_te$V7, by=list(inter_te$V1, inter_te$V2, inter_te$V3), FUN=sum, na.rm=TRUE)
dim(aggte)
# 4,731,106       4   however this 4.7 is not unique TE because it contains match_part as well
sum(aggte$x)
# 8,594,324,451
sum(inter_te$V7)
#8,594,324,451

sum(agg$Group.3-agg$Group.2)
#12,097,462,561
#8594324451/12097462561
#0.71  TE region covered

#to find unique entries of TE we need to merge it with Tae.Chinese_Spring.refSeqv2.1.TE_annotation.txt
temerged <-  merge(inter_te, te , by=c(1,2,3), all.x = T)
dim(temerged)
#6,408,983      10
#now aggregate based on V4.y (TE ID)
#now filter for complete TE
#merged <- merged[merged$V6.y=="complete",]
#dim(merged)
# 2,437,596      11
agg2 <- aggregate(temerged$V7, by=list(temerged$V4.y), FUN=sum, na.rm=TRUE)
dim(agg2)
#3,997,530       2 :This is total number of unique TE
#1109706       2  :This is total number of unique complete TE
dim(agg2[agg2$x==0,])
#699420      2  17.5% of TE not covered at all
#699420/3997530

#to find number of TE fully covered we need more information:
temerged$len <- temerged$V3 - temerged$V2 + 1 
agg3 <- aggregate(temerged$len, by=list(temerged$V4.y), FUN=mean, na.rm=TRUE)
dim(agg3)
#3,997,530       2
#reorder agg3 based on agg2
row.names(agg3) <- agg3$Group.1
agg3 <- agg3[agg2$Group.1, ]
length(which((agg3$x-3)<= agg2$x))
#2,476,922
#2476922/3997530
#61.96% fully covered



#now we want separate positions of TE in core based on different TE types
levels(as.factor(temerged$V5.y))
#"DHH" "DTA" "DTC" "DTH" "DTM" "DTT" "DTX" "DXX" "RIX" "RLC" "RLG" "RLX" "SIX" "XXX"
agg4 <- aggregate(temerged$V7, by=list(temerged$V5.y, temerged$V1), FUN=sum, na.rm=TRUE)
#use output of agg4 to fill the supplementary table 1 in the manuscript
agg4[agg4$Group.1=="DTH",]
#now extract 
stat1 <- temerged[,c(1,2,3,7,9)]
stat1 <- stat1[stat1$V7!=0,]
stat1 <- stat1[order(stat1$V1, stat1$V2),]
stat1$million <- floor(stat1$V2/1000000)
agg5 <- aggregate(stat1$V7, by=list(stat1$V1, stat1$million), FUN=sum, na.rm=TRUE)
max(agg5$x)
#1044299

agg5[agg5$x==1044299,]
#2A     530 1027102
dev.off();
plot.new()
par(mar=c(0,2,0,0), mfrow = c(23,1), lwd=0.8)
plot(agg5[agg5$Group.1=="1A","Group.2"],agg5[agg5$Group.1=="1A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,1000000)); abline(v=213, col="green"); mtext(side = 2, text = "1A", line = 0.5)
plot(agg5[agg5$Group.1=="1B","Group.2"],agg5[agg5$Group.1=="1B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,1000000)); abline(v=240, col="green"); mtext(side = 2, text = "1B", line = 0.5)
plot(agg5[agg5$Group.1=="1D","Group.2"],agg5[agg5$Group.1=="1D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,1000000)); abline(v=175, col="green"); mtext(side = 2, text = "1D", line = 0.5)
plot(agg5[agg5$Group.1=="2A","Group.2"],agg5[agg5$Group.1=="2A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,1000000)); abline(v=335, col="green"); mtext(side = 2, text = "2A", line = 0.5)
plot(agg5[agg5$Group.1=="2B","Group.2"],agg5[agg5$Group.1=="2B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,1000000)); abline(v=358, col="green"); mtext(side = 2, text = "2B", line = 0.5)
plot(agg5[agg5$Group.1=="2D","Group.2"],agg5[agg5$Group.1=="2D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,1000000)); abline(v=273, col="green"); mtext(side = 2, text = "2D", line = 0.5)
plot(agg5[agg5$Group.1=="3A","Group.2"],agg5[agg5$Group.1=="3A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,1000000)); abline(v=330, col="green");  mtext(side = 2, text = "3A", line = 0.5)
plot(agg5[agg5$Group.1=="3B","Group.2"],agg5[agg5$Group.1=="3B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,1000000)); abline(v=355, col="green"); mtext(side = 2, text = "3B", line = 0.5)
plot(agg5[agg5$Group.1=="3D","Group.2"],agg5[agg5$Group.1=="3D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,1000000)); abline(v=245, col="green"); mtext(side = 2, text = "3D", line = 0.5)
plot(agg5[agg5$Group.1=="4A","Group.2"],agg5[agg5$Group.1=="4A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,1000000)); abline(v=313, col="green"); mtext(side = 2, text = "4A", line = 0.5)
plot(agg5[agg5$Group.1=="4B","Group.2"],agg5[agg5$Group.1=="4B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,1000000)); abline(v=300, col="green"); mtext(side = 2, text = "4B", line = 0.5)
plot(agg5[agg5$Group.1=="4D","Group.2"],agg5[agg5$Group.1=="4D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,1000000)); abline(v=205, col="green");mtext(side = 2, text = "4D", line = 0.5)
plot(agg5[agg5$Group.1=="5A","Group.2"],agg5[agg5$Group.1=="5A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,1000000)); abline(v=253, col="green"); mtext(side = 2, text = "5A", line = 0.5)
plot(agg5[agg5$Group.1=="5B","Group.2"],agg5[agg5$Group.1=="5B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,1000000)); abline(v=200, col="green"); mtext(side = 2, text = "5B", line = 0.5)
plot(agg5[agg5$Group.1=="5D","Group.2"],agg5[agg5$Group.1=="5D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,1000000)); abline(v=190, col="green"); mtext(side = 2, text = "5D", line = 0.5)
plot(agg5[agg5$Group.1=="6A","Group.2"],agg5[agg5$Group.1=="6A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,1000000)); abline(v=290, col="green"); mtext(side = 2, text = "6A", line = 0.5)
plot(agg5[agg5$Group.1=="6B","Group.2"],agg5[agg5$Group.1=="6B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,1000000)); abline(v=333, col="green"); mtext(side = 2, text = "6B", line = 0.5)
plot(agg5[agg5$Group.1=="6D","Group.2"],agg5[agg5$Group.1=="6D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,1000000)); abline(v=235, col="green"); mtext(side = 2, text = "6D", line = 0.5)
plot(agg5[agg5$Group.1=="7A","Group.2"],agg5[agg5$Group.1=="7A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,1000000)); abline(v=365, col="green"); mtext(side = 2, text = "7A", line = 0.5)
plot(agg5[agg5$Group.1=="7B","Group.2"],agg5[agg5$Group.1=="7B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,1000000)); abline(v=310, col="green"); mtext(side = 2, text = "7B", line = 0.5)
plot(agg5[agg5$Group.1=="7D","Group.2"],agg5[agg5$Group.1=="7D","x"],type="l",  yaxt="n", xlim=c(0,900), ylim=c(0,1000000)); abline(v=340, col="green"); mtext(side = 2, text = "7D", line = 0.5)
mtext(side = 1, text = "Chromosome Position(1Mbp)", line = 2.5)


#now plot each TE separately
stat2<- stat1[stat1$V5.y=="RLG",]
agg6 <- aggregate(stat2$V7, by=list(stat2$V1, stat2$million), FUN=sum, na.rm=TRUE)

dev.off();
plot.new()
par(mar=c(0,2,0,0), mfrow = c(23,1))
plot(agg6[agg6$Group.1=="1A","Group.2"],agg6[agg6$Group.1=="1A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,1000000)); abline(v=215, col="green"); mtext(side = 2, text = "1A", line = 0.5)
plot(agg6[agg6$Group.1=="1B","Group.2"],agg6[agg6$Group.1=="1B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,1000000)); abline(v=240, col="green"); mtext(side = 2, text = "1B", line = 0.5)
plot(agg6[agg6$Group.1=="1D","Group.2"],agg6[agg6$Group.1=="1D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,1000000)); abline(v=175, col="green"); mtext(side = 2, text = "1D", line = 0.5)
plot(agg6[agg6$Group.1=="2A","Group.2"],agg6[agg6$Group.1=="2A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,1000000)); abline(v=335, col="green"); mtext(side = 2, text = "2A", line = 0.5)
plot(agg6[agg6$Group.1=="2B","Group.2"],agg6[agg6$Group.1=="2B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,1000000)); abline(v=360, col="green"); mtext(side = 2, text = "2B", line = 0.5)
plot(agg6[agg6$Group.1=="2D","Group.2"],agg6[agg6$Group.1=="2D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,1000000)); abline(v=275, col="green"); mtext(side = 2, text = "2D", line = 0.5)
plot(agg6[agg6$Group.1=="3A","Group.2"],agg6[agg6$Group.1=="3A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,1000000)); abline(v=330, col="green");  mtext(side = 2, text = "3A", line = 0.5)
plot(agg6[agg6$Group.1=="3B","Group.2"],agg6[agg6$Group.1=="3B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,1000000)); abline(v=355, col="green"); mtext(side = 2, text = "3B", line = 0.5)
plot(agg6[agg6$Group.1=="3D","Group.2"],agg6[agg6$Group.1=="3D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,1000000)); abline(v=245, col="green"); mtext(side = 2, text = "3D", line = 0.5)
plot(agg6[agg6$Group.1=="4A","Group.2"],agg6[agg6$Group.1=="4A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,1000000)); abline(v=250, col="green"); mtext(side = 2, text = "4A", line = 0.5)
plot(agg6[agg6$Group.1=="4B","Group.2"],agg6[agg6$Group.1=="4B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,1000000)); abline(v=300, col="green"); mtext(side = 2, text = "4B", line = 0.5)
plot(agg6[agg6$Group.1=="4D","Group.2"],agg6[agg6$Group.1=="4D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,1000000)); abline(v=205, col="green");mtext(side = 2, text = "4D", line = 0.5)
plot(agg6[agg6$Group.1=="5A","Group.2"],agg6[agg6$Group.1=="5A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,1000000)); abline(v=253, col="green"); mtext(side = 2, text = "5A", line = 0.5)
plot(agg6[agg6$Group.1=="5B","Group.2"],agg6[agg6$Group.1=="5B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,1000000)); abline(v=200, col="green"); mtext(side = 2, text = "5B", line = 0.5)
plot(agg6[agg6$Group.1=="5D","Group.2"],agg6[agg6$Group.1=="5D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,1000000)); abline(v=190, col="green"); mtext(side = 2, text = "5D", line = 0.5)
plot(agg6[agg6$Group.1=="6A","Group.2"],agg6[agg6$Group.1=="6A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,1000000)); abline(v=290, col="green"); mtext(side = 2, text = "6A", line = 0.5)
plot(agg6[agg6$Group.1=="6B","Group.2"],agg6[agg6$Group.1=="6B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,1000000)); abline(v=333, col="green"); mtext(side = 2, text = "6B", line = 0.5)
plot(agg6[agg6$Group.1=="6D","Group.2"],agg6[agg6$Group.1=="6D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,1000000)); abline(v=235, col="green"); mtext(side = 2, text = "6D", line = 0.5)
plot(agg6[agg6$Group.1=="7A","Group.2"],agg6[agg6$Group.1=="7A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,1000000)); abline(v=368, col="green"); mtext(side = 2, text = "7A", line = 0.5)
plot(agg6[agg6$Group.1=="7B","Group.2"],agg6[agg6$Group.1=="7B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,1000000)); abline(v=305, col="green"); mtext(side = 2, text = "7B", line = 0.5)
plot(agg6[agg6$Group.1=="7D","Group.2"],agg6[agg6$Group.1=="7D","x"],type="l",  yaxt="n", xlim=c(0,900), ylim=c(0,1000000)); abline(v=345, col="green"); mtext(side = 2, text = "7D", line = 0.5)
mtext(side = 1, text = "Chromosome Position(1Mbp)", line = 2.5)

#find correlation with geneics region or aggm
gene_RLG <- merge(aggm, agg6, by=c(1,2))
cor(gene_RLG$x.x, gene_RLG$x.y)
#-0.2062

stat2<- stat1[stat1$V5.y=="RLC",]
agg7 <- aggregate(stat2$V7, by=list(stat2$V1, stat2$million), FUN=sum, na.rm=TRUE)
dev.off();
plot.new()

par(mar=c(0,2,0,0), mfrow = c(23,1))
plot(agg7[agg7$Group.1=="1A","Group.2"],agg7[agg7$Group.1=="1A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,300000)); abline(v=215, col="green"); mtext(side = 2, text = "1A", line = 0.5)
plot(agg7[agg7$Group.1=="1B","Group.2"],agg7[agg7$Group.1=="1B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,300000)); abline(v=240, col="green"); mtext(side = 2, text = "1B", line = 0.5)
plot(agg7[agg7$Group.1=="1D","Group.2"],agg7[agg7$Group.1=="1D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,300000)); abline(v=175, col="green"); mtext(side = 2, text = "1D", line = 0.5)
plot(agg7[agg7$Group.1=="2A","Group.2"],agg7[agg7$Group.1=="2A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,300000)); abline(v=335, col="green"); mtext(side = 2, text = "2A", line = 0.5)
plot(agg7[agg7$Group.1=="2B","Group.2"],agg7[agg7$Group.1=="2B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,300000)); abline(v=360, col="green"); mtext(side = 2, text = "2B", line = 0.5)
plot(agg7[agg7$Group.1=="2D","Group.2"],agg7[agg7$Group.1=="2D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,300000)); abline(v=275, col="green"); mtext(side = 2, text = "2D", line = 0.5)
plot(agg7[agg7$Group.1=="3A","Group.2"],agg7[agg7$Group.1=="3A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,300000)); abline(v=330, col="green");  mtext(side = 2, text = "3A", line = 0.5)
plot(agg7[agg7$Group.1=="3B","Group.2"],agg7[agg7$Group.1=="3B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,300000)); abline(v=355, col="green"); mtext(side = 2, text = "3B", line = 0.5)
plot(agg7[agg7$Group.1=="3D","Group.2"],agg7[agg7$Group.1=="3D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,300000)); abline(v=245, col="green"); mtext(side = 2, text = "3D", line = 0.5)
plot(agg7[agg7$Group.1=="4A","Group.2"],agg7[agg7$Group.1=="4A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,300000)); abline(v=250, col="green"); mtext(side = 2, text = "4A", line = 0.5)
plot(agg7[agg7$Group.1=="4B","Group.2"],agg7[agg7$Group.1=="4B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,300000)); abline(v=300, col="green"); mtext(side = 2, text = "4B", line = 0.5)
plot(agg7[agg7$Group.1=="4D","Group.2"],agg7[agg7$Group.1=="4D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,300000)); abline(v=205, col="green");mtext(side = 2, text = "4D", line = 0.5)
plot(agg7[agg7$Group.1=="5A","Group.2"],agg7[agg7$Group.1=="5A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,300000)); abline(v=253, col="green"); mtext(side = 2, text = "5A", line = 0.5)
plot(agg7[agg7$Group.1=="5B","Group.2"],agg7[agg7$Group.1=="5B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,300000)); abline(v=200, col="green"); mtext(side = 2, text = "5B", line = 0.5)
plot(agg7[agg7$Group.1=="5D","Group.2"],agg7[agg7$Group.1=="5D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,300000)); abline(v=190, col="green"); mtext(side = 2, text = "5D", line = 0.5)
plot(agg7[agg7$Group.1=="6A","Group.2"],agg7[agg7$Group.1=="6A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,300000)); abline(v=290, col="green"); mtext(side = 2, text = "6A", line = 0.5)
plot(agg7[agg7$Group.1=="6B","Group.2"],agg7[agg7$Group.1=="6B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,300000)); abline(v=333, col="green"); mtext(side = 2, text = "6B", line = 0.5)
plot(agg7[agg7$Group.1=="6D","Group.2"],agg7[agg7$Group.1=="6D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,300000)); abline(v=235, col="green"); mtext(side = 2, text = "6D", line = 0.5)
plot(agg7[agg7$Group.1=="7A","Group.2"],agg7[agg7$Group.1=="7A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,300000)); abline(v=368, col="green"); mtext(side = 2, text = "7A", line = 0.5)
plot(agg7[agg7$Group.1=="7B","Group.2"],agg7[agg7$Group.1=="7B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,300000)); abline(v=305, col="green"); mtext(side = 2, text = "7B", line = 0.5)
plot(agg7[agg7$Group.1=="7D","Group.2"],agg7[agg7$Group.1=="7D","x"],type="l",  yaxt="n", xlim=c(0,900), ylim=c(0,300000)); abline(v=345, col="green"); mtext(side = 2, text = "7D", line = 0.5)
mtext(side = 1, text = "Chromosome Position(1Mbp)", line = 2.5)
gene_RLC <- merge(aggm, agg7, by=c(1,2))
cor(gene_RLC$x.x, gene_RLC$x.y)
#0.2099358

#now combination of LTR retrotransposons : RLC+RLG+RLX
#supplementary Figure 3
stat2<- stat1[stat1$V5.y%in% c("RLC","RLG","RLX"),]
agg7 <- aggregate(stat2$V7, by=list(stat2$V1, stat2$million), FUN=sum, na.rm=TRUE)
max(agg7$x)
#915719
#528,285 complete retrotrasposon
ymax <- 1000000
dev.off();
plot.new()

par(mar=c(0,2,0,0), mfrow = c(23,1))
plot(agg7[agg7$Group.1=="1A","Group.2"],agg7[agg7$Group.1=="1A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=213, col="green"); mtext(side = 2, text = "1A", line = 0.5)
plot(agg7[agg7$Group.1=="1B","Group.2"],agg7[agg7$Group.1=="1B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=240, col="green"); mtext(side = 2, text = "1B", line = 0.5)
plot(agg7[agg7$Group.1=="1D","Group.2"],agg7[agg7$Group.1=="1D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=175, col="green"); mtext(side = 2, text = "1D", line = 0.5)
plot(agg7[agg7$Group.1=="2A","Group.2"],agg7[agg7$Group.1=="2A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=335, col="green"); mtext(side = 2, text = "2A", line = 0.5)
plot(agg7[agg7$Group.1=="2B","Group.2"],agg7[agg7$Group.1=="2B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=358, col="green"); mtext(side = 2, text = "2B", line = 0.5)
plot(agg7[agg7$Group.1=="2D","Group.2"],agg7[agg7$Group.1=="2D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=273, col="green"); mtext(side = 2, text = "2D", line = 0.5)
plot(agg7[agg7$Group.1=="3A","Group.2"],agg7[agg7$Group.1=="3A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=330, col="green");  mtext(side = 2, text = "3A", line = 0.5)
plot(agg7[agg7$Group.1=="3B","Group.2"],agg7[agg7$Group.1=="3B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=355, col="green"); mtext(side = 2, text = "3B", line = 0.5)
plot(agg7[agg7$Group.1=="3D","Group.2"],agg7[agg7$Group.1=="3D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=245, col="green"); mtext(side = 2, text = "3D", line = 0.5)
plot(agg7[agg7$Group.1=="4A","Group.2"],agg7[agg7$Group.1=="4A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=313, col="green"); mtext(side = 2, text = "4A", line = 0.5)
plot(agg7[agg7$Group.1=="4B","Group.2"],agg7[agg7$Group.1=="4B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=300, col="green"); mtext(side = 2, text = "4B", line = 0.5)
plot(agg7[agg7$Group.1=="4D","Group.2"],agg7[agg7$Group.1=="4D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=205, col="green");mtext(side = 2, text = "4D", line = 0.5)
plot(agg7[agg7$Group.1=="5A","Group.2"],agg7[agg7$Group.1=="5A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=253, col="green"); mtext(side = 2, text = "5A", line = 0.5)
plot(agg7[agg7$Group.1=="5B","Group.2"],agg7[agg7$Group.1=="5B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=200, col="green"); mtext(side = 2, text = "5B", line = 0.5)
plot(agg7[agg7$Group.1=="5D","Group.2"],agg7[agg7$Group.1=="5D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=190, col="green"); mtext(side = 2, text = "5D", line = 0.5)
plot(agg7[agg7$Group.1=="6A","Group.2"],agg7[agg7$Group.1=="6A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=290, col="green"); mtext(side = 2, text = "6A", line = 0.5)
plot(agg7[agg7$Group.1=="6B","Group.2"],agg7[agg7$Group.1=="6B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=333, col="green"); mtext(side = 2, text = "6B", line = 0.5)
plot(agg7[agg7$Group.1=="6D","Group.2"],agg7[agg7$Group.1=="6D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=235, col="green"); mtext(side = 2, text = "6D", line = 0.5)
plot(agg7[agg7$Group.1=="7A","Group.2"],agg7[agg7$Group.1=="7A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=365, col="green"); mtext(side = 2, text = "7A", line = 0.5)
plot(agg7[agg7$Group.1=="7B","Group.2"],agg7[agg7$Group.1=="7B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=310, col="green"); mtext(side = 2, text = "7B", line = 0.5)
plot(agg7[agg7$Group.1=="7D","Group.2"],agg7[agg7$Group.1=="7D","x"],type="l",  yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=340, col="green"); mtext(side = 2, text = "7D", line = 0.5)
mtext(side = 1, text = "Chromosome Position(1Mbp)", line = 2.5)
gene_RL <- merge(aggm, agg7, by=c(1,2))
cor(gene_RL$x.x, gene_RL$x.y)
#-0.1352366

#Supplementary Figure 4
stat2<- stat1[stat1$V5.y=="DTC",]
agg8 <- aggregate(stat2$V7, by=list(stat2$V1, stat2$million), FUN=sum, na.rm=TRUE)
max(agg8$x)
#369579
#156,602 complete CACTA transposon
agg8[agg8$x==369579,]
# 4D     411 369579
ymax <- 400000
#ymax <- 200000
dev.off();
plot.new()
par(mar=c(0,2,0,0), mfrow = c(23,1))
plot(agg8[agg8$Group.1=="1A","Group.2"],agg8[agg8$Group.1=="1A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=215, col="green"); mtext(side = 2, text = "1A", line = 0.5)
plot(agg8[agg8$Group.1=="1B","Group.2"],agg8[agg8$Group.1=="1B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=240, col="green"); mtext(side = 2, text = "1B", line = 0.5)
plot(agg8[agg8$Group.1=="1D","Group.2"],agg8[agg8$Group.1=="1D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=175, col="green"); mtext(side = 2, text = "1D", line = 0.5)
plot(agg8[agg8$Group.1=="2A","Group.2"],agg8[agg8$Group.1=="2A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=335, col="green"); mtext(side = 2, text = "2A", line = 0.5)
plot(agg8[agg8$Group.1=="2B","Group.2"],agg8[agg8$Group.1=="2B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=360, col="green"); mtext(side = 2, text = "2B", line = 0.5)
plot(agg8[agg8$Group.1=="2D","Group.2"],agg8[agg8$Group.1=="2D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=275, col="green"); mtext(side = 2, text = "2D", line = 0.5)
plot(agg8[agg8$Group.1=="3A","Group.2"],agg8[agg8$Group.1=="3A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=330, col="green");  mtext(side = 2, text = "3A", line = 0.5)
plot(agg8[agg8$Group.1=="3B","Group.2"],agg8[agg8$Group.1=="3B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=355, col="green"); mtext(side = 2, text = "3B", line = 0.5)
plot(agg8[agg8$Group.1=="3D","Group.2"],agg8[agg8$Group.1=="3D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=245, col="green"); mtext(side = 2, text = "3D", line = 0.5)
plot(agg8[agg8$Group.1=="4A","Group.2"],agg8[agg8$Group.1=="4A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=250, col="green"); mtext(side = 2, text = "4A", line = 0.5)
plot(agg8[agg8$Group.1=="4B","Group.2"],agg8[agg8$Group.1=="4B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=300, col="green"); mtext(side = 2, text = "4B", line = 0.5)
plot(agg8[agg8$Group.1=="4D","Group.2"],agg8[agg8$Group.1=="4D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=205, col="green");mtext(side = 2, text = "4D", line = 0.5)
plot(agg8[agg8$Group.1=="5A","Group.2"],agg8[agg8$Group.1=="5A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=253, col="green"); mtext(side = 2, text = "5A", line = 0.5)
plot(agg8[agg8$Group.1=="5B","Group.2"],agg8[agg8$Group.1=="5B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=200, col="green"); mtext(side = 2, text = "5B", line = 0.5)
plot(agg8[agg8$Group.1=="5D","Group.2"],agg8[agg8$Group.1=="5D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=190, col="green"); mtext(side = 2, text = "5D", line = 0.5)
plot(agg8[agg8$Group.1=="6A","Group.2"],agg8[agg8$Group.1=="6A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=290, col="green"); mtext(side = 2, text = "6A", line = 0.5)
plot(agg8[agg8$Group.1=="6B","Group.2"],agg8[agg8$Group.1=="6B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=333, col="green"); mtext(side = 2, text = "6B", line = 0.5)
plot(agg8[agg8$Group.1=="6D","Group.2"],agg8[agg8$Group.1=="6D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=235, col="green"); mtext(side = 2, text = "6D", line = 0.5)
plot(agg8[agg8$Group.1=="7A","Group.2"],agg8[agg8$Group.1=="7A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=368, col="green"); mtext(side = 2, text = "7A", line = 0.5)
plot(agg8[agg8$Group.1=="7B","Group.2"],agg8[agg8$Group.1=="7B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=305, col="green"); mtext(side = 2, text = "7B", line = 0.5)
plot(agg8[agg8$Group.1=="7D","Group.2"],agg8[agg8$Group.1=="7D","x"],type="l",  yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=345, col="green"); mtext(side = 2, text = "7D", line = 0.5)
mtext(side = 1, text = "Chromosome Position(1Mbp)", line = 2.5)
gene_DTC <- merge(aggm, agg8, by=c(1,2))
cor(gene_DTC$x.x, gene_DTC$x.y)
#-0.1246951


stat2<- stat1[stat1$V5.y=="RIX",]
agg9 <- aggregate(stat2$V7, by=list(stat2$V1, stat2$million), FUN=sum, na.rm=TRUE)
max(agg9$x)
#53188
agg9[agg9$x==53188,]
# 2B      26 53188
ymax <- 100000

dev.off();
plot.new()
par(mar=c(0,2,0,0), mfrow = c(23,1))
plot(agg9[agg9$Group.1=="1A","Group.2"],agg9[agg9$Group.1=="1A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=215, col="green"); mtext(side = 2, text = "1A", line = 0.5)
plot(agg9[agg9$Group.1=="1B","Group.2"],agg9[agg9$Group.1=="1B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=240, col="green"); mtext(side = 2, text = "1B", line = 0.5)
plot(agg9[agg9$Group.1=="1D","Group.2"],agg9[agg9$Group.1=="1D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=175, col="green"); mtext(side = 2, text = "1D", line = 0.5)
plot(agg9[agg9$Group.1=="2A","Group.2"],agg9[agg9$Group.1=="2A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=335, col="green"); mtext(side = 2, text = "2A", line = 0.5)
plot(agg9[agg9$Group.1=="2B","Group.2"],agg9[agg9$Group.1=="2B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=360, col="green"); mtext(side = 2, text = "2B", line = 0.5)
plot(agg9[agg9$Group.1=="2D","Group.2"],agg9[agg9$Group.1=="2D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=275, col="green"); mtext(side = 2, text = "2D", line = 0.5)
plot(agg9[agg9$Group.1=="3A","Group.2"],agg9[agg9$Group.1=="3A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=330, col="green");  mtext(side = 2, text = "3A", line = 0.5)
plot(agg9[agg9$Group.1=="3B","Group.2"],agg9[agg9$Group.1=="3B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=355, col="green"); mtext(side = 2, text = "3B", line = 0.5)
plot(agg9[agg9$Group.1=="3D","Group.2"],agg9[agg9$Group.1=="3D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=245, col="green"); mtext(side = 2, text = "3D", line = 0.5)
plot(agg9[agg9$Group.1=="4A","Group.2"],agg9[agg9$Group.1=="4A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=250, col="green"); mtext(side = 2, text = "4A", line = 0.5)
plot(agg9[agg9$Group.1=="4B","Group.2"],agg9[agg9$Group.1=="4B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=300, col="green"); mtext(side = 2, text = "4B", line = 0.5)
plot(agg9[agg9$Group.1=="4D","Group.2"],agg9[agg9$Group.1=="4D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=205, col="green");mtext(side = 2, text = "4D", line = 0.5)
plot(agg9[agg9$Group.1=="5A","Group.2"],agg9[agg9$Group.1=="5A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=253, col="green"); mtext(side = 2, text = "5A", line = 0.5)
plot(agg9[agg9$Group.1=="5B","Group.2"],agg9[agg9$Group.1=="5B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=200, col="green"); mtext(side = 2, text = "5B", line = 0.5)
plot(agg9[agg9$Group.1=="5D","Group.2"],agg9[agg9$Group.1=="5D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=190, col="green"); mtext(side = 2, text = "5D", line = 0.5)
plot(agg9[agg9$Group.1=="6A","Group.2"],agg9[agg9$Group.1=="6A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=290, col="green"); mtext(side = 2, text = "6A", line = 0.5)
plot(agg9[agg9$Group.1=="6B","Group.2"],agg9[agg9$Group.1=="6B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=333, col="green"); mtext(side = 2, text = "6B", line = 0.5)
plot(agg9[agg9$Group.1=="6D","Group.2"],agg9[agg9$Group.1=="6D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=235, col="green"); mtext(side = 2, text = "6D", line = 0.5)
plot(agg9[agg9$Group.1=="7A","Group.2"],agg9[agg9$Group.1=="7A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=368, col="green"); mtext(side = 2, text = "7A", line = 0.5)
plot(agg9[agg9$Group.1=="7B","Group.2"],agg9[agg9$Group.1=="7B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=305, col="green"); mtext(side = 2, text = "7B", line = 0.5)
plot(agg9[agg9$Group.1=="7D","Group.2"],agg9[agg9$Group.1=="7D","x"],type="l",  yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=345, col="green"); mtext(side = 2, text = "7D", line = 0.5)
mtext(side = 1, text = "Chromosome Position(1Mbp)", line = 2.5)

gene_RIX <- merge(aggm, agg9, by=c(1,2))
cor(gene_RIX$x.x, gene_RIX$x.y)
#0.2469125

stat2<- stat1[stat1$V5.y=="SIX",]
agg10 <- aggregate(stat2$V7, by=list(stat2$V1, stat2$million), FUN=sum, na.rm=TRUE)
gene_SIX <- merge(aggm, agg10, by=c(1,2))
cor(gene_SIX$x.x, gene_SIX$x.y)
#0.05328177

stat2<- stat1[stat1$V5.y=="DHH",]
agg11 <- aggregate(stat2$V7, by=list(stat2$V1, stat2$million), FUN=sum, na.rm=TRUE)
gene_DHH <- merge(aggm, agg11, by=c(1,2))
cor(gene_DHH$x.x, gene_DHH$x.y)
#0.030

stat2<- stat1[stat1$V5.y=="DTM",]
agg12 <- aggregate(stat2$V7, by=list(stat2$V1, stat2$million), FUN=sum, na.rm=TRUE)
gene_DTM <- merge(aggm, agg12, by=c(1,2))
cor(gene_DTM$x.x, gene_DTM$x.y)
#0.1838

stat2<- stat1[stat1$V5.y=="DTA",]
agg13 <- aggregate(stat2$V7, by=list(stat2$V1, stat2$million), FUN=sum, na.rm=TRUE)
gene_DTA <- merge(aggm, agg13, by=c(1,2))
cor(gene_DTA$x.x, gene_DTA$x.y)
#-0.05553094

stat2<- stat1[stat1$V5.y=="DTH",]
agg14 <- aggregate(stat2$V7, by=list(stat2$V1, stat2$million), FUN=sum, na.rm=TRUE)
gene_DTH <- merge(aggm, agg14, by=c(1,2))
cor(gene_DTH$x.x, gene_DTH$x.y)
#0.288926
stat2<- stat1[stat1$V5.y=="DTT",]
agg15 <- aggregate(stat2$V7, by=list(stat2$V1, stat2$million), FUN=sum, na.rm=TRUE)
max(agg15$x)
#10,159
ymax <- 10000
dev.off();
plot.new()
par(mar=c(0,2,0,0), mfrow = c(23,1))
plot(agg15[agg15$Group.1=="1A","Group.2"],agg15[agg15$Group.1=="1A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=215, col="green"); mtext(side = 2, text = "1A", line = 0.5)
plot(agg15[agg15$Group.1=="1B","Group.2"],agg15[agg15$Group.1=="1B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=240, col="green"); mtext(side = 2, text = "1B", line = 0.5)
plot(agg15[agg15$Group.1=="1D","Group.2"],agg15[agg15$Group.1=="1D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=175, col="green"); mtext(side = 2, text = "1D", line = 0.5)
plot(agg15[agg15$Group.1=="2A","Group.2"],agg15[agg15$Group.1=="2A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=335, col="green"); mtext(side = 2, text = "2A", line = 0.5)
plot(agg15[agg15$Group.1=="2B","Group.2"],agg15[agg15$Group.1=="2B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=360, col="green"); mtext(side = 2, text = "2B", line = 0.5)
plot(agg15[agg15$Group.1=="2D","Group.2"],agg15[agg15$Group.1=="2D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=275, col="green"); mtext(side = 2, text = "2D", line = 0.5)
plot(agg15[agg15$Group.1=="3A","Group.2"],agg15[agg15$Group.1=="3A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=330, col="green");  mtext(side = 2, text = "3A", line = 0.5)
plot(agg15[agg15$Group.1=="3B","Group.2"],agg15[agg15$Group.1=="3B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=355, col="green"); mtext(side = 2, text = "3B", line = 0.5)
plot(agg15[agg15$Group.1=="3D","Group.2"],agg15[agg15$Group.1=="3D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=245, col="green"); mtext(side = 2, text = "3D", line = 0.5)
plot(agg15[agg15$Group.1=="4A","Group.2"],agg15[agg15$Group.1=="4A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=250, col="green"); mtext(side = 2, text = "4A", line = 0.5)
plot(agg15[agg15$Group.1=="4B","Group.2"],agg15[agg15$Group.1=="4B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=300, col="green"); mtext(side = 2, text = "4B", line = 0.5)
plot(agg15[agg15$Group.1=="4D","Group.2"],agg15[agg15$Group.1=="4D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=205, col="green");mtext(side = 2, text = "4D", line = 0.5)
plot(agg15[agg15$Group.1=="5A","Group.2"],agg15[agg15$Group.1=="5A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=253, col="green"); mtext(side = 2, text = "5A", line = 0.5)
plot(agg15[agg15$Group.1=="5B","Group.2"],agg15[agg15$Group.1=="5B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=200, col="green"); mtext(side = 2, text = "5B", line = 0.5)
plot(agg15[agg15$Group.1=="5D","Group.2"],agg15[agg15$Group.1=="5D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=190, col="green"); mtext(side = 2, text = "5D", line = 0.5)
plot(agg15[agg15$Group.1=="6A","Group.2"],agg15[agg15$Group.1=="6A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=290, col="green"); mtext(side = 2, text = "6A", line = 0.5)
plot(agg15[agg15$Group.1=="6B","Group.2"],agg15[agg15$Group.1=="6B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=333, col="green"); mtext(side = 2, text = "6B", line = 0.5)
plot(agg15[agg15$Group.1=="6D","Group.2"],agg15[agg15$Group.1=="6D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=235, col="green"); mtext(side = 2, text = "6D", line = 0.5)
plot(agg15[agg15$Group.1=="7A","Group.2"],agg15[agg15$Group.1=="7A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=368, col="green"); mtext(side = 2, text = "7A", line = 0.5)
plot(agg15[agg15$Group.1=="7B","Group.2"],agg15[agg15$Group.1=="7B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=305, col="green"); mtext(side = 2, text = "7B", line = 0.5)
plot(agg15[agg15$Group.1=="7D","Group.2"],agg15[agg15$Group.1=="7D","x"],type="l",  yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=345, col="green"); mtext(side = 2, text = "7D", line = 0.5)
mtext(side = 1, text = "Chromosome Position(1Mbp)", line = 2.5)




gene_DTT <- merge(aggm, agg15, by=c(1,2))
cor(gene_DTT$x.x, gene_DTT$x.y)
#0.3028991

#now sum of DNA transposan
stat2<- stat1[stat1$V5.y%in%c("DTM","DTC","DTA","DTH","DTT","DTX", "DHH"),]
agg16 <- aggregate(stat2$V7, by=list(stat2$V1, stat2$million), FUN=sum, na.rm=TRUE)
max(agg16$x)
#367997
ymax <- 400000
dev.off();
plot.new()
par(mar=c(0,2,0,0), mfrow = c(23,1))
plot(agg16[agg16$Group.1=="1A","Group.2"],agg16[agg16$Group.1=="1A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=215, col="green"); mtext(side = 2, text = "1A", line = 0.5)
plot(agg16[agg16$Group.1=="1B","Group.2"],agg16[agg16$Group.1=="1B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=240, col="green"); mtext(side = 2, text = "1B", line = 0.5)
plot(agg16[agg16$Group.1=="1D","Group.2"],agg16[agg16$Group.1=="1D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=175, col="green"); mtext(side = 2, text = "1D", line = 0.5)
plot(agg16[agg16$Group.1=="2A","Group.2"],agg16[agg16$Group.1=="2A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=335, col="green"); mtext(side = 2, text = "2A", line = 0.5)
plot(agg16[agg16$Group.1=="2B","Group.2"],agg16[agg16$Group.1=="2B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=360, col="green"); mtext(side = 2, text = "2B", line = 0.5)
plot(agg16[agg16$Group.1=="2D","Group.2"],agg16[agg16$Group.1=="2D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=275, col="green"); mtext(side = 2, text = "2D", line = 0.5)
plot(agg16[agg16$Group.1=="3A","Group.2"],agg16[agg16$Group.1=="3A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=330, col="green");  mtext(side = 2, text = "3A", line = 0.5)
plot(agg16[agg16$Group.1=="3B","Group.2"],agg16[agg16$Group.1=="3B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=355, col="green"); mtext(side = 2, text = "3B", line = 0.5)
plot(agg16[agg16$Group.1=="3D","Group.2"],agg16[agg16$Group.1=="3D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=245, col="green"); mtext(side = 2, text = "3D", line = 0.5)
plot(agg16[agg16$Group.1=="4A","Group.2"],agg16[agg16$Group.1=="4A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=250, col="green"); mtext(side = 2, text = "4A", line = 0.5)
plot(agg16[agg16$Group.1=="4B","Group.2"],agg16[agg16$Group.1=="4B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=300, col="green"); mtext(side = 2, text = "4B", line = 0.5)
plot(agg16[agg16$Group.1=="4D","Group.2"],agg16[agg16$Group.1=="4D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=205, col="green");mtext(side = 2, text = "4D", line = 0.5)
plot(agg16[agg16$Group.1=="5A","Group.2"],agg16[agg16$Group.1=="5A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=253, col="green"); mtext(side = 2, text = "5A", line = 0.5)
plot(agg16[agg16$Group.1=="5B","Group.2"],agg16[agg16$Group.1=="5B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=200, col="green"); mtext(side = 2, text = "5B", line = 0.5)
plot(agg16[agg16$Group.1=="5D","Group.2"],agg16[agg16$Group.1=="5D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=190, col="green"); mtext(side = 2, text = "5D", line = 0.5)
plot(agg16[agg16$Group.1=="6A","Group.2"],agg16[agg16$Group.1=="6A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=290, col="green"); mtext(side = 2, text = "6A", line = 0.5)
plot(agg16[agg16$Group.1=="6B","Group.2"],agg16[agg16$Group.1=="6B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=333, col="green"); mtext(side = 2, text = "6B", line = 0.5)
plot(agg16[agg16$Group.1=="6D","Group.2"],agg16[agg16$Group.1=="6D","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=235, col="green"); mtext(side = 2, text = "6D", line = 0.5)
plot(agg16[agg16$Group.1=="7A","Group.2"],agg16[agg16$Group.1=="7A","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=368, col="green"); mtext(side = 2, text = "7A", line = 0.5)
plot(agg16[agg16$Group.1=="7B","Group.2"],agg16[agg16$Group.1=="7B","x"],type="l", xaxt="n", yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=305, col="green"); mtext(side = 2, text = "7B", line = 0.5)
plot(agg16[agg16$Group.1=="7D","Group.2"],agg16[agg16$Group.1=="7D","x"],type="l",  yaxt="n", xlim=c(0,900), ylim=c(0,ymax)); abline(v=345, col="green"); mtext(side = 2, text = "7D", line = 0.5)
mtext(side = 1, text = "Chromosome Position(1Mbp)", line = 2.5)


gene_DT <- merge(aggm, agg16, by=c(1,2))
cor(gene_DT$x.x, gene_DT$x.y)
#-0.0232
###############################################################
#below code is for generating Figure 4d(heatmap) of manuscript
dfcnv <- read.table("cnv_table.txt", header=F, sep="\t")
rownames(dfcnv) <- dfcnv$V1
dfcnv <- dfcnv[,-1]
#dfcnv <- dfcnv[,-ncol(dfcnv)]
colnames(dfcnv) <- c("1A","1B","1D","2A","2B","2D","3A","3B","3D","4A","4B","4D","5A","5B","5D","6A","6B","6D","7A","7B","7D")
dfcnv[,] <- lapply(dfcnv[,], as.character)
dfcnv[,] <- lapply(dfcnv[,], strsplit, ":", fixed=T)
cnv_vec <- sapply(as.vector(as.matrix(dfcnv)), "[", 2)
cnv_mat <- matrix(cnv_vec, nrow = nrow(dfcnv), ncol = ncol(dfcnv))
dfcnv2 <- as.data.frame(cnv_mat)
rownames(dfcnv2) <- rownames(dfcnv) 
colnames(dfcnv2) <- colnames(dfcnv)
dfcnv2[,] <- lapply(dfcnv2[,], as.numeric)
heatmap(as.matrix(dfcnv2), scale='none' )
points(0.26,0.9, col = "black", pch = 19, cex = 1)#top left
points(0.71,0, col = "black", pch = 7, cex = 1)#bottom   right
points(0.28,0.9, col = "black", pch = 7, cex = 1)#bottom   right
points(0.42,0.9, col = "black", pch = 17, cex = 1)

dev.off();
###############################################################

#bedtools intersect -u -f 0.9 -a  iwgsc_refseqv2.1_annotation_200916_HC+.bed -b cnv_700_sorted.bed > genes_cnv_90.bed
dfcnv <- read.table("genes_cnv.bed", sep="\t", header=F)
dim(dfcnv)
#591   3
dfcnv$cc <- rep(nrow(dfcnv),1)
agg <- aggregate(cc~V1+V2+V3, dfcnv, length)
dim(agg)
#584 4
df <- read.table("iwgsc_refseqv2.1_annotation_200916_HC.gen", sep="\t", header=F)
write.table(df[,c(2,3,4,1)], file="iwgsc_refseqv2.1_annotation_200916_HC.txt", sep="\t", quote = F,  row.names = F, col.names=F)
#copy number variations intersect with TE
gcnvte <- read.table("gene_cnv_te.txt", sep="\t", header=F)
dim(gcnvte)
#646     12
length(which(gcnvte$V12==0))
#495 genes have no intersection with TE, we extract 1k flanking region
gcnvte_zero <- gcnvte[gcnvte$V12==0,]
lst <- list(); 
for(i in 1:nrow(gcnvte_zero)){
  j <- (i*2) -1
  lst[[j]] <- list(gcnvte_zero[i,1],gcnvte_zero[i,2]-1000,gcnvte_zero[i,2], gcnvte_zero[i,4])
  lst[[j+1]] <- list(gcnvte_zero[i,1],gcnvte_zero[i,3],gcnvte_zero[i,3]+1000, gcnvte_zero[i,4])
}
df <- do.call("rbind", lst)
df<-as.data.frame(df)
df$V1 <- as.character(df$V1)
df$V4 <- as.character(df$V4)
df$V2 <- as.numeric(df$V2)
df$V3 <- as.numeric(df$V3)
write.table(df, file="full_genes_flank.txt", sep="\t", quote = F,  row.names = F, col.names=F)

gcnvte <- gcnvte[gcnvte$V12>0,]
dim(gcnvte)
#151   12
gcnvte$cc <- rep(1, nrow(gcnvte))
agg <- aggregate(cc~., gcnvte[,c(10, ncol(gcnvte))], length)
agg
'
1  DTC 48
2  DTH  3
3  DTM  9
4  DTT 38
5  DTX 14
6  RIX 11
7  RLC  8
8  RLG  8
9  RLX  5
10 XXX  7
'
sum(agg$cc)-7
#151-7
aggregate(cc~., gcnvte[gcnvte$V9=="complete",c(8, ncol(gcnvte))], length)
'V8 cc
1 DTC  7
2 DTH  1
3 DTM  3
4 DTT 29
5 DTX  7
6 RLC  2
7 RLX  1'
#Now investigate promoter of all genes
promte <- read.table("promoter_te.txt", sep="\t", header=F)
promte <- promte[promte$V11>0,]
promte$cc <- rep(1, nrow(promte))
agg <- aggregate(cc~., promte[,c(9, ncol(promte))], length)
agg
'
    V9    cc
1  DHH    27
2  DTA   102
3  DTC 30763
4  DTH  6748
5  DTM  9746
6  DTT 16483
7  DTX 12299
8  DXX   649
9  RIX 12319
10 RLC 15915
11 RLG 21020
12 RLX  5084
13 SIX   490
14 XXX 25813
'
cnvpromte <- read.table("591_cnv_gene_promoter_te.txt", sep="\t", header=F)
cnvpromte <- cnvpromte[cnvpromte$V11>0,]
cnvpromte$cc <- rep(1, nrow(cnvpromte))
#see how many genes promoters have intersection with TE
agg <- aggregate(cc~., cnvpromte[,c(4, ncol(cnvpromte))], length)
dim(agg)
#412
agg <- aggregate(cc~., cnvpromte[,c(9, ncol(cnvpromte))], length)
agg
'
1  DTC 177
2  DTH  34
3  DTM  61
4  DTT  71
5  DTX  65
6  DXX   1
7  RIX  47
8  RLC  66
9  RLG  80
10 RLX  22
11 SIX   2
12 XXX 118
'
sum(agg$cc)
#744
###############################################

dfcnvte <- read.table("cnv_te.bed", sep="\t", header=F)
dim(dfcnvte)
#113,313     10

# Now perform the same analysis for 591 genes with cnv
dfcnv <- read.table("cnv_sorted.bed", sep="\t", header=F)
dfcnv$ID <- paste("CNV",dfcnv$V1, dfcnv$V2, sep="-")
dfte <- read.table("Tae.Chinese_Spring.refSeqv2.1.TE_annotation.txt", sep="\t", header=F)
#now see how length of TE are based on Type
dfte$len <- dfte$V3 - dfte$V2 + 1
aggregate(len~V5, dfte, median)
'
1  DHH  299.0
2  DTA  330.0
3  DTC  477.0
4  DTH  230.5
5  DTM  233.0
6  DTT   95.0
7  DTX  154.0
8  DXX  449.0
9  RIX  340.0
10 RLC 2092.0
11 RLG 2409.0
12 RLX 1236.0
13 SIX  172.0
14 XXX  234.0
'
#generally retrotransposon are longer than DNA transposon

colnames(dfte)[4] <- "ID"  # make column names equal
dfcnv <- rbind(dfcnv, dfte[, c(1,2,3,4)])
#now double rows of dfcnv
lst <- list(); 
for(i in 1:nrow(dfcnv)){
  j <- (i*2) -1
  lst[[j]] <- list(dfcnv[i,1],dfcnv[i,2],dfcnv[i,4])
  lst[[j+1]] <- list(dfcnv[i,1],dfcnv[i,3],dfcnv[i,4])
}
df <- do.call("rbind", lst)
df<-as.data.frame(df)
df$V1 <- as.character(df$V1)
df$V3 <- as.character(df$V3)
df$V2 <- as.numeric(df$V2)
dim(df)
# 9617936       3
class(df)
#now sort based on position
df <- df[order(df$V1,df$V2),]

#Now do the same logic on dfte to see if TE elements are overlapped
lst <- list(); 
for(i in 1:nrow(dfte)){
  j <- (i*2) -1
  lst[[j]] <- list(dfte[i,1],dfte[i,2],dfte[i,4])
  lst[[j+1]] <- list(dfte[i,1],dfte[i,3],dfte[i,4])
}

dfte_flat <- do.call("rbind", lst)
dfte_flat<-as.data.frame(dfte_flat)
dfte_flat$V1 <- as.character(dfte_flat$V1)
dfte_flat$V3 <- as.character(dfte_flat$V3)
dfte_flat$V2 <- as.numeric(dfte_flat$V2)
dfte_flat <- dfte_flat[order(dfte_flat$V1,dfte_flat$V2),]
dfte_flat$V4 <- c(dfte_flat$V3[2:nrow(dfte_flat)],"")
#now remove even rows
dfte_flat <- dfte_flat[seq(1,nrow(dfte_flat),2),]
all(dfte_flat$V3==dfte_flat$V4)
#FALSE
length(which(dfte_flat$V3!=dfte_flat$V4))
#1479  only 1479 TE elemnts have overlap with each other, and they are because of part_of that is not really an overlap
embedded <- dfte_flat[which(dfte_flat$V3!=dfte_flat$V4),3]
dfte_flat2 <- dfte_flat[dfte_flat$V3%in%embedded,]
#Now l

#now look at cnv regions only (first regions fully overlapped by TE)
dfcnvte <- read.table("cnv_1k_te_09.txt", sep="\t", header=F)
dfcnvte <- dfcnvte[dfcnvte$V10>0,]
#55,119     3
dfcnvte$cc <- rep(nrow(dfcnvte),1)
agg <- aggregate(cc~V1+V2+V3, dfcnvte, length)
dim(agg)
#55088     4
agg <- aggregate(cc~., dfcnvte[,c(8, ncol(dfcnvte))], length)
#below table goes to supplementary table 5
'
    V8    cc
1  DHH     1
2  DTA     1
3  DTC 13,674
4  DTH    12
5  DTM   151
6  DTT     3
7  DTX    43
8  RIX   452
9  RLC 13,215
10 RLG 24,092
11 RLX  3,424
12 XXX    51
'
sum(agg$cc)-51
#now look at cnv regions only that fully or partially overlapped with TE
dfcnvte <- read.table("cnv_1k_te.txt", sep="\t", header=F)
dfcnvte <- dfcnvte[dfcnvte$V10>0,]
# 109009     3
dfcnvte$cc <- rep(nrow(dfcnvte),1)
agg <- aggregate(cc~V1+V2+V3, dfcnvte, length)
dim(agg)
#73431     4

#now for cnv that map to multiple TE we pick the most frequent one
agg2 <- aggregate(V8~V1+V2+V3, dfcnvte , FUN=function(v){return(names(which.max(table(v))))})
dim(agg2)
#73431     4
agg2$cc <- rep(nrow(agg2),1)
agg3 <- aggregate(cc~V8, agg2[,c('V8','cc')], length)
agg3
'V8    cc
1  DHH    36
2  DTA    22
3  DTC 20387
4  DTH   216
5  DTM   951
6  DTT   623
7  DTX   472
8  DXX    60
9  RIX  1439
10 RLC 16106
11 RLG 28495
12 RLX  3736
13 SIX     8
14 XXX   880'
sum(agg3$cc)


#############################################################################################
#circlize coverage graphs Figure 5 of manuscript
#First we need to re-format dfcovg as: chr start end value
dim(dfcovg)
covlst <- list()
i <- 1
for(r in 1:nrow(dfcovg)){
  chr <- rownames(dfcovg)[r]
  for(c in 1:ncol(dfcovg)){
    start <- (c-1)*1000000 +1
    end <- c*1000000
    start <- as.integer(start)
    end <- as.integer(end)
    value <- dfcovg[r,c]
    covlst[[i]] <- list(chr, start, end, value)
    
    i <- i+1
  }
}
dffullcovg <- do.call("rbind", covlst)
dffullcovg <- as.data.frame(dffullcovg)
dim(dffullcovg)
#21000,4
colnames(dffullcovg) <- c("chr","start","end","value")
dffullcovg$chr <- as.character(dffullcovg$chr)
dffullcovg$start <- as.integer(dffullcovg$start)
dffullcovg$end <- as.integer(dffullcovg$end)
dffullcovg$value <- as.numeric(dffullcovg$value)

#now for intersection of geneic regions and core
dfgcovg <- aggm
dfgcovg$start <- (dfgcovg$Group.2 * 1000000) + 1
dfgcovg$end <- (dfgcovg$Group.2 + 1) * 1000000
dfgcovg <- dfgcovg[,c("Group.1","start","end","x")]
colnames(dfgcovg) <- c("chr","start","end","value")
dfgcovg$chr <- as.character(dfgcovg$chr)
dfgcovg$start <- as.integer(dfgcovg$start)
dfgcovg$end <- as.integer(dfgcovg$end)
dfgcovg$value <- as.numeric(dfgcovg$value)


#now for retrotransposons agg7
agg7$start <- (agg7$Group.2 * 1000000) + 1
agg7$end <- (agg7$Group.2 + 1) * 1000000
agg7 <- agg7[,c("Group.1","start","end","x")]
colnames(agg7) <- c("chr","start","end","value")
agg7$chr <- as.character(agg7$chr)
agg7$start <- as.integer(agg7$start)
agg7$end <- as.integer(agg7$end)
agg7$value <- as.numeric(agg7$value)

#now for CACTA 
agg8$start <- (agg8$Group.2 * 1000000) + 1
agg8$end <- (agg8$Group.2 + 1) * 1000000
agg8 <- agg8[,c("Group.1","start","end","x")]
colnames(agg8) <- c("chr","start","end","value")
agg8$chr <- as.character(agg8$chr)
agg8$start <- as.integer(agg8$start)
agg8$end <- as.integer(agg8$end)
agg8$value <- as.numeric(agg8$value)


max(as.numeric(dffullcovg[(dffullcovg$chr=='1A' & dffullcovg$value>0),3]))
#[1] 5.99e+08
max(as.numeric(dffullcovg[(dffullcovg$chr=='2A' & dffullcovg$value>0),3]))
#[1] 7.88e+08
max(as.numeric(dffullcovg[(dffullcovg$chr=='3A' & dffullcovg$value>0),3]))
#[1] 7.55e+08
max(as.numeric(dffullcovg[(dffullcovg$chr=='4A' & dffullcovg$value>0),3]))
#[1] 7.55e+08
max(as.numeric(dffullcovg[(dffullcovg$chr=='5A' & dffullcovg$value>0),3]))
#[1] 7.14e+08
max(as.numeric(dffullcovg[(dffullcovg$chr=='6A' & dffullcovg$value>0),3]))
#[1] 6.23e+08
max(as.numeric(dffullcovg[(dffullcovg$chr=='7A' & dffullcovg$value>0),3]))
#[1] 7.45e+08

max(as.numeric(dffullcovg[(dffullcovg$chr=='1B' & dffullcovg$value>0),3]))
#[1]  7.01e+08
max(as.numeric(dffullcovg[(dffullcovg$chr=='2B' & dffullcovg$value>0),3]))
#[1] 8.13e+08
max(as.numeric(dffullcovg[(dffullcovg$chr=='3B' & dffullcovg$value>0),3]))
#[1] 8.52e+08
max(as.numeric(dffullcovg[(dffullcovg$chr=='4B' & dffullcovg$value>0),3]))
#[1] 6.74e+08
max(as.numeric(dffullcovg[(dffullcovg$chr=='5B' & dffullcovg$value>0),3]))
#[1] 7.15e+08
max(as.numeric(dffullcovg[(dffullcovg$chr=='6B' & dffullcovg$value>0),3]))
#[1] 7.32e+08
max(as.numeric(dffullcovg[(dffullcovg$chr=='7B' & dffullcovg$value>0),3]))
#[1] 7.64e+08
max(as.numeric(dffullcovg[(dffullcovg$chr=='1D' & dffullcovg$value>0),3]))
#[1]  4.99e+08
max(as.numeric(dffullcovg[(dffullcovg$chr=='2D' & dffullcovg$value>0),3]))
#[1] 6.57e+08
max(as.numeric(dffullcovg[(dffullcovg$chr=='3D' & dffullcovg$value>0),3]))
#[1] 6.2e+08
max(as.numeric(dffullcovg[(dffullcovg$chr=='4D' & dffullcovg$value>0),3]))
#[1] 5.19e+08
max(as.numeric(dffullcovg[(dffullcovg$chr=='5D' & dffullcovg$value>0),3]))
#[1] 5.7e+08
max(as.numeric(dffullcovg[(dffullcovg$chr=='6D' & dffullcovg$value>0),3]))
#[1] 4.96e+08
max(as.numeric(dffullcovg[(dffullcovg$chr=='7D' & dffullcovg$value>0),3]))
#[1] 6.43e+08

library(circlize)
circos.clear()
circos.par("start.degree" = 90, gap.after = c(rep(2,20), 10))


#genomics application




circos.initializeWithIdeogram()
circos.genomicTrackPlotRegion(bed, panel.fun = function(region, value, ...) {
  circos.genomicPoints(region, value, ...)
})
class(bed[,5])

dfA_covg <- dffullcovg[dffullcovg$chr %in% c("1A","2A","3A","4A","5A","6A","7A"),]
dfA_covg <- dfA_covg[dfA_covg$value>0,]
dfA_covg$value <- dfA_covg$value/1000000;
dfA_covg$chr <- paste("chr", dfA_covg$chr, sep="")
max(dfA_covg$value)
#0.93627
dfgA_covg <- dfgcovg[dfgcovg$chr %in% c("1A","2A","3A","4A","5A","6A","7A"),]
dfgA_covg$value <- dfgA_covg$value/1000000;
dfgA_covg$chr <- paste("chr", dfgA_covg$chr, sep="")
max(dfgA_covg$value)
#0.131636
retroA_covg <- agg7[agg7$chr %in% c("1A","2A","3A","4A","5A","6A","7A"),]
retroA_covg$value <- retroA_covg$value/1000000;
retroA_covg$chr <- paste("chr", retroA_covg$chr, sep="")
max(retroA_covg$value)
#0.915719

dtcA_covg <- agg8[agg8$chr %in% c("1A","2A","3A","4A","5A","6A","7A"),]
dtcA_covg$value <- dtcA_covg$value/1000000;
dtcA_covg$chr <- paste("chr", dtcA_covg$chr, sep="")
max(dtcA_covg$value)
#0.279875

centromere_A <- c(213000000, 335000000, 330000000, 313000000, 253000000, 290000000, 365000000)
names(centromere_A) <- c("chr1A","chr2A","chr3A","chr4A","chr5A","chr6A","chr7A")
centromere_B <- c(24000000, 358000000, 355000000, 300000000, 200000000, 333000000, 310000000)
names(centromere_B) <- c("chr1B","chr2B","chr3B","chr4B","chr5B","chr6B","chr7B")
centromere_D <- c(175000000, 273000000, 245000000, 205000000, 190000000, 235000000, 340000000)
names(centromere_D) <- c("chr1D","chr2D","chr3D","chr4D","chr5D","chr6D","chr7D")
#repeat below code 3 times
genome <- "D"
chrs <- paste(c(1:7),genome,sep="")
chrs2 <- paste("chr",chrs, sep="")
if (genome=="A"){
  centromere <- c(213000000, 335000000, 330000000, 313000000, 253000000, 290000000, 365000000)
}else if (genome=="B"){
  centromere <- c(240000000, 358000000, 355000000, 300000000, 200000000, 333000000, 310000000)
}else{
  centromere <- c(175000000, 273000000, 245000000, 205000000, 190000000, 235000000, 340000000)
}
names(centromere) <- chrs2

dfA_covg <- dffullcovg[dffullcovg$chr %in% chrs,]
dfA_covg <- dfA_covg[dfA_covg$value>0,]
dfA_covg$value <- dfA_covg$value/1000000;
dfA_covg$chr <- paste("chr", dfA_covg$chr, sep="")
max(dfA_covg$value)
#0.936277
dfgA_covg <- dfgcovg[dfgcovg$chr %in% chrs,]
dfgA_covg$value <- dfgA_covg$value/1000000;
dfgA_covg$chr <- paste("chr", dfgA_covg$chr, sep="")
max(dfgA_covg$value)
#0.131636
retroA_covg <- agg7[agg7$chr %in% chrs,]
retroA_covg$value <- retroA_covg$value/1000000;
retroA_covg$chr <- paste("chr", retroA_covg$chr, sep="")
max(retroA_covg$value)
#0.915719

dtcA_covg <- agg8[agg8$chr %in% chrs,]
dtcA_covg$value <- dtcA_covg$value/1000000;
dtcA_covg$chr <- paste("chr", dtcA_covg$chr, sep="")
max(dtcA_covg$value)  
dev.off()
circos.clear()
circos.par("start.degree" = 90, "track.height" = 0.1)
#dim(cytoband)
tdf <- cytoband[which(cytoband[,1] %in% chrs2),]
circos.initializeWithIdeogram(tdf)
#class(tdf[,3])

circos.genomicTrack(dfA_covg, ylim = c(0, 1), panel.fun = function(region, value, ...) {
  circos.genomicLines(region, value,lwd=0.1, type='l', col="black")
  if(CELL_META$sector.index == "chr2B"){
    circos.rect(xleft = 100000000, ybottom = 0, xright = 527000000, ytop = 1, border="red")
  }
  x_centre = centromere[CELL_META$sector.index]
  circos.lines(c(x_centre, x_centre), c(0,1.2), col="green", lty=2, lwd=0.5)
  
})
circos.genomicTrack(retroA_covg, ylim = c(0, 1), panel.fun = function(region, value, ...) {
  circos.genomicLines(region, value,lwd=0.1, type='l', col="purple")
  x_centre = centromere[CELL_META$sector.index]
  circos.lines(c(x_centre, x_centre), c(0,1.2), col="green", lty=2, lwd=0.5)
  
})


circos.genomicTrack(dtcA_covg, ylim = c(0, 1), panel.fun = function(region, value, ...) {
  circos.genomicLines(region, value,lwd=0.1, type="l", col="orange")
  x_centre = centromere[CELL_META$sector.index]
  circos.lines(c(x_centre, x_centre), c(0,1.2), col="green", lty=2, lwd=0.5)
  
})
circos.genomicTrack(dfgA_covg, ylim = c(0, 0.15), panel.fun = function(region, value, ...) {
  circos.genomicLines(region, value,lwd=0.1, type="l", col="brown")
  x_centre = centromere[CELL_META$sector.index]
  circos.lines(c(x_centre, x_centre), c(0,1.2), col="green", lty=2, lwd=0.5)
  
})

legend("center", # or use x=0.5, y=0.5 for precise control
       legend = c("ALL", "LTR-retrotransposon", "CACTA", "Genes"),
       fill = c("black", "purple", "orange", "brown"), # use colors corresponding to your plot elements
       cex = 0.8, # adjust size if needed
       bty = "n") # removes the box around the legend

###########################
