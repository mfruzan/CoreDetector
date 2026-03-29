#code for Figure 1
chrlens <- c(598660471, 700547350, 498638509, 787782082,  812755788, 656544405, 754128162, 851934019, 619618552, 
              754227511, 673810255, 518332611, 713360525, 714697677, 569951140, 622669697, 731188232, 495380293, 744491536, 764072961, 642921167)
chrnames <- c("1A","1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D", "4A", "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D", "7A", "7B","7D")
dfcovg <- read.table("coverage.txt", header=F, row.names = 1, fill=TRUE)
rownames(dfcovg) <- substr(rownames(dfcovg),2,3)
dfcovg <- dfcovg[chrnames,]
sums <- rowSums(dfcovg)

ratios <- sums/chrlens
matdata <- c(ratios, 1-ratios)
matratio <- matrix(matdata, nrow=2, ncol=21, byrow=T)
colnames(matratio) <- chrnames
par(cex=0.8)
barplot(matratio, besidae=F)
total_ratio <- c(0.72, 0.28)
my_labels <- c("Core", " ")
cols <- c("Black", "grey")
pie(total_ratio, labels=my_labels, col=cols)
dev.off();
plot.new()

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
#this code is for visualizing gaps in wheat
chromspecies <- c("CS","spelt","mace","julius","norin61","lancer","mattis","jagger","landmark","arina","stanley", "paragon")
species <- c("CS","spelt","mace","julius","norin61","lancer","mattis","jagger","landmark","arina","stanley", "paragon","weebil","robigus","claire","cadenza")
gaps <- read.table("gaps.txt", sep="\t", header=F)
dim(gaps)
#13,265,923  6 
colnames(gaps) <- c("hsp", "specie", "chrom", "start", "realstart", "len")
length(unique(gaps$hsp))
#630,493  hsp's that contain gaps larger than 100bp
#aggregate based on chromosome, start postion and number of cultivars related to that gap
gap_agg1 <- aggregate(specie~hsp+start, gaps , FUN=function(v){return(length(unique(v)))})
gap_agg1 <- gap_agg1[order(gap_agg1$hsp, gap_agg1$start),]
max(gap_agg1$specie)
#15
uins <- gap_agg1[which(gap_agg1$specie==15)  ,]
dim(uins)
#580,696
merged <- merge(gaps, uins, by.x=c(1,4), by.y=c(1,2))
dim(merged)
#8,710,440       7
merged10k <- merged[merged$len>10000,]
length(which(merged10k$len>10000))/15
#4819 number of unique insertions larger than 10K
#now for each gap alignment we find the missing one that would be insertion specie
agg_func <- function(x){
  #cat(paste(length(x[,3]),","))
  if(dim(x)[1]==15){
     spc <- species[which(!species%in%x[,3])]
     #find chr and position of CS
     vv <- data.frame()
     if (spc!="CS"){
       vv <- x[x[,3]=="CS", 4:6]
     }else{
       #stanley is closest in terms of CNV to CS
       vv <- x[x[,3]=="stanley", 4:6]
     }
     return (list(spc, vv[1,1], vv[1,2], vv[1,3]))
  }else
    return (list(NA,NA,NA,NA))
}
res <- by(merged10k, merged10k[,c(1,2)], FUN=agg_func, simplify=FALSE)
res <- unlist(as.list(res))
length(res)
#19288
insdf <- data.frame(matrix(res, ncol=4, byrow =T))
dim(insdf)
#4822 4
insdf <- insdf[-which(rowSums(is.na(insdf))==4),]
dim(insdf)
#4818
insdf$X3 <- as.numeric(insdf$X3)
insdf$X4 <- as.numeric(insdf$X4)
max(insdf$X4)
#99198
insdf[which(insdf$X4==99198),]
colnames(insdf) <- c("specie","chrom","start","len")
#lancer 1B 139451389 99198
#see if chinese spring has unique insertion larger than 10K
insdf[insdf$specie=="CS",]
#3722     CS    7B 181434262 17128
#3839     CS    1A  49406678 11611
sum_ins <- aggregate(len~specie+chrom, insdf , FUN=length)
sum_ins <- sum_ins[sum_ins$specie%in%chromspecies,]
min(sum_ins$len)
#1
max(sum_ins$len)
#199
tab <-  xtabs(len ~ specie + chrom, data = sum_ins)

temat <- as.matrix(tab)
heatmap(temat, scale='none' )

#now for uniue insertion larger than 20k now visulize it in circlize for 12 cultivars
merged20k <- merged[merged$len>20000,]
dim(merged20k)[1]/15
#512
#
merged20k <- merged20k[merged20k$specie.x%in%chromspecies, ]
dim(merged20k)[1]/11
gap_agg2 <- aggregate(specie~hsp+start, gaps , FUN=function(v){return(length(v))})
max(gap_agg2$specie)
#15
max(gaps$len)
#99198
gaps[gaps$hsp==92740,]
#it shows that lancer has insertion of size 99198
#to find unique insertions, we look for entries having 15  


gaps_12 <- gaps[which(gaps$specie%in%chromspecies), ]
dim(gaps)
#13,265,923
sum_gaps <- aggregate(len~specie+chrom, gaps_12 , FUN=length)
min(sum_gaps$len)
#13,461
max(sum_gaps$len)
#57,988
tab <-  xtabs(len ~ specie + chrom, data = sum_gaps)

temat <- as.matrix(tab)
heatmap(temat, scale='none' )

gaps_12_10k <- gaps_12[which(gaps_12$len>10000),]
dim(gaps_12_10k)
#103,339    6
sum_gaps <- aggregate(len~specie+chrom, gaps_12_10k , FUN=length)
min(sum_gaps$len)
#58
max(sum_gaps$len)
#1173
tab <-  xtabs(len ~ specie + chrom, data = sum_gaps)

temat <- as.matrix(tab)
heatmap(temat, scale='none' )


#gaps larger that 10K
length(which(gaps_12$len>10000))
#94,665
length(which(gaps$len>30000))
#2142
#####################################################################################################
#bedtools intersect -wao -a /scratchdata1/users/a1195806/mario/references/CS2.1/whole/iwgsc_refseqv2.1_annotation_200916_HC.bed  -b /scratchdata1/users/a1195806/mario/projects/msa/wheat/minimap_out_sm_15/msa.bed > /scratchdata1/users/a1195806/mario/projects/msa/wheat/minimap_out_sm_15/intersect.bed
inter <- read.table("intersect.bed", sep="\t", header=F)
dim(inter)
#129,091      7
#aggregate based on gene region (v1,v2,v3) and return sum of core coverage for that region
agg <- aggregate(inter$V7, by=list(inter$V1, inter$V2, inter$V3), FUN=sum, na.rm=TRUE)
dim(agg)
# 105463      3
dim(agg[agg$x==0,])
#16977   3  16% of genes not covered at all
sum(agg$x)
# 287846655

sum(agg$Group.3-agg$Group.2)
#368986606
#78% genic region covered 

length(which((agg$Group.3-agg$Group.2 -2)<= agg$x))
#64346
#0.598646   61% of gene are fully covered
#now draw the graph of gene area
agg$million <- floor(agg$Group.2/1000000)
aggm <- aggregate(agg$x, by=list(agg$Group.1, agg$million), FUN=sum, na.rm=TRUE)
max(aggm$x)
#149664
aggm[aggm$x==max(aggm$x),]
#7B     592 149664
ymax <- 200000 
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
merged <- merge(notcovered, gen , by.x=c(1,2), by.y=c(2,3), all.x = T)
dim(merged)
#17022
length(unique(merged$V1))
write.table(merged$V1, file="notfoundgenes.txt", quote=F)


#############################################################################################
te <- read.table("Tae.Chinese_Spring.refSeqv2.1.TE_annotation.txt", sep="\t", header=F)
dim(te)
#4731260       6
tab <- table(te$V1,te$V5)
addmargins(tab)
tab2 <- table(te$V5, te$V6)
tab2
te$V2 <- te$V2 -1
te$V2 <- as.integer(te$V2) #this will avoid scientific notation
#te$V1 <- substr(te$V1,4,5)
write.table(te[,1:3], file="core/Tae.Chinese_Spring.refSeqv2.1.TE_annotation.bed", sep="\t", quote = F,  row.names = F, col.names=F)
#bedtools intersect -wao -a /scratchdata1/users/a1195806/mario/references/CS2.1/whole/Tae.Chinese_Spring.refSeqv2.1.TE_annotation.bed  -b /scratchdata1/users/a1195806/mario/projects/msa/wheat/minimap_out_sm_15/msa.bed > /scratchdata1/users/a1195806/mario/projects/msa/wheat/minimap_out_sm_15/intersect_te.bed
inter_te <- read.table("intersect_te.bed", sep="\t", header=F)
dim(inter_te)
#6,284,637       7
class(inter_te$V3)
#aggregate based on TE region and return sum of core coverage of that region
agg <- aggregate(inter_te$V7, by=list(inter_te$V1, inter_te$V2, inter_te$V3), FUN=sum, na.rm=TRUE)
dim(agg)
# 4,731,106       4   however this 4.7 is not unique TE because it contains match_part as well
max(inter_te$V2)
#0-based
sum(agg$x)
# 8,246,891,905
sum(inter_te$V7)
#8,246,891,905

sum(agg$Group.3-agg$Group.2)
#12,102,193,667
#0.6814378 TE region covered 

#to find unique entries of TE we need to merge it with Tae.Chinese_Spring.refSeqv2.1.TE_annotation.txt
#to match coordinates of bed and gff we need to decrease te$V2 by 1
merged <-  merge(inter_te, te , by=c(1,2,3), all.x = T)
dim(merged)
#6,285,027      10
#now aggregate based on V4.y (TE ID)
#now filter for complete TE
#merged <- merged[merged$V6.y=="complete",]
#dim(merged)
# 2,437,596      11
agg2 <- aggregate(merged$V7, by=list(merged$V4.y), FUN=sum, na.rm=TRUE)
dim(agg2)
#3,997,530       2 :This is total number of unique TE
#1109706       2  :This is total number of unique complete TE
dim(agg2[agg2$x==0,])
#863,904      2  21.6% of TE not covered at all
#213449      2   19.2% of complete TE not covered at all

#to find number of TE fully covered we need more information:
merged$len <- merged$V3 - merged$V2 + 1 
agg3 <- aggregate(merged$len, by=list(merged$V4.y), FUN=mean, na.rm=TRUE)
dim(agg3)
#3,997,530       2
#1,109,706       2
#reorder agg3 based on agg2
row.names(agg3) <- agg3$Group.1
agg3 <- agg3[agg2$Group.1, ]
length(which((agg3$x-3)<= agg2$x))
#2,360,925

#59.06% fully covered
#522,002
#47% of complete TE fully covered



#now we want separate positions of TE in core based on different TE types
levels(as.factor(merged$V5.y))
#"DHH" "DTA" "DTC" "DTH" "DTM" "DTT" "DTX" "DXX" "RIX" "RLC" "RLG" "RLX" "SIX" "XXX"
agg4 <- aggregate(merged$V7, by=list(merged$V5.y, merged$V1), FUN=sum, na.rm=TRUE)
#use output of agg4 to fill the supplementary table 1 in the manuscript
agg4[agg4$Group.1=="RLG",]
#now extract 
stat1 <- merged[,c(1,2,3,7,9)]
stat1 <- stat1[stat1$V7!=0,]
stat1 <- stat1[order(stat1$V1, stat1$V2),]
stat1$million <- floor(stat1$V2/1000000)
agg5 <- aggregate(stat1$V7, by=list(stat1$V1, stat1$million), FUN=sum, na.rm=TRUE)
max(agg5$x)
#1027102
#586648  complete TE
agg5[agg5$x==1027102,]
#2A     530 1027102
agg5[agg5$x==586648,]
#4B     185 586648
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
#-0.1209303
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
#0.26023

#now combination of LTR retrotransposons : RLC+RLG+RLX

stat2<- stat1[stat1$V5.y%in% c("RLC","RLG","RLX"),]
agg7 <- aggregate(stat2$V7, by=list(stat2$V1, stat2$million), FUN=sum, na.rm=TRUE)
max(agg7$x)
#902632
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
#-0.0449

stat2<- stat1[stat1$V5.y=="DTC",]
agg8 <- aggregate(stat2$V7, by=list(stat2$V1, stat2$million), FUN=sum, na.rm=TRUE)
max(agg8$x)
#361203
#156,602 complete CACTA transposon
agg8[agg8$x==361203,]
#3D     193 361203
ymax <- 400000
ymax <- 200000
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
#0.063


stat2<- stat1[stat1$V5.y=="RIX",]
agg9 <- aggregate(stat2$V7, by=list(stat2$V1, stat2$million), FUN=sum, na.rm=TRUE)
max(agg9$x)
#52976
agg9[agg9$x==52976,]
# 2B      26 52976
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
#0.25

stat2<- stat1[stat1$V5.y=="SIX",]
agg10 <- aggregate(stat2$V7, by=list(stat2$V1, stat2$million), FUN=sum, na.rm=TRUE)
gene_SIX <- merge(aggm, agg10, by=c(1,2))
cor(gene_SIX$x.x, gene_SIX$x.y)
#0.054

stat2<- stat1[stat1$V5.y=="DHH",]
agg11 <- aggregate(stat2$V7, by=list(stat2$V1, stat2$million), FUN=sum, na.rm=TRUE)
gene_DHH <- merge(aggm, agg11, by=c(1,2))
cor(gene_DHH$x.x, gene_DHH$x.y)
#0.030

stat2<- stat1[stat1$V5.y=="DTM",]
agg12 <- aggregate(stat2$V7, by=list(stat2$V1, stat2$million), FUN=sum, na.rm=TRUE)
gene_DTM <- merge(aggm, agg12, by=c(1,2))
cor(gene_DTM$x.x, gene_DTM$x.y)
#0.1889

stat2<- stat1[stat1$V5.y=="DTA",]
agg13 <- aggregate(stat2$V7, by=list(stat2$V1, stat2$million), FUN=sum, na.rm=TRUE)
gene_DTA <- merge(aggm, agg13, by=c(1,2))
cor(gene_DTA$x.x, gene_DTA$x.y)
#-0.05305163

stat2<- stat1[stat1$V5.y=="DTH",]
agg14 <- aggregate(stat2$V7, by=list(stat2$V1, stat2$million), FUN=sum, na.rm=TRUE)
gene_DTH <- merge(aggm, agg14, by=c(1,2))
cor(gene_DTH$x.x, gene_DTH$x.y)
#0.2868945
stat2<- stat1[stat1$V5.y=="DTT",]
agg15 <- aggregate(stat2$V7, by=list(stat2$V1, stat2$million), FUN=sum, na.rm=TRUE)
max(agg15$x)
#9323
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
#0.305503

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
#now look at 10k upstream TSS of all genes: Figure 4 of manuscript
df10kte <- read.table("10K_TSS_te.txt", sep="\t", header=F)
dim(df10kte)
#545,047     11
#now convert coordinates in respect to TSS as 0-10,000

corconvert <- function(row){
  c2 <- 0
  c3 <- 0
  vec <- as.numeric(row[c(1,2,4,5)])

  
if(row[3]=='-'){
#  cat(class(row[2]), ":" ,class(row[7]))
  c2 <- vec[1] - vec[3];
  c3 <- vec[2] - vec[3];
}else{
  c2 <- vec[4] - vec[2] 
  c3 <- vec[4] - vec[1] 
}
  return(c(c2,c3));
}

dfout <- apply(df10kte[,c(2,3,5,6,7)],1, corconvert)
dfout <- t(dfout)
dim(dfout)
# 545047  2
df10kte[,2:3] <- dfout
max(df10kte$V2)
#9998
min(df10kte$V2)
#0
all(df10kte$V2<df10kte$V3)
#TRUE

lst <- split(df10kte[,c(1,2,3,12)], df10kte$V12)
length(lst)
#14
class(lst[[1]])
#each element of list is dataframe
sum_up <- function(df){
  #returns a vector of length 10K
  vec <- rep(0,10000)
  
  for(i in 1:nrow(df)){
    vec[df[i,2]:df[i,3]] <- vec[df[i,2]:df[i,3]]+1
  }
  return(vec)
  
}
conv <- lapply(lst, FUN=sum_up)
class(conv)
#"list"
length(conv)
#14
class(conv[[2]])
# "numeric"
names(conv)
#"DHH" "DTA" "DTC" "DTH" "DTM" "DTT" "DTX" "DXX" "RIX" "RLC" "RLG" "RLX" "SIX" "XXX"
#red, green,blue,yellow,cyan,magenta,black,gray,orange,purple,brown,pink,teal

dev.off();
plot(0, 0, type = "l", xlim = c(0,9999),ylim = c(0,4000),      xlab = "Position upstream TSS", ylab = "Number of TE" )
lapply(conv, max)
lines(conv[["DTM"]], type='l',col='red')
lines(conv[["DTH"]], type='l',col='green')
lines(conv[["DTT"]], type='l',col='blue')
lines(conv[["DTA"]], type='l',col='cyan')
lines(conv[["RIX"]], type='l',col='magenta')
lines(conv[["SIX"]], type='l',col='orange')
lines(conv[["DTX"]], type='l',col='black')
legend("topright", c("Mutator", "Harbinger", "Mariner", "hAT", "LINE", "SINE", "Unclassified DT"), fill=c("red", "green","blue","cyan","magenta", "orange","black"), bty="n",cex=0.75, ncol = 2)

dev.off();
plot(0, 0, type = "l", xlim = c(0,9999),ylim = c(0,30000),      xlab = "Position upstream TSS", ylab = "Number of TE" )


lines(conv[["DTC"]], type='l',col='gray')
lines(conv[["DHH"]], type='l',col='orange')
lines(conv[["RLC"]], type='l',col='pink')
lines(conv[["RLG"]], type='l',col='#1B9E77')

legend("topleft", c("CACTA", "Helitron", "Copia", "Gypsy"), fill=c("grey", "orange","pink","#1B9E77"), bty="n",cex=0.75, ncol = 2)

dfcnv <- read.table("genes_cnv_90.bed", sep="\t", header=F)
dim(dfcnv)
#591 5

df10kte <- df10kte[df10kte$V4 %in%  dfcnv$V5, ]
dim(df10kte)
#2838  13
lst <- split(df10kte[,c(1,2,3,12)], df10kte$V12)
length(lst)
#14
conv <- lapply(lst, FUN=sum_up)
lapply(conv, max)

dev.off();
plot(0, 0, type = "l", xlim = c(0,9999),ylim = c(0,25),      xlab = "Position upstream TSS", ylab = "Number of TE" )
lapply(conv, max)
lines(conv[["DTM"]], type='l',col='red')
lines(conv[["DTH"]], type='l',col='green')
lines(conv[["DTT"]], type='l',col='blue')
lines(conv[["DTA"]], type='l',col='cyan')
lines(conv[["RIX"]], type='l',col='magenta')
lines(conv[["SIX"]], type='l',col='orange')
lines(conv[["DTX"]], type='l',col='black')
legend("topright", c("Mutator", "Harbinger", "Mariner", "hAT", "LINE", "SINE", "Unclassified DT"), fill=c("red", "green","blue","cyan","magenta", "orange","black"), bty="n",cex=0.75, ncol = 2)

dev.off();
plot(0, 0, type = "l", xlim = c(0,9999),ylim = c(0,130),      xlab = "Position upstream TSS", ylab = "Number of TE" )
lines(conv[["DTC"]], type='l',col='gray')
lines(conv[["DHH"]], type='l',col='orange')
lines(conv[["RLC"]], type='l',col='pink')
lines(conv[["RLG"]], type='l',col='#1B9E77')

legend("topleft", c("CACTA", "Helitron", "Copia", "Gypsy"), fill=c("grey", "orange","pink","#1B9E77"), bty="n",cex=0.75, ncol = 2)
######################################################################
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
####################################################################
#below code is for generating figure 2(heatmap) of manuscript
agg4 <- agg4[agg4$Group.1!="XXX",]
te_symbols <- c("Copia", "Gypsy", "Line", "Sine", "Unclassified LTR", "Helitron","Mutator", "hAT","CACTA","PIF/Harrbinger","Tc1/Mariner", "Unclassified TIR", "Unkonwn Transposon")
names(te_symbols) <- c("RLC","RLG","RIX","SIX","RLX","DHH","DTM","DTA","DTC","DTH","DTT", "DTX", "DXX")
agg4_high <- agg4[agg4$Group.1 %in% c("RLC","RLG","DTC", "RLX", "RIX") ,]
agg4_low <- agg4[!agg4$Group.1 %in% c("RLC","RLG","DTC", "RLX", "RIX") ,]
tab <-  xtabs(x ~ Group.1 + Group.2, data = agg4_low)
temat <- as.matrix(tab)
row.names(temat)<- te_symbols[row.names(temat)]
heatmap(temat, scale='none' )
min(temat)
mean(temat)
max(temat)
dev.off()
tab2 <-  xtabs(x ~ Group.1 + Group.2, data = agg4_high)
temat2 <- as.matrix(tab2)
row.names(temat2)<- te_symbols[row.names(temat2)]
heatmap(temat2, scale='none' )
min(temat2)
max(temat2)
####################################################################
#below code is for generating figure 3(heatmap) of manuscript
dfcnv <- read.table("cnv.txt", header=F, sep="\t")
rownames(dfcnv) <- dfcnv$V1
dfcnv <- dfcnv[,-1]
dfcnv <- dfcnv[,-ncol(dfcnv)]
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
min(dfcnv2)
max(dfcnv2)
#legend(x="topright", legend=c("min", "ave", "max"),  fill=colorRampPalette(brewer.pal(8, "Oranges"))(3))

#############################################################################################
#circlize coverage graphs
#First we need to re-format dfcov as: chr start end value
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
#genomics application

dfA_covg <- dffullcovg[dffullcovg$chr %in% c("1A","2A","3A","4A","5A","6A","7A"),]
dfA_covg <- dfA_covg[dfA_covg$value>0,]
dfA_covg$value <- dfA_covg$value/1000000;
dfA_covg$chr <- paste("chr", dfA_covg$chr, sep="")
max(dfA_covg$value)
#0.937
dfgA_covg <- dfgcovg[dfgcovg$chr %in% c("1A","2A","3A","4A","5A","6A","7A"),]
dfgA_covg$value <- dfgA_covg$value/1000000;
dfgA_covg$chr <- paste("chr", dfgA_covg$chr, sep="")
max(dfgA_covg$value)
#0.1324
retroA_covg <- agg7[agg7$chr %in% c("1A","2A","3A","4A","5A","6A","7A"),]
retroA_covg$value <- retroA_covg$value/1000000;
retroA_covg$chr <- paste("chr", retroA_covg$chr, sep="")
max(retroA_covg$value)
#0.902

dtcA_covg <- agg8[agg8$chr %in% c("1A","2A","3A","4A","5A","6A","7A"),]
dtcA_covg$value <- dtcA_covg$value/1000000;
dtcA_covg$chr <- paste("chr", dtcA_covg$chr, sep="")
max(dtcA_covg$value)
#0.26

centromere_A <- c(213000000, 335000000, 330000000, 313000000, 253000000, 290000000, 365000000)
names(centromere_A) <- c("chr1A","chr2A","chr3A","chr4A","chr5A","chr6A","chr7A")
centromere_B <- c(24000000, 358000000, 355000000, 300000000, 200000000, 333000000, 310000000)
names(centromere_B) <- c("chr1B","chr2B","chr3B","chr4B","chr5B","chr6B","chr7B")
centromere_D <- c(175000000, 273000000, 245000000, 205000000, 190000000, 235000000, 340000000)
names(centromere_D) <- c("chr1D","chr2D","chr3D","chr4D","chr5D","chr6D","chr7D")
#repeat below code 3 times
  genome <- "B"
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
  #0.937
  dfgA_covg <- dfgcovg[dfgcovg$chr %in% chrs,]
  dfgA_covg$value <- dfgA_covg$value/1000000;
  dfgA_covg$chr <- paste("chr", dfgA_covg$chr, sep="")
  max(dfgA_covg$value)
  #0.1324
  retroA_covg <- agg7[agg7$chr %in% chrs,]
  retroA_covg$value <- retroA_covg$value/1000000;
  retroA_covg$chr <- paste("chr", retroA_covg$chr, sep="")
  max(retroA_covg$value)
  #0.902
  
  dtcA_covg <- agg8[agg8$chr %in% chrs,]
  dtcA_covg$value <- dtcA_covg$value/1000000;
  dtcA_covg$chr <- paste("chr", dtcA_covg$chr, sep="")
  max(dtcA_covg$value)  
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
