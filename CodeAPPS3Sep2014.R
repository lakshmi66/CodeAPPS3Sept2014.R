library(polysat)

# Data manipulation -----

treedatamaster = read.table("/Users/Lakshmi/Desktop/Analysis/TreeDataFiles/TreeData.MasterJun2014.txt", header = TRUE)
treedatamaster$NewSite = paste(treedatamaster$Site, treedatamaster$Plot, sep = "")
treedatamaster$TreeID = paste(treedatamaster$NewSite, treedatamaster$X16th, treedatamaster$Tree.num, sep = ".")
row.names(treedatamaster) = treedatamaster$TreeID

analysisa = read.table("/Users/Lakshmi/Desktop/Labwork/AnalysisChecklistA.txt",fill=T,header=T)
analysisa$NewSite = paste(analysisa$Site, analysisa$Plot, sep = "")
analysisa$TreeID = paste(analysisa$NewSite,analysisa$X16th,analysisa$Tree.num, sep = ".")
row.names(analysisa) = analysisa$TreeID
analysisa$SampleNo1 

treeswithids = grep("[0-9]",analysisa$SampleNo1) ## vector of trees that have SampleIDs
analysisb = analysisa[treeswithids,]
twidsrownames = analysisa[treeswithids,"TreeID"] ## vector of TreeIDs for trees with SampleIDs
datamastera = treedatamaster[twidsrownames,]
datamastera$SampleID = analysisb$SampleNo1

datamasterb = datamastera[-which(datamastera$DBH<10),]
datamasterc = datamasterb[-which(datamasterb$SampleID==424|datamasterb$SampleID=="c583"|datamasterb$SampleID==520|datamasterb$SampleID=="c335"|datamasterb$SampleID=="c338"|datamasterb$SampleID=="c314"|datamasterb$SampleID=="c336"|datamasterb$SampleID=="c776"|datamasterb$SampleID=="c412"|datamasterb$SampleID==741|datamasterb$SampleID=="c442"|datamasterb$SampleID=="c447"|datamasterb$SampleID=="c325"|datamasterb$SampleID==540),]

bbadata = datamasterc[which(datamasterc$NewSite=="BBA"),]
bbbdata = datamasterc[which(datamasterc$NewSite=="BBB"),]
hrwadata = datamasterc[which(datamasterc$NewSite=="HRWA"),]
hrwbdata = datamasterc[which(datamasterc$NewSite=="HRWB"),]
rnpdata = datamasterc[which(datamasterc$NewSite=="RNPA"),]
pcdata = datamasterc[which(datamasterc$NewSite=="PCA"),]

# Importing and subsetting GeneMapper data -----

cleandata = read.GeneMapper("/Users/Lakshmi/Desktop/Analysis/mergedcleanGT.txt")
Usatnts(cleandata)
Usatnts(cleandata) <- c(4,3,2,4,4,2)

ploidycol = rep(6,length(Samples(cleandata)))

ploidymat = as.matrix(cbind(ploidycol,ploidycol,ploidycol,ploidycol,ploidycol,ploidycol))

colnames(ploidymat) = c("RW28","SEQ18D73","SEQ8E8","RW39","RW56","RWDI11") 

Ploidies(cleandata) = ploidymat

Genotype(cleandata, "1095", "RW56") = 205 ##updating fixed sample

# Making distance histograms of genetic distances by plot ----

dismatbb = meandistance.matrix(cleandata, samples = as.character(bbadata$SampleID), distmetric = Bruvo.distance)
hist(as.vector(dismatbb),breaks=30)

dismatbbb = meandistance.matrix(cleandata, samples = as.character(bbbdata$SampleID), distmetric = Bruvo.distance)
hist(as.vector(dismatbbb),breaks=30)

dismathrwa = meandistance.matrix(cleandata, samples = as.character(hrwadata$SampleID), distmetric = Bruvo.distance)
hist(as.vector(dismathrwa),breaks=30)

dismathrwb = meandistance.matrix(cleandata, samples = as.character(hrwbdata$SampleID), distmetric = Bruvo.distance)
hist(as.vector(dismathrwb),breaks=30)

dismatrnp = meandistance.matrix(cleandata, samples = as.character(rnpdata$SampleID), distmetric = Bruvo.distance)
hist(as.vector(dismatrnp),breaks=30)

dismatpc = meandistance.matrix(cleandata, samples = as.character(pcdata$SampleID), distmetric = Bruvo.distance)
hist(as.vector(dismatpc),breaks=30)

# Figure 1 -----
par$mar
postscript(file="figure1.rev2.eps")
par(mfrow=c(3,2),mar=c(4,4,4,2),oma=c(0,0,0,2))
hist(as.vector(dismatbb),breaks=30,xlab="Bruvo genetic distance",main="A",xlim=c(0,1),las=1)
hist(as.vector(dismatbbb),breaks=30,xlab="Bruvo genetic distance",main="B",xlim=c(0,1),las=1)
hist(as.vector(dismathrwa),breaks=30,xlab="Bruvo genetic distance",main="C",xlim=c(0,1),las=1)
hist(as.vector(dismathrwb),breaks=30,xlab="Bruvo genetic distance",main="D",xlim=c(0,1),las=1)
hist(as.vector(dismatrnp),breaks=30,xlab="Bruvo genetic distance",main="E",xlim=c(0,1),las=1)
hist(as.vector(dismatpc),breaks=30,xlab="Bruvo genetic distance",main="F",xlim=c(0,1),las=1)
dev.off()

# Making distance histograms of genetic distances by plot with Lynch distances -----

dismatbbL = meandistance.matrix(cleandata, samples = as.character(bbadata$SampleID), distmetric = Lynch.distance)
hist(as.vector(dismatbbL),breaks=30)

dismatbbbL = meandistance.matrix(cleandata, samples = as.character(bbbdata$SampleID), distmetric = Lynch.distance)
hist(as.vector(dismatbbbL),breaks=30)

dismathrwaL = meandistance.matrix(cleandata, samples = as.character(hrwadata$SampleID), distmetric = Lynch.distance)
hist(as.vector(dismathrwaL),breaks=30)

dismathrwbL = meandistance.matrix(cleandata, samples = as.character(hrwbdata$SampleID), distmetric = Lynch.distance)
hist(as.vector(dismathrwbL),breaks=30)

dismatrnpL = meandistance.matrix(cleandata, samples = as.character(rnpdata$SampleID), distmetric = Lynch.distance)
hist(as.vector(dismatrnpL),breaks=30)

dismatpcL = meandistance.matrix(cleandata, samples = as.character(pcdata$SampleID), distmetric = Lynch.distance)
hist(as.vector(dismatpcL),breaks=30)

par(mfrow=c(3,2))
par(mar=rep(2,4))

# Test samples with duplicates -----

dupvec = as.character(c(520,1048,1049,562,1050,1051,575,1052,1053,623,1054,1055,684,1056,1057,533,1058,1059,699,1060,1061,317,1069,319,318,321,320,323,322,327,326,339,338,342,343,344,345,346,347,354,355,386,387,379,378,426,1070,428,1071,439,1072,446,1073,459,1066,501,1075,385,384,513,1076,540,1077,551,1078,552,1079,565,1080,579,1081,587,1082,611,1083,617,1084,619,1085,628,1086,646,1087,671,1089,686,1090,692,1091,702,1067,707,1092,716,1062,724,1093,725,1094,729,1095,738,1068,742,1096,744,1064,746,1097))

testplate = as.character(c(562,1050,1051,575,1052,1053,623,1054,1055,684,1056,1057,533,1058,1059,699,1060,1061))

dismattestplate = meandistance.matrix(cleandata, samples = testplate, distmetric = Bruvo.distance)
testplateclones = assignClones(dismattestplate, threshold = 0.2)

plate1bb = as.character(c(317,1069,319,318,321,320,323,322,327,326,339,338,342,343,344,345,346,347,354,355,386,387,379,378,426,1070,428,1071,439,1072,446,1073,459,1066,501,1075,385,384))

dismatplate1bb = meandistance.matrix(cleandata, samples = plate1bb, distmetric = Bruvo.distance)
plate1bbclones = assignClones(dismatplate1bb, threshold = 0.2)

plate1hrw = as.character(c(513,1076,540,1077,551,1078,552,1079,565,1080,579,1081,587,1082,611,1083,617,1084,619,1085,628,1086,646,1087,671,1089,686,1090,692,1091))
dismatplate1hrw = meandistance.matrix(cleandata, samples = plate1hrw, distmetric = Bruvo.distance)
plate1hrwclones = assignClones(dismatplate1hrw, threshold = 0.2)

plate1rnp = as.character(c(702,1067,707,1092,716,1062,1063,724,1093,725,1094,729,1095,738,1068,742,1096,744,1064,1065,746,1097))
dismatplate1rnp = meandistance.matrix(cleandata, samples = plate1rnp, distmetric = Bruvo.distance)
plate1rnpclones = assignClones(dismatplate1rnp, threshold = 0.2)

basal = as.character(c(1108,492,1110,610,645,683,703,709,721,731,1236,1237,1238,1239,1240,1242,1243,1244,1245,1247,1250))

basalmatches = as.character(c("c220","c291","c320","c358","c372","c391","c408","c420","c443","c458","c684","c700","c720","c732","c744","c784","c615","c618","1143","c628","c655"))

viewGenotypes(cleandata, samples = c("729","1095"))

dismatbasal = meandistance.matrix(cleandata, samples = c(basal,basalmatches), distmetric = Bruvo.distance)
basalclones = assignClones(dismatbasal, threshold = 0.2)

forkedtrees = as.character(c(740,741,"c416","c417","c428","c429","c441","c442","c446","c447","c459","c460","c461","c473","c474","c552","c553","c598","c599","c604","c605","c618","c619","c629","c630","c678","c679"))
dismatforked = meandistance.matrix(cleandata, samples = forkedtrees, distmetric = Bruvo.distance)
forkedclones = assignClones(dismatforked, threshold = 0.2)

viewGenotypes(cleandata, samples = c("c678","c679"))

# Eliminating clones within plots ------

distbba = meandistance.matrix(cleandata, samples = as.character(bbadata$SampleID), distmetric = Bruvo.distance)

bbaclones = assignClones(distbba,threshold=0.2)
bbaincl = names(bbaclones[-which(duplicated(bbaclones))])

distbbb = meandistance.matrix(cleandata, samples = as.character(bbbdata$SampleID), distmetric = Bruvo.distance)

bbbclones = assignClones(distbbb,threshold=0.2)
bbbincl = names(bbbclones[-which(duplicated(bbbclones))])

disthrwa = meandistance.matrix(cleandata, samples = as.character(hrwadata$SampleID), distmetric = Bruvo.distance)

hrwaclones = assignClones(disthrwa,threshold=0.2)
hrwaincl = names(hrwaclones[-which(duplicated(hrwaclones))])

disthrwb = meandistance.matrix(cleandata, samples = as.character(hrwbdata$SampleID), distmetric = Bruvo.distance)

hrwbclones = assignClones(disthrwb,threshold=0.2)
hrwbincl = names(hrwbclones[-which(duplicated(hrwbclones))])

distrnp = meandistance.matrix(cleandata, samples = as.character(rnpdata$SampleID), distmetric = Bruvo.distance)

rnpclones = assignClones(distrnp,threshold=0.2)
rnpincl = names(rnpclones[-which(duplicated(rnpclones))])

distpc = meandistance.matrix(cleandata, samples = as.character(pcdata$SampleID), distmetric = Bruvo.distance)

pcclones = assignClones(distpc,threshold=0.2)
pcincl = names(pcclones[-which(duplicated(pcclones))])

# Combining samples from different plots to summarize allele counts (Table 1) ------

incl = c(bbaincl,bbbincl,hrwaincl,hrwbincl,rnpincl,pcincl) ##450 samples
length(incl)
Genotype(cleandata, sample = incl[446], locus = "RW28")

#RW28

genomat = matrix(2,nrow=450,ncol=7)
genomat[,1] = incl
tail(genomat)
i=449
for(i in 1:length(incl)){
  for(j in 2:7){
    genomat[i,j] = Genotype(cleandata, sample = incl[i], locus = "RW28")[j-1]
  }
}

genomat28 = as.data.frame(genomat)
head(genomat28)

allelecountsrw28 = rep(7,450)
for(i in 1:length(incl)){
  allelecountsrw28[i] = length(which(!is.na(genomat28[i,]) & genomat28[i,]!="-9"))-1
}

class(allelecountsrw28)
length(which(allelecountsrw28==1))
length(which(allelecountsrw28==2))
length(which(allelecountsrw28==3))
length(which(allelecountsrw28==4))
length(which(allelecountsrw28==0))

#RW39

genomatrw39 = matrix(2,nrow=450,ncol=7)
genomatrw39[,1] = incl
tail(genomat)
for(i in 1:length(incl)){
  for(j in 2:7){
    genomatrw39[i,j] = Genotype(cleandata, sample = incl[i], locus = "RW39")[j-1]
  }
}
Genotype(cleandata, sample = incl[450], locus = "RW39")
genomatrw39 = as.data.frame(genomatrw39)
tail(genomatrw39)

allelecountsrw39 = rep(7,450)
for(i in 1:length(incl)){
  allelecountsrw39[i] = length(which(!is.na(genomatrw39[i,]) & genomatrw39[i,]!="-9"))-1
}

class(allelecountsrw28)
length(which(allelecountsrw39==1))
length(which(allelecountsrw39==2))
length(which(allelecountsrw39==3))
length(which(allelecountsrw39==4))
length(which(allelecountsrw39==5))
length(which(allelecountsrw39==6))
length(which(allelecountsrw39==0))
2+9+52+138+149+97+3

# SEQ18D73
genomatseq18d73 = matrix(2,nrow=450,ncol=7)
genomatseq18d73[,1] = incl

for(i in 1:length(incl)){
  for(j in 2:7){
    genomatseq18d73[i,j] = Genotype(cleandata, sample = incl[i], locus = "SEQ18D73")[j-1]
  }
}
Genotype(cleandata, sample = incl[1], locus = "SEQ18D73")
genomatseq18d73 = as.data.frame(genomatseq18d73)
tail(genomatseq18d73)

allelecountsseq18d73 = rep(7,450)
for(i in 1:length(incl)){
  allelecountsseq18d73[i] = length(which(!is.na(genomatseq18d73[i,]) & genomatseq18d73[i,]!="-9"))-1
}

length(which(allelecountsseq18d73==1))
length(which(allelecountsseq18d73==2))
length(which(allelecountsseq18d73==3))
length(which(allelecountsseq18d73==4))
length(which(allelecountsseq18d73==5))
length(which(allelecountsseq18d73==6))
length(which(allelecountsseq18d73==0))
136+180+88+15+5+26

# SEQ8E8
genomatseq8e8 = matrix(2,nrow=450,ncol=7)
genomatseq8e8[,1] = incl

for(i in 1:length(incl)){
  for(j in 2:7){
    genomatseq8e8[i,j] = Genotype(cleandata, sample = incl[i], locus = "SEQ8E8")[j-1]
  }
}
Genotype(cleandata, sample = incl[1], locus = "SEQ8E8")
genomatseq8e8 = as.data.frame(genomatseq8e8)
head(genomatseq8e8)

allelecountsseq8e8 = rep(7,450)
for(i in 1:length(incl)){
  allelecountsseq8e8[i] = length(which(!is.na(genomatseq8e8[i,]) & genomatseq8e8[i,]!="-9"))-1
}

length(which(allelecountsseq8e8==1))
length(which(allelecountsseq8e8==2))
length(which(allelecountsseq8e8==3))
length(which(allelecountsseq8e8==4))
length(which(allelecountsseq8e8==5))
length(which(allelecountsseq8e8==6))
length(which(allelecountsseq8e8==0))
246+65+9+2+128

# RW56
genomatrw56 = matrix(2,nrow=450,ncol=7)
genomatrw56[,1] = incl

for(i in 1:length(incl)){
  for(j in 2:7){
    genomatrw56[i,j] = Genotype(cleandata, sample = incl[i], locus = "RW56")[j-1]
  }
}
Genotype(cleandata, sample = incl[1], locus = "RW56")
genomatrw56 = as.data.frame(genomatrw56)
tail(genomatrw56)

allelecountsrw56 = rep(7,450)
for(i in 1:length(incl)){
  allelecountsrw56[i] = length(which(!is.na(genomatrw56[i,]) & genomatrw56[i,]!="-9"))-1
}

length(which(allelecountsrw56==1))
length(which(allelecountsrw56==2))
length(which(allelecountsrw56==3))
length(which(allelecountsrw56==4))
length(which(allelecountsrw56==5))
length(which(allelecountsrw56==6))
length(which(allelecountsrw56==0))

# RWDI11
genomatrwdi11 = matrix(2,nrow=450,ncol=7)
genomatrwdi11[,1] = incl

for(i in 1:length(incl)){
  for(j in 2:7){
    genomatrwdi11[i,j] = Genotype(cleandata, sample = incl[i], locus = "RWDI11")[j-1]
  }
}
Genotype(cleandata, sample = incl[1], locus = "RWDI11")
genomatrwdi11 = as.data.frame(genomatrwdi11)
head(genomatrwdi11)

allelecountsrwdi11 = rep(7,450)
for(i in 1:length(incl)){
  allelecountsrwdi11[i] = length(which(!is.na(genomatrwdi11[i,]) & genomatrwdi11[i,]!="-9"))-1
}

length(which(allelecountsrwdi11==1))
length(which(allelecountsrwdi11==2))
length(which(allelecountsrwdi11==3))
length(which(allelecountsrwdi11==4))
length(which(allelecountsrwdi11==5))
length(which(allelecountsrwdi11==6))
length(which(allelecountsrwdi11==0))

# Counting different types of duplicate samples -----

duptable = read.delim("/Users/Lakshmi/Desktop/Analysis/duplicatelist.txt", header = TRUE)

sampleswithdups = c(as.character(duptable$sample1),as.character(duptable$sample2)) 
dups2 = unique(sampleswithdups) #making a list of all samples

Genotype(cleandata, "1095", "RW56") = 205 ##updating fixed sample

dismatdups2 = meandistance.matrix(cleandata, samples = dups2, distmetric = Bruvo.distance)  

dist = rep(2,88)
for(i in 1:length(dist)){
  dist[i] = dismatdups2[as.character(duptable$sample1[i]),as.character(duptable$sample2[i])]
}

duptable$dist = dist
  
duptableno709 = duptable[-which(duptable$sample1==709),]

duptableno709nofolcam = duptableno709[-which(duptableno709$type=="fol.cam"),]

duptableno709nofolcam$type = factor(duptableno709nofolcam$type,levels=c("fol.epi","fol.bas","epi.bas","bas.cam"))
plot(dist~type,data=duptableno709nofolcam)

anova = aov(dist~type, data = duptableno709nofolcam)
summary(anova)
plot(residuals(anova)~fitted(anova))

lm1 = lm(dist~type, data = duptableno709nofolcam)
summary(lm1)
TukeyHSD(anova)
shapiro.test(residuals(anova))

dim(duptableno709nofolcam)
duptableno709nofolcam$type
duptab3 = duptableno709nofolcam
tail(duptab3)

folvec = rep(2,86)
for(i in 1:length(folvec)){
  folvec[i] = ifelse(duptab3[i,2]=="fol.epi"|duptab3[i,2]=="fol.bas",1,0)
}
duptab3$fol = folvec

epivec = rep(2,86)
for(i in 1:length(epivec)){
  epivec[i] = ifelse(duptab3[i,2]=="fol.epi"|duptab3[i,2]=="epi.bas",1,0)
}
duptab3$epi = epivec

basvec = rep(2,86)
for(i in 1:length(basvec)){
  basvec[i] = ifelse(duptab3[i,2]=="fol.bas"|duptab3[i,2]=="epi.bas"|duptab3[i,2]=="bas.cam",1,0)
}
duptab3$bas = basvec

camvec = rep(2,86)
for(i in 1:length(camvec)){
  camvec[i] = ifelse(duptab3[i,2]=="bas.cam",1,0)
}
duptab3$cam = camvec

lm2 = lm(dist ~ fol + epi + bas + cam, data=duptab3)
summary(lm2)

duptab3[which(duptab3$type=="bas.cam"&duptab3$dist>0),]

viewGenotypes(cleandata, samples = c("1250","c655"))

# dup table with no zeros

duptablenozeroes = duptable2[-which(duptable2$dist==0),]
duptablenozeroes[which(duptablenozeroes$type=="fol.cam"),]
duptablenozeroes$type = factor(duptablenozeroes$type,levels=c("fol.epi","fol.bas","epi.bas","bas.cam"))
plot(dist~type,data=duptablenozeroes)
anova2 = aov(dist~type,data=duptablenozeroes)
summary(anova)

viewGenotypes(cleandata, samples = c("729","1095"))

# counting zero and nonzero distances in each category
dim(duptablenofolcam[which(duptablenofolcam$dist==0 & duptablenofolcam$type=="bas.cam"),])

fol.epi = c(17,5)
fol.bas = c(9,7)
epi.bas = c(17,12)
bas.cam = c(6,14)

chisqtab = rbind(fol.epi,fol.bas,epi.bas,bas.cam)

test = chisq.test(chisqtab)

library(ggplot2)
library(labeling)

# Looking for relationship between genetic distance and proportion of loci amplified in same plate for paired samples ------

duptab = read.delim("/Users/Lakshmi/Desktop/Analysis/duplicatelist2.txt", header = TRUE)

sampleswithdups = c(as.character(duptab$sample1),as.character(duptab$sample2)) 
dups3 = unique(sampleswithdups) #making a list of all samples

dismatdups3 = meandistance.matrix(cleandata, samples = dups3, distmetric = Bruvo.distance)  

dist = rep(2,88)
for(i in 1:length(dist)){
  dist[i] = dismatdups3[as.character(duptab$sample1[i]),as.character(duptab$sample2[i])]
}

duptab$dist = dist

## eliminating the outlier

duptab4 = duptab[-which(duptab$sample1==709 | duptab$type=="fol.cam"),]
duptab4$type = factor(duptab4$type,levels=c("fol.epi","fol.bas","epi.bas","bas.cam"))

lm1 = lm(dist~proplocionsameplate+type, data = duptab4)
summary(lm1)
plot(residuals(lm1)~duptab4$proplocionsameplate)

##add plot with symbols

library(ggplot2)
p = ggplot(duptab4,aes(proplocionsameplate,dist,shape=type))
p + geom_point(position = position_jitter(w=0.01,h=0)) + scale_shape_manual(values=c(1,2,3,4)) + theme_bw()

##hard to look at.  barplot instead

proplocibin = rep(11,dim(duptab4)[1])
for(i in 1:length(proplocibin)){
  proplocibin[i] = ifelse(0<=duptab4[i,6]&duptab4[i,6]<0.1,"<0.1",ifelse(0.1<=duptab4[i,6]&duptab4[i,6]<0.2,"0.1-0.2",ifelse(0.2<=duptab4[i,6]&duptab4[i,6]<0.3,"0.2-0.3",ifelse(0.3<=duptab4[i,6]&duptab4[i,6]<0.4,"0.3-0.4",ifelse(0.4<=duptab4[i,6]&duptab4[i,6]<0.5,"0.4-0.5",ifelse(0.5<=duptab4[i,6]&duptab4[i,6]<0.6,"0.5-0.6",ifelse(0.6<=duptab4[i,6]&duptab4[i,6]<0.7,"0.6-0.7",ifelse(0.7<=duptab4[i,6]&duptab4[i,6]<0.8,"0.7-0.8",ifelse(0.8<=duptab4[i,6]&duptab4[i,6]<0.9,"0.8-0.9",ifelse(duptab4[i,6]>=0.9,"0.9-1",NA))))))))))
}

duptab4$proplocibin = proplocibin

qplot(factor(proplocibin), data=duptab4, geom="bar", fill=factor(type),)

## duplicate samples without bas.cam

duptab = read.delim("/Users/Lakshmi/Desktop/Analysis/duplicatelist3.txt", header = TRUE)

sampleswithdups = c(as.character(duptab$sample1),as.character(duptab$sample2)) 
dups3 = unique(sampleswithdups) #making a list of all samples

dismatdups3 = meandistance.matrix(cleandata, samples = dups3, distmetric = Bruvo.distance)  

dist = rep(2,67)
for(i in 1:length(dist)){
  dist[i] = dismatdups3[as.character(duptab$sample1[i]),as.character(duptab$sample2[i])]
}

duptab$dist = dist

## eliminating the outlier

duptab$type = factor(duptab$type,levels=c("fol.epi","fol.bas","epi.bas","bas.cam"))

plot(dist~proplocionsameplate, data=duptab)

lm2 = lm(dist~proplocionsameplate, data=duptab)
summary(lm2)

# Figure 3 -----

##dotplot with zeroes
?tiff
tiff(file="figure3b.tiff",res=1200,units="in",height=6,width=6)
p = ggplot(duptab3,aes(factor(type),dist))
p + geom_point(alpha = 0.3, position = position_jitter(w=0.10,h=0)) + theme_bw() + xlab("Tissue types") + ylab("Bruvo genetic distance")
dev.off()


##dotplot without zeroes

p = ggplot(duptablenozeroes,aes(factor(type),dist),label=Name)

p + geom_text(aes(label = "o",size=0.05),position = position_jitter(w=0.05,h=0))

##revised version scaling dots by number of observations-original data in "duptab3"

folepitab = duptab3[which(duptab3$type=="fol.epi"),]
folepidist = unique(folepitab$dist)
countfolepi = rep(100,length(folepidist))
for(i in 1:length(countfolepi)){
  countfolepi[i] = dim(folepitab[which(folepitab$dist==folepidist[i]),])[1]
}
folepitab = cbind(rep("foliage.epicormic",length(folepidist)),folepidist,countfolepi)
colnames(folepitab) = c("type","dist","count")

folbastab = duptab3[which(duptab3$type=="fol.bas"),]
folbasdist = unique(folbastab$dist)
countfolbas = rep(100,length(folbasdist))
for(i in 1:length(countfolbas)){
  countfolbas[i] = dim(folbastab[which(folbastab$dist==folbasdist[i]),])[1]
}
folbastab = cbind(rep("foliage.basal",length(folbasdist)),folbasdist,countfolbas)
colnames(folbastab) = c("type","dist","count")

epibastab = duptab3[which(duptab3$type=="epi.bas"),]
epibasdist = unique(epibastab$dist)
countepibas = rep(100,length(epibasdist))
for(i in 1:length(countepibas)){
  countepibas[i] = dim(epibastab[which(epibastab$dist==epibasdist[i]),])[1]
}
epibastab = cbind(rep("epicormic.basal",length(epibasdist)),epibasdist,countepibas)
colnames(epibastab) = c("type","dist","count")

bascamtab = duptab3[which(duptab3$type=="bas.cam"),]
bascamdist = unique(bascamtab$dist)
countbascam = rep(100,length(bascamdist))
for(i in 1:length(countbascam)){
  countbascam[i] = dim(bascamtab[which(bascamtab$dist==bascamdist[i]),])[1]
}
bascamtab = cbind(rep("basal.cambium",length(bascamdist)),bascamdist,countbascam)
colnames(bascamtab) = c("type","dist","count")

fig3tab = as.data.frame(rbind(folepitab,folbastab,epibastab,bascamtab))
fig3tab$count = as.numeric(as.character(fig3tab$count))
fig3tab$dist = as.numeric(as.character(fig3tab$dist))

postscript(file="figure3rev2.eps")
par(mar=c(4,4,4,4))

p = ggplot(fig3tab, aes(factor(type),dist))
p + geom_point(aes(size = count)) + theme_bw() + labs(size="# pairs", x= "Tissue type", y="Bruvo genetic distance") + scale_y_continuous(limits=c(0,0.2), breaks=seq(0,0.2,0.01))

dev.off()
ggsave("figure3rev2.eps",dpi=1200)

# Calculating He based on allele frequencies -----

rw28freqtab  = read.table("/Users/Lakshmi/Desktop/Analysis/allelefreqrw28.txt",sep = "\t",header = TRUE)
herw28 = 1 - sum(rw28freqtab^6)  

rw39freqtab  = read.table("/Users/Lakshmi/Desktop/Analysis/allelefreqrw39.txt",sep = "\t",header = TRUE)
herw39 = 1 - sum(rw39freqtab^6)

seq18d73freqtab  = read.table("/Users/Lakshmi/Desktop/Analysis/allelefreqseq18d73.txt",sep = "\t",header = TRUE)
heseq18d73 = 1 - sum(seq18d73freqtab^6)

seq8e8freqtab = read.table("/Users/Lakshmi/Desktop/Analysis/allelefreqseq8e8.txt",sep = "\t",header = TRUE)
heseq8e8 = 1 - sum(seq8e8freqtab^6)

rw56freqtab = read.table("/Users/Lakshmi/Desktop/Analysis/allelefreqrw56.txt",sep = "\t",header = TRUE)
herw56 = 1 - sum(rw56freqtab^6)

rwdi11freqtab  = read.table("/Users/Lakshmi/Desktop/Analysis/allelefreqrwdi11.txt",sep = "\t",header = TRUE)
herwdi11 = 1 - sum(rwdi11freqtab^6)

# Testing what happens with increasing numbers of null alleles -----

## run before for loop
library(plyr)
seq18d73freqtab  = read.table("/Users/Lakshmi/Desktop/Analysis/allelefreqseq18d73.txt",sep = "\t",header = TRUE) 
seq18d73freq = as.numeric(seq18d73freqtab)
seq18d73namevec = as.character(names(seq18d73freqtab)) 
seq18d73namevec2 = unlist(strsplit(seq18d73namevec,"[.]"))[1:length(seq18d73namevec)*2]
rw39freqtab  = read.table("/Users/Lakshmi/Desktop/Analysis/allelefreqrw39.txt",sep = "\t",header = TRUE) 
rw39freq = as.numeric(rw39freqtab)
rw39namevec = as.character(names(rw39freqtab)) 
rw39namevec2 = unlist(strsplit(rw39namevec,"[.]"))[1:length(rw39namevec)*2]
rw28freqtab  = read.table("/Users/Lakshmi/Desktop/Analysis/allelefreqrw28.txt",sep = "\t",header = TRUE) 
rw28freq = as.numeric(rw28freqtab)
rw28namevec = as.character(names(rw28freqtab)) 
rw28namevec2 = unlist(strsplit(rw28namevec,"[.]"))[1:length(rw28namevec)*2]
seq8e8freqtab = read.table("/Users/Lakshmi/Desktop/Analysis/allelefreqseq8e8.txt",sep = "\t",header = TRUE) 
seq8e8freq = as.numeric(seq8e8freqtab)
seq8e8namevec = as.character(names(seq8e8freqtab)) 
seq8e8namevec2 = unlist(strsplit(seq8e8namevec,"[.]"))[1:length(seq8e8namevec)*2]
rw56freqtab  = read.table("/Users/Lakshmi/Desktop/Analysis/allelefreqrw56.txt",sep = "\t",header = TRUE) 
rw56freq = as.numeric(rw56freqtab)
rw56namevec = as.character(names(rw56freqtab)) 
rw56namevec2 = unlist(strsplit(rw56namevec,"[.]"))[1:length(rw56namevec)*2]
rwdi11freqtab  = read.table("/Users/Lakshmi/Desktop/Analysis/allelefreqrwdi11.txt",sep = "\t",header = TRUE) 
rwdi11freq = as.numeric(rwdi11freqtab)
rwdi11namevec = as.character(names(rwdi11freqtab)) 
rwdi11namevec2 = unlist(strsplit(rwdi11namevec,"[.]"))[1:length(rwdi11namevec)*2]

# Testing: 0 null alleles -----
matchvec = rep(190,94)
n=1
for(i in 1:94){
  ##simulating seq18d73
  seq18d73samplevec = sample(seq18d73namevec2,1092,replace=TRUE,prob=seq18d73freq)
  dfseq18d73 = matrix(seq18d73samplevec,nrow=182,ncol=6)
  dfseq18d73b = (apply(dfseq18d73,1,unique)) #eliminating duplicate alleles 
  dfseq18d73c = ldply(dfseq18d73b, rbind) #converting from list to data frame
  matseq18d73 = as.matrix(dfseq18d73c)
  dfseq18d73d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfseq18d73)[1]){
    testrow = matseq18d73[i,]
    dfseq18d73d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("SEQ18D73",182)
  dfseq18d73e = as.data.frame(dfseq18d73d)
  colnames(dfseq18d73e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfseq18d73f = cbind(Sample.Name,Marker,dfseq18d73e)
  ##simulating RW39
  rw39samplevec = sample(rw39namevec2,1092,replace=TRUE,prob=rw39freq)
  dfrw39 = matrix(rw39samplevec,nrow=182,ncol=6)
  dfrw39b = (apply(dfrw39,1,unique)) #eliminating duplicate alleles 
  dfrw39c = ldply(dfrw39b, rbind) #converting from list to data frame
  matrw39 = as.matrix(dfrw39c)
  dfrw39d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrw39)[1]){
    testrow = matrw39[i,]
    dfrw39d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RW39",182)
  dfrw39e = as.data.frame(dfrw39d)
  colnames(dfrw39e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrw39f = cbind(Sample.Name,Marker,dfrw39e)
  ##simulating RW28
  rw28samplevec = sample(rw28namevec2,1092,replace=TRUE,prob=rw28freq)
  dfrw28 = matrix(rw28samplevec,nrow=182,ncol=6)
  dfrw28b = (apply(dfrw28,1,unique)) #eliminating duplicate alleles 
  dfrw28c = ldply(dfrw28b, rbind) #converting from list to data frame
  matrw28 = as.matrix(dfrw28c)
  dfrw28d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrw28)[1]){
    testrow = matrw28[i,]
    dfrw28d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RW28",182)
  dfrw28e = as.data.frame(dfrw28d)
  colnames(dfrw28e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrw28f = cbind(Sample.Name,Marker,dfrw28e)
  ##simulating SEQ8E8
  seq8e8samplevec = sample(seq8e8namevec2,1092,replace=TRUE,prob=seq8e8freq)
  dfseq8e8 = matrix(seq8e8samplevec,nrow=182,ncol=6)
  dfseq8e8b = (apply(dfseq8e8,1,unique)) #eliminating duplicate alleles 
  dfseq8e8c = ldply(dfseq8e8b, rbind) #converting from list to data frame
  matseq8e8 = as.matrix(dfseq8e8c)
  dfseq8e8d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfseq8e8)[1]){
    testrow = matseq8e8[i,]
    dfseq8e8d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("SEQ8E8",182)
  dfseq8e8e = as.data.frame(dfseq8e8d)
  colnames(dfseq8e8e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfseq8e8f = cbind(Sample.Name,Marker,dfseq8e8e)
  ##simulating RW56
  rw56samplevec = sample(rw56namevec2,1092,replace=TRUE,prob=rw56freq)
  dfrw56 = matrix(rw56samplevec,nrow=182,ncol=6)
  dfrw56b = (apply(dfrw56,1,unique)) #eliminating duplicate alleles 
  dfrw56c = ldply(dfrw56b, rbind) #converting from list to data frame
  matrw56 = as.matrix(dfrw56c)
  dfrw56d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrw56)[1]){
    testrow = matrw56[i,]
    dfrw56d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RW56",182)
  dfrw56e = as.data.frame(dfrw56d)
  colnames(dfrw56e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrw56f = cbind(Sample.Name,Marker,dfrw56e)
  ##simulating RWDI11
  rwdi11samplevec = sample(rwdi11namevec2,1092,replace=TRUE,prob=rwdi11freq)
  dfrwdi11 = matrix(rwdi11samplevec,nrow=182,ncol=6)
  dfrwdi11b = (apply(dfrwdi11,1,unique)) #eliminating duplicate alleles 
  dfrwdi11c = ldply(dfrwdi11b, rbind) #converting from list to data frame
  matrwdi11 = as.matrix(dfrwdi11c)
  dfrwdi11d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrwdi11)[1]){
    testrow = matrwdi11[i,]
    dfrwdi11d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RWDI11",182)
  dfrwdi11e = as.data.frame(dfrwdi11d)
  colnames(dfrwdi11e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrwdi11f = cbind(Sample.Name,Marker,dfrwdi11e)
  ##combining markers, exporting and importing files
  allmarkers=rbind(dfseq18d73f,dfrw39f,dfrw28f,dfseq8e8f,dfrw56f,dfrwdi11f)
  write.table(allmarkers,"/Users/Lakshmi/Desktop/Analysis/simdata.txt",sep = "\t",quote = FALSE, row.names=FALSE)
  simdata = read.GeneMapper("/Users/Lakshmi/Desktop/Analysis/simdata.txt")
  Usatnts(simdata) <- c(3,4,4,2,4,2)
  ploidycol = rep(6,length(Samples(simdata)))
  ploidymat = as.matrix(cbind(ploidycol,ploidycol,ploidycol,ploidycol,ploidycol,ploidycol))
  colnames(ploidymat) = c("SEQ18D73","RW39","RW28","SEQ8E8","RW56","RWDI11") 
  Ploidies(simdata) = ploidymat
  simdist = meandistance.matrix(simdata,distmetric = Bruvo.distance)
  simclones = assignClones(simdist,threshold=0.2)
  matches = length(simclones)-length(unique(simclones))
  matchvec[n] = matches
  n = n+1
}

# Testing: 1 null allele per sample -----
matchvec = rep(190,100)
n=1
i=1
##change to i in 1:100
for(i in 1:100){
  ##simulating seq18d73
  seq18d73samplevec = sample(seq18d73namevec2,1092,replace=TRUE,prob=seq18d73freq)
  dfseq18d73 = matrix(seq18d73samplevec,nrow=182,ncol=6)
  dfseq18d73b = (apply(dfseq18d73,1,unique)) #eliminating duplicate alleles 
  dfseq18d73c = ldply(dfseq18d73b, rbind) #converting from list to data frame
  matseq18d73 = as.matrix(dfseq18d73c)
  dfseq18d73d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfseq18d73)[1]){
    testrow = matseq18d73[i,]
    dfseq18d73d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("SEQ18D73",182)
  dfseq18d73e = as.data.frame(dfseq18d73d)
  colnames(dfseq18d73e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfseq18d73f = cbind(Sample.Name,Marker,dfseq18d73e)
  ##simulating RW39
  rw39samplevec = sample(rw39namevec2,1092,replace=TRUE,prob=rw39freq)
  dfrw39 = matrix(rw39samplevec,nrow=182,ncol=6)
  dfrw39b = (apply(dfrw39,1,unique)) #eliminating duplicate alleles 
  dfrw39c = ldply(dfrw39b, rbind) #converting from list to data frame
  matrw39 = as.matrix(dfrw39c)
  dfrw39d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrw39)[1]){
    testrow = matrw39[i,]
    dfrw39d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RW39",182)
  dfrw39e = as.data.frame(dfrw39d)
  colnames(dfrw39e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrw39f = cbind(Sample.Name,Marker,dfrw39e)
  ##simulating RW28
  rw28samplevec = sample(rw28namevec2,1092,replace=TRUE,prob=rw28freq)
  dfrw28 = matrix(rw28samplevec,nrow=182,ncol=6)
  dfrw28b = (apply(dfrw28,1,unique)) #eliminating duplicate alleles 
  dfrw28c = ldply(dfrw28b, rbind) #converting from list to data frame
  matrw28 = as.matrix(dfrw28c)
  dfrw28d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrw28)[1]){
    testrow = matrw28[i,]
    dfrw28d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RW28",182)
  dfrw28e = as.data.frame(dfrw28d)
  colnames(dfrw28e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrw28f = cbind(Sample.Name,Marker,dfrw28e)
  ##simulating SEQ8E8
  seq8e8samplevec = sample(seq8e8namevec2,1092,replace=TRUE,prob=seq8e8freq)
  dfseq8e8 = matrix(seq8e8samplevec,nrow=182,ncol=6)
  dfseq8e8b = (apply(dfseq8e8,1,unique)) #eliminating duplicate alleles 
  dfseq8e8c = ldply(dfseq8e8b, rbind) #converting from list to data frame
  matseq8e8 = as.matrix(dfseq8e8c)
  dfseq8e8d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfseq8e8)[1]){
    testrow = matseq8e8[i,]
    dfseq8e8d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("SEQ8E8",182)
  dfseq8e8e = as.data.frame(dfseq8e8d)
  colnames(dfseq8e8e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfseq8e8f = cbind(Sample.Name,Marker,dfseq8e8e)
  ##simulating RW56
  rw56samplevec = sample(rw56namevec2,1092,replace=TRUE,prob=rw56freq)
  dfrw56 = matrix(rw56samplevec,nrow=182,ncol=6)
  dfrw56b = (apply(dfrw56,1,unique)) #eliminating duplicate alleles 
  dfrw56c = ldply(dfrw56b, rbind) #converting from list to data frame
  matrw56 = as.matrix(dfrw56c)
  dfrw56d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrw56)[1]){
    testrow = matrw56[i,]
    dfrw56d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RW56",182)
  dfrw56e = as.data.frame(dfrw56d)
  colnames(dfrw56e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrw56f = cbind(Sample.Name,Marker,dfrw56e)
  ##simulating RWDI11
  rwdi11samplevec = sample(rwdi11namevec2,1092,replace=TRUE,prob=rwdi11freq)
  dfrwdi11 = matrix(rwdi11samplevec,nrow=182,ncol=6)
  dfrwdi11b = (apply(dfrwdi11,1,unique)) #eliminating duplicate alleles 
  dfrwdi11c = ldply(dfrwdi11b, rbind) #converting from list to data frame
  matrwdi11 = as.matrix(dfrwdi11c)
  dfrwdi11d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrwdi11)[1]){
    testrow = matrwdi11[i,]
    dfrwdi11d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RWDI11",182)
  dfrwdi11e = as.data.frame(dfrwdi11d)
  colnames(dfrwdi11e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrwdi11f = cbind(Sample.Name,Marker,dfrwdi11e)
  ##combining markers, exporting and importing files
  allmarkers=rbind(dfseq18d73f,dfrw39f,dfrw28f,dfseq8e8f,dfrw56f,dfrwdi11f)
  #deleting one allele per individual
  rowsel = rep(c(0,182,364,546,728,910),30)
  rowsel[181:182] = c(0,182)
  levelsa = c(levels(allmarkers$Allele.1),levels(allmarkers$Allele.2),levels(allmarkers$Allele.3),levels(allmarkers$Allele.4),levels(allmarkers$Allele.5),levels(allmarkers$Allele.6))
  levelsb=unique(levelsa)
  allmarkers$Allele.1 = factor(allmarkers$Allele.1,levels=levelsb)
  allmarkers$Allele.2 = factor(allmarkers$Allele.2,levels=levelsb)
  allmarkers$Allele.3 = factor(allmarkers$Allele.3,levels=levelsb)
  allmarkers$Allele.4 = factor(allmarkers$Allele.4,levels=levelsb)
  allmarkers$Allele.5 = factor(allmarkers$Allele.5,levels=levelsb)
  allmarkers$Allele.6 = factor(allmarkers$Allele.6,levels=levelsb)
  allmarkers2 = allmarkers
  for(i in 1:182){
    allmarkers2[i+rowsel[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  write.table(allmarkers2,"/Users/Lakshmi/Desktop/Analysis/simdata.txt",sep = "\t",quote = FALSE, row.names=FALSE)
  simdata = read.GeneMapper("/Users/Lakshmi/Desktop/Analysis/simdata.txt")
  Usatnts(simdata) <- c(3,4,4,2,4,2)
  ploidycol = rep(6,length(Samples(simdata)))
  ploidymat = as.matrix(cbind(ploidycol,ploidycol,ploidycol,ploidycol,ploidycol,ploidycol))
  colnames(ploidymat) = c("SEQ18D73","RW39","RW28","SEQ8E8","RW56","RWDI11") 
  Ploidies(simdata) = ploidymat
  simdist = meandistance.matrix(simdata,distmetric = Bruvo.distance)
  simclones = assignClones(simdist,threshold=0.2)
  matches = length(simclones)-length(unique(simclones))
  matchvec[n] = matches
  n = n+1
}

# Testing: 2 null alleles per sample -----
matchvec = rep(190,100)
n=1
i=1
for(i in 1:100){
  ##simulating seq18d73
  seq18d73samplevec = sample(seq18d73namevec2,1092,replace=TRUE,prob=seq18d73freq)
  dfseq18d73 = matrix(seq18d73samplevec,nrow=182,ncol=6)
  dfseq18d73b = (apply(dfseq18d73,1,unique)) #eliminating duplicate alleles 
  dfseq18d73c = ldply(dfseq18d73b, rbind) #converting from list to data frame
  matseq18d73 = as.matrix(dfseq18d73c)
  dfseq18d73d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfseq18d73)[1]){
    testrow = matseq18d73[i,]
    dfseq18d73d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("SEQ18D73",182)
  dfseq18d73e = as.data.frame(dfseq18d73d)
  colnames(dfseq18d73e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfseq18d73f = cbind(Sample.Name,Marker,dfseq18d73e)
  ##simulating RW39
  rw39samplevec = sample(rw39namevec2,1092,replace=TRUE,prob=rw39freq)
  dfrw39 = matrix(rw39samplevec,nrow=182,ncol=6)
  dfrw39b = (apply(dfrw39,1,unique)) #eliminating duplicate alleles 
  dfrw39c = ldply(dfrw39b, rbind) #converting from list to data frame
  matrw39 = as.matrix(dfrw39c)
  dfrw39d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrw39)[1]){
    testrow = matrw39[i,]
    dfrw39d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RW39",182)
  dfrw39e = as.data.frame(dfrw39d)
  colnames(dfrw39e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrw39f = cbind(Sample.Name,Marker,dfrw39e)
  ##simulating RW28
  rw28samplevec = sample(rw28namevec2,1092,replace=TRUE,prob=rw28freq)
  dfrw28 = matrix(rw28samplevec,nrow=182,ncol=6)
  dfrw28b = (apply(dfrw28,1,unique)) #eliminating duplicate alleles 
  dfrw28c = ldply(dfrw28b, rbind) #converting from list to data frame
  matrw28 = as.matrix(dfrw28c)
  dfrw28d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrw28)[1]){
    testrow = matrw28[i,]
    dfrw28d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RW28",182)
  dfrw28e = as.data.frame(dfrw28d)
  colnames(dfrw28e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrw28f = cbind(Sample.Name,Marker,dfrw28e)
  ##simulating SEQ8E8
  seq8e8samplevec = sample(seq8e8namevec2,1092,replace=TRUE,prob=seq8e8freq)
  dfseq8e8 = matrix(seq8e8samplevec,nrow=182,ncol=6)
  dfseq8e8b = (apply(dfseq8e8,1,unique)) #eliminating duplicate alleles 
  dfseq8e8c = ldply(dfseq8e8b, rbind) #converting from list to data frame
  matseq8e8 = as.matrix(dfseq8e8c)
  dfseq8e8d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfseq8e8)[1]){
    testrow = matseq8e8[i,]
    dfseq8e8d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("SEQ8E8",182)
  dfseq8e8e = as.data.frame(dfseq8e8d)
  colnames(dfseq8e8e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfseq8e8f = cbind(Sample.Name,Marker,dfseq8e8e)
  ##simulating RW56
  rw56samplevec = sample(rw56namevec2,1092,replace=TRUE,prob=rw56freq)
  dfrw56 = matrix(rw56samplevec,nrow=182,ncol=6)
  dfrw56b = (apply(dfrw56,1,unique)) #eliminating duplicate alleles 
  dfrw56c = ldply(dfrw56b, rbind) #converting from list to data frame
  matrw56 = as.matrix(dfrw56c)
  dfrw56d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrw56)[1]){
    testrow = matrw56[i,]
    dfrw56d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RW56",182)
  dfrw56e = as.data.frame(dfrw56d)
  colnames(dfrw56e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrw56f = cbind(Sample.Name,Marker,dfrw56e)
  ##simulating RWDI11
  rwdi11samplevec = sample(rwdi11namevec2,1092,replace=TRUE,prob=rwdi11freq)
  dfrwdi11 = matrix(rwdi11samplevec,nrow=182,ncol=6)
  dfrwdi11b = (apply(dfrwdi11,1,unique)) #eliminating duplicate alleles 
  dfrwdi11c = ldply(dfrwdi11b, rbind) #converting from list to data frame
  matrwdi11 = as.matrix(dfrwdi11c)
  dfrwdi11d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrwdi11)[1]){
    testrow = matrwdi11[i,]
    dfrwdi11d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RWDI11",182)
  dfrwdi11e = as.data.frame(dfrwdi11d)
  colnames(dfrwdi11e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrwdi11f = cbind(Sample.Name,Marker,dfrwdi11e)
  ##combining markers, exporting and importing files
  allmarkers=rbind(dfseq18d73f,dfrw39f,dfrw28f,dfseq8e8f,dfrw56f,dfrwdi11f)
  #deleting one allele per individual
  rowsel = rep(c(0,182,364,546,728,910),30)
  rowsel[181:182] = c(0,182)
  levelsa = c(levels(allmarkers$Allele.1),levels(allmarkers$Allele.2),levels(allmarkers$Allele.3),levels(allmarkers$Allele.4),levels(allmarkers$Allele.5),levels(allmarkers$Allele.6))
  levelsb=unique(levelsa)
  allmarkers$Allele.1 = factor(allmarkers$Allele.1,levels=levelsb)
  allmarkers$Allele.2 = factor(allmarkers$Allele.2,levels=levelsb)
  allmarkers$Allele.3 = factor(allmarkers$Allele.3,levels=levelsb)
  allmarkers$Allele.4 = factor(allmarkers$Allele.4,levels=levelsb)
  allmarkers$Allele.5 = factor(allmarkers$Allele.5,levels=levelsb)
  allmarkers$Allele.6 = factor(allmarkers$Allele.6,levels=levelsb)
  allmarkers2 = allmarkers
  ## randomly eliminate first allele
  for(i in 1:182){
    allmarkers2[i+rowsel[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  rowsel2 = rep(c(182,364,546,728,910,0),30)
  rowsel2[181:182] = c(182,364)
  ##random eliminate second allele
  for(i in 1:182){
    allmarkers2[i+rowsel2[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel2[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel2[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  write.table(allmarkers2,"/Users/Lakshmi/Desktop/Analysis/simdata.txt",sep = "\t",quote = FALSE, row.names=FALSE)
  simdata = read.GeneMapper("/Users/Lakshmi/Desktop/Analysis/simdata.txt")
  Usatnts(simdata) <- c(3,4,4,2,4,2)
  ploidycol = rep(6,length(Samples(simdata)))
  ploidymat = as.matrix(cbind(ploidycol,ploidycol,ploidycol,ploidycol,ploidycol,ploidycol))
  colnames(ploidymat) = c("SEQ18D73","RW39","RW28","SEQ8E8","RW56","RWDI11") 
  Ploidies(simdata) = ploidymat
  simdist = meandistance.matrix(simdata,distmetric = Bruvo.distance)
  simclones = assignClones(simdist,threshold=0.2)
  matches = length(simclones)-length(unique(simclones))
  matchvec[n] = matches
  n = n+1
}

# Testing: 3 null alleles per sample ----- no repeats
matchvec = rep(190,100)
n=1
i=1
##change to i in 1:100
for(i in 1:100){
  ##simulating seq18d73
  seq18d73samplevec = sample(seq18d73namevec2,1092,replace=TRUE,prob=seq18d73freq)
  dfseq18d73 = matrix(seq18d73samplevec,nrow=182,ncol=6)
  dfseq18d73b = (apply(dfseq18d73,1,unique)) #eliminating duplicate alleles 
  dfseq18d73c = ldply(dfseq18d73b, rbind) #converting from list to data frame
  matseq18d73 = as.matrix(dfseq18d73c)
  dfseq18d73d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfseq18d73)[1]){
    testrow = matseq18d73[i,]
    dfseq18d73d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("SEQ18D73",182)
  dfseq18d73e = as.data.frame(dfseq18d73d)
  colnames(dfseq18d73e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfseq18d73f = cbind(Sample.Name,Marker,dfseq18d73e)
  ##simulating RW39
  rw39samplevec = sample(rw39namevec2,1092,replace=TRUE,prob=rw39freq)
  dfrw39 = matrix(rw39samplevec,nrow=182,ncol=6)
  dfrw39b = (apply(dfrw39,1,unique)) #eliminating duplicate alleles 
  dfrw39c = ldply(dfrw39b, rbind) #converting from list to data frame
  matrw39 = as.matrix(dfrw39c)
  dfrw39d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrw39)[1]){
    testrow = matrw39[i,]
    dfrw39d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RW39",182)
  dfrw39e = as.data.frame(dfrw39d)
  colnames(dfrw39e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrw39f = cbind(Sample.Name,Marker,dfrw39e)
  ##simulating RW28
  rw28samplevec = sample(rw28namevec2,1092,replace=TRUE,prob=rw28freq)
  dfrw28 = matrix(rw28samplevec,nrow=182,ncol=6)
  dfrw28b = (apply(dfrw28,1,unique)) #eliminating duplicate alleles 
  dfrw28c = ldply(dfrw28b, rbind) #converting from list to data frame
  matrw28 = as.matrix(dfrw28c)
  dfrw28d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrw28)[1]){
    testrow = matrw28[i,]
    dfrw28d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RW28",182)
  dfrw28e = as.data.frame(dfrw28d)
  colnames(dfrw28e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrw28f = cbind(Sample.Name,Marker,dfrw28e)
  ##simulating SEQ8E8
  seq8e8samplevec = sample(seq8e8namevec2,1092,replace=TRUE,prob=seq8e8freq)
  dfseq8e8 = matrix(seq8e8samplevec,nrow=182,ncol=6)
  dfseq8e8b = (apply(dfseq8e8,1,unique)) #eliminating duplicate alleles 
  dfseq8e8c = ldply(dfseq8e8b, rbind) #converting from list to data frame
  matseq8e8 = as.matrix(dfseq8e8c)
  dfseq8e8d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfseq8e8)[1]){
    testrow = matseq8e8[i,]
    dfseq8e8d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("SEQ8E8",182)
  dfseq8e8e = as.data.frame(dfseq8e8d)
  colnames(dfseq8e8e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfseq8e8f = cbind(Sample.Name,Marker,dfseq8e8e)
  ##simulating RW56
  rw56samplevec = sample(rw56namevec2,1092,replace=TRUE,prob=rw56freq)
  dfrw56 = matrix(rw56samplevec,nrow=182,ncol=6)
  dfrw56b = (apply(dfrw56,1,unique)) #eliminating duplicate alleles 
  dfrw56c = ldply(dfrw56b, rbind) #converting from list to data frame
  matrw56 = as.matrix(dfrw56c)
  dfrw56d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrw56)[1]){
    testrow = matrw56[i,]
    dfrw56d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RW56",182)
  dfrw56e = as.data.frame(dfrw56d)
  colnames(dfrw56e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrw56f = cbind(Sample.Name,Marker,dfrw56e)
  ##simulating RWDI11
  rwdi11samplevec = sample(rwdi11namevec2,1092,replace=TRUE,prob=rwdi11freq)
  dfrwdi11 = matrix(rwdi11samplevec,nrow=182,ncol=6)
  dfrwdi11b = (apply(dfrwdi11,1,unique)) #eliminating duplicate alleles 
  dfrwdi11c = ldply(dfrwdi11b, rbind) #converting from list to data frame
  matrwdi11 = as.matrix(dfrwdi11c)
  dfrwdi11d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrwdi11)[1]){
    testrow = matrwdi11[i,]
    dfrwdi11d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RWDI11",182)
  dfrwdi11e = as.data.frame(dfrwdi11d)
  colnames(dfrwdi11e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrwdi11f = cbind(Sample.Name,Marker,dfrwdi11e)
  ##combining markers, exporting and importing files
  allmarkers=rbind(dfseq18d73f,dfrw39f,dfrw28f,dfseq8e8f,dfrw56f,dfrwdi11f)
  #deleting one allele per individual
  levelsa = c(levels(allmarkers$Allele.1),levels(allmarkers$Allele.2),levels(allmarkers$Allele.3),levels(allmarkers$Allele.4),levels(allmarkers$Allele.5),levels(allmarkers$Allele.6))
  levelsb=unique(levelsa)
  allmarkers$Allele.1 = factor(allmarkers$Allele.1,levels=levelsb)
  allmarkers$Allele.2 = factor(allmarkers$Allele.2,levels=levelsb)
  allmarkers$Allele.3 = factor(allmarkers$Allele.3,levels=levelsb)
  allmarkers$Allele.4 = factor(allmarkers$Allele.4,levels=levelsb)
  allmarkers$Allele.5 = factor(allmarkers$Allele.5,levels=levelsb)
  allmarkers$Allele.6 = factor(allmarkers$Allele.6,levels=levelsb)
  allmarkers2 = allmarkers
  ## randomly eliminate first allele
  rowsel = rep(c(0,182,364,546,728,910),30)
  rowsel[181:182] = c(0,182)
  for(i in 1:182){
    allmarkers2[i+rowsel[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##random eliminate second allele
  rowsel2 = rep(c(182,364,546,728,910,0),30)
  rowsel2[181:182] = c(182,364)
  for(i in 1:182){
    allmarkers2[i+rowsel2[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel2[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel2[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##random eliminate third allele
  rowsel3 = rep(c(364,546,728,910,0,182),30)
  rowsel3[181:182] = c(364,546)
  for(i in 1:182){
    allmarkers2[i+rowsel3[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel3[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel3[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  write.table(allmarkers2,"/Users/Lakshmi/Desktop/Analysis/simdata.txt",sep = "\t",quote = FALSE, row.names=FALSE)
  simdata = read.GeneMapper("/Users/Lakshmi/Desktop/Analysis/simdata.txt")
  Usatnts(simdata) <- c(3,4,4,2,4,2)
  ploidycol = rep(6,length(Samples(simdata)))
  ploidymat = as.matrix(cbind(ploidycol,ploidycol,ploidycol,ploidycol,ploidycol,ploidycol))
  colnames(ploidymat) = c("SEQ18D73","RW39","RW28","SEQ8E8","RW56","RWDI11") 
  Ploidies(simdata) = ploidymat
  simdist = meandistance.matrix(simdata,distmetric = Bruvo.distance)
  simclones = assignClones(simdist,threshold=0.2)
  matches = length(simclones)-length(unique(simclones))
  matchvec[n] = matches
  n = n+1
}

# Testing: 4 null alleles per sample -----
matchvec = rep(190,99)
n=1
i=1
##change to i in 1:100
for(i in 1:99){
  ##simulating seq18d73
  seq18d73samplevec = sample(seq18d73namevec2,1092,replace=TRUE,prob=seq18d73freq)
  dfseq18d73 = matrix(seq18d73samplevec,nrow=182,ncol=6)
  dfseq18d73b = (apply(dfseq18d73,1,unique)) #eliminating duplicate alleles 
  dfseq18d73c = ldply(dfseq18d73b, rbind) #converting from list to data frame
  matseq18d73 = as.matrix(dfseq18d73c)
  dfseq18d73d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfseq18d73)[1]){
    testrow = matseq18d73[i,]
    dfseq18d73d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("SEQ18D73",182)
  dfseq18d73e = as.data.frame(dfseq18d73d)
  colnames(dfseq18d73e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfseq18d73f = cbind(Sample.Name,Marker,dfseq18d73e)
  ##simulating RW39
  rw39samplevec = sample(rw39namevec2,1092,replace=TRUE,prob=rw39freq)
  dfrw39 = matrix(rw39samplevec,nrow=182,ncol=6)
  dfrw39b = (apply(dfrw39,1,unique)) #eliminating duplicate alleles 
  dfrw39c = ldply(dfrw39b, rbind) #converting from list to data frame
  matrw39 = as.matrix(dfrw39c)
  dfrw39d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrw39)[1]){
    testrow = matrw39[i,]
    dfrw39d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RW39",182)
  dfrw39e = as.data.frame(dfrw39d)
  colnames(dfrw39e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrw39f = cbind(Sample.Name,Marker,dfrw39e)
  ##simulating RW28
  rw28samplevec = sample(rw28namevec2,1092,replace=TRUE,prob=rw28freq)
  dfrw28 = matrix(rw28samplevec,nrow=182,ncol=6)
  dfrw28b = (apply(dfrw28,1,unique)) #eliminating duplicate alleles 
  dfrw28c = ldply(dfrw28b, rbind) #converting from list to data frame
  matrw28 = as.matrix(dfrw28c)
  dfrw28d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrw28)[1]){
    testrow = matrw28[i,]
    dfrw28d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RW28",182)
  dfrw28e = as.data.frame(dfrw28d)
  colnames(dfrw28e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrw28f = cbind(Sample.Name,Marker,dfrw28e)
  ##simulating SEQ8E8
  seq8e8samplevec = sample(seq8e8namevec2,1092,replace=TRUE,prob=seq8e8freq)
  dfseq8e8 = matrix(seq8e8samplevec,nrow=182,ncol=6)
  dfseq8e8b = (apply(dfseq8e8,1,unique)) #eliminating duplicate alleles 
  dfseq8e8c = ldply(dfseq8e8b, rbind) #converting from list to data frame
  matseq8e8 = as.matrix(dfseq8e8c)
  dfseq8e8d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfseq8e8)[1]){
    testrow = matseq8e8[i,]
    dfseq8e8d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("SEQ8E8",182)
  dfseq8e8e = as.data.frame(dfseq8e8d)
  colnames(dfseq8e8e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfseq8e8f = cbind(Sample.Name,Marker,dfseq8e8e)
  ##simulating RW56
  rw56samplevec = sample(rw56namevec2,1092,replace=TRUE,prob=rw56freq)
  dfrw56 = matrix(rw56samplevec,nrow=182,ncol=6)
  dfrw56b = (apply(dfrw56,1,unique)) #eliminating duplicate alleles 
  dfrw56c = ldply(dfrw56b, rbind) #converting from list to data frame
  matrw56 = as.matrix(dfrw56c)
  dfrw56d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrw56)[1]){
    testrow = matrw56[i,]
    dfrw56d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RW56",182)
  dfrw56e = as.data.frame(dfrw56d)
  colnames(dfrw56e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrw56f = cbind(Sample.Name,Marker,dfrw56e)
  ##simulating RWDI11
  rwdi11samplevec = sample(rwdi11namevec2,1092,replace=TRUE,prob=rwdi11freq)
  dfrwdi11 = matrix(rwdi11samplevec,nrow=182,ncol=6)
  dfrwdi11b = (apply(dfrwdi11,1,unique)) #eliminating duplicate alleles 
  dfrwdi11c = ldply(dfrwdi11b, rbind) #converting from list to data frame
  matrwdi11 = as.matrix(dfrwdi11c)
  dfrwdi11d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrwdi11)[1]){
    testrow = matrwdi11[i,]
    dfrwdi11d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RWDI11",182)
  dfrwdi11e = as.data.frame(dfrwdi11d)
  colnames(dfrwdi11e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrwdi11f = cbind(Sample.Name,Marker,dfrwdi11e)
  ##combining markers, exporting and importing files
  allmarkers=rbind(dfseq18d73f,dfrw39f,dfrw28f,dfseq8e8f,dfrw56f,dfrwdi11f)
  #deleting one allele per individual
  levelsa = c(levels(allmarkers$Allele.1),levels(allmarkers$Allele.2),levels(allmarkers$Allele.3),levels(allmarkers$Allele.4),levels(allmarkers$Allele.5),levels(allmarkers$Allele.6))
  levelsb=unique(levelsa)
  allmarkers$Allele.1 = factor(allmarkers$Allele.1,levels=levelsb)
  allmarkers$Allele.2 = factor(allmarkers$Allele.2,levels=levelsb)
  allmarkers$Allele.3 = factor(allmarkers$Allele.3,levels=levelsb)
  allmarkers$Allele.4 = factor(allmarkers$Allele.4,levels=levelsb)
  allmarkers$Allele.5 = factor(allmarkers$Allele.5,levels=levelsb)
  allmarkers$Allele.6 = factor(allmarkers$Allele.6,levels=levelsb)
  allmarkers2 = allmarkers
  ## randomly eliminate first allele
  rowsel = rep(c(0,182,364,546,728,910),30)
  rowsel[181:182] = c(0,182)
  for(i in 1:182){
    allmarkers2[i+rowsel[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##random eliminate second allele
  rowsel2 = rep(c(182,364,546,728,910,0),30)
  rowsel2[181:182] = c(182,364)
  for(i in 1:182){
    allmarkers2[i+rowsel2[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel2[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel2[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##random eliminate third allele
  rowsel3 = rep(c(364,546,728,910,0,182),30)
  rowsel3[181:182] = c(364,546)
  for(i in 1:182){
    allmarkers2[i+rowsel3[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel3[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel3[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate fourth allele
  rowsel4 = rep(c(546,728,910,0,182,364),30)
  rowsel4[181:182] = c(546,728)
  for(i in 1:182){
    allmarkers2[i+rowsel4[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel4[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel4[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  write.table(allmarkers2,"/Users/Lakshmi/Desktop/Analysis/simdata.txt",sep = "\t",quote = FALSE, row.names=FALSE)
  simdata = read.GeneMapper("/Users/Lakshmi/Desktop/Analysis/simdata.txt")
  Usatnts(simdata) <- c(3,4,4,2,4,2)
  ploidycol = rep(6,length(Samples(simdata)))
  ploidymat = as.matrix(cbind(ploidycol,ploidycol,ploidycol,ploidycol,ploidycol,ploidycol))
  colnames(ploidymat) = c("SEQ18D73","RW39","RW28","SEQ8E8","RW56","RWDI11") 
  Ploidies(simdata) = ploidymat
  simdist = meandistance.matrix(simdata,distmetric = Bruvo.distance)
  simclones = assignClones(simdist,threshold=0.2)
  matches = length(simclones)-length(unique(simclones))
  matchvec[n] = matches
  n = n+1
}

# Testing: 5 null alleles per sample -----
matchvec = rep(190,99)
n=1
i=1

for(i in 1:99){
  ##simulating seq18d73
  seq18d73samplevec = sample(seq18d73namevec2,1092,replace=TRUE,prob=seq18d73freq)
  dfseq18d73 = matrix(seq18d73samplevec,nrow=182,ncol=6)
  dfseq18d73b = (apply(dfseq18d73,1,unique)) #eliminating duplicate alleles 
  dfseq18d73c = ldply(dfseq18d73b, rbind) #converting from list to data frame
  matseq18d73 = as.matrix(dfseq18d73c)
  dfseq18d73d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfseq18d73)[1]){
    testrow = matseq18d73[i,]
    dfseq18d73d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("SEQ18D73",182)
  dfseq18d73e = as.data.frame(dfseq18d73d)
  colnames(dfseq18d73e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfseq18d73f = cbind(Sample.Name,Marker,dfseq18d73e)
  ##simulating RW39
  rw39samplevec = sample(rw39namevec2,1092,replace=TRUE,prob=rw39freq)
  dfrw39 = matrix(rw39samplevec,nrow=182,ncol=6)
  dfrw39b = (apply(dfrw39,1,unique)) #eliminating duplicate alleles 
  dfrw39c = ldply(dfrw39b, rbind) #converting from list to data frame
  matrw39 = as.matrix(dfrw39c)
  dfrw39d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrw39)[1]){
    testrow = matrw39[i,]
    dfrw39d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RW39",182)
  dfrw39e = as.data.frame(dfrw39d)
  colnames(dfrw39e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrw39f = cbind(Sample.Name,Marker,dfrw39e)
  ##simulating RW28
  rw28samplevec = sample(rw28namevec2,1092,replace=TRUE,prob=rw28freq)
  dfrw28 = matrix(rw28samplevec,nrow=182,ncol=6)
  dfrw28b = (apply(dfrw28,1,unique)) #eliminating duplicate alleles 
  dfrw28c = ldply(dfrw28b, rbind) #converting from list to data frame
  matrw28 = as.matrix(dfrw28c)
  dfrw28d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrw28)[1]){
    testrow = matrw28[i,]
    dfrw28d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RW28",182)
  dfrw28e = as.data.frame(dfrw28d)
  colnames(dfrw28e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrw28f = cbind(Sample.Name,Marker,dfrw28e)
  ##simulating SEQ8E8
  seq8e8samplevec = sample(seq8e8namevec2,1092,replace=TRUE,prob=seq8e8freq)
  dfseq8e8 = matrix(seq8e8samplevec,nrow=182,ncol=6)
  dfseq8e8b = (apply(dfseq8e8,1,unique)) #eliminating duplicate alleles 
  dfseq8e8c = ldply(dfseq8e8b, rbind) #converting from list to data frame
  matseq8e8 = as.matrix(dfseq8e8c)
  dfseq8e8d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfseq8e8)[1]){
    testrow = matseq8e8[i,]
    dfseq8e8d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("SEQ8E8",182)
  dfseq8e8e = as.data.frame(dfseq8e8d)
  colnames(dfseq8e8e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfseq8e8f = cbind(Sample.Name,Marker,dfseq8e8e)
  ##simulating RW56
  rw56samplevec = sample(rw56namevec2,1092,replace=TRUE,prob=rw56freq)
  dfrw56 = matrix(rw56samplevec,nrow=182,ncol=6)
  dfrw56b = (apply(dfrw56,1,unique)) #eliminating duplicate alleles 
  dfrw56c = ldply(dfrw56b, rbind) #converting from list to data frame
  matrw56 = as.matrix(dfrw56c)
  dfrw56d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrw56)[1]){
    testrow = matrw56[i,]
    dfrw56d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RW56",182)
  dfrw56e = as.data.frame(dfrw56d)
  colnames(dfrw56e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrw56f = cbind(Sample.Name,Marker,dfrw56e)
  ##simulating RWDI11
  rwdi11samplevec = sample(rwdi11namevec2,1092,replace=TRUE,prob=rwdi11freq)
  dfrwdi11 = matrix(rwdi11samplevec,nrow=182,ncol=6)
  dfrwdi11b = (apply(dfrwdi11,1,unique)) #eliminating duplicate alleles 
  dfrwdi11c = ldply(dfrwdi11b, rbind) #converting from list to data frame
  matrwdi11 = as.matrix(dfrwdi11c)
  dfrwdi11d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrwdi11)[1]){
    testrow = matrwdi11[i,]
    dfrwdi11d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RWDI11",182)
  dfrwdi11e = as.data.frame(dfrwdi11d)
  colnames(dfrwdi11e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrwdi11f = cbind(Sample.Name,Marker,dfrwdi11e)
  ##combining markers, exporting and importing files
  allmarkers=rbind(dfseq18d73f,dfrw39f,dfrw28f,dfseq8e8f,dfrw56f,dfrwdi11f)
  #deleting one allele per individual
  levelsa = c(levels(allmarkers$Allele.1),levels(allmarkers$Allele.2),levels(allmarkers$Allele.3),levels(allmarkers$Allele.4),levels(allmarkers$Allele.5),levels(allmarkers$Allele.6))
  levelsb=unique(levelsa)
  allmarkers$Allele.1 = factor(allmarkers$Allele.1,levels=levelsb)
  allmarkers$Allele.2 = factor(allmarkers$Allele.2,levels=levelsb)
  allmarkers$Allele.3 = factor(allmarkers$Allele.3,levels=levelsb)
  allmarkers$Allele.4 = factor(allmarkers$Allele.4,levels=levelsb)
  allmarkers$Allele.5 = factor(allmarkers$Allele.5,levels=levelsb)
  allmarkers$Allele.6 = factor(allmarkers$Allele.6,levels=levelsb)
  allmarkers2 = allmarkers
  ## randomly eliminate first allele
  rowsel = rep(c(0,182,364,546,728,910),30)
  rowsel[181:182] = c(0,182)
  for(i in 1:182){
    allmarkers2[i+rowsel[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate second allele
  rowsel2 = rep(c(182,364,546,728,910,0),30)
  rowsel2[181:182] = c(182,364)
  for(i in 1:182){
    allmarkers2[i+rowsel2[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel2[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel2[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate third allele
  rowsel3 = rep(c(364,546,728,910,0,182),30)
  rowsel3[181:182] = c(364,546)
  for(i in 1:182){
    allmarkers2[i+rowsel3[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel3[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel3[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate fourth allele
  rowsel4 = rep(c(546,728,910,0,182,364),30)
  rowsel4[181:182] = c(546,728)
  for(i in 1:182){
    allmarkers2[i+rowsel4[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel4[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel4[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate fifth allele
  rowsel5 = rep(c(728,910,0,182,364,546),30)
  rowsel5[181:182] = c(728,910)
  for(i in 1:182){
    allmarkers2[i+rowsel5[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel5[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel5[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  write.table(allmarkers2,"/Users/Lakshmi/Desktop/Analysis/simdata.txt",sep = "\t",quote = FALSE, row.names=FALSE)
  simdata = read.GeneMapper("/Users/Lakshmi/Desktop/Analysis/simdata.txt")
  Usatnts(simdata) <- c(3,4,4,2,4,2)
  ploidycol = rep(6,length(Samples(simdata)))
  ploidymat = as.matrix(cbind(ploidycol,ploidycol,ploidycol,ploidycol,ploidycol,ploidycol))
  colnames(ploidymat) = c("SEQ18D73","RW39","RW28","SEQ8E8","RW56","RWDI11") 
  Ploidies(simdata) = ploidymat
  simdist = meandistance.matrix(simdata,distmetric = Bruvo.distance)
  simclones = assignClones(simdist,threshold=0.2)
  matches = length(simclones)-length(unique(simclones))
  matchvec[n] = matches
  n = n+1
}

# Testing: 6 null alleles per sample -----
matchvec = rep(190,99)
n=1
i=1

for(i in 1:1){
  ##simulating seq18d73
  seq18d73samplevec = sample(seq18d73namevec2,1092,replace=TRUE,prob=seq18d73freq)
  dfseq18d73 = matrix(seq18d73samplevec,nrow=182,ncol=6)
  dfseq18d73b = (apply(dfseq18d73,1,unique)) #eliminating duplicate alleles 
  dfseq18d73c = ldply(dfseq18d73b, rbind) #converting from list to data frame
  matseq18d73 = as.matrix(dfseq18d73c)
  dfseq18d73d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfseq18d73)[1]){
    testrow = matseq18d73[i,]
    dfseq18d73d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("SEQ18D73",182)
  dfseq18d73e = as.data.frame(dfseq18d73d)
  colnames(dfseq18d73e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfseq18d73f = cbind(Sample.Name,Marker,dfseq18d73e)
  ##simulating RW39
  rw39samplevec = sample(rw39namevec2,1092,replace=TRUE,prob=rw39freq)
  dfrw39 = matrix(rw39samplevec,nrow=182,ncol=6)
  dfrw39b = (apply(dfrw39,1,unique)) #eliminating duplicate alleles 
  dfrw39c = ldply(dfrw39b, rbind) #converting from list to data frame
  matrw39 = as.matrix(dfrw39c)
  dfrw39d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrw39)[1]){
    testrow = matrw39[i,]
    dfrw39d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RW39",182)
  dfrw39e = as.data.frame(dfrw39d)
  colnames(dfrw39e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrw39f = cbind(Sample.Name,Marker,dfrw39e)
  ##simulating RW28
  rw28samplevec = sample(rw28namevec2,1092,replace=TRUE,prob=rw28freq)
  dfrw28 = matrix(rw28samplevec,nrow=182,ncol=6)
  dfrw28b = (apply(dfrw28,1,unique)) #eliminating duplicate alleles 
  dfrw28c = ldply(dfrw28b, rbind) #converting from list to data frame
  matrw28 = as.matrix(dfrw28c)
  dfrw28d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrw28)[1]){
    testrow = matrw28[i,]
    dfrw28d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RW28",182)
  dfrw28e = as.data.frame(dfrw28d)
  colnames(dfrw28e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrw28f = cbind(Sample.Name,Marker,dfrw28e)
  ##simulating SEQ8E8
  seq8e8samplevec = sample(seq8e8namevec2,1092,replace=TRUE,prob=seq8e8freq)
  dfseq8e8 = matrix(seq8e8samplevec,nrow=182,ncol=6)
  dfseq8e8b = (apply(dfseq8e8,1,unique)) #eliminating duplicate alleles 
  dfseq8e8c = ldply(dfseq8e8b, rbind) #converting from list to data frame
  matseq8e8 = as.matrix(dfseq8e8c)
  dfseq8e8d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfseq8e8)[1]){
    testrow = matseq8e8[i,]
    dfseq8e8d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("SEQ8E8",182)
  dfseq8e8e = as.data.frame(dfseq8e8d)
  colnames(dfseq8e8e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfseq8e8f = cbind(Sample.Name,Marker,dfseq8e8e)
  ##simulating RW56
  rw56samplevec = sample(rw56namevec2,1092,replace=TRUE,prob=rw56freq)
  dfrw56 = matrix(rw56samplevec,nrow=182,ncol=6)
  dfrw56b = (apply(dfrw56,1,unique)) #eliminating duplicate alleles 
  dfrw56c = ldply(dfrw56b, rbind) #converting from list to data frame
  matrw56 = as.matrix(dfrw56c)
  dfrw56d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrw56)[1]){
    testrow = matrw56[i,]
    dfrw56d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RW56",182)
  dfrw56e = as.data.frame(dfrw56d)
  colnames(dfrw56e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrw56f = cbind(Sample.Name,Marker,dfrw56e)
  ##simulating RWDI11
  rwdi11samplevec = sample(rwdi11namevec2,1092,replace=TRUE,prob=rwdi11freq)
  dfrwdi11 = matrix(rwdi11samplevec,nrow=182,ncol=6)
  dfrwdi11b = (apply(dfrwdi11,1,unique)) #eliminating duplicate alleles 
  dfrwdi11c = ldply(dfrwdi11b, rbind) #converting from list to data frame
  matrwdi11 = as.matrix(dfrwdi11c)
  dfrwdi11d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrwdi11)[1]){
    testrow = matrwdi11[i,]
    dfrwdi11d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RWDI11",182)
  dfrwdi11e = as.data.frame(dfrwdi11d)
  colnames(dfrwdi11e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrwdi11f = cbind(Sample.Name,Marker,dfrwdi11e)
  ##combining markers, exporting and importing files
  allmarkers=rbind(dfseq18d73f,dfrw39f,dfrw28f,dfseq8e8f,dfrw56f,dfrwdi11f)
  #deleting one allele per individual
  levelsa = c(levels(allmarkers$Allele.1),levels(allmarkers$Allele.2),levels(allmarkers$Allele.3),levels(allmarkers$Allele.4),levels(allmarkers$Allele.5),levels(allmarkers$Allele.6))
  levelsb=unique(levelsa)
  allmarkers$Allele.1 = factor(allmarkers$Allele.1,levels=levelsb)
  allmarkers$Allele.2 = factor(allmarkers$Allele.2,levels=levelsb)
  allmarkers$Allele.3 = factor(allmarkers$Allele.3,levels=levelsb)
  allmarkers$Allele.4 = factor(allmarkers$Allele.4,levels=levelsb)
  allmarkers$Allele.5 = factor(allmarkers$Allele.5,levels=levelsb)
  allmarkers$Allele.6 = factor(allmarkers$Allele.6,levels=levelsb)
  allmarkers2 = allmarkers
  ## randomly eliminate first allele
  rowsel = rep(c(0,182,364,546,728,910),30)
  rowsel[181:182] = c(0,182)
  for(i in 1:182){
    allmarkers2[i+rowsel[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate second allele
  rowsel2 = rep(c(182,364,546,728,910,0),30)
  rowsel2[181:182] = c(182,364)
  for(i in 1:182){
    allmarkers2[i+rowsel2[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel2[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel2[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate third allele
  rowsel3 = rep(c(364,546,728,910,0,182),30)
  rowsel3[181:182] = c(364,546)
  for(i in 1:182){
    allmarkers2[i+rowsel3[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel3[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel3[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate fourth allele
  rowsel4 = rep(c(546,728,910,0,182,364),30)
  rowsel4[181:182] = c(546,728)
  for(i in 1:182){
    allmarkers2[i+rowsel4[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel4[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel4[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate fifth allele
  rowsel5 = rep(c(728,910,0,182,364,546),30)
  rowsel5[181:182] = c(728,910)
  for(i in 1:182){
    allmarkers2[i+rowsel5[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel5[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel5[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate sixth allele
  rowsel6 = rep(c(910,0,182,364,546,728),30)
  rowsel6[181:182] = c(910,0)
  for(i in 1:182){
    allmarkers2[i+rowsel6[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel6[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel6[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  write.table(allmarkers2,"/Users/Lakshmi/Desktop/Analysis/simdata.txt",sep = "\t",quote = FALSE, row.names=FALSE)
  simdata = read.GeneMapper("/Users/Lakshmi/Desktop/Analysis/simdata.txt")
  Usatnts(simdata) <- c(3,4,4,2,4,2)
  ploidycol = rep(6,length(Samples(simdata)))
  ploidymat = as.matrix(cbind(ploidycol,ploidycol,ploidycol,ploidycol,ploidycol,ploidycol))
  colnames(ploidymat) = c("SEQ18D73","RW39","RW28","SEQ8E8","RW56","RWDI11") 
  Ploidies(simdata) = ploidymat
  simdist = meandistance.matrix(simdata,distmetric = Bruvo.distance)
  simclones = assignClones(simdist,threshold=0.2)
  matches = length(simclones)-length(unique(simclones))
  matchvec[n] = matches
  n = n+1
}

# Testing: 7 null alleles per sample -----
matchvec = rep(190,99)
n=1
i=1

for(i in 1:99){
  ##simulating seq18d73
  seq18d73samplevec = sample(seq18d73namevec2,1092,replace=TRUE,prob=seq18d73freq)
  dfseq18d73 = matrix(seq18d73samplevec,nrow=182,ncol=6)
  dfseq18d73b = (apply(dfseq18d73,1,unique)) #eliminating duplicate alleles 
  dfseq18d73c = ldply(dfseq18d73b, rbind) #converting from list to data frame
  matseq18d73 = as.matrix(dfseq18d73c)
  dfseq18d73d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfseq18d73)[1]){
    testrow = matseq18d73[i,]
    dfseq18d73d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("SEQ18D73",182)
  dfseq18d73e = as.data.frame(dfseq18d73d)
  colnames(dfseq18d73e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfseq18d73f = cbind(Sample.Name,Marker,dfseq18d73e)
  ##simulating RW39
  rw39samplevec = sample(rw39namevec2,1092,replace=TRUE,prob=rw39freq)
  dfrw39 = matrix(rw39samplevec,nrow=182,ncol=6)
  dfrw39b = (apply(dfrw39,1,unique)) #eliminating duplicate alleles 
  dfrw39c = ldply(dfrw39b, rbind) #converting from list to data frame
  matrw39 = as.matrix(dfrw39c)
  dfrw39d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrw39)[1]){
    testrow = matrw39[i,]
    dfrw39d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RW39",182)
  dfrw39e = as.data.frame(dfrw39d)
  colnames(dfrw39e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrw39f = cbind(Sample.Name,Marker,dfrw39e)
  ##simulating RW28
  rw28samplevec = sample(rw28namevec2,1092,replace=TRUE,prob=rw28freq)
  dfrw28 = matrix(rw28samplevec,nrow=182,ncol=6)
  dfrw28b = (apply(dfrw28,1,unique)) #eliminating duplicate alleles 
  dfrw28c = ldply(dfrw28b, rbind) #converting from list to data frame
  matrw28 = as.matrix(dfrw28c)
  dfrw28d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrw28)[1]){
    testrow = matrw28[i,]
    dfrw28d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RW28",182)
  dfrw28e = as.data.frame(dfrw28d)
  colnames(dfrw28e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrw28f = cbind(Sample.Name,Marker,dfrw28e)
  ##simulating SEQ8E8
  seq8e8samplevec = sample(seq8e8namevec2,1092,replace=TRUE,prob=seq8e8freq)
  dfseq8e8 = matrix(seq8e8samplevec,nrow=182,ncol=6)
  dfseq8e8b = (apply(dfseq8e8,1,unique)) #eliminating duplicate alleles 
  dfseq8e8c = ldply(dfseq8e8b, rbind) #converting from list to data frame
  matseq8e8 = as.matrix(dfseq8e8c)
  dfseq8e8d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfseq8e8)[1]){
    testrow = matseq8e8[i,]
    dfseq8e8d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("SEQ8E8",182)
  dfseq8e8e = as.data.frame(dfseq8e8d)
  colnames(dfseq8e8e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfseq8e8f = cbind(Sample.Name,Marker,dfseq8e8e)
  ##simulating RW56
  rw56samplevec = sample(rw56namevec2,1092,replace=TRUE,prob=rw56freq)
  dfrw56 = matrix(rw56samplevec,nrow=182,ncol=6)
  dfrw56b = (apply(dfrw56,1,unique)) #eliminating duplicate alleles 
  dfrw56c = ldply(dfrw56b, rbind) #converting from list to data frame
  matrw56 = as.matrix(dfrw56c)
  dfrw56d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrw56)[1]){
    testrow = matrw56[i,]
    dfrw56d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RW56",182)
  dfrw56e = as.data.frame(dfrw56d)
  colnames(dfrw56e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrw56f = cbind(Sample.Name,Marker,dfrw56e)
  ##simulating RWDI11
  rwdi11samplevec = sample(rwdi11namevec2,1092,replace=TRUE,prob=rwdi11freq)
  dfrwdi11 = matrix(rwdi11samplevec,nrow=182,ncol=6)
  dfrwdi11b = (apply(dfrwdi11,1,unique)) #eliminating duplicate alleles 
  dfrwdi11c = ldply(dfrwdi11b, rbind) #converting from list to data frame
  matrwdi11 = as.matrix(dfrwdi11c)
  dfrwdi11d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrwdi11)[1]){
    testrow = matrwdi11[i,]
    dfrwdi11d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RWDI11",182)
  dfrwdi11e = as.data.frame(dfrwdi11d)
  colnames(dfrwdi11e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrwdi11f = cbind(Sample.Name,Marker,dfrwdi11e)
  ##combining markers, exporting and importing files
  allmarkers=rbind(dfseq18d73f,dfrw39f,dfrw28f,dfseq8e8f,dfrw56f,dfrwdi11f)
  #deleting one allele per individual
  levelsa = c(levels(allmarkers$Allele.1),levels(allmarkers$Allele.2),levels(allmarkers$Allele.3),levels(allmarkers$Allele.4),levels(allmarkers$Allele.5),levels(allmarkers$Allele.6))
  levelsb=unique(levelsa)
  allmarkers$Allele.1 = factor(allmarkers$Allele.1,levels=levelsb)
  allmarkers$Allele.2 = factor(allmarkers$Allele.2,levels=levelsb)
  allmarkers$Allele.3 = factor(allmarkers$Allele.3,levels=levelsb)
  allmarkers$Allele.4 = factor(allmarkers$Allele.4,levels=levelsb)
  allmarkers$Allele.5 = factor(allmarkers$Allele.5,levels=levelsb)
  allmarkers$Allele.6 = factor(allmarkers$Allele.6,levels=levelsb)
  allmarkers2 = allmarkers
  ## randomly eliminate first allele
  rowsel = rep(c(0,182,364,546,728,910),30)
  rowsel[181:182] = c(0,182)
  for(i in 1:182){
    allmarkers2[i+rowsel[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate second allele
  rowsel2 = rep(c(182,364,546,728,910,0),30)
  rowsel2[181:182] = c(182,364)
  for(i in 1:182){
    allmarkers2[i+rowsel2[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel2[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel2[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate third allele
  rowsel3 = rep(c(364,546,728,910,0,182),30)
  rowsel3[181:182] = c(364,546)
  for(i in 1:182){
    allmarkers2[i+rowsel3[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel3[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel3[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate fourth allele
  rowsel4 = rep(c(546,728,910,0,182,364),30)
  rowsel4[181:182] = c(546,728)
  for(i in 1:182){
    allmarkers2[i+rowsel4[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel4[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel4[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate fifth allele
  rowsel5 = rep(c(728,910,0,182,364,546),30)
  rowsel5[181:182] = c(728,910)
  for(i in 1:182){
    allmarkers2[i+rowsel5[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel5[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel5[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate sixth allele
  rowsel6 = rep(c(910,0,182,364,546,728),30)
  rowsel6[181:182] = c(910,0)
  for(i in 1:182){
    allmarkers2[i+rowsel6[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel6[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel6[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ## randomly eliminate seventh allele
  rowsel7 = rep(c(0,182,364,546,728,910),30)
  rowsel7[181:182] = c(0,182)
  for(i in 1:182){
    allmarkers2[i+rowsel7[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel7[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel7[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  write.table(allmarkers2,"/Users/Lakshmi/Desktop/Analysis/simdata.txt",sep = "\t",quote = FALSE, row.names=FALSE)
  simdata = read.GeneMapper("/Users/Lakshmi/Desktop/Analysis/simdata.txt")
  Usatnts(simdata) <- c(3,4,4,2,4,2)
  ploidycol = rep(6,length(Samples(simdata)))
  ploidymat = as.matrix(cbind(ploidycol,ploidycol,ploidycol,ploidycol,ploidycol,ploidycol))
  colnames(ploidymat) = c("SEQ18D73","RW39","RW28","SEQ8E8","RW56","RWDI11") 
  Ploidies(simdata) = ploidymat
  simdist = meandistance.matrix(simdata,distmetric = Bruvo.distance)
  simclones = assignClones(simdist,threshold=0.2)
  matches = length(simclones)-length(unique(simclones))
  matchvec[n] = matches
  n = n+1
}

# Testing: 8 null alleles per sample -----
matchvec = rep(190,99)
n=1
i=1

for(i in 1:99){
  ##simulating seq18d73
  seq18d73samplevec = sample(seq18d73namevec2,1092,replace=TRUE,prob=seq18d73freq)
  dfseq18d73 = matrix(seq18d73samplevec,nrow=182,ncol=6)
  dfseq18d73b = (apply(dfseq18d73,1,unique)) #eliminating duplicate alleles 
  dfseq18d73c = ldply(dfseq18d73b, rbind) #converting from list to data frame
  matseq18d73 = as.matrix(dfseq18d73c)
  dfseq18d73d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfseq18d73)[1]){
    testrow = matseq18d73[i,]
    dfseq18d73d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("SEQ18D73",182)
  dfseq18d73e = as.data.frame(dfseq18d73d)
  colnames(dfseq18d73e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfseq18d73f = cbind(Sample.Name,Marker,dfseq18d73e)
  ##simulating RW39
  rw39samplevec = sample(rw39namevec2,1092,replace=TRUE,prob=rw39freq)
  dfrw39 = matrix(rw39samplevec,nrow=182,ncol=6)
  dfrw39b = (apply(dfrw39,1,unique)) #eliminating duplicate alleles 
  dfrw39c = ldply(dfrw39b, rbind) #converting from list to data frame
  matrw39 = as.matrix(dfrw39c)
  dfrw39d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrw39)[1]){
    testrow = matrw39[i,]
    dfrw39d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RW39",182)
  dfrw39e = as.data.frame(dfrw39d)
  colnames(dfrw39e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrw39f = cbind(Sample.Name,Marker,dfrw39e)
  ##simulating RW28
  rw28samplevec = sample(rw28namevec2,1092,replace=TRUE,prob=rw28freq)
  dfrw28 = matrix(rw28samplevec,nrow=182,ncol=6)
  dfrw28b = (apply(dfrw28,1,unique)) #eliminating duplicate alleles 
  dfrw28c = ldply(dfrw28b, rbind) #converting from list to data frame
  matrw28 = as.matrix(dfrw28c)
  dfrw28d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrw28)[1]){
    testrow = matrw28[i,]
    dfrw28d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RW28",182)
  dfrw28e = as.data.frame(dfrw28d)
  colnames(dfrw28e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrw28f = cbind(Sample.Name,Marker,dfrw28e)
  ##simulating SEQ8E8
  seq8e8samplevec = sample(seq8e8namevec2,1092,replace=TRUE,prob=seq8e8freq)
  dfseq8e8 = matrix(seq8e8samplevec,nrow=182,ncol=6)
  dfseq8e8b = (apply(dfseq8e8,1,unique)) #eliminating duplicate alleles 
  dfseq8e8c = ldply(dfseq8e8b, rbind) #converting from list to data frame
  matseq8e8 = as.matrix(dfseq8e8c)
  dfseq8e8d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfseq8e8)[1]){
    testrow = matseq8e8[i,]
    dfseq8e8d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("SEQ8E8",182)
  dfseq8e8e = as.data.frame(dfseq8e8d)
  colnames(dfseq8e8e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfseq8e8f = cbind(Sample.Name,Marker,dfseq8e8e)
  ##simulating RW56
  rw56samplevec = sample(rw56namevec2,1092,replace=TRUE,prob=rw56freq)
  dfrw56 = matrix(rw56samplevec,nrow=182,ncol=6)
  dfrw56b = (apply(dfrw56,1,unique)) #eliminating duplicate alleles 
  dfrw56c = ldply(dfrw56b, rbind) #converting from list to data frame
  matrw56 = as.matrix(dfrw56c)
  dfrw56d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrw56)[1]){
    testrow = matrw56[i,]
    dfrw56d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RW56",182)
  dfrw56e = as.data.frame(dfrw56d)
  colnames(dfrw56e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrw56f = cbind(Sample.Name,Marker,dfrw56e)
  ##simulating RWDI11
  rwdi11samplevec = sample(rwdi11namevec2,1092,replace=TRUE,prob=rwdi11freq)
  dfrwdi11 = matrix(rwdi11samplevec,nrow=182,ncol=6)
  dfrwdi11b = (apply(dfrwdi11,1,unique)) #eliminating duplicate alleles 
  dfrwdi11c = ldply(dfrwdi11b, rbind) #converting from list to data frame
  matrwdi11 = as.matrix(dfrwdi11c)
  dfrwdi11d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrwdi11)[1]){
    testrow = matrwdi11[i,]
    dfrwdi11d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RWDI11",182)
  dfrwdi11e = as.data.frame(dfrwdi11d)
  colnames(dfrwdi11e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrwdi11f = cbind(Sample.Name,Marker,dfrwdi11e)
  ##combining markers, exporting and importing files
  allmarkers=rbind(dfseq18d73f,dfrw39f,dfrw28f,dfseq8e8f,dfrw56f,dfrwdi11f)
  #deleting one allele per individual
  levelsa = c(levels(allmarkers$Allele.1),levels(allmarkers$Allele.2),levels(allmarkers$Allele.3),levels(allmarkers$Allele.4),levels(allmarkers$Allele.5),levels(allmarkers$Allele.6))
  levelsb=unique(levelsa)
  allmarkers$Allele.1 = factor(allmarkers$Allele.1,levels=levelsb)
  allmarkers$Allele.2 = factor(allmarkers$Allele.2,levels=levelsb)
  allmarkers$Allele.3 = factor(allmarkers$Allele.3,levels=levelsb)
  allmarkers$Allele.4 = factor(allmarkers$Allele.4,levels=levelsb)
  allmarkers$Allele.5 = factor(allmarkers$Allele.5,levels=levelsb)
  allmarkers$Allele.6 = factor(allmarkers$Allele.6,levels=levelsb)
  allmarkers2 = allmarkers
  ## randomly eliminate first allele
  rowsel = rep(c(0,182,364,546,728,910),30)
  rowsel[181:182] = c(0,182)
  for(i in 1:182){
    allmarkers2[i+rowsel[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate second allele
  rowsel2 = rep(c(182,364,546,728,910,0),30)
  rowsel2[181:182] = c(182,364)
  for(i in 1:182){
    allmarkers2[i+rowsel2[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel2[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel2[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate third allele
  rowsel3 = rep(c(364,546,728,910,0,182),30)
  rowsel3[181:182] = c(364,546)
  for(i in 1:182){
    allmarkers2[i+rowsel3[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel3[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel3[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate fourth allele
  rowsel4 = rep(c(546,728,910,0,182,364),30)
  rowsel4[181:182] = c(546,728)
  for(i in 1:182){
    allmarkers2[i+rowsel4[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel4[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel4[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate fifth allele
  rowsel5 = rep(c(728,910,0,182,364,546),30)
  rowsel5[181:182] = c(728,910)
  for(i in 1:182){
    allmarkers2[i+rowsel5[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel5[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel5[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate sixth allele
  rowsel6 = rep(c(910,0,182,364,546,728),30)
  rowsel6[181:182] = c(910,0)
  for(i in 1:182){
    allmarkers2[i+rowsel6[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel6[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel6[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ## randomly eliminate seventh allele
  rowsel7 = rep(c(0,182,364,546,728,910),30)
  rowsel7[181:182] = c(0,182)
  for(i in 1:182){
    allmarkers2[i+rowsel7[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel7[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel7[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate 8th allele
  rowsel8 = rep(c(182,364,546,728,910,0),30)
  rowsel8[181:182] = c(182,364)
  for(i in 1:182){
    allmarkers2[i+rowsel8[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel8[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel8[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  write.table(allmarkers2,"/Users/Lakshmi/Desktop/Analysis/simdata.txt",sep = "\t",quote = FALSE, row.names=FALSE)
  simdata = read.GeneMapper("/Users/Lakshmi/Desktop/Analysis/simdata.txt")
  Usatnts(simdata) <- c(3,4,4,2,4,2)
  ploidycol = rep(6,length(Samples(simdata)))
  ploidymat = as.matrix(cbind(ploidycol,ploidycol,ploidycol,ploidycol,ploidycol,ploidycol))
  colnames(ploidymat) = c("SEQ18D73","RW39","RW28","SEQ8E8","RW56","RWDI11") 
  Ploidies(simdata) = ploidymat
  simdist = meandistance.matrix(simdata,distmetric = Bruvo.distance)
  simclones = assignClones(simdist,threshold=0.2)
  matches = length(simclones)-length(unique(simclones))
  matchvec[n] = matches
  n = n+1
}

# Testing: 9 null alleles per sample -----
matchvec = rep(190,99)
n=1
i=1

for(i in 1:99){
  ##simulating seq18d73
  seq18d73samplevec = sample(seq18d73namevec2,1092,replace=TRUE,prob=seq18d73freq)
  dfseq18d73 = matrix(seq18d73samplevec,nrow=182,ncol=6)
  dfseq18d73b = (apply(dfseq18d73,1,unique)) #eliminating duplicate alleles 
  dfseq18d73c = ldply(dfseq18d73b, rbind) #converting from list to data frame
  matseq18d73 = as.matrix(dfseq18d73c)
  dfseq18d73d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfseq18d73)[1]){
    testrow = matseq18d73[i,]
    dfseq18d73d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("SEQ18D73",182)
  dfseq18d73e = as.data.frame(dfseq18d73d)
  colnames(dfseq18d73e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfseq18d73f = cbind(Sample.Name,Marker,dfseq18d73e)
  ##simulating RW39
  rw39samplevec = sample(rw39namevec2,1092,replace=TRUE,prob=rw39freq)
  dfrw39 = matrix(rw39samplevec,nrow=182,ncol=6)
  dfrw39b = (apply(dfrw39,1,unique)) #eliminating duplicate alleles 
  dfrw39c = ldply(dfrw39b, rbind) #converting from list to data frame
  matrw39 = as.matrix(dfrw39c)
  dfrw39d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrw39)[1]){
    testrow = matrw39[i,]
    dfrw39d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RW39",182)
  dfrw39e = as.data.frame(dfrw39d)
  colnames(dfrw39e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrw39f = cbind(Sample.Name,Marker,dfrw39e)
  ##simulating RW28
  rw28samplevec = sample(rw28namevec2,1092,replace=TRUE,prob=rw28freq)
  dfrw28 = matrix(rw28samplevec,nrow=182,ncol=6)
  dfrw28b = (apply(dfrw28,1,unique)) #eliminating duplicate alleles 
  dfrw28c = ldply(dfrw28b, rbind) #converting from list to data frame
  matrw28 = as.matrix(dfrw28c)
  dfrw28d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrw28)[1]){
    testrow = matrw28[i,]
    dfrw28d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RW28",182)
  dfrw28e = as.data.frame(dfrw28d)
  colnames(dfrw28e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrw28f = cbind(Sample.Name,Marker,dfrw28e)
  ##simulating SEQ8E8
  seq8e8samplevec = sample(seq8e8namevec2,1092,replace=TRUE,prob=seq8e8freq)
  dfseq8e8 = matrix(seq8e8samplevec,nrow=182,ncol=6)
  dfseq8e8b = (apply(dfseq8e8,1,unique)) #eliminating duplicate alleles 
  dfseq8e8c = ldply(dfseq8e8b, rbind) #converting from list to data frame
  matseq8e8 = as.matrix(dfseq8e8c)
  dfseq8e8d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfseq8e8)[1]){
    testrow = matseq8e8[i,]
    dfseq8e8d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("SEQ8E8",182)
  dfseq8e8e = as.data.frame(dfseq8e8d)
  colnames(dfseq8e8e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfseq8e8f = cbind(Sample.Name,Marker,dfseq8e8e)
  ##simulating RW56
  rw56samplevec = sample(rw56namevec2,1092,replace=TRUE,prob=rw56freq)
  dfrw56 = matrix(rw56samplevec,nrow=182,ncol=6)
  dfrw56b = (apply(dfrw56,1,unique)) #eliminating duplicate alleles 
  dfrw56c = ldply(dfrw56b, rbind) #converting from list to data frame
  matrw56 = as.matrix(dfrw56c)
  dfrw56d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrw56)[1]){
    testrow = matrw56[i,]
    dfrw56d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RW56",182)
  dfrw56e = as.data.frame(dfrw56d)
  colnames(dfrw56e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrw56f = cbind(Sample.Name,Marker,dfrw56e)
  ##simulating RWDI11
  rwdi11samplevec = sample(rwdi11namevec2,1092,replace=TRUE,prob=rwdi11freq)
  dfrwdi11 = matrix(rwdi11samplevec,nrow=182,ncol=6)
  dfrwdi11b = (apply(dfrwdi11,1,unique)) #eliminating duplicate alleles 
  dfrwdi11c = ldply(dfrwdi11b, rbind) #converting from list to data frame
  matrwdi11 = as.matrix(dfrwdi11c)
  dfrwdi11d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrwdi11)[1]){
    testrow = matrwdi11[i,]
    dfrwdi11d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RWDI11",182)
  dfrwdi11e = as.data.frame(dfrwdi11d)
  colnames(dfrwdi11e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrwdi11f = cbind(Sample.Name,Marker,dfrwdi11e)
  ##combining markers, exporting and importing files
  allmarkers=rbind(dfseq18d73f,dfrw39f,dfrw28f,dfseq8e8f,dfrw56f,dfrwdi11f)
  #deleting one allele per individual
  levelsa = c(levels(allmarkers$Allele.1),levels(allmarkers$Allele.2),levels(allmarkers$Allele.3),levels(allmarkers$Allele.4),levels(allmarkers$Allele.5),levels(allmarkers$Allele.6))
  levelsb=unique(levelsa)
  allmarkers$Allele.1 = factor(allmarkers$Allele.1,levels=levelsb)
  allmarkers$Allele.2 = factor(allmarkers$Allele.2,levels=levelsb)
  allmarkers$Allele.3 = factor(allmarkers$Allele.3,levels=levelsb)
  allmarkers$Allele.4 = factor(allmarkers$Allele.4,levels=levelsb)
  allmarkers$Allele.5 = factor(allmarkers$Allele.5,levels=levelsb)
  allmarkers$Allele.6 = factor(allmarkers$Allele.6,levels=levelsb)
  allmarkers2 = allmarkers
  ## randomly eliminate first allele
  rowsel = rep(c(0,182,364,546,728,910),30)
  rowsel[181:182] = c(0,182)
  for(i in 1:182){
    allmarkers2[i+rowsel[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate second allele
  rowsel2 = rep(c(182,364,546,728,910,0),30)
  rowsel2[181:182] = c(182,364)
  for(i in 1:182){
    allmarkers2[i+rowsel2[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel2[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel2[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate third allele
  rowsel3 = rep(c(364,546,728,910,0,182),30)
  rowsel3[181:182] = c(364,546)
  for(i in 1:182){
    allmarkers2[i+rowsel3[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel3[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel3[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate fourth allele
  rowsel4 = rep(c(546,728,910,0,182,364),30)
  rowsel4[181:182] = c(546,728)
  for(i in 1:182){
    allmarkers2[i+rowsel4[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel4[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel4[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate fifth allele
  rowsel5 = rep(c(728,910,0,182,364,546),30)
  rowsel5[181:182] = c(728,910)
  for(i in 1:182){
    allmarkers2[i+rowsel5[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel5[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel5[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate sixth allele
  rowsel6 = rep(c(910,0,182,364,546,728),30)
  rowsel6[181:182] = c(910,0)
  for(i in 1:182){
    allmarkers2[i+rowsel6[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel6[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel6[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ## randomly eliminate seventh allele
  rowsel7 = rep(c(0,182,364,546,728,910),30)
  rowsel7[181:182] = c(0,182)
  for(i in 1:182){
    allmarkers2[i+rowsel7[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel7[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel7[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate 8th allele
  rowsel8 = rep(c(182,364,546,728,910,0),30)
  rowsel8[181:182] = c(182,364)
  for(i in 1:182){
    allmarkers2[i+rowsel8[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel8[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel8[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate 9th allele (like 3rd)
  rowsel9 = rep(c(364,546,728,910,0,182),30)
  rowsel9[181:182] = c(364,546)
  for(i in 1:182){
    allmarkers2[i+rowsel9[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel9[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel9[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  write.table(allmarkers2,"/Users/Lakshmi/Desktop/Analysis/simdata.txt",sep = "\t",quote = FALSE, row.names=FALSE)
  simdata = read.GeneMapper("/Users/Lakshmi/Desktop/Analysis/simdata.txt")
  Usatnts(simdata) <- c(3,4,4,2,4,2)
  ploidycol = rep(6,length(Samples(simdata)))
  ploidymat = as.matrix(cbind(ploidycol,ploidycol,ploidycol,ploidycol,ploidycol,ploidycol))
  colnames(ploidymat) = c("SEQ18D73","RW39","RW28","SEQ8E8","RW56","RWDI11") 
  Ploidies(simdata) = ploidymat
  simdist = meandistance.matrix(simdata,distmetric = Bruvo.distance)
  simclones = assignClones(simdist,threshold=0.2)
  matches = length(simclones)-length(unique(simclones))
  matchvec[n] = matches
  n = n+1
}

# Testing: 10 null alleles per sample -----
matchvec = rep(190,99)
n=1
i=1

for(i in 1:99){
  ##simulating seq18d73
  seq18d73samplevec = sample(seq18d73namevec2,1092,replace=TRUE,prob=seq18d73freq)
  dfseq18d73 = matrix(seq18d73samplevec,nrow=182,ncol=6)
  dfseq18d73b = (apply(dfseq18d73,1,unique)) #eliminating duplicate alleles 
  dfseq18d73c = ldply(dfseq18d73b, rbind) #converting from list to data frame
  matseq18d73 = as.matrix(dfseq18d73c)
  dfseq18d73d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfseq18d73)[1]){
    testrow = matseq18d73[i,]
    dfseq18d73d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("SEQ18D73",182)
  dfseq18d73e = as.data.frame(dfseq18d73d)
  colnames(dfseq18d73e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfseq18d73f = cbind(Sample.Name,Marker,dfseq18d73e)
  ##simulating RW39
  rw39samplevec = sample(rw39namevec2,1092,replace=TRUE,prob=rw39freq)
  dfrw39 = matrix(rw39samplevec,nrow=182,ncol=6)
  dfrw39b = (apply(dfrw39,1,unique)) #eliminating duplicate alleles 
  dfrw39c = ldply(dfrw39b, rbind) #converting from list to data frame
  matrw39 = as.matrix(dfrw39c)
  dfrw39d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrw39)[1]){
    testrow = matrw39[i,]
    dfrw39d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RW39",182)
  dfrw39e = as.data.frame(dfrw39d)
  colnames(dfrw39e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrw39f = cbind(Sample.Name,Marker,dfrw39e)
  ##simulating RW28
  rw28samplevec = sample(rw28namevec2,1092,replace=TRUE,prob=rw28freq)
  dfrw28 = matrix(rw28samplevec,nrow=182,ncol=6)
  dfrw28b = (apply(dfrw28,1,unique)) #eliminating duplicate alleles 
  dfrw28c = ldply(dfrw28b, rbind) #converting from list to data frame
  matrw28 = as.matrix(dfrw28c)
  dfrw28d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrw28)[1]){
    testrow = matrw28[i,]
    dfrw28d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RW28",182)
  dfrw28e = as.data.frame(dfrw28d)
  colnames(dfrw28e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrw28f = cbind(Sample.Name,Marker,dfrw28e)
  ##simulating SEQ8E8
  seq8e8samplevec = sample(seq8e8namevec2,1092,replace=TRUE,prob=seq8e8freq)
  dfseq8e8 = matrix(seq8e8samplevec,nrow=182,ncol=6)
  dfseq8e8b = (apply(dfseq8e8,1,unique)) #eliminating duplicate alleles 
  dfseq8e8c = ldply(dfseq8e8b, rbind) #converting from list to data frame
  matseq8e8 = as.matrix(dfseq8e8c)
  dfseq8e8d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfseq8e8)[1]){
    testrow = matseq8e8[i,]
    dfseq8e8d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("SEQ8E8",182)
  dfseq8e8e = as.data.frame(dfseq8e8d)
  colnames(dfseq8e8e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfseq8e8f = cbind(Sample.Name,Marker,dfseq8e8e)
  ##simulating RW56
  rw56samplevec = sample(rw56namevec2,1092,replace=TRUE,prob=rw56freq)
  dfrw56 = matrix(rw56samplevec,nrow=182,ncol=6)
  dfrw56b = (apply(dfrw56,1,unique)) #eliminating duplicate alleles 
  dfrw56c = ldply(dfrw56b, rbind) #converting from list to data frame
  matrw56 = as.matrix(dfrw56c)
  dfrw56d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrw56)[1]){
    testrow = matrw56[i,]
    dfrw56d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RW56",182)
  dfrw56e = as.data.frame(dfrw56d)
  colnames(dfrw56e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrw56f = cbind(Sample.Name,Marker,dfrw56e)
  ##simulating RWDI11
  rwdi11samplevec = sample(rwdi11namevec2,1092,replace=TRUE,prob=rwdi11freq)
  dfrwdi11 = matrix(rwdi11samplevec,nrow=182,ncol=6)
  dfrwdi11b = (apply(dfrwdi11,1,unique)) #eliminating duplicate alleles 
  dfrwdi11c = ldply(dfrwdi11b, rbind) #converting from list to data frame
  matrwdi11 = as.matrix(dfrwdi11c)
  dfrwdi11d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrwdi11)[1]){
    testrow = matrwdi11[i,]
    dfrwdi11d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RWDI11",182)
  dfrwdi11e = as.data.frame(dfrwdi11d)
  colnames(dfrwdi11e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrwdi11f = cbind(Sample.Name,Marker,dfrwdi11e)
  ##combining markers, exporting and importing files
  allmarkers=rbind(dfseq18d73f,dfrw39f,dfrw28f,dfseq8e8f,dfrw56f,dfrwdi11f)
  #deleting one allele per individual
  levelsa = c(levels(allmarkers$Allele.1),levels(allmarkers$Allele.2),levels(allmarkers$Allele.3),levels(allmarkers$Allele.4),levels(allmarkers$Allele.5),levels(allmarkers$Allele.6))
  levelsb=unique(levelsa)
  allmarkers$Allele.1 = factor(allmarkers$Allele.1,levels=levelsb)
  allmarkers$Allele.2 = factor(allmarkers$Allele.2,levels=levelsb)
  allmarkers$Allele.3 = factor(allmarkers$Allele.3,levels=levelsb)
  allmarkers$Allele.4 = factor(allmarkers$Allele.4,levels=levelsb)
  allmarkers$Allele.5 = factor(allmarkers$Allele.5,levels=levelsb)
  allmarkers$Allele.6 = factor(allmarkers$Allele.6,levels=levelsb)
  allmarkers2 = allmarkers
  ## randomly eliminate first allele
  rowsel = rep(c(0,182,364,546,728,910),30)
  rowsel[181:182] = c(0,182)
  for(i in 1:182){
    allmarkers2[i+rowsel[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate second allele
  rowsel2 = rep(c(182,364,546,728,910,0),30)
  rowsel2[181:182] = c(182,364)
  for(i in 1:182){
    allmarkers2[i+rowsel2[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel2[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel2[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate third allele
  rowsel3 = rep(c(364,546,728,910,0,182),30)
  rowsel3[181:182] = c(364,546)
  for(i in 1:182){
    allmarkers2[i+rowsel3[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel3[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel3[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate fourth allele
  rowsel4 = rep(c(546,728,910,0,182,364),30)
  rowsel4[181:182] = c(546,728)
  for(i in 1:182){
    allmarkers2[i+rowsel4[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel4[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel4[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate fifth allele
  rowsel5 = rep(c(728,910,0,182,364,546),30)
  rowsel5[181:182] = c(728,910)
  for(i in 1:182){
    allmarkers2[i+rowsel5[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel5[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel5[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate sixth allele
  rowsel6 = rep(c(910,0,182,364,546,728),30)
  rowsel6[181:182] = c(910,0)
  for(i in 1:182){
    allmarkers2[i+rowsel6[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel6[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel6[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ## randomly eliminate seventh allele
  rowsel7 = rep(c(0,182,364,546,728,910),30)
  rowsel7[181:182] = c(0,182)
  for(i in 1:182){
    allmarkers2[i+rowsel7[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel7[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel7[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate 8th allele
  rowsel8 = rep(c(182,364,546,728,910,0),30)
  rowsel8[181:182] = c(182,364)
  for(i in 1:182){
    allmarkers2[i+rowsel8[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel8[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel8[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate 9th allele (like 3rd)
  rowsel9 = rep(c(364,546,728,910,0,182),30)
  rowsel9[181:182] = c(364,546)
  for(i in 1:182){
    allmarkers2[i+rowsel9[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel9[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel9[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate 10th allele (like 4th)
  rowsel10 = rep(c(546,728,910,0,182,364),30)
  rowsel10[181:182] = c(546,728)
  for(i in 1:182){
    allmarkers2[i+rowsel10[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel10[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel10[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  write.table(allmarkers2,"/Users/Lakshmi/Desktop/Analysis/simdata.txt",sep = "\t",quote = FALSE, row.names=FALSE)
  simdata = read.GeneMapper("/Users/Lakshmi/Desktop/Analysis/simdata.txt")
  Usatnts(simdata) <- c(3,4,4,2,4,2)
  ploidycol = rep(6,length(Samples(simdata)))
  ploidymat = as.matrix(cbind(ploidycol,ploidycol,ploidycol,ploidycol,ploidycol,ploidycol))
  colnames(ploidymat) = c("SEQ18D73","RW39","RW28","SEQ8E8","RW56","RWDI11") 
  Ploidies(simdata) = ploidymat
  simdist = meandistance.matrix(simdata,distmetric = Bruvo.distance)
  simclones = assignClones(simdist,threshold=0.2)
  matches = length(simclones)-length(unique(simclones))
  matchvec[n] = matches
  n = n+1
}

# Testing: 14 null alleles per sample -----
matchvec14 = rep(190,100)
n=1
i=1

for(i in 1:100){
  ##simulating seq18d73
  seq18d73samplevec = sample(seq18d73namevec2,1092,replace=TRUE,prob=seq18d73freq)
  dfseq18d73 = matrix(seq18d73samplevec,nrow=182,ncol=6)
  dfseq18d73b = (apply(dfseq18d73,1,unique)) #eliminating duplicate alleles 
  dfseq18d73c = ldply(dfseq18d73b, rbind) #converting from list to data frame
  matseq18d73 = as.matrix(dfseq18d73c)
  dfseq18d73d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfseq18d73)[1]){
    testrow = matseq18d73[i,]
    dfseq18d73d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("SEQ18D73",182)
  dfseq18d73e = as.data.frame(dfseq18d73d)
  colnames(dfseq18d73e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfseq18d73f = cbind(Sample.Name,Marker,dfseq18d73e)
  ##simulating RW39
  rw39samplevec = sample(rw39namevec2,1092,replace=TRUE,prob=rw39freq)
  dfrw39 = matrix(rw39samplevec,nrow=182,ncol=6)
  dfrw39b = (apply(dfrw39,1,unique)) #eliminating duplicate alleles 
  dfrw39c = ldply(dfrw39b, rbind) #converting from list to data frame
  matrw39 = as.matrix(dfrw39c)
  dfrw39d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrw39)[1]){
    testrow = matrw39[i,]
    dfrw39d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RW39",182)
  dfrw39e = as.data.frame(dfrw39d)
  colnames(dfrw39e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrw39f = cbind(Sample.Name,Marker,dfrw39e)
  ##simulating RW28
  rw28samplevec = sample(rw28namevec2,1092,replace=TRUE,prob=rw28freq)
  dfrw28 = matrix(rw28samplevec,nrow=182,ncol=6)
  dfrw28b = (apply(dfrw28,1,unique)) #eliminating duplicate alleles 
  dfrw28c = ldply(dfrw28b, rbind) #converting from list to data frame
  matrw28 = as.matrix(dfrw28c)
  dfrw28d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrw28)[1]){
    testrow = matrw28[i,]
    dfrw28d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RW28",182)
  dfrw28e = as.data.frame(dfrw28d)
  colnames(dfrw28e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrw28f = cbind(Sample.Name,Marker,dfrw28e)
  ##simulating SEQ8E8
  seq8e8samplevec = sample(seq8e8namevec2,1092,replace=TRUE,prob=seq8e8freq)
  dfseq8e8 = matrix(seq8e8samplevec,nrow=182,ncol=6)
  dfseq8e8b = (apply(dfseq8e8,1,unique)) #eliminating duplicate alleles 
  dfseq8e8c = ldply(dfseq8e8b, rbind) #converting from list to data frame
  matseq8e8 = as.matrix(dfseq8e8c)
  dfseq8e8d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfseq8e8)[1]){
    testrow = matseq8e8[i,]
    dfseq8e8d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("SEQ8E8",182)
  dfseq8e8e = as.data.frame(dfseq8e8d)
  colnames(dfseq8e8e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfseq8e8f = cbind(Sample.Name,Marker,dfseq8e8e)
  ##simulating RW56
  rw56samplevec = sample(rw56namevec2,1092,replace=TRUE,prob=rw56freq)
  dfrw56 = matrix(rw56samplevec,nrow=182,ncol=6)
  dfrw56b = (apply(dfrw56,1,unique)) #eliminating duplicate alleles 
  dfrw56c = ldply(dfrw56b, rbind) #converting from list to data frame
  matrw56 = as.matrix(dfrw56c)
  dfrw56d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrw56)[1]){
    testrow = matrw56[i,]
    dfrw56d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RW56",182)
  dfrw56e = as.data.frame(dfrw56d)
  colnames(dfrw56e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrw56f = cbind(Sample.Name,Marker,dfrw56e)
  ##simulating RWDI11
  rwdi11samplevec = sample(rwdi11namevec2,1092,replace=TRUE,prob=rwdi11freq)
  dfrwdi11 = matrix(rwdi11samplevec,nrow=182,ncol=6)
  dfrwdi11b = (apply(dfrwdi11,1,unique)) #eliminating duplicate alleles 
  dfrwdi11c = ldply(dfrwdi11b, rbind) #converting from list to data frame
  matrwdi11 = as.matrix(dfrwdi11c)
  dfrwdi11d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrwdi11)[1]){
    testrow = matrwdi11[i,]
    dfrwdi11d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RWDI11",182)
  dfrwdi11e = as.data.frame(dfrwdi11d)
  colnames(dfrwdi11e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrwdi11f = cbind(Sample.Name,Marker,dfrwdi11e)
  ##combining markers, exporting and importing files
  allmarkers=rbind(dfseq18d73f,dfrw39f,dfrw28f,dfseq8e8f,dfrw56f,dfrwdi11f)
  #deleting one allele per individual
  levelsa = c(levels(allmarkers$Allele.1),levels(allmarkers$Allele.2),levels(allmarkers$Allele.3),levels(allmarkers$Allele.4),levels(allmarkers$Allele.5),levels(allmarkers$Allele.6))
  levelsb=unique(levelsa)
  allmarkers$Allele.1 = factor(allmarkers$Allele.1,levels=levelsb)
  allmarkers$Allele.2 = factor(allmarkers$Allele.2,levels=levelsb)
  allmarkers$Allele.3 = factor(allmarkers$Allele.3,levels=levelsb)
  allmarkers$Allele.4 = factor(allmarkers$Allele.4,levels=levelsb)
  allmarkers$Allele.5 = factor(allmarkers$Allele.5,levels=levelsb)
  allmarkers$Allele.6 = factor(allmarkers$Allele.6,levels=levelsb)
  allmarkers2 = allmarkers
  ## randomly eliminate first allele
  rowsel = rep(c(0,182,364,546,728,910),30)
  rowsel[181:182] = c(0,182)
  for(i in 1:182){
    allmarkers2[i+rowsel[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate second allele
  rowsel2 = rep(c(182,364,546,728,910,0),30)
  rowsel2[181:182] = c(182,364)
  for(i in 1:182){
    allmarkers2[i+rowsel2[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel2[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel2[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate third allele
  rowsel3 = rep(c(364,546,728,910,0,182),30)
  rowsel3[181:182] = c(364,546)
  for(i in 1:182){
    allmarkers2[i+rowsel3[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel3[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel3[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate fourth allele
  rowsel4 = rep(c(546,728,910,0,182,364),30)
  rowsel4[181:182] = c(546,728)
  for(i in 1:182){
    allmarkers2[i+rowsel4[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel4[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel4[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate fifth allele
  rowsel5 = rep(c(728,910,0,182,364,546),30)
  rowsel5[181:182] = c(728,910)
  for(i in 1:182){
    allmarkers2[i+rowsel5[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel5[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel5[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate sixth allele
  rowsel6 = rep(c(910,0,182,364,546,728),30)
  rowsel6[181:182] = c(910,0)
  for(i in 1:182){
    allmarkers2[i+rowsel6[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel6[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel6[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ## randomly eliminate seventh allele
  rowsel7 = rep(c(0,182,364,546,728,910),30)
  rowsel7[181:182] = c(0,182)
  for(i in 1:182){
    allmarkers2[i+rowsel7[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel7[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel7[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate 8th allele
  rowsel8 = rep(c(182,364,546,728,910,0),30)
  rowsel8[181:182] = c(182,364)
  for(i in 1:182){
    allmarkers2[i+rowsel8[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel8[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel8[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate 9th allele (like 3rd)
  rowsel9 = rep(c(364,546,728,910,0,182),30)
  rowsel9[181:182] = c(364,546)
  for(i in 1:182){
    allmarkers2[i+rowsel9[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel9[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel9[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate 10th allele (like 4th)
  rowsel10 = rep(c(546,728,910,0,182,364),30)
  rowsel10[181:182] = c(546,728)
  for(i in 1:182){
    allmarkers2[i+rowsel10[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel10[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel10[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate 11th allele (like 5th)
  for(i in 1:182){
    allmarkers2[i+rowsel5[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel5[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel5[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate 12th allele (like 6th)
  for(i in 1:182){
    allmarkers2[i+rowsel6[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel6[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel6[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate 13th allele (like 1st)
  for(i in 1:182){
    allmarkers2[i+rowsel[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate 14th allele (like 2nd)
  for(i in 1:182){
    allmarkers2[i+rowsel2[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel2[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel2[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  write.table(allmarkers2,"/Users/Lakshmi/Desktop/Analysis/simdata.txt",sep = "\t",quote = FALSE, row.names=FALSE)
  simdata = read.GeneMapper("/Users/Lakshmi/Desktop/Analysis/simdata.txt")
  Usatnts(simdata) <- c(3,4,4,2,4,2)
  ploidycol = rep(6,length(Samples(simdata)))
  ploidymat = as.matrix(cbind(ploidycol,ploidycol,ploidycol,ploidycol,ploidycol,ploidycol))
  colnames(ploidymat) = c("SEQ18D73","RW39","RW28","SEQ8E8","RW56","RWDI11") 
  Ploidies(simdata) = ploidymat
  simdist = meandistance.matrix(simdata,distmetric = Bruvo.distance)
  simclones = assignClones(simdist,threshold=0.2)
  matches = length(simclones)-length(unique(simclones))
  matchvec14[n] = matches
  n = n+1
}


# Testing: 22 null alleles per sample with writing data files-----

namevec1 = rep("/Users/Lakshmi/Desktop/Analysis/simdata",100)
namevec2 = seq(1:100)
namevec3 = rep(".txt",100)
namevec4 = paste(namevec1,namevec2,namevec3,sep="")

matchvec22 = rep(190,99)
n=1
i=1

for(i in 1:99){
  ##simulating seq18d73
  seq18d73samplevec = sample(seq18d73namevec2,1092,replace=TRUE,prob=seq18d73freq)
  dfseq18d73 = matrix(seq18d73samplevec,nrow=182,ncol=6)
  dfseq18d73b = (apply(dfseq18d73,1,unique)) #eliminating duplicate alleles 
  dfseq18d73c = ldply(dfseq18d73b, rbind) #converting from list to data frame
  matseq18d73 = as.matrix(dfseq18d73c)
  dfseq18d73d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfseq18d73)[1]){
    testrow = matseq18d73[i,]
    dfseq18d73d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("SEQ18D73",182)
  dfseq18d73e = as.data.frame(dfseq18d73d)
  colnames(dfseq18d73e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfseq18d73f = cbind(Sample.Name,Marker,dfseq18d73e)
  ##simulating RW39
  rw39samplevec = sample(rw39namevec2,1092,replace=TRUE,prob=rw39freq)
  dfrw39 = matrix(rw39samplevec,nrow=182,ncol=6)
  dfrw39b = (apply(dfrw39,1,unique)) #eliminating duplicate alleles 
  dfrw39c = ldply(dfrw39b, rbind) #converting from list to data frame
  matrw39 = as.matrix(dfrw39c)
  dfrw39d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrw39)[1]){
    testrow = matrw39[i,]
    dfrw39d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RW39",182)
  dfrw39e = as.data.frame(dfrw39d)
  colnames(dfrw39e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrw39f = cbind(Sample.Name,Marker,dfrw39e)
  ##simulating RW28
  rw28samplevec = sample(rw28namevec2,1092,replace=TRUE,prob=rw28freq)
  dfrw28 = matrix(rw28samplevec,nrow=182,ncol=6)
  dfrw28b = (apply(dfrw28,1,unique)) #eliminating duplicate alleles 
  dfrw28c = ldply(dfrw28b, rbind) #converting from list to data frame
  matrw28 = as.matrix(dfrw28c)
  dfrw28d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrw28)[1]){
    testrow = matrw28[i,]
    dfrw28d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RW28",182)
  dfrw28e = as.data.frame(dfrw28d)
  colnames(dfrw28e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrw28f = cbind(Sample.Name,Marker,dfrw28e)
  ##simulating SEQ8E8
  seq8e8samplevec = sample(seq8e8namevec2,1092,replace=TRUE,prob=seq8e8freq)
  dfseq8e8 = matrix(seq8e8samplevec,nrow=182,ncol=6)
  dfseq8e8b = (apply(dfseq8e8,1,unique)) #eliminating duplicate alleles 
  dfseq8e8c = ldply(dfseq8e8b, rbind) #converting from list to data frame
  matseq8e8 = as.matrix(dfseq8e8c)
  dfseq8e8d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfseq8e8)[1]){
    testrow = matseq8e8[i,]
    dfseq8e8d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("SEQ8E8",182)
  dfseq8e8e = as.data.frame(dfseq8e8d)
  colnames(dfseq8e8e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfseq8e8f = cbind(Sample.Name,Marker,dfseq8e8e)
  ##simulating RW56
  rw56samplevec = sample(rw56namevec2,1092,replace=TRUE,prob=rw56freq)
  dfrw56 = matrix(rw56samplevec,nrow=182,ncol=6)
  dfrw56b = (apply(dfrw56,1,unique)) #eliminating duplicate alleles 
  dfrw56c = ldply(dfrw56b, rbind) #converting from list to data frame
  matrw56 = as.matrix(dfrw56c)
  dfrw56d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrw56)[1]){
    testrow = matrw56[i,]
    dfrw56d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RW56",182)
  dfrw56e = as.data.frame(dfrw56d)
  colnames(dfrw56e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrw56f = cbind(Sample.Name,Marker,dfrw56e)
  ##simulating RWDI11
  rwdi11samplevec = sample(rwdi11namevec2,1092,replace=TRUE,prob=rwdi11freq)
  dfrwdi11 = matrix(rwdi11samplevec,nrow=182,ncol=6)
  dfrwdi11b = (apply(dfrwdi11,1,unique)) #eliminating duplicate alleles 
  dfrwdi11c = ldply(dfrwdi11b, rbind) #converting from list to data frame
  matrwdi11 = as.matrix(dfrwdi11c)
  dfrwdi11d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrwdi11)[1]){
    testrow = matrwdi11[i,]
    dfrwdi11d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RWDI11",182)
  dfrwdi11e = as.data.frame(dfrwdi11d)
  colnames(dfrwdi11e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrwdi11f = cbind(Sample.Name,Marker,dfrwdi11e)
  ##combining markers, exporting and importing files
  allmarkers=rbind(dfseq18d73f,dfrw39f,dfrw28f,dfseq8e8f,dfrw56f,dfrwdi11f)
  #deleting one allele per individual
  levelsa = c(levels(allmarkers$Allele.1),levels(allmarkers$Allele.2),levels(allmarkers$Allele.3),levels(allmarkers$Allele.4),levels(allmarkers$Allele.5),levels(allmarkers$Allele.6))
  levelsb=unique(levelsa)
  allmarkers$Allele.1 = factor(allmarkers$Allele.1,levels=levelsb)
  allmarkers$Allele.2 = factor(allmarkers$Allele.2,levels=levelsb)
  allmarkers$Allele.3 = factor(allmarkers$Allele.3,levels=levelsb)
  allmarkers$Allele.4 = factor(allmarkers$Allele.4,levels=levelsb)
  allmarkers$Allele.5 = factor(allmarkers$Allele.5,levels=levelsb)
  allmarkers$Allele.6 = factor(allmarkers$Allele.6,levels=levelsb)
  allmarkers2 = allmarkers
  ## randomly eliminate first allele
  rowsel = rep(c(0,182,364,546,728,910),30)
  rowsel[181:182] = c(0,182)
  for(i in 1:182){
    allmarkers2[i+rowsel[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate second allele
  rowsel2 = rep(c(182,364,546,728,910,0),30)
  rowsel2[181:182] = c(182,364)
  for(i in 1:182){
    allmarkers2[i+rowsel2[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel2[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel2[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate third allele
  rowsel3 = rep(c(364,546,728,910,0,182),30)
  rowsel3[181:182] = c(364,546)
  for(i in 1:182){
    allmarkers2[i+rowsel3[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel3[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel3[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate fourth allele
  rowsel4 = rep(c(546,728,910,0,182,364),30)
  rowsel4[181:182] = c(546,728)
  for(i in 1:182){
    allmarkers2[i+rowsel4[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel4[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel4[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate fifth allele
  rowsel5 = rep(c(728,910,0,182,364,546),30)
  rowsel5[181:182] = c(728,910)
  for(i in 1:182){
    allmarkers2[i+rowsel5[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel5[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel5[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate sixth allele
  rowsel6 = rep(c(910,0,182,364,546,728),30)
  rowsel6[181:182] = c(910,0)
  for(i in 1:182){
    allmarkers2[i+rowsel6[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel6[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel6[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ## randomly eliminate seventh allele
  rowsel7 = rep(c(0,182,364,546,728,910),30)
  rowsel7[181:182] = c(0,182)
  for(i in 1:182){
    allmarkers2[i+rowsel7[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel7[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel7[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate 8th allele
  rowsel8 = rep(c(182,364,546,728,910,0),30)
  rowsel8[181:182] = c(182,364)
  for(i in 1:182){
    allmarkers2[i+rowsel8[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel8[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel8[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate 9th allele (like 3rd)
  rowsel9 = rep(c(364,546,728,910,0,182),30)
  rowsel9[181:182] = c(364,546)
  for(i in 1:182){
    allmarkers2[i+rowsel9[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel9[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel9[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate 10th allele (like 4th)
  rowsel10 = rep(c(546,728,910,0,182,364),30)
  rowsel10[181:182] = c(546,728)
  for(i in 1:182){
    allmarkers2[i+rowsel10[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel10[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel10[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate 11th allele (like 5th)
  for(i in 1:182){
    allmarkers2[i+rowsel5[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel5[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel5[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate 12th allele (like 6th)
  for(i in 1:182){
    allmarkers2[i+rowsel6[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel6[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel6[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate 13th allele (like 1st)
  for(i in 1:182){
    allmarkers2[i+rowsel[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate 14th allele (like 2nd)
  for(i in 1:182){
    allmarkers2[i+rowsel2[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel2[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel2[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate 15th allele (like 3rd)
  for(i in 1:182){
    allmarkers2[i+rowsel3[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel3[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel3[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate 16th allele (like 4th)
  for(i in 1:182){
    allmarkers2[i+rowsel4[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel4[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel4[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate 17th allele (like 5th)
  for(i in 1:182){
    allmarkers2[i+rowsel5[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel5[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel5[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate 18th allele (like 6th)
  for(i in 1:182){
    allmarkers2[i+rowsel6[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel6[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel6[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate 19th allele (like 1st)
  for(i in 1:182){
    allmarkers2[i+rowsel[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate 20th allele (like 2nd)
  for(i in 1:182){
    allmarkers2[i+rowsel2[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel2[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel2[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate 21st allele (like 3rd)
  for(i in 1:182){
    allmarkers2[i+rowsel3[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel3[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel3[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate 22nd allele (like 4th)
  for(i in 1:182){
    allmarkers2[i+rowsel4[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel4[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel4[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  write.table(allmarkers2,namevec4[n],sep = "\t",quote = FALSE, row.names=FALSE)
  simdata = read.GeneMapper("/Users/Lakshmi/Desktop/Analysis/simdata.txt")
  Usatnts(simdata) <- c(3,4,4,2,4,2)
  ploidycol = rep(6,length(Samples(simdata)))
  ploidymat = as.matrix(cbind(ploidycol,ploidycol,ploidycol,ploidycol,ploidycol,ploidycol))
  colnames(ploidymat) = c("SEQ18D73","RW39","RW28","SEQ8E8","RW56","RWDI11") 
  Ploidies(simdata) = ploidymat
  simdist = meandistance.matrix(simdata,distmetric = Bruvo.distance)
  simclones = assignClones(simdist,threshold=0.2)
  matches = length(simclones)-length(unique(simclones))
  matchvec22[n] = matches
  n = n+1
}

## Testing 23 null alleles -----

matchvec23 = rep(190,100)
n=1
i=1

for(i in 1:100){
  ##simulating seq18d73
  seq18d73samplevec = sample(seq18d73namevec2,1092,replace=TRUE,prob=seq18d73freq)
  dfseq18d73 = matrix(seq18d73samplevec,nrow=182,ncol=6)
  dfseq18d73b = (apply(dfseq18d73,1,unique)) #eliminating duplicate alleles 
  dfseq18d73c = ldply(dfseq18d73b, rbind) #converting from list to data frame
  matseq18d73 = as.matrix(dfseq18d73c)
  dfseq18d73d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfseq18d73)[1]){
    testrow = matseq18d73[i,]
    dfseq18d73d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("SEQ18D73",182)
  dfseq18d73e = as.data.frame(dfseq18d73d)
  colnames(dfseq18d73e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfseq18d73f = cbind(Sample.Name,Marker,dfseq18d73e)
  ##simulating RW39
  rw39samplevec = sample(rw39namevec2,1092,replace=TRUE,prob=rw39freq)
  dfrw39 = matrix(rw39samplevec,nrow=182,ncol=6)
  dfrw39b = (apply(dfrw39,1,unique)) #eliminating duplicate alleles 
  dfrw39c = ldply(dfrw39b, rbind) #converting from list to data frame
  matrw39 = as.matrix(dfrw39c)
  dfrw39d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrw39)[1]){
    testrow = matrw39[i,]
    dfrw39d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RW39",182)
  dfrw39e = as.data.frame(dfrw39d)
  colnames(dfrw39e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrw39f = cbind(Sample.Name,Marker,dfrw39e)
  ##simulating RW28
  rw28samplevec = sample(rw28namevec2,1092,replace=TRUE,prob=rw28freq)
  dfrw28 = matrix(rw28samplevec,nrow=182,ncol=6)
  dfrw28b = (apply(dfrw28,1,unique)) #eliminating duplicate alleles 
  dfrw28c = ldply(dfrw28b, rbind) #converting from list to data frame
  matrw28 = as.matrix(dfrw28c)
  dfrw28d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrw28)[1]){
    testrow = matrw28[i,]
    dfrw28d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RW28",182)
  dfrw28e = as.data.frame(dfrw28d)
  colnames(dfrw28e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrw28f = cbind(Sample.Name,Marker,dfrw28e)
  ##simulating SEQ8E8
  seq8e8samplevec = sample(seq8e8namevec2,1092,replace=TRUE,prob=seq8e8freq)
  dfseq8e8 = matrix(seq8e8samplevec,nrow=182,ncol=6)
  dfseq8e8b = (apply(dfseq8e8,1,unique)) #eliminating duplicate alleles 
  dfseq8e8c = ldply(dfseq8e8b, rbind) #converting from list to data frame
  matseq8e8 = as.matrix(dfseq8e8c)
  dfseq8e8d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfseq8e8)[1]){
    testrow = matseq8e8[i,]
    dfseq8e8d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("SEQ8E8",182)
  dfseq8e8e = as.data.frame(dfseq8e8d)
  colnames(dfseq8e8e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfseq8e8f = cbind(Sample.Name,Marker,dfseq8e8e)
  ##simulating RW56
  rw56samplevec = sample(rw56namevec2,1092,replace=TRUE,prob=rw56freq)
  dfrw56 = matrix(rw56samplevec,nrow=182,ncol=6)
  dfrw56b = (apply(dfrw56,1,unique)) #eliminating duplicate alleles 
  dfrw56c = ldply(dfrw56b, rbind) #converting from list to data frame
  matrw56 = as.matrix(dfrw56c)
  dfrw56d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrw56)[1]){
    testrow = matrw56[i,]
    dfrw56d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RW56",182)
  dfrw56e = as.data.frame(dfrw56d)
  colnames(dfrw56e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrw56f = cbind(Sample.Name,Marker,dfrw56e)
  ##simulating RWDI11
  rwdi11samplevec = sample(rwdi11namevec2,1092,replace=TRUE,prob=rwdi11freq)
  dfrwdi11 = matrix(rwdi11samplevec,nrow=182,ncol=6)
  dfrwdi11b = (apply(dfrwdi11,1,unique)) #eliminating duplicate alleles 
  dfrwdi11c = ldply(dfrwdi11b, rbind) #converting from list to data frame
  matrwdi11 = as.matrix(dfrwdi11c)
  dfrwdi11d = matrix(,nrow=182,ncol=6) #ordering alleles by size
  for(i in 1:dim(dfrwdi11)[1]){
    testrow = matrwdi11[i,]
    dfrwdi11d[i,] = testrow[order(testrow)]
  }
  Sample.Name = c(1:182)
  Marker = rep("RWDI11",182)
  dfrwdi11e = as.data.frame(dfrwdi11d)
  colnames(dfrwdi11e) = c("Allele.1","Allele.2","Allele.3","Allele.4","Allele.5","Allele.6")  
  dfrwdi11f = cbind(Sample.Name,Marker,dfrwdi11e)
  ##combining markers, exporting and importing files
  allmarkers=rbind(dfseq18d73f,dfrw39f,dfrw28f,dfseq8e8f,dfrw56f,dfrwdi11f)
  #deleting one allele per individual
  levelsa = c(levels(allmarkers$Allele.1),levels(allmarkers$Allele.2),levels(allmarkers$Allele.3),levels(allmarkers$Allele.4),levels(allmarkers$Allele.5),levels(allmarkers$Allele.6))
  levelsb=unique(levelsa)
  allmarkers$Allele.1 = factor(allmarkers$Allele.1,levels=levelsb)
  allmarkers$Allele.2 = factor(allmarkers$Allele.2,levels=levelsb)
  allmarkers$Allele.3 = factor(allmarkers$Allele.3,levels=levelsb)
  allmarkers$Allele.4 = factor(allmarkers$Allele.4,levels=levelsb)
  allmarkers$Allele.5 = factor(allmarkers$Allele.5,levels=levelsb)
  allmarkers$Allele.6 = factor(allmarkers$Allele.6,levels=levelsb)
  allmarkers2 = allmarkers
  ## randomly eliminate first allele
  rowsel = rep(c(0,182,364,546,728,910),30)
  rowsel[181:182] = c(0,182)
  for(i in 1:182){
    allmarkers2[i+rowsel[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate second allele
  rowsel2 = rep(c(182,364,546,728,910,0),30)
  rowsel2[181:182] = c(182,364)
  for(i in 1:182){
    allmarkers2[i+rowsel2[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel2[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel2[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate third allele
  rowsel3 = rep(c(364,546,728,910,0,182),30)
  rowsel3[181:182] = c(364,546)
  for(i in 1:182){
    allmarkers2[i+rowsel3[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel3[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel3[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate fourth allele
  rowsel4 = rep(c(546,728,910,0,182,364),30)
  rowsel4[181:182] = c(546,728)
  for(i in 1:182){
    allmarkers2[i+rowsel4[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel4[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel4[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate fifth allele
  rowsel5 = rep(c(728,910,0,182,364,546),30)
  rowsel5[181:182] = c(728,910)
  for(i in 1:182){
    allmarkers2[i+rowsel5[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel5[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel5[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate sixth allele
  rowsel6 = rep(c(910,0,182,364,546,728),30)
  rowsel6[181:182] = c(910,0)
  for(i in 1:182){
    allmarkers2[i+rowsel6[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel6[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel6[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ## randomly eliminate seventh allele
  rowsel7 = rep(c(0,182,364,546,728,910),30)
  rowsel7[181:182] = c(0,182)
  for(i in 1:182){
    allmarkers2[i+rowsel7[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel7[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel7[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate 8th allele
  rowsel8 = rep(c(182,364,546,728,910,0),30)
  rowsel8[181:182] = c(182,364)
  for(i in 1:182){
    allmarkers2[i+rowsel8[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel8[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel8[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate 9th allele (like 3rd)
  rowsel9 = rep(c(364,546,728,910,0,182),30)
  rowsel9[181:182] = c(364,546)
  for(i in 1:182){
    allmarkers2[i+rowsel9[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel9[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel9[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate 10th allele (like 4th)
  rowsel10 = rep(c(546,728,910,0,182,364),30)
  rowsel10[181:182] = c(546,728)
  for(i in 1:182){
    allmarkers2[i+rowsel10[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel10[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel10[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate 11th allele (like 5th)
  for(i in 1:182){
    allmarkers2[i+rowsel5[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel5[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel5[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate 12th allele (like 6th)
  for(i in 1:182){
    allmarkers2[i+rowsel6[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel6[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel6[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate 13th allele (like 1st)
  for(i in 1:182){
    allmarkers2[i+rowsel[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate 14th allele (like 2nd)
  for(i in 1:182){
    allmarkers2[i+rowsel2[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel2[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel2[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate 15th allele (like 3rd)
  for(i in 1:182){
    allmarkers2[i+rowsel3[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel3[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel3[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate 16th allele (like 4th)
  for(i in 1:182){
    allmarkers2[i+rowsel4[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel4[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel4[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate 17th allele (like 5th)
  for(i in 1:182){
    allmarkers2[i+rowsel5[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel5[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel5[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate 18th allele (like 6th)
  for(i in 1:182){
    allmarkers2[i+rowsel6[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel6[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel6[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate 19th allele (like 1st)
  for(i in 1:182){
    allmarkers2[i+rowsel[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate 20th allele (like 2nd)
  for(i in 1:182){
    allmarkers2[i+rowsel2[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel2[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel2[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate 21st allele (like 3rd)
  for(i in 1:182){
    allmarkers2[i+rowsel3[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel3[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel3[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate 22nd allele (like 4th)
  for(i in 1:182){
    allmarkers2[i+rowsel4[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel4[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel4[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  ##randomly eliminate 23rd allele (like 5th)
  for(i in 1:182){
    allmarkers2[i+rowsel5[i],ifelse(length(which(!is.na(allmarkers2[i+rowsel5[i],])))>3,sample(3:length(which(!is.na(allmarkers2[i+rowsel5[i],]))),1),8)] = NA
  }
  ## reorder alleles so NAs are last
  for(i in 1:dim(allmarkers2)[1]){
    alleles = allmarkers2[i,3:8]
    allmarkers2[i,3:8] = alleles[order(alleles)]
  }
  write.table(allmarkers2,namevec4[n],sep = "\t",quote = FALSE, row.names=FALSE)
  simdata = read.GeneMapper("/Users/Lakshmi/Desktop/Analysis/simdata.txt")
  Usatnts(simdata) <- c(3,4,4,2,4,2)
  ploidycol = rep(6,length(Samples(simdata)))
  ploidymat = as.matrix(cbind(ploidycol,ploidycol,ploidycol,ploidycol,ploidycol,ploidycol))
  colnames(ploidymat) = c("SEQ18D73","RW39","RW28","SEQ8E8","RW56","RWDI11") 
  Ploidies(simdata) = ploidymat
  simdist = meandistance.matrix(simdata,distmetric = Bruvo.distance)
  simclones = assignClones(simdist,threshold=0.2)
  matches = length(simclones)-length(unique(simclones))
  matchvec23[n] = matches
  n = n+1
}

# Figure 2 -----

nulldata = read.table("nullalleletrials.txt",header=TRUE)

postscript(file="figure2.eps")
par(mar=c(4,4,4,4))
plot(matchcount~Trial,data=nulldata,type="l",ylab="",xlab="Rounds of deletions",lty=3)
lines(actualdel~Trial,data=nulldata,type="l",lty=1)
lines(simcount~Trial,data=nulldata,type="l",lty=2)
legend(1,170,c("Average deletions","Percent simulations with false positives","Total false positives in 100 simulations"),lty=c(1,2,3))
dev.off()

##Figure 2 revision

nulldata = read.table("nullalleletrials.txt",header=TRUE)
nulldata$percdel = (nulldata$actualdel/30)*100
nulldata$percsim = (nulldata$simcount/100)*100
nulldata$percfp = (nulldata$matchcount/(100*181))*100

postscript(file="figure2rev4.eps")
par(mar=c(4,5,4,4))
plot(percsim~Trial,data=nulldata,type="l",ylab="Percent of max",xlab="Rounds of deletions",lty=2,las=1,lwd=3,cex.axis=1.5,cex.lab=1.5)
lines(percdel~Trial,data=nulldata,type="l",lwd=3,lty=1)
lines(percfp~Trial,data=nulldata,type="l",lty=3,lwd=3)
legend(1,75,c("Average deletions","Simulations with false positives","Number of false positives"),lty=c(1,2,3),lwd=3,cex=1.5)
dev.off()
