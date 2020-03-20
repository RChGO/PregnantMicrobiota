setwd('D:/puyuan/Project/PMD/2018.PMD3/2018.PMD3/2018.Oct24.PMD3_qiime2_lish')
library("fdrtool")
library(ggplot2)
library(ggrepel)
library(cowplot)
library(xlsx)
wtest2<-function(tab,map,group,p=F,sample='a'){
  if (p==T){  
    fil<-plyr::count(map[,sample])
    fil$x<-as.character(fil$x)
    fil<-fil$x[which(fil$freq==2)]
    map<-map[which(map[,sample]%in%fil),]
    map<-map[order(map[,group],map[,sample]),]
  }
  #row.names(tab)<-gsub('.*s__','',row.names(tab))
  tab<-tab[,map[,1]]
  g<-unique(map[,group])
  table<-data.frame(ID=row.names(tab),g1=NA,g2=NA,pvalue=NA,stringsAsFactors = F)
  table$g1<-apply(tab,1,function(x){mean(x[map[,group]==g[1]])})
  table$g2<-apply(tab,1,function(x){mean(x[map[,group]==g[2]])})
  if(p==T){
    table$pvalue<-apply(tab,1,function(x){
      wilcox.test(as.numeric(x[map[,group]==g[1]]),as.numeric(x[map[,group]==g[2]]),paired = T)$p.value
    })
  }else{
    table$pvalue<-apply(tab,1,function(x){wilcox.test(as.numeric(x[map[,group]==g[1]]),as.numeric(x[map[,group]==g[2]]))$p.value})
  }
  colnames(table)[2:3]<-g
  return(table)
}
species<-function(ii){
  ii<-as.character(ii)
  #ii[intersect(grep('.*[|,.;]s__.*[|,.;]t__..*',ii),grep('[|,.;]t__$',ii,invert=T))]<-as.character(lapply(ii[intersect(grep('.*[|,.;]s__.*[|,.;]t__..*',ii),grep('[|,.;]t__$',ii,invert=T))],function(x){gsub('.*[|,.;]t','t',x)}))
  #ii[intersect(grep('.*[|,.;]g__.*[|,.;]s__..*',ii),grep('s__[|,.;]t__',ii,invert=T))]<-as.character(lapply(ii[intersect(grep('.*[|,.;]g__.*[|,.;]s__..*',ii),grep('s__[|,.;]t__',ii,invert=T))],function(x){gsub('.*[|,.;]s','s',x)}))
  ii[intersect(grep('.*[|,.;]f__.*[|,.;]g__..*',ii),grep('g__[|,.;]s__',ii,invert=T))]<-as.character(lapply(ii[intersect(grep('.*[|,.;]f__.*[|,.;]g__..*',ii),grep('g__[|,.;]s__',ii,invert=T))],function(x){gsub('.*[|,.;]g','g',x)}))
  ii[intersect(grep('.*[|,.;]o__.*[|,.;]f__..*',ii),grep('f__[|,.;]g__',ii,invert=T))]<-as.character(lapply(ii[intersect(grep('.*[|,.;]o__.*[|,.;]f__..*',ii),grep('f__[|,.;]g__',ii,invert=T))],function(x){gsub('.*[|,.;]f','f',x)}))
  ii[intersect(grep('.*[|,.;]c__.*[|,.;]o__..*',ii),grep('o__[|,.;]f__',ii,invert=T))]<-as.character(lapply(ii[intersect(grep('.*[|,.;]c__.*[|,.;]o__..*',ii),grep('o__[|,.;]f__',ii,invert=T))],function(x){gsub('.*[|,.;]o','o',x)}))
  ii[intersect(grep('.*[|,.;]p__.*[|,.;]c__..*',ii),grep('c__[|,.;]o__',ii,invert=T))]<-as.character(lapply(ii[intersect(grep('.*[|,.;]p__.*[|,.;]c__..*',ii),grep('c__[|,.;]o__',ii,invert=T))],function(x){gsub('.*[|,.;]c','c',x)}))
  ii[intersect(grep('k__.*[|,.;]p__..*',ii),grep('k__[|,.;]p__',ii,invert=T))]<-as.character(lapply(ii[intersect(grep('k__.*[|,.;]p__..*',ii),grep('k__[|,.;]p__',ii,invert=T))],function(x){gsub('.*[|,.;]p','p',x)}))
  return(ii)
}

#wb <- createWorkbook()
#P <- c()
tax1 = 'data/data_generation/qiime2/s06/rep-seqs.tax/taxonomy.tsv'
tax = 'data/mapping.file/mapping.Feature'
#tax = 'analysis/4.abundance.sig/diff_sig/moduleLD.annote'
mapall = 'data/mapping.file/mapping11.file' 

inmap = 'data/mapping.file/mapping13.file'
infile = 'data/profile/Feature.txt'

#Group = 'Age_group'
#Group = 'GWG'
#Group = 'Thalassemia'
#Group ='Hepatopathy'
#Group ='Thyroid_disease'
#Group = 'GDM'
#Group = 'Native_place'
Group='BMI_Group'
cutoff = -log10(0.01)
FD = 5
func_xlim = seq(-7,24,6) 
Test = 'pvalue'

########
prof <- read.csv(infile,sep='\t',row.names = 1,check.names = F,stringsAsFactors = F)
map_all <- read.csv(mapall,sep = '\t',check.names = F,stringsAsFactors = F)
prof <- prof[,map_all[,1]]
ff <-rowMeans(prof)
prof <- prof[ff>0.0001,]
#####################
tax1 <- read.csv(tax1,sep = '\t',check.names = F,stringsAsFactors = F)
#tax1<- data.frame(ID=)
tax1$Taxon <- unlist(lapply(tax1$Taxon,function(x){gsub(" ",'',x)}))
tax1$Taxon <- species(tax1$Taxon)
colnames(tax1)[1] <- 'ID'
tax1 = tax1[,-3]
colnames(tax1) <- c('ID','Taxonomy')

tax <- read.csv(tax,sep = '\t',check.names = F,stringsAsFactors = F)
colnames(tax) <- c('ID','Taxon')
##########################
map <- read.csv(inmap,sep = '\t',check.names = F,stringsAsFactors = F)
map[,Group] = as.character(map[,Group])
map <- map[order(map[,Group]),]
prof <- prof[,map[,1]]
#prof <- prof[rowSums(prof)!=0,]
mea = data.frame(ID=row.names(prof),Mean=rowMeans(prof))

tt <- wtest2(prof,map,Group)
#tt$qvalue <- NA
#tt$qvalue[!is.na(tt$pvalue)] = fdrtool(tt$pvalue[!is.na(tt$pvalue)],statistic="pvalue")$qval
tt$Fold_Change <- tt[,2]/tt[,3]
tt$tranF <- log2(tt$Fold_Change)
tt$tranQ <- -log10(tt[,Test])

tt$Group <- colnames(tt)[2]
tt$Group[tt$tranF<0] <- colnames(tt)[3]
tt$Group[tt$tranQ < cutoff] <- 'PASS'
tt <- merge(tt,mea,by = 'ID',all.x = T)
tt$Mean[tt$Group==colnames(tt)[2]] <- tt[tt$Group==colnames(tt)[2],2]
tt$Mean[tt$Group==colnames(tt)[3]] <- tt[tt$Group==colnames(tt)[3],3]
tt$Group = factor(tt$Group,c(colnames(tt)[2],colnames(tt)[3],'PASS'))
tt$tranF[tt$tranF>=FD] <- FD
tt$tranF[tt$tranF<= -FD] <- -FD
#############################
#tt$Taxon = tt$ID
tt <- merge(tt,tax,by = 'ID',all.x = T)
tt <- merge(tt,tax1,by = 'ID',all.x = T)
#tt$Taxon[abs(tt$tranQ) < cutoff] <- NA
dat = data.frame(FeatureID=tt$Taxon,Taxonomy=tt$Taxonomy,G1=tt[,2],G2=tt[,3],'wilcox pvalue'=tt$pvalue,'Fold change'=tt$Fold_Change,check.names = F)
colnames(dat)[3:4] <- colnames(tt)[2:3]
#write.table(dat,paste('analysis/4.abundance.sig/diff_sig/table/coreFeature.',Group,'.txt',sep='')
#            ,sep='\t',row.names = F,quote = F)

#sheet  <- createSheet(wb, sheetName= paste(Group,'CoreFeature',sep='.'))
#addDataFrame(x = dat, sheet ,
#             row.names = FALSE)
#saveWorkbook(wb, "analysis/4.abundance.sig/diff_sig/all_test.volcano.xlsx") 
#####################
mycolor<-c('red','steelblue','Forest green')

p1<-ggplot(tt,
       aes(x=tranF,y=tranQ,color=Group)
       #aes(x=FD*tranF,y=tranQ,color=Group)
       )+  
  geom_point(aes(size=Mean,fill=Group),color='black',shape=21,alpha=0.6)+
  geom_text_repel(aes(label=Taxon),size=2.5)+
  geom_vline(xintercept=0,linetype="dashed")+
  geom_hline(yintercept=cutoff,linetype="dashed")+
  #geom_point(aes(size=Mean),shape=1)+
  theme_minimal()+
  #theme_bw()+
  theme(
    #panel.grid =element_blank(),
    #panel.border = element_blank(),
    axis.line = element_line(
      #size=1,
      colour = "black"),
    axis.text=element_text(
      family="Arial",
      #face="bold",
      colour="black",
      size=10
    ),
    axis.title =element_text(
      family="Arial",
      #face="bold",
      colour="black",
      size=12
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )+
  scale_fill_manual(values = mycolor)+
  scale_color_manual(values = mycolor)+
  #scale_fill_brewer(palette = 'Paired')+
  #scale_x_continuous(trans='sqrt')+
  guides(label=F)+
  #xlim(-6,6)+
  xlim(-FD,FD)+
  xlab('log2 ( Fold change )')+
  #xlab(paste(as.character(FD),'*log2 ( Fold change )',sep=''))+
  ylab('-log10 ( wilcox pvalue)')
#dev.off()
#p1

###############################################################################################
#tax = 'data/data_generation/qiime2/s06/rep-seqs.tax/taxonomy.tsv'
#tax = 'data/mapping.file/mapping.Feature'
tax = 'analysis/4.abundance.sig/diff_sig/moduleLD.annote'
#mapall = 'data/mapping.file/mapping11.file'  

#inmap = 'data/mapping.file/mapping12.file'
infile = 'data/profile/kegg.module.LD.prof'
#Group = 'Age_group'
#Group = 'GWG'
#Group = 'Thalassemia'
#Group ='Hepatopathy'
#Group ='Thyroid_disease'
#Group = 'GDM'
#Group = 'Native_place'
cutoff = -log10(0.01)
FD = 20
Test = 'pvalue'

########
prof <- read.csv(infile,sep='\t',row.names = 1,check.names = F,stringsAsFactors = F)
map_all <- read.csv(mapall,sep = '\t',check.names = F,stringsAsFactors = F)
prof <- prof[,map_all[,1]]
ff <-rowMeans(prof)
prof <- prof[ff>0.0001,]
#####################
#tax <- read.csv(tax,sep = '\t',check.names = F,stringsAsFactors = F)
#tax <- data.frame(ID=)
#tax$Taxon <- unlist(lapply(tax$Taxon,function(x){gsub(" ",'',x)}))
#tax$Taxon <- species(tax$Taxon)
#colnames(tax)[1] <- 'ID'
#tax = tax[,-3]

tax <- read.csv(tax,sep = '\t',check.names = F,stringsAsFactors = F)
colnames(tax) <- c('ID','Taxon')
##########################
map <- read.csv(inmap,sep = '\t',check.names = F,stringsAsFactors = F)
map[,Group] = as.character(map[,Group])
prof <- prof[,map[,1]]
#prof <- prof[rowSums(prof)!=0,]
mea = data.frame(ID=row.names(prof),Mean=rowMeans(prof))

tt <- wtest2(prof,map,Group)
#tt$qvalue <- NA
#tt$qvalue[!is.na(tt$pvalue)] = fdrtool(tt$pvalue[!is.na(tt$pvalue)],statistic="pvalue")$qval
tt$Fold_Change <- tt[,2]/tt[,3]
tt$tranF <- log2(tt$Fold_Change)
tt$tranQ <- -log10(tt[,Test])

tt$Group <- colnames(tt)[2]
tt$Group[tt$tranF<0] <- colnames(tt)[3]
tt$Group[tt$tranQ < cutoff] <- 'PASS'
tt <- merge(tt,mea,by = 'ID',all.x = T)
tt$Mean[tt$Group==colnames(tt)[2]] <- tt[tt$Group==colnames(tt)[2],2]
tt$Mean[tt$Group==colnames(tt)[3]] <- tt[tt$Group==colnames(tt)[3],3]
tt$Group = factor(tt$Group,c(colnames(tt)[2],colnames(tt)[3],'PASS'))
tt$tranF[tt$tranF>=FD] <- FD
tt$tranF[tt$tranF<= -FD] <- -FD
#############################
#tt$Taxon = tt$ID
tt <- merge(tt,tax,by = 'ID',all.x = T)
#tt$Taxon[abs(tt$tranQ) < cutoff] <- NA
dat = data.frame(ModuleID=tt$ID,Module_annot=tt$Taxon,G1=tt[,2],G2=tt[,3],'wilcox pvalue'=tt$pvalue,'Fold change'=tt$Fold_Change,check.names = F)
colnames(dat)[3:4] <- colnames(tt)[2:3]
#write.table(tt,paste('analysis/4.abundance.sig/diff_sig/table/moduleLD.',Group,'.txt',sep=''),
#            sep='\t',row.names = F,quote = F)

sheet  <- createSheet(wb, sheetName= paste(Group,'ModuleLD',sep='.'))
addDataFrame(x = dat, sheet ,
             row.names = FALSE)
saveWorkbook(wb, "analysis/4.abundance.sig/diff_sig/all_test.volcano.xlsx") 
#####################
mycolor<-c('red','steelblue','Forest green')

p2<-ggplot(tt,
       #aes(x=tranF,y=tranQ,color=Group)
       aes(x=FD*tranF,y=tranQ,color=Group)
)+  
  geom_point(aes(size=Mean,fill=Group),color='black',shape=21,alpha=0.6)+
  geom_text_repel(aes(label=Taxon),size=2.5)+
  geom_vline(xintercept=0,linetype="dashed")+
  geom_hline(yintercept=cutoff,linetype="dashed")+
  #geom_point(aes(size=Mean),shape=1)+
  theme_minimal()+
  #theme_bw()+
  theme(
    #panel.grid =element_blank(),
    #panel.border = element_blank(),
    axis.line = element_line(
      #size=1,
      colour = "black"),
    axis.text=element_text(
      family="Arial",
      #face="bold",
      colour="black",
      size=10
    ),
    axis.title =element_text(
      family="Arial",
      #face="bold",
      colour="black",
      size=12
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
    #panel.grid.minor.y = element_line(  
    #  linetype = 'dashed',
    #  colour='grey'
    #),
  )+
  scale_fill_manual(values = mycolor)+
  scale_color_manual(values = mycolor)+
  #scale_fill_brewer(palette = 'Paired')+
  #scale_x_continuous(trans='sqrt')+
  guides(label=F)+
  #xlim()+
  scale_x_continuous(breaks = func_xlim)+
  #xlim(-FD,FD)+
  #xlab('log2 ( Fold change )')+
  xlab(paste(as.character(FD),'*log2 ( Fold change )',sep=''))+
  ylab('-log10 ( wilcox pvalue)')
#dev.off()
#p2
#print(summary(20*tt$tranF))

#P[[Group]]<-plot_grid(p1,p2,labels = Group)
#P[[Group]]

#plot_grid(plotlist = P,ncol = 1)
