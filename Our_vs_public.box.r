setwd('D:/puyuan/Project/PMD/2018.PMD3/2018.PMD3/2018.Oct24.PMD3_qiime2_lish/')
library(ggplot2)
library(ggpubr)
library(cowplot)
library(xlsx)
library(reshape2)
library(RColorBrewer)

species<-function(ii){
  ii<-as.character(ii)
  #ii[intersect(grep('.*[|,.;]s__.*[|,.;]t__..*',ii),grep('[|,.;]t__$',ii,invert=T))]<-as.character(lapply(ii[intersect(grep('.*[|,.;]s__.*[|,.;]t__..*',ii),grep('[|,.;]t__$',ii,invert=T))],function(x){gsub('.*[|,.;]t','t',x)}))
  #ii[intersect(grep('.*[|,.;]g__.*[|,.;]s__..*',ii),grep('s__[|,.;]t__',ii,invert=T))]<-as.character(lapply(ii[intersect(grep('.*[|,.;]g__.*[|,.;]s__..*',ii),grep('s__[|,.;]t__',ii,invert=T))],function(x){gsub('.*[|,.;]s','s',x)}))
  #ii[intersect(grep('.*[|,.;]f__.*[|,.;]g__..*',ii),grep('g__[|,.;]s__',ii,invert=T))]<-as.character(lapply(ii[intersect(grep('.*[|,.;]f__.*[|,.;]g__..*',ii),grep('g__[|,.;]s__',ii,invert=T))],function(x){gsub('.*[|,.;]g','g',x)}))
  ii[intersect(grep('.*[|,.;]o__.*[|,.;]f__..*',ii),grep('f__[|,.;]g__',ii,invert=T))]<-as.character(lapply(ii[intersect(grep('.*[|,.;]o__.*[|,.;]f__..*',ii),grep('f__[|,.;]g__',ii,invert=T))],function(x){gsub('.*[|,.;]f','f',x)}))
  ii[intersect(grep('.*[|,.;]c__.*[|,.;]o__..*',ii),grep('o__[|,.;]f__',ii,invert=T))]<-as.character(lapply(ii[intersect(grep('.*[|,.;]c__.*[|,.;]o__..*',ii),grep('o__[|,.;]f__',ii,invert=T))],function(x){gsub('.*[|,.;]o','o',x)}))
  ii[intersect(grep('.*[|,.;]p__.*[|,.;]c__..*',ii),grep('c__[|,.;]o__',ii,invert=T))]<-as.character(lapply(ii[intersect(grep('.*[|,.;]p__.*[|,.;]c__..*',ii),grep('c__[|,.;]o__',ii,invert=T))],function(x){gsub('.*[|,.;]c','c',x)}))
  ii[intersect(grep('k__.*[|,.;]p__..*',ii),grep('k__[|,.;]p__',ii,invert=T))]<-as.character(lapply(ii[intersect(grep('k__.*[|,.;]p__..*',ii),grep('k__[|,.;]p__',ii,invert=T))],function(x){gsub('.*[|,.;]p','p',x)}))
  return(ii)
}

inguang = 'D:/puyuan/Project/Guangdong/2019.Guangdong/2019.Jan3.Guangdong/data/profile/profile/L6.txt'
setd = c('p__OD1;c__;o__;f__;g__','Unassigned;__;__;__;__;__','k__Bacteria;__;__;__;__;__','o__RF39;f__;g__','o__Clostridiales;__;__')
Na_color = 'analysis/8.PMD_vs_Guangdong/PMDvsGuangdong.Name_color'
Na_color<- read.csv(Na_color,sep='\t',stringsAsFactors = F,check.names = F)

prof <- read.csv('data/profile/L6.txt',sep = '\t',check.names = F,row.names = 1)
map<- read.csv('./data/mapping.file/mapping3.file',sep='\t',stringsAsFactors = F,check.names = F)
map2 <- read.csv('./data/mapping.file/mapping5.file',sep='\t',stringsAsFactors = F,check.names = F)
map[,colnames(map2)[2]]<- map2[,2] 
prof <- prof[,map[,1]]
prof <- prof[which(rowMeans(prof)>=0.005),]
#prof <- prof[order(-rowMeans(prof)),]
prof <- prof[order(-as.numeric(apply(prof,1,median))),]
pp <- data.frame(ID=row.names(prof),prof,Group='PMD',stringsAsFactors = F)
pp <- melt(pp,id=c('ID','Group'))
pp$Name <- species(pp$ID)
pp$GG <- gsub('.*;p__','',pp$ID)
pp$GG <- gsub(';.*','',pp$GG)
#pp <- merge(pp,Na_color,by = 'Name')
pp$Name <- factor(pp$Name,unique(pp$Name))


prof_g <- read.csv(inguang,sep = '\t',check.names = F,row.names = 1)
prof_g <- prof_g[which(rowMeans(prof_g)>=0.005),]
#prof_g <- prof_g[order(-rowMeans(prof_g)),]
prof_g <- prof_g[order(-as.numeric(apply(prof_g,1,median))),]
pg <- data.frame(ID=row.names(prof_g),prof_g,Group='Guangdong')
pg <- melt(pg,id=c('ID','Group'))
pg$Name <- species(pg$ID)
pg$GG <- gsub('.*;p__','',pg$ID)
pg$GG <- gsub(';.*','',pg$GG)
pg$GG[pg$GG=="k__Bacteria"] <- "k__Bacteria;Others" 
pg <- pg[!pg$Name%in%setd,]
#pg <- merge(pg,Na_color,by = 'Name')
pg$Name <- factor(pg$Name,unique(pg$Name))

Color = data.frame(GG=sort(union(unique(pg$GG),unique(pp$GG))),Color=brewer.pal(8,'Paired')[c(3,6,5,1,4,2)],stringsAsFactors = F)
pp = merge(pp,Color,by = 'GG',all.x=T)
pg = merge(pg,Color,by = 'GG',all.x=T)
pp$GG <- factor(pp$GG,unique(pp$GG)) 
pg$GG <-  factor(pg$GG,unique(pg$GG)) 

a = unique(c(as.character(pp$ID),as.character(pg$ID)))
#write.table(a,'analysis/8.PMD_vs_Guangdong/PMDvsGuangdong1_new.rename',sep='\t',row.names = F,quote = F)
a = read.csv('analysis/8.PMD_vs_Guangdong/PMDvsGuangdong1_new.rename.f',sep='\t',stringsAsFactors = F,check.names = F,header = F)

pp = merge(pp,a,by.x = 'ID',by.y='V1',all.x = T)
pp = pp[order(as.numeric(pp$Name)),]
pp$V2 = factor(pp$V2,unique(pp$V2)) 

pg = merge(pg,a,by.x = 'ID',by.y='V1',all.x = T)
pg = pg[order(as.numeric(pg$Name)),]
pg$V2 = factor(pg$V2,unique(pg$V2))

D = pp

p1 <- ggplot(data = D,
       aes(x=V2,
           y=value)) +
  geom_boxplot(aes(fill=GG),size = 0.2, outlier.size = 0.1, alpha = 0.8)+
  scale_y_sqrt(limits=c(0,1))+
                   
  labs(title="", x="", y="Relative Abundance")+
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title=element_text(
      family="Helvetica",
      face="bold",
      colour="black",
      size=6
    ),
    axis.text.y=element_text(
      family="Helvetica",
      colour="black",
      size=6
    ),
    axis.text.x=element_text(
      family="Helvetica",
      colour="black",
      size=6,
      vjust=1,
      hjust=1,
      angle = 45
    ),
    axis.line=element_blank(),   
    strip.text.y =element_text(  
      family="Helvetica",
      colour="black",
      angle=180,
      vjust=0.5,
      hjust=1,
      size=10
    ),
    strip.background=element_rect(  
      size=0.3,
      colour = 'black'
    ),
    panel.border=element_rect(  
      fill=NA,
      colour="black",
      size=0.3
    ),
    legend.title=element_text(
      family="Helvetica",
      colour="black",
      size=10
    ),
    legend.text=element_text(
      family="Helvetica",
      colour="black",
      size=8
    )
  ) +
  guides(fill = guide_legend("",ncol=1), color = F) + 
  scale_fill_manual(values = unique(D$Color))#+


plot_grid(p1,p2,nrow = 2)
