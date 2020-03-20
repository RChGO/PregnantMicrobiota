
library(vegan) 
library(ape)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)

vector_fd<-function(tab,sam){  
  ms<-max(sam[,1]^2+sam[,2]^2)
  mt<-max(tab[,1]^2+tab[,2]^2)
  return(sqrt(ms)/sqrt(mt))
}

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

#inf='D:/puyuan/Project/RA/2018.RA_tonsil/2018.Oct28.MAG.RA_tonsil/data/profile/20181112.GeneSet.kegg.module.LD.prof'
#inp='E:/puyuan/Project/RA/2018.RA_tonsil/2018.Oct28.MAG.RA_tonsil/data/mapping.file/mapping.file'
inf='D:/puyuan/Project/PMD/2018.PMD3/2018.PMD3/2018.Oct24.PMD3_qiime2_lish/data/profile/L6.txt'
inp='D:/puyuan/Project/PMD/2018.PMD3/2018.PMD3/2018.Oct24.PMD3_qiime2_lish/data/mapping.file/mapping3.file'
inc <- "Gestation_age_G2"
Axis1 <- 1
Axis2 <- 2
xx <- 1
yy <- 1

prof <- read.csv(file=inf,sep="\t",row.names=1,check.names = F,stringsAsFactors = F)#,skip = 1)  
map <- read.csv(file=inp,sep="\t",check.names = F,stringsAsFactors = F)
prof <- apply(prof,2,function(x){x/sum(x)})
#map = map[map$GA!=0,]
#map = map[map$GA!=1,]
map <- data.frame(ID=map[,1],Group=map[,inc],stringsAsFactors = F)
map<- map[order(map[,1]),]
prof <- t(prof[,map[,1]])
prof <- prof[,colSums(prof)!=0]
prof2 <- sqrt(prof)

if(length(which(plyr::count(map$Group)$freq>3))<1){
  cat("\nWarning: There are not enough samples!\n\n")
  q()
}

ord <- capscale(prof2~Group,map,distance = 'bray') #dbRDA分析,'canberra'距离相比'bray'往往更明�?
#eig <- summary(eigenvals(ord))[[1]]
eig <- summary(eigenvals(ord))



p<-plot(ord,choices = c(Axis1,Axis2))
site.scores <- as.data.frame(scores(p, "site")) 
Axis_name<-colnames(site.scores)
site.scores$lab <- row.names(site.scores)
site.scores <- merge(site.scores,map[,c(1,2)],by.x = 'lab',by.y = 'ID',on='left')
#site.scores$z <- NA 
colnames(site.scores)[2:3]<- c('A1','A2')

species.scores <- as.data.frame(scores(p, "species")) 
or<-species.scores[,1]^2+species.scores[,2]^2
species.scores <-species.scores[head(order(or,decreasing = T),6),]

species.scores$lab <- species(row.names(species.scores))
chang <- vector_fd(species.scores[,c(1,2)],site.scores[,c(2,3)])*0.1
colnames(species.scores)[1:2]<- c('A1','A2')

fig<-ggplot(data = site.scores) +
  #geom_text_repel(data=species.scores, aes(x=A1*chang, y=A2*chang, label=lab),
  #                family="Helvetica", fontface="italic", size=4, check_overlap=TRUE) +
  geom_text(data=species.scores, aes(x=A1*chang, y=A2*chang, label=lab),
                  family="Helvetica", fontface="italic", size=4) +
  #geom_segment(data=species.scores, aes(x=0, y=0, xend=A1*chang, yend=A2*chang),
  #             arrow = arrow(length = unit(0.3, "cm")), size=0.8, alpha=0.5)+
  geom_hline(yintercept = 0,color='grey')+
  geom_vline(xintercept = 0,color='grey')+
  stat_ellipse(aes(x=A1, y=A2,fill=as.character(Group)),size=1, geom="polygon", level=0.8, alpha=0.3) +
  geom_point(size=2.5,aes(x=A1,y=A2,fill=as.character(Group)),alpha=0.6,color='black',shape=21)+
  theme_bw()+
  guides(fill=guide_legend(title = 'Group'),color=guide_legend(title = 'Group'))+
  #scale_color_manual(values = c('#00468B','#ED0000'))+
  #scale_fill_manual(values = c('#00468B','#ED0000'))+
  #coord_cartesian(ylim = c(-5, 5),xlim = c(-5, 5))+
  labs(title="dbRDA", 
       x=paste(Axis_name[1]," (",as.character(round(eig[2,Axis1],4)*100),'%)'), 
       y=paste(Axis_name[2]," (",as.character(round(eig[2,Axis2],4)*100),'%)')
  )+
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5,
                              size=20),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_text(size=12)
  )+
  xlim(min(site.scores[,'A1'],species.scores[,'A1']*chang, 0)*xx, max(site.scores[,'A1'],species.scores[,'A1']*chang)*xx)+
  ylim(min(site.scores[,'A2'],species.scores[,'A2']*chang, 0)*yy, max(site.scores[,'A2'],species.scores[,'A2']*chang)*yy)

fig

#######################################################################################################
rda_summ = summary(ord)
site.scores = rda_summ$sites

library(shapes)
library(ggplot2)
library(reshape2)

#A <- as.matrix(read.table("species-site.txt",head=T))
#B <- as.matrix(read.table("composition.site.txt",head=T))
#B <- B[rownames(A),]
map <- read.csv(file='D:/puyuan/Project/PMD/2018.PMD3/2018.PMD3/2018.Oct24.PMD3_qiime2_lish/data/mapping.file/mapping3.file',
                sep="\t",check.names = F,stringsAsFactors = F)
#map<- map[!is.na(map$`16S_sampleID`),]
AA <- site.scores
AA <- AA[map[,1],]


BB <- site.scores
BB <- BB[map[,1],]
#
A = scale(AA)
B = scale(BB)
####
p <- protest(B,A,symmetric = TRUE) #procrustes(B,A)
p
pt <- data.frame(rda1=p$Yrot[,1], rda2=p$Yrot[,2],
                 xrda1=p$X[,1], xrda2=p$X[,2])
pt$SampleID <- row.names(pt)
pt <- merge(pt,map,by.x = 'SampleID',by.y = '#SampleID',all.x = T )
colnames(pt)[ncol(pt)] = 'Group'


g<-c("1-12W","13-16W","17-20W","21-24W","25-28W","29-32W","33-36W","37-45W")
for(x in 1:length(g)){
  pt$Group[pt$Group == as.character(x-1)]<-g[x]
}
pt$Group = as.character(pt$Group)

ggplot(pt) +
  geom_segment(aes(x=rda1,y=rda2,xend=xrda1,yend=xrda2),color='grey60',arrow=arrow(length=unit(0.2,"cm")),size=0.01)+
  geom_point(aes(x=rda1, y=rda2,fill=Group),color="black",shape=21,alpha=0.6,size=2.5) +
  geom_point(aes(x=xrda1, y=xrda2,fill=Group),color="black",shape=22,alpha=0.6,size=2.5) +
  #geom_hline(yintercept = 0,color='grey')+
  #geom_vline(xintercept = 0,color='grey')+
  theme_bw()+
  guides(fill=guide_legend(title = 'Group'))+
  #scale_color_manual(values = c('#00468B','#ED0000'))+
  #scale_fill_manual(values = c('#00468B','#ED0000'))+
  #coord_cartesian(ylim = c(-5, 5),xlim = c(-5, 5))+
  labs(title="Procrustes errors", 
       x='Dimension1', 
       y='Dimension2')+
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5,
                              size=16),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_text(size=12)
  )
  
p$signif
p$permutations

#####

ans <- procOPA(A,B) 
plotshapes(A,B,joinline=1:nrow(A))
plotshapes(ans$Ahat,ans$Bhat,joinline=1:nrow(A))

