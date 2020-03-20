setwd('D:/puyuan/Project/PMD/2018.PMD3/2018.PMD3/2018.Oct24.PMD3_qiime2_lish/')
library(ggplot2)
library(cowplot)
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
corr_1 = function(matr){
  cor_r = matrix(rep(0,ncol(matr)*ncol(matr)),ncol(matr),ncol(matr))
  row.names(cor_r) = colnames(matr)
  colnames(cor_r) = colnames(matr)
  cor_p = matrix(rep(0,ncol(matr)*ncol(matr)),ncol(matr),ncol(matr))
  row.names(cor_p) = colnames(matr)
  colnames(cor_p) = colnames(matr)
  for(i in 1:ncol(matr)){
    for(j in 1:ncol(matr)){
      if (length(which(matr[,i]!=0))>1 & length(which(matr[,j]!=0))>1){
        cor_1 = cor.test(matr[,i],matr[,j],method = "spearman")
        cor_r[i,j] = cor_1$estimate
        cor_p[i,j] = cor_1$p.value
      }else {
        cor_r[i,j]<-0
        cor_p[i,j]<-1
      }
    }
  }
  k = list(p = cor_p,r = cor_r)
  return(k)
}
cor_filter_diff <- function(r,p,rlim=0.3,plim=0.05){
  #r <- corr$r
  #p <- corr$p
  r_f <- r
  r_f[p>plim]<-NA
  r_f[r_f<rlim & r_f> -rlim]<-NA
  fr <- apply(r_f,1,function(x){length(which(is.na(x)==F))})
  fc <- apply(r_f,2,function(x){length(which(is.na(x)==F))})
  r <- r[fr!=1,fc!=1]
  p <- p[fr!=1,fc!=1]
  return(list(r=r,p=p))
}
reshape_dat<-function(matr){
  #matr<-p
  matr<-data.frame(ID=row.names(matr),matr,stringsAsFactors = F,check.names = F)
  matr_n<-reshape2::melt(matr,id='ID')
  i=1
  num<-c()
  for (x in 1:ncol(matr)){
    for (y in (1+(x-1)*nrow(matr)):(x+(x-1)*nrow(matr))){
      num<-c(num,y)
    }
  }
  matr_n<-matr_n[-num,]
  return(matr_n)
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
mycolor = c("#5893c8","#e27933","#a5a5a5","#f4ba1b","#416cb4","#6ea647","#245a8a","#984823","#636363","#947226","#254173","#426630","#79a8d3","#e89259",
            "#b8b7b7","#f5c63b","#6587c4","#88b966","#3476b5","#c95e17","#848485","#c49611","#325697","#588638","#9abdde","#ecad81","#c9c9c8","#f6d368",
            "#8ba5d3","#a4ca8b","#1e4c75","#803c1f","#525252","#7b5f26","#1e3661","#37552a","#89b2d9","#e9a06d","#c0c0bf","#f5cc50","#7896cb","#95c179",
            "#2a689f","#b05321","#737474","#ad831e","#2a4c85","#4d7634")

prof <- read.csv('data/profile/L6.txt',sep = '\t',check.names = F,row.names = 1)
ma1 <- read.csv('data/mapping.file/mapping4.file',check.names = F,sep='\t',stringsAsFactors = F)
ma2 <- read.csv('data/mapping.file/mapping3.file',check.names = F,sep='\t',stringsAsFactors = F)
ma3 <- read.csv('data/mapping.file/mapping5.file',check.names = F,sep='\t',stringsAsFactors = F)
ma3$Gestation_age_G2 <- ma2$Gestation_age_G2
map <- ma3
map <- merge(map,ma1,by='#SampleID',all.x = T)

#map$Gestation_age[map$Gestation_age<=10]<-10
prof <- prof[,map[,1]]
row.names(prof) <- species(row.names(prof))

#f <- apply(prof,1,function(x){length(which(x>0.001))})
f <- rowMeans(prof)
prof <- prof[f>0.005,]

group<- sort(unique(map$Gestation_age_G2))
tt <- c()
for(x in 1:3){
  for(y in x:3){
    if (x != y){
      m <- map[map$Gestation_age_G1%in%c(group[x],group[y]),]
      p<- prof[,m[,1]]
      nam <- paste(as.character(group[x]),as.character(group[y]),sep=' vs ') 
      tt[[nam]] <- wtest2(p,m,'Gestation_age_G1')
    }
  }
}
#View(tt$`0 vs 1`)

#a<-which(tt$`0 vs 1`$`0`/tt$`0 vs 1`$`1`>2 | tt$`0 vs 1`$`0`/tt$`0 vs 1`$`1`< -1/2)
#b<-which(tt$`0 vs 2`$`0`/tt$`0 vs 2`$`2`>2 | tt$`0 vs 2`$`0`/tt$`0 vs 2`$`2`< -1/2)
#c<-which(tt$`1 vs 2`$`1`/tt$`1 vs 2`$`2`>2 | tt$`1 vs 2`$`1`/tt$`1 vs 2`$`2`< -1/2)
#f <- sort(unique(c(a,b,c)))
#prof_fd = prof[f,]
prof_fd=prof

corr <- corr_1(t(prof_fd))
corr_f <- cor_filter_diff(corr$r,corr$p,0,0.01)
#corr_fdat <- two_matr_rp(corr_f$r,corr_f$p)
corr_fr <- reshape_dat(corr_f$r)
corr_fp <- reshape_dat(corr_f$p)
corr_fr$pvalue <- corr_fp$value
corr_fr$lab <- ''
corr_fr$lab[corr_fr$pvalue<=0.01]<-'+'
#corr_fr <- corr_fr[corr_fr$lab == '+',]

#corr_fdat <- corr_fr[which(corr_fr$lab=='+'),]
corr_fdat <- corr_fr
corr_fdat$variable <- as.character(corr_fdat$variable)
corr_fdat$edge_size <- corr_fdat$value
corr_fdat$edge_size[corr_fdat$lab!='+']<-0
#corr_fdat$ID <- species(corr_fdat$ID)
#corr_fdat$variable <- species(corr_fdat$variable)
#write.table(corr_fdat,'analysis/5.correlation/genus_selfCorrelation_trajectories1.edge',sep='\t',row.names = F,quote = F)

node = data.frame(ID=unique(c(as.character(corr_fdat$ID),as.character(corr_fdat$variable))),stringsAsFactors = F)
for(x in group){
  p <- prof[node$ID,map[,3]==x]
  node[,as.character(x)]=rowMeans(p)
}
node$all <- rowMeans(prof[node$ID,])

#corr_fdat$edge_color <- abs(corr_fdat$value)
#corr_fdat$edge_color[which(corr_fdat$value<0)]<- -corr_fdat$edge_color[which(corr_fdat$value<0)]
#write.table(node,'analysis/5.correlation/genus_selfCorrelation_trajectories1.node',sep='\t',row.names = F,quote = F)

######################################
tax <- read.csv('analysis/5.correlation/genus_selfCorrelation_trajectories1附件/tax.mapping',check.names = F,sep='\t',stringsAsFactors = F)
node <- merge(node,tax,all.x = T)

#
GG <- 0:7
#GA2 <- data.frame(ID=corr_fdat$ID,variable=corr_fdat$variable,stringsAsFactors = F)
GA2 <- corr_fdat
for(i in GG){
  prof_fd=prof[,map[which(map$Gestation_age_G2==i),1]]
  
  corr <- corr_1(t(prof_fd))
  corr_f <- cor_filter_diff(corr$r,corr$p,0,1)
  #corr_fdat <- two_matr_rp(corr_f$r,corr_f$p)
  corr_fr <- reshape_dat(corr_f$r)
  corr_fp <- reshape_dat(corr_f$p)
  corr_fr$pvalue <- corr_fp$value
  corr_fr$lab <- ''
  corr_fr$lab[corr_fr$pvalue<=0.01]<-'+'
  #corr_fr <- corr_fr[corr_fr$lab == '+',]
  
  #corr_fdat <- corr_fr[which(corr_fr$lab=='+'),]
  corr_fdat <- corr_fr
  corr_fdat$variable <- as.character(corr_fdat$variable)
  corr_fdat$edge_size <- corr_fdat$value
  corr_fdat$edge_size[corr_fdat$lab!='+']<-0
  
  colnames(corr_fdat)[3:6] <- paste(colnames(corr_fdat)[3:6],as.character(i),sep='')
  GA2 = cbind.data.frame(GA2,corr_fdat[,3:6])
}

#write.table(GA2,'analysis/5.correlation/genus_selfCorrelation_trajectories1附件/genus_selfCorrelation_trajectories1.GA2.edge',sep='\t',row.names = F,quote = F)
#write.table(node,'analysis/5.correlation/genus_selfCorrelation_trajectories1附件/genus_selfCorrelation_trajectories1.GA2.node',sep='\t',row.names = F,quote = F)
#===================

V1 = 'value'
L1 = 'lab'
V2 = 'value7'
L2 = 'lab7'

length(which(GA2[,L1]=='+'))
length(which(GA2[,L2]=='+'))
length(which(GA2[,V1]< 0 & GA2[,L1]=='+' & GA2[,V2] < 0 & GA2[,L2]=='+')) + length(which(GA2[,V1]> 0 & GA2[,L1]=='+' & GA2[,V2] > 0 & GA2[,L2]=='+'))
#smooth#####################################################################
p <- prof[node$ID,]
#p <- prof
plo <- c()
for(x in 1:nrow(p)){
  #x=3
  pp <- data.frame(value=as.numeric(p[x,]),week=map[,4],ID=row.names(p)[x])
  pp <- pp[which(pp$week>50),]
 plo[[x]] <- 
    ggplot(pp,aes(x=week,y=value))+
      scale_x_continuous(breaks =seq(60,300,60))+
      #geom_point()+
      #geom_smooth(method = 'gam')+
      geom_smooth()+
      labs(title=row.names(p)[x],y='Relative Abundance',x='Gestational age')+
      theme_grey()+
     theme(
      plot.title = element_text(hjust=0.5,size = 10),
      panel.grid.minor = element_blank(),
      axis.text.x = element_blank(),
       axis.ticks.x = element_blank(),
       axis.title = element_blank()
   )
}
plo[[x+1]] <- 
  ggplot(pp,aes(x=week,y=c(rep(seq(0.001,0.01,0.001),150),0)[1:1477]))+
  scale_x_continuous(breaks =seq(60,300,60))+
  #geom_point()+
  #geom_smooth(method = 'gam')+
  geom_smooth()+
  labs(title='Key',y='Relative Abundance',x='Gestational week')+
  theme_grey()+
  theme(
      plot.title = element_text(hjust=0.5,size = 10),
      panel.grid.minor = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title = element_text(size=10)
  )

plot_grid(plotlist = plo)
#glm(value~week,data = pp)






