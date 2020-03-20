#!/usr/bin/env Rscript

# -------------- NOTE: ---------------
# The last column should be the class.
# ------------------------------------
library(argparse)
library(parallel)

parser <- ArgumentParser()
parser$add_argument("-i", help = "input profile")
parser$add_argument("-p", help = "input index.prof")
parser$add_argument("-m", help = "input mapping of index.prof")
parser$add_argument("-c", help = "input  group and variable type of mapping, e.g. Group,Type")
#parser$add_argument("-f", help = "input group of medicine mapping")
parser$add_argument("-o", help = "input outfile name")
args <- parser$parse_args()
inf <- file.path(args$i)
inp <- file.path(args$p)
inpp <- file.path(args$m)
inc <-  file.path(args$c)
ino <- file.path(args$o)

corr_1 = function(matr){
  cor_r = matrix(rep(0,ncol(matr)*ncol(matr)),ncol(matr),ncol(matr))
  row.names(cor_r) = colnames(matr)
  colnames(cor_r) = colnames(matr)
  cor_p = matrix(rep(0,ncol(matr)*ncol(matr)),ncol(matr),ncol(matr))
  row.names(cor_p) = colnames(matr)
  colnames(cor_p) = colnames(matr)
  for(i in 1:ncol(matr)){
    for(j in 1:ncol(matr)){
      cor_1 = cor.test(matr[,i],matr[,j],method = "pearson")
      cor_r[i,j] = cor_1$estimate
      cor_p[i,j] = cor_1$p.value
    }
  }
  k = list(p = cor_p,r = cor_r)
  return(k)
}
cor_filter_diff <- function(r,p,rlim=0.3,plim=0.01){
  #r <- a$`Urine test`$r
  #p <- a$`Urine test`$p
  r_f <- r
  r_f[p>plim]<-NA
  r_f[abs(r_f)<rlim]<-NA
  fr <- apply(r_f,1,function(x){length(which(is.na(x)==F))})
  fc <- apply(r_f,2,function(x){length(which(is.na(x)==F))})
  r <- r[fr!=1,fc!=1]
  p <- p[fr!=1,fc!=1]
  #n <- c()
  #for(x in 1:nrow(r_f)){
  #  for(y in x:nrow(r_f)){
  #   if(x!=y & !is.na(r_f[x,y]) &abs(r_f[x,y])>=rlim){
  #     n <- c(n,colnames(r_f)[y])
  #   } 
  #  }
  #}
  #f = setdiff(colnames(r_f),n)
  n<-c()
  for(x in 1:nrow(r_f)){
    for(y in x:nrow(r_f)){
      if(x != y & identical(is.na(r_f[x,]),is.na(r_f[y,]))){
        n <-c(n,y)
      }
    }
  }
  f = setdiff(1:nrow(r_f),unique(n))
  f = row.names(r_f)[f]
  return(list(r=r,p=p,f=f))
}
Dereduntance<-function(m,p,incc){
  #m =map ;p = phen ; incc = 'Group'
  p <- p[,m[,1]]
  g <- unique(m[,incc])
  a<- c() ; b<-c() ; fn<-c()
  for(x in 1:length(g)){
    mf <- m [m[,incc]==g[x],]
    pf <- p[,mf[,1]]
    a[[g[x]]]<- corr_1(pf)
    b[[g[x]]]<- cor_filter_diff(a[[g[x]]]$r,a[[g[x]]]$p)
    fn <- c(fn,b[[g[x]]]$f)
  }
  return(m[m[,1]%in%fn,])
}
library(ggplot2)
library(vegan)
#inf = '/share/data2/guorc/Project/PMD/2018.PMD3/2018.Sep21.PMD3_qiime2_lish/data/medicine/medicine1.prof'
#inp = '/share/data2/guorc/Project/PMD/2018.PMD3/2018.Sep21.PMD3_qiime2_lish/data/medicine/medicine.mapping'
#inpp ='/share/data2/guorc/Project/PMD/2018.PMD3/2018.Sep21.PMD3_qiime2_lish/data/profile/otus.prof'


inp = 'D:/puyuan/Project/PMD/2018.PMD3/2018.PMD3/2018.Oct24.PMD3_qiime2_lish/data/medicine/medicine7.prof'
inpp = 'D:/puyuan/Project/PMD/2018.PMD3/2018.PMD3/2018.Oct24.PMD3_qiime2_lish/data/medicine/medicine.mapping_lish'
inf ='D:/puyuan/Project/PMD/2018.PMD3/2018.PMD3/2018.Oct24.PMD3_qiime2_lish/data/profile/kegg.ko.prof'
inc = 'Variable group,Type'

phen = read.csv(inp,sep='\t',row.names = 1,check.names = F,stringsAsFactors = F)
map <- read.csv(inpp,sep='\t',check.names = F,stringsAsFactors = F)
prof = read.csv(inf,sep='\t',row.names = 1,check.names = F,stringsAsFactors = F)
prof = t(prof)
inc = strsplit(inc,',')[[1]]
f = intersect(colnames(phen),map[,1])
phen = phen[,f]
map  = map[map[,1]%in%f,]

for(x in 1:nrow(map)){
  if(map[x,inc[2]]=='Discrete'){
    phen[,map[x,1]] = as.character(phen[,map[x,1]])
  }else{
    phen[,map[x,1]] = as.numeric(phen[,map[x,1]])  
  }
}


#ff=as.numeric(strsplit(inff,',')[[1]])
g <- unique(map[,inc[1]])[3]
#g <- g[ff]
#map <- map[map$Group%in%g,]
#map <- map[!map$Group%in%c('Immune Parameters','Liver Function','Renal Function'),]

#map <- Dereduntance(map,phen,'Group')

#for (x in 1:ncol(phen)){
#  if(class(phen[,x])=='character'){
#    print(colnames(phen)[x])
#    }
#}

a <- c() ; sample_num <-c();bianliang_num<-c()
tab <- data.frame(medicine=g,r2=NA,p=NA,sample_num=NA,variable_num=NA,R2_adj=NA,stringsAsFactors = F)

for (x in 1:length(g)){
  # x=1  
  m <- map[map[,inc[1]]==g[x],]
  ph <- as.data.frame(phen[,m[,1]])
  row.names(ph)=row.names(phen) ;colnames(ph)= m[,1]
  
  n <- apply(ph,1,function(x){length(which(is.na(x)==T))})
  #mm <- apply(ph,2,function(x){length(which(is.na(x)==T))})
  #print(paste(g[x],length(which(n>0))))
 
  ph <- as.data.frame(ph[n==0,])
  row.names(ph)=row.names(phen)[n==0] ;colnames(ph)= m[,1]
  
 # for(i in 1:ncol(ph)){
 #   if(m$Type[i]=='Discrete'){
 #     ph[,i] = as.character(ph[,i])
 #   }else{
 #     ph[,i] = as.numeric(ph[,i])
 #   }
 # }
  
  
  n <- apply(ph,2,function(x){length(unique(x))})
  ph <- ph[,which(n>1)]
  p <- prof[row.names(ph),]
  tab[x,4] <- nrow(ph)
  tab[x,5] <- ncol(ph)
 
  a[[g[x]]] <- adonis(p~. ,data = ph)
  tab[x,2] <- 1-a[[g[x]]]$aov.tab$R2[length(a[[g[x]]]$aov.tab$R2)-1]
  print(tab[x,])
}
tab[nrow(tab)+1,]<-NA
tab[nrow(tab),1]<-'all'

pp <- phen[,map[,1]]
n <- apply(pp,1,function(x){length(which(is.na(x)==T))})
pp <- pp[n==0,]
pro <- prof[row.names(pp),]
tab[nrow(tab),4]<- nrow(pp) ; tab[nrow(tab),5]<-ncol(pp)
a <- adonis(pro~. ,data = pp)
tab$r2[nrow(tab)] = 1-a$aov.tab$R2[length(a$aov.tab$R2)-1]

write.table(tab,ino,
            row.names = F,sep='\t',quote=F)

q()

library(vegan)
tab<-read.csv('E:/puyuan/Project/PMD/2018.PMD3/2018.PMD3/2018.Oct24.PMD3_qiime2_lish/analysis/2.beta_div/compare/adonis/adonis_combine.feature',sep = '\t',stringsAsFactors = F)
for(x in 1:nrow(tab)){
  tab$R2_adj[x] <- RsquareAdj(tab$r2[x],tab$sample_num[x],tab$variable_num[x])
}
write.table(tab,'E:/puyuan/Project/PMD/2018.PMD3/2018.PMD3/2018.Oct24.PMD3_qiime2_lish/analysis/2.beta_div/compare/adonis/adonis_combine.feature2',
            sep='\t',row.names = F,quote = F)

tab<-tab[tab$medicine!='all',]
tab <- tab[order(tab$r2),]
tab$medicine <- factor(tab$medicine,levels = tab$medicine)
ggplot(data = tab,
       aes(x=medicine,
           y=r2,fill=medicine
           )) +
  geom_bar(position='dodge',stat='identity',color='grey20',alpha=0.8,width = 0.8)+
  coord_flip()+                      
  labs(title="", x="", y="R2")+
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title=element_text(
      family="Helvetica",
      face="bold",
      colour="black",
      size=12
    ),
    axis.text=element_text(
      family="Helvetica",
      colour="black",
      size=10
    ),
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
    legend.title=element_text(
      family="Helvetica",
      colour="black",
      size=18
    ),
    legend.text=element_text(
      family="Helvetica",
      colour="black",
      size=16
    )
  ) +
  guides(fill = F) + 
  scale_color_manual(values = c('#4da968','#604c90','#7b92bd','#6cbbc9','#00572f'))+
  scale_y_continuous(position = 'top')#,breaks = c(0.001,0.005,0.01,0.10,0.15,0.20))
