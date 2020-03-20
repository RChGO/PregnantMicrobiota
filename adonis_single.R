#!/usr/bin/env Rscript

# -------------- NOTE: ---------------
# The last column should be the class.
# ------------------------------------
library(argparse)

parser <- ArgumentParser()
parser$add_argument("-i", help = "input profile")
parser$add_argument("-p", help = "input mapping file")
parser$add_argument("-c", help = "input column name of group, default all")
parser$add_argument("-t", help = "input type of group, default Discrete",choices=c('Discrete','Continuous'),default='Dscrete')
parser$add_argument("-o", help = "input outfile name")
args <- parser$parse_args()
infile <- file.path(args$i)
inm <- file.path(args$p)
inc <- file.path(args$c)
ino <- file.path(args$o)
intype <- file.path(args$t)

library(vegan)
#infile <- 'E:/puyuan/Project/PMD/2018.PMD3/2018.PMD3/2018.Oct24.PMD3_qiime2_lish/data/profile/Feature.txt'
#inm <- 'E:/puyuan/Project/PMD/2018.PMD3/2018.PMD3/2018.Oct24.PMD3_qiime2_lish/data/medicine/medicine6.prof'
#inc <- 'Gestation_age_G1'
#intype <- 'Discrete'
#====================================
prof <- read.csv(infile,sep='\t',row.names = 1,check.names = F,stringsAsFactors = F)
map <- read.csv(inm,sep='\t',check.names = F,stringsAsFactors = F)
prof <- t(prof[,map[,1]])
f = as.numeric(apply(map,1,function(x){length(which(is.na(x)==T))}))
map <- map[f==0,]


if(length(inc)==0){
  inc <- colnames(map)[-1]
}
map <- data.frame(ID=map[,1],map[,inc],stringsAsFactors = F)
colnames(map)[-1]<-inc

if(intype=='Discrete'){
  for(x in 2:ncol(map)){
    map[,x]<- as.character(map[,x])
  }
}else{
  for(x in 2:ncol(map)){
    map[,x]<- as.numeric(map[,x])
  }
}
tab <- data.frame(ID = colnames(map)[-1],R2=NA,Pvalue=NA,Sample=NA,stringsAsFactors = F)
for(x in 2:ncol(map)){
  m <- data.frame(ID=map[,1],Group=map[,x],stringsAsFactors = F)
  m <- m[!is.na(m$Group),]
  p <- prof[m[,1],]
  ado <- adonis(p~Group,data=m)
  tab[x-1,'Pvalue'] <- ado$aov.tab$`Pr(>F)`[1]
  tab[x-1,'R2'] <- ado$aov.tab$R2[1]
  tab[x-1,'Sample'] <- nrow(m)
}

#print(tab)
write.table(tab,ino,sep='\t',quote = F,row.names = F)
q()



#plot===============================

tab <- read.csv('analysis/2.beta_div/compare/adonis/adonis_split_otus.all',sep = '\t',check.names = F,stringsAsFactors = F)
tab$lab <- ''
tab$lab[tab$Pvalue<0.05]<-'*'
#tab$lab <- paste(tab$medicine,tab$lab)
tab <- tab[order(tab$R2),]
#tab$lab <- factor(tab$lab,levels = tab$lab) 
tab$ID <- factor(tab$ID,levels = tab$ID) 
library(ggplot2)
ggplot(data = tab,
       aes(x=ID,
           y=R2,fill=lab)) +
  geom_bar(position='dodge',stat='identity',color='grey20',width = 0.8,size=0.1)+
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
    )) +
  guides(fill = F) + 
  scale_fill_manual(values = c('#00bfc4','#f8766d'))+
  scale_y_continuous(position = 'top')#,breaks = c(0.001,0.005,0.01,0.10,0.15,0.20))
