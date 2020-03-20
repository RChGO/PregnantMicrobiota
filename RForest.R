library(argparse)

parser <- ArgumentParser()

parser$add_argument("-i", help = "input otufile or Kegg with the descriptors")
parser$add_argument("-p", help = "input file with the classified")
parser$add_argument("-c", help = "input column of map, default 2",default=2)
parser$add_argument("-n", help = "input row of sample")
parser$add_argument("-o", help = "input outfile name")

args <- parser$parse_args()

infile <- file.path(args$i)
inmap <- file.path(args$p)
inc <- as.numeric(file.path(args$c))
num <- as.numeric(file.path(args$n))
outfile <- file.path(args$o)


library(randomForest)
library(pROC)

wtest2<-function(tab,map,group,p=F,sample='a'){
  if (p==T){   #配对检验需要
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

#infile = 'data/profile/coreFeature.txt'
#inmap = 'data/mapping.file/mapping6.file'
#inc = 2
#num = 1
#outfile = 'analysis/6.classification/single/mapping6'

prof <- read.csv(infile,sep='\t',row.names = 1,check.names = F,stringsAsFactors = F)
map <- read.csv(inmap,sep = '\t',check.names = F,stringsAsFactors = F)
map = data.frame(ID=map[,1],Group=map[,inc],check.names = F,stringsAsFactors = F)
prof <- prof[,map[,1]]

tt <- wtest2(prof,map,'Group')
prof <- prof[tt$ID[tt$pvalue<0.01],]

prof <- as.data.frame(t(prof))
prof[,'Group'] <- as.factor(map$Group)
colnames(prof)<-make.names(colnames(prof))

train <- prof[-num,]
test <- prof[num,]

fit <- randomForest(Group~ .,data=train,ntree=2000,
               proximity=TRUE,importance=T) 

test_true <- test$Group
test <- subset(test,select = -c(Group))
test_pred <- predict(fit,type="prob",newdata = test)
#p_true <- as.numeric(factor(percent20_testdata$pclass))
test_pred <- as.numeric(test_pred[,2])
pred <- data.frame(ID=row.names(test),true=as.character(test_true),pred=test_pred,stringsAsFactors = F)
write.table(pred,paste(outfile,pred$ID,sep='.'),sep='\t',row.names = F,quote = F,col.names = FALSE)
#p_myroc <- roc(test_true,test_pred,percent=TRUE,print.auc=TRUE,auc.polygon=TRUE,auc.polygon.col=NA,plot = T)
#p_auc=as.character(lapply(as.numeric(ci.auc(p_myroc)), function(x){round(x,2)}))
#ci_se =ci.se(p_myroc)
#plot.roc(p_myroc)
#plot(ci_se,type = 'shape',col='grey90')
#plot(ci_se)
#text(20,20,paste('AUC: ',p_auc[2],'%','\n(95% CI: ',p_auc[1],'% ~ ',p_auc[3],'%)',sep = ''))