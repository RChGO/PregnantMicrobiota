library(ggrepel)
library(ggplot2)
library(RColorBrewer)
species<-function(ii){
  ii<-as.character(ii)
  ii[intersect(grep('.*[|,.;]s__.*[|,.;]t__..*',ii),grep('[|,.;]t__$',ii,invert=T))]<-as.character(lapply(ii[intersect(grep('.*[|,.;]s__.*[|,.;]t__..*',ii),grep('[|,.;]t__$',ii,invert=T))],function(x){gsub('.*[|,.;]t','t',x)}))
  ii[intersect(grep('.*[|,.;]g__.*[|,.;]s__..*',ii),grep('s__[|,.;]t__',ii,invert=T))]<-as.character(lapply(ii[intersect(grep('.*[|,.;]g__.*[|,.;]s__..*',ii),grep('s__[|,.;]t__',ii,invert=T))],function(x){gsub('.*[|,.;]s','s',x)}))
  ii[intersect(grep('.*[|,.;]f__.*[|,.;]g__..*',ii),grep('g__[|,.;]s__',ii,invert=T))]<-as.character(lapply(ii[intersect(grep('.*[|,.;]f__.*[|,.;]g__..*',ii),grep('g__[|,.;]s__',ii,invert=T))],function(x){gsub('.*[|,.;]g','g',x)}))
  ii[intersect(grep('.*[|,.;]o__.*[|,.;]f__..*',ii),grep('f__[|,.;]g__',ii,invert=T))]<-as.character(lapply(ii[intersect(grep('.*[|,.;]o__.*[|,.;]f__..*',ii),grep('f__[|,.;]g__',ii,invert=T))],function(x){gsub('.*[|,.;]f','f',x)}))
  ii[intersect(grep('.*[|,.;]c__.*[|,.;]o__..*',ii),grep('o__[|,.;]f__',ii,invert=T))]<-as.character(lapply(ii[intersect(grep('.*[|,.;]c__.*[|,.;]o__..*',ii),grep('o__[|,.;]f__',ii,invert=T))],function(x){gsub('.*[|,.;]o','o',x)}))
  ii[intersect(grep('.*[|,.;]p__.*[|,.;]c__..*',ii),grep('c__[|,.;]o__',ii,invert=T))]<-as.character(lapply(ii[intersect(grep('.*[|,.;]p__.*[|,.;]c__..*',ii),grep('c__[|,.;]o__',ii,invert=T))],function(x){gsub('.*[|,.;]c','c',x)}))
  ii[intersect(grep('k__.*[|,.;]p__..*',ii),grep('k__[|,.;]p__',ii,invert=T))]<-as.character(lapply(ii[intersect(grep('k__.*[|,.;]p__..*',ii),grep('k__[|,.;]p__',ii,invert=T))],function(x){gsub('.*[|,.;]p','p',x)}))
  return(ii)
}

infile = "E:/puyuan/Project/PMD/2018.PMD3/2018.PMD3/2018.Oct24.PMD3_qiime2_lish/analysis/2.beta_div/plots/pcoa/mapping3.Gestation_age_G2.enterotype_pcoa.sample.dat"
ine='E:/puyuan/Project/PMD/2018.PMD3/2018.PMD3/2018.Oct24.PMD3_qiime2_lish/analysis/2.beta_div/plots/pcoa/mapping3.Gestation_age_G2.enterotype_pcoa.eig.dat'
inp = 'E:/puyuan/Project/PMD/2018.PMD3/2018.PMD3/2018.Oct24.PMD3_qiime2_lish/analysis/2.beta_div/plots/pcoa/mapping3.Gestation_age_G2.enterotype_pcoa.importan_bacteria.dat'

samp=read.csv(infile, header=T, check.names = F, sep="\t",row.names=1,stringsAsFactors = F)
eig = read.csv(ine, header=T, check.names = F, sep="\t",stringsAsFactors = F)
pc = read.csv(inp, header=T, check.names = F, sep="\t",row.names=1,stringsAsFactors = F)
pc$ID = species(row.names(pc)) 
pc <- pc[1:8,]
cutoff = 0.4
pc$A1<-pc$A1*cutoff ; pc$A2=pc$A2*cutoff
pc$ID<- gsub('^g__','',pc$ID)
p<-ggplot(samp) +
  xlab("") +
  ylab("") +
  #geom_vline(xintercept = 0) +
  #geom_hline(yintercept = 0) +
  #p+geom_point(aes(x=d[,1], y=d[,2], color=Type), size=2, shape=19) +
  geom_point(aes(x=samp[,1], y=samp[,2], color=as.character(samp[,3])), size=0.8) +
  #stat_ellipse(aes(x=samp[,1], y=samp[,2], color=as.character(samp[,3])),fill=NA,size=1, geom="polygon", level=0.8, alpha=0.3) +
  geom_text_repel(data=pc, aes(x=A1, y=A2, label=ID),
                  size=4, check_overlap=TRUE,fontface="italic") +
  geom_segment(data=pc, aes(x=0, y=0, xend=A1, yend=A2),
               arrow = arrow(length = unit(0.3, "cm")), size=0.5, alpha=0.5)+
  #geom_segment(data=track, aes(x=x1, y=y1, xend=x2, yend=y2),
  #             arrow = arrow(length = unit(0.4, "cm")), size=1, alpha=0.8) +
  #geom_point(data=points, aes(x=x1, y=y1, color=factor(Type)),
  #           size=6, shape=19, alpha=0.7) +
  #geom_text(data=points, aes(x=x1, y=y1, label=Class),
  #          family="Helvetica", size=6, check_overlap=TRUE) +
  scale_color_manual(values=brewer.pal(9,"Set2"))+
  #  scale_fill_manual(values=brewer.pal(9,"Set1"))+
  #scale_shape_manual(values=1:20)+
  guides(color=F
  #       fill=guide_legend(colnames(data1)[0+1]),
  #       shape=guide_legend(colnames(data1)[0+1]) 
  ) +
  #theme_classic()+
  theme_minimal()+
  theme(
        #panel.grid = element_blank(),
        axis.text = element_blank(),
        #axis.ticks = element_blank(),
        #axis.line = element_blank(),
        axis.title = element_text(
          size=13
        ),
        legend.position=c(0.9,0.7),
        text=element_text(
          family="Helvetica",
          face="bold",
          colour="black",
          size=9
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
  )+
  xlab(paste("PCoA",1," (",round(eig$eig[1]*100,0),"%)",sep=""))+
  ylab(paste("PCoA",2," (",round(eig$eig[2]*100,0),"%)",sep=""))
#xlim(min(d[,1], 0)*xx, max(d[,1])*xx)+ylim(min(d[,2], 0)*yy, max(d[,2])*yy)
#xlim(min(d[,1], 0)*3, max(d[,1])*3)+ylim(min(d[,2], 0)*3, max(d[,2])*3)

p



plot(x=samp[,1],y=samp[,2])
z=kde2d(x=samp[,1],y=samp[,2])
contour(z, col = 'red', drawlabel=FALSE, main="Density estimation :cont Plot")

persp(z,col = "lightblue",
      theta = 0, 
      phi = 70,
      #theta = 160, 
      #phi = 40,
      expand = 0.5,
      ltheta = 120, 
      shade = 0.75,
      xlab=paste("PCoA",1," (",round(eig$eig[1]*100,0),"%)",sep=""),
        ylab=paste("PCoA",2," (",round(eig$eig[2]*100,0),"%)",sep=""),
      zlab = "Density"
)
