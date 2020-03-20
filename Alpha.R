setwd('E:/puyuan/Project/PMD/2018.PMD3/2018.PMD3/2018.Oct24.PMD3_qiime2_lish/')
library('ggplot2');
library('RColorBrewer');
library('ggsci');
library('ggpubr');
library(cowplot)
D <- read.table(file=".//analysis/1.alpha_div/faith_pd_vector.txt", sep='\t', head=T,stringsAsFactors = F);
colnames(D)[2] <- 'Alpha'
map1 <- read.csv('./data/mapping.file/mapping3.file',sep='\t',stringsAsFactors = F,check.names = F)
map2 <- read.csv('./data/mapping.file/mapping5.file',sep='\t',stringsAsFactors = F,check.names = F)
D <- merge(D,map1[,c(1,2)],by.x = 'X',by.y ='#SampleID',all.y = T)
D <- merge(D,map2[,c(1,2)],by.x = 'X',by.y ='#SampleID',all.y = T)

D$Gestation_age_G1<- as.character(D$Gestation_age_G1)
D$Gestation_age_G2<- as.character(D$Gestation_age_G2)
#pdf(".//analysis/1.alpha_div/boxplot/Gestation_age_G2.Gestation_age_G2.shannon.pdf", width=9, height=4);
#png(".//analysis/1.alpha_div/boxplot/Gestation_age_G2.Gestation_age_G2.shannon.png", width=9, height=4, units = "in", res=100);
#svg(".//analysis/1.alpha_div/boxplot/Gestation_age_G2.Gestation_age_G2.shannon.svg", width=9, height=4);
#mypal = pal_npg("nrc", alpha = 0.7)(9)
mypal = pal_lancet("lanonc", alpha = 0.7)(9)
mycolor = colorRampPalette(mypal)
#my_compare=list(c("",""),c("",""))

g<-c("1-12W","13-16W","17-20W","21-24W","25-28W","29-32W","33-36W","37-45W")
for(x in 1:length(g)){
  D$Gestation_age_G2[D$Gestation_age_G2 == as.character(x-1)]<-g[x]
}



D$Gestation_age_G2 <-factor(D$Gestation_age_G2,levels = sort(unique(D$Gestation_age_G2)))

p4<-
  ggplot(data = D, aes(x = Gestation_age_G2, y = Alpha)) +
#geom_violin(aes(fill=factor(Group)), size = 1) +
#geom_boxplot(size = 1, width=0.15, outlier.size = 0.4) +
#geom_jitter(aes(color=factor(Group)), width = 0.3, size=2, alpha = 0.6) +
geom_boxplot( aes(fill=Gestation_age_G1),size = 0.5, width=0.7, outlier.size = 0.4, alpha = 0.8) +
#stat_compare_means(aes(label=..p.format..), method='anova') +
  stat_compare_means(method='anova') +
#facet_grid(~Methods, scales = "free", space = "free_x") +
  theme_classic() +
  theme(
        #panel.background=element_blank(),
        #axis.line=element_line(
        #                       colour="black",
        #                       size = 1),
        #legend.box.background=element_blank(),
        #legend.box.background=element_blank(),
        axis.title=element_text(
                                family="Helvetica",
                                face="bold",
                                colour="black",
                                size=16
                                ),
        axis.text=element_text(
                               family="Helvetica",
                               colour="black",
                               size=14
                               ),
        axis.text.x=element_text(
                                 angle=90,
                                 hjust=1,
                                 vjust=0.5
                                 ),
        strip.text=element_text(
                                family="Helvetica",
                                colour="black",
                                size=14
                                )

        ) +
guides(color=FALSE) +
#guides(fill=FALSE) +
#scale_fill_discrete(h = c(0, 360) + 270, breaks = rev(as.vector(unique(D$Group)))) +
#scale_fill_brewer(palette='Paired', breaks = rev(as.vector(unique(D$Group)))) +
#scale_fill_manual(values = mycolor(9), breaks = rev(as.vector(unique(D$Group)))) +
#scale_color_manual(values = mycolor(9), breaks = rev(as.vector(unique(D$Group)))) +
#scale_fill_lancet() +
#scale_color_lancet() +
#labs(color="Group") +
xlab('') +
#ylab('Shannon Diversity');
#ylab('Observed otus');
#ylab('Evenness Diversity');
ylab('Phylogenetic Diversity ')
  #dev.off()

#plot_grid(plotlist = p)
plot_grid(p1,p2,p3,p4,nrow=2)
bartlett.test(D$Alpha ~ D$Group)
a<-aov(D$Alpha ~ D$Group)
TukeyHSD(a)
summary(a)
kruskal.test(D$Alpha ~ D$Group)
