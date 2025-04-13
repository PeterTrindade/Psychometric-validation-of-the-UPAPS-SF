#### (0) Packages ####
rm(list=ls())

pacotes=c(
  "dplyr","ggplot2", 'biostatUZH','irr', 'devtools', 'factoextra', 
  "foreign", "lavaan", "GPArotation", "ROCR", "FactoMineR", "burt",
  "tidySEM", "reshape2","knitr","kableExtra", "pROC", "ggpubr",
  "nlme","lmtest","fastDummies","msm","lmeInfo","jtools", "lme4", 
  'lsmeans', 'multcompView', 'multcomp', "epiR", "DescTools", "patchwork",
  "GDAtools"
)
if(sum(as.numeric(!pacotes %in% installed.packages()))!= 0){
  instalador=pacotes[!pacotes %in%installed.packages()]
  for(i in 1:length(instalador)){
    install.packages(instalador,dependencies=T)
    break()}
  sapply(pacotes,require,character=T) 
} else {
  sapply(pacotes, require,character=T) 
}

rm(list=ls())

#### (1) Reliability ####
rm(list=ls())
#----Importing data-------
df = read.csv2("UPAPS Short data.csv")
set.seed(360)
#----Pre processing-----

#1. Filter subsets by observer
#EX: df_observer2 = subset(df_sheep_subset, Observer == "2")


table(df$Observer)
df$Order=paste(df$Video.YT,df$Litter, df$Observer, df$Study)
df=df[order(df$Order),]

df_observer1 = subset(df, Observer == "1" & Study == "training")
df_observer2 = subset(df, Observer == "2"& Study == "training")
df_observer3 = subset(df, Observer == "3"& Study == "training")
df_observer4 = subset(df, Observer == "4"& Study == "training")
df_observer5 = subset(df, Observer == "5"& Study == "training")
df_observer6 = subset(df, Observer == "6"& Study == "training")

#-----Intra-observer reliability-----


#Looping

#----Observer 1------

table(df_observer1$Head.down[df_observer1$Phase == "1"],
      df_observer1$Head.down[df_observer1$Phase == "2"])

table(df_observer1$Interaction[df_observer1$Phase == "1"],
      df_observer1$Interaction[df_observer1$Phase == "2"])

table(df_observer1$Activity[df_observer1$Phase == "1"],
      df_observer1$Activity[df_observer1$Phase == "2"])

table(df_observer1$Sits.with.difficulty[df_observer1$Phase == "1"],
      df_observer1$Sits.with.difficulty[df_observer1$Phase == "2"])

table(df_observer1$Wags.tail[df_observer1$Phase == "1"],
      df_observer1$Wags.tail[df_observer1$Phase == "2"])



df_selec1=df_observer1[,c('Total',
                          #'Head.down',
                          #'Interaction',
                          #'Activity',
                          'Sits.with.difficulty',
                          'Wags.tail'
                          
)]

?confIntKappa
KappaEstimated=lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df_observer1$Phase=="1"],
    x[df_observer1$Phase=="2"]),nrow=20,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$kappa,dig=2))

KappaLower = lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df_observer1$Phase=="1"],
    x[df_observer1$Phase=="2"]),nrow=20,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$boot.quant[[1]], dig=2))

KappaUpper=lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df_observer1$Phase=="1"],
    x[df_observer1$Phase=="2"]),nrow=20,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$boot.quant[[2]], dig=2))

ICC=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df_observer1$Phase=="1"],
    x[df_observer1$Phase=="2"]),nrow=20,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$value,dig=2))

ICC_CI_Lower=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df_observer1$Phase=="1"],
    x[df_observer1$Phase=="2"]),nrow=20,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$lbound,dig=2))

ICC_CI_Upper=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df_observer1$Phase=="1"],
    x[df_observer1$Phase=="2"]),nrow=20,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$ubound,dig=2))

ICC_Pvalue=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df_observer1$Phase=="1"],
    x[df_observer1$Phase=="2"]),nrow=20,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$p.value,dig=40))


col0=c(colnames(df_selec1))
row0=c("KappaEstimated",
       "KappaLower",
       "KappaUpper",
       "ICC",
       "ICC_CI_Lower",
       "ICC_CI_Upper",
       "ICC_Pvalue")

table=matrix(cbind(KappaEstimated,
                   KappaLower,
                   KappaUpper,
                   ICC,
                   ICC_CI_Lower,
                   ICC_CI_Upper,
                   ICC_Pvalue),
             nrow=length(col0),ncol=length(row0),dimnames=list(col0,row0))

write.csv2(table,"(1) Intra-observer reliability - observer 1.csv")



#----Observer 2------
table(df_observer2$Head.down[df_observer2$Phase == "1"],
      df_observer2$Head.down[df_observer2$Phase == "2"])

table(df_observer2$Interaction[df_observer2$Phase == "1"],
      df_observer2$Interaction[df_observer2$Phase == "2"])

table(df_observer2$Activity[df_observer2$Phase == "1"],
      df_observer2$Activity[df_observer2$Phase == "2"])

table(df_observer2$Sits.with.difficulty[df_observer2$Phase == "1"],
      df_observer2$Sits.with.difficulty[df_observer2$Phase == "2"])

table(df_observer2$Wags.tail[df_observer2$Phase == "1"],
      df_observer2$Wags.tail[df_observer2$Phase == "2"])



df_selec1=df_observer2[,c('Total',
                          'Head.down',
                          #'Interaction',
                          'Activity',
                          'Sits.with.difficulty',
                          'Wags.tail'
                          
)]

KappaEstimated=lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df_observer2$Phase=="1"],
    x[df_observer2$Phase=="2"]),nrow=20,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$kappa,dig=2))

KappaLower = lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df_observer2$Phase=="1"],
    x[df_observer2$Phase=="2"]),nrow=20,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$boot.quant[[1]], dig=2))

KappaUpper=lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df_observer2$Phase=="1"],
    x[df_observer2$Phase=="2"]),nrow=20,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$boot.quant[[2]], dig=2))

ICC=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df_observer2$Phase=="1"],
    x[df_observer2$Phase=="2"]),nrow=20,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$value,dig=2))

ICC_CI_Lower=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df_observer2$Phase=="1"],
    x[df_observer2$Phase=="2"]),nrow=20,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$lbound,dig=2))

ICC_CI_Upper=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df_observer2$Phase=="1"],
    x[df_observer2$Phase=="2"]),nrow=20,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$ubound,dig=2))

ICC_Pvalue=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df_observer2$Phase=="1"],
    x[df_observer2$Phase=="2"]),nrow=20,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$p.value,dig=40))


col0=c(colnames(df_selec1))
row0=c("KappaEstimated",
       "KappaLower",
       "KappaUpper",
       "ICC",
       "ICC_CI_Lower",
       "ICC_CI_Upper",
       "ICC_Pvalue")

table=matrix(cbind(KappaEstimated,
                   KappaLower,
                   KappaUpper,
                   ICC,
                   ICC_CI_Lower,
                   ICC_CI_Upper,
                   ICC_Pvalue),
             nrow=length(col0),ncol=length(row0),dimnames=list(col0,row0))

write.csv2(table,"(1) Intra-observer reliability - observer 2.csv")



#----Observer 3------
table(df_observer3$Head.down[df_observer3$Phase == "1"],
      df_observer3$Head.down[df_observer3$Phase == "2"])

table(df_observer3$Interaction[df_observer3$Phase == "1"],
      df_observer3$Interaction[df_observer3$Phase == "2"])

table(df_observer3$Activity[df_observer3$Phase == "1"],
      df_observer3$Activity[df_observer3$Phase == "2"])

table(df_observer3$Sits.with.difficulty[df_observer3$Phase == "1"],
      df_observer3$Sits.with.difficulty[df_observer3$Phase == "2"])

table(df_observer3$Wags.tail[df_observer3$Phase == "1"],
      df_observer3$Wags.tail[df_observer3$Phase == "2"])





df_selec1=df_observer3[,c('Total',
                          'Head.down',
                          #'Interaction',
                          #'Activity',
                          #'Sits.with.difficulty',
                          'Wags.tail'
                          
)]

KappaEstimated=lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df_observer3$Phase=="1"],
    x[df_observer3$Phase=="2"]),nrow=20,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$kappa,dig=2))

KappaLower = lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df_observer3$Phase=="1"],
    x[df_observer3$Phase=="2"]),nrow=20,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$boot.quant[[1]], dig=2))

KappaUpper=lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df_observer3$Phase=="1"],
    x[df_observer3$Phase=="2"]),nrow=20,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$boot.quant[[2]], dig=2))

ICC=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df_observer3$Phase=="1"],
    x[df_observer3$Phase=="2"]),nrow=20,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$value,dig=2))

ICC_CI_Lower=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df_observer3$Phase=="1"],
    x[df_observer3$Phase=="2"]),nrow=20,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$lbound,dig=2))

ICC_CI_Upper=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df_observer3$Phase=="1"],
    x[df_observer3$Phase=="2"]),nrow=20,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$ubound,dig=2))

ICC_Pvalue=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df_observer3$Phase=="1"],
    x[df_observer3$Phase=="2"]),nrow=20,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$p.value,dig=40))


col0=c(colnames(df_selec1))
row0=c("KappaEstimated",
       "KappaLower",
       "KappaUpper",
       "ICC",
       "ICC_CI_Lower",
       "ICC_CI_Upper",
       "ICC_Pvalue")

table=matrix(cbind(KappaEstimated,
                   KappaLower,
                   KappaUpper,
                   ICC,
                   ICC_CI_Lower,
                   ICC_CI_Upper,
                   ICC_Pvalue),
             nrow=length(col0),ncol=length(row0),dimnames=list(col0,row0))

write.csv2(table,"(1) Intra-observer reliability - observer 3.csv")





#----Observer 4------
table(df_observer4$Head.down[df_observer4$Phase == "1"],
      df_observer4$Head.down[df_observer4$Phase == "2"])

table(df_observer4$Interaction[df_observer4$Phase == "1"],
      df_observer4$Interaction[df_observer4$Phase == "2"])

table(df_observer4$Activity[df_observer4$Phase == "1"],
      df_observer4$Activity[df_observer4$Phase == "2"])


table(df_observer4$Sits.with.difficulty[df_observer4$Phase == "1"],
      df_observer4$Sits.with.difficulty[df_observer4$Phase == "2"])

table(df_observer4$Wags.tail[df_observer4$Phase == "1"],
      df_observer4$Wags.tail[df_observer4$Phase == "2"])




df_selec1=df_observer4[,c('Total',
                          #'Head.down',
                          #'Interaction',
                          #'Activity',
                          'Sits.with.difficulty',
                          'Wags.tail'
                          
)]

KappaEstimated=lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df_observer4$Phase=="1"],
    x[df_observer4$Phase=="2"]),nrow=20,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$kappa,dig=2))

KappaLower = lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df_observer4$Phase=="1"],
    x[df_observer4$Phase=="2"]),nrow=20,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$boot.quant[[1]], dig=2))

KappaUpper=lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df_observer4$Phase=="1"],
    x[df_observer4$Phase=="2"]),nrow=20,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$boot.quant[[2]], dig=2))

ICC=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df_observer4$Phase=="1"],
    x[df_observer4$Phase=="2"]),nrow=20,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$value,dig=2))

ICC_CI_Lower=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df_observer4$Phase=="1"],
    x[df_observer4$Phase=="2"]),nrow=20,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$lbound,dig=2))

ICC_CI_Upper=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df_observer4$Phase=="1"],
    x[df_observer4$Phase=="2"]),nrow=20,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$ubound,dig=2))

ICC_Pvalue=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df_observer4$Phase=="1"],
    x[df_observer4$Phase=="2"]),nrow=20,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$p.value,dig=40))


col0=c(colnames(df_selec1))
row0=c("KappaEstimated",
       "KappaLower",
       "KappaUpper",
       "ICC",
       "ICC_CI_Lower",
       "ICC_CI_Upper",
       "ICC_Pvalue")

table=matrix(cbind(KappaEstimated,
                   KappaLower,
                   KappaUpper,
                   ICC,
                   ICC_CI_Lower,
                   ICC_CI_Upper,
                   ICC_Pvalue),
             nrow=length(col0),ncol=length(row0),dimnames=list(col0,row0))

write.csv2(table,"(1) Intra-observer reliability - observer 4.csv")



#----Observer 5------
table(df_observer5$Head.down[df_observer5$Phase == "1"],
      df_observer5$Head.down[df_observer5$Phase == "2"])

table(df_observer5$Interaction[df_observer5$Phase == "1"],
      df_observer5$Interaction[df_observer5$Phase == "2"])

table(df_observer5$Activity[df_observer5$Phase == "1"],
      df_observer5$Activity[df_observer5$Phase == "2"])

table(df_observer5$Sits.with.difficulty[df_observer5$Phase == "1"],
      df_observer5$Sits.with.difficulty[df_observer5$Phase == "2"])

table(df_observer5$Wags.tail[df_observer5$Phase == "1"],
      df_observer5$Wags.tail[df_observer5$Phase == "2"])



df_selec1=df_observer5[,c('Total',
                          #'Head.down',
                          'Interaction',
                          'Activity',
                          #'Sits.with.difficulty',
                          'Wags.tail'
                          
)]

KappaEstimated=lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df_observer5$Phase=="1"],
    x[df_observer5$Phase=="2"]),nrow=20,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$kappa,dig=2))

KappaLower = lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df_observer5$Phase=="1"],
    x[df_observer5$Phase=="2"]),nrow=20,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$boot.quant[[1]], dig=2))

KappaUpper=lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df_observer5$Phase=="1"],
    x[df_observer5$Phase=="2"]),nrow=20,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$boot.quant[[2]], dig=2))

ICC=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df_observer5$Phase=="1"],
    x[df_observer5$Phase=="2"]),nrow=20,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$value,dig=2))

ICC_CI_Lower=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df_observer5$Phase=="1"],
    x[df_observer5$Phase=="2"]),nrow=20,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$lbound,dig=2))

ICC_CI_Upper=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df_observer5$Phase=="1"],
    x[df_observer5$Phase=="2"]),nrow=20,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$ubound,dig=2))

ICC_Pvalue=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df_observer5$Phase=="1"],
    x[df_observer5$Phase=="2"]),nrow=20,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$p.value,dig=40))


col0=c(colnames(df_selec1))
row0=c("KappaEstimated",
       "KappaLower",
       "KappaUpper",
       "ICC",
       "ICC_CI_Lower",
       "ICC_CI_Upper",
       "ICC_Pvalue")

table=matrix(cbind(KappaEstimated,
                   KappaLower,
                   KappaUpper,
                   ICC,
                   ICC_CI_Lower,
                   ICC_CI_Upper,
                   ICC_Pvalue),
             nrow=length(col0),ncol=length(row0),dimnames=list(col0,row0))

write.csv2(table,"(1) Intra-observer reliability - observer 5.csv")




#----Observer 6------
table(df_observer6$Head.down[df_observer6$Phase == "1"],
      df_observer6$Head.down[df_observer6$Phase == "2"])

table(df_observer6$Interaction[df_observer6$Phase == "1"],
      df_observer6$Interaction[df_observer6$Phase == "2"])

table(df_observer6$Activity[df_observer6$Phase == "1"],
      df_observer6$Activity[df_observer6$Phase == "2"])

table(df_observer6$Sits.with.difficulty[df_observer6$Phase == "1"],
      df_observer6$Sits.with.difficulty[df_observer6$Phase == "2"])

table(df_observer6$Wags.tail[df_observer6$Phase == "1"],
      df_observer6$Wags.tail[df_observer6$Phase == "2"])

df_selec1=df_observer6[,c('Total',
                          'Head.down',
                          #'Interaction',
                          'Activity',
                          #'Sits.with.difficulty',
                          'Wags.tail'
                          
)]

KappaEstimated=lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df_observer6$Phase=="1"],
    x[df_observer6$Phase=="2"]),nrow=20,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$kappa,dig=2))

KappaLower = lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df_observer6$Phase=="1"],
    x[df_observer6$Phase=="2"]),nrow=20,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$boot.quant[[1]], dig=2))

KappaUpper=lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df_observer6$Phase=="1"],
    x[df_observer6$Phase=="2"]),nrow=20,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$boot.quant[[2]], dig=2))

ICC=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df_observer6$Phase=="1"],
    x[df_observer6$Phase=="2"]),nrow=20,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$value,dig=2))

ICC_CI_Lower=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df_observer6$Phase=="1"],
    x[df_observer6$Phase=="2"]),nrow=20,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$lbound,dig=2))

ICC_CI_Upper=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df_observer6$Phase=="1"],
    x[df_observer6$Phase=="2"]),nrow=20,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$ubound,dig=2))

ICC_Pvalue=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df_observer6$Phase=="1"],
    x[df_observer6$Phase=="2"]),nrow=20,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$p.value,dig=40))


col0=c(colnames(df_selec1))
row0=c("KappaEstimated",
       "KappaLower",
       "KappaUpper",
       "ICC",
       "ICC_CI_Lower",
       "ICC_CI_Upper",
       "ICC_Pvalue")

table=matrix(cbind(KappaEstimated,
                   KappaLower,
                   KappaUpper,
                   ICC,
                   ICC_CI_Lower,
                   ICC_CI_Upper,
                   ICC_Pvalue),
             nrow=length(col0),ncol=length(row0),dimnames=list(col0,row0))

write.csv2(table,"(1) Intra-observer reliability - observer 6.csv")



#----Figure----

df_figure = read.csv2("(1) Reliability figure.csv", dec=".")

tiff("(1) Figure 1 version 2. Reliability.tiff",width=6.5,height=5,units='in',res=600,
     compression = c("lzw"),
     family="sans")
#png("(1) Figure 1 version 2. Reliability.png", width = 6.5, height = 5, units = 'in', res=300, pointsize=12)
df_figure %>%
  ggplot(aes(x=Variables, y=Estimates, fill=Variables))+
  stat_summary(fun.y=mean, geom="point", shape=18, color="coral", size=3)+
  theme_classic()+
  theme(legend.position="none")+
  scale_y_continuous(breaks=c(0, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1), limits=c(0, 1))+
  labs(x="")+
  annotate(geom="rect", xmin=0, xmax=7, ymin=0, ymax=0.2, fill="brown", alpha=0.2)+
  annotate(geom="rect", xmin=0, xmax=7, ymin=0.21, ymax=0.4, fill="coral", alpha=0.2)+
  annotate(geom="rect", xmin=0, xmax=7, ymin=0.41, ymax=0.6, fill="gold", alpha=0.2)+
  annotate(geom="rect", xmin=0, xmax=7, ymin=0.61, ymax=0.8, fill="limegreen", alpha=0.2)+
  annotate(geom="rect", xmin=0, xmax=7, ymin=0.81, ymax=1, fill="skyblue", alpha=0.2)+
  geom_boxplot(width=.30, fill="gray50")+
  geom_point(size=2, position=position_jitter(width=.02),
             shape=21, fill="gray70")+
  annotate(geom="text", x=0.3, y=0.05, label="Poor")+
  annotate(geom="text", x=0.55, y=0.25, label="Reasonable")+
  annotate(geom="text", x=0.45, y=0.45, label="Moderate")+
  annotate(geom="text", x=0.3, y=0.65, label="Good")+
  annotate(geom="text", x=0.5, y=0.85, label="Very good")
#scale_fill_viridis_d(alpha=.8)+
#scale_color_viridis_d(alpha=.8)
dev.off()

#-----Inter-observer reliability-----
table(df$Observer)
df$Order=paste(df$Video.YT,df$Litter, df$Observer, df$Study)
df=df[order(df$Order),]

df_selec1=df[,c('Total',
                'Head.down',
                'Interaction',
                'Activity',
                'Sits.with.difficulty',
                'Wags.tail'
                
)]

#Looping

#----Observer 1x2------
KappaEstimated=lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df$Observer=="1"],
    x[df$Observer=="2"]),nrow=203,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$kappa,dig=2))

KappaLower = lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df$Observer=="1"],
    x[df$Observer=="2"]),nrow=203,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$boot.quant[[1]], dig=2))

KappaUpper=lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df$Observer=="1"],
    x[df$Observer=="2"]),nrow=203,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$boot.quant[[2]], dig=2))

ICC=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="1"],
    x[df$Observer=="2"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$value,dig=2))

ICC_CI_Lower=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="1"],
    x[df$Observer=="2"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$lbound,dig=2))

ICC_CI_Upper=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="1"],
    x[df$Observer=="2"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$ubound,dig=2))

ICC_Pvalue=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="1"],
    x[df$Observer=="2"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$p.value,dig=40))


col0=c(colnames(df_selec1))
row0=c("KappaEstimated",
       "KappaLower",
       "KappaUpper",
       "ICC",
       "ICC_CI_Lower",
       "ICC_CI_Upper",
       "ICC_Pvalue")

table=matrix(cbind(KappaEstimated,
                   KappaLower,
                   KappaUpper,
                   ICC,
                   ICC_CI_Lower,
                   ICC_CI_Upper,
                   ICC_Pvalue),
             nrow=length(col0),ncol=length(row0),dimnames=list(col0,row0))

write.csv2(table,"(2) Inter-observer reliability - observer 1x2.csv")





#----Observer 1x3------
KappaEstimated=lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df$Observer=="1"],
    x[df$Observer=="3"]),nrow=203,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$kappa,dig=2))

KappaLower = lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df$Observer=="1"],
    x[df$Observer=="3"]),nrow=203,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$boot.quant[[1]], dig=2))

KappaUpper=lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df$Observer=="1"],
    x[df$Observer=="3"]),nrow=203,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$boot.quant[[2]], dig=2))

ICC=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="1"],
    x[df$Observer=="3"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$value,dig=2))

ICC_CI_Lower=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="1"],
    x[df$Observer=="3"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$lbound,dig=2))

ICC_CI_Upper=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="1"],
    x[df$Observer=="3"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$ubound,dig=2))

ICC_Pvalue=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="1"],
    x[df$Observer=="3"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$p.value,dig=40))


col0=c(colnames(df_selec1))
row0=c("KappaEstimated",
       "KappaLower",
       "KappaUpper",
       "ICC",
       "ICC_CI_Lower",
       "ICC_CI_Upper",
       "ICC_Pvalue")

table=matrix(cbind(KappaEstimated,
                   KappaLower,
                   KappaUpper,
                   ICC,
                   ICC_CI_Lower,
                   ICC_CI_Upper,
                   ICC_Pvalue),
             nrow=length(col0),ncol=length(row0),dimnames=list(col0,row0))

write.csv2(table,"(2) Inter-observer reliability - observer 1x3.csv")





#----Observer 1x4------
KappaEstimated=lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df$Observer=="1"],
    x[df$Observer=="4"]),nrow=203,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$kappa,dig=2))

KappaLower = lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df$Observer=="1"],
    x[df$Observer=="4"]),nrow=203,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$boot.quant[[1]], dig=2))

KappaUpper=lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df$Observer=="1"],
    x[df$Observer=="4"]),nrow=203,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$boot.quant[[2]], dig=2))

ICC=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="1"],
    x[df$Observer=="4"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$value,dig=2))

ICC_CI_Lower=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="1"],
    x[df$Observer=="4"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$lbound,dig=2))

ICC_CI_Upper=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="1"],
    x[df$Observer=="4"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$ubound,dig=2))

ICC_Pvalue=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="1"],
    x[df$Observer=="4"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$p.value,dig=40))


col0=c(colnames(df_selec1))
row0=c("KappaEstimated",
       "KappaLower",
       "KappaUpper",
       "ICC",
       "ICC_CI_Lower",
       "ICC_CI_Upper",
       "ICC_Pvalue")

table=matrix(cbind(KappaEstimated,
                   KappaLower,
                   KappaUpper,
                   ICC,
                   ICC_CI_Lower,
                   ICC_CI_Upper,
                   ICC_Pvalue),
             nrow=length(col0),ncol=length(row0),dimnames=list(col0,row0))

write.csv2(table,"(2) Inter-observer reliability - observer 1x4.csv")





#----Observer 1x5------
KappaEstimated=lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df$Observer=="1"],
    x[df$Observer=="5"]),nrow=203,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$kappa,dig=2))

KappaLower = lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df$Observer=="1"],
    x[df$Observer=="5"]),nrow=203,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$boot.quant[[1]], dig=2))

KappaUpper=lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df$Observer=="1"],
    x[df$Observer=="5"]),nrow=203,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$boot.quant[[2]], dig=2))

ICC=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="1"],
    x[df$Observer=="5"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$value,dig=2))

ICC_CI_Lower=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="1"],
    x[df$Observer=="5"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$lbound,dig=2))

ICC_CI_Upper=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="1"],
    x[df$Observer=="5"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$ubound,dig=2))

ICC_Pvalue=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="1"],
    x[df$Observer=="5"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$p.value,dig=40))


col0=c(colnames(df_selec1))
row0=c("KappaEstimated",
       "KappaLower",
       "KappaUpper",
       "ICC",
       "ICC_CI_Lower",
       "ICC_CI_Upper",
       "ICC_Pvalue")

table=matrix(cbind(KappaEstimated,
                   KappaLower,
                   KappaUpper,
                   ICC,
                   ICC_CI_Lower,
                   ICC_CI_Upper,
                   ICC_Pvalue),
             nrow=length(col0),ncol=length(row0),dimnames=list(col0,row0))

write.csv2(table,"(2) Inter-observer reliability - observer 1x5.csv")





#----Observer 1x6------
KappaEstimated=lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df$Observer=="1"],
    x[df$Observer=="6"]),nrow=203,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$kappa,dig=2))

KappaLower = lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df$Observer=="1"],
    x[df$Observer=="6"]),nrow=203,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$boot.quant[[1]], dig=2))

KappaUpper=lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df$Observer=="1"],
    x[df$Observer=="6"]),nrow=203,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$boot.quant[[2]], dig=2))

ICC=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="1"],
    x[df$Observer=="6"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$value,dig=2))

ICC_CI_Lower=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="1"],
    x[df$Observer=="6"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$lbound,dig=2))

ICC_CI_Upper=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="1"],
    x[df$Observer=="6"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$ubound,dig=2))

ICC_Pvalue=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="1"],
    x[df$Observer=="6"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$p.value,dig=40))


col0=c(colnames(df_selec1))
row0=c("KappaEstimated",
       "KappaLower",
       "KappaUpper",
       "ICC",
       "ICC_CI_Lower",
       "ICC_CI_Upper",
       "ICC_Pvalue")

table=matrix(cbind(KappaEstimated,
                   KappaLower,
                   KappaUpper,
                   ICC,
                   ICC_CI_Lower,
                   ICC_CI_Upper,
                   ICC_Pvalue),
             nrow=length(col0),ncol=length(row0),dimnames=list(col0,row0))

write.csv2(table,"(2) Inter-observer reliability - observer 1x6.csv")





#----Observer 2x3------
KappaEstimated=lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df$Observer=="2"],
    x[df$Observer=="3"]),nrow=203,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$kappa,dig=2))

KappaLower = lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df$Observer=="2"],
    x[df$Observer=="3"]),nrow=203,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$boot.quant[[1]], dig=2))

KappaUpper=lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df$Observer=="2"],
    x[df$Observer=="3"]),nrow=203,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$boot.quant[[2]], dig=2))

ICC=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="2"],
    x[df$Observer=="3"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$value,dig=2))

ICC_CI_Lower=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="2"],
    x[df$Observer=="3"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$lbound,dig=2))

ICC_CI_Upper=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="2"],
    x[df$Observer=="3"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$ubound,dig=2))

ICC_Pvalue=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="2"],
    x[df$Observer=="3"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$p.value,dig=40))


col0=c(colnames(df_selec1))
row0=c("KappaEstimated",
       "KappaLower",
       "KappaUpper",
       "ICC",
       "ICC_CI_Lower",
       "ICC_CI_Upper",
       "ICC_Pvalue")

table=matrix(cbind(KappaEstimated,
                   KappaLower,
                   KappaUpper,
                   ICC,
                   ICC_CI_Lower,
                   ICC_CI_Upper,
                   ICC_Pvalue),
             nrow=length(col0),ncol=length(row0),dimnames=list(col0,row0))

write.csv2(table,"(2) Inter-observer reliability - observer 2x3.csv")





#----Observer 2x4------
KappaEstimated=lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df$Observer=="2"],
    x[df$Observer=="4"]),nrow=203,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$kappa,dig=2))

KappaLower = lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df$Observer=="2"],
    x[df$Observer=="4"]),nrow=203,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$boot.quant[[1]], dig=2))

KappaUpper=lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df$Observer=="2"],
    x[df$Observer=="4"]),nrow=203,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$boot.quant[[2]], dig=2))

ICC=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="2"],
    x[df$Observer=="4"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$value,dig=2))

ICC_CI_Lower=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="2"],
    x[df$Observer=="4"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$lbound,dig=2))

ICC_CI_Upper=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="2"],
    x[df$Observer=="4"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$ubound,dig=2))

ICC_Pvalue=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="2"],
    x[df$Observer=="4"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$p.value,dig=40))


col0=c(colnames(df_selec1))
row0=c("KappaEstimated",
       "KappaLower",
       "KappaUpper",
       "ICC",
       "ICC_CI_Lower",
       "ICC_CI_Upper",
       "ICC_Pvalue")

table=matrix(cbind(KappaEstimated,
                   KappaLower,
                   KappaUpper,
                   ICC,
                   ICC_CI_Lower,
                   ICC_CI_Upper,
                   ICC_Pvalue),
             nrow=length(col0),ncol=length(row0),dimnames=list(col0,row0))

write.csv2(table,"(2) Inter-observer reliability - observer 2x4.csv")





#----Observer 2x5------
KappaEstimated=lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df$Observer=="2"],
    x[df$Observer=="5"]),nrow=203,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$kappa,dig=2))

KappaLower = lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df$Observer=="2"],
    x[df$Observer=="5"]),nrow=203,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$boot.quant[[1]], dig=2))

KappaUpper=lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df$Observer=="2"],
    x[df$Observer=="5"]),nrow=203,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$boot.quant[[2]], dig=2))

ICC=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="2"],
    x[df$Observer=="5"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$value,dig=2))

ICC_CI_Lower=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="2"],
    x[df$Observer=="5"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$lbound,dig=2))

ICC_CI_Upper=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="2"],
    x[df$Observer=="5"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$ubound,dig=2))

ICC_Pvalue=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="2"],
    x[df$Observer=="5"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$p.value,dig=40))


col0=c(colnames(df_selec1))
row0=c("KappaEstimated",
       "KappaLower",
       "KappaUpper",
       "ICC",
       "ICC_CI_Lower",
       "ICC_CI_Upper",
       "ICC_Pvalue")

table=matrix(cbind(KappaEstimated,
                   KappaLower,
                   KappaUpper,
                   ICC,
                   ICC_CI_Lower,
                   ICC_CI_Upper,
                   ICC_Pvalue),
             nrow=length(col0),ncol=length(row0),dimnames=list(col0,row0))

write.csv2(table,"(2) Inter-observer reliability - observer 2x5.csv")





#----Observer 2x6------
KappaEstimated=lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df$Observer=="2"],
    x[df$Observer=="6"]),nrow=203,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$kappa,dig=2))

KappaLower = lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df$Observer=="2"],
    x[df$Observer=="6"]),nrow=203,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$boot.quant[[1]], dig=2))

KappaUpper=lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df$Observer=="2"],
    x[df$Observer=="6"]),nrow=203,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$boot.quant[[2]], dig=2))

ICC=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="2"],
    x[df$Observer=="6"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$value,dig=2))

ICC_CI_Lower=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="2"],
    x[df$Observer=="6"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$lbound,dig=2))

ICC_CI_Upper=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="2"],
    x[df$Observer=="6"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$ubound,dig=2))

ICC_Pvalue=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="2"],
    x[df$Observer=="6"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$p.value,dig=40))


col0=c(colnames(df_selec1))
row0=c("KappaEstimated",
       "KappaLower",
       "KappaUpper",
       "ICC",
       "ICC_CI_Lower",
       "ICC_CI_Upper",
       "ICC_Pvalue")

table=matrix(cbind(KappaEstimated,
                   KappaLower,
                   KappaUpper,
                   ICC,
                   ICC_CI_Lower,
                   ICC_CI_Upper,
                   ICC_Pvalue),
             nrow=length(col0),ncol=length(row0),dimnames=list(col0,row0))

write.csv2(table,"(2) Inter-observer reliability - observer 2x6.csv")





#----Observer 3x4------
KappaEstimated=lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df$Observer=="3"],
    x[df$Observer=="4"]),nrow=203,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$kappa,dig=2))

KappaLower = lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df$Observer=="3"],
    x[df$Observer=="4"]),nrow=203,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$boot.quant[[1]], dig=2))

KappaUpper=lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df$Observer=="3"],
    x[df$Observer=="4"]),nrow=203,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$boot.quant[[2]], dig=2))

ICC=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="3"],
    x[df$Observer=="4"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$value,dig=2))

ICC_CI_Lower=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="3"],
    x[df$Observer=="4"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$lbound,dig=2))

ICC_CI_Upper=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="3"],
    x[df$Observer=="4"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$ubound,dig=2))

ICC_Pvalue=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="3"],
    x[df$Observer=="4"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$p.value,dig=40))


col0=c(colnames(df_selec1))
row0=c("KappaEstimated",
       "KappaLower",
       "KappaUpper",
       "ICC",
       "ICC_CI_Lower",
       "ICC_CI_Upper",
       "ICC_Pvalue")

table=matrix(cbind(KappaEstimated,
                   KappaLower,
                   KappaUpper,
                   ICC,
                   ICC_CI_Lower,
                   ICC_CI_Upper,
                   ICC_Pvalue),
             nrow=length(col0),ncol=length(row0),dimnames=list(col0,row0))

write.csv2(table,"(2) Inter-observer reliability - observer 3x4.csv")





#----Observer 3x5------
KappaEstimated=lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df$Observer=="3"],
    x[df$Observer=="5"]),nrow=203,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$kappa,dig=2))

KappaLower = lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df$Observer=="3"],
    x[df$Observer=="5"]),nrow=203,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$boot.quant[[1]], dig=2))

KappaUpper=lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df$Observer=="3"],
    x[df$Observer=="5"]),nrow=203,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$boot.quant[[2]], dig=2))

ICC=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="3"],
    x[df$Observer=="5"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$value,dig=2))

ICC_CI_Lower=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="3"],
    x[df$Observer=="5"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$lbound,dig=2))

ICC_CI_Upper=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="3"],
    x[df$Observer=="5"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$ubound,dig=2))

ICC_Pvalue=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="3"],
    x[df$Observer=="5"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$p.value,dig=40))


col0=c(colnames(df_selec1))
row0=c("KappaEstimated",
       "KappaLower",
       "KappaUpper",
       "ICC",
       "ICC_CI_Lower",
       "ICC_CI_Upper",
       "ICC_Pvalue")

table=matrix(cbind(KappaEstimated,
                   KappaLower,
                   KappaUpper,
                   ICC,
                   ICC_CI_Lower,
                   ICC_CI_Upper,
                   ICC_Pvalue),
             nrow=length(col0),ncol=length(row0),dimnames=list(col0,row0))

write.csv2(table,"(2) Inter-observer reliability - observer 3x5.csv")





#----Observer 3x6------
KappaEstimated=lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df$Observer=="3"],
    x[df$Observer=="6"]),nrow=203,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$kappa,dig=2))

KappaLower = lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df$Observer=="3"],
    x[df$Observer=="6"]),nrow=203,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$boot.quant[[1]], dig=2))

KappaUpper=lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df$Observer=="3"],
    x[df$Observer=="6"]),nrow=203,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$boot.quant[[2]], dig=2))

ICC=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="3"],
    x[df$Observer=="6"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$value,dig=2))

ICC_CI_Lower=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="3"],
    x[df$Observer=="6"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$lbound,dig=2))

ICC_CI_Upper=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="3"],
    x[df$Observer=="6"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$ubound,dig=2))

ICC_Pvalue=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="3"],
    x[df$Observer=="6"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$p.value,dig=40))


col0=c(colnames(df_selec1))
row0=c("KappaEstimated",
       "KappaLower",
       "KappaUpper",
       "ICC",
       "ICC_CI_Lower",
       "ICC_CI_Upper",
       "ICC_Pvalue")

table=matrix(cbind(KappaEstimated,
                   KappaLower,
                   KappaUpper,
                   ICC,
                   ICC_CI_Lower,
                   ICC_CI_Upper,
                   ICC_Pvalue),
             nrow=length(col0),ncol=length(row0),dimnames=list(col0,row0))

write.csv2(table,"(2) Inter-observer reliability - observer 3x6.csv")





#----Observer 4x5------
KappaEstimated=lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df$Observer=="4"],
    x[df$Observer=="5"]),nrow=203,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$kappa,dig=2))

KappaLower = lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df$Observer=="4"],
    x[df$Observer=="5"]),nrow=203,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$boot.quant[[1]], dig=2))

KappaUpper=lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df$Observer=="4"],
    x[df$Observer=="5"]),nrow=203,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$boot.quant[[2]], dig=2))

ICC=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="4"],
    x[df$Observer=="5"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$value,dig=2))

ICC_CI_Lower=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="4"],
    x[df$Observer=="5"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$lbound,dig=2))

ICC_CI_Upper=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="4"],
    x[df$Observer=="5"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$ubound,dig=2))

ICC_Pvalue=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="4"],
    x[df$Observer=="5"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$p.value,dig=40))


col0=c(colnames(df_selec1))
row0=c("KappaEstimated",
       "KappaLower",
       "KappaUpper",
       "ICC",
       "ICC_CI_Lower",
       "ICC_CI_Upper",
       "ICC_Pvalue")

table=matrix(cbind(KappaEstimated,
                   KappaLower,
                   KappaUpper,
                   ICC,
                   ICC_CI_Lower,
                   ICC_CI_Upper,
                   ICC_Pvalue),
             nrow=length(col0),ncol=length(row0),dimnames=list(col0,row0))

write.csv2(table,"(2) Inter-observer reliability - observer 4x5.csv")





#----Observer 4x6------
KappaEstimated=lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df$Observer=="4"],
    x[df$Observer=="6"]),nrow=203,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$kappa,dig=2))

KappaLower = lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df$Observer=="4"],
    x[df$Observer=="6"]),nrow=203,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$boot.quant[[1]], dig=2))

KappaUpper=lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df$Observer=="4"],
    x[df$Observer=="6"]),nrow=203,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$boot.quant[[2]], dig=2))

ICC=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="4"],
    x[df$Observer=="6"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$value,dig=2))

ICC_CI_Lower=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="4"],
    x[df$Observer=="6"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$lbound,dig=2))

ICC_CI_Upper=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="4"],
    x[df$Observer=="6"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$ubound,dig=2))

ICC_Pvalue=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="4"],
    x[df$Observer=="6"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$p.value,dig=40))


col0=c(colnames(df_selec1))
row0=c("KappaEstimated",
       "KappaLower",
       "KappaUpper",
       "ICC",
       "ICC_CI_Lower",
       "ICC_CI_Upper",
       "ICC_Pvalue")

table=matrix(cbind(KappaEstimated,
                   KappaLower,
                   KappaUpper,
                   ICC,
                   ICC_CI_Lower,
                   ICC_CI_Upper,
                   ICC_Pvalue),
             nrow=length(col0),ncol=length(row0),dimnames=list(col0,row0))

write.csv2(table,"(2) Inter-observer reliability - observer 4x6.csv")






#----Observer 5x6------
KappaEstimated=lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df$Observer=="5"],
    x[df$Observer=="6"]),nrow=203,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$kappa,dig=2))

KappaLower = lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df$Observer=="5"],
    x[df$Observer=="6"]),nrow=203,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$boot.quant[[1]], dig=2))

KappaUpper=lapply(df_selec1,function(x)round(
  biostatUZH::confIntKappa(matrix(c(
    x[df$Observer=="5"],
    x[df$Observer=="6"]),nrow=203,ncol=2),type="not Cohen",weights="squared"[1],
    m = 1001,conf.level=.95)$boot.quant[[2]], dig=2))

ICC=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="5"],
    x[df$Observer=="6"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$value,dig=2))

ICC_CI_Lower=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="5"],
    x[df$Observer=="6"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$lbound,dig=2))

ICC_CI_Upper=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="5"],
    x[df$Observer=="6"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$ubound,dig=2))

ICC_Pvalue=lapply(df_selec1,function(x)round(
  irr::icc(matrix(c(
    x[df$Observer=="5"],
    x[df$Observer=="6"]),nrow=203,ncol=2),model="twoway",unit="average",
    type="consistency",conf.level=0.95)$p.value,dig=40))


col0=c(colnames(df_selec1))
row0=c("KappaEstimated",
       "KappaLower",
       "KappaUpper",
       "ICC",
       "ICC_CI_Lower",
       "ICC_CI_Upper",
       "ICC_Pvalue")

table=matrix(cbind(KappaEstimated,
                   KappaLower,
                   KappaUpper,
                   ICC,
                   ICC_CI_Lower,
                   ICC_CI_Upper,
                   ICC_Pvalue),
             nrow=length(col0),ncol=length(row0),dimnames=list(col0,row0))

write.csv2(table,"(2) Inter-observer reliability - observer 5x6.csv")





#### (3) Distribution of scores ####
#----Importing data-------
rm(list=ls())

df = read.csv2("UPAPS Short data.csv")
set.seed(360)

df$Moment[df$Moment == "1"] = "0-1"
df$Moment = factor(df$Moment,
                   levels=c("-24", "0-1", "24"),
                   labels=c("24h before\n surgery",
                            "Immediately \nafter surgery",
                            "24h after\n surgery"))

df$Treatment[df$Treatment == "Iodine"] = "C"
df$Treatment[df$Treatment == "C "] = "C"
df$Treatment = factor(df$Treatment,
                      levels=c("C", "CF"),
                      labels=c("Control",
                               "Flunixin"))

df = df[df$Study != "final",]

#----Control------
df_C=filter(df,Treatment=='Control' &
              Study=='main')

#Head.down
df_proportion_Head.down <- df_C %>%
  group_by(Head.down,Moment) %>%
  summarise(Count = n()) %>%
  group_by(Moment) %>%
  mutate(Proporcao = Count / sum(Count)*100)

df_proportion_Head.down$Head.down = factor(df_proportion_Head.down$Head.down,
                                           levels=c(1, 0))

a = df_proportion_Head.down %>%
  ggplot(aes(x=Moment, y=Proporcao, fill=Head.down))+
  geom_bar(stat="identity", width=0.7,color='black',size=.8)+
  labs(title="A. Head down", fill="Score", y="Percentage",
       x="")+
  scale_y_continuous(n.breaks = 5, limits = c(0,105))+
  scale_fill_viridis_d(alpha=0.50,option = 'D', end=.7)+
  theme_classic()+theme(panel.grid.minor = element_blank(),
                        panel.grid.major = element_blank(),
                        legend.position = 'none')+
  annotate(geom="text", x=1, y=95, label="9.03%", size=3)+
  annotate(geom="text", x=1, y=45, label="90.97%", size=3)+
  annotate(geom="text", x=2, y=60, label="82.64%", size=3)+
  annotate(geom="text", x=2, y=8, label="17.36%", size=3)+
  annotate(geom="text", x=3, y=105, label="6.94%", size=3)+
  annotate(geom="text", x=3, y=46, label="93.06%", size=3);a

#Interaction
df_proportion_Interaction <- df_C %>%
  group_by(Interaction,Moment) %>%
  summarise(Count = n()) %>%
  group_by(Moment) %>%
  mutate(Proporcao = Count / sum(Count)*100)

df_proportion_Interaction$Interaction = factor(df_proportion_Interaction$Interaction,
                                               levels=c(1, 0))

b = df_proportion_Interaction %>%
  ggplot(aes(x=Moment, y=Proporcao, fill=Interaction))+
  geom_bar(stat="identity", width=0.7,color='black',size=.8)+
  labs(title="B. Interaction", fill="Score", y="Percentage",
       x="")+
  scale_y_continuous(n.breaks = 5, limits = c(0,105))+
  scale_fill_viridis_d(alpha=0.50,option = 'D', end=.7)+
  theme_classic()+theme(panel.grid.minor = element_blank(),
                        panel.grid.major = element_blank(),
                        legend.position = 'none')+
  annotate(geom="text", x=1, y=105, label="6.94%", size=3)+
  annotate(geom="text", x=1, y=45, label="93.06%", size=3)+
  annotate(geom="text", x=2, y=55, label="75.00%", size=3)+
  annotate(geom="text", x=2, y=12, label="25.00%", size=3)+
  annotate(geom="text", x=3, y=92, label="17.36%", size=3)+
  annotate(geom="text", x=3, y=41, label="82.64%", size=3);b

#Activity
df_proportion_Activity <- df_C %>%
  group_by(Activity,Moment) %>%
  summarise(Count = n()) %>%
  group_by(Moment) %>%
  mutate(Proporcao = Count / sum(Count)*100)

df_proportion_Activity$Activity = factor(df_proportion_Activity$Activity,
                                         levels=c(1, 0))

c = df_proportion_Activity %>%
  ggplot(aes(x=Moment, y=Proporcao, fill=Activity))+
  geom_bar(stat="identity", width=0.7,color='black',size=.8)+
  labs(title="C. Activity", fill="Score", y="Percentage",
       x="")+
  scale_y_continuous(n.breaks = 5, limits = c(0,105))+
  scale_fill_viridis_d(alpha=0.50,option = 'D', end=.7)+
  theme_classic()+theme(panel.grid.minor = element_blank(),
                        panel.grid.major = element_blank(),
                        legend.position = 'none')+
  annotate(geom="text", x=1, y=93, label="14.58%", size=3)+
  annotate(geom="text", x=1, y=42, label="85.42%", size=3)+
  annotate(geom="text", x=2, y=58, label="77.08%", size=3)+
  annotate(geom="text", x=2, y=10, label="22.92%", size=3)+
  annotate(geom="text", x=3, y=90, label="21.53%", size=3)+
  annotate(geom="text", x=3, y=39, label="78.47%", size=3);c

#Sits.with.difficulty
df_proportion_Sits.with.difficulty <- df_C %>%
  group_by(Sits.with.difficulty,Moment) %>%
  summarise(Count = n()) %>%
  group_by(Moment) %>%
  mutate(Proporcao = Count / sum(Count)*100)

df_proportion_Sits.with.difficulty$Sits.with.difficulty = factor(df_proportion_Sits.with.difficulty$Sits.with.difficulty,
                                                                 levels=c(1, 0))

d = df_proportion_Sits.with.difficulty %>%
  ggplot(aes(x=Moment, y=Proporcao, fill=Sits.with.difficulty))+
  geom_bar(stat="identity", width=0.7,color='black',size=.8)+
  labs(title="D. Sits with difficulty", fill="Score", y="Percentage",
       x="")+
  scale_y_continuous(n.breaks = 5, limits = c(0,105))+
  scale_fill_viridis_d(alpha=0.50,option = 'D', end=.7)+
  theme_classic()+theme(panel.grid.minor = element_blank(),
                        panel.grid.major = element_blank(),
                        legend.position = 'none')+
  annotate(geom="text", x=1, y=90, label="18.75%", size=3)+
  annotate(geom="text", x=1, y=40, label="81.25%", size=3)+
  annotate(geom="text", x=2, y=70, label="60.42%", size=3)+
  annotate(geom="text", x=2, y=20, label="39.58%", size=3)+
  annotate(geom="text", x=3, y=83, label="34.03%", size=3)+
  annotate(geom="text", x=3, y=33, label="65.97%", size=3);d

#Wags.tail
df_proportion_Wags.tail <- df_C %>%
  group_by(Wags.tail,Moment) %>%
  summarise(Count = n()) %>%
  group_by(Moment) %>%
  mutate(Proporcao = Count / sum(Count)*100)

df_proportion_Wags.tail$Wags.tail = factor(df_proportion_Wags.tail$Wags.tail,
                                           levels=c(1, 0))

e = df_proportion_Wags.tail %>%
  ggplot(aes(x=Moment, y=Proporcao, fill=Wags.tail))+
  geom_bar(stat="identity", width=0.7,color='black',size=.8)+
  labs(title="E. Wags tail", fill="Score", y="Percentage",
       x="")+
  scale_y_continuous(n.breaks = 5, limits = c(0,105))+
  scale_fill_viridis_d(alpha=0.50,option = 'D', end=.7)+
  theme_classic()+theme(panel.grid.minor = element_blank(),
                        panel.grid.major = element_blank(),
                        legend.position = 'none')+
  annotate(geom="text", x=1, y=82, label="42.36%", size=3)+
  annotate(geom="text", x=1, y=31, label="57.64%", size=3)+
  annotate(geom="text", x=2, y=75, label="45.14%", size=3)+
  annotate(geom="text", x=2, y=24, label="54.86%", size=3)+
  annotate(geom="text", x=3, y=60, label="79.86%", size=3)+
  annotate(geom="text", x=3, y=10, label="20.14%", size=3);e


legend <- cowplot::get_legend(e)
legend_plot <- as_ggplot(legend)


tiff("(3) Figure 1. Distribution of scores - Control2.tiff",width=6,height=8.5,units='in',res=600,
     compression = c("lzw"),
     family="sans")
png("(3) Figure 1. Distribution of scores - Control2.png", width = 6, height = 8.5,
    units = 'in', res=300, pointsize=12)
a+b+c+d+e+legend+plot_layout(ncol=2)
dev.off()


#----Flunixin-----
df_CF=filter(df,Treatment=='Flunixin' &
               Study=='main')

#Head.down
df_proportion_Head.down <- df_CF %>%
  group_by(Head.down,Moment) %>%
  summarise(Count = n()) %>%
  group_by(Moment) %>%
  mutate(Proporcao = Count / sum(Count)*100)

df_proportion_Head.down$Head.down = factor(df_proportion_Head.down$Head.down,
                                           levels=c(1, 0))

a = df_proportion_Head.down %>%
  ggplot(aes(x=Moment, y=Proporcao, fill=Head.down))+
  geom_bar(stat="identity", width=0.7,color='black',size=.8)+
  labs(title="A. Head down", fill="Score", y="Percentage",
       x="")+
  scale_y_continuous(n.breaks = 5, limits = c(0,105))+
  scale_fill_viridis_d(alpha=0.50,option = 'D', end=.7)+
  theme_classic()+theme(panel.grid.minor = element_blank(),
                        panel.grid.major = element_blank(),
                        legend.position = 'none')+
  annotate(geom="text", x=1, y=94, label="13.19%", size=3)+
  annotate(geom="text", x=1, y=43, label="86.81%", size=3)+
  annotate(geom="text", x=2, y=78, label="63.77%", size=3)+
  annotate(geom="text", x=2, y=18, label="36.23%", size=3)+
  annotate(geom="text", x=3, y=96, label="9.72%", size=3)+
  annotate(geom="text", x=3, y=45, label="90.28%", size=3);a

#Interaction
df_proportion_Interaction <- df_CF %>%
  group_by(Interaction,Moment) %>%
  summarise(Count = n()) %>%
  group_by(Moment) %>%
  mutate(Proporcao = Count / sum(Count)*100)

df_proportion_Interaction$Interaction = factor(df_proportion_Interaction$Interaction,
                                               levels=c(1, 0))

b = df_proportion_Interaction %>%
  ggplot(aes(x=Moment, y=Proporcao, fill=Interaction))+
  geom_bar(stat="identity", width=0.7,color='black',size=.8)+
  labs(title="B. Interaction", fill="Score", y="Percentage",
       x="")+
  scale_y_continuous(n.breaks = 5, limits = c(0,105))+
  scale_fill_viridis_d(alpha=0.50,option = 'D', end=.7)+
  theme_classic()+theme(panel.grid.minor = element_blank(),
                        panel.grid.major = element_blank(),
                        legend.position = 'none')+
  annotate(geom="text", x=1, y=90, label="20.13%", size=3)+
  annotate(geom="text", x=1, y=40, label="79.86%", size=3)+
  annotate(geom="text", x=2, y=77, label="55.80%", size=3)+
  annotate(geom="text", x=2, y=22, label="44.20%", size=3)+
  annotate(geom="text", x=3, y=93, label="15.97%", size=3)+
  annotate(geom="text", x=3, y=42, label="84.03%", size=3);b

#Activity
df_proportion_Activity <- df_CF %>%
  group_by(Activity,Moment) %>%
  summarise(Count = n()) %>%
  group_by(Moment) %>%
  mutate(Proporcao = Count / sum(Count)*100)

df_proportion_Activity$Activity = factor(df_proportion_Activity$Activity,
                                         levels=c(1, 0))

c = df_proportion_Activity %>%
  ggplot(aes(x=Moment, y=Proporcao, fill=Activity))+
  geom_bar(stat="identity", width=0.7,color='black',size=.8)+
  labs(title="C. Activity", fill="Score", y="Percentage",
       x="")+
  scale_y_continuous(n.breaks = 5, limits = c(0,105))+
  scale_fill_viridis_d(alpha=0.50,option = 'D', end=.7)+
  theme_classic()+theme(panel.grid.minor = element_blank(),
                        panel.grid.major = element_blank(),
                        legend.position = 'none')+
  annotate(geom="text", x=1, y=86, label="27.78%", size=3)+
  annotate(geom="text", x=1, y=36, label="72.22%", size=3)+
  annotate(geom="text", x=2, y=76, label="53.62%", size=3)+
  annotate(geom="text", x=2, y=22, label="46.38%", size=3)+
  annotate(geom="text", x=3, y=91, label="18.75%", size=3)+
  annotate(geom="text", x=3, y=40, label="81.25%", size=3);c

#Sits.with.difficulty
df_proportion_Sits.with.difficulty <- df_CF %>%
  group_by(Sits.with.difficulty,Moment) %>%
  summarise(Count = n()) %>%
  group_by(Moment) %>%
  mutate(Proporcao = Count / sum(Count)*100)

df_proportion_Sits.with.difficulty$Sits.with.difficulty = factor(df_proportion_Sits.with.difficulty$Sits.with.difficulty,
                                                                 levels=c(1, 0))

d = df_proportion_Sits.with.difficulty %>%
  ggplot(aes(x=Moment, y=Proporcao, fill=Sits.with.difficulty))+
  geom_bar(stat="identity", width=0.7,color='black',size=.8)+
  labs(title="D. Sits with difficulty", fill="Score", y="Percentage",
       x="")+
  scale_y_continuous(n.breaks = 5, limits = c(0,105))+
  scale_fill_viridis_d(alpha=0.50,option = 'D', end=.7)+
  theme_classic()+theme(panel.grid.minor = element_blank(),
                        panel.grid.major = element_blank(),
                        legend.position = 'none')+
  annotate(geom="text", x=1, y=89, label="24.31%", size=3)+
  annotate(geom="text", x=1, y=38, label="75.69%", size=3)+
  annotate(geom="text", x=2, y=88, label="26.09%", size=3)+
  annotate(geom="text", x=2, y=36, label="73.91%", size=3)+
  annotate(geom="text", x=3, y=88, label="27.78%", size=3)+
  annotate(geom="text", x=3, y=36, label="72.22%", size=3);d

#Wags.tail
df_proportion_Wags.tail <- df_CF %>%
  group_by(Wags.tail,Moment) %>%
  summarise(Count = n()) %>%
  group_by(Moment) %>%
  mutate(Proporcao = Count / sum(Count)*100)

df_proportion_Wags.tail$Wags.tail = factor(df_proportion_Wags.tail$Wags.tail,
                                           levels=c(1, 0))

e = df_proportion_Wags.tail %>%
  ggplot(aes(x=Moment, y=Proporcao, fill=Wags.tail))+
  geom_bar(stat="identity", width=0.7,color='black',size=.8)+
  labs(title="E. Wags tail", fill="Score", y="Percentage",
       x="")+
  scale_y_continuous(n.breaks = 5, limits = c(0,105))+
  scale_fill_viridis_d(alpha=0.50,option = 'D', end=.7)+
  theme_classic()+theme(panel.grid.minor = element_blank(),
                        panel.grid.major = element_blank(),
                        legend.position = 'none')+
  annotate(geom="text", x=1, y=75, label="52.08%", size=3)+
  annotate(geom="text", x=1, y=24, label="47.92%", size=3)+
  annotate(geom="text", x=2, y=76, label="47.10%", size=3)+
  annotate(geom="text", x=2, y=26, label="52.90%", size=3)+
  annotate(geom="text", x=3, y=50, label="77.78%", size=3)+
  annotate(geom="text", x=3, y=10, label="22.22%", size=3);e


legend <- cowplot::get_legend(e)
legend_plot <- as_ggplot(legend)


tiff("(3) Figure 2. Distribution of scores - Flunixin2.tiff",width=6,height=8.5,units='in',res=600,
     compression = c("lzw"),
     family="sans")
#png("(3) Figure 2. Distribution of scores - Flunixin2.png", width = 6, height = 8.5,
#units = 'in', res=300, pointsize=12)
a+b+c+d+e+legend+plot_layout(ncol=2)
dev.off()


#### (4) Multiple association and dimensional structure ####
#-----Importing data----
rm(list=ls())

df = read.csv2("UPAPS Short data.csv")
set.seed(360)

df$Moment[df$Moment == "1"] = "0-1"
df$Moment = factor(df$Moment,
                   levels=c("-24", "0-1", "24"),
                   labels=c("24h before\n surgery",
                            "Immediately \nafter surgery",
                            "24h after\n surgery"))

df$Treatment[df$Treatment == "Iodine"] = "C"
df$Treatment[df$Treatment == "C "] = "C"
df$Treatment = factor(df$Treatment,
                      levels=c("C", "CF"),
                      labels=c("Control",
                               "Flunixin"))

df = df[df$Study == "main",]

# Parallel analysis
df = df[df$Treatment=="Control",]
horn= psych::fa.parallel(df[,c('Head.down',
                               'Interaction',
                               'Activity',
                               'Sits.with.difficulty',
                               'Wags.tail')],fa="fa",n.iter=1001,cor='tet',show.legend=T)

obs = data.frame(horn$fa.values)
obs$type = c('Observed Data')
obs$num = c(row.names(obs))
obs$num = as.numeric(obs$num)
colnames(obs) = c('eigenvalue', 'type', 'num')

#----Parallel analysis (scree plot)-----

percentile = apply(horn$values,2,function(x) quantile(x,.95))
min = as.numeric(nrow(obs))
min = (4*min) - (min-1)
max = as.numeric(nrow(obs))
max = 4*max
percentile1 = percentile[min:max]

sim = data.frame(percentile1)
sim$type = c('Simulated Data')
sim$num = c(row.names(obs))
sim$num = as.numeric(sim$num)
colnames(sim) = c('eigenvalue', 'type', 'num')

eigendat = rbind(obs,sim)

p = ggplot(eigendat, aes(x=num, y=eigenvalue, shape=type, color=type)) +
  geom_line()+
  geom_point(size=2)+
  scale_y_continuous(name='Eigenvalues')+
  scale_x_continuous(name='Factor number', breaks=min(eigendat$num):max(eigendat$num))+
  scale_shape_manual(values=c(16,1)) +
  labs(shape="", color="", tag="A")+
  geom_vline(xintercept = horn$nfact, linetype = 'dashed')+
  scale_color_viridis_d(end=.9, begin=.3)+
  theme_classic();p



tiff("(4) Figure 1. Horns Parallel Analysis.tiff",width=5,height=3,units='in',res=300,
     compression = c("lzw"),
     family="sans")
png("(4) Figure 1. Horns Parallel Analysis.png", width = 5, height = 3,
    units = 'in', res=300, pointsize=12)
p
dev.off()

#----MCA-----

df$Head.down = factor(df$Head.down)
df$Interaction = factor(df$Interaction)
df$Activity = factor(df$Activity)
df$Sits.with.difficulty = factor(df$Sits.with.difficulty)
df$Wags.tail = factor(df$Wags.tail)

mca1=MCA(df[,c(15:19)],graph=T,method="Burt")

contrib = mca1$var[[2]][,1:3]
cos2 = mca1$var[[3]][,1:3]
mca_result = data.frame(contrib, cos2)
#write.csv2(mca_result, "(4) MCA result.csv")
mca1$var
mca1$eig
#write.csv2(mca1$eig, "(4) MCA eigenvalues and variance of dimensions.csv")

biplot1 = fviz_mca_biplot(mca1,
                          axes=c(1,2),
                          geom.ind="point",
                          #pointshape=21,
                          pointsize=3,
                          labelsize=5,
                          geom.var=c("arrow","text"),
                          fill.ind=df$Moment,
                          col.ind=df$Moment,
                          #palette=c("#22A884FF","#E16462FF","#0D0887FF","#F0F921FF"),
                          addEllipses=F,
                          ellipse.type="convex",
                          legend.title="Timepoint",
                          repel=T, 
                          xlab="Factor 1 (Eig: 0.26; Var: 77.73%)",
                          ylab="Factor 2 (Eig: 0.04; Var: 11.88%)",
                          col.var="gray10",
                          title=NULL,
                          ellipse.level=0.95)+
  labs(tag="A")+
  scale_x_continuous(n.breaks = 10)+
  scale_y_continuous(n.breaks = 10)+
  theme_classic()+theme(axis.text=element_text(size=10),
                        axis.title=element_text(size=12),
                        legend.text=element_text(size=12),
                        legend.title = element_text(size=12),
                        panel.grid.minor = element_blank(),
                        panel.grid.major = element_blank())+
  scale_shape_manual(values=c(19,19,19,19))+
  scale_fill_viridis_d(begin=0.10, end=0.98)+
  scale_color_viridis_d(begin=0.10, end=0.98); biplot1



tiff("(4) Figure 2. MCA Perceptual map.tiff",width=7,height=7,units='in',res=300,
     compression = c("lzw"),
     family="sans")
#png("(4) Figure 2. MCA Perceptual map.png", width = 7, height = 7,
#units = 'in', res=300, pointsize=12)
biplot
dev.off()

table.burt=burt(df[,c(15:19)]) 
chi.square.burt=chisq.test(table.burt,simulate.p.value=TRUE);chi.square.burt
chi.square.burt$expected
chi.square.burt$residuals
e=chi.square.burt$stdres
upper.tri(e)
e[upper.tri(e)] <- NA
diag(e)<-NA
write.csv2(e,"(4) Residual.csv")

#----Exploratory factor analysis------
df <- df %>%
  mutate(across(c(Head.down, Interaction, Activity, Sits.with.difficulty, Wags.tail),
                ~ as.numeric(as.character(.))))

EFA=psych::fa(df[,c('Head.down',
                    'Interaction',
                    'Activity',
                    'Sits.with.difficulty',
                    'Wags.tail')],1,rotate="oblimin",fm='ml',cor = 'tet')

EFA$loadings

EFA$e.values

write.csv2(EFA$loadings,"(4) EFA_Oblimin_Loading.csv")
write.csv2(EFA$values,"(4) EFA_eigean.csv")

#One dimension
scale.model <- 'D1=~ Head.down+Interaction+Activity+Sits.with.difficulty+Wags.tail'
onefac8items_a <- cfa(scale.model, data=df,std.lv=TRUE) 
one_=summary(onefac8items_a, fit.measures=TRUE, standardized=TRUE)
write.csv2(one_, "(4) One_CFA.csv")

one_




#----Refinement: Parallel analysis (scree plot)-----
horn= psych::fa.parallel(df[,c('Head.down',
                               'Interaction',
                               'Activity',
                               'Sits.with.difficulty')],fa="fa",n.iter=1001,cor='tet',show.legend=T)


obs = data.frame(horn$fa.values)
obs$type = c('Observed Data')
obs$num = c(row.names(obs))
obs$num = as.numeric(obs$num)
colnames(obs) = c('eigenvalue', 'type', 'num')

#Calculate quantiles for eigenvalues, but only store those from simulated CF model in percentile1
percentile = apply(horn$values,2,function(x) quantile(x,.95))
min = as.numeric(nrow(obs))
min = (4*min) - (min-1)
max = as.numeric(nrow(obs))
max = 4*max
percentile1 = percentile[min:max]

#Create data frame called &amp;amp;amp;amp;amp;quot;sim&amp;amp;amp;amp;amp;quot; with simulated eigenvalue data
sim = data.frame(percentile1)
sim$type = c('Simulated Data')
sim$num = c(row.names(obs))
sim$num = as.numeric(sim$num)
colnames(sim) = c('eigenvalue', 'type', 'num')

#Merge the two data frames (obs and sim) together into data frame called eigendat
eigendat = rbind(obs,sim)

p2 = ggplot(eigendat, aes(x=num, y=eigenvalue, shape=type, color=type)) +
  geom_line()+
  geom_point(size=2)+
  scale_y_continuous(name='Eigenvalues')+
  scale_x_continuous(name='Factor number', breaks=min(eigendat$num):max(eigendat$num))+
  scale_shape_manual(values=c(16,1)) +
  labs(shape="", color="", tag="B")+
  geom_vline(xintercept = horn$nfact, linetype = 'dashed')+
  scale_color_viridis_d(end=.9, begin=.3)+
  theme_classic();p2



tiff("(4) Figure 1. Horns Parallel Analysis A and B.tiff",width=7.5,height=3,units='in',res=300,
     compression = c("lzw"),
     family="sans")
png("(4) Figure 1. Horns Parallel Analysis A and B.png", width = 7.55, height = 3,
    units = 'in', res=300, pointsize=12)
p+p2+plot_layout(guides="collect")
dev.off()
#----Refinement: MCA-----
df$Head.down = factor(df$Head.down)
df$Interaction = factor(df$Interaction)
df$Activity = factor(df$Activity)
df$Sits.with.difficulty = factor(df$Sits.with.difficulty)
df$Wags.tail = factor(df$Wags.tail)

mca2=MCA(df[,c(15:18)],graph=T,method="Burt")

contrib = mca2$var[[2]][,1:3]
cos2 = mca2$var[[3]][,1:3]
mca_result2 = data.frame(contrib, cos2)
write.csv2(mca_result2, "(4) MCA result version 2.csv")

mca2$eig

biplot2 = fviz_mca_biplot(mca2,
                          axes=c(1,2),
                          geom.ind="point",
                          #pointshape=21,
                          pointsize=3,
                          labelsize=5,
                          geom.var=c("arrow","text"),
                          fill.ind=df$Moment,
                          col.ind=df$Moment,
                          #palette=c("#22A884FF","#E16462FF","#0D0887FF","#F0F921FF"),
                          addEllipses=F,
                          ellipse.type="convex",
                          legend.title="Timepoint",
                          repel=T, 
                          xlab="Factor 1 (Eig: 0.41; Var: 88.04%)",
                          ylab="Factor 2 (Eig: 0.04; Var: 8.72%)",
                          col.var="gray10",
                          title=NULL,
                          ellipse.level=0.95)+
  scale_x_continuous(n.breaks = 10)+
  scale_y_continuous(n.breaks = 10)+
  labs(tag="B")+
  theme_classic()+theme(axis.text=element_text(size=10),
                        axis.title=element_text(size=12),
                        legend.text=element_text(size=12),
                        legend.title = element_text(size=12),
                        panel.grid.minor = element_blank(),
                        panel.grid.major = element_blank())+
  scale_shape_manual(values=c(19,19,19,19))+
  scale_fill_viridis_d(begin=0.10, end=0.98)+
  scale_color_viridis_d(begin=0.10, end=0.98); biplot2



tiff("(4) Figure 2. MCA Perceptual map version 2.tiff",width=7,height=8,units='in',res=300,
     compression = c("lzw"),
     family="sans")
png("(4) Figure 2. MCA Perceptual map version 2.png", width = 7, height = 8,
    units = 'in', res=300, pointsize=12)
biplot1+biplot2+plot_layout(guides="collect", nrow=2)
dev.off()


table.burt=burt(df[,c(15:19)]) 
chi.square.burt=chisq.test(table.burt,simulate.p.value=TRUE);chi.square.burt
chi.square.burt$expected
chi.square.burt$residuals
e=chi.square.burt$stdres
upper.tri(e)
e[upper.tri(e)] <- NA
diag(e)<-NA
write.csv2(e,"(4) Residual.csv")

#----Refinement: Exploratory factor analysis------
df <- df %>%
  mutate(across(c(Head.down, Interaction, Activity, Sits.with.difficulty, Wags.tail),
                ~ as.numeric(as.character(.))))

EFA=psych::fa(df[,c('Head.down',
                    'Interaction',
                    'Activity',
                    'Sits.with.difficulty')],1,rotate="oblimin",fm='ml',cor = 'tet')

EFA$loadings

EFA$e.values

write.csv2(EFA$loadings,"(4) EFA_Oblimin_Loading.csv")
write.csv2(EFA$values,"(4) EFA_eigean.csv")

#One dimension
scale.model <- 'D1=~ Head.down+Interaction+Activity+Sits.with.difficulty'
onefac8items_a <- cfa(scale.model, data=df,std.lv=TRUE) 
one_=summary(onefac8items_a, fit.measures=TRUE, standardized=TRUE)
write.csv2(one_, "(4) One_CFA.csv")

one_

#Two dimension correlated - not used
scale.modelc <- 'D1=~ Head.down+Interaction+Activity
                 D2=~ Sits.with.difficulty'

twofac8items_c <- cfa(scale.modelc, data=df,std.lv=TRUE) 
two_c=summary(twofac8items_c, fit.measures=TRUE, standardized=TRUE)
write.csv2(two_c$PE,"Correlated_CFA.csv")

two_c


#### (5) Construct validity & responsiveness ####
#----Importing data----
rm(list=ls())

df = read.csv2("UPAPS Short data.csv")
set.seed(360)

df$Moment[df$Moment == "1"] = "0-1"
df$Moment = factor(df$Moment,
                   levels=c("-24", "0-1", "24"),
                   labels=c("24h before\n surgery",
                            "Immediately \nafter surgery",
                            "24h after\n surgery"))

df$Treatment[df$Treatment == "Iodine"] = "C"
df$Treatment[df$Treatment == "C "] = "C"
df$Treatment = factor(df$Treatment,
                      levels=c("C", "CF"),
                      labels=c("Control",
                               "Flunixin"))

df = df[df$Study == "main",]

df$Experience = case_when(df$Observer == 1 ~ "Some experience",
                          df$Observer == 2 ~ "Little to no experience",
                          df$Observer == 3 ~ "Extensive experience",
                          df$Observer == 4 ~ "Extensive experience",
                          df$Observer == 5 ~ "Little to no experience",
                          df$Observer == 6 ~ "Little to no experience")

df$Experience = factor(df$Experience, levels=c("Little to no experience",
                                               "Some experience",
                                               "Extensive experience"),
                       labels=c("Little to\n no experience",
                                "Some \n experience",
                                "Extensive \n experience"))

df$Gender = case_when(df$Observer == 1 ~ "Male",
                      df$Observer == 2 ~ "Female",
                      df$Observer == 3 ~ "Male",
                      df$Observer == 4 ~ "Female",
                      df$Observer == 5 ~ "Female",
                      df$Observer == 6 ~ "Male")
df$Observer = factor(df$Observer)

df %>% ggplot(aes(y=Total, x=Observer, fill=Observer))+
  geom_boxplot()+
  theme_classic()+
  scale_fill_viridis_d()

df %>% ggplot(aes(y=Total, x=Gender, fill=Gender))+
  geom_boxplot()+
  theme_classic()+
  scale_fill_viridis_d()

df %>% ggplot(aes(y=Total, x=Experience, fill=Experience))+
  geom_boxplot()+
  theme_classic()+
  scale_fill_viridis_d()


wilcox.test(df$Total~df$Experience)

test = kruskal.test(Total~Experience, data=df); test
dunn.test::dunn.test(x=df$Total, g=df$Experience, method="bonferroni")

wilcox.test(Total~Gender, data=df)



#----Checking distribution & setting functions-----------------

hist(df$Total)
mean(df$Total)
median(df$Total)
var(df$Total)
df$Litter = factor(df$Litter)

mod1 = glmer.nb(Total ~ Moment * Treatment +
                  Experience +
                  (1 | Observer) +
                  (1 | Litter),
                data = df)

mod2 = glmer(Total ~ Moment * Treatment +
               Experience +
               (1 | Observer) +
               (1 | Litter),
             family = "poisson",
             data = df)


mod3 = glmer(Total ~ Moment * Treatment +
               Experience +
               (1 | Observer) +
               (1 | Litter),
             family = "poisson",
             control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),
             data = df)

summary(mod1)
summary(mod2)
summary(mod3)


BIC(mod1, mod2, mod3)
AIC(mod1, mod2, mod3)

#Conclusion: Poisson has lower BIC, and will not use observer*moment interaction

#Descriptive statistics function

InSummary = function(x){summarise(df,
                                  Median=round(median(na.omit(x)),2),
                                  FirstQuartile=round(quantile(na.omit(x),type=5,.25),2),
                                  ThirdQuartile=round(quantile(na.omit(x),type=5,.75),2))}



#----Extracting metrics-----

## 1. Total Score
# 1.1 Total Score effects on the model
summary(mod2)
FitModel = mod2

# 1.2 Post hoc test
func_fitmodel(mod2)

cld(lsmeans(glmer(Total ~ Moment * Treatment +
                    Experience +
                    (1 | Observer) +
                    (1 | Litter),
                  family = "poisson",
                  data = df),
            ~Moment*Treatment,p.adjust="bonferroni"),Letters="abcdef",alpha=.05,type="response")

cld(lsmeans(glmer(Total ~ Moment * Treatment +
                    Experience +
                    (1 | Observer) +
                    (1 | Litter),
                  family = "poisson",
                  data = df),
            ~Treatment,p.adjust="bonferroni"),Letters="abcdef",alpha=.05,type="response")

cld(lsmeans(glmer(Total ~ Moment * Treatment +
                    Experience +
                    (1 | Observer) +
                    (1 | Litter),
                  family = "poisson",
                  data = df),
            ~Experience,p.adjust="bonferroni"),Letters="abcdef",alpha=.05,type="response")


# 1.3 Descriptives

for (treatment in c("Control", "Flunixin")) {
  print(treatment)
  for (moment in c("24h before\n surgery",
                   "Immediately \nafter surgery",
                   "24h after\n surgery")) {
    
    print(moment)
    print(InSummary(df$Total[df$Treatment == treatment & df$Moment == moment]))
  }}

## 4. Binary variables

# 4.1 Descriptives

#Head down
for (treatment in c("Control", "Flunixin")) {
  print(treatment)
  for (moment in c("24h before\n surgery",
                   "Immediately \nafter surgery",
                   "24h after\n surgery")) {
    
    print(moment)
    print(InSummary(df$Head.down[df$Treatment == treatment & df$Moment == moment]))
  }}

#Interaction
for (treatment in c("Control", "Flunixin")) {
  print(treatment)
  for (moment in c("24h before\n surgery",
                   "Immediately \nafter surgery",
                   "24h after\n surgery")) {
    
    print(moment)
    print(InSummary(df$Interaction[df$Treatment == treatment & df$Moment == moment]))
  }}

#Activity
for (treatment in c("Control", "Flunixin")) {
  print(treatment)
  for (moment in c("24h before\n surgery",
                   "Immediately \nafter surgery",
                   "24h after\n surgery")) {
    
    print(moment)
    print(InSummary(df$Activity[df$Treatment == treatment & df$Moment == moment]))
  }}

#Sits with difficulty
for (treatment in c("Control", "Flunixin")) {
  print(treatment)
  for (moment in c("24h before\n surgery",
                   "Immediately \nafter surgery",
                   "24h after\n surgery")) {
    
    print(moment)
    print(InSummary(df$Sits.with.difficulty[df$Treatment == treatment & df$Moment == moment]))
  }}

#Wags tail
for (treatment in c("Control", "Flunixin")) {
  print(treatment)
  for (moment in c("24h before\n surgery",
                   "Immediately \nafter surgery",
                   "24h after\n surgery")) {
    
    print(moment)
    print(InSummary(df$Wags.tail[df$Treatment == treatment & df$Moment == moment]))
  }}


# 4.2 Effects on the model and post-hoc tests
lapply(df[,c('Head.down',
             'Interaction',
             'Activity',
             'Sits.with.difficulty',
             'Wags.tail')],function(x)
               list(
                 summary(glmer(x~Moment * Treatment +
                                 Experience +
                                 (1 | Observer) +
                                 (1 | Litter),family='binomial',df)),
                 cld(lsmeans(glmer(x~Moment * Treatment +
                                     Experience +
                                     (1 | Observer) +
                                     (1 | Litter),family='binomial', df),
                             ~Moment*Treatment, p.adjust="bonferroni"),
                     Letters="abcdef",alpha=.05,type="response")))

cld(lsmeans(glmer(Total ~ Moment * Treatment +
                    Experience +
                    (1 | Observer) +
                    (1 | Litter),
                  family = "poisson",
                  data = df),
            ~Experience,p.adjust="bonferroni"),Letters="abcdef",alpha=.05,type="response")


#----Boxplots-----


# 1 - Moments

tiff("(3) BoxPlot 1 - Short UPAPS and Moments tiff.tiff",width=6,height=5,units='in',res=300,
     compression = c("lzw"),
     family="sans")
#png("(3) BoxPlot 1 - UPAPS-SF and Moments.png", width=6, height=6, units='in', res=300)

moments = ggplot(df,aes(y=Total,x=Treatment,fill=Moment))+
  geom_boxplot(show.legend=T)+
  #geom_jitter(width=.23,size=1.5,alpha=.9,color='gray50',fill='gray90',shape=21)+
  stat_summary(fun=mean,geom="point",size=4,shape=18,color="coral", alpha=0.9,position=position_dodge(0.75),show.legend=F)+
  theme_classic()+
  xlab("") +
  ylab(expression(paste("UPAPS-SF total sum")))+
  labs(tag="A", fill="Timepoint")+
  scale_colour_viridis_d(begin=0.10, end=0.9, alpha=.7)+
  scale_fill_viridis_d(begin=0.10, end=0.9, alpha=.7)+
  annotate("text",label="d",x=0.75,y=2.2,colour="black")+
  annotate("text",label="a",x=1,y=5.2,colour="black")+
  annotate("text",label="c",x=1.25,y=3.2,colour="black")+
  annotate("text",label="c",x=1.75,y=5.2,colour="black")+
  annotate("text",label="b",x=2,y=5.2,colour="black")+
  annotate("text",label="c",x=2.25,y=3.2,colour="black"); moments
dev.off()

# Experience
experience = ggplot(df,aes(y=Total,x=Experience, fill=Experience))+
  geom_boxplot(show.legend=F, width=.5)+
  #geom_jitter(width=.23,size=1.5,alpha=.9,color='gray50',fill='gray90',shape=21)+
  stat_summary(fun=mean,geom="point",size=4,shape=18,color="coral", alpha=0.9,position=position_dodge(0.75),show.legend=F)+
  theme_classic()+
  xlab("") +
  ylab(expression(paste("UPAPS-SF total sum")))+
  labs(tag="A")+
  scale_colour_viridis_d(begin=0.10, end=0.9, alpha=.7)+
  scale_fill_viridis_d(alpha=.7, option="E")+
  annotate("text",label="a",x=1,y=5.2,colour="black")+
  annotate("text",label="ab",x=2,y=5.2,colour="black")+
  annotate("text",label="b",x=3,y=3.2,colour="black"); experience


library(patchwork)
tiff("(9) BoxPlots - Timepoints treatment experience suplementar.tiff",width=5,height=8,units='in',res=300,
     compression = c("lzw"),
     family="sans")
#png("(9) BoxPlots - Timepoints treatment experience.png", width=5, height=8, units='in', res=300)
moments+experience+plot_layout(ncol=1)
dev.off()


#----Smoothline plots-------

#Observers

resp_treatment = df %>%
  ggplot(aes(x=as.numeric(Moment),y=Total,color=Treatment))+
  geom_smooth(method="loess",formula=y~x,se=T,alpha=.3)+
  scale_color_viridis_d(begin=.1, end=.9, alpha=.8)+
  labs(x="",
       y="UPAPS-SF total sum",
       tag="A")+
  scale_x_discrete(limits=c("24h before\n surgery",
                            "Immediately \nafter surgery",
                            "24h after\n surgery"), labels=c("24h before\n surgery",
                                                             "Immediately \nafter surgery",
                                                             "24h after\n surgery"))+
  scale_y_continuous(limits=c(0, 5), n.breaks=6)+
  theme_classic(); resp_treatment

resp_experience = df %>%
  ggplot(aes(x=as.numeric(Moment),y=Total,color=Experience))+
  geom_smooth(method="loess",formula=y~x,se=T,alpha=.3)+
  scale_color_viridis_d(begin=.1, end=.9, alpha=.8, option="A")+
  labs(x="",
       y="UPAPS-SF total sum",
       tag="C")+
  scale_x_discrete(limits=c("24h before\n surgery",
                            "Immediately \nafter surgery",
                            "24h after\n surgery"), labels=c("24h before\n surgery",
                                                             "Immediately \nafter surgery",
                                                             "24h after\n surgery"))+
  scale_y_continuous(limits=c(0, 5), n.breaks=6)+
  theme_classic(); resp_experience


df$Litter = factor(df$Litter)

resp_pig = df %>%
  ggplot(aes(x=as.numeric(Moment),y=Total,color=Litter))+
  geom_smooth(method="loess",formula=y~x,se=T,alpha=.1)+
  scale_color_viridis_d(alpha=.8)+
  labs(x="",
       y="UPAPS-SF total sum",
       tag="A",
       color="Pig")+
  scale_x_discrete(limits=c("24h before\n surgery",
                            "Immediately \nafter surgery",
                            "24h after\n surgery"), labels=c("24h before\n surgery",
                                                             "Immediately \nafter surgery",
                                                             "24h after\n surgery"))+
  scale_y_continuous(limits=c(-1, 6.2),
                     breaks=c(0, 1, 2, 3, 4, 5))+
  theme_classic(); resp_pig

resp_gender = df %>%
  ggplot(aes(x=as.numeric(Moment),y=Total,color=Gender))+
  geom_smooth(method="loess",formula=y~x,se=T,alpha=.3)+
  scale_color_viridis_d(end=.9, alpha=.8, option="E")+
  labs(x="",
       y="UPAPS-SF total sum",
       tag="E")+
  scale_x_discrete(limits=c("24h before\n surgery",
                            "Immediately \nafter surgery",
                            "24h after\n surgery"), labels=c("24h before\n surgery",
                                                             "Immediately \nafter surgery",
                                                             "24h after\n surgery"))+
  scale_y_continuous(limits=c(0, 5), n.breaks=6)+
  theme_classic(); resp_gender




tiff("(9) Smoothlines1.tiff",width=5,height=7,units='in',res=300,
     compression = c("lzw"),
     family="sans")
#png("(9) Smoothlines1.png", width=5, height=7, units='in', res=300)
resp_treatment+resp_experience+resp_gender+plot_layout(ncol=1)
dev.off()

tiff("(9) Smoothlines2.tiff",width=6,height=5,units='in',res=300,
     compression = c("lzw"),
     family="sans")
#png("(9) Smoothlines2.png", width=6, height=5, units='in', res=300)
resp_pig
dev.off()


#----Refinement: Extracting metrics--------------
df$Total.new = df$Head.down + df$Interaction + df$Activity + df$Sits.with.difficulty

wilcox.test(df$Total~df$Experience)

test = kruskal.test(Total.new~Experience, data=df); test
dunn.test::dunn.test(x=df$Total, g=df$Experience, method="bonferroni")

wilcox.test(Total.new~Gender, data=df)

mod2_r = glmer(Total.new ~ Moment * Treatment +
                 Experience +
                 (1 | Observer) +
                 (1 | Litter),
               family = "poisson",
               data = df)

cld(lsmeans(glmer(Total.new ~ Moment * Treatment +
                    Experience +
                    (1 | Observer) +
                    (1 | Litter),
                  family = "poisson",
                  data = df),
            ~Moment*Treatment,p.adjust="bonferroni"),Letters="abcdef",alpha=.05,type="response")

cld(lsmeans(glmer(Total.new ~ Moment * Treatment +
                    Experience +
                    (1 | Observer) +
                    (1 | Litter),
                  family = "poisson",
                  data = df),
            ~Experience,p.adjust="bonferroni"),Letters="abcdef",alpha=.05,type="response")



summary(mod2_r)

#----Refinement: Boxplots---------
# 1 - Moments
moments2 = ggplot(df,aes(y=Total.new,x=Treatment,fill=Moment))+
  geom_boxplot(show.legend=T)+
  #geom_jitter(width=.23,size=1.5,alpha=.9,color='gray50',fill='gray90',shape=21)+
  stat_summary(fun=mean,geom="point",size=4,shape=18,color="coral", alpha=0.9,position=position_dodge(0.75),show.legend=F)+
  theme_classic()+
  xlab("") +
  ylab(expression(paste("UPAPS-SF total sum")))+
  labs(tag="B", fill="Timepoint")+
  scale_colour_viridis_d(begin=0.10, end=0.9, alpha=.7)+
  scale_fill_viridis_d(begin=0.10, end=0.9, alpha=.7)+
  annotate("text",label="d",x=0.75,y=2.2,colour="black")+
  annotate("text",label="a",x=1,y=4.2,colour="black")+
  annotate("text",label="c",x=1.25,y=2.2,colour="black")+
  annotate("text",label="c",x=1.75,y=4.2,colour="black")+
  annotate("text",label="b",x=2,y=4.2,colour="black")+
  annotate("text",label="cd",x=2.25,y=2.2,colour="black"); moments2


# Experience
experience2 = ggplot(df,aes(y=Total.new,x=Experience, fill=Experience))+
  geom_boxplot(show.legend=F, width=.5)+
  #geom_jitter(width=.23,size=1.5,alpha=.9,color='gray50',fill='gray90',shape=21)+
  stat_summary(fun=mean,geom="point",size=4,shape=18,color="coral", alpha=0.9,position=position_dodge(0.75),show.legend=F)+
  theme_classic()+
  xlab("") +
  ylab(expression(paste("UPAPS-SF total sum")))+
  labs(tag="B")+
  scale_colour_viridis_d(begin=0.10, end=0.9, alpha=.7)+
  scale_fill_viridis_d(alpha=.7, option="E")+
  annotate("text",label="a",x=1,y=4.2,colour="black")+
  annotate("text",label="a",x=2,y=4.2,colour="black")+
  annotate("text",label="b",x=3,y=4.2,colour="black"); experience2

tiff("(9) BoxPlot 1 - Short UPAPS and Moments tiff version 2.tiff",width=7.5,height=4,units='in',res=300,
     compression = c("lzw"),
     family="sans")
#png("(9) BoxPlot 1 - Short UPAPS and Moments tiff version 2.png", width=7.5, height=4, units='in', res=300)
moments+moments2+plot_layout(guides="collect")
dev.off()

tiff("(9) BoxPlot 2 - Short UPAPS and Experience tiff version 2.tiff",width=7,height=4,units='in',res=300,
     compression = c("lzw"),
     family="sans")
#png("(9) BoxPlot 2 - Short UPAPS and Experience tiff version 2.png", width=7, height=4, units='in', res=300)
experience+experience2
dev.off()


#----Refinement: Smoothline plots-----------------
resp_treatment2 = df %>%
  ggplot(aes(x=as.numeric(Moment),y=Total.new,color=Treatment))+
  geom_smooth(method="loess",formula=y~x,se=T,alpha=.3)+
  scale_color_viridis_d(begin=.1, end=.9, alpha=.8)+
  labs(x="",
       y="UPAPS-SF total sum",
       tag="B")+
  scale_x_discrete(limits=c("24h before\n surgery",
                            "Immediately \nafter surgery",
                            "24h after\n surgery"), labels=c("24h before\n surgery",
                                                             "Immediately \nafter surgery",
                                                             "24h after\n surgery"))+
  scale_y_continuous(limits=c(0, 4), n.breaks=6)+
  theme_classic(); resp_treatment2

resp_experience2 = df %>%
  ggplot(aes(x=as.numeric(Moment),y=Total.new,color=Experience))+
  geom_smooth(method="loess",formula=y~x,se=T,alpha=.3)+
  scale_color_viridis_d(begin=.1, end=.9, alpha=.8, option="A")+
  labs(x="",
       y="UPAPS-SF total sum",
       tag="D")+
  scale_x_discrete(limits=c("24h before\n surgery",
                            "Immediately \nafter surgery",
                            "24h after\n surgery"), labels=c("24h before\n surgery",
                                                             "Immediately \nafter surgery",
                                                             "24h after\n surgery"))+
  scale_y_continuous(limits=c(0, 4), n.breaks=6)+
  theme_classic(); resp_experience2


df$Litter = factor(df$Litter)

resp_pig2 = df %>%
  ggplot(aes(x=as.numeric(Moment),y=Total.new,color=Litter))+
  geom_smooth(method="loess",formula=y~x,se=T,alpha=.1)+
  scale_color_viridis_d(alpha=.8)+
  labs(x="",
       y="UPAPS-SF total sum",
       tag="B",
       color="Pig")+
  scale_x_discrete(limits=c("24h before\n surgery",
                            "Immediately \nafter surgery",
                            "24h after\n surgery"), labels=c("24h before\n surgery",
                                                             "Immediately \nafter surgery",
                                                             "24h after\n surgery"))+
  scale_y_continuous(limits=c(-1, 5.2),
                     breaks=c(0, 1, 2, 3, 4))+
  theme_classic(); resp_pig2

resp_gender2 = df %>%
  ggplot(aes(x=as.numeric(Moment),y=Total.new,color=Gender))+
  geom_smooth(method="loess",formula=y~x,se=T,alpha=.3)+
  scale_color_viridis_d(end=.9, alpha=.8, option="E")+
  labs(x="",
       y="UPAPS-SF total sum",
       tag="F")+
  scale_x_discrete(limits=c("24h before\n surgery",
                            "Immediately \nafter surgery",
                            "24h after\n surgery"), labels=c("24h before\n surgery",
                                                             "Immediately \nafter surgery",
                                                             "24h after\n surgery"))+
  scale_y_continuous(limits=c(0, 4), n.breaks=6)+
  theme_classic(); resp_gender2


tiff("(9) Smoothlines1 version 2.tiff",width=7,height=7,units='in',res=300,
     compression = c("lzw"),
     family="sans")
#png("(9) Smoothlines1 version 2.png", width=7, height=7, units='in', res=300)
resp_treatment+resp_treatment2+
  resp_experience+resp_experience2+
  resp_gender+resp_gender2+plot_layout(ncol=2, guides="collect")
dev.off()


tiff("(9) Smoothlines1 pigs version 2.tiff",width=7,height=8,units='in',res=300,
     compression = c("lzw"),
     family="sans")
#png("(9) Smoothlines1 pigs version 2.png", width=7, height=8, units='in', res=300)
resp_pig+resp_pig2+plot_layout(guides="collect", nrow=2)
dev.off()

#Extracting median

for (treatment in c("Control", "Flunixin")) {
  print(treatment)
  for (moment in c("24h before\n surgery",
                   "Immediately \nafter surgery",
                   "24h after\n surgery")) {
    
    print(moment)
    print(InSummary(df$Total.new[df$Treatment == treatment & df$Moment == moment]))
  }}


#### (6) Item-total correlation ####
#----Importing data----
rm(list=ls())

df = read.csv2("UPAPS Short data.csv")
set.seed(360)

df$Moment[df$Moment == "1"] = "0-1"
df$Moment = factor(df$Moment,
                   levels=c("-24", "0-1", "24"),
                   labels=c("24h before\n surgery",
                            "Immediately \nafter surgery",
                            "24h after\n surgery"))

df$Treatment[df$Treatment == "Iodine"] = "C"
df$Treatment[df$Treatment == "C "] = "C"
df$Treatment = factor(df$Treatment,
                      levels=c("C", "CF"),
                      labels=c("Control",
                               "Flunixin"))

df = df[df$Study == "main",]
df = df[df$Treatment=="Control",]

#----Item-Total---- 



dfIT=data.frame(df[,c('Head.down',
                      'Interaction',
                      'Activity',
                      'Sits.with.difficulty',
                      'Wags.tail')]);str(dfIT) 

dfIT[length(dfIT)+1]=df['Total']-dfIT[1];names(dfIT)[length(dfIT)]=paste0("Exc_",colnames(dfIT)[[1]])
dfIT[length(dfIT)+1]=df['Total']-dfIT[2];names(dfIT)[length(dfIT)]=paste0("Exc_",colnames(dfIT)[[2]])
dfIT[length(dfIT)+1]=df['Total']-dfIT[3];names(dfIT)[length(dfIT)]=paste0("Exc_",colnames(dfIT)[[3]])
dfIT[length(dfIT)+1]=df['Total']-dfIT[4];names(dfIT)[length(dfIT)]=paste0("Exc_",colnames(dfIT)[[4]])
dfIT[length(dfIT)+1]=df['Total']-dfIT[5];names(dfIT)[length(dfIT)]=paste0("Exc_",colnames(dfIT)[[5]])

col0=c(colnames(dfIT[,c(1:5)]), colnames(dfIT[,c(1:5)]))
row0=c("rpb",'P-value')

library("psych")
table=matrix(c(round(corr.test(dfIT[1],dfIT[6],method="spearman")$r[[1]],dig=2),
               round(corr.test(dfIT[2],dfIT[7],method="spearman")$r[[1]],dig=2),
               round(corr.test(dfIT[3],dfIT[8],method="spearman")$r[[1]],dig=2),
               round(corr.test(dfIT[4],dfIT[9],method="spearman")$r[[1]],dig=2),
               round(corr.test(dfIT[5],dfIT[10],method="spearman")$r[[1]],dig=2),
               round(corr.test(df[15],df[20],method="spearman")$r[[1]],dig=2),
               round(corr.test(df[16],df[20],method="spearman")$r[[1]],dig=2),
               round(corr.test(df[17],df[20],method="spearman")$r[[1]],dig=2),
               round(corr.test(df[18],df[20],method="spearman")$r[[1]],dig=2),
               round(corr.test(df[19],df[20],method="spearman")$r[[1]],dig=2),
               
               round(corr.test(dfIT[1],dfIT[6],method="spearman")$p[[1]],dig=100),
               round(corr.test(dfIT[2],dfIT[7],method="spearman")$p[[1]],dig=100),
               round(corr.test(dfIT[3],dfIT[8],method="spearman")$p[[1]],dig=100),
               round(corr.test(dfIT[4],dfIT[9],method="spearman")$p[[1]],dig=100),
               round(corr.test(dfIT[5],dfIT[10],method="spearman")$p[[1]],dig=100),
               round(corr.test(df[15],df[20],method="spearman")$p[[1]],dig=100),
               round(corr.test(df[16],df[20],method="spearman")$p[[1]],dig=100),
               round(corr.test(df[17],df[20],method="spearman")$p[[1]],dig=100),
               round(corr.test(df[18],df[20],method="spearman")$p[[1]],dig=100),
               round(corr.test(df[19],df[20],method="spearman")$p[[1]],dig=100)
),
nrow=length(col0),ncol=length(row0),dimnames=list(col0,row0))

write.csv2(table,"(6) Item Total correlation.csv")

#----Refinement: Item-Total---- 

df$Total.new = df$Head.down + df$Interaction + df$Activity + df$Sits.with.difficulty


dfIT=data.frame(df[,c('Head.down',
                      'Interaction',
                      'Activity',
                      'Sits.with.difficulty')]);str(dfIT) 

dfIT[length(dfIT)+1]=df['Total.new']-dfIT[1];names(dfIT)[length(dfIT)]=paste0("Exc_",colnames(dfIT)[[1]])
dfIT[length(dfIT)+1]=df['Total.new']-dfIT[2];names(dfIT)[length(dfIT)]=paste0("Exc_",colnames(dfIT)[[2]])
dfIT[length(dfIT)+1]=df['Total.new']-dfIT[3];names(dfIT)[length(dfIT)]=paste0("Exc_",colnames(dfIT)[[3]])
dfIT[length(dfIT)+1]=df['Total.new']-dfIT[4];names(dfIT)[length(dfIT)]=paste0("Exc_",colnames(dfIT)[[4]])

col0=c(colnames(dfIT[,c(1:4)]), colnames(dfIT[,c(1:4)]))
row0=c("Rho",'P-value')

library("psych")

table=matrix(c(round(corr.test(dfIT[1],dfIT[5],method="spearman")$r[[1]],dig=2),
               round(corr.test(dfIT[2],dfIT[6],method="spearman")$r[[1]],dig=2),
               round(corr.test(dfIT[3],dfIT[7],method="spearman")$r[[1]],dig=2),
               round(corr.test(dfIT[4],dfIT[8],method="spearman")$r[[1]],dig=2),
               round(corr.test(df[15],df[23],method="spearman")$r[[1]],dig=2),
               round(corr.test(df[16],df[23],method="spearman")$r[[1]],dig=2),
               round(corr.test(df[17],df[23],method="spearman")$r[[1]],dig=2),
               round(corr.test(df[18],df[23],method="spearman")$r[[1]],dig=2),
               
               round(corr.test(dfIT[1],dfIT[5],method="spearman")$p[[1]],dig=100),
               round(corr.test(dfIT[2],dfIT[6],method="spearman")$p[[1]],dig=100),
               round(corr.test(dfIT[3],dfIT[7],method="spearman")$p[[1]],dig=100),
               round(corr.test(dfIT[4],dfIT[8],method="spearman")$p[[1]],dig=100),
               round(corr.test(df[15],df[23],method="spearman")$p[[1]],dig=100),
               round(corr.test(df[16],df[23],method="spearman")$p[[1]],dig=100),
               round(corr.test(df[17],df[23],method="spearman")$p[[1]],dig=100),
               round(corr.test(df[18],df[23],method="spearman")$p[[1]],dig=100)
),
nrow=length(col0),ncol=length(row0),dimnames=list(col0,row0))

write.csv2(table,"(6) Item Total correlation refined.csv")

#### (7) Internal consistency ####
#----Importing data----
rm(list=ls())

df = read.csv2("UPAPS Short data.csv")
set.seed(360)

df$Moment[df$Moment == "1"] = "0-1"
df$Moment = factor(df$Moment,
                   levels=c("-24", "0-1", "24"),
                   labels=c("24h before\n surgery",
                            "Immediately \nafter surgery",
                            "24h after\n surgery"))

df$Treatment[df$Treatment == "Iodine"] = "C"
df$Treatment[df$Treatment == "C "] = "C"
df$Treatment = factor(df$Treatment,
                      levels=c("C", "CF"),
                      labels=c("Control",
                               "Flunixin"))

df = df[df$Study == "main",]
df = df[df$Treatment=="Control",]

#----Alpha----
alpha=data.frame(df[,c('Head.down',
                       'Interaction',
                       'Activity',
                       'Sits.with.difficulty',
                       'Wags.tail')]);str(alpha) 

col0=c('All',colnames(alpha))
row0=c("GM")

alfa.table=matrix(c(
  psy::cronbach(alpha)$alpha,
  psy::cronbach(alpha[,-1])$alpha,
  psy::cronbach(alpha[,-2])$alpha,
  psy::cronbach(alpha[,-3])$alpha,
  psy::cronbach(alpha[,-4])$alpha,
  psy::cronbach(alpha[,-5])$alpha
),
nrow=length(col0),ncol=length(row0),dimnames=list(col0,row0))

write.csv2(round(alfa.table,dig=2),"(5) Alpha_Internal_Consistency.csv")


#----Omega----

col0=c('All',colnames(alpha))
row0=c("GM")

omega=matrix(c(
  psych::omega(alpha,fm='pc',nfactors=1,covar=F,rotation='none', poly=T)[[4]],
  psych::omega(alpha[,-1],fm='pc',nfactors=1,covar=F,rotation='none', poly=T)[[4]],
  psych::omega(alpha[,-2],fm='pc',nfactors=1,covar=F,rotation='none', poly=T)[[4]],
  psych::omega(alpha[,-3],fm='pc',nfactors=1,covar=F,rotation='none', poly=T)[[4]],
  psych::omega(alpha[,-4],fm='pc',nfactors=1,covar=F,rotation='none', poly=T)[[4]],
  psych::omega(alpha[,-5],fm='pc',nfactors=1,covar=F,rotation='none', poly=T)[[4]]
),
nrow=length(col0),ncol=length(row0),dimnames=list(col0,row0))

write.csv2(round(omega,2),"(5) Omega_Internal_Consistency1.csv")

#----Refinement: Alpha----
alpha=data.frame(df[,c('Head.down',
                       'Interaction',
                       'Activity',
                       'Sits.with.difficulty')]);str(alpha) 

col0=c('All',colnames(alpha))
row0=c("GM")

alfa.table=matrix(c(
  psy::cronbach(alpha)$alpha,
  psy::cronbach(alpha[,-1])$alpha,
  psy::cronbach(alpha[,-2])$alpha,
  psy::cronbach(alpha[,-3])$alpha,
  psy::cronbach(alpha[,-4])$alpha
),
nrow=length(col0),ncol=length(row0),dimnames=list(col0,row0))

write.csv2(round(alfa.table,dig=2),"(5) Alpha_Internal_Consistency refined.csv")


#----Refinement: Omega----

col0=c('All',colnames(alpha))
row0=c("GM")

omega=matrix(c(
  psych::omega(alpha,fm='pc',nfactors=1,covar=F,rotation='none', poly=T)[[4]],
  psych::omega(alpha[,-1],fm='pc',nfactors=1,covar=F,rotation='none', poly=T)[[4]],
  psych::omega(alpha[,-2],fm='pc',nfactors=1,covar=F,rotation='none', poly=T)[[4]],
  psych::omega(alpha[,-3],fm='pc',nfactors=1,covar=F,rotation='none', poly=T)[[4]],
  psych::omega(alpha[,-4],fm='pc',nfactors=1,covar=F,rotation='none', poly=T)[[4]]
),
nrow=length(col0),ncol=length(row0),dimnames=list(col0,row0))

write.csv2(round(omega,2),"(5) Omega_Internal_Consistency1 refined.csv")

#### (8) Optimal cut-off point ####
#----Importing data----
rm(list=ls())

df = read.csv2("UPAPS Short data.csv")
set.seed(360)

df$Moment[df$Moment == "1"] = "0-1"
df$Moment = factor(df$Moment,
                   levels=c("-24", "0-1", "24"),
                   labels=c("24h before\n surgery",
                            "Immediately \nafter surgery",
                            "24h after\n surgery"))

df$Treatment[df$Treatment == "Iodine"] = "C"
df$Treatment[df$Treatment == "C "] = "C"
df$Treatment = factor(df$Treatment,
                      levels=c("C", "CF"),
                      labels=c("Control",
                               "Flunixin"))

df = df[df$Study == "main",]
df = df[df$Treatment=="Control",]

df = df[df$Moment == "24h before\n surgery"| df$Moment == "Immediately \nafter surgery",]

df$Condition = ifelse(df$Moment == "24h before\n surgery", 0, 1)

#----ROC------

roccurve=roc(Condition~Total,data=df,plot=T,
             algorithm=2,smooth=F,boot.n=1001,boot.stratified=T,ci.auc=T,auc=T)

roccurve$auc
ci.auc(roccurve)

roc_coords = ci.coords(roccurve,x="best",
                       input=c("threshold", "specificity", "sensitivity"),
                       ret=c("threshold", "specificity", "sensitivity"),
                       best.method="youden",
                       best.policy = "random",
                       conf.level=0.95, boot.n=1001,
                       boot.stratified=TRUE)

roc_coords


#----Plots-----

#ROC UPAPS-SF
predicoes <- prediction(predictions = df$Total,labels = df$Condition) #atualizar variáveis aqui
dados_curva_roc <- performance(predicoes, measure = "sens") 
sensitividade <- (performance(predicoes, measure = "sens"))@y.values[[1]] 
especificidade <- (performance(predicoes, measure = "spec"))@y.values[[1]]
cutoffs <- dados_curva_roc@x.values[[1]] 
dados_plotagem <- cbind.data.frame(cutoffs, especificidade, sensitividade)

p=dados_plotagem %>%
  ggplot(aes(x = cutoffs, y = especificidade)) +
  ggtitle("A")+
  geom_vline(xintercept = 1.5,color='#f9a242ff',size=.8,linetype = "dashed")+
  annotate("rect",xmin=1.5,xmax=2.5,ymin=0,ymax=1,alpha=.2,fill="gray")+
  annotate("text",x=1.6,y=.3,label="Optimal cut-off point",angle=90,color='#f9a242ff',size=4)+
  annotate("text",x=2.6,y=.4,label="Diagnostic Uncertainly Zone",angle=90,color='gray',size=4)+
  #annotate("text",x=1,y=.03,label="True Negative",angle=0,color='black',size=3)+
  #annotate("text",x=5.5,y=.03,label="True Positive",angle=0,color='black',size=3)+
  annotate("text",x=4.5,y=.9,label="Specificity",angle=0,color='#95D840FF',size=4)+
  annotate("text",x=.7,y=.9,label="Sensitivity",angle=0,color='#440154FF',size=4)+
  geom_line(size = 1,color = "#95D840FF") +
  #geom_point(color = "#95D840FF",size = 1.9) +
  geom_line(aes(x = cutoffs, y = sensitividade),size = 1,color = "#440154FF") +
  #geom_point(aes(x = cutoffs, y = sensitividade),color = "#440154FF",size = 1.9) +
  labs(x = "UPAPS-SF total sum",y = "Specificity and Sensitivity") +
  scale_color_manual("Legend:",values = c("#95D840FF", "#440154FF")) +
  scale_x_continuous(n.breaks = 5)+
  scale_y_continuous(n.breaks = 10,
                     labels=c("0.0", 0.1, 0.2, 0.3, 0.4, 0.5,
                              0.6, 0.7, 0.8, 0.9, "1.0"))+
  theme_classic()+theme(axis.text=element_text(size=10),
                        axis.title=element_text(size=12),
                        legend.text=element_text(size=12),
                        legend.title = element_text(size=12),
                        legend.position = "bottom");p

obj=roc(df$Condition,df$Total, ci=T, plot=F) 
ciobj=ci.se(obj,specificities=seq(0,1,l=25))
dat.ci=data.frame(x=as.numeric(rownames(ciobj)),lower=ciobj[,1],upper=ciobj[,3])


a = ggroc(obj,color='steelblue',size=1.2)+ 
  theme_minimal()+ 
  geom_abline(slope=1,intercept=1,linetype="dashed",alpha=.7,color="grey")+ 
  coord_equal()+ 
  geom_ribbon(data=dat.ci,aes(x=x,ymin=lower,ymax=upper),fill="steelblue",alpha=.2)+
  ggtitle("B")+ 
  xlab("1 – Specificity") +
  ylab("Sensitivity")+
  scale_y_continuous(n.breaks=10,
                     labels=c("0.0", 0.1, 0.2, 0.3, 0.4, 0.5,
                              0.6, 0.7, 0.8, 0.9, "1.0"))+
  scale_x_reverse(n.breaks=10,
                  labels=c("0.0", 0.1, 0.2, 0.3, 0.4, 0.5,
                           0.6, 0.7, 0.8, 0.9, "1.0"))+
  theme_classic()+theme(axis.text=element_text(size=10),
                        axis.title=element_text(size=12),
                        legend.text=element_text(size=12),
                        legend.title = element_text(size=12))+
  annotate(geom="text", x=0.2, y=0.1, label="AUC: 91.63\n(88.33 - 94.92)");a

auc(roccurve)
ci.auc(roccurve)

tiff("(7) ROC and 2 line.tiff",width=10,height=7,units='in',res=300,
     compression = c("lzw"),
     family="sans")
#png("(7) ROC and 2 line.png", width=10, height=7, units='in', res=300)
p+a#+plot_layout(nrow=2, ncol=1)
dev.off()


#----Refinement: ROC-----
df$Total.new = df$Head.down + df$Interaction + df$Activity + df$Sits.with.difficulty

roccurve2=roc(Condition~Total.new,data=df,plot=T,
              algorithm=2,smooth=F,boot.n=1001,boot.stratified=T,ci.auc=T,auc=T)

roccurve2$auc
ci.auc(roccurve2)

roc_coords2 = ci.coords(roccurve2,x="best",
                        input=c("threshold", "specificity", "sensitivity"),
                        ret=c("threshold", "specificity", "sensitivity"),
                        best.method="youden",
                        best.policy = "random",
                        conf.level=0.95, boot.n=1001,
                        boot.stratified=TRUE)

roc_coords2

#----Refinement: Plots-----

#ROC UPAPS-SF
predicoes <- prediction(predictions = df$Total.new,labels = df$Condition)
dados_curva_roc <- performance(predicoes, measure = "sens") 
sensitividade <- (performance(predicoes, measure = "sens"))@y.values[[1]] 
especificidade <- (performance(predicoes, measure = "spec"))@y.values[[1]]
cutoffs <- dados_curva_roc@x.values[[1]] 
dados_plotagem <- cbind.data.frame(cutoffs, especificidade, sensitividade)

p2=dados_plotagem %>%
  ggplot(aes(x = cutoffs, y = especificidade)) +
  ggtitle("C")+
  geom_vline(xintercept = 1.5,color='#f9a242ff',size=.8,linetype = "dashed")+
  annotate("rect",xmin=1.5,xmax=2.5,ymin=0,ymax=1,alpha=.2,fill="gray")+
  annotate("text",x=1.6,y=.3,label="Optimal cut-off point",angle=90,color='#f9a242ff',size=4)+
  annotate("text",x=2.6,y=.4,label="Diagnostic Uncertainly Zone",angle=90,color='gray',size=4)+
  #annotate("text",x=1,y=.03,label="True Negative",angle=0,color='black',size=3)+
  #annotate("text",x=5.5,y=.03,label="True Positive",angle=0,color='black',size=3)+
  annotate("text",x=3.5,y=.9,label="Specificity",angle=0,color='#95D840FF',size=4)+
  annotate("text",x=.7,y=.9,label="Sensitivity",angle=0,color='#440154FF',size=4)+
  geom_line(size = 1,color = "#95D840FF") +
  #geom_point(color = "#95D840FF",size = 1.9) +
  geom_line(aes(x = cutoffs, y = sensitividade),size = 1,color = "#440154FF") +
  #geom_point(aes(x = cutoffs, y = sensitividade),color = "#440154FF",size = 1.9) +
  labs(x = "UPAPS-SF total sum",y = "Specificity and Sensitivity") +
  scale_color_manual("Legend:",values = c("#95D840FF", "#440154FF")) +
  scale_x_continuous(n.breaks = 5)+
  scale_y_continuous(n.breaks = 10,
                     labels=c("0.0", 0.1, 0.2, 0.3, 0.4, 0.5,
                              0.6, 0.7, 0.8, 0.9, "1.0"))+
  theme_classic()+theme(axis.text=element_text(size=10),
                        axis.title=element_text(size=12),
                        legend.text=element_text(size=12),
                        legend.title = element_text(size=12),
                        legend.position = "bottom");p2

obj=roc(df$Condition,df$Total.new, ci=T, plot=F) 
ciobj=ci.se(obj,specificities=seq(0,1,l=25))
dat.ci=data.frame(x=as.numeric(rownames(ciobj)),lower=ciobj[,1],upper=ciobj[,3])


a2 = ggroc(obj,color='steelblue',size=1.2)+ 
  theme_minimal()+ 
  geom_abline(slope=1,intercept=1,linetype="dashed",alpha=.7,color="grey")+ 
  coord_equal()+ 
  geom_ribbon(data=dat.ci,aes(x=x,ymin=lower,ymax=upper),fill="steelblue",alpha=.2)+
  ggtitle("D")+ 
  xlab("1 – Specificity") +
  ylab("Sensitivity")+
  scale_y_continuous(n.breaks=10,
                     labels=c("0.0", 0.1, 0.2, 0.3, 0.4, 0.5,
                              0.6, 0.7, 0.8, 0.9, "1.0"))+
  scale_x_reverse(n.breaks=10,
                  labels=c("0.0", 0.1, 0.2, 0.3, 0.4, 0.5,
                           0.6, 0.7, 0.8, 0.9, "1.0"))+
  theme_classic()+theme(axis.text=element_text(size=10),
                        axis.title=element_text(size=12),
                        legend.text=element_text(size=12),
                        legend.title = element_text(size=12))+
  annotate(geom="text", x=0.2, y=0.1, label="AUC: 92.96\n(89.96 - 95.96)");a2

auc(roccurve2)
ci.auc(roccurve2)

tiff("(7) ROC and 2 line version 2.tiff",width=7.5,height=7.5,units='in',res=300,
     compression = c("lzw"),
     family="sans")
#png("(7) ROC and 2 line version 2.png", width=7.5, height=7.5, units='in', res=300)
p+a+p2+a2+plot_layout(nrow=2, ncol=2)
dev.off()


#### (9) Sensitivity and specifity ####
#----Importing data----
rm(list=ls())

df = read.csv2("UPAPS Short data.csv")
set.seed(360)

df$Moment[df$Moment == "1"] = "0-1"
df$Moment = factor(df$Moment,
                   levels=c("-24", "0-1", "24"),
                   labels=c("24h before\n surgery",
                            "Immediately \nafter surgery",
                            "24h after\n surgery"))

df$Treatment[df$Treatment == "Iodine"] = "C"
df$Treatment[df$Treatment == "C "] = "C"
df$Treatment = factor(df$Treatment,
                      levels=c("C", "CF"),
                      labels=c("Control",
                               "Flunixin"))

df = df[df$Study == "main",]
df = df[df$Treatment=="Control",]

df = df[df$Moment == "24h before\n surgery"| df$Moment == "Immediately \nafter surgery",]

df$Condition = ifelse(df$Moment == "24h before\n surgery", 0, 1)

df = df[df$Total != 2,]
#----Specificity and sensitivity-------------
df_select = df[,c(15:20, 23)]

df_select$Total[df_select$Total < 2] <- 0
df_select$Total[df_select$Total > 2] <- 1


table(df_select$Total, df_select$Condition)

Sp_estimate=lapply(df_select, function(x) round(
  epi.tests(table(df_select$Condition,x),conf.level=.95)$detail$est[
    epi.tests(table(df_select$Condition,x),conf.level=.95)$detail$statistic=="sp"],dig=4)*100)

Sp_low=lapply(df_select, function(x) round(
  epi.tests(table(df_select$Condition,x),conf.level=.95)$detail$lower[
    epi.tests(table(df_select$Condition,x),conf.level=.95)$detail$statistic=="sp"],dig=4)*100)

Sp_high=lapply(df_select, function(x) round(
  epi.tests(table(df_select$Condition,x),conf.level=.95)$detail$upper[
    epi.tests(table(df_select$Condition,x),conf.level=.95)$detail$statistic=="sp"],dig=4)*100)

Se_estimate=lapply(df_select, function(x) round(
  epi.tests(table(df_select$Condition,x),conf.level=.95)$detail$est[
    epi.tests(table(df_select$Condition,x),conf.level=.95)$detail$statistic=="se"],dig=4)*100)

Se_low=lapply(df_select, function(x) round(
  epi.tests(table(df_select$Condition,x),conf.level=.95)$detail$lower[
    epi.tests(table(df_select$Condition,x),conf.level=.95)$detail$statistic=="se"],dig=4)*100)

Se_high=lapply(df_select, function(x) round(
  epi.tests(table(df_select$Condition,x),conf.level=.95)$detail$upper[
    epi.tests(table(df_select$Condition,x),conf.level=.95)$detail$statistic=="se"],dig=4)*100)

col0=c(colnames(df_select))
row0=c("Sp_estimate",
       "Sp_low",
       "Sp_high",
       "Se_estimate",
       "Se_low",
       "Se_high")

table=matrix(cbind(Sp_estimate,
                   Sp_low,
                   Sp_high,
                   Se_estimate,
                   Se_low,
                   Se_high),
             ncol=length(row0), nrow=length(col0),dimnames=list(col0,row0))

write.csv2(table,"(8) Sp and Se.csv")


#----Refinement: Specificity and sensitivity-------------
df$Total.new = df$Head.down + df$Interaction + df$Activity + df$Sits.with.difficulty
df = df[df$Total.new != 2,]


df_select = df[,c(15:18, 24, 23)]

df_select$Total.new[df_select$Total.new < 2] <- 0
df_select$Total.new[df_select$Total.new > 2] <- 1


table(df_select$Total.new, df_select$Condition)

Sp_estimate=lapply(df_select, function(x) round(
  epi.tests(table(df_select$Condition,x),conf.level=.95)$detail$est[
    epi.tests(table(df_select$Condition,x),conf.level=.95)$detail$statistic=="sp"],dig=4)*100)

Sp_low=lapply(df_select, function(x) round(
  epi.tests(table(df_select$Condition,x),conf.level=.95)$detail$lower[
    epi.tests(table(df_select$Condition,x),conf.level=.95)$detail$statistic=="sp"],dig=4)*100)

Sp_high=lapply(df_select, function(x) round(
  epi.tests(table(df_select$Condition,x),conf.level=.95)$detail$upper[
    epi.tests(table(df_select$Condition,x),conf.level=.95)$detail$statistic=="sp"],dig=4)*100)

Se_estimate=lapply(df_select, function(x) round(
  epi.tests(table(df_select$Condition,x),conf.level=.95)$detail$est[
    epi.tests(table(df_select$Condition,x),conf.level=.95)$detail$statistic=="se"],dig=4)*100)

Se_low=lapply(df_select, function(x) round(
  epi.tests(table(df_select$Condition,x),conf.level=.95)$detail$lower[
    epi.tests(table(df_select$Condition,x),conf.level=.95)$detail$statistic=="se"],dig=4)*100)

Se_high=lapply(df_select, function(x) round(
  epi.tests(table(df_select$Condition,x),conf.level=.95)$detail$upper[
    epi.tests(table(df_select$Condition,x),conf.level=.95)$detail$statistic=="se"],dig=4)*100)

col0=c(colnames(df_select))
row0=c("Sp_estimate",
       "Sp_low",
       "Sp_high",
       "Se_estimate",
       "Se_low",
       "Se_high")

table=matrix(cbind(Sp_estimate,
                   Sp_low,
                   Sp_high,
                   Se_estimate,
                   Se_low,
                   Se_high),
             ncol=length(row0), nrow=length(col0),dimnames=list(col0,row0))

write.csv2(table,"(8) Sp and Se refined.csv")


#### (10) Criterion validity ####
#----Importing data----
rm(list=ls())

df = read.csv2("UPAPS Short data.csv")
set.seed(360)

df$Moment[df$Moment == "1"] = "0-1"
df$Moment = factor(df$Moment,
                   levels=c("-24", "0-1", "24"),
                   labels=c("24h before\n surgery",
                            "Immediately \nafter surgery",
                            "24h after\n surgery"))

df$Treatment[df$Treatment == "Iodine"] = "C"
df$Treatment[df$Treatment == "C "] = "C"
df$Treatment = factor(df$Treatment,
                      levels=c("C", "CF"),
                      labels=c("Control",
                               "Flunixin"))

df = df[df$Study == "main",]
df$Litter = factor(df$Litter)
df$ID = paste(df$Litter, df$Moment)

df$Total.new = df$Head.down + df$Interaction + df$Activity + df$Sits.with.difficulty

#-----Predictive criterion validity---------
#Negative at M1
#Before refinement
nrow(df[df$Moment == "24h before\n surgery"
        & df$Total < 2
        & df$Observer == 1,]) / 48 * 100

nrow(df[df$Moment == "24h before\n surgery"
        & df$Total < 2
        & df$Observer == 2,]) / 48 * 100

nrow(df[df$Moment == "24h before\n surgery"
        & df$Total < 2
        & df$Observer == 3,]) / 48 * 100

nrow(df[df$Moment == "24h before\n surgery"
        & df$Total < 2
        & df$Observer == 4,]) / 48 * 100

nrow(df[df$Moment == "24h before\n surgery"
        & df$Total < 2
        & df$Observer == 5,]) / 48 * 100

nrow(df[df$Moment == "24h before\n surgery"
        & df$Total < 2
        & df$Observer == 6,]) / 48 * 100

#After refinement
nrow(df[df$Moment == "24h before\n surgery"
        & df$Total.new < 2
        & df$Observer == 1,]) / 48 * 100

nrow(df[df$Moment == "24h before\n surgery"
        & df$Total.new < 2
        & df$Observer == 2,]) / 48 * 100

nrow(df[df$Moment == "24h before\n surgery"
        & df$Total.new < 2
        & df$Observer == 3,]) / 48 * 100

nrow(df[df$Moment == "24h before\n surgery"
        & df$Total.new < 2
        & df$Observer == 4,]) / 48 * 100

nrow(df[df$Moment == "24h before\n surgery"
        & df$Total.new < 2
        & df$Observer == 5,]) / 48 * 100

nrow(df[df$Moment == "24h before\n surgery"
        & df$Total.new < 2
        & df$Observer == 6,]) / 48 * 100







#Positive at M2

#Before refinement
nrow(df[df$Moment == "Immediately \nafter surgery"
        & df$Total > 2
        & df$Observer == 1
        & df$Treatment == "Control",]) / 24 * 100

nrow(df[df$Moment == "Immediately \nafter surgery"
        & df$Total > 2
        & df$Observer == 2
        & df$Treatment == "Control",]) / 24 * 100

nrow(df[df$Moment == "Immediately \nafter surgery"
        & df$Total > 2
        & df$Observer == 3
        & df$Treatment == "Control",]) / 24 * 100

nrow(df[df$Moment == "Immediately \nafter surgery"
        & df$Total > 2
        & df$Observer == 4
        & df$Treatment == "Control",]) / 24 * 100

nrow(df[df$Moment == "Immediately \nafter surgery"
        & df$Total > 2
        & df$Observer == 5
        & df$Treatment == "Control",]) / 24 * 100

nrow(df[df$Moment == "Immediately \nafter surgery"
        & df$Total > 2
        & df$Observer == 6
        & df$Treatment == "Control",]) / 24 * 100

#After refinement
nrow(df[df$Moment == "Immediately \nafter surgery"
        & df$Total.new > 2
        & df$Observer == 1
        & df$Treatment == "Control",]) / 24 * 100

nrow(df[df$Moment == "Immediately \nafter surgery"
        & df$Total.new > 2
        & df$Observer == 2
        & df$Treatment == "Control",]) / 24 * 100

nrow(df[df$Moment == "Immediately \nafter surgery"
        & df$Total.new > 2
        & df$Observer == 3
        & df$Treatment == "Control",]) / 24 * 100

nrow(df[df$Moment == "Immediately \nafter surgery"
        & df$Total.new > 2
        & df$Observer == 4
        & df$Treatment == "Control",]) / 24 * 100

nrow(df[df$Moment == "Immediately \nafter surgery"
        & df$Total.new > 2
        & df$Observer == 5
        & df$Treatment == "Control",]) / 24 * 100

nrow(df[df$Moment == "Immediately \nafter surgery"
        & df$Total.new > 2
        & df$Observer == 6
        & df$Treatment == "Control",]) / 24 * 100



table(df$Treatment, df$Moment)






#Diagnostic uncertainty zone
nrow(df[df$Total == 2
        & df$Observer == 1,]) / 143 * 100

nrow(df[df$Total == 2 
        & df$Observer == 2,]) / 143 * 100

nrow(df[df$Total == 2
        & df$Observer == 3,]) / 143 * 100

nrow(df[df$Total == 2
        & df$Observer == 4,]) / 143 * 100

nrow(df[df$Total == 2
        & df$Observer == 5,]) / 143 * 100

nrow(df[df$Total == 2
        & df$Observer == 6,]) / 143 * 100

#After refinement
nrow(df[df$Total.new == 2
        & df$Observer == 1,]) / 143 * 100

nrow(df[df$Total.new == 2
        & df$Observer == 2,]) / 143 * 100

nrow(df[df$Total.new == 2
        & df$Observer == 3,]) / 143 * 100

nrow(df[df$Total.new == 2
        & df$Observer == 4,]) / 143 * 100

nrow(df[df$Total.new == 2
        & df$Observer == 5,]) / 143 * 100

nrow(df[df$Total.new == 2
        & df$Observer == 6,]) / 143 * 100


hist(df$Total.new)

hist(df2$Total.Pain.Score)


