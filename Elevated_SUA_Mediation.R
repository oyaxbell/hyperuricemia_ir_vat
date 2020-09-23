#######################################################################
###ENSANUT-NHANES joint cohort analysis. SIGMA confirmatory analysis###
#######################################################################
#1. Package library####
library(alr3); library(tidyverse); library(dplyr); library(mice); library(readxl); 
library(corrplot); library(nortest); library(lmtest); library(dunn.test);library(mediation);
library(bestNormalize); library(boot); library(simpleboot); library(randtests);
library(pROC); library(OptimalCutpoints); library(caret); library(ggpubr); 
library(cowplot); library(gridExtra); library(GmAMisc); library(interactions);
library(nhanesA); library(qpcR); library(gvlma); library(ggdag); library(lme4); 
library(lmerTest); library(knitr); library(sjPlot); library(broom); library(stargazer); 
library(effects)

#2. Additional functions module####

Chsq <- function(x){
  x <- matrix(x,byrow =TRUE,nrow=1)
  return(chisq.test(x)$p.value)
}

ttest <- function(x,y){
  return(t.test(x,y)$p.value)
}

ktest<-function(x,y){
  return(kruskal.test(x~y)$p.value)
}

categorize<-function(x){
  c1<-c()
  name<-deparse(substitute(x))
  for (value in x){
    if (value==1){
      cat=1
    }else{
      cat=0
    }
    c1<-c(c1, cat)
  }
  return(name<-c1)
}

categorizeinv<-function(x){
  c1<-c()
  name<-deparse(substitute(x))
  for (value in x){
    if (value==0){
      cat=0
    }else{
      cat=1
    }
    c1<-c(c1, cat)
  }
  return(name<-c1)
}

categ<-function(x, lim){
  category<-c()
  for(value in x){
    if(value<lim){
      cat=0
    }else{
      cat=1
    }
    category<-c(category,cat)
  }
  return(category)
}

inv.categ<-function(x, lim){
  category<-c()
  for(value in x){
    if(value>lim){
      cat=0
    }else{
      cat=1
    }
    category<-c(category,cat)
  }
  return(category)
}

cap<-function(x){
  qnt <- quantile(x, probs=c(.25, .75), na.rm = T)
  caps <- quantile(x, probs=c(.05, .95), na.rm = T)
  H <- 1.5 * IQR(x, na.rm = T)
  x[x < (qnt[1] - H)] <- caps[1]
  x[x > (qnt[2] + H)] <- caps[2]
  return(x)
}

#3. NHANES/ENSANUT JOINT COHORT####
##3.1 NHANES import and database management####
DEMO.2004<-nhanes("DEMO_C")%>%dplyr::select(SEQN, RIAGENDR, RIDAGEYR, RIDRETH1)

lab10.2004<-nhanes("L10_C")%>%dplyr::select(SEQN, LBXGH)
lab10AM.2004<-nhanes("L10AM_C")%>%dplyr::select(SEQN, LBXGLU, LBXIN)
lab13ct.2004<-nhanes("l13_c")%>%dplyr::select(SEQN, LBXTC)
lab13hdl.2004<-nhanes("l13_c")%>%dplyr::select(SEQN, LBXHDD)
lab13AM.2004<-nhanes("L13AM_C")%>%dplyr::select(SEQN, LBXTR)
lab18.2004<-nhanes("L40_C")%>%dplyr::select(SEQN, LBXSUA)

labs.2004<-merge(lab10.2004, 
                 merge(lab10AM.2004, 
                       merge(lab13ct.2004, 
                             merge(lab13hdl.2004, 
                                   merge(lab18.2004, lab13AM.2004, by="SEQN", all.x=T), 
                                   by="SEQN", all.x=T), 
                             by="SEQN", all.x=T), 
                       by="SEQN", all.x=T), 
                 by="SEQN", all.x=T)

BMX.2004<-nhanes("BMX_C")%>%dplyr::select(SEQN, BMXWT, BMXHT, BMXWAIST)

DXA.2004<-nhanesDXA(2003)%>%dplyr::select(SEQN, DXXTRFAT)

DEMO.2012<-nhanes("DEMO_G")%>%dplyr::select(SEQN, RIAGENDR, RIDAGEYR, RIDRETH1)

lab10.2012<-nhanes("GHB_G")%>%dplyr::select(SEQN, LBXGH)
lab10AM.2012<-nhanes("GLU_G")%>%dplyr::select(SEQN, LBXGLU, LBXIN)
lab13ct.2012<-nhanes("TCHOL_G")%>%dplyr::select(SEQN, LBXTC)
lab13hdl.2012<-nhanes("HDL_G")%>%dplyr::select(SEQN, LBDHDD)
lab13AM.2012<-nhanes("TRIGLY_G")%>%dplyr::select(SEQN, LBXTR)
lab18.2012<-nhanes("BIOPRO_G")%>%dplyr::select(SEQN, LBXSUA)

labs.2012<-merge(lab10.2012, 
                 merge(lab10AM.2012, 
                       merge(lab13ct.2012, 
                             merge(lab13hdl.2012, 
                                   merge(lab18.2012, lab13AM.2012, by="SEQN", all.x=T), 
                                   by="SEQN", all.x=T), 
                             by="SEQN", all.x=T), 
                       by="SEQN", all.x=T), 
                 by="SEQN", all.x=T)

BMX.2012<-nhanes("BMX_G")%>%dplyr::select(SEQN, BMXWT, BMXHT, BMXWAIST)

data_names<-c("ID", "sex", "age", "ethnicity", "glycohemoglobin", "glucose", "insulin",
              "tcholesterol", "cHDL", "uric_acid", "TGL", "weight", "height", "waist")

au.nhanes.2012<-merge(DEMO.2012, 
                      merge(labs.2012, BMX.2012, by="SEQN", all.x=T), 
                      by="SEQN", all.x=T)
au.nhanes.2012.1<-au.nhanes.2012[!duplicated(au.nhanes.2012$SEQN),]
colnames(au.nhanes.2012.1)<-data_names
au.nhanes.2012.1$dataset<-2

au.nhanes.2004<-merge(DEMO.2004, 
                      merge(labs.2004, BMX.2004, by="SEQN", all.x=T), 
                      by="SEQN", all.x=T)
au.nhanes.2004.1<-au.nhanes.2004[!duplicated(au.nhanes.2004$SEQN),]
colnames(au.nhanes.2004.1)<-data_names
au.nhanes.2004.1$dataset<-1

au.nhanes<-rbind(au.nhanes.2004.1, au.nhanes.2012.1)

au.nhanes1<-au.nhanes[!duplicated(au.nhanes$ID),]
au.nhanes2<-au.nhanes1[!is.na(au.nhanes1$uric_acid),]

md.pattern(au.nhanes2)
set.seed(123)
imp<-mice(au.nhanes2, m=1, maxit=1)
d.au.nhanes<-complete(imp, "long")
md.pattern(d.au.nhanes)

write.csv(d.au.nhanes, "d.au.nhanes.csv") 

nhanes.homa<-read_excel("NHANES_ENSANUT_HOMA.xls")
homas.nhanes<-nhanes.homa%>%dplyr::select(ID, HOMA2IR, HOMA2B, HOMA2S)

d.au.nhanes1<-merge(d.au.nhanes, homas.nhanes, by="ID", all.x=T)
d.au.nhanes1$BMI<-as.numeric(d.au.nhanes1$weight/((d.au.nhanes1$height)/100)^2)

d.au.nhanes1$sex<-categorize(d.au.nhanes1$sex)

d.au.nhanes2<-d.au.nhanes1%>%dplyr::select(ID, sex, age, ethnicity, glycohemoglobin, 
                                           glucose, insulin, tcholesterol, cHDL, uric_acid, TGL, 
                                           weight, height, waist, HOMA2IR, HOMA2B, HOMA2S, BMI,
                                           dataset)


ensanuthecha<-read_excel("~/DOCUMENTOS/UNAM/TESIS/Hiperuricemia/ensanuthecha.xlsx")
au.ensanut<-ensanuthecha %>% dplyr::select(id, sexo.x, edad.x, etnia, valor_HB1AC, valor_GLU_SUERO, 
                                           valor_INSULINA, valor_COLEST, valor_COL_HDL,
                                           AC_UR_MG_DL, valor_TRIG, peso, prom_talla, 
                                           prom_cintura, HOMA2IR, HOMA2B, HOMA2S, imc2)

data_names1<-c("ID", "sex", "age", "ethnicity", "glycohemoglobin", "glucose", "insulin",
               "tcholesterol", "cHDL", "uric_acid", "TGL", "weight", "height", "waist",
               "HOMA2IR", "HOMA2B", "HOMA2S", "BMI")
colnames(au.ensanut)<-data_names1

md.pattern(au.ensanut)
set.seed(123)
imp<-mice(au.ensanut, m=1, maxit=1)
au.ensanut.1<-complete(imp, "long")
md.pattern(au.ensanut.1)

au.ensanut1<-au.ensanut.1%>%dplyr::select(ID, sex, age, ethnicity, glycohemoglobin, 
                                          glucose, insulin, tcholesterol, cHDL, uric_acid, TGL, 
                                          weight, height, waist, HOMA2IR, HOMA2B, HOMA2S, BMI)
au.ensanut1$dataset<-0

nhanes.ensanut<-as.data.frame(rbind(au.ensanut1, d.au.nhanes2))
md.pattern(nhanes.ensanut)

write.csv(nhanes.ensanut, "nhanes_ensanut.csv") 

##3.2 Variables####
data<-nhanes.ensanut 

data[, c(1:18)]<-sapply(data[, c(1:18)], as.numeric)
data$METSIR<-as.numeric((log(2*data$glucose+data$TGL)*data$BMI)/log(data$cHDL))
data$WHR<-as.numeric(data$waist)/as.numeric(data$height)
data$METSVF<-(4.466+0.011*(log(data$METSIR)^3)+3.239*(log(data$WHR)^3)+0.319*(data$sex)
              +0.594*(log(data$age)))

md.pattern(data)
set.seed(123)
imp<-mice(data, m=1, maxit=1)
data1<-complete(imp, "long")
md.pattern(data1)

ID<-data1$ID
sex<-data1$sex
age<-data1$age
ethnicity<-data1$ethnicity
glycohemoglobin<-data1$glycohemoglobin
glucose<-data1$glucose
insulin<-data1$insulin
tcholesterol<-data1$tcholesterol
cHDL<-data1$cHDL
uric_acid<-data1$uric_acid
TGL<-data1$TGL
weight<-data1$weight
height<-data1$height
waist<-data1$waist
HOMA2IR<-data1$HOMA2IR
HOMA2B<-data1$HOMA2B
HOMA2S<-data1$HOMA2S
BMI<-data1$BMI
METSIR<-data1$METSIR
METSVF<-data1$METSVF
WHR<-data1$WHR

data1$HOMA2IRcat<-categ(HOMA2IR, 2.5)
HOMA2IRcat<-data1$HOMA2IRcat

data1$METSVFcat<-categ(METSVF, 7.18)
METSVFcat<-data1$METSVFcat

data1$UAcat<-categ(uric_acid, 5.5)
UAcat<-data1$UAcat

c1<-c()
for(value in BMI){
  if(value<18.54){
    cat=0
  }else{
    if(value<24.9){
      cat=1
    }else{
      if(value<29.9){
        cat=2
      }else{
        if(value<40){
          cat=3
        }else{
          cat=4
        }
      }
    }
  }
  c1<-c(c1,cat)
}
data1$BMIcat<-c1
BMIcat<-data1$BMIcat

diab<-categ(glycohemoglobin, 6.5)+categ(glucose, 126)
data1$diabetes<-categorize(diab)
diabetes<-data1$diabetes

METSVF.IR<-2*METSVFcat+HOMA2IRcat
data1$METSVF.IR<-METSVF.IR

data1$homa2irhealthy<-factor(data1$HOMA2IRcat, labels=c("Non resistant", 
                                                        "Resistant"))
data1$hiperuricemia<-factor(data1$UAcat, labels=c("Good", "Sick"))
data1$fat<-factor(data1$METSVFcat, labels=c("Non obese", "Obese"))
data1$imc<-factor(data1$BMIcat, labels=c("Underweight", "Normal", "Overweight",
                                         "Obese", "Morbidly obese"))
data1$diabetic<-factor(data1$diabetes, labels=c("Non diabetic", "Diabetic"))
data1$sex1<-factor(data1$sex, labels=c("Female", "Male"))

write.csv(data1, "nhanes_ensanut_cat.csv") 

##3.3 Variable analysis####

log_insulin<-log(insulin)

log_h2ir<-log(HOMA2IR)
sqrt_h2b<-sqrt(HOMA2B)
log_h2s<-log(HOMA2S)
log_ua<-log(uric_acid)

qnt <- quantile(log_h2ir, probs=c(.25, .75), na.rm = T)
caps <- quantile(log_h2ir, probs=c(.05, .95), na.rm = T)
H <- 1.5 * IQR(log_h2ir, na.rm = T)
log_h2ir[log_h2ir < (qnt[1] - H)] <- caps[1]
log_h2ir[log_h2ir > (qnt[2] + H)] <- caps[2]

qnt <- quantile(log_ua, probs=c(.25, .75), na.rm = T)
caps <- quantile(log_ua, probs=c(.05, .95), na.rm = T)
H <- 1.5 * IQR(log_h2ir, na.rm = T)
log_ua[log_ua < (qnt[1] - H)] <- caps[1]
log_ua[log_ua > (qnt[2] + H)] <- caps[2]

o<-orderNorm(METSVF)
o_METSVF<-o$x.t

df1<-cbind(uric_acid, log_h2ir, log_h2s, sqrt_h2b, o_METSVF, log_insulin)
res1<-cor.mtest(df1, conf.level = .95)
colnames(df1)<-c("Uric\nacid", "HOMA2\nIR", "HOMA2\n%S", "HOMA2\n%B", "METS-\nVF", 
                 "Insulin")
cor1<-cor(df1, method="spearman", use = "complete.obs")
corrplot.mixed(cor1,lower.col = "black", number.cex = 1, upper = "circle",
               p.mat=res1$p, sig.level=.05)


##3.4 Regression Analysis. Graphs for best model fit.####
#1

r.1<-lmer(log_h2ir~o_METSVF+ethnicity+(1|data1$dataset))
summary(r.1)

r.2<-lmer(log_h2ir~poly(o_METSVF,3)+factor(ethnicity)+(1|data1$dataset))
summary(r.2)

g.r2<-ggplot(data1, aes(x=o_METSVF, y=log_h2ir))+
  geom_point(alpha=0.5)+
  stat_smooth(method = "lm", formula = y ~ poly(x,3), size = 1)+
  ylab("Ln HOMA2-IR")+
  xlab("ON METS-VF")+
  theme_classic()
print(g.r2)

#2
r.3<-lm(log_ua~o_METSVF+ethnicity)
summary(r.3)

r.4<-lmer(log_ua~poly(o_METSVF,3)+factor(ethnicity)+(1|data1$dataset))
summary(r.4)

g.r4<-ggplot(data1, aes(x=o_METSVF, y=log_ua))+
  geom_point(alpha=0.5)+
  stat_smooth(method = "lm", formula = y ~ poly(x,3), size = 1)+
  ylab("Ln uric acid")+
  xlab("ON METS-VF")+
  theme_classic()
print(g.r4)

#3
r.5<-lm(log_h2ir~log_ua+ethnicity+age+sex)
summary(r.5)

r.6<-lmer(log_h2ir~poly(log_ua,2)+factor(ethnicity)+age+sex+(1|data1$dataset))
summary(r.6)

g.r6<-ggplot(data1, aes(x=log_ua, y=log_h2ir))+
  geom_point(alpha=0.5)+
  stat_smooth(method = "lm", formula = y ~ poly(x,2), size = 1)+
  ylab("Ln HOMA2-IR")+
  xlab("Ln uric acid")+
  theme_classic()
print(g.r6)

#4

r.7<-lm(o_METSVF~log_h2ir+ethnicity)
summary(r.7)

r.8<-lmer(o_METSVF~poly(log_h2ir,3)+factor(ethnicity)+(1|data1$dataset))
summary(r.8)

g.r8<-ggplot(data1, aes(x=log_h2ir, y=o_METSVF))+
  geom_point(alpha=0.5)+
  stat_smooth(method = "lm", formula = y ~ poly(x,3), size = 1)+
  ylab("ON METS-VF")+
  xlab("Ln HOMA2-IR")+
  theme_classic()
print(g.r8)

#5

r.9<-lmer(log_ua~log_h2ir+factor(ethnicity)+age+sex+(1|data1$dataset))
summary(r.9)

r.10<-lmer(log_ua~poly(log_h2ir,2)+factor(ethnicity)+age+sex+(1|data1$dataset))
summary(r.10)

g.r9<-ggplot(data1, aes(x=log_h2ir, y=log_ua))+
  geom_point(alpha=0.5)+
  stat_smooth(method = "lm", formula = y ~ poly(x,1), size = 1)+
  ylab("Ln uric acid")+
  xlab("Ln HOMA2-IR")+
  theme_classic()
print(g.r9)

#6

r.11<-lm(o_METSVF~log_ua+ethnicity)
summary(r.11)

r.12<-lmer(o_METSVF~poly(log_ua, 3)+factor(ethnicity)+(1|data1$dataset))
summary(r.12)


g.r12<-ggplot(data1, aes(x=log_ua, y=o_METSVF))+
  geom_point(alpha=0.5)+
  stat_smooth(method = "lm", formula = y ~ poly(x,3), size = 1)+
  xlab("Ln uric acid")+
  ylab("ON METS-VF")+
  theme_classic()
print(g.r12)

reg.1<-list(r.2, r.4, r.6, r.8, r.9, r.12)
regresiones<-reg.sum(reg.1)

reg1<-tab_model(r.2, r.4, r.6, r.8, r.8, r.12)

#7

r.1.1<-lm(log_h2ir~o_METSVF+log_ua+ethnicity)
summary(r.1.1)

r.1.2<-lmer(log_h2ir~poly(o_METSVF,3)+log_ua+factor(ethnicity)+(1|data1$dataset))
summary(r.1.2)

#8

r.1.3<-lm(log_ua~o_METSVF+log_h2ir+ethnicity)
summary(r.1.3)

r.1.4<-lmer(log_ua~poly(o_METSVF,3)+log_h2ir+factor(ethnicity)+(1|data1$dataset))
summary(r.1.4)

#9

r.1.5<-lmer(log_h2ir~log_ua+o_METSVF+factor(ethnicity)+(1|data1$dataset))
summary(r.1.5)

r.1.6<-lm(log_h2ir~poly(log_ua, 3)+o_METSVF+ethnicity)
summary(r.1.6)

#10

r.1.7<-lm(o_METSVF~log_h2ir+ethnicity+log_ua)
summary(r.1.7)

r.1.8<-lmer(o_METSVF~poly(log_h2ir,3)+factor(ethnicity)+log_ua+(1|data1$dataset))
summary(r.1.8)

#11

r.1.9<-lmer(log_ua~log_h2ir+factor(ethnicity)+o_METSVF+(1|data1$dataset))
summary(r.1.9)

r.1.10<-lm(log_ua~poly(log_h2ir, 3)+ethnicity+o_METSVF)
summary(r.1.10)

#12

r.1.11<-lm(o_METSVF~log_ua+ethnicity+log_h2ir)
summary(r.1.11)

r.1.12<-lmer(o_METSVF~poly(log_ua,3)+ethnicity+log_h2ir+(1|data1$dataset))
summary(r.1.12). Gra

r.arrange<-ggarrange(g.r2, g.r4, g.r6, g.r8, g.r9, g.r12, 
                     labels=c("A", "B", "C", "D", "E", "F"), nrow=2, ncol=3)
reg2<-tab_model(r.1.2, r.1.4, r.1.5, r.1.8, r.1.9, r.1.12, collapse.ci = T)

##3.5 Mediation analysis####

detach("package:lmerTest", unload = TRUE)

#C~A - 1
set.seed(123)
ac.1<-lmer(log_ua~log_h2ir+age+sex+factor(ethnicity)+(1|data1$dataset))
summary(ac.1)

#B~A+C - 1
bac.1<-lmer(o_METSVF~log_h2ir+log_ua+factor(ethnicity)+(1|data1$dataset))
summary(bac.1)

#Mediation - 1
med.1<- mediate(ac.1, bac.1, treat = "log_h2ir", 
                mediator = "log_ua", sims = 1000, boot.ci.type = "perc", boot = F)
summary(med.1)

#C~A - 2
set.seed(123)
ac.4<-lmer(log_ua~o_METSVF+factor(ethnicity)+(1|data1$dataset))
summary(ac.4)

#B~A+C - 2
bac.4<-lmer(log_h2ir~o_METSVF+log_ua+factor(ethnicity)+(1|data1$dataset))
summary(bac.4)

#Mediation - 2
med.4<-mediate(ac.4, bac.4, treat = "o_METSVF", 
               mediator = "log_ua", sims = 1000, boot.ci.type = "perc", boot = F)
summary(med.4)

##3.6 Fatty acids (AdipoIR) analysis####

###3.6.1 Database management
FFA.2004<-nhanes("SSFA_C")%>%dplyr::select(SEQN, SSPM1.N)
colnames(FFA.2004)<-c("ID", "palmitate")
FFA.2012<-nhanes("FAS_G")%>%dplyr::select(SEQN, LBXPM1)
colnames(FFA.2012)<-c("ID", "palmitate")

FFA<-rbind(FFA.2004, FFA.2012)

d.au.nhanes3<-merge(data1, FFA, by="ID")
d.au.nhanes3<-d.au.nhanes3[!duplicated(d.au.nhanes3$ID),]

d.nhanes.ffa<-d.au.nhanes3%>%dplyr::select(ID, uric_acid, palmitate, glucose, insulin, HOMA2S, 
                                           HOMA2IR, HOMA2B, sex, age, METSVF, palmitate, ethnicity, 
                                           dataset, METSVF.IR)

d.nhanes.ffa$insulinSI<-d.nhanes.ffa$insulin*0.144

d.nhanes.ffa$adipoIR<-(d.nhanes.ffa$insulinSI*d.nhanes.ffa$palmitate)*4

ffa.homa<-read_excel("~/DOCUMENTOS/UNAM/TESIS/Hiperuricemia/ffahoma.xlsx")%>%
  dplyr::select(HOMA2IR, HOMA2B)

md.pattern(d.nhanes.ffa)
set.seed(123)
imp<-mice(d.nhanes.ffa, m=1, maxit=1)
data0<-complete(imp, "long")
md.pattern(data0)

write.csv(data0, "nhanesffa.csv") #Final database

###3.6.2 Variables and analysis

fa.uric_acid<-cap(log(data0$uric_acid))
fa.palmitate<-log(data0$palmitate)
fa.adipoIR<-cap(log(data0$adipoIR))
fa.METSVF<-data0$METSVF
fa.sex<-data0$sex
fa.age<-data0$age
fa.eth<-data0$ethnicity
fa.insulin<-cap(log(data0$insulin))
fa.HOMA2IR<-log(data0$HOMA2IR)
fa.HOMA2B<-log(data0$HOMA2B)
fa.dataset<-data0$dataset

o.f<-orderNorm(fa.METSVF)
fa.o_METSVF<-o.f$x.t

data0$log_adipoIR<-log(data0$adipoIR, base=10)
data0$fa.adipoIRcat<-categ(data0$log_adipoIR, 4.6)
data0$UAcat<-categ(data0$uric_acid, 5.5)
data0$METSVFcat<-categ(data0$METSVF, 7.18)
data0$HOMA2IRcat<-categ(data0$HOMA2IR, 2.5)

df2<-cbind(fa.uric_acid, fa.HOMA2IR, fa.adipoIR, fa.HOMA2B, fa.o_METSVF, fa.insulin)
res2<-cor.mtest(df2, conf.level = .95)
colnames(df2)<-c("Uric\nacid","HOMA2\nIR","adipo\nIR", "HOMA2\n%B", "METS-\nVF", "Insulin")
cor2<-cor(df2, method="spearman", use = "complete.obs")
corrplot.mixed(cor2,lower.col = "black", number.cex = 1, upper = "circle",
               p.mat=res2$p, sig.level=.05)

###3.6.3 Fatty Acid Regression analysis
library(lmerTest)

#1
fa.1<-lmer(fa.HOMA2IR~poly(fa.o_METSVF,3)+factor(fa.eth)+(1|fa.dataset))
summary(fa.1)

fa.2<-lmer(fa.adipoIR~poly(fa.o_METSVF,3)+factor(fa.eth)+(1|fa.dataset))
summary(fa.2)

g.fa2<-ggplot(data=data0, aes(x=fa.o_METSVF, y=fa.adipoIR))+
  geom_point(alpha=0.5)+
  stat_smooth(method = "lm", formula = y ~ poly(x,3), size = 1)+
  ylab("Ln adipoIR")+
  xlab("ON METS-VF")+
  theme_classic()
print(g.fa2)

#2
fa.3<-lm(fa.uric_acid~fa.o_METSVF+fa.eth)
summary(fa.3)

fa.4<-lmer(fa.uric_acid~poly(fa.o_METSVF,3)+factor(fa.eth)+(1|fa.dataset))
summary(fa.4)

g.fa4<-ggplot(data0, aes(x=fa.o_METSVF, y=fa.uric_acid))+
  geom_point(alpha=0.5)+
  stat_smooth(method = "lm", formula = y ~ poly(x,3), size = 1)+
  ylab("Ln uric acid")+
  xlab("ON METS-VF")+
  theme_classic()
print(g.fa4)

#3
fa.5<-lmer(fa.HOMA2IR~poly(fa.uric_acid,2)+factor(fa.eth)+fa.age+fa.sex+(1|fa.dataset))
summary(fa.5)

fa.6<-lmer(fa.adipoIR~poly(fa.uric_acid,2)+factor(fa.eth)+fa.age+fa.sex+(1|fa.dataset))
summary(fa.6)

g.fa6<-ggplot(data0, aes(x=fa.uric_acid, y=fa.adipoIR))+
  geom_point(alpha=0.5)+
  stat_smooth(method = "lm", formula = y ~ poly(x,2), size = 1)+
  ylab("Ln adipoIR")+
  xlab("Ln uric acid")+
  theme_classic()
print(g.fa6)

#4

fa.7<-lmer(fa.o_METSVF~poly(fa.HOMA2IR,3)+factor(fa.eth)+(1|fa.dataset))
summary(fa.7)

fa.8<-lmer(fa.o_METSVF~poly(fa.adipoIR,3)+factor(fa.eth)+(1|fa.dataset))
summary(fa.8)

g.fa8<-ggplot(data0, aes(x=fa.adipoIR, y=fa.o_METSVF))+
  geom_point(alpha=0.5)+
  stat_smooth(method = "lm", formula = y ~ poly(x,3), size = 1)+
  ylab("ON METS-VF")+
  xlab("Ln adipoIR")+
  theme_classic()
print(g.fa8)

#5

fa.9<-lmer(fa.uric_acid~fa.HOMA2IR+factor(fa.eth)+fa.age+fa.sex+(1|fa.dataset))
summary(fa.9)

fa.10<-lmer(fa.uric_acid~poly(fa.adipoIR,3)+factor(fa.eth)+fa.age+fa.sex+(1|fa.dataset))
summary(fa.10)

g.fa10<-ggplot(data0, aes(x=fa.adipoIR, y=fa.uric_acid))+
  geom_point(alpha=0.5)+
  stat_smooth(method = "lm", formula = y ~ poly(x,3), size = 1)+
  ylab("Ln uric acid")+
  xlab("Ln adipoIR")+
  theme_classic()
print(g.fa10)

#6

fa.11<-lm(fa.o_METSVF~fa.uric_acid+fa.eth)
summary(fa.11)

fa.12<-lmer(fa.o_METSVF~poly(fa.uric_acid,3)+factor(fa.eth)+(1|fa.dataset))
summary(fa.12)

g.fa12<-ggplot(data0, aes(x=fa.uric_acid, y=fa.o_METSVF))+
  geom_point(alpha=0.5)+
  stat_smooth(method = "lm", formula = y ~ poly(x,3), size = 1)+
  xlab("Ln uric acid")+
  ylab("ON METS-VF")+
  theme_classic()
print(g.fa12)

reg.3<-list(fa.2, fa.4, fa.6, fa.8, fa.10, fa.12)
regresiones.fa<-reg.sum(reg.3)

#7

fa.1.1<-lmer(fa.HOMA2IR~poly(fa.o_METSVF,3)+fa.uric_acid+factor(fa.eth)+(1|fa.dataset))
summary(fa.1.1)

fa.1.2<-lmer(fa.adipoIR~poly(fa.o_METSVF,3)+fa.uric_acid+factor(fa.eth)+(1|fa.dataset))
summary(fa.1.2)

#8

fa.1.3<-lmer(fa.uric_acid~poly(fa.o_METSVF,3)+fa.HOMA2IR+factor(fa.eth)+(1|fa.dataset))
summary(fa.1.3)

fa.1.4<-lmer(fa.uric_acid~poly(fa.o_METSVF,3)+fa.adipoIR+factor(fa.eth)+(1|fa.dataset))
summary(fa.1.4)

#9

fa.1.5<-lmer(fa.adipoIR~poly(fa.uric_acid,2)+fa.o_METSVF+factor(fa.eth)+(1|fa.dataset))
summary(fa.1.5)

fa.1.6<-lmer(fa.adipoIR~poly(fa.uric_acid,2)+fa.o_METSVF+factor(fa.eth)+(1|fa.dataset))
summary(fa.1.6)

#10

fa.1.7<-lmer(fa.o_METSVF~poly(fa.HOMA2IR,3)+factor(fa.eth)+fa.uric_acid+(1|fa.dataset))
summary(fa.1.7)

fa.1.8<-lmer(fa.o_METSVF~poly(fa.adipoIR,3)+factor(fa.eth)+fa.uric_acid+(1|fa.dataset))
summary(fa.1.8)

#11

fa.1.9<-lmer(fa.uric_acid~fa.HOMA2IR+factor(fa.eth)+fa.o_METSVF+(1|fa.dataset))
summary(fa.1.9)

fa.1.10<-lmer(fa.uric_acid~poly(fa.adipoIR,2)+factor(fa.eth)+fa.o_METSVF+(1|fa.dataset))
summary(fa.1.10)

#12

fa.1.11<-lmer(fa.o_METSVF~poly(fa.uric_acid,3)+factor(fa.eth)+fa.HOMA2IR+(1|fa.dataset))
summary(fa.1.11)

fa.1.12<-lmer(fa.o_METSVF~poly(fa.uric_acid,3)+factor(fa.eth)+fa.adipoIR+(1|fa.dataset))
summary(fa.1.12)

r.arrange<-ggarrange(g.fa2, g.fa4, g.fa6, g.fa8, g.fa10, g.fa12, 
                     labels=c("A", "B", "C", "D", "E", "F"), nrow=2, ncol=3)

###3.6.4 Fatty Acid Mediation analysis

detach("package:lmerTest", unload = TRUE)

#C~A - 3
set.seed(123)
ac.fa<-lmer(fa.uric_acid~fa.adipoIR+fa.age+fa.sex+factor(fa.eth)+(1|fa.dataset))
summary(ac.fa)

#B~A+C - 3
bac.fa<-lmer(fa.o_METSVF~fa.adipoIR+fa.uric_acid+factor(fa.eth)+(1|fa.dataset))
summary(bac.fa)

#Mediation - 3
med.fa<-mediate(ac.fa, bac.fa, treat = "fa.adipoIR", 
                mediator = "fa.uric_acid", sims = 1000, boot.ci.type = "perc", boot = F)
summary(med.fa)

#C~A - 4
set.seed(123)
ac.fa1<-lmer(fa.uric_acid~fa.o_METSVF+factor(fa.eth)+(1|fa.dataset))
summary(ac.fa1)

#B~A+C - 4
bac.fa1<-lmer(fa.adipoIR~fa.o_METSVF+fa.uric_acid+factor(fa.eth)+(1|fa.dataset))
summary(bac.fa1)

#Mediation - 4
med.fa1<-mediate(ac.fa1, bac.fa1, treat = "fa.o_METSVF", 
                 mediator = "fa.uric_acid", sims = 1000, boot.ci.type = "perc", boot = F)
summary(med.fa1)

#C~A - 5(HOMA2-IR)
set.seed(123)
ac.fa2<-lmer(fa.uric_acid~fa.HOMA2IR+fa.age+fa.sex+factor(fa.eth)+(1|fa.dataset))
summary(ac.fa2)

#B~A+C -5(HOMA2-IR)
bac.fa2<-lmer(fa.o_METSVF~fa.HOMA2IR+fa.uric_acid+factor(fa.eth)+(1|fa.dataset))
summary(bac.fa2)

#Mediation - 5(HOMA2-IR)
med.fa2<-mediate(ac.fa2, bac.fa2, treat = "fa.HOMA2IR", 
                 mediator = "fa.uric_acid", sims = 1000, boot.ci.type = "perc", boot = F)
summary(med.fa2)

#C~A - 6(HOMA2IR)
set.seed(123)
ac.fa3<-lmer(fa.uric_acid~fa.o_METSVF+factor(fa.eth)+(1|fa.dataset))
summary(ac.fa3)

#B~A+C - 6(HOMA2-IR)
bac.fa3<-lmer(fa.HOMA2IR~fa.o_METSVF+fa.uric_acid+factor(fa.eth)+(1|fa.dataset))
summary(bac.fa3)

#Mediation - 6(HOMA2-IR)
med.fa3<-mediate(ac.fa3, bac.fa3, treat = "fa.o_METSVF", 
                 mediator = "fa.uric_acid", sims = 1000, boot.ci.type = "perc", boot = F)
summary(med.fa3)

###3.6.5 Fatty Acid population table
data0$UAcat<-categ(data0$uric_acid, 5.5)

palmh<-data0%>%filter(UAcat==0)%>%dplyr::select(palmitate)
palmnh<-data0%>%filter(UAcat==1)%>%dplyr::select(palmitate)

t.test(palmh, palmnh)

##3.7 Logistic regressions####
##Logistic regression 1
rlog.1<-glmer(METSVFcat~HOMA2IRcat+factor(ethnicity)+(1|dataset), 
              data=data1, family="binomial")
summary(rlog.1)
exp(coef(rlog.1)$dataset) 

#Logistic regression 2
rlog.2<-glmer(METSVFcat~UAcat+factor(ethnicity)+(1|dataset), 
              data=data1, family="binomial")
summary(rlog.2)
exp(coef(rlog.2)$dataset)

#Logistic regression 3
rlog.3<-glmer(UAcat~HOMA2IRcat+age+sex+factor(ethnicity)+(1|dataset), 
              data=data1, family="binomial")
summary(rlog.3)
exp(coef(rlog.3)$dataset)

#Logistic regression 4
rlog.4<-glmer(UAcat~METSVF.IR+factor(ethnicity)+(1|dataset), 
              data=data1, family="binomial")
summary(rlog.4)
exp(coef(rlog.4)$dataset)

##Logistic regression 1.1: NHANES only
rlog.1.1<-glmer(METSVFcat~HOMA2IRcat+factor(ethnicity)+(1|dataset), 
                data=data0, family="binomial")
summary(rlog.1.1)
exp(coef(rlog.1.1)$dataset)

#Logistic regression 2.1: NHANES only
rlog.2.1<-glmer(METSVFcat~UAcat+factor(ethnicity)+(1|dataset), 
                data=data0, family="binomial")
summary(rlog.2.1)
exp(coef(rlog.2.1)$dataset)

#Logistic regression 3.1: NHANES only
rlog.3.1<-glmer(UAcat~HOMA2IRcat+age+sex+factor(ethnicity)+(1|dataset), 
                data=data0, family="binomial")
summary(rlog.3.1)
exp(coef(rlog.3.1)$dataset)

#Logistic regression 4.1: NHANES only
rlog.4.1<-glmer(UAcat~METSVF.IR+factor(ethnicity)+(1|dataset), 
                data=data0, family="binomial")
summary(rlog.4.1)
exp(coef(rlog.4.1)$dataset)

#Logistic regression 5
rlog.5<-glmer(METSVFcat~fa.adipoIRcat+factor(ethnicity)+(1|dataset), 
              data=data0, family="binomial")
summary(rlog.5)
exp(coef(rlog.5)$dataset) 

#Logistic regression 6
rlog.6<-glmer(METSVFcat~UAcat+factor(ethnicity)+(1|dataset), 
              data=data0, family="binomial")
summary(rlog.6)
exp(coef(rlog.6)$dataset)

#Logistic regression 7
rlog.7<-glmer(UAcat~fa.adipoIRcat+age+sex+factor(ethnicity)+(1|dataset), 
              data=data0, family="binomial")
summary(rlog.7)
exp(coef(rlog.7)$dataset)

#Logistic regression 8
rlog.8<-glmer(UAcat~METSVF.IR+factor(ethnicity)+(1|dataset), data=data0, family="binomial")
summary(rlog.8)
exp(coef(rlog.8)$dataset)

##3.8 Logistic Mediation analysis ####

detach("package:lmerTest", unload = TRUE)

#C~A - 7
set.seed(123)
ac.fa4<-glmer(UAcat~fa.adipoIRcat+age+sex+factor(ethnicity)+(1|dataset), 
              data=data0, family="binomial")
summary(ac.fa4)

#B~A+C - 7
bac.fa4<-glmer(METSVFcat~fa.adipoIRcat+UAcat+factor(ethnicity)+(1|dataset), data=data0,
               family="binomial")
summary(bac.fa4)

#Mediation - 7
med.fa4<-mediate(ac.fa4, bac.fa4, treat = "fa.adipoIRcat", 
                 mediator = "UAcat", sims = 1000, boot.ci.type = "perc", boot = F)
summary(med.fa4)

#C~A - 8
set.seed(123)
ac.fa5<-glmer(UAcat~METSVFcat+age+sex+factor(ethnicity)+(1|dataset),
              data=data0, family="binomial")
summary(ac.fa5)

#B~A+C - 8
bac.fa5<-glmer(fa.adipoIRcat~METSVFcat+UAcat+factor(ethnicity)+(1|dataset),
               data=data0, family="binomial")
summary(bac.fa5)

#Mediation - 8
med.fa5<-mediate(ac.fa5, bac.fa5, treat = "METSVFcat", 
                 mediator = "UAcat", sims = 1000, boot.ci.type = "perc", boot = F)
summary(med.fa5)

#C~A - 9
set.seed(123)
ac.fa6<-glmer(UAcat~HOMA2IRcat+age+sex+factor(ethnicity)+(1|dataset),
              data=data1, family="binomial")
summary(ac.fa6)

#B~A+C - 9
bac.fa6<-glmer(METSVFcat~HOMA2IRcat+UAcat+factor(ethnicity)+(1|dataset), 
               data=data1, family="binomial")
summary(bac.fa6)

#Mediation - 9
med.fa6<-mediate(ac.fa6, bac.fa6, treat = "HOMA2IRcat", 
                 mediator = "UAcat", sims = 1000, boot.ci.type = "perc", boot = F)
summary(med.fa6)

#C~A - 10
set.seed(123)
ac.fa7<-glmer(UAcat~METSVFcat+factor(ethnicity)+(1|dataset), 
              data=data1, family="binomial")
summary(ac.fa7)

#B~A+C - 10
bac.fa7<-glmer(HOMA2IRcat~METSVFcat+UAcat+factor(ethnicity)+(1|dataset),
               data=data1, family="binomial")
summary(bac.fa7)

#Mediation - 10
med.fa7<-mediate(ac.fa7, bac.fa7, treat = "METSVFcat", 
                 mediator = "UAcat", sims = 1000, boot.ci.type = "perc", boot = F)
summary(med.fa7)


##3.9 ROC curve analysis####
roc.1<-roc(HOMA2IRcat~uric_acid, data=data1, ci=TRUE)
plot(roc.1)

m.1 <- optimal.cutpoints(X = "uric_acid", status = "homa2irhealthy", 
                         tag.healthy = "Non resistant", methods = "Youden", 
                         data = data1, pop.prev = NULL, 
                         control = control.cutpoints(), ci.fit = TRUE, 
                         conf.level = 0.95, trace = FALSE, categorical.cov="sex")
summary(m.1)
plot(m.1)

table(data1$homa2irhealthy, data1$hiperuricemia)

roc.3<-roc(UAcat~METSVF, data=data1, ci=TRUE)
plot(roc.3)
m.3 <- optimal.cutpoints(X = "uric_acid", status = "fat", 
                         tag.healthy = "Non obese", methods = "Youden", 
                         data = data1, pop.prev = NULL, 
                         control = control.cutpoints(), ci.fit = TRUE, 
                         conf.level = 0.95, trace = FALSE, categorical.cov = "sex")
summary(m.3)

roc.4<-roc(fa.adipoIRcat~uric_acid, data=data0, ci=TRUE)
plot(roc.4)

data0$adipoirhealthy<-factor(data0$fa.adipoIRcat, labels=c("Non resistant", "Resistant"))

m.4 <- optimal.cutpoints(X = "uric_acid", status = "adipoirhealthy", 
                         tag.healthy = "Non resistant", methods = "Youden", 
                         data = data0, pop.prev = NULL, 
                         control = control.cutpoints(), ci.fit = TRUE, 
                         conf.level = 0.95, trace = FALSE, categorical.cov = "sex")
summary(m.4)

##3.10 DAG####
hir_metsvf_dag <- dagify(fa.METSVF ~ fa.uric_acid,
                         fa.METSVF ~ fa.HOMA2IR + fa.uric_acid,
                         fa.uric_acid ~ fa.HOMA2IR,
                         fa.HOMA2IR ~ fa.adipoIR,
                         fa.adipoIR~fa.uric_acid,
                         fa.METSVF~fa.adipoIR+fa.uric_acid,
                         labels = c("fa.METSVF" = "Visceral\nobesity", 
                                    "fa.HOMA2IR" = "Peripheric\ninsulin\nresistance",
                                    "fa.uric_acid" = "Elevated SUA",
                                    "fa.adipoIR" = "Adipose\ninsulin\nresistance"),
                         exposure = "fa.HOMA2IR",
                         outcome = "fa.METSVF")

dag1<-ggdag(hir_metsvf_dag, text = FALSE, use_labels = "label")+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

print(dag1)

ggsave(dag1, file="Supfig_2.png", height=15, width=20, units= "cm", dpi=300)

##3.11 Boxplots####

#Sex, SUA, BMI
auc3<-ggplot(data1, aes(y=uric_acid, x=fat, fill=factor(BMIcat)))+
  geom_boxplot(outlier.colour="gray", outlier.shape=16,
               outlier.size=1, notch=FALSE)+
  geom_hline(yintercept=5.8, linetype="dashed", color ="red")+
  geom_hline(yintercept=4.6, linetype="dashed", color="red")+
  geom_text(x=1.5, y=5.9, label="M")+
  geom_text(x=1.5, y=4.7, label="F")+
  facet_wrap(~sex1, labeller=label_value)+
  ylab("Uric acid [mg/dL]")+
  xlab("Visceral obesity [METS-VF>7.18]")+
  scale_fill_manual(name="BMI categories",
                    breaks=c("0", "1", "2", "3", "4"),
                    labels=c("Underweight", "Normal", "Overweight",
                             "Obese", "Morbidly obese"),
                    values=c("#0800FF", "#3BC516","#FAEF02","#FF7800","#F90B00"))+
  theme_classic()
print(auc3)

#SUA, BMI, IR
ri.imc<-ggplot(data1, aes(y=uric_acid, x=fat, fill=factor(BMIcat)))+
  geom_boxplot(outlier.colour="gray", outlier.shape=16,
               outlier.size=1, notch=FALSE)+
  geom_hline(yintercept=5.8, linetype="dashed", color ="red")+
  geom_hline(yintercept=4.6, linetype="dashed", color="red")+
  geom_text(x=1.5, y=5.9, label="M")+
  geom_text(x=1.5, y=4.7, label="F")+
  ylab("Uric acid [mg/dL]")+
  xlab("Visceral obesity [METS-VF>7.18]")+
  scale_fill_manual(name="BMI categories",
                    breaks=c("0", "1", "2", "3", "4"),
                    labels=c("Underweight", "Normal", "Overweight",
                             "Obese", "Morbidly obese"),
                    values=c("#0800FF", "#3BC516","#FAEF02","#FF7800","#F90B00"))+
  facet_wrap(~homa2irhealthy, labeller=label_value)+
  theme_classic()
print(ri.imc)

imcgrid<-ggarrange(auc3, ri.imc, common.legend = T, ncol=2, labels=c("A", "B"))

ggsave(imcgrid, file="Figure_1.png", height=15, width=30, units= "cm", dpi=300)

#METSVFIR, SUA
my_comparisons<-list(c("0", "1"), c("0", "2"), c("0", "3"), c("1", "2"), c("1", "3"),
                     c("2", "3"))
auc<-ggplot(data1, aes(y=uric_acid, x=factor(METSVF.IR), fill=factor(METSVF.IR)))+
  geom_boxplot(outlier.colour="gray", outlier.shape=16,
               outlier.size=2, notch=FALSE)+
  geom_hline(yintercept=5.8, linetype="dashed", color ="red")+
  geom_hline(yintercept=4.6, linetype="dashed", color="red")+
  geom_text(x=1.5, y=5.9, label="M")+
  geom_text(x=1.5, y=4.7, label="F")+
  facet_wrap(~sex1, labeller=label_value)+
  ylab("Uric acid [mg/dL]")+
  xlab("Joint variable (METS-VF and HOMA2-IR categorization)")+
  scale_fill_discrete(name="Joint variable\n(METS-VF and HOMA2-IR\ncategorization)",
                      breaks=c("0", "1", "2", "3"),
                      labels=c("Non resistant/non obese", "Resistant/non obese", 
                               "Non resistant/obese", "Resistant/obese"))+
  stat_compare_means(comparisons = my_comparisons) + 
  stat_compare_means(label.y = 1, label.x=1.5)+
  theme_classic()
print(auc)

ggsave(auc, file="Figure_2.png", height=15, width=30, units= "cm", dpi=300)

##3.12 Population tables####

crow_names<-c("Age [years]", "Glucose [mg/dL]", "METSVF", "BMI [kg/m^2]", "cHDL [mg/dL]", 
              "Total cholesterol [mg/dL]", "WHR")
mrow_names<-c("Insulin [uU/mL]", "HOMA2-%B", "HOMA2-IR", "HOMA2-%S", "TGL [mg/dL]")

p.engen<-data1%>%dplyr::select(age, glucose,  METSVF, BMI, cHDL, tcholesterol, WHR)

p.engen1<-cbind(as.matrix(sapply(p.engen, mean, na.rm=TRUE)), 
                as.matrix(sapply(p.engen, sd, na.rm=TRUE)))
p.engen2<-data.frame(p.engen1[,1], p.engen1[,2])

p.engen.med<-data1%>%dplyr::select(insulin, HOMA2B, HOMA2S, HOMA2IR, TGL)
p.engen.med1<-cbind(as.matrix(sapply(p.engen.med, median, na.rm=TRUE)), 
                    as.matrix(sapply(p.engen.med, IQR, na.rm=TRUE)))
p.engen.med2<-data.frame(p.engen.med1[,1], p.engen.med1[,2])

colnames(p.engen2)<-colnames(p.engen.med2)<-c("Mean", "Stdev")

prom.h<-data.frame(data1 %>% filter(UAcat==1))

en.vch<-prom.h %>%dplyr:: select(age, glucose, METSVF, BMI, cHDL, tcholesterol, WHR)

p.enh<-cbind(as.matrix(sapply(en.vch, mean, na.rm=TRUE)), 
             as.matrix(sapply(en.vch, sd, na.rm=TRUE)))
p.enh1<-data.frame(p.enh[,1], p.enh[,2])

p.enh.med<-prom.h%>%dplyr::select(insulin, HOMA2B, HOMA2S, HOMA2IR, TGL)
p.enh.med1<-cbind(as.matrix(sapply(p.enh.med, median, na.rm=TRUE)), 
                  as.matrix(sapply(p.enh.med, IQR, na.rm=TRUE)))
p.enh.med2<-data.frame(p.enh.med1[,1], p.enh.med1[,2])

colnames(p.enh1)<-colnames(p.enh.med2)<-c("Mean_hua", "Stdev_hua")

prom.nh<-data.frame(data1 %>% filter(UAcat==0))

en.vcnh<-prom.nh %>% dplyr:: select(age, glucose, METSVF, BMI, cHDL, tcholesterol, WHR)

p.ennh<-cbind(as.matrix(sapply(en.vcnh, mean, na.rm=TRUE)), 
              as.matrix(sapply(en.vcnh, sd, na.rm=TRUE)))
p.ennh1<-data.frame(p.ennh[,1], p.ennh[,2])

p.ennh.med<-prom.nh%>%dplyr::select(insulin, HOMA2B, HOMA2S, HOMA2IR, TGL)
p.ennh.med1<-cbind(as.matrix(sapply(p.ennh.med, median, na.rm=TRUE)), 
                   as.matrix(sapply(p.ennh.med, IQR, na.rm=TRUE)))
p.ennh.med2<-data.frame(p.ennh.med1[,1], p.ennh.med1[,2])

colnames(p.ennh1)<-colnames(p.ennh.med2)<-c("Mean_nhua", "Stdev_nhua")

rownames(p.engen2)<-rownames(p.enh1)<-rownames(p.ennh1)<-crow_names
rownames(p.engen.med2)<-rownames(p.enh.med2)<-rownames(p.ennh.med2)<-mrow_names

p_vs<-as.vector(mapply(ttest, en.vch, en.vcnh))
rs<-cbind(crow_names, p_vs)
p.envc<-cbind(p.engen2, p.enh1, p.ennh1, p_vs)

p.enh.med$hua<-1
p.ennh.med$hua<-0
ktdf<-rbind(p.enh.med, p.ennh.med)
p_vm<-as.vector(sapply(FUN=ktest, X=ktdf, y=ktdf$hua))
p_vs<-p_vm[-length(p_vm)]
rs<-cbind(mrow_names, p_vs)
p.menvc<-cbind(p.engen.med2, p.enh.med2, p.ennh.med2, p_vs)

p.envc<-rbind(p.envc, p.menvc)

write.csv(p.envc, "nhp_envc.csv") 

en.gen<-NULL
en.gen<-as.matrix(data1 %>% dplyr::select(diabetes, HOMA2IRcat, METSVFcat, sex)
                  %>% sapply(sum, en.gen))

en.vnch<-NULL
en.vnch<-as.matrix(prom.h %>% dplyr::select(diabetes, HOMA2IRcat, METSVFcat, sex) 
                   %>% sapply(sum, en.vnch))

en.vncnh<-NULL
en.vncnh<-as.matrix(prom.nh %>% dplyr::select(diabetes, HOMA2IRcat, METSVFcat, sex)
                    %>% sapply(sum, en.vncnh))

en.chi.d<-NULL
en.chi.d<-as.table(matrix(c(en.vnch[1,1], en.vnch[2,1], en.vnch[3,1], en.vnch[4,1], 
                            en.vncnh[1,1], en.vncnh[2,1], en.vncnh[3,1], en.vncnh[4,1]), 
                          nrow=4))
names.en.chi<-c("Diabetics [n, %]", "IR [HOMA2-IR>2.5]", "Viscerally obese [METS-VF>7.18]", 
                "Female sex")
names.en.chi2<-c("hiperuricemia", "no_hiperuricemia")
rownames(en.chi.d)<-names.en.chi
colnames(en.chi.d)<-names.en.chi2

Ps <- as.vector(apply(en.chi.d,1,Chsq))
r.t <- cbind(rownames(en.chi.d),Ps)

en.chi.gen<-c(en.gen[1,1], en.gen[2,1], en.gen[3,1], en.gen[4,1])

en.chi.d1<-data.frame(en.chi.d[,1], en.chi.d[,2])
en.chi.d1<-cbind(en.chi.d1, en.chi.gen)
colnames(en.chi.d1)<-c(names.en.chi2, "general")
en.chi.p<-en.chi.d1%>%mutate(percenthua=(hiperuricemia/(hiperuricemia+no_hiperuricemia)*100),
                             percentnhua=(no_hiperuricemia/(hiperuricemia+no_hiperuricemia)*100),
                             percentgen=(general)/nrow(data1)*100)
rownames(en.chi.p)<-names.en.chi

en.chi.p<-cbind(en.chi.p, Ps)

en.chi.p<-en.chi.p[c(3,6,1,4,2,5,7)]
colnames(en.chi.p)<-colnames(p.envc)

pop.table<-rbind(p.envc, en.chi.p)

knitr::kable(pop.table, "latex", digits=c(2, 2, 2, 2, 2, 2, 3))

write.csv(en.chi.p, "nh_chi_p.csv") 

#O: Mexican
#1: Mexican American 	
#2: Other Hispanic 		
#3: Non-Hispanic White 	
#4: Non-Hispanic Black 		
#5: Other Race - Including Multi-Racial




#4. ENSANUT cohort####
rm(list = ls(all.names = TRUE))
#Function reload

Chsq <- function(x){
  x <- matrix(x,byrow =TRUE,nrow=1)
  return(chisq.test(x)$p.value)
}

ttest <- function(x,y){
  return(t.test(x,y)$p.value)
}

ktest<-function(x,y){
  return(kruskal.test(x~y)$p.value)
}

categorize<-function(x){
  c1<-c()
  name<-deparse(substitute(x))
  for (value in x){
    if (value==1){
      cat=1
    }else{
      cat=0
    }
    c1<-c(c1, cat)
  }
  return(name<-c1)
}

categorizeinv<-function(x){
  c1<-c()
  name<-deparse(substitute(x))
  for (value in x){
    if (value==0){
      cat=0
    }else{
      cat=1
    }
    c1<-c(c1, cat)
  }
  return(name<-c1)
}

categ<-function(x, lim){
  category<-c()
  for(value in x){
    if(value<lim){
      cat=0
    }else{
      cat=1
    }
    category<-c(category,cat)
  }
  return(category)
}

inv.categ<-function(x, lim){
  category<-c()
  for(value in x){
    if(value>lim){
      cat=0
    }else{
      cat=1
    }
    category<-c(category,cat)
  }
  return(category)
}

cap<-function(x){
  qnt <- quantile(x, probs=c(.25, .75), na.rm = T)
  caps <- quantile(x, probs=c(.05, .95), na.rm = T)
  H <- 1.5 * IQR(x, na.rm = T)
  x[x < (qnt[1] - H)] <- caps[1]
  x[x > (qnt[2] + H)] <- caps[2]
  return(x)
}

##4.1 Database management/variables####
ensanuthecha <- read_excel("ENSANUT.xlsx")
tabvar<-ensanuthecha %>% dplyr::select(id, valor_GLU_SUERO, AC_UR_MG_DL, prom_cintura, 
                                       prom_talla, sexo.x, edad.x, valor_INSULINA, imc2, 
                                       valor_TRIG, valor_COL_HDL, HOMA2B, HOMA2IR, HOMA2S, 
                                       valor_COLEST, valor_HB1AC)
tabvar[, c(1:16)]<-sapply(tabvar[, c(1:16)], as.numeric)

md.pattern(tabvar)
set.seed(123)
imp<-mice(tabvar, m=1, maxit=1)
tabvar1<-complete(imp, "long")
md.pattern(tabvar1)

tabvar1<-subset(tabvar1, tabvar1$AC_UR_MG_DL!="14.9")

tabvar1$METSIR<-log(2*tabvar1$valor_GLU_SUERO+tabvar1$valor_TRIG)*tabvar1$imc2/
  log(tabvar1$valor_COL_HDL)
tabvar1$whr<-(tabvar1$prom_cintura)/as.numeric(tabvar1$prom_talla)
tabvar1$METSVF<-(4.466+0.011*(log(tabvar1$METSIR)^3)+3.239*(log(tabvar1$whr)^3)
                 +0.319*(tabvar1$sexo.x)+0.594*(log(tabvar1$edad.x)))

glucosa<-tabvar1$valor_GLU_SUERO
insulina<-tabvar1$valor_INSULINA
acidourico<-tabvar1$AC_UR_MG_DL
HDL<-tabvar1$valor_COL_HDL
TGL<-tabvar1$valor_TRIG
IMC<-tabvar1$imc2
colesterol<-tabvar1$valor_COLEST
homa2ir<-tabvar1$HOMA2IR
homa2s<-tabvar1$HOMA2S
homa2b<-tabvar1$HOMA2B
METSIR<-tabvar1$METSIR
sexo<-tabvar1$sexo.x
edad<-tabvar1$edad.x
whr<-tabvar1$whr
cintura<-tabvar1$prom_cintura
talla<-tabvar1$prom_talla
METSVF<-tabvar1$METSVF
hb1ac<-tabvar1$valor_HB1AC

log_acidourico<-cap(log(acidourico))
log_homa2ir<-cap(log(homa2ir))
log_homa2s<-cap(log(homa2s))
log_homa2b<-cap(log(homa2b))
y<-yeojohnson(METSVF)
log_mets<-y$x.t
log_insulina<-cap(log(insulina))

tabvar1$homa2ircat<-categ(homa2ir, 2.5)
homa2ircat<-tabvar1$homa2ircat

tabvar1$metsvfcat<-categ(METSVF, 7.18)
metsvfcat<-tabvar1$metsvfcat

tabvar1$aucat<-categ(acidourico, 5.5)
aucat<-tabvar1$aucat

diab1<-categ(hb1ac, 6.5)+categ(glucosa, 126)
tabvar1$diabetes<-categorize(diab1)
diabetes<-tabvar1$diabetes

c1<-c()
for(value in IMC){
  if(value<18.54){
    cat=0
  }else{
    if(value<24.9){
      cat=1
    }else{
      if(value<29.9){
        cat=2
      }else{
        if(value<40){
          cat=3
        }else{
          cat=4
        }
      }
    }
  }
  c1<-c(c1,cat)
}
tabvar1$imccat<-c1
imccat<-tabvar1$imccat

#Joint variable -> 2*Visceral obesity+IR (x3):
###0: Non resistant/non obese
###1: Resistant/non obese
###2: Non resistant/obese
###3: Resistant/obese

metsvfir<-2*metsvfcat+homa2ircat
tabvar1$metsvfir<-metsvfir

##4.2 Variable analysis####

dunn.test(acidourico, metsvfir, kw=TRUE, label=TRUE)
dunn.test(acidourico, imccat, kw=TRUE, label=TRUE)

wilcox.test(acidourico~homa2ircat, data=tabvar1)
wilcox.test(acidourico~metsvfcat, data=tabvar1)

wilcox.test(acidourico~sexo, data=tabvar1)

#Correlation plot

df1<-cbind(log_acidourico, log_homa2s, log_homa2ir, log_homa2b, log_mets, log_insulina)
res1<-cor.mtest(df1, conf.level = .95)
colnames(df1)<-c("Ácido\núrico", "HOMA2-S", "HOMA2-IR", "HOMA2-B", "METS-VF", "Insulina")
cor1<-cor(df1, method="spearman", use = "complete.obs")
corrplot.mixed(cor1,lower.col = "black", number.cex = 1, upper = "circle",
                         p.mat=res1$p, sig.level=.05)

##4.3 Regression analysis####

#1
regr1<-lm(log_mets~log_homa2ir)
summary(regr1)
confint(regr1)

r1<-lm(log_mets~poly(log_homa2ir,2))
summary(r1)
confint(r1)

bptest(regr1, studentize = FALSE)

shapiro.test(regr1$residuals)
ad.test(regr1$residuals)
hist(regr1$residuals)

gvlma(regr1)
gvlma(r1)

BIC(regr1)-BIC(r1)

#2
regr2<-lm(log_mets~log_acidourico)
summary(regr2)
confint(regr2)

r2<-lm(log_mets~poly(log_acidourico,2))
summary(r2)
confint(r2)

bptest(regr2)

shapiro.test(regr2$residuals)
ad.test(regr2$residuals)
hist(regr2$residuals)

BIC(regr2)-BIC(r2)
gvlma(regr2)
gvlma(r2)

#3
regr3<-lm(log_acidourico~log_homa2ir+edad+sexo)
summary(regr3)
confint(regr3)

r3<-lm(log_acidourico~poly(log_homa2ir,2)+edad+sexo)
summary(r3)
confint(r3)

bptest(regr3)

shapiro.test(regr3$residuals)
ad.test(regr3$residuals)
hist(regr3$residuals)

BIC(regr3)-BIC(r3)

#4
regrev1<-lm(log_homa2ir~log_mets)
summary(regrev1)
confint(regrev1)

rev1<-lm(log_homa2ir~poly(log_mets,2))
summary(rev1)
confint(rev1)

bptest(rev1)

shapiro.test(rev1$residuals)
ad.test(rev1$residuals)
hist(rev1$residuals)

BIC(regrev1)-BIC(rev1)

#5
regrev2<-lm(log_homa2ir~log_acidourico+edad+sexo)
summary(regrev2)
confint(regrev2)

rev2<-lm(log_homa2ir~poly(log_acidourico,2)+edad+sexo)
summary(rev2)
confint(rev2)

bptest(regrev2)

shapiro.test(regrev2$residuals)
ad.test(regrev2$residuals)
hist(regrev2$residuals)

BIC(regrev2)-BIC(rev2)

#6
regrev3<-lm(log_acidourico~log_mets)
summary(regrev3)
confint(regrev3)

rev3<-lm(log_acidourico~poly(log_mets,2))
summary(rev3)
confint(rev3)

bptest(regrev3)

shapiro.test(rev3$residuals)
ad.test(rev3$residuals)
hist(rev3$residuals)

BIC(regrev3)-BIC(rev3)

###4.3.1 Multiple regression analysis
#1
regrm<-lm(log_mets~log_homa2ir+log_acidourico*sexo)
summary(regrm)
confint(regrm)

regrm1<-lm(log_mets~poly(log_homa2ir,2)+log_acidourico)
summary(regrm1)
confint(regrm1)

bptest(regrm1)

shapiro.test(regrm1$residuals)
ad.test(regrm1$residuals)
hist(regrm1$residuals)

BIC(regrm)-BIC(regrm1)

#2
regrml2<-lm(log_mets~log_acidourico+log_homa2ir)
summary(regrml2)
confint(regrml2)

regrm2<-lm(log_mets~poly(log_acidourico,2)+log_homa2ir)
summary(regrm2)
confint(regrm2)

bptest(regrml2)

shapiro.test(regrml2$residuals)
ad.test(regrml2$residuals)
hist(regrml2$residuals)

BIC(regrm2)-BIC(regrml2)

#3
regrml3<-lm(log_acidourico~log_homa2ir+log_mets)
summary(regrml3)
confint(regrml3)

regrm3<-lm(log_acidourico~poly(log_homa2ir,2)+log_mets)
summary(regrm3)
confint(regrm3)

bptest(regrml3)

shapiro.test(regrml3$residuals)
ad.test(regrml3$residuals)
hist(regrml3$residuals)

BIC(regrm3)-BIC(regrml3)

#4
regrevm<-lm(log_homa2ir~log_mets+log_acidourico)
summary(regrevm)
confint(regrevm)

regrevm1<-lm(log_homa2ir~poly(log_mets,2)+log_acidourico)
summary(regrevm1)
confint(regrevm1)

bptest(regrevm1)

shapiro.test(regrevm1$residuals)
ad.test(regrevm1$residuals)
hist(regrevm1$residuals)

BIC(regrevm1)-BIC(regrevm)

#5
regrevml2<-lm(log_homa2ir~log_acidourico+log_mets)
summary(regrevml2)
confint(regrevml2)

regrevm2<-lm(log_homa2ir~poly(log_acidourico,2)+log_mets)
summary(regrevm2)
confint(regrevm2)

bptest(regrevml2)

shapiro.test(regrevml2$residuals)
ad.test(regrevml2$residuals)
hist(regrevml2$residuals)

BIC(regrevm2)-BIC(regrevml2)

#6
regrevml3<-lm(log_acidourico~log_mets+log_homa2ir)
summary(regrevml3)
confint(regrevml3)

regrevm3<-lm(log_acidourico~poly(log_mets,2)+log_homa2ir)
summary(regrevm3)
confint(regrevm3)

bptest(regrevm3)

shapiro.test(regrevm3$residuals)
ad.test(regrevm3$residuals)
hist(regrevm3$residuals)

BIC(regrevm3)-BIC(regrevml3)

##4.4 Mediation analysis####
#C~A - 1
set.seed(123)
ac_1<-lm(log_acidourico~log_homa2ir+sexo+edad+diabetes)
summary(ac_1)

#B~A+C - 1
bac_1<-lm(log_mets~log_acidourico+sexo+log_homa2ir+diabetes)
summary(bac_1)

#Mediation - 1
med_1<- mediate(ac_1, bac_1, treat = "log_homa2ir", 
                mediator = "log_acidourico", sims = 1000, boot.ci.type = "perc", boot = T)
summary(med_1)

test.modmed(med_1, covariates.1 = list(diabetes=0),
            covariates.2 = list(diabetes=1), sims = 1000)

#C~A - 1.1
set.seed(123)
ac_1.1<-lm(log_acidourico~log_homa2ir+edad+sexo)
summary(ac_1.1)

#B~A+C - 1.1
bac_1.1<-lm(log_mets~log_homa2ir+log_acidourico)
summary(bac_1.1)

#Mediation 1.1
med_1.1<-mediate(ac_1.1, bac_1.1, treat = "log_homa2ir", 
                 mediator = "log_acidourico", sims = 1000, boot.ci.type = "perc", boot = T)
summary(med_1.1)

#C~A - 2
set.seed(123)
ac_4<-lm(log_acidourico~log_mets+diabetes)
summary(ac_4)

#B~A+C - 2
bac_4<-lm(log_homa2ir~log_mets+log_acidourico+diabetes)
summary(bac_4)

#Mediation - 2
med_4<-mediate(ac_4, bac_4, treat = "log_mets", 
               mediator = "log_acidourico", sims = 1000, boot.ci.type = "perc", boot = T)
summary(med_4)

#C~A - 2.1
set.seed(123)
ac_4.1<-lm(log_acidourico~log_mets)
summary(ac_4.1)

#B~A+C - 2.1
bac_4.1<-lm(log_homa2ir~log_mets+log_acidourico)
summary(bac_4.1)

#Mediation - 2.1
med_4.1<-mediate(ac_4.1, bac_4.1, treat = "log_mets", 
                 mediator = "log_acidourico", sims = 1000, boot.ci.type = "perc", boot = T)
summary(med_4.1)

##4.5 Logistic regressions####
#1
rlog1<-glm(metsvfcat~homa2ircat, data=tabvar2, family="binomial")
summary(rlog1)
confint(rlog1)
exp(coef(rlog1)) 
exp(confint(rlog1))

#2
rlog2<-glm(metsvfcat~aucat, data=tabvar2, family="binomial")
summary(rlog2)
confint(rlog2)
exp(coef(rlog2))
exp(confint(rlog2))

#3
rlog3<-glm(aucat~homa2ircat+edad+sexo, data=tabvar2, family="binomial")
summary(rlog3)
confint(rlog3)
exp(coef(rlog3))
exp(confint(rlog3))

#4
rlog4<-glm(metsvfir~aucat, data=tabvar2, family="binomial")
summary(rlog4)
confint(rlog4)
exp(coef(rlog4))
exp(confint(rlog4))

##4.6 Logistic Mediation analysis####
#C~A - 3
set.seed(123)
aclog<-glm(aucat~homa2ircat+edad+sexo+diabetes, family="binomial")
summary(aclog)

#B~A+C - 3
baclog<-glm(metsvfcat~homa2ircat+aucat+diabetes, family="binomial")
summary(baclog)

#Mediation - 3
med_log<- mediate(aclog, baclog, treat = "homa2ircat", 
                  mediator = "aucat", sims = 1000, boot.ci.type = "perc", boot = T)
summary(med_log)

#C~A - 4
set.seed(123)
aclog2<-glm(aucat~metsvfcat+diabetes, family="binomial")
summary(aclog2)

#B~A+C - 4
baclog2<-glm(homa2ircat~metsvfcat+aucat+diabetes, family="binomial")
summary(baclog2)

#Mediation - 4
med_log2<- mediate(aclog2, baclog2, treat = "metsvfcat", 
                   mediator = "aucat", sims = 1000, boot.ci.type = "perc", boot = T)
summary(med_log2)

#C~A - 3.1
set.seed(123)
aclog.1<-glm(aucat~homa2ircat+edad+sexo, family="binomial")
summary(aclog.1)

#B~A+C - 3.1
baclog.1<-glm(metsvfcat~homa2ircat+aucat, family="binomial")
summary(baclog.1)

#Mediation - 3.1
med_log.1<- mediate(aclog.1, baclog.1, treat = "homa2ircat", 
                    mediator = "aucat", sims = 1000, boot.ci.type = "perc", boot = T)
summary(med_log.1)

#C~A - 4.1
set.seed(123)
aclog2.1<-glm(aucat~metsvfcat, family="binomial")
summary(aclog2.1)

#B~A+C - 4.1
baclog2.1<-glm(homa2ircat~metsvfcat+aucat, family="binomial")
summary(baclog2.1)

#Mediation - 4.1
med_log2.1<- mediate(aclog2.1, baclog2.1, treat = "metsvfcat", 
                     mediator = "aucat", sims = 1000, boot.ci.type = "perc", boot = T)
summary(med_log2.1)

##4.7 ROC curve analysis####
roc1<-roc(homa2ircat~acidourico, data=tabvar2, ci=TRUE)
plot(roc1)
category4<-NULL
for(value in HOMA2IR){
  if(value<2.5){
    cat="Good"
  }else{
    cat="Sick"
  }
  category4<-c(category4,cat)
}
tabvar2$homa2irhealthy<-category4
homa2irhealthy<-tabvar2$homa2irhealthy

m1 <- optimal.cutpoints(X = "acidourico", status = "homa2irhealthy", 
                        tag.healthy = "Good", methods = "Youden", 
                        data = tabvar2, pop.prev = NULL, 
                        control = control.cutpoints(), ci.fit = TRUE, 
                        conf.level = 0.95, trace = FALSE, categorical.cov="sexo1")
summary(m1)
plot(m1)

category5<-NULL
for(value in acidourico){
  if(value<5.5){
    cat="Good"
  }else{
    cat="Sick"
  }
  category5<-c(category5,cat)
}
tabvar2$hiperuricemia<-category5
hiperuricemia<-tabvar2$hiperuricemia
table(tabvar2$homa2irhealthy, hiperuricemia)

caret::confusionMatrix(factor(tabvar2$homa2irhealthy), factor(hiperuricemia))

category7<-NULL
for(value in METSVF){
  if(value<7.18){
    cat="Good"
  }else{
    cat="Sick"
  }
  category7<-c(category7,cat)
}
tabvar2$fat<-category7
fat<-tabvar2$fat
table(tabvar2$auhealthy, fat)

roc3<-roc(aucat~METSVF, data=tabvar2, ci=TRUE)
plot(roc3)
m3 <- optimal.cutpoints(X = "acidourico", status = "fat", 
                        tag.healthy = "Good", methods = "Youden", 
                        data = tabvar2, pop.prev = NULL, 
                        control = control.cutpoints(), ci.fit = TRUE, 
                        conf.level = 0.95, trace = FALSE, categorical.cov="sexo1")
summary(m3)

caret::confusionMatrix(factor(tabvar2$auhealthy), factor(fat))

##4.8 Population tables####
p.engen<-tabvar2%>%dplyr::select(edad, glucosa, insulina, homa2b, sensinsulina, 
                                 HOMA2IR, METSVF, IMC, HDL, TGL, colesterol, whr)

p.engen1<-cbind(as.matrix(sapply(p.engen, mean, na.rm=TRUE)), 
                as.matrix(sapply(p.engen, sd, na.rm=TRUE)))
p.engen2<-data.frame(p.engen1[,1], p.engen1[,2])
colnames(p.engen2)<-c("Mean", "Stdev")

prom.h<-data.frame(tabvar2 %>% filter(aucat==1))

en.vch<-prom.h %>%dplyr:: select(edad, glucosa, insulina, homa2b, sensinsulina, 
                                 HOMA2IR, METSVF, IMC, HDL, TGL, colesterol, whr)

p.enh<-cbind(as.matrix(sapply(en.vch, mean, na.rm=TRUE)), 
             as.matrix(sapply(en.vch, sd, na.rm=TRUE)))
p.enh1<-data.frame(p.enh[,1], p.enh[,2])
colnames(p.enh1)<-c("Media_hua", "Stdev_hua")

prom.nh<-data.frame(tabvar2 %>% filter(aucat==0))

en.vcnh<-prom.nh %>%dplyr:: select(edad, glucosa, insulina, homa2b, sensinsulina, 
                                   HOMA2IR, METSVF, IMC, HDL, TGL, colesterol, whr)

p.ennh<-cbind(as.matrix(sapply(en.vcnh, mean, na.rm=TRUE)), 
              as.matrix(sapply(en.vcnh, sd, na.rm=TRUE)))
p.ennh1<-data.frame(p.ennh[,1], p.ennh[,2])
colnames(p.ennh1)<-c("Media_nhua", "Stdev_nhua")

p_vs <- as.vector(mapply(ttest, en.vch, en.vcnh))
rs <- cbind(rownames(p.ennh1), p_vs)

p.envc<-cbind(p.engen2, p.enh1, p.ennh1, p_vs)
write.csv(p.envc, "p_envc.csv") 

en.vnch<-NULL
en.vnch<-as.matrix(prom.h %>% dplyr::select(diabetes, homa2ircat, metsvfcat, sexo, 
                                            .imp) %>% sapply(sum, en.vnch))
en.vnch<-rbind(en.vnch, (en.vnch[5,1]-en.vnch[4,1]))

en.vncnh<-NULL
en.vncnh<-as.matrix(prom.nh %>% dplyr::select(diabetes, homa2ircat, metsvfcat,
                                              sexo, .imp) %>% sapply(sum, en.vncnh))
en.vncnh<-rbind(en.vncnh, (en.vncnh[5,1]-en.vncnh[4,1]))

en.chi.d<-NULL
en.chi.d<-as.table(matrix(c(en.vnch[1,1], en.vnch[2,1], en.vnch[3,1], en.vnch[4,1], 
                            en.vnch[6,1], en.vncnh[1,1], en.vncnh[2,1], 
                            en.vncnh[3,1], en.vncnh[4,1], en.vncnh[6,1]), nrow=5))
names.en.chi<-c("diabetes", "resistentes h2ir", "obesos mets", "hombres", "mujeres")
names.en.chi2<-c("hiperuricemia", "no_hiperuricemia")
rownames(en.chi.d)<-names.en.chi
colnames(en.chi.d)<-names.en.chi2

Ps <- as.vector(apply(en.chi.d,1,Chsq))
r.t <- cbind(rownames(en.chi.d),Ps)

en.chi.d1<-data.frame(en.chi.d[,1], en.chi.d[,2])
colnames(en.chi.d1)<-names.en.chi2
en.chi.p<-en.chi.d1%>%mutate(percenthua=(hiperuricemia/(no_hiperuricemia+hiperuricemia)*100),
                             percentnhua=(no_hiperuricemia/(no_hiperuricemia+hiperuricemia))*100,
                             general=hiperuricemia+no_hiperuricemia)
rownames(en.chi.p)<-names.en.chi

en.chi.p<-cbind(en.chi.p, Ps)
write.csv(en.chi.p, "en_chi_p.csv") 

#5. SIGMA cohort####
rm(list = ls(all.names = TRUE))
#Function reload

Chsq <- function(x){
  x <- matrix(x,byrow =TRUE,nrow=1)
  return(chisq.test(x)$p.value)
}

ttest <- function(x,y){
  return(t.test(x,y)$p.value)
}

ktest<-function(x,y){
  return(kruskal.test(x~y)$p.value)
}

categorize<-function(x){
  c1<-c()
  name<-deparse(substitute(x))
  for (value in x){
    if (value==1){
      cat=1
    }else{
      cat=0
    }
    c1<-c(c1, cat)
  }
  return(name<-c1)
}

categorizeinv<-function(x){
  c1<-c()
  name<-deparse(substitute(x))
  for (value in x){
    if (value==0){
      cat=0
    }else{
      cat=1
    }
    c1<-c(c1, cat)
  }
  return(name<-c1)
}

categ<-function(x, lim){
  category<-c()
  for(value in x){
    if(value<lim){
      cat=0
    }else{
      cat=1
    }
    category<-c(category,cat)
  }
  return(category)
}

inv.categ<-function(x, lim){
  category<-c()
  for(value in x){
    if(value>lim){
      cat=0
    }else{
      cat=1
    }
    category<-c(category,cat)
  }
  return(category)
}

cap<-function(x){
  qnt <- quantile(x, probs=c(.25, .75), na.rm = T)
  caps <- quantile(x, probs=c(.05, .95), na.rm = T)
  H <- 1.5 * IQR(x, na.rm = T)
  x[x < (qnt[1] - H)] <- caps[1]
  x[x > (qnt[2] + H)] <- caps[2]
  return(x)
}

##5.1 Database management####

data<-read_sav("SIGMA.sav")
homas<-read_excel("SIGMA_HOMA.xls")
data2<-data %>% dplyr::select(sexo, edad, Cintura, promedio_peso, promed_talla, au_b, 
                              masa_VAT, vol_VAT, M_kgFFM, M_kgpeso, METS_VF, METS_IR,
                              adiponectina, promed_imc, ins_b, DM, glu_b, tgb, ct_b, hdl_b, leptina)
data2[, c(1:20)]<-sapply(data2[, c(1:20)], as.numeric)
data2$whr2<-as.numeric((data2$Cintura)/(data2$promed_talla))

md.pattern(data2)
set.seed(123)
imp<-mice(data2, m=1, maxit=1)
data3<-complete(imp, "long")
md.pattern(data3)

homas<-homas %>% dplyr::select(`HOMA2%B`,`HOMA2%S`, HOMA2IR)

data3<-cbind(homas, data3)

data3<-data3%>%rename(homa2b2=`HOMA2%B`,homa2s2=`HOMA2%S`, homa2ir2=HOMA2IR)

md.pattern(data3)
set.seed(123)
imp<-mice(data3, m=1, maxit=1)
data3<-complete(imp, "long")
md.pattern(data3)

adiponectina<-data3$adiponectina
acidourico1<-data3$au_b
masa_vat<-data3$masa_VAT
vol_vat<-data3$vol_VAT
valorM<-data3$M_kgpeso
metsvf<-data3$METS_VF
promimc<-data3$promed_imc
insulina1<-data3$ins_b
glucosa1<-data3$glu_b
Mffm<-data3$M_kgFFM
diabetes1<-data3$DM
edad2<-data3$edad
sexo2<-data3$sexo
tgb<-data3$tgb
chdl<-data3$hdl_b
colt<-data3$ct_b
whr2<-data3$whr2
homa2b2<-data3$homa2b2
homa2ir2<-data3$homa2ir2
homa2s2<-data3$homa2s2

log_adiponectina<-log(adiponectina)
log_acidourico1<-sqrt(acidourico1)
o<-orderNorm(masa_vat)
log_masavat<-o$x.t
log_volvat<-log(vol_vat)
log_valorM<-log(valorM)
log_metsvf<-log(metsvf)
log_insulina<-log(insulina1)
r_Mffm<-sqrt(Mffm)
r_valorM<-sqrt(valorM)

data3$valorMcat<-inv.categ(valorM, 4.7)
valorMcat<-data3$valorMcat

data3$acidouricocat<-categ(acidourico1, 5.5)
acidouricocat<-data3$acidouricocat

data3$h2ir<-categ(homa2ir2, 2.5)
h2ir<-data3$h2ir

data3$vatcat<-categ(masa_vat, 1000)
vatcat<-data3$vatcat

category3<-c()
for(value in promimc){
  if(value<18.54){
    cat=0
  }else{
    if(value<24.9){
      cat=1
    }else{
      if(value<29.9){
        cat=2
      }else{
        if(value<40){
          cat=3
        }else{
          cat=4
        }
      }
    }
  }
  category3<-c(category3,cat)
}
data3$imcat<-category3
imcat<-data3$imcat

data3$catmets<-categ(metsvf, 7.18)
catmets<-data3$catmets

c1<-c()
for (value in diabetes1){
  if (value!=0){
    cat=1
  }else{
    cat=0
  }
  c1<-c(c1, cat)
}

data3$diabetescat<-c1
diabetescat<-data3$diabetescat

##5.2 Correlation plot####
df2<-cbind(log_acidourico1, log_adiponectina, log_masavat, r_valorM, log_metsvf, log_insulina,
           Mffm, data3$homa2ir2)
res2<-cor.mtest(df2, conf.level = .95)
colnames(df2)<-c("Uric\nacid", "Adiponectin", "VAT mass", "M value", "METS-VF", "Insulin", "FFM",
                 "homa2ir")
cor2<-cor(df2, method="spearman", use = "complete.obs")
corrplot.mixed(cor2,lower.col = "black", number.cex = 1, upper = "circle",
               p.mat=res2$p, sig.level=.05)

##5.3 Mediation analyses####
#1
#C~A
polir<-r_valorM^3
set.seed(123)
ac1<-lm(log_acidourico1~polir+edad2+sexo2+diabetes1)
summary(ac1)

#B~A+C
bac1<-lm(log_masavat~polir+log_acidourico1+edad2+sexo2+diabetes1)
summary(bac1)

#Mediation
mediation1<- mediate(ac1, bac1, treat = "polir", 
                     mediator = "log_acidourico1", sims = 1000, boot.ci.type = "perc", 
                     boot = T)
summary(mediation1)

test.modmed(mediation1, covariates.1 = list(diabetes1=0),
            covariates.2 = list(diabetes1=1), sims = 1000)

#C~A*
set.seed(123)
ac1.1<-lm(log_acidourico1~polir+edad2+sexo2)
summary(ac1.1)

#B~A+C*
bac1.1<-lm(log_masavat~polir+log_acidourico1+edad2+sexo2)
summary(bac1.1)

#Mediation*
mediation1.1<- mediate(ac1.1, bac1.1, treat = "polir", 
                       mediator = "log_acidourico1", sims = 1000, boot.ci.type = "perc", 
                       boot = T)
summary(mediation1.1)

#2

#C~A
set.seed(123)
ac2<-lm(log_acidourico1~log_masavat+edad2+sexo2+diabetes1)
summary(ac2)

#B~A+C
bac2<-lm(r_valorM~log_masavat+log_acidourico1+edad2+sexo2+diabetes1)
summary(bac2)

#Mediation
mediation2<- mediate(ac2, bac2, treat = "log_masavat", 
                     mediator = "log_acidourico1", sims = 1000, boot.ci.type = "perc", 
                     boot = T)
summary(mediation2)

test.modmed(mediation2, covariates.1 = list(diabetes1=0),
            covariates.2 = list(diabetes1=1), sims = 1000)

#C~A*
polyvat<-log_masavat^3

set.seed(123)
ac2.1<-lm(log_acidourico1~polyvat+edad2+sexo2)
summary(ac2.1)

#B~A+C*
bac2.1<-lm(r_valorM~polyvat+log_acidourico1+edad2+sexo2)
summary(bac2.1)

#Mediation*
mediation2.1<- mediate(ac2.1, bac2.1, treat = "polyvat", 
                       mediator = "log_acidourico1", sims = 1000, boot.ci.type = "perc", 
                       boot = T)
summary(mediation2.1)

#3

#C~A
set.seed(123)
ac3<-lm(log_adiponectina~log_acidourico1+edad2+sexo2+diabetes1)
summary(ac3)

#B~A+C
bac3<-lm(log_masavat~log_acidourico1+log_adiponectina+edad2+sexo2+diabetes1)
summary(bac3)

#Mediation
mediation3<-mediate(ac3, bac3, treat = "log_acidourico1", 
                    mediator = "log_adiponectina", sims = 1000, boot.ci.type = "perc", 
                    boot = T)
summary(mediation3)

test.modmed(mediation3, covariates.1 = list(diabetes1=0),
            covariates.2 = list(diabetes1=1), sims = 1000)

#C~A*
set.seed(123)
ac3.1<-lm(log_adiponectina~log_acidourico1+edad2+sexo2)
summary(ac3.1)

#B~A+C*
bac3.1<-lm(log_masavat~log_acidourico1+log_adiponectina+edad2+sexo2)
summary(bac3.1)

#Mediation*
mediation3.1<-mediate(ac3.1, bac3.1, treat = "log_acidourico1", 
                      mediator = "log_adiponectina", sims = 1000, boot.ci.type = "perc", 
                      boot = T)
summary(mediation3.1)

#4

#C~A
set.seed(123)
ac4<-lm(log_adiponectina~log_masavat+edad2+sexo2+diabetes1)
summary(ac4)

#B~A+C
bac4<-lm(log_acidourico1~log_masavat+log_adiponectina+edad2+sexo2+diabetes1)
summary(bac4)

#Mediation
mediation4<-mediate(ac4, bac4, treat = "log_masavat", 
                    mediator = "log_adiponectina", sims = 1000, boot.ci.type = "perc", 
                    boot = T)
summary(mediation4)

test.modmed(mediation4, covariates.1 = list(diabetes1=0),
            covariates.2 = list(diabetes1=1), sims = 1000)
#C~A*
set.seed(123)
ac4.1<-lm(log_adiponectina~log_masavat+edad2+sexo2)
summary(ac4.1)

#B~A+C*
bac4.1<-lm(log_acidourico1~log_masavat+log_adiponectina+edad2+sexo2)
summary(bac4.1)

#Mediation*
mediation4.1<-mediate(ac4.1, bac4.1, treat = "log_masavat", 
                      mediator = "log_adiponectina", sims = 1000, boot.ci.type = "perc", 
                      boot = T)
summary(mediation4.1)

#5

#C~A
set.seed(123)
ac5<-lm(log_adiponectina~r_valorM+edad2+sexo2+diabetes1)
summary(ac5)

#B~A+C
bac5<-lm(log_masavat~r_valorM+log_adiponectina+edad2+sexo2+diabetes1)
summary(bac5)

#Mediation
mediation5<-mediate(ac5, bac5, treat = "r_valorM", 
                    mediator = "log_adiponectina", sims = 1000, boot.ci.type = "perc", 
                    boot = T)
summary(mediation5)

test.modmed(mediation5, covariates.1 = list(diabetes1=0),
            covariates.2 = list(diabetes1=1), sims = 1000)

#C~A*
set.seed(123)
ac5.1<-lm(log_adiponectina~r_valorM+edad2+sexo2)
summary(ac5.1)

#B~A+C*
bac5.1<-lm(log_masavat~r_valorM+log_adiponectina+edad2+sexo2)
summary(bac5.1)

#Mediation
mediation5.1<-mediate(ac5.1, bac5.1, treat = "r_valorM", 
                      mediator = "log_adiponectina", sims = 1000, boot.ci.type = "perc", 
                      boot = T)
summary(mediation5.1)

#6

#C~A
set.seed(123)
ac6<-lm(log_adiponectina~log_masavat+sexo2+edad2+diabetes1)
summary(ac6)

#B~A+C
bac6<-lm(r_valorM~log_masavat+log_adiponectina+sexo2+edad2+diabetes1)
summary(bac6)

#Mediation
mediation6<-mediate(ac6, bac6, treat = "log_masavat", 
                    mediator = "log_adiponectina", sims = 1000, boot.ci.type = "perc", 
                    boot = T)
summary(mediation6)

test.modmed(mediation6, covariates.1 = list(diabetes1=0),
            covariates.2 = list(diabetes1=1), sims = 1000)

#C~A*
set.seed(123)
ac6.1<-lm(log_adiponectina~log_masavat+sexo2+edad2)
summary(ac6.1)

#B~A+C*
bac6.1<-lm(r_valorM~log_masavat+log_adiponectina+sexo2+edad2)
summary(bac6.1)

#Mediation
mediation6.1<-mediate(ac6.1, bac6.1, treat = "log_masavat", 
                      mediator = "log_adiponectina", sims = 1000, boot.ci.type = "perc", 
                      boot = T)
summary(mediation6.1)

#7

#C~A
set.seed(123)
ac7<-lm(log_adiponectina~r_valorM+edad2+sexo2+diabetes1)
summary(ac7)

#B~A+C
bac7<-lm(log_acidourico1~r_valorM+log_adiponectina+edad2+sexo2+diabetes1)
summary(bac7)

#Mediation
mediation7<-mediate(ac7, bac7, treat = "r_valorM", 
                    mediator = "log_adiponectina", sims = 1000, boot.ci.type = "perc", 
                    boot = T)
summary(mediation7)

test.modmed(mediation7, covariates.1 = list(diabetes1=0),
            covariates.2 = list(diabetes1=1), sims = 1000)

#C~A
set.seed(123)
ac7.1<-lm(log_adiponectina~r_valorM+edad2+sexo2)
summary(ac7.1)

#B~A+C
bac7.1<-lm(log_acidourico1~r_valorM+log_adiponectina+edad2+sexo2)
summary(bac7.1)

#Mediation
mediation7.1<-mediate(ac7.1, bac7.1, treat = "r_valorM", 
                      mediator = "log_adiponectina", sims = 1000, boot.ci.type = "perc", 
                      boot = T)
summary(mediation7.1)

#8

#C~A
set.seed(123)
ac8<-lm(log_adiponectina~log_acidourico1+edad2+sexo2+diabetes1)
summary(ac8)

#B~A+C
bac8<-lm(r_valorM~log_acidourico1+log_adiponectina+edad2+sexo2+diabetes1)
summary(bac8)

#Mediation
mediation8<-mediate(ac8, bac8, treat = "log_acidourico1", 
                    mediator = "log_adiponectina", sims = 1000, boot.ci.type = "perc", 
                    boot = T)
summary(mediation8)

test.modmed(mediation8, covariates.1 = list(diabetes1=0),
            covariates.2 = list(diabetes1=1), sims = 1000)

#C~A*
set.seed(123)
ac8.1<-lm(log_adiponectina~log_acidourico1+edad2+sexo2)
summary(ac8.1)

#B~A+C*
bac8.1<-lm(r_valorM~log_acidourico1+log_adiponectina+edad2+sexo2)
summary(bac8.1)

#Mediation
mediation8.1<-mediate(ac8.1, bac8.1, treat = "log_acidourico1", 
                      mediator = "log_adiponectina", sims = 1000, boot.ci.type = "perc", 
                      boot = T)
summary(mediation8.1)

#9
mediador<-log_acidourico1/log_adiponectina
log_mediador<-log(mediador)

#C~A
set.seed(123)
ac9<-lm(log_mediador~r_valorM+edad2+sexo2+diabetes1)
summary(ac9)

#B~A+C
bac9<-lm(log_masavat~r_valorM+log_mediador+edad2+sexo2+diabetes1)
summary(bac9)

#Mediation
mediation9<-mediate(ac9, bac9, treat = "r_valorM", 
                    mediator = "log_mediador", sims = 1000, boot.ci.type = "perc", 
                    boot = T)
summary(mediation9)

test.modmed(mediation9, covariates.1 = list(diabetes1=0),
            covariates.2 = list(diabetes1=1), sims = 1000)

#C~A*
set.seed(123)
ac9.1<-lm(log_mediador~r_valorM+edad2+sexo2)
summary(ac9.1)

#B~A+C*
bac9.1<-lm(log_masavat~r_valorM+log_mediador+edad2+sexo2)
summary(bac9.1)

#Mediation
mediation9.1<-mediate(ac9.1, bac9.1, treat = "r_valorM", 
                      mediator = "log_mediador", sims = 1000, boot.ci.type = "perc", 
                      boot = T)
summary(mediation9.1)

#10
#C~A
set.seed(123)
ac10<-lm(log_mediador~log_masavat+edad2+sexo2+diabetes1)
summary(ac10)

#B~A+C
bac10<-lm(r_valorM~log_masavat+log_mediador+edad2+sexo2+diabetes1)
summary(bac10)

#Mediation
mediation10<-mediate(ac10, bac10, treat = "log_masavat", 
                     mediator = "log_mediador", sims = 1000, boot.ci.type = "perc", 
                     boot = T)
summary(mediation10)

test.modmed(mediation10, covariates.1 = list(diabetes1=0),
            covariates.2 = list(diabetes1=1), sims = 1000)

#C~A*
set.seed(123)
ac10.1<-lm(log_mediador~log_masavat+edad2+sexo2)
summary(ac10.1)

#B~A+C*
bac10.1<-lm(r_valorM~log_masavat+log_mediador+edad2+sexo2)
summary(bac10.1)

#Mediation
mediation10.1<-mediate(ac10.1, bac10.1, treat = "log_masavat", 
                       mediator = "log_mediador", sims = 1000, boot.ci.type = "perc", 
                       boot = T)
summary(mediation10.1)

#10 (multimed)

data4<-as.data.frame(cbind(sexo2, edad2, diabetes1, log_masavat, r_valorM, 
                           log_acidourico1, log_adiponectina))
r_multimed<-lm(log_masavat~r_valorM+log_acidourico1+log_adiponectina, data=data4)
summary(r_multimed)

xnames<-c("sexo2","edad2","diabetes1")

mediation11<-multimed(outcome="r_valorM", med.main="log_acidourico1", 
                      med.alt="log_adiponectina", treat="log_masavat", covariates=xnames,
                      sims=100, data=data4)
summary(mediation11)

mediation12<-multimed(outcome="log_masavat", med.main="log_acidourico1", 
                      med.alt="log_adiponectina", treat="r_valorM", covariates=xnames,
                      sims=100, data=data4)
summary(mediation12)


##5.4 Population tables####
row_namesc<-c("Age [years]", "Glucose [mg/dL]", "METSVF", "BMI [kg/m^2]", "cHDL [mg/dL]", 
              "Total cholesterol [mg/dL]", "WHR", "VAT mass[g]")
row_namesm<-c("Insulin [uU/mL]","Mvalue [mg/min/kg]", "HOMA2-IR", "HOMA2-%B", "HOMA2-%S", 
              "TGL [mg/dL]", "Adiponectin [ug/mL]", "Leptin [ug/mL]")

p.gen<-data3%>%dplyr::select(edad, glu_b, METS_VF, promed_imc, hdl_b, ct_b, whr2, masa_VAT)

p.gen1<-cbind(as.matrix(sapply(p.gen, mean, na.rm=TRUE)), 
              as.matrix(sapply(p.gen, sd, na.rm=TRUE)))
p.gen2<-data.frame(p.gen1[,1], p.gen1[,2])

p.gen.med<-data3%>%dplyr::select(ins_b, M_kgpeso, homa2ir2, homa2b2, homa2s2,  tgb, 
                                 adiponectina, leptina)
p.gen.med1<-cbind(as.matrix(sapply(p.gen.med, median, na.rm=TRUE)), 
                  as.matrix(sapply(p.gen.med, IQR, na.rm=TRUE)))
p.gen.med2<-data.frame(p.gen.med1[,1], p.gen.med1[,2])

colnames(p.gen2)<-colnames(p.gen.med2)<-c("Mean", "Stdev")


promediosh<-data.frame(data3 %>% filter(acidouricocat==1))

vcph<-promediosh %>%dplyr:: select(edad, glu_b, METS_VF, promed_imc, hdl_b, ct_b, whr2, masa_VAT)

p.hua<-cbind(as.matrix(sapply(vcph, mean, na.rm=TRUE)), 
             as.matrix(sapply(vcph, sd, na.rm=TRUE)))
p.hua1<-data.frame(p.hua[,1], p.hua[,2])

p.hua.med<-promediosh%>%dplyr::select(ins_b, M_kgpeso, homa2ir2, homa2b2, homa2s2,  tgb, 
                                      adiponectina, leptina)
p.hua.med1<-cbind(as.matrix(sapply(p.hua.med, median, na.rm=TRUE)), 
                  as.matrix(sapply(p.hua.med, IQR, na.rm=TRUE)))
p.hua.med2<-data.frame(p.hua.med1[,1], p.hua.med1[,2])

colnames(p.hua1)<-colnames(p.hua.med2)<-c("Media_hua", "Stdev_hua")


promediosnh<-data.frame(data3 %>% filter(acidouricocat==0))

vcpnh<-promediosnh %>%dplyr:: select(edad, glu_b, METS_VF, promed_imc, hdl_b, ct_b, whr2, masa_VAT)

p.nhua<-cbind(as.matrix(sapply(vcpnh, mean, na.rm=TRUE)), 
              as.matrix(sapply(vcpnh, sd, na.rm=TRUE)))
p.nhua1<-data.frame(p.nhua[,1], p.nhua[,2])

p.nhua.med<-promediosnh%>%dplyr::select(ins_b, M_kgpeso, homa2ir2, homa2b2, homa2s2,  tgb, 
                                        adiponectina, leptina)
p.nhua.med1<-cbind(as.matrix(sapply(p.nhua.med, median, na.rm=TRUE)), 
                   as.matrix(sapply(p.nhua.med, IQR, na.rm=TRUE)))
p.nhua.med2<-data.frame(p.nhua.med1[,1], p.nhua.med1[,2])

colnames(p.nhua1)<-colnames(p.nhua.med2)<-c("Media_nhua", "Stdev_nhua")

rownames(p.hua1)<-rownames(p.nhua1)<-rownames(p.gen2)<-row_namesc
rownames(p.hua.med2)<-rownames(p.nhua.med2)<-rownames(p.gen.med2)<-row_namesm

p_values <- as.vector(mapply(ttest, vcph, vcpnh))
results <- cbind(row_namesc, p_values)
p.vc<-cbind(p.gen2, p.hua1, p.nhua1, p_values)

p.hua.med$hua<-1
p.nhua.med$hua<-0
ktdf<-rbind(p.hua.med, p.nhua.med)
p_vm<-as.vector(sapply(FUN=ktest, X=ktdf, y=ktdf$hua))
p_values<-p_vm[-length(p_vm)]
rs<-cbind(row_namesm, p_values)
p.vm<-cbind(p.gen.med2, p.hua.med2, p.nhua.med2, p_values)

p.envc<-rbind(p.vc, p.vm)

gen<-NULL
gen<-as.matrix((data3 %>% dplyr::select(diabetescat, vatcat, catmets, valorMcat, 
                                        h2ir, sexo) %>%
                  sapply(sum, gen)))

vncph<-NULL
vncph<-as.matrix((promediosh %>% dplyr::select(diabetescat, vatcat, catmets, valorMcat, 
                                               h2ir, sexo) %>%
                    sapply(sum, vncph)))

vncpnh<-NULL
vncpnh<-as.matrix(promediosnh %>% dplyr::select(diabetescat, vatcat, catmets, valorMcat, 
                                                h2ir, sexo) %>%
                    sapply(sum, vncpnh))

chi.d<-NULL
chi.d<-as.table(matrix(c(gen[1,1], gen[2,1], gen[3,1], gen[4,1], gen[5,1], gen[6,1],
                         vncph[1,1],vncph[2,1],vncph[3,1], vncph[4,1], vncph[5,1],vncph[6,1], 
                         vncpnh[1,1],vncpnh[2,1],vncpnh[3,1],vncpnh[4,1],vncpnh[5,1],vncpnh[6,1]),
                       nrow=6))
names.chi<-c("Diabetics", "Visceral obese [DXA1000g]", "Visceral obese [METS-VF7.18]", 
             "IR [Mvalue4.6]",  "IR [HOMA2-IR2.5", "Male sex")
names.chi2<-c("gen","hiperuricemia", "no_hiperuricemia")
rownames(chi.d)<-names.chi
colnames(chi.d)<-names.chi2

P_Values <- as.vector(apply(chi.d,1,Chsq))
result <- cbind(names.chi,P_Values)

chi.d1<-data.frame(chi.d[,1], chi.d[,2], chi.d[,3])
colnames(chi.d1)<-names.chi2
chi.p<-chi.d1%>%mutate(percenthua=(hiperuricemia/(hiperuricemia+no_hiperuricemia)*100),
                       percentnhua=(no_hiperuricemia/(hiperuricemia+no_hiperuricemia)*100),
                       general=gen/nrow(data3)*100)
rownames(chi.p)<-names.chi

chi.p<-cbind(chi.p, P_Values)
chi.p<-chi.p[c(1,6,2,4,3,5,7)]
colnames(chi.p)<-colnames(p.envc)

write.csv(chi.p, "chi_p.csv") 

tabla.sigma<-rbind(p.envc, chi.p)
knitr::kable(tabla.sigma, "latex", digits=c(2, 2, 2, 2, 2, 2, 3))