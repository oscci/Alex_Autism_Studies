##AlldataAnalyse: CFA group comparison

###Script collates data from autistic and non-autistic adults, then runs confirmatory factor analyses to test for separable core language and social understanding factors,
###and group differences in scores on these factors

library(psych)
library(rcompanion)
library(MVN)
library(lavaan)
library(semPlot)
library(dplyr)

##Specify directory
dir<-""

##==============================================================================================================================================
##=======================================
## COLLATING ALL DATA
##=======================================

##Read in processed data for tasks (this data incorporates that from autistic and non-autistic people)

impdata<-read.csv(paste0(dir,"Implicature.csv"))
vocabdata<-read.csv(paste0(dir,"Vocab.csv"))
grammardata<-read.csv(paste0(dir,"Grammar.csv"))
inferencingdata<-read.csv(paste0(dir,"Inferencing.csv"))
backchanneldata<-read.csv(paste0(dir,"Backchannel.csv"))
TOMdata<-read.csv(paste0(dir,"TOManimations.csv"))
awkwarddata<-read.csv(paste0(dir,"AwkwardDialogues.csv"))

##Read in further variables
demographicsAut<-read.csv(paste0(dir,"DemographicsAutistic.csv"),stringsAsFactors = F)
demographicsNonAut<-read.csv(paste0(dir,"DemographicsNonautistic.csv"),stringsAsFactors = F)
IQdata<-read.csv(paste0(dir,"IQ.csv")) 
ADOSdata<-read.csv(paste0(dir,"ADOScodes.csv"),stringsAsFactors=F)
CCSRdata<-read.csv(paste0(dir,"CCSR.csv"),stringsAsFactors=F)

##Make full data frame
alldata<-data.frame(impdata$ID,impdata$Group,
                    impdata$imp_total,
                    backchanneldata$total,awkwarddata$factor_scores,TOMdata$factor_scores,
                    grammardata$total,vocabdata$total,
                    inferencingdata$total,TOMdata$raw_CONTtotal,impdata$exp_acc)

colnames(alldata)<-c("ID","Group","implicature","fillers","awkwardconvo","TOManimations","grammar","vocab","inferencing","TOMcontrol","impcontrol")

alldata$gender<-NA

##Get gender variable for non-autistic adults
for(i in 1:120){
  myID<-impdata$ID[i]
  myID<-as.numeric(as.character(myID))
  
  if(myID %in% demographicsNonAut$ID){
    alldata$gender[i]<-demographicsNonAut$gender[which(demographicsNonAut$ID==myID)]
  }
}

##Get gender for autistic adults
for(i in 121:194){
  
  myID<-impdata$ID[i]
  myID<-as.numeric(as.character(myID))
  
  if(myID %in% demographicsAut$ID){
    alldata$gender[i]<-demographicsAut$gender[which(demographicsAut$ID==myID)]
  }
}

#===================================================
# Also get IQ, ADOS and CC-SR data for 
# autistic participants
#===================================================

alldata$IQ<-NA
alldata$ADOSscore<-NA
alldata$ADOSclass<-NA
alldata$ADOSzscore<-NA
alldata$CCSRpragZscore<-NA
alldata$CCSRstructZscore<-NA

for(i in 121:length(alldata$ID)){
  
  myID<-alldata$ID[i]
  
  if(myID %in% ADOSdata$ID){
    alldata$ADOSscore[i]<-ADOSdata$TOTAL[which(ADOSdata$ID==myID)]
    alldata$ADOSclass[i]<-ADOSdata$Classification[which(ADOSdata$ID==myID)]
    alldata$ADOSzscore[i]<-ADOSdata$Zscore[which(ADOSdata$ID==myID)]
  }
  if(myID %in% CCSRdata$ID){
    alldata$CCSRpragZscore[i]<-CCSRdata$PragZscore[which(CCSRdata$ID==myID)]
    alldata$CCSRstructZscore[i]<-CCSRdata$StructZscore[which(CCSRdata$ID==myID)]
  }
  
  if(myID %in% IQdata$ID){
    alldata$IQ[i]<-IQdata$total[which(IQdata$ID==myID)]
  }
}

#====================================
# Carry out exclusions
#====================================

#Check IDs all lined up. They are!
IDs<-data.frame(impdata$ID,backchanneldata$ID,awkwarddata$ID,TOMdata$ID,
                inferencingdata$ID,grammardata$ID,vocabdata$ID)

##Need to check inferencing totals for outliers
inf.desc.stats<-describe(alldata$inferencing,IQR=T,quant=c(.25,.75))
find.outliers<-inf.desc.stats$Q0.25-inf.desc.stats$IQR*2.2
outliers<-which(alldata$inferencing<find.outliers)
inf.outliers<-rep(0,194)
inf.outliers[outliers]<-1

##Recode univariate outliers as NA prior to analysis below
##(vocab, fillers and TOManimations have no outliers)
alldata.outliers<-data.frame(impdata$outlier,grammardata$outlier,awkwarddata$outlier,inf.outliers)
alldata$implicature[which(alldata.outliers$impdata.outlier==1)]<-NA
alldata$impcontrol[which(alldata.outliers$impdata.outlier==1)]<-NA
alldata$grammar[which(alldata.outliers$grammardata.outlier==1)]<-NA
alldata$inferencing[which(alldata.outliers$inf.outliers==1)]<-NA
alldata$awkwardconvo[which(alldata.outliers$awkwarddata.outlier==1)]<-NA

#find autistic participants who did not have clinical diagnoses
exclude<-c(507125,507148,551069)
exclude.indices<-numeric()
for(i in 1:3){
  exclude.indices<-c(exclude.indices,which(alldata$ID==exclude[i]))
}

rawdata<-alldata # keep record of all data, including excluded participants
alldata<-alldata[-exclude.indices,] # main data-set excluding 3 participants; autistic group includes anyone reporting a diagnosis
alldata1<-alldata[-which(alldata$ADOSclass=="non spectrum"),] # secondary data-set; autistic group only includes
#those meeting ADOS-2 criteria

alldata.autism<-alldata[alldata$Group=="Autistic",] #drop control participants for analysis of variability in autistic sample

##=======================================
## MAKE TABLE OF DESCRIPTIVE STATS
## and also store info re exclusions
##=======================================

all.desc.stats<-psych::describeBy(alldata[,c(3,11,4:6,10,7:9)],group=alldata$Group)
desc.stats.autistic<-select(all.desc.stats[[1]],"n","mean","sd","min","max","skew","kurtosis")
desc.stats.autistic<-print(desc.stats.autistic,signif = 3)
desc.stats.nonautistic<-dplyr::select(all.desc.stats[[2]],"n","mean","sd","min","max","skew","kurtosis")
desc.stats.nonautistic<-print(desc.stats.nonautistic,signif = 3)

##Calculate Cohen's d for all variables:
##In columns 1:3, place estimates and CIs for all autistic people compared to non-autistic
##In columns 4:6, place estimates and CIs for autistic people meeting ADOS-2 criteria compared to non-autistic people
cohen.d<-data.frame(rep(NA,7),rep(NA,7),rep(NA,7),rep(NA,7),rep(NA,7),rep(NA,7))
row<-1
for(i in 3:8){
  cohen.d[row,1]<-effsize::cohen.d(d=alldata[,i],f=alldata$Group,na.rm=T)$estimate
  cohen.d[row,3:2]<-effsize::cohen.d(d=alldata[,i],f=alldata$Group,na.rm=T)$conf.int
  cohen.d[row,4]<-effsize::cohen.d(d=alldata1[,i],f=alldata1$Group,na.rm=T)$estimate
  cohen.d[row,5:6]<-effsize::cohen.d(d=alldata1[,i],f=alldata1$Group,na.rm=T)$conf.int
  row<-row+1
}

exclusions<-data.frame(c("implicature","Awkward Dialogues","Inferencing","Grammar"),rep(NA,4),rep(NA,4))
exclusions[1,2]<-paste(impdata$imp_total[which(impdata$outlier==1)],collapse=" ")
exclusions[4,2]<-paste(grammardata$total[which(grammardata$outlier==1)],collapse=" ")
exclusions[3,2]<-paste(inferencingdata$total[which(inf.outliers==1)],collapse=" ")
exclusions[2,2]<-paste(round(awkwarddata$factor_scores[which(awkwarddata$outlier==1)],2),collapse=" ")
exclusions[1,3]<-paste(which(impdata$outlier==1),collapse=" ")
exclusions[4,3]<-paste(which(grammardata$outlier==1),collapse=" ")
exclusions[3,3]<-paste(which(inf.outliers==1),collapse=" ")
exclusions[2,3]<-paste(which(awkwarddata$outlier==1),collapse=" ")

##==============================================================================================================================================
##=======================================
## CONFIRMATORY FACTOR ANALYSES
##=======================================

#========================================
# PREPARE DATA
#========================================

##In selecting estimator for lavaan models, note that MLR estimator 
#allows calculation of robust standard errors and test statistic for incomplete data. 

alldata$Group<-factor(alldata$Group,levels=c("Non-autistic","Autistic"))

##Inspect non-transformed correlations
cor.matrix.notrans<-cor(alldata[,3:8],method="spearman",use="pairwise.complete.obs")

##Check normality of variables.
normplots<-mvn(na.omit(alldata[,3:8]), mvnTest = "royston", univariatePlot = "histogram")
normtest.nontrans<-mvn(na.omit(alldata[,3:8]), multivariateOutlierMethod = "adj",showOutliers=T)

##Note that several cases are identified as multivariate outliers.
##With transformations, none are identified, suggesting
##analysis should proceed with transformed data.

##Let's take a look at the outliers in the non-transformed data
mvoutliers.nontrans.find<-normtest.nontrans$multivariateOutliers$Observation
mvoutliers.nontrans<-na.omit(alldata)[mvoutliers.nontrans.find,]

##Do we want to transform?
transform<-1

if(transform==1){
  ##Add an arbitrary constant to TOM and awkward conversations data, as the
  ##transformation requires there be no negative numbers
  alldata$TOManimations<-alldata$TOManimations+3
  alldata$awkwardconvo<-alldata$awkwardconvo+3
  
  ##Use the Tukey transform -- this establishes the best power to which to raise each variable,
  ##so that it is more normally distributed. 
  for(i in 3:8){
    alldata[,i]<-transformTukey(alldata[,i],plotit=FALSE)
  }
  
  ##We need to change the scales of some variables so that they are more similar
  alldata[,c(3,4)]<-alldata[,c(3,4)]/100
  alldata[,c(7)]<-alldata[,c(7)]/10000
  
  #Now inspect transformed data
  cor.matrix.trans<-cor(alldata[,3:8],method="pearson",use="complete.obs")
  normplots.trans<-mvn(na.omit(alldata[,3:8]), mvnTest = "royston", univariatePlot = "histogram")
  normtest.trans<-mvn(na.omit(alldata[,3:8]), multivariateOutlierMethod = "adj",showOutliers=T)
  mvoutliers.trans.find<-as.numeric(normtest.trans$multivariateOutliers$Observation)
}

#========================================
# TEST HYPOTHESES
#========================================

## HYPOTHESIS 1: Does a two-factor or one-factor model better account for the data across autistic and non-autistic individuals?

##Specify models

twofactor<-'
ONE=~implicature+fillers+awkwardconvo+TOManimations
TWO=~vocab+grammar
'

onefactor<-'
ONE=~implicature+fillers+vocab+grammar+awkwardconvo+TOManimations
'

##Now run two models and perform chi-difference test to assess whether one or two factor model fits better across whole sample

twofit<-cfa(model=twofactor,data=alldata,meanstructure=T,std.lv=T,missing="fiml")
two.fit.stats<-as.data.frame(fitMeasures(twofit,c("cfi","rmsea","rmsea.ci.lower","rmsea.ci.upper")))
twofit.result<-standardizedSolution(twofit)
fscores<-data.frame(lavPredict(twofit)) #extract factor scores for analysis later

onefit<-cfa(model=onefactor,data=alldata,meanstructure=T,std.lv=T,missing="fiml") 
one.fit.stats<-as.data.frame(fitMeasures(onefit,c("cfi","rmsea","rmsea.ci.lower","rmsea.ci.upper")))

chidiff<-anova(onefit,twofit)

## HYPOTHESIS 2: Are there significant group differences at the level of the factors?

##Now add Group as a covariate to the model, and compare models where path between Group (i.e. autistic or non-autistic) and the factors
##is allowed to vary and where it isn't

twofactor.covariate<-'
social=~implicature+fillers+awkwardconvo+TOManimations
core=~vocab+grammar
social+core~Group
'
twofactor.covariate0<-'
social=~implicature+fillers+awkwardconvo+TOManimations
core=~vocab+grammar
social+core~0*Group
'

twofit.covariate<-cfa(model=twofactor.covariate,data=alldata,meanstructure=T,std.lv=T,missing="fiml",estimator="MLR")
two.fit.stats.covariate<-as.data.frame(fitMeasures(twofit.covariate,c("cfi.robust","rmsea.robust","rmsea.ci.lower.robust","rmsea.ci.upper.robust")))

twofit.covariate0<-cfa(model=twofactor.covariate0,data=alldata,meanstructure=T,std.lv=T,missing="fiml",estimator="MLR")
two.fit.stats.covariate0<-as.data.frame(fitMeasures(twofit.covariate0,c("cfi.robust","rmsea.robust","rmsea.ci.lower.robust","rmsea.ci.upper.robust")))

chidiff.covariate<-anova(twofit.covariate,twofit.covariate0)

##Investigate whether the residuals for particular test(s) are contributing
##to problems with model fit
resid.two<-residuals(twofit.covariate,type="cor")$cov

# Looks like we need paths between Group, and imp and vocab

twofactor.2<-'
core=~vocab+grammar
social=~implicature+fillers+awkwardconvo+TOManimations
core+social+implicature+vocab~Group
'
twofit.2<-cfa(model=twofactor.2,data=alldata,meanstructure=T,std.lv=T,estimator="MLR",missing="fiml")
two.fit.stats.2<-as.data.frame(fitMeasures(twofit.2,c("cfi.robust","rmsea.robust","rmsea.ci.lower.robust","rmsea.ci.upper.robust")))
twofit.result.2<-standardizedSolution(twofit.2)

#========================================
# MAKE PLOTS
#========================================

# Make path diagram of factor analysis
x<-c(1.5,2.25,-2.25,-3.5,-1.5,-0.75,0,1.875,-1.25)
y<-c(-1,-1,-1,-0.5,-1,-1,-1,0,0)
lay<-matrix(c(x,y),ncol=2)
semPlot::semPaths(twofit.2,"std",intercepts=F,layout=lay,nodeLabels = c("vocb","grmr","impl","Diag","fillrs","awk","ToM","CLang","Social"))

# Now make scatter plot showing factor scores in two groups

alldata$gender[alldata$gender=="Female "]<-"Female"
for.scatter<-which(alldata$gender=="Male"|alldata$gender=="Female")
pd <- position_dodge(0.5) # move them .05 to the left and right
alldata$gender<-factor(alldata$gender)
plot<-ggplot(fscores[for.scatter,], aes(x=ONE, y=TWO)) + 
  geom_point(aes(shape=alldata$Group[for.scatter],color=alldata$gender[for.scatter]),size=4) + theme_bw() + 
  scale_shape_manual(values=c(23, 19)) + labs(shape="Group") +
  scale_color_manual(name="Gender",breaks = c("Female", "Male"),values=c("orangered","black")) +
  xlab("Social Understanding") + ylab("Core Language")

plot +  theme(text = element_text(size=15),axis.text.x = element_text(size=15), axis.text.y = element_text(size=15))

##==============================================================================================================================================
##========================================
## EXPLORATORY ANALYSIS: 
## DOMAIN-GENERAL INFERENCING FACTOR
##========================================

# This analysis tests whether the group differences we observed are better explained by a specific difficulty with
# social inferencing or a domain-general inferencing difference

#Specify new two-factor model in which "social understanding" tests are replaced by tests without a "theory of mind" component
twofactor.2.inf<-'
core=~vocab+grammar
inf=~implicature+inferencing+TOMcontrol
core+inf+implicature+vocab~Group
'

twofit.2.inf<-cfa(model=twofactor.2.inf,data=alldata,meanstructure=T,std.lv=T,estimator="MLR",missing="fiml")
two.fit.stats.2.inf<-as.data.frame(fitMeasures(twofit.2.inf,c("cfi.robust","rmsea.robust","rmsea.ci.lower.robust","rmsea.ci.upper.robust")))
twofit.result.2.inf<-standardizedSolution(twofit.2.inf)
# This tests the fit of a very similar CFA to that defined above, including group as a covariate.
# Difference: the "social" factor has been redefined to not involve theory of mind

# Now extract factor scores from the social factor and the inferencing factor, and compute the correlation
# between them
twofactor.2.inf.nocov<-'
core=~vocab+grammar
inf=~inferencing+TOMcontrol+implicature
'
twofit.2.inf.nocov<-cfa(model=twofactor.2.inf.nocov,data=alldata,meanstructure=T,std.lv=T,estimator="MLR",missing="fiml")
inf.factorscores<-data.frame(lavPredict(twofit.2.inf.nocov)[,2])
scores.cor<-cor(fscores$ONE,inf.factorscores,use="complete.obs")
social.IQ.cor<-cor(fscores$ONE,alldata$IQ,use="complete.obs")

##==============================================================================================================================================
#========================================
# VARIABILITY IN AUTISTIC GROUP
#========================================

#add factor scores to the data frame that only includes autistic participants
alldata.autism$social<-fscores$ONE[alldata$Group=="Autistic"]
alldata.autism$core<-fscores$TWO[alldata$Group=="Autistic"]

# Convert factor scores to z-scores based on data distribution in the control sample
alldata.autism$social<-(alldata.autism$social-.24)/.7
alldata.autism$core<-(alldata.autism$core-.06)/1.3  

alldata.autism.small<-data.frame(alldata.autism$core,alldata.autism$social,
                                 alldata.autism$ADOSzscore,
                                 alldata.autism$CCSRpragZscore,alldata.autism$CCSRstructZscore,
                                 alldata.autism$IQ)
colnames(alldata.autism.small)<-c("core","social","ADOSzscore","CCSRpragZscore","CCSRstructZscore","IQ")
desc.stats.comm.measures.autistic.group<-psych::describe(alldata.autism.small)
desc.stats.comm.measures.autistic.group<-data.frame(select(desc.stats.comm.measures.autistic.group,"n","mean","sd","min","max","skew","kurtosis"))

# HYPOTHESIS 3
# Run multivariate test: social understanding and core language factor scores, predicting CCSR prag scale + ADOS score
multivar.reg.globalcomm<-lm(cbind(alldata.autism.small$CCSRpragZscore,alldata.autism.small$ADOSzscore) ~ alldata.autism.small$core + alldata.autism.small$social,data=alldata.autism.small)
multivar.reg.globalcomm.results<-heplots::etasq(car::Anova(multivar.reg.globalcomm),anova=T,partial=F)

# Exploratory Analysis


# Test for differences in mean z-scores on CCSR, ADOS and factors
t.test.ADOS.CCSR<-t.test(na.omit(alldata.autism.small$ADOSzscore-alldata.autism.small$CCSRpragZscore),mu=0)
t.test.socialfac.CCSR<-t.test(alldata.autism.small$social,alldata.autism.small$CCSRpragZscore)

alldata.autism$gender[which(alldata.autism$gender=="NB")]<-NA

# Get desc stats in autistic group by gender
desc.stats.gender<-psych::describeBy(cbind(alldata.autism.small),group=alldata.autism$gender)
desc.stats.women<-dplyr::select(desc.stats.gender[[1]],"n","mean","sd")
desc.stats.men<-dplyr::select(desc.stats.gender[[2]],"n","mean","sd")

# Compute Cohen's d for gender differences on communication measures
cohen.d.gender<-data.frame(rep(NA,6),rep(NA,6),rep(NA,6))
row<-1
for(i in 1:6){
  cohen.d.gender[row,1]<-effsize::cohen.d(d=alldata.autism.small[,i],f=alldata.autism$gender,na.rm=T)$estimate
  cohen.d.gender[row,2:3]<-effsize::cohen.d(d=alldata.autism.small[,i],f=alldata.autism$gender,na.rm=T)$conf.int
  row<-row+1
}
desc.stats.gender<-cbind(desc.stats.men,desc.stats.women,cohen.d.gender)

# Plot z-scores for ADOS and CC_SR pragmatic scale against each other; scatter plot with gender coded by colour
plot1<-ggplot(alldata.autism.small, aes(x=ADOSzscore, y=CCSRpragZscore, fill=gender)) + 
  theme_bw() + 
  scale_x_continuous(name="ADOS-2 z-score") + 
  scale_y_continuous(name="CC-SR pragmatic language z-score") +
  annotate("rect", fill = "red", alpha = 0.2, xmin = -Inf, xmax = -2, ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "red", alpha = 0.2, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -2) +
  geom_point(show.legend = T, shape=21, color="black",alpha=0.5,size=5) +
  scale_fill_manual(values=c("yellow", "steelblue", "grey"),labels = c("Female", "Male"))

plot1 + guides(fill = guide_legend(override.aes = list(size=10))) +
  theme(text = element_text(size=15),axis.text.x = element_text(size=15), axis.text.y = element_text(size=15))

# Plot social understanding factor score against CC-SR pragmatic z-score; scatter plot with gender coded by colour
plot2<-ggplot(alldata.autism.small, aes(x=social, y=CCSRpragZscore, fill=gender)) + 
  theme_bw() +
  scale_x_continuous(name="Social understanding z-score") +
  scale_y_continuous(name="CC-SR pragmatic language z-score") +
  annotate("rect", fill = "red", alpha = 0.2, xmin = -Inf, xmax = -1, ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "red", alpha = 0.2, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -2) +
  geom_point(show.legend = T, shape=21, color="black",alpha=0.5,size=5) +  
  scale_fill_manual(values=c("yellow", "steelblue", "grey"),labels = c("Female", "Male", "Non-binary"))

plot2 + guides(fill = guide_legend(override.aes = list(size=10))) +
  theme(text = element_text(size=15),axis.text.x = element_text(size=15), axis.text.y = element_text(size=15))