##AlldataAnalyse: CFA group comparison

###Script collates data from autistic and non-autistic adults, then runs confirmatory factor analyses 
###to assess relationships between tests, and group differences in scores on the tests

library(MVN)
library(psych)
library(rcompanion)
library(lavaan)
library(semPlot)
library(dplyr)
library(eeptools)
library(yarrr)

##Specify directory
dir<-"C:\\Users\\alexa\\OneDrive\\Desktop\\CFA\\"

##==============================================================================================================================================
##=======================================
## COLLATING ALL DATA
##=======================================
 
##Read in processed data for tasks (these files incorporate data from autistic and non-autistic people)

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

##Get Gender variable for non-autistic adults
alldata$Gender<-NA
for(i in 1:120){
  myID<-impdata$ID[i]
  myID<-as.numeric(as.character(myID))
  
  if(myID %in% demographicsNonAut$ID){
    alldata$Gender[i]<-demographicsNonAut$gender[which(demographicsNonAut$ID==myID)]
  }
}

##Get Gender for autistic adults
for(i in 121:194){
  
  myID<-impdata$ID[i]
  myID<-as.numeric(as.character(myID))
  
  if(myID %in% demographicsAut$ID){
    alldata$Gender[i]<-demographicsAut$gender[which(demographicsAut$ID==myID)]
  }
}

#===================================================
# Also get IQ, ADOS and CC-SR data for 
# autistic participants
#===================================================

alldata$IQ<-NA
alldata$ADOSscore<-NA
alldata$ADOSclass<-NA
alldata$CCSRpragRawscore<-NA

for(i in 121:length(alldata$ID)){
  
  myID<-alldata$ID[i]

  if(myID %in% ADOSdata$ID){
    alldata$ADOSscore[i]<-ADOSdata$TOTAL[which(ADOSdata$ID==myID)]
    alldata$ADOSclass[i]<-ADOSdata$Classification[which(ADOSdata$ID==myID)]
  }
  if(myID %in% CCSRdata$ID){
    alldata$CCSRpragRawscore[i]<-CCSRdata$Pragcorrected[which(CCSRdata$ID==myID)]
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
##Note that outliers can be inspected below under exclusions
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
alldata1<-alldata[c(which(alldata$Group=="Non-autistic"),which(alldata$ADOSclass=="autism spectrum"),which(alldata$ADOSclass=="autism")),] # secondary data-set; autistic group only includes
#those meeting ADOS-2 criteria

##=======================================
## MAKE TABLE OF DESCRIPTIVE STATS
## and also store info re exclusions
##=======================================

all.desc.stats<-psych::describeBy(alldata[,c(3,11,4:6,10,7:9)],group=alldata$Group)
#print(all.desc.stats,signif=3)

##Calculate Cohen's d for all variables:
##In columns 1:3, place estimates and CIs for all autistic people compared to non-autistic
##In columns 4:6, place estimates and CIs for autistic people meeting ADOS-2 criteria compared to non-autistic people
cohen.d<-data.frame(rep(NA,8),rep(NA,8),rep(NA,8),rep(NA,8),rep(NA,8),rep(NA,8))
row<-1
for(i in 3:10){
cohen.d[row,1]<-effsize::cohen.d(d=alldata[,i],f=alldata$Group,na.rm=T)$estimate
cohen.d[row,2:3]<-effsize::cohen.d(d=alldata[,i],f=alldata$Group,na.rm=T)$conf.int
cohen.d[row,4]<-effsize::cohen.d(d=alldata1[,i],f=alldata1$Group,na.rm=T)$estimate
cohen.d[row,5:6]<-effsize::cohen.d(d=alldata1[,i],f=alldata1$Group,na.rm=T)$conf.int
row<-row+1
}
row.names(cohen.d)<-colnames(alldata)[3:10]
colnames(cohen.d)<-c("Cohen.d.full.sample","Lower95CI.1","Upper95CI.1","Cohen.d.AutMeetADOS","Lower95CI.1","Upper95CI.1")

exclusions<-data.frame(c("implicature","Awkward Dialogues","Inferencing","Grammar"),rep(NA,4),rep(NA,4))
colnames(exclusions)<-c("Test","Outlying.value","Participant")
#Collect outlying values, and the 
exclusions[1,2]<-paste(impdata$imp_total[which(impdata$outlier==1)],collapse=" ")
exclusions[4,2]<-paste(grammardata$total[which(grammardata$outlier==1)],collapse=" ")
exclusions[3,2]<-paste(inferencingdata$total[which(inf.outliers==1)],collapse=" ")
exclusions[2,2]<-paste(round(awkwarddata$factor_scores[which(awkwarddata$outlier==1)],2),collapse=" ")
exclusions[1,3]<-paste(impdata$ID[which(impdata$outlier==1)],collapse=" ")
exclusions[4,3]<-paste(grammardata$ID[which(grammardata$outlier==1)],collapse=" ")
exclusions[3,3]<-paste(inferencingdata$ID[which(inf.outliers==1)],collapse=" ")
exclusions[2,3]<-paste(awkwarddata$ID[which(awkwarddata$outlier==1)],collapse=" ")

##==============================================================================================================================================
##=======================================
## CONFIRMATORY FACTOR ANALYSES
##=======================================

#========================================
# PREPARE DATA
#========================================

##In selecting estimator for lavaan models, note that MLR estimator 
#allows calculation of robust standard errors and test statistic for incomplete data. 

##Inspect non-transformed correlations
cor.matrix.notrans<-cor(alldata[,3:9],method="pearson",use="pairwise.complete.obs")

##Check normality of variables.
normplots<-mvn(na.omit(alldata[,3:9]), mvnTest = "royston", univariatePlot = "histogram")
normtest.nontrans<-mvn(na.omit(alldata[,3:9]), multivariateOutlierMethod = "adj",showOutliers=T)

##Note that several cases are identified as multivariate outliers.
##With transformations, none are identified, suggesting
##analysis should proceed with transformed data.

##Let's take a look at the outliers in the non-transformed data
mvoutliers.nontrans.find<-as.numeric(normtest.nontrans$multivariateOutliers$Observation)
mvoutliers.nontrans<-alldata[mvoutliers.nontrans.find,]

##Do we want to transform?
transform<-1

if(transform==1){
  ##Add an arbitrary constant to TOM and awkward conversations data, as the
  ##transformation requires there be no negative numbers
  alldata$TOManimations<-alldata$TOManimations+3
  alldata$awkwardconvo<-alldata$awkwardconvo+3
  
  ##Use the Tukey transform -- this establishes the best power to which to raise each variable,
  ##so that it is more normally distributed. 
  for(i in 3:9){
    alldata[,i]<-transformTukey(alldata[,i],plotit=FALSE)
  }

  ##We need to change the scales of some variables so that they are more similar
  alldata[,c(3,4)]<-alldata[,c(3,4)]/100
  alldata[,c(7,9)]<-alldata[,c(7,9)]/100000

  #Now inspect transformed data
  cor.matrix.trans<-cor(alldata[,3:9],method="pearson",use="pairwise.complete.obs")
  normplots.trans<-mvn(na.omit(alldata[,3:9]), mvnTest = "royston", univariatePlot = "histogram")
  normtest.trans<-mvn(na.omit(alldata[,3:9]), multivariateOutlierMethod = "adj",showOutliers=T)
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

##Now run two models

twofit<-cfa(model=twofactor,data=alldata,meanstructure=T,std.lv=T,missing="fiml",estimator="MLR")
two.fit.stats<-as.data.frame(fitMeasures(twofit,c("cfi","rmsea","rmsea.ci.lower","rmsea.ci.upper")))
#There is a warning hear about negative variances (i.e. there is a Heywood case)
#This is an indication that there is model mis-specification
#Inspecting correlations above, looks like this is due to grammar clustering
#with factor ONE

onefit<-cfa(model=onefactor,data=alldata,meanstructure=T,std.lv=T,missing="fiml",estimator="MLR") 
one.fit.stats<-as.data.frame(fitMeasures(onefit,c("cfi","rmsea","rmsea.ci.lower","rmsea.ci.upper")))
#One factor model is weak

#####Define new model without vocab, as this variable seems distinct from other language tests

revised<-'
COMM=~implicature+fillers+grammar+awkwardconvo+TOManimations+inferencing
'
fit.revised<-cfa(model=revised,data=alldata,meanstructure=T,std.lv=T,missing="fiml",estimator="MLR") 
fit.revised.stats<-as.data.frame(fitMeasures(fit.revised,c("cfi","rmsea","rmsea.ci.lower","rmsea.ci.upper")))
##Good fit stats!

##Store factor scores for this model
alldata<-data.frame(alldata,lavPredict(fit.revised))

## HYPOTHESIS 2: Are there significant group differences at the level of the factors?

##Now add Group as a covariate to the model, and compare models where path between Group (i.e. autistic or non-autistic) and the factors
##is allowed to vary and where it isn't

with.Group.effect<-'
ONE=~implicature+fillers+awkwardconvo+TOManimations+grammar+inferencing
ONE~Group
'
no.Group.effect<-'
ONE=~implicature+fillers+awkwardconvo+TOManimations+grammar+inferencing
ONE~0*Group
'

fit.with.Group<-cfa(model=with.Group.effect,data=alldata,meanstructure=T,std.lv=T,missing="fiml",estimator="MLR")
Group.fit.stats<-as.data.frame(fitMeasures(fit.with.Group,c("cfi.robust","rmsea.robust","rmsea.ci.lower.robust","rmsea.ci.upper.robust")))

fit.no.Group<-cfa(model=no.Group.effect,data=alldata,meanstructure=T,std.lv=T,missing="fiml",estimator="MLR")
no.Group.fit.stats<-as.data.frame(fitMeasures(fit.no.Group,c("cfi.robust","rmsea.robust","rmsea.ci.lower.robust","rmsea.ci.upper.robust")))

chidiff.group.effect<-anova(fit.with.Group,fit.no.Group)
##Clearly there are group differences, as model with the group effect fixed
##to zero is very weak. However, model allowing a group difference could have
##better fit .... so investigate whether the residuals for particular test(s) are contributing
##to problems with model fit

resid<-residuals(fit.with.Group,type="cor")$cov
##Looks like we need paths between Group and imp too

with.Group.effectV2<-'
ONE=~implicature+fillers+awkwardconvo+TOManimations+grammar+inferencing
ONE+implicature~Group
'
revised.Group.fit<-cfa(model=with.Group.effectV2,data=alldata,meanstructure=T,std.lv=T,estimator="MLR",missing="fiml")
revised.Group.fit.stats<-as.data.frame(fitMeasures(revised.Group.fit,c("cfi.robust","rmsea.robust","rmsea.ci.lower.robust","rmsea.ci.upper.robust")))
revised.Group.result<-standardizedSolution(revised.Group.fit)

##Test for gender differences (there is a group by gender interaction)
alldata$Gender[alldata$Gender=="Female "]<-"Female"
alldata$Gender[alldata$Gender==""]<-NA
Gender.diffs<-summary(aov(alldata$COMM~alldata$Gender*alldata$Group,data=alldata))

## HYPOTHESIS 3: Are factor scores linked to self-reported and observed difficulties in
## social interaction in autistic individuals?

variability.communication<-'
COMM=~implicature+fillers+awkwardconvo+TOManimations+grammar+inferencing
ADOSscore+CCSRpragRawscore~COMM
'

fit.variability.communication<-cfa(model=variability.communication,data=alldata[alldata$Group=="Autistic",],meanstructure=T,std.lv=T,estimator="MLR",missing="fiml")
variability.communication.fit.stats<-as.data.frame(fitMeasures(fit.variability.communication,c("cfi.robust","rmsea.robust","rmsea.ci.lower.robust","rmsea.ci.upper.robust")))
variability.communication.results<-standardizedSolution(fit.variability.communication)
#YES! Moderate pathways of about .4

##Can also check extent to which links between factor scores and other measures
##remain when correcting for general intelligence -- these regressions allow
##extraction of IQ-corrected residuals
alldata$ADOSresid<-resid(lm(alldata$ADOSscore ~ alldata$IQ, na.action = na.exclude))
alldata$CCSRresid<-resid(lm(alldata$CCSRpragRawscore ~ alldata$IQ, na.action = na.exclude))

#========================================
# MAKE PLOTS
#========================================

##Make path diagram of factor analysis comparing groups
x<-c(-2,-1.5,-1,-0.5,0,0.5,-2,-0.75)
y<-c(0,0,0,0,0,0,0.6,1)
lay<-matrix(c(x,y),ncol=2)
semPlot::semPaths(revised.Group.fit,"std",intercepts=F,layout=lay,nodeLabels = c("imp","flr","awk","TOM","grm","inf","Diag","COMM"))

# Make path diagram of factor analysis exploring variability in autistic group
x<-c(-2,-1.5,-1,-0.5,0,0.5,-1.5,0,-0.75)
y<-c(0,0,0,0,0,0,-2,-2,-1)
lay<-matrix(c(x,y),ncol=2)
semPlot::semPaths(fit.variability.communication,"std",intercepts=F,layout=lay,nodeLabels = c("imp","flr","awk","TOM","grm","inf","ADOS","CCSR","COMM"))

##Now make pirate plot showing factor scores in two groups by Gender

##Standardise factor scores to have SD of 1 across the sample
##to allow easier interpretation of group differences in pirate plot
alldata$COMM<-alldata$COMM/sd(alldata$COMM)
pirateplot(COMM ~ Gender*Group, data=alldata,
           ylab="Receptive communication factor score",
           yaxt.y=seq(-2,3,1),pal="gray",point.o=1,inf.method="ci")

##==============================================================================================================================================

##Finally, check whether the factor analysis measured communication skills
##in a similar way across groups through measurement invariance testing.
##We need to set up models in which we progressing fix parameters. Note that
##scores on the implicature test showed a group difference even accounting
##for factor level differences, so we need to leave the intercept on that
##test to vary.

# Model with no constraints
fit.group.noconstraints<-cfa(model=revised,data=alldata,group="Group",
                                std.lv=T,missing="fiml",estimator="MLR")
fit.stats.no.constraints<-t(data.frame(fitMeasures(fit.group.noconstraints,c("cfi.robust","rmsea.robust"))))

# Model constraining loadings to be the same across groups
fit.group.loadings<-cfa(model=revised,data=alldata,group="Group",group.equal=c("loadings"),
                        std.lv=T,missing="fiml",estimator="MLR")
fit.stats.loadings<-t(data.frame(fitMeasures(fit.group.loadings,c("cfi.robust","rmsea.robust"))))

# Model constraining intercepts to be the same across groups
fit.group.intercepts<-cfa(model=revised,data=alldata,group="Group",group.equal=c("loadings","intercepts"),group.partial=c("implicature~1"),
                       std.lv=T,missing="fiml",estimator="MLR")
fit.stats.intercepts<-t(data.frame(fitMeasures(fit.group.intercepts,c("cfi.robust","rmsea.robust"))))

# Model constraining residuals to be the same across groups
fit.group.residuals<-cfa(model=revised,data=alldata,group="Group",group.equal=c("loadings","intercepts","residuals"),group.partial=c("implicature~1"),
                         std.lv=T,missing="fiml",estimator="MLR")
fit.stats.residuals<-t(data.frame(fitMeasures(fit.group.residuals,c("cfi.robust","rmsea.robust"))))

# Now test for difference between models with and without constraints
anova.test.measInv.metric<-as.data.frame(anova(fit.group.noconstraints,fit.group.loadings))
anova.test.measInv.scalar<-as.data.frame(anova(fit.group.loadings,fit.group.intercepts))
anova.test.measInv.strict<-as.data.frame(anova(fit.group.intercepts,fit.group.residuals))

#summary(fit.group.residuals,fit.measures=T)
##If we look at the most constrained model, we see that fit stats are high, indicating
##that the model fits equally well in both groups. The only parameters left to vary
##are means and variances of the factors. Autistic group has lower mean and bigger variance
