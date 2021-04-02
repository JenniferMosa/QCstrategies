## Final script QC manuscript
## 25-MAR-2021
## Jennifer Monereo 
## Input: FreeSurfer output before editing (before), FreeSurfer output after editing (after), Demographics (IDS)
plotsfolder<-'C:/Users/Jennifer\ Monereo/Dropbox/Shared/01.2_MCIsample/Manuscripts/Plots/'

## 1. Sample characteristics ################################################################################################################################################
#############################################################################################################################################################################
## Table of demographics and differences between groups (wilcox & chisquared)
variables = c( "Age", "SEX", "bmi", "N_Education_3cat", "MMSEscore", "NDiabetesWHO2")
var_cont = c("Age", "bmi", "MMSEscore" )
var_cat = c("SEX",  "NEducation3cat", "NDiabetesWHO2" )
demographics <- data.frame(Demographic=variables, NonMCI=NA, MCI=NA, WillnPval=NA)

for (j in var_cont) {
  print(j) 
  sum<-ddply(IDS,~CogAdd,summarise,mean=mean(get(j)),sd=sd(get(j)))
  demographics$NonMCI[demographics$Demographic==j] <-paste0(sum[1,2],"(", sum[1,3], ")")
  demographics$MCI[demographics$Demographic==j] <-paste0(sum[2,2],"(", sum[2,3], ")")
  temp<-wilcox.test(get(j)~ IDS$CogAdd, data = IDS)
  demographics$WillnPval[demographics$Demographic==j] <-temp$p.value
  }

for (j in var_cat) {
  print(j)
  sum<-table(IDS$CogAdd, IDS[[j]])
  demographics$NonMCI[demographics$Demographic==j] <-sum[1,2]
  demographics$MCI[demographics$Demographic==j] <-sum[2,2]
  temp<-chisq.test(table(IDS$CogAdd, IDS[[j]]))
  demographics$WillnPval[demographics$Demographic==j] <-temp$p.value
}

write_xlsx(format_headers = F, demographics, paste0(plotsfolder,"demographics.xlsx") ,col_names = T)
rm(list=setdiff(ls(), c("before", "after", "IDS", "plotsfolder")))


## 2. Creation of datasets ##################################################################################################################################################
#############################################################################################################################################################################
## Merge
BEFORE<-cbind(IDS, before)
AFTER<-cbind(IDS, after)

# Change QOala T level tags
BEFORE$Recommendation<-as.numeric(gsub("Exclude", 1, gsub("Include", 0, BEFORE$Recommendation)))
AFTER$Recommendation<-as.numeric(gsub("Exclude", 1, gsub("Include", 0, AFTER$Recommendation)))

## NO QC  ___________________________________________________________________________________________________________________________________________________________________
NoQC<- BEFORE

## Manual QC (+ edits) = Gold standard ______________________________________________________________________________________________________________________________________
GoldStandard <- AFTER[!grepl(3, AFTER$Scoreagreement),]

## Manual QC (+ exclusion)  _________________________________________________________________________________________________________________________________________________
Manual_exclusion <- BEFORE[!grepl(3, AFTER$Scoreagreement),]
Manual_exclusion <- Manual_exclusion[!grepl(2, Manual_exclusion$Scoreagreement),]

bla1<-AFTER[grepl(3, AFTER$Scoreagreement),][1] ## Create vector with excluded ID
bla2<-AFTER[grepl(2, AFTER$Scoreagreement),][1] ## Create vector with edited ID

Excl_visual_V<- rbind(bla1, bla2 ) ## Create vector which the tagged IDs by visual inspection

## Automated MRIQC (+ exclusions) ___________________________________________________________________________________________________________________________________________
excl_MRIQC <- BEFORE[!grepl(1, BEFORE$MRIQCpred),]

excl_MRIQC_V <- BEFORE[grepl(1, BEFORE$MRIQCpred),][1] ## Create vector with excluded ID

## Automated  EN (+ exclusions) _____________________________________________________________________________________________________________________________________________
Q_L <- quantile(BEFORE$FSLH_en, probs=c(.25, .75), na.rm = FALSE)
Q_R <- quantile(BEFORE$FSRH_en, probs=c(.25, .75), na.rm = FALSE)
iqr_L <- IQR(BEFORE$FSLH_en, na.rm = T)
iqr_R <- IQR(BEFORE$FSRH_en, na.rm = T)
low_L<- Q_L[1]-1.5*iqr_L # Left lower Range
low_R<- Q_R[1]-1.5*iqr_R # Right lower Range

excl_EN<- subset(BEFORE, BEFORE$FSLH_en > low_L & BEFORE$FSRH_en > low_R) # Only exclude lower range

excl_EN_V<- subset(BEFORE, BEFORE$FSLH_en < low_L | BEFORE$FSRH_en < low_R)[1] ## Create vector with excluded ID

## Automated CNR (+ exclusions) _____________________________________________________________________________________________________________________________________________
Q_L <- quantile(BEFORE$FSLH_cnr, probs=c(.25, .75), na.rm = FALSE)
Q_R <- quantile(BEFORE$FSRH_cnr, probs=c(.25, .75), na.rm = FALSE)
iqr_L <- IQR(BEFORE$FSLH_cnr, na.rm = T)
iqr_R <- IQR(BEFORE$FSRH_cnr, na.rm = T)
low_L<- Q_L[1]-1.5*iqr_L # Lower Range
low_R<- Q_R[1]-1.5*iqr_R # Lower Range

excl_CNR<- subset(BEFORE, BEFORE$FSLH_cnr > low_L & BEFORE$FSRH_cnr > low_R) # Only exclude lower range

excl_CNR_V<- subset(BEFORE, BEFORE$FSLH_cnr < low_L | BEFORE$FSRH_cnr < low_R)[1] ## Create vector with excluded ID


## Automated  Morpho outliers GLOBAL (+ exclusions) _________________________________________________________________________________________________________________________
names <- data.frame(df=c("subcort","hippo","thick","area"),
                    exclude1=c("EstimatedTotalIntraCranialVol","LH_Wholehippocampus","LH_MeanThickness_thickness", "LH_WhiteSurfArea_area"),
                    exclude2=c("MaskVol","RH_Wholehippocampus","LH_MeanThickness_thickness", "RH_WhiteSurfArea_area"))
for(i in 1:4){
  excl_MORPHOglobal <- BEFORE
  Regions <- c(names$exclude1[i], names$exclude2[i])
  for(Region in Regions){
    Q <- quantile(BEFORE[,Region], probs=c(.25, .75), na.rm = FALSE)
    iqr <- IQR(BEFORE[,Region], na.rm = T)
    up <-  Q[1]+1.5*iqr # Upper Range  
    low<- Q[1]-1.5*iqr # Lower Range
    Outliers <- subset(BEFORE, BEFORE[,Region] < low | BEFORE[,Region] > up)$RandomID 
    assign(paste0("excl_MG_V_", names$df[i]), Outliers) ## Create vector with excluded ID
    excl_MORPHOglobal[excl_MORPHOglobal$RandomID %in% Outliers,Region] <- NA
    excl_MORPHOglobal <-na.omit(excl_MORPHOglobal)
  }
  assign(paste0("excl_MG_", names$df[i]), excl_MORPHOglobal)
}

## Automated Qoala-T (+ exclusions) _________________________________________________________________________________________________________________________________________
excl_Qoala <- BEFORE[!grepl(1, BEFORE$Recommendation),]
excl_Qoala_V <- BEFORE[grepl(1, BEFORE$Recommendation),][1] ## Create vector with excluded ID

## Semi automated MRIQC (+ edits) ___________________________________________________________________________________________________________________________________________
MRIQCbef <- BEFORE[grepl(0, BEFORE$MRIQCpred),] # take from the BEFORE df all the MRIQC scores of 0
MRIQCaft <- AFTER[grepl(1, AFTER$MRIQCpred),] # take from the AFTER df all the MRIQC scores of 1
MRIQCaft <- MRIQCaft[!grepl(3, MRIQCaft$Scoreagreement),] ## Remove manual exclusions
semi_MRIQC <- rbind(MRIQCbef, MRIQCaft) ## Merge
rm(MRIQCbef, MRIQCaft)

## Semi automated EN (+ edits) ______________________________________________________________________________________________________________________________________________
Q_L <- quantile(BEFORE$FSLH_en, probs=c(.25, .75), na.rm = FALSE)
Q_R <- quantile(BEFORE$FSRH_en, probs=c(.25, .75), na.rm = FALSE)
iqr_L <- IQR(BEFORE$FSLH_en, na.rm = T)
iqr_R <- IQR(BEFORE$FSRH_en, na.rm = T)
up_L <-  Q_L[1]+1.5*iqr_L # Upper Range  
up_R <-  Q_R[1]+1.5*iqr_R # Upper Range  
low_L<- Q_L[1]-1.5*iqr_L # Lower Range
low_R<- Q_R[1]-1.5*iqr_R # Lower Range
enbef<- subset(BEFORE, BEFORE$FSLH_en > low_L & BEFORE$FSRH_en > low_R)## Lower
enaft<- subset(AFTER, AFTER$FSLH_en <= low_L | AFTER$FSRH_en <= low_R)## Upper
count(enaft$Scoreagreement)
enaft <- enaft[!grepl(3, enaft$Scoreagreement),] ## Remove manual exclusions
semi_EN <- rbind(enbef, enaft) ## Merge
semi_EN <-na.omit(semi_EN)
rm(enbef, enaft)

## Semi automated  CNR (+ edits) ____________________________________________________________________________________________________________________________________________
Q_L <- quantile(BEFORE$FSLH_cnr, probs=c(.25, .75), na.rm = FALSE)
Q_R <- quantile(BEFORE$FSRH_cnr, probs=c(.25, .75), na.rm = FALSE)
iqr_L <- IQR(BEFORE$FSLH_cnr, na.rm = T)
iqr_R <- IQR(BEFORE$FSRH_cnr, na.rm = T)
up_L <-  Q_L[1]+1.5*iqr_L # Upper Range  
up_R <-  Q_R[1]+1.5*iqr_R # Upper Range  
low_L<- Q_L[1]-1.5*iqr_L # Lower Range
low_R<- Q_R[1]-1.5*iqr_R # Lower Range
cnrbef1<- subset(BEFORE, BEFORE$FSLH_cnr > low_L & BEFORE$FSRH_cnr > low_R) # Lower
cnraft1<- subset(AFTER, AFTER$FSLH_cnr < low_L | AFTER$FSRH_cnr < low_R) # upper
count(cnraft1$Scoreagreement)
cnraft1 <- cnraft1[!grepl(3, cnraft1$Scoreagreement),] ## Remove manual exclusions
semi_CNR <- rbind(cnrbef1, cnraft1) ## Merge
semi_CNR <-na.omit(semi_CNR)
rm(cnrbef1, cnraft1)

## Semi automated Morpho outliers GLOBAL (+ edits) __________________________________________________________________________________________________________________________
## SC
remaining<-AFTER[ !(AFTER$RandomID %in% excl_MG_subcort$RandomID), ] ## Take cases based in the IDs missing from the exclusion based on morphological
count(remaining$Scoreagreement)
remaining <- remaining[!grepl(3, remaining$Scoreagreement),] ## Remove manual exclusions
semi_MG_subcort <- rbind(excl_MG_subcort, remaining) ## Merge
semi_MG_subcort <-na.omit(semi_MG_subcort)
one<-remaining
## HC
remaining<-AFTER[ !(AFTER$RandomID %in% excl_MG_hippo$RandomID), ]
count(remaining$Scoreagreement)
remaining <- remaining[!grepl(3, remaining$Scoreagreement),] ## Remove manual exclusions
semi_MG_hippo <- rbind(excl_MG_hippo, remaining) ## Merge
semi_MG_hippo <-na.omit(semi_MG_hippo)
two<-remaining

## TH
remaining<-AFTER[ !(AFTER$RandomID %in% excl_MG_thick$RandomID), ]
count(remaining$Scoreagreement)
remaining <- remaining[!grepl(3, remaining$Scoreagreement),] ## Remove manual exclusions
semi_MG_thick <- rbind(excl_MG_thick, remaining) ## Merge
semi_MG_thick <-na.omit(semi_MG_thick)
three<-remaining

## AR
remaining<-AFTER[ !(AFTER$RandomID %in% excl_MG_area$RandomID), ]
count(remaining$Scoreagreement)
remaining <- remaining[!grepl(3, remaining$Scoreagreement),] ## Remove manual exclusions
semi_MG_area <- rbind(excl_MG_area, remaining) ## Merge
semi_MG_area <-na.omit(semi_MG_area)
four<-remaining

bla<-rbind(one, two, three, four)
bla<-unique(bla)
count(bla$Scoreagreement)

## Semi automated Qoala (+ edits) ___________________________________________________________________________________________________________________________________________
Qooala.bef <- BEFORE[grepl(0, BEFORE$Recommendation),] # take from the BEFORE df all the MRIQC scores of 0
Qoala.aft <- AFTER[grepl(1, AFTER$Recommendation),] # take from the AFTER df all the MRIQC scores of 1
count(Qoala.aft$Scoreagreement)
Qoala.aft <- Qoala.aft[!grepl(3, Qoala.aft$Scoreagreement),] # take from the AFTER df all the MRIQC scores of 1
semi_Qoala <- rbind(Qooala.bef, Qoala.aft) ## Merge
semi_Qoala <-na.omit(semi_Qoala)
rm(Qooala.bef, Qoala.aft)



##  3. Explore the flagged cases by each QC strategy ########################################################################################################################
#############################################################################################################################################################################

## Merge the morphological vectors
excl_MG_V<-as.data.frame(c(excl_MG_V_area, excl_MG_V_hippo, excl_MG_V_subcort, excl_MG_V_thick))
excl_MG_V<-as.data.frame(excl_MG_V[,1][!duplicated(excl_MG_V[,1])])

## Add dummies
dataframes = c("excl_MG_V", "excl_CNR_V", "excl_EN_V", "excl_MRIQC_V", "excl_Qoala_V", "Excl_visual_V")
for (dataframe in dataframes) {
  temp<-cbind(get(dataframe), 1)
  colnames(temp)<- c("RandomID", dataframe)
  assign(dataframe, temp)
}

# Merge dummies
list<-BEFORE[1]
bla<-full_join(list, excl_MG_V, by = "RandomID")
bla0<-full_join(bla, excl_CNR_V, by = "RandomID")
bla1<-full_join(bla0, excl_EN_V, by = "RandomID")
bla2<-full_join(bla1, excl_MRIQC_V, by = "RandomID")
bla3<-full_join(bla2, excl_Qoala_V, by = "RandomID")
flagged_DF<-full_join(bla3, Excl_visual_V, by = "RandomID")

# Code NA as 0
flagged_DF[is.na(flagged_DF)]=0


## Are more excluded cases in the MCI group?  _______________________________________________________________________________________________________________________________
## Add MCI status
dataframe<-left_join(flagged_DF, IDS[c(1, 7, 8, 13)], by = "RandomID")

vector = c("excl_MG_V", "excl_CNR_V", "excl_EN_V", "excl_MRIQC_V", "excl_Qoala_V", "Excl_visual_V")
comparisons <- data.frame(QC=vector, MCIyes=NA, MCIno=NA, pval = NA)

for (i in vector) {
  temp<-chisq.test(dataframe$CogAdd, dataframe[[i]])
  tab<-print(table(dataframe$CogAdd, dataframe[[i]]))
  comparisons$pval[comparisons$QC == i]<-temp[3]
  comparisons$MCIyes[comparisons$QC == i]<-tab[2,2]
  comparisons$MCIno[comparisons$QC == i]<-tab[1,2]
  assign("comparisons_MCI", comparisons)
}

write_xlsx(comparisons_MCI,paste0(plotsfolder, "comparisons_MCI.xlsx"), col_names = T)

## Do flagged IDS overlap between QC strategies? ____________________________________________________________________________________________________________________________

## Remove ID and add 0 when NA
flagged_DF<-flagged_DF[,2:ncol(flagged_DF)]

## Table for number of overlapping
empty_df = flagged_DF[FALSE,]
for (c in 1:6) {
  for (r in 1:6) {
    empty_df[r,c]<-sum(ifelse(flagged_DF[c] ==1 & flagged_DF[r] ==1, 1, 0))
  }}
rownames(empty_df)<- colnames(empty_df)
overlap<-empty_df
overlap$names<-rownames(overlap)
write_xlsx(overlap,paste0(plotsfolder,"overlaptable.xlsx"), col_names = T)


## Plot an overly complicated venn diagram

png("Venn.png",width = 9, height = 9, units = 'in', res = 300,bg = "transparent")
venn(flagged_DF, ilabels = TRUE, zcolor = "style", small=0.7)
dev.off()

## Table for dice coefficient
dice_df = flagged_DF[FALSE,]
for (c in 1:6) {
  for (r in 1:6) {
    dice_df[r,c]<-(2*sum(ifelse(flagged_DF[c] ==1 & flagged_DF[r] ==1, 1, 0)))/((sum(ifelse(flagged_DF[c] ==1 , 1, 0)))+(sum(ifelse(flagged_DF[r] ==1 , 1, 0))))
  }}

rownames(dice_df)<- colnames(dice_df)
dice_df$names<-rownames(dice_df)
write_xlsx(dice_df, paste0(plotsfolder,"dicetab.xlsx"), col_names = T)

## Clean stuff
rm(list=setdiff(ls(), c("before", "after", "IDS", "plotsfolder", 
                        "excl_CNR", "excl_EN", "excl_MG_subcort", "excl_MG_hippo", "excl_MG_thick", "excl_MG_area", "excl_MRIQC", "excl_Qoala",
                        "semi_CNR", "semi_EN", "semi_MG_subcort", "semi_MG_hippo", "semi_MG_thick", "semi_MG_area", "semi_MRIQC", "semi_Qoala",
                        "NoQC", "Manual_exclusion", "GoldStandard")))



## 4. Manual editing effects ################################################################################################################################################
#############################################################################################################################################################################

## Calculate percentage of change __________________________________________________________________________________________________________________________________________
percentages<-((after-before)/before)*100

## Create vectors for edited people
rowsvector<-rowSums(percentages, na.rm = T)
percentages<-cbind(rowsvector, IDS, percentages)

## Remove all participants with no editions (i.e. rowsvector = to 0)
edited<-percentages[(!percentages$rowsvector==0),]
editedPercent<-edited[c(54:ncol(edited))] ## Preserve only the aseg/aparc output
rm(percentages)

# Make the average change per structure
MeanPercent<-as.data.frame(sapply(editedPercent, mean, na.rm = T))    

MeanPercent<-setDT(MeanPercent, keep.rownames = TRUE)[] ## Make row name into column
names(MeanPercent)[names(MeanPercent) == "rn"] <- "region" ## Change colname
names(MeanPercent)[names(MeanPercent) == "sapply(editedPercentabs, mean, na.rm = T)"] <- "PercentChange" ## Change colname

## What are the structures that change the most?
MaxAndMin <- data.frame(region=colnames(editedPercent[c(1:ncol(editedPercent))]), max=NA, min=NA)

for (i in c(1:ncol(editedPercent))){
  temp<- max(editedPercent[i])
  MaxAndMin$max[i] <- temp
  temp<- min(editedPercent[i])
  MaxAndMin$min[i] <- temp
}

## Wilcoxon signed rank test on paired samples & Effect Size ________________________________________________________________________________________________________________

before.ch<-cbind(IDS$RandomID, before) ## Add ID to before
names(before.ch)[names(before.ch) == "IDS$RandomID"] <- "RandomID" ## Change the ID col name

after.ch<-cbind(IDS$RandomID, after) ## Add ID to after
names(after.ch)[names(after.ch) == "IDS$RandomID"] <- "RandomID" ## Change the ID colname

before.ch<- before.ch[before.ch$RandomID %in% edited$RandomID,] ## Reduce the DF to only those that are present in the edited dataframe
after.ch<- after.ch[after.ch$RandomID %in% edited$RandomID,]

before.ch$group <- c("before") ## add group
after.ch$group <- c("after")

grouped.ch<-rbind(before.ch, after.ch) ## merge before and after into a grouped DF

grouped.ch$group <- as.factor(grouped.ch$group) ## Make grouped$group a factor

## Paired wilcoxon rank test, FDR corrected
will <- data.frame(region=colnames(grouped.ch[c(20:217)]), pval=NA)
for (i in c(20:217)){
  temp<- pairwise_wilcox_test(data = grouped.ch, reformulate(response=colnames(grouped.ch)[i],termlabels =  "group"), paired = TRUE, p.adjust.method = "fdr")
  will$pval[i-19] <- temp$p
}

## Effect size 
effsIZ <- data.frame(region=colnames(grouped.ch[c(20:217)]), r=NA)
for (i in c(20:217)){
  temp <- wilcox_effsize(data = grouped.ch, reformulate(response=colnames(grouped.ch)[i], termlabels =  "group"), paired = TRUE, alternative = "two.sided")
  effsIZ$r[i-19] <- temp$effsize
}

## Merge
res.will.ef<- merge(will, effsIZ, by = "region")
All<-merge(MeanPercent, res.will.ef, by = "region")
names(All)[names(All) == "region"] <- "REGION"

# Add significancestars
All$sig <- ifelse(All$pval>.05, NA, 
                  ifelse(All$pval<=0.05 & All$pval>0.01, "*",
                         ifelse(All$pval<=0.01 & All$pval>0.001, "**",
                                ifelse(All$pval<=0.001, "***", 0))))   

names(All)[names(All) == "sapply(editedPercent, mean, na.rm = T)"] <- "PercChange" ## Change colname

## Separate by type of region & make sure the region names are GGSEG friendly
SUBC <- All %>% select(-contains("Estimated"))
R2_SC <- SUBC[-grep("thickness", SUBC$REGION), ]
R2_SC <- R2_SC[-grep("area", R2_SC$REGION), ]
R2_SC <- R2_SC[-grep("ippocampaltai", R2_SC$REGION), ]
R2_SC <- R2_SC[-grep("molecularlayer", R2_SC$REGION), ]
R2_SC <- R2_SC[-grep("_fimbria", R2_SC$REGION), ]
R2_SC <- R2_SC[-grep("subiculum", R2_SC$REGION), ]
R2_SC <- R2_SC[-grep(	"_CA", R2_SC$REGION), ]
R2_SC <- R2_SC[-grep(	"GCMLDG", R2_SC$REGION), ]
R2_SC <- R2_SC[-grep(	"_HATA", R2_SC$REGION), ]
R2_SC <- R2_SC[-grep(	"hippocampalfissure", R2_SC$REGION), ]

AR_sig <- All[grep("_area", All$REGION), ]

TH_sig <- All[grep("_thickness", All$REGION), ]

## Sup fig legends _______________________________________________________________________________________________________________________________________________________
ggseg(colour="black", size=.7, mapping=aes(fill=region)) +
  theme_void()+
  theme(legend.position = "bottom")+
  scale_color_viridis_d(na.value="light gray")

ggseg(atlas="aseg", colour="black", size=.7,mapping=aes(fill=region))+
  theme_void()+
  theme(legend.position = "bottom")+
  scale_fill_manual(values=as.vector(polychrome()), na.value="light gray")

## Figure 4 ______________________________________________________________________________________________________________________________________________________________

## %
ggseg(.data= R2_SC, atlas = "aseg",colour = "gray", mapping=aes(fill= PercChange)) +
  labs(title=" Percentage of change, subcortical volumes", fill="%")+
  scale_fill_gradient2(high="#E69F00",low="#56B4E9",na.value="light gray", limits=c(-7.3, 9), breaks=seq(-7.3, 9,by=5))

ggseg(.data=TH_sig, 
      mapping=aes(fill=PercChange), position="stacked", colour="gray") +
  labs(title="Percentage of change, cortical thickness", fill="%") +
  scale_fill_gradient2(high="#E69F00",low="#56B4E9",na.value="light gray", limits=c(-2.5, 3.5), breaks=seq(-2, 3,by=1))

ggseg(.data=AR_sig, colour="gray", 
      mapping=aes(fill=PercChange), position="stacked") +
  labs(title="Percentage of change, cortical area", fill="%") +
  scale_fill_gradient2(high="#E69F00",low="#56B4E9",na.value="light gray", limits=c(-2.5, 3.5), breaks=seq(-2, 3,by=1))

## Effect size 
# Subcortical
ggseg(.data= R2_SC, atlas = "aseg",colour = "gray", 
      mapping=aes(fill= r)) +
  labs(title="          Effect size, subcortical volume", fill="r")+
  scale_fill_continuous(low="white", high="#CC79A7", na.value="light gray", limits=c(0, 1), breaks=seq(0,1,by=.4)) +
  theme(legend.position = "bottom")

# Cortical thickness
ggseg(.data=TH_sig, colour="gray", position="stacked",
      mapping=aes(fill=r)) +
  labs(title="Effect size, cortical thickness", fill="r") +
  scale_fill_continuous(low="white", high="#CC79A7", na.value="light gray", limits=c(0, 1), breaks=seq(0,1,by=.4)) +
  theme(legend.position = "bottom")

# Cortical area
ggseg(.data=AR_sig, colour="gray", position="stacked",
      mapping=aes(fill=r)) +
  labs(title="Effect size, cortical area", fill="r") +
  scale_fill_continuous(low="white", high="#CC79A7", na.value="light gray", limits=c(0, 1), breaks=seq(0,1,by=.4)) +
  theme(legend.position = "bottom")

## Sup fig: P-value _______________________________________________________________________________________________________________________________________________________
# Subcortical
ggseg(.data= R2_SC, atlas = "aseg",colour = "gray", 
      mapping=aes(fill= pval)) +
  labs(title="          P-value, subcortical structures", fill="P-value")+
  scale_fill_continuous(low="#009E73", high="white", na.value="light gray", limits=c(0, 0.05), breaks=seq(0,0.05,by=.05)) +
  theme(legend.position = "bottom")

# Cortical thickness
ggseg(.data=TH_sig, colour="gray", position="stacked",
      mapping=aes(fill=pval)) +
  labs(title="P-value, cortical thickness", fill="P-value") +
  scale_fill_continuous(low="#009E73", high="white", na.value="light gray", limits=c(0, 0.05), breaks=seq(0,0.05,by=.05)) +
  theme(legend.position = "bottom")

# Cortical area
ggseg(.data=AR_sig, colour="gray", position="stacked",
      mapping=aes(fill=pval)) +
  labs(title="P-value, cortical area", fill="P-value") +
  scale_fill_continuous(low="#009E73", high="white", na.value="light gray", limits=c(0, 0.05), breaks=seq(0,0.05,by=.05)) +
  theme(legend.position = "bottom")


## 5. QC strategies, consequences in a regression analysis  #################################################################################################################
#############################################################################################################################################################################

## Do  a test case linear model _____________________________________________________________________________________________________________________________________________

DATAFRAMES=c("NoQC", "GoldStandard", 
             "excl_CNR", "excl_EN", "excl_MG_subcort", "excl_MG_hippo", "excl_MG_thick", "excl_MG_area","excl_MRIQC", "excl_Qoala",
             "semi_CNR", "semi_EN", "semi_MG_subcort", "semi_MG_hippo", "semi_MG_thick", "semi_MG_area","semi_MRIQC", "semi_Qoala", 
             "Manual_exclusion")

for (DATAFRAME in DATAFRAMES) {
  temp1 <- data.frame(Region=colnames(get(DATAFRAME))[c(53:ncol(get(DATAFRAME)))],Beta=NA,SE=NA,Pval=NA, R2=NA, adj.R2=NA)
  for (i in c(43:ncol(get(DATAFRAME)))) {
    templm<-summary(lm(get(colnames(get(DATAFRAME)[i])) ~  Age + SEX , data = get(DATAFRAME))) ; model = "Model1"
    #templm<-summary(lm(get(colnames(get(DATAFRAME)[i])) ~  Age + SEX +bmi, data = get(DATAFRAME))); model = "Model2_BMI"
    temp1$Beta[temp1$Region==colnames(get(DATAFRAME)[i])] <- templm$coefficients[2,1]
    temp1$SE[temp1$Region==colnames(get(DATAFRAME)[i])] <- templm$coefficients[2,2]
    temp1$Pval[temp1$Region==colnames(get(DATAFRAME)[i])] <- templm$coefficients[2,4]
    temp1$R2[temp1$Region==colnames(get(DATAFRAME)[i])] <- templm$r.squared
    temp1$adj.R2[temp1$Region==colnames(get(DATAFRAME)[i])] <- templm$adj.r.squared
    assign(paste0("LR_", DATAFRAME), temp1)
  }}

## Extract R2 ______________________________________________________________________________________________________________________________________________________________
## Select the R2 obtained by each QC dataset and place it into the same DF
R2 <- data.frame(region=LR_NoQC$Region)

for (DATAFRAME in DATAFRAMES) {
  temp<-get(paste0("LR_", DATAFRAME))[5]
  R2[ , ncol(R2) + 1] <-temp
  colnames(R2)[ncol(R2)] <- DATAFRAME
  }

## Change names to fit the paper
colnames(R2)
colnames(R2) <- c("Region", "Non-QC", "Visual-edit", 
                  "Auto-CNR", "Auto-EN", "Auto-SC_Morphological", "Auto-HC_Morphological", "Auto-TH_Morphological", "Auto-AR_Morphological", "Auto-MRIQC", "Auto-Qoala", 
                  "Semi-CNR", "Semi-EN", "Semi-SC_Morphological", "Semi-HC_Morphological", "Semi-TH_Morphological", "Semi-AR_Morphological", "Semi-MRIQC" , "Semi-Qoala", 
                  "Visual-excl" ) 


## Separate by morphological groups 
R2_TH <- R2[grep("thickness", R2$Region), ] # Thickness
R2_AR<- R2[grep("area", R2$Region), ] # Area
R2_SC <- R2[-grep("thickness", R2$Region), ] # Subcortical
R2_SC<- R2_SC[-grep("area", R2_SC$Region), ]
R2_SC <- R2_SC[-grep("ppocampaltail", R2_SC$Region), ]
R2_SC <- R2_SC[-grep("subiculum", R2_SC$Region), ]
R2_SC <- R2_SC[-grep("hippocampalfissure", R2_SC$Region), ]
R2_SC <- R2_SC[-grep("molecularlayer", R2_SC$Region), ]
R2_SC <- R2_SC[-grep(	"_CA", R2_SC$Region), ]
R2_SC <- R2_SC[-grep(	"GCMLDG", R2_SC$Region), ]
R2_SC <- R2_SC[-grep(	"_HATA", R2_SC$Region), ]
R2_SC <- R2_SC[-grep(	"fimbria", R2_SC$Region), ]
a1 <- R2[grep("ppocampaltail", R2$Region), ] # Hippocampal subfields
a2 <- R2[grep("subiculum", R2$Region), ]
a3 <- R2[grep("hippocampalfissure", R2$Region), ]
a4 <- R2[grep("molecularlayer", R2$Region), ]
a5 <- R2[grep("_CA", R2$Region), ]
a6 <- R2[grep("GCMLDG", R2$Region), ]
a7 <- R2[grep("_HATA", R2$Region), ]
a9 <- R2[grep("fimbria", R2$Region), ]
R2_HC<-rbind(a1, a2, a3, a4, a5, a6, a7, a9)
rm(a1, a2, a3, a4, a5, a6, a7, a9)

## Remove the non matching morphological outliers columns
drops <- c( 'Auto-HC_Morphological',"Auto-TH_Morphological", "Auto-AR_Morphological", 
            'Semi-HC_Morphological',"Semi-TH_Morphological", "Semi-AR_Morphological")
R2_SC <-R2_SC[ , !(names(R2_SC) %in% drops)]
drops <- c("Auto-SC_Morphological", "Auto-TH_Morphological", "Auto-AR_Morphological", 
           "Semi-SC_Morphological", "Semi-TH_Morphological", "Semi-AR_Morphological")
R2_HC <-R2_HC[ , !(names(R2_HC) %in% drops)]
drops <- c("Auto-SC_Morphological", 'Auto-HC_Morphological', "Auto-AR_Morphological", 
           "Semi-SC_Morphological", 'Semi-HC_Morphological', "Semi-AR_Morphological")
R2_TH <-R2_TH[ , !(names(R2_TH) %in% drops)]
drops <- c("Auto-SC_Morphological", 'Auto-HC_Morphological',"Auto-TH_Morphological",  
           "Semi-SC_Morphological", 'Semi-HC_Morphological',"Semi-TH_Morphological")
R2_AR <-R2_AR[ , !(names(R2_AR) %in% drops)]

## Change colnames
names(R2_SC)[names(R2_SC) == "Auto-SC_Morphological"] <- "Auto-Morphological"
names(R2_HC)[names(R2_HC) == "Auto-HC_Morphological"] <- "Auto-Morphological"
names(R2_TH)[names(R2_TH) == "Auto-TH_Morphological"] <- "Auto-Morphological"
names(R2_AR)[names(R2_AR) == "Auto-AR_Morphological"] <- "Auto-Morphological"
names(R2_SC)[names(R2_SC) == "Semi-SC_Morphological"] <- "Semi-Morphological"
names(R2_HC)[names(R2_HC) == "Semi-HC_Morphological"] <- "Semi-Morphological"
names(R2_TH)[names(R2_TH) == "Semi-TH_Morphological"] <- "Semi-Morphological"
names(R2_AR)[names(R2_AR) == "Semi-AR_Morphological"] <- "Semi-Morphological"

## Merge back to create the final R2 distributions
R2_all <- rbind(R2_AR, R2_TH, R2_HC, R2_SC)
R2_all<-as.data.frame(R2_all)

## Export for supplementary table
write_xlsx(R2_all,paste0(plotsfolder, "R2distribution.xlsx"), col_names = T)

## Calculate delta R2  for all brain regions _______________________________________________________________________________________________________________________________
DeltaR2_all <- data.frame(region=R2_all$Region)

for (i in colnames(R2_all[3:ncol(R2_all)])) {
  print(i)
  temp<- (R2_all[i] - R2_all$`Non-QC`) ## remove the initial R2 to calculate the change each QC strategy produced
  DeltaR2_all[ , ncol(DeltaR2_all) + 1] <-temp }

write.table(x = DeltaR2_all, file =  paste0(plotsfolder, "DeltaMeans_", model, ".txt"), row.names = F, col.names = T)  ## Save the R2 to compare between models

## Delta R2 means and CI for all morphological estimates
bla<-melt(DeltaR2_all)
allmean_tot <- bla %>%
  group_by(variable) %>%
  dplyr::summarise(n=n(),
                   mean=mean(value),
                   sd=sd(value)) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))

allmean_tot<-allmean_tot[c(1, 3:6)]
colnames(allmean_tot)<- c("QC", "mean", "sd", "se", "ic")

## Total change in R2, by type of region

#Remove first column & make rest numeric
R2_SC1<-R2_SC[2:ncol(R2_SC)]
R2_HC1<-R2_HC[2:ncol(R2_HC)]
R2_TH1<-R2_TH[2:ncol(R2_TH)]
R2_AR1<-R2_AR[2:ncol(R2_AR)]

R2_SC1[] <- lapply(R2_SC1, function(x) as.numeric(as.character(x)))
R2_HC1[] <- lapply(R2_HC1, function(x) as.numeric(as.character(x)))
R2_TH1[] <- lapply(R2_TH1, function(x) as.numeric(as.character(x)))
R2_AR1[] <- lapply(R2_AR1, function(x) as.numeric(as.character(x)))


# Create empty dataframe
R2s=c("R2_SC1", "R2_HC1", "R2_TH1", "R2_AR1")

for(index in 1:4){
  R2=R2s[index]
  r2name=R2
  R2 <- as.data.frame(get(R2))
  i =2:ncol(R2)
  temp_tot <- data.frame(matrix(ncol = 12, nrow = nrow(R2)))
  x <- colnames(R2[2:13])
  colnames(temp_tot) <- x
  temp_tot$`Visual-edit`<-(R2$`Visual-edit`-R2$`Non-QC`)  ## Calculate the change in R2 total
  temp_tot$`Auto-CNR`<-(R2$`Auto-CNR`-R2$`Non-QC`)  
  temp_tot$`Semi-CNR`<-(R2$`Semi-CNR`-R2$`Non-QC`)  
  temp_tot$`Auto-EN`<-(R2$`Auto-EN`-R2$`Non-QC`)  
  temp_tot$`Semi-EN`<-(R2$`Semi-EN`-R2$`Non-QC`)  
  temp_tot$`Auto-Morphological`<-(R2$`Auto-Morphological`-R2$`Non-QC`)  
  temp_tot$`Semi-Morphological`<-(R2$`Semi-Morphological`-R2$`Non-QC`)  
  temp_tot$`Auto-MRIQC`<-(R2$`Auto-MRIQC`-R2$`Non-QC`)  
  temp_tot$`Semi-MRIQC`<-(R2$`Semi-MRIQC`-R2$`Non-QC`)  
  temp_tot$`Auto-Qoala`<-(R2$`Auto-Qoala`-R2$`Non-QC`)  
  temp_tot$`Semi-Qoala`<-(R2$`Semi-Qoala`-R2$`Non-QC`)
  temp_tot$`Visual-excl`<-(R2$`Visual-excl`-R2$`Non-QC`)
  assign(paste0("Delta_",r2name) , temp_tot)}

## Delta R2 means and CI for type of region
df = c("Delta_R2_AR1", "Delta_R2_TH1", "Delta_R2_HC1", "Delta_R2_SC1")

for (i in df) {
  bla<-melt(get(i))
  mean_temp <- bla %>%
    group_by(variable) %>%
    dplyr::summarise(n=n(),
                     mean=mean(value),
                     sd=sd(value)) %>%
    mutate( se=sd/sqrt(n))  %>%
    mutate( ic=se * qt((1-0.05)/2 + .5, n-1))
    mean_temp<-mean_temp[c(1, 3:6)]
    colnames(mean_temp)<- c("QC", "mean", "sd", "se", "ic")
    assign(paste0("allmean_", i), mean_temp)}



## Figure 5 _______________________________________________________________________________________________________________________________________________________________
## Mergve dfs
AR_1<-allmean_Delta_R2_AR1[,c(1:2, 5)]; TH_1<-allmean_Delta_R2_TH1[,c(1:2, 5)]; SC_1<-allmean_Delta_R2_SC1[,c(1:2, 5)]; HC_1<-allmean_Delta_R2_HC1[,c(1:2, 5)]
mean<-allmean_tot[,c(1:2, 5)]
AR_1$Group<-"Area" ;TH_1$Group<-"Thickness"; SC_1$Group<-"Subcortical"; HC_1$Group<-"Hippocampal subfields"; mean$Group<-"Average"
meanTotal_type1<-rbind(AR_1, TH_1,SC_1, HC_1, mean)

colnames(meanTotal_type1) <- c("QC", "value", "ic", "region")
meanTotal_type1$value1<- formatC(meanTotal_type1$value, format = "e", digits = 2) ## Make a col value with only 3 decimals 
meanTotal_type1$QC<-as.factor(meanTotal_type1$QC)

# Change to regions to for paper
meanTotal_type1$region <- gsub('Subcortical', 'Subcortical volumes', meanTotal_type1$region)
meanTotal_type1$region <- gsub('Hippocampal subfields', 'Hippocampal subfields volumes', meanTotal_type1$region)
meanTotal_type1$region <- gsub('Area', 'Cortical area', meanTotal_type1$region)
meanTotal_type1$region <- gsub('Thickness', 'Cortical thickness', meanTotal_type1$region)
meanTotal_type1$region <- gsub('Average', 'Average across all types of brain measures', meanTotal_type1$region)

# Order the regions manually
meanTotal_type1$region <- ordered(meanTotal_type1$region, levels = c("Average across all types of brain measures","Cortical thickness", "Cortical area", "Subcortical volumes",
                                                                     "Hippocampal subfields volumes"))
# Order the QC manually
meanTotal_type1$QC <- ordered(meanTotal_type1$QC, levels = c("Auto-Morphological","Auto-Qoala","Semi-Morphological","Semi-Qoala","Auto-CNR", "Semi-CNR",
                                                             "Semi-MRIQC","Semi-EN", "Auto-MRIQC", "Visual-edit", "Auto-EN", "Visual-excl"))
# Plot it
plot<-
  ggplot(meanTotal_type1) +
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  geom_errorbar( aes(x=QC, ymin=value-ic, ymax=value+ic), width=0.4, colour="black", alpha=0.9, size=0.5) +
  theme_bw()+
  geom_point( aes(x=QC, y=value, size=5, fill=value), colour="black",pch=22, size=3) +
  coord_flip() +
  facet_wrap(~region, ncol=1, scale="free_y") + 
  scale_fill_gradientn(
    colors=c("#D55E00","white","#009E73"),
    values=scales::rescale(c(min(meanTotal_type1$value),0, -min(meanTotal_type1$value))),
    limits=c(min(meanTotal_type1$value),-min(meanTotal_type1$value)),) +
  theme(
    legend.position = "right",
    legend.title = element_text(size =12),
    legend.text = element_text(size = 10),
    panel.spacing = unit(0.2, "lines"),
    strip.text.x = element_text(size = 12),
    axis.text.x = element_text(size=10),
    axis.text.y = element_text(size=8),
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12)) +
  labs(fill = (expression(Delta~R^2))) +
  xlab("Qc strategy") +
  ylab(expression("Delta R squared "(Delta~R^2)))

ggsave(plot, filename = paste0(plotsfolder, "Figure5.png"), width = 190, height = 200, units = "mm")


## Supplementary figure  ________________________________________________________________________________________________________________________________________________
## QQ plot of R2 between models 

deltaR2_model1 <-fread(paste0(plotsfolder, "DeltaMeans_Model1.txt"), sep=" ", sep2="auto", dec=".", quote="\"",
                       header="auto",na.strings = "NA")

deltaR2_model2 <-fread(paste0(plotsfolder, "DeltaMeans_Model2_BMI.txt"), sep=" ", sep2="auto", dec=".", quote="\"",
                       header="auto",na.strings = "NA")


deltaR2_model2<-deltaR2_model2[,2:ncol(deltaR2_model2)]
deltaR2_model1<-deltaR2_model1[,2:ncol(deltaR2_model1)]
deltaR2_model2<-melt(deltaR2_model2)
deltaR2_model1<-melt(deltaR2_model1)

qq<-ggplot() +
  geom_point(aes(y=deltaR2_model1$value,x=deltaR2_model2$value, color= deltaR2_model1$variable), alpha = 0.5)+
  geom_abline(aes(slope = 1, intercept = 0), linetype = 2) +
  xlab(expression("Original model + BMI "(Delta~R^2))) + ylab(expression("Original model "(Delta~R^2))) +
  xlim(-0.08, 0.09) + ylim(-0.08, 0.09) +
  theme(
    legend.position = "right",
    legend.title = element_text(size =12),
    legend.text = element_text(size = 10),
    strip.text.x = element_text(size = 12),
    axis.text.x = element_text(size=10),
    axis.text.y = element_text(size=10),
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12)) +
  labs(color = (expression("QC strategy"))) 

ggsave(qq, filename = paste0(plotsfolder, "model1vs2bmi.png"), width = 140, height = 80, units = "mm")


## 7. Visual inspection agreement ##########################################################################################################################################
#############################################################################################################################################################################
ManualQC<-read.xlsx("C:/Users/Jennifer\ Monereo/Dropbox/Shared/01.2_MCIsample/data/All_Manual_QC.xlsx") 
ManualQC[ManualQC==""] <- NA
ManualQC[ManualQC=="Na"] <- NA
ManualQC <- mutate_all(ManualQC, function(x) as.numeric(as.character(x)))
ManualQC <- ManualQC[complete.cases(ManualQC),] 

ManQC<-ManualQC[-c(1,6 )]

BarDF <- data.frame(Rater=rep(c("Accorded rating", "Rater1, no SOP","Rater2, SOP","Rater1, SOP"),each=4),
                    Score=as.factor(rep(0:3,4)),Count=c(as.vector(table(ManQC$Score_agreement)),
                                                        as.vector(table(ManQC$Score_J_TP1)),0,
                                                        as.vector(table(ManQC$Score_M_TP2)),
                                                        as.vector(table(ManQC$Score_J_TP2))))

BarDF$Score <- gsub('0', '0=Perfect quality', BarDF$Score)
BarDF$Score <- gsub('1', '1=Good quality', BarDF$Score)
BarDF$Score <- gsub('2', '2=Needs manual editing', BarDF$Score)
BarDF$Score <- gsub('3', '3=Unfixable, exclude', BarDF$Score)

## Figure 3 Rater counts _________________________________________________________________________________________________________________________________________________
cbPalette <- c(  "#009E73", "#56B4E9","#CC79A7",   "#D55E00", "#999999", "#E69F00" )
sapply(BarDF, class)
BarDF$Rater<-as.factor(BarDF$Rater)
levels(BarDF$Rater)
BarDF$Rater <- ordered(BarDF$Rater, levels = c("Rater1, no SOP","Rater1, SOP","Rater2, SOP","Accorded rating"))

RaterCounts<-ggplot(BarDF,aes(x=Rater,y=Count)) + 
  theme_bw() + 
  geom_bar(position="dodge",stat="identity",aes(fill=Score)) +
  scale_fill_manual(values=cbPalette)+
  labs( x="\n Rater", y="Count \n") +
  theme(title =element_text(size=8),
        legend.title=element_text(size=7), 
        legend.text=element_text(size=7),
        axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        plot.background = element_blank(),
        axis.line = element_line(c("black")),
        panel.background = element_rect(colour = NA, fill = 'transparent'))

# Save
ggsave(RaterCounts, filename = paste0(plotsfolder, "Figure3.png"), width = 140, height = 80, units = "mm")

## How many discordant cases between reviewers? 
count(!ManualQC$Score_J_TP2 == ManualQC$Score_M_TP2)
count((!ManualQC$Score_J_TP2 == ManualQC$Score_M_TP2) & (!ManualQC$Score_J_TP2 == ManualQC$Score_M_TP2 +1) & (!ManualQC$Score_J_TP2 == ManualQC$Score_M_TP2 -1 ))
