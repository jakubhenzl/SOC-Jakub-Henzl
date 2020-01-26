
##Nactení klinickych dat pacientu

#Nastaveni pracovniho adresare na slozku, kde mate ulozene soubory donor.tsv a nature.csv (napr.: "F:/Jakub/SOC")
setwd("F:/...")

#Nactení souboru donor.tsv a nature.csv
donor<-read.delim("donor.tsv")
compl_meta <- read.csv2("nature.csv")




##Statisticka analyza vstupnich dat

#Nacteni knihoven potrebnych pro statistickou analyzu
library(plyr)
library(dplyr)
library(survival)
library(ggplot2 )
library(survminer)

#spojeni tabulek donor a nature (=compl_meta) podle identifikacniho cisla pacientu
donor <- donor[order(donor$submitted_donor_id),]
data<-merge(donor,compl_meta, by.x="submitted_donor_id", by.y="Case")

#Odstraneni pacientu s neznamym stavem IGHV nebo s IGHV jinym nez mutovane nebo nemutovane
data<-filter(data,IGHV.status=="MUT"|IGHV.status=="UNMUT")

#Vytvoreni histogramu veku pacientu pøi diagnoze
hist(donor$donor_age_at_diagnosis,
     main = "Histogram vìku pacientù pøi diagnóze",
     xlab = "Vìk",
     ylab = "Poèet pacientù",
     xlim = c(20,90),
     ylim = c(0,200),
     col="red"
)

#Prevedeni sloupce zivotniho stavu pacienta na numericke hodnoty (mrtvy=1, zivy=0)
data$donor_vital_status <- as.character(data$donor_vital_status)
data$donor_vital_status[data$donor_vital_status == "deceased"] <- "1"
data$donor_vital_status[data$donor_vital_status == "alive"] <- "0"
data$donor_vital_status <- as.numeric(data$donor_vital_status)

#Vytvoreni Kaplan-Meierovych krivek preziti skupiny s mutovanym a nemutovanym IGHV (porovnani celkove doby preziti)
survtime<-Surv(time=data$donor_survival_time/365, event=data$donor_vital_status)
fit_temp<-survfit(survtime~IGHV.status,data=data)
ggsurvplot(fit_temp,data=data,pval=TRUE,xlab = "Roky",ylab= "Pravdìpodobnost pøežití")

#Porovnani doby do lecby mezi skupinou s mutovanym a nemutovanym IGHV
treatment <- rep(1,length(data$donor_relapse_interval))
treatment[is.na(data$donor_relapse_interval)] <- 0
treatment_interval <- data$donor_relapse_interval
treatment_interval[is.na(treatment_interval)] <- data[is.na(treatment_interval),"donor_interval_of_last_followup"]
treatment <- data.frame(is_treatment = treatment, treatment_interval = treatment_interval)
treatment$IGHV.status <- data$IGHV.status
treatment_obj<-Surv(time=treatment$treatment_interval/365, event=treatment$is_treatment)
treatment_fit<-survfit(treatment_obj~IGHV.status,data=treatment)
ggsurvplot(treatment_fit,data=treatment,fun="event",pval=TRUE,xlab="Roky",ylab="Pravdìpodobnost léèby")

#Prevedeni sloupce se statusem IGHV na faktor
data$IGHV.status<-as.factor(as.character(data$IGHV.status))

