
##Shlukova analyza Ensemble
##Ensemble clustering

#Tvorba funkci pro vypocet vzdalenosti (binomial,mahalanobis)

#Funkce binomial
#custom binomial function
binomial <- function(x) {
  vegan::vegdist(x, method = "binomial")
}
#Funkce mahalanobis
#custom mahalanobis function
mahalanobis <- function(x) {
  vegan::vegdist(x, method = "mahalanobis")
}

#Nacteni knihovny diceR urcene pro provedeni shlukove analyzy Ensemble
library(diceR)


#Shlukova analaza Ensemble pouzita na matici s drahovym mutacnim skorem genu asociovanych s rakovinou
cc_cancer<-consensus_cluster(pathmutscore_can, nk = 2:7, p.item = 0.9, reps = 5,
                             algorithms = c("hc", "km","pam"),hc.method = "ward.D2",distance = c("binomial","mahalanobis"))

#Vyhodnoceni shlukove analyzy pomoci funkce consensus_evaluate
cc_cancer_evaluate<-consensus_evaluate(pathmutscore_cancer,cc_cancer)

#Vypocet hodnoty PAC
apply(cc_cancer_evaluate$pac[,-1], 1, sum)

#Urèeni cisla shluku pomoci funkce consensus_combine z peti opakovani kazdeho algoritmu
cc_cancer_combine<- consensus_combine(cc_cancer, element = "class")

#Urceni vysledneho cisla shluku pomoci funkce k_modes pro shlukovani do dvou a do ctyr shluku
cancer_ensemble2<-k_modes(cc_cancer_combine$`2`)
cancer_ensemble4<-k_modes(cc_cancer_combine$`5`)

#Porovnani casu do prvni lecby pri shlukove analyze Ensemble do dvou shluku s geny asociovanymi s rakovinou 
data_cancer<-data_order
data_cancer$cluster<-cancer_ensemble2
treatment <- rep(1,length(data_cancer$donor_relapse_interval))
treatment[is.na(data_cancer$donor_relapse_interval)] <- 0
treatment_interval <- data_cancer$donor_relapse_interval
treatment_interval[is.na(treatment_interval)] <- data_cancer[is.na(treatment_interval),"donor_interval_of_last_followup"]
treatment <- data.frame(is_treatment = treatment, treatment_interval = treatment_interval)
treatment$cluster <- data_cancer$cluster
treatment_obj<-Surv(time=treatment$treatment_interval/365, event=treatment$is_treatment)
treatment_fit<-survfit(treatment_obj~cluster,data=treatment)
ggsurvplot(treatment_fit,data=treatment,legend.title="Legenda",palette = c("blue","red"),fun="event",pval=TRUE,pval.coord = c(20, 0.1),xlab="Roky",ylab="Pravdìpodobnost léèby")

#Porovnani casu do prvni lecby pri shlukove analyze Ensemble do peti shluku s geny asociovanymi s rakovinou 
data_cancer<-data_order
data_cancer$cluster<-cancer_ensemble5
treatment <- rep(1,length(data_cancer$donor_relapse_interval))
treatment[is.na(data_cancer$donor_relapse_interval)] <- 0
treatment_interval <- data_cancer$donor_relapse_interval
treatment_interval[is.na(treatment_interval)] <- data_cancer[is.na(treatment_interval),"donor_interval_of_last_followup"]
treatment <- data.frame(is_treatment = treatment, treatment_interval = treatment_interval)
treatment$cluster <- data_cancer$cluster
treatment_obj<-Surv(time=treatment$treatment_interval/365, event=treatment$is_treatment)
treatment_fit<-survfit(treatment_obj~cluster,data=treatment)
ggsurvplot(treatment_fit,data=treatment,legend.title="Legenda",palette = c("blue","red","green","orange","grey"),fun="event",pval=TRUE,pval.coord = c(20, 0.1),xlab="Roky",ylab="Pravdìpodobnost léèby")




#Shlukova analyza Ensemble pouzita na matici s drahovym mutacnim skorem vsech genu
cc_allgenes<-consensus_cluster(pathmutscore_all, nk = 2:7, p.item = 0.9, reps = 5,
                               algorithms = c("hc","km","pam"),hc.method = "ward.D2",distance = c("binomial","mahalanobis"))
#Vyhodnoceni shlukove analyzy pomoci funkce consensus_evaluate
cc_allgenes_evaluate<-consensus_evaluate(pathmutscore_allgenes,cc_allgenes)

#Vypocet hodnoty PAC
apply(cc_allgenes_evaluate$pac[,-1], 1, sum)

#Urceni cisla shluku pomoci funkce consensus_combine z peti opakovani kazdeho algoritmu
cc_allgenes_combine<- consensus_combine(cc_allgenes, element = "class")

#Urceni vysledneho cisla shluku pomoci funkce k_modes pro shlukovani do dvou a do ctyr shluku
allgenes_ensemble2<-k_modes(cc_allgenes_combine$`5`)
allgenes_ensemble4<-k_modes(cc_allgenes_combine$`6`)

#Porovnani casu do prvni lecby pri shlukove analyze Ensemble do peti shluku se vsemi geny 
data_allgenes<-data_order
data_allgenes$cluster<-allgenes_ensemble5
treatment <- rep(1,length(data_allgenes$donor_relapse_interval))
treatment[is.na(data_allgenes$donor_relapse_interval)] <- 0
treatment_interval <- data_allgenes$donor_relapse_interval
treatment_interval[is.na(treatment_interval)] <- data_allgenes[is.na(treatment_interval),"donor_interval_of_last_followup"]
treatment <- data.frame(is_treatment = treatment, treatment_interval = treatment_interval)
treatment$cluster <- data_allgenes$cluster
treatment_obj<-Surv(time=treatment$treatment_interval/365, event=treatment$is_treatment)
treatment_fit<-survfit(treatment_obj~cluster,data=treatment)
ggsurvplot(treatment_fit,data=treatment,legend.title="Legenda",palette=c("blue","red","green","orange","grey"),fun="event",pval=TRUE,pval.coord = c(20, 0.1),xlab="Roky",ylab="Pravdìpodobnost léèby")

#Porovnani casu do prvni lecby pri shlukove analyze Ensemble do sesti shluku se vsemi geny 
data_allgenes<-data_order
data_allgenes$cluster<-allgenes_ensemble6
treatment <- rep(1,length(data_allgenes$donor_relapse_interval))
treatment[is.na(data_allgenes$donor_relapse_interval)] <- 0
treatment_interval <- data_allgenes$donor_relapse_interval
treatment_interval[is.na(treatment_interval)] <- data_allgenes[is.na(treatment_interval),"donor_interval_of_last_followup"]
treatment <- data.frame(is_treatment = treatment, treatment_interval = treatment_interval)
treatment$cluster <- data_allgenes$cluster
treatment_obj<-Surv(time=treatment$treatment_interval/365, event=treatment$is_treatment)
treatment_fit<-survfit(treatment_obj~cluster,data=treatment)
ggsurvplot(treatment_fit,data=treatment,legend.title="Legenda",palette=c("blue","red","green","orange","grey","purple"),fun="event",pval=TRUE,pval.coord = c(20, 0.1),xlab="Roky",ylab="Pravdìpodobnost léèby")
