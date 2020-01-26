

##Hierarchicka shlukova analyza SAMBAR
##Hierarchical clustering SAMBAR

#Nastaveni pracovniho adresare na slozku, kde mate ulzeny soubor gene_symbols  (napø.: "F:/Jakub/SOC/R")
setwd("F:/...")

#Nacteni souboru gene_symbols.txt a prevedeni do formy vektoru
gene_symbols<-read.delim("gene_symbols.txt")
gene_symbols<-as.vector(gene_symbols)

#Nacteni knihovny SAMBAR
#Loading SAMBAR package
library(SAMBAR)

#Hierarchicka shlukova analyza pomoci funkce sambar s geny asociovanymi s rakovinou (k=2:7)
cancerass_sambar<-sambar(sign_order,signatureset="F:/Jakub/SOC/R/pathway.txt",cangenes=SAMBAR::genes,kmin=2,kmax=7)

#Zobrazeni vysledku funkce sambar s geny asociovanymi s rakovinou
for (i in 1:6) {print(table(cancerass_sambar[[i]]))}

#Prevedeni vysledneho listu do podoby tabulky
cancerass_sambar<-as.data.frame(cancerass_sambar)


#Hierarchicka shlukova analyza pomoci funkce sambar se vsemi geny (vektor gene_symbols, k=2:7) 
allgenes_sambar<-sambar(sign_order,signatureset="F:/Jakub/SOC/R/pathway.txt",cangenes=gene_symbols$Approved.symbol,kmin=2,kmax=7)

#Zobrazeni vysledku funkce sambar se vsemi geny
for (i in 1:6) {print(table(allgenes_sambar[[i]]))}

#Prevedeni vysledneho listu do podoby tabulky
allgenes_sambar<-as.data.frame(allgenes_sambar)


#Hierarchicka shlukova analyza pomoci funkce sambar s geny asociovanymi s rakovinou (k=8:25)
cancerass_sambar<-sambar(sign_order,signatureset="F:/Jakub/SOC/R/pathway.txt",cangenes=SAMBAR::genes,kmin=8,kmax=25)

#Zobrazeni vysledku funkce sambar s geny asociovanymi s rakovinou
for (i in 1:18) {print(table(cancerass_sambar[[i]]))}

#Porovnánani casu do prvni lecby mezi skupinami ze shlukove analyzy SAMBAR do 25 shluku
data_cancluster<-data_order
data_cancluster$cluster<-cancerass_sambar$X25
data_cancluster<-filter(data_cancluster,cluster==1|cluster==2|cluster==3|cluster==4|cluster==7|cluster==8|cluster==11
                        |cluster==22)


treatment <- rep(1,length(data_cancluster$donor_relapse_interval))
treatment[is.na(data_cancluster$donor_relapse_interval)] <- 0
treatment_interval <- data_cancluster$donor_relapse_interval
treatment_interval[is.na(treatment_interval)] <- data_cancluster[is.na(treatment_interval),"donor_interval_of_last_followup"]
treatment <- data.frame(is_treatment = treatment, treatment_interval = treatment_interval)
treatment$cluster <- data_cancluster$cluster
treatment_obj<-Surv(time=treatment$treatment_interval/365, event=treatment$is_treatment)
treatment_fit<-survfit(treatment_obj~cluster,data=treatment)
ggsurvplot(treatment_fit,data=treatment,fun="event",pval=TRUE,palette=c("blue","red","green","orange","grey","purple","brown","pink","yellow"),pval.coord = c(20, 0.1),,xlab="Roky",legend.title = "Legenda",ylab="Pravdìpodobnost léèby")

#Porovnani p-hodnot mezi jednotlivymi shluky
pairwise_survdiff(Surv(treatment_interval,is_treatment)~cluster,data=treatment,p.adjust.method = "fdr")







