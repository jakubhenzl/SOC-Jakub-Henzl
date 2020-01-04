
#SOC skript-Jakub Henzl




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




##Nacteni dat z celoexomoveho sekvenovani nasich pacientu

#Nastaveni pracovniho adresare na slozku, kde mate ulozene soubory vsech pacientu z celoexomoveho sekvenovani  (napr.: "F:/Jakub/SOC/newfilter") 
setwd("F:/...")

#Nacteni nazve souboru z celoexomoveho sekvenovani do listu files
files = list.files(pattern="*.extracted")

#Nacteni souboru z celoexomoveho sekvenovani do listu s nazvem myfiles
myfiles = lapply(files, read.delim, header=FALSE,
                 
                 col.names=c("variation","consequences","impact", "symbol"),
                 
                 colClasses=c("factor","factor","factor","factor"))

#Odstraneni z nazvu souboru cast _extracted_gnomAD
names<-as.factor(sub("_extracted_gnomAD","",files))
names(myfiles) <- names




##Vytvoreni matice s genovým mutacnim skorem nasich pacientu

#Odstraneni identickych transkriptu
#Remove identical transcripts
var_symbol<-lapply(myfiles, function(x) {paste(x$variation,x$symbol)})
for (i in 1:length(myfiles)) {myfiles[[i]]$var_symbol<-var_symbol[[i]]}
unique_genes<-lapply(myfiles, function(x) { x[!duplicated(x$var_symbol),] })

#Vytvorení vektoru s nazvy vsech genu, u kterych se u nasich pacientu vyskytovaly muatce (=vektor all_genes)
df_ncol <- max(sapply(unique_genes,dim)[1,])
df_nrow <- length(unique_genes)
df_zero <- data.frame(matrix(NA, nrow=df_nrow, ncol=df_ncol))
rownames(df_zero) <- names
j=1
for (i in unique_genes){
  df_zero[j,]<-c(as.character(i$symbol),rep("NA",(df_ncol-length(i$symbol))))
  j=j+1
}
gen_mat <- as.matrix(df_zero)
all_genes <- unique(unlist(c(gen_mat)))
all_genes <- all_genes[order(all_genes)]
all_genes <- all_genes[which(all_genes != "NA")]

#Vytvoreni prazdne matice pro dosazeni genoveho mutacniho skore
sign_mat <- matrix(0, nrow(gen_mat), length(all_genes))
row.names(sign_mat) <- row.names(gen_mat)
colnames(sign_mat) <- all_genes

#Prevedeni teto matice do podoby tabulky
#Transform this matrix into data frame
sign_mat<-as.data.frame(sign_mat)

#Vypocet genového mutaniho skore
#Calculation of gene-level mutation score
u<-nrow(sign_mat)
w=1
for (i in 1:u) {
  quantity<-table(unique_genes[[i]]$symbol)
  
  genes<-names(quantity)
  quantity<-as.vector(quantity)
  y=1
  
  v<-length(all_genes)
  for (z in 1:v){
    if(identical(genes[y],all_genes[z])){sign_mat[w,z]=quantity[y];y=y+1} else {sign_mat[w,z]=0}}
  w=w+1
} 




##Propojeni tabulky s genovym mutacnim skorem(tabulka sign_mat) s tabulkou s klinickymi daty pacientu(tabulka data) 

#Vytvoreni sloupce s identifikacnim cislem v tabulce data, ktere bude odpovidat nazvum sloupcu v tabulce s genovym mutacnim skorem
library(stringr)
number<-str_pad(data$submitted_donor_id,3,pad="0")
number<-paste(number,"ND",sep="")
data$id<-number

#Upraveni poètu pacientu a jejich poradi v tabulce data, tak aby odpovidali poradi pacientum v tabulce sign_mat
data_mod<-data[data$id %in% rownames(sign_mat),] 
data_order<-data_mod[order(data_mod$id),]

#Upraveni poctu pacientu a jejich poradi v tabulce sign_mat, tak aby odpovidali poradi pacientum v tabulce data
sign_mat_mod<-sign_mat[rownames(sign_mat) %in% data$id,]  
sign_order<-sign_mat_mod[order(rownames(sign_mat_mod)),]  




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

#Hierarchicka shlukova analyza pomoci funkce sambar s geny asociovanymi s rakovinou
cancerass_sambar<-sambar(sign_order,signatureset="F:/Jakub/SOC/R/pathway.txt",cangenes=SAMBAR::genes,kmin=2,kmax=7)

#Zobrazeni vysledku funkce sambar s geny asociovanymi s rakovinou
for (i in 1:6) {print(table(cancerass_sambar[[i]]))}

#Prevedeni vysledneho listu do podoby tabulky
cancerass_sambar<-as.data.frame(cancerass_sambar)


#Hierarchicka shlukova analyza pomoci funkce sambar se vsemi geny (vektor gene_symbols) 
hallgenes_sambar<-sambar(sign_order,signatureset="F:/Jakub/SOC/R/pathway.txt",cangenes=gene_symbols$Approved.symbol,kmin=2,kmax=7)

#Zobrazeni vysledku funkce sambar se vsemi geny
for (i in 1:6) {print(table(cancerass_sambar[[i]]))}

#Prevedeni vysledneho listu do podoby tabulky
allgenes_sambar<-as.data.frame(allgenes_sambar)

#Porovnani casu do lecby u pacientu zarazenych do dvou velkych shluku (n1=99 a n2=309) pri analyze se vsemi geny do peti shluku
data_allgenescluster<-data[data$id%in%rownames(allgenes_sambar),]
data_allgenescluster$cluster<-allgenes_sambar$X5
data_allgenescluster_main<-filter(data_allgenescluster,cluster==1|cluster==2)
treatment <- rep(1,length(data_allgenescluster_main$donor_relapse_interval))
treatment[is.na(data_allgenescluster_main$donor_relapse_interval)] <- 0
treatment_interval <- data_allgenescluster_main$donor_relapse_interval
treatment_interval[is.na(treatment_interval)] <- data_allgenescluster_main[is.na(treatment_interval),"donor_interval_of_last_followup"]
treatment <- data.frame(is_treatment = treatment, treatment_interval = treatment_interval)
treatment$cluster <- data_allgenescluster_main$cluster
treatment_obj<-Surv(time=treatment$treatment_interval/365)
treatment_fit<-survfit(treatment_obj~cluster,data=treatment)
ggsurvplot(treatment_fit,data=treatment,fun="event",pval=TRUE,pval.coord = c(20, 0.1),xlab="Roky",ylab="Pravdìpodobnost léèby", 
           legend.title = "Legenda")




##Vypocet drahoveho mutacniho skore
##Calculation of pathway mutation scores

#Nacteni funkci potrebnych pro vypocet drahoveho mutacniho skore (vzdy uplna cesta k mistu uloziste funkce ve vasem pocitaci napr: source('F:/Jakub/SOC/R/path_mut_score.R'))
source('F:/.../path_mut_score.R')
source('F:/.../path_mut_score_allgenes.R')
source('F:/.../convertgmt_mod.R')

#Vypocteni drahoveho mutacniho skore genu asociovanych s rakovinou
convert_cancerass<-convertgmt_mod(signature="pathway.txt",cagenes=SAMBAR::genes)
pathmutscore_cancer<-path_mut_score(sign_order,convert_cancerass)

#Vypocteni drahoveho mutacniho skore vsech genu
convert_allgenes<-convertgmt_mod(signature="pathway.txt",cagenes=gene_symbols$Approved.symbol)
pathmutscore_allgenes<-path_mut_score_allgenes(sign_order,convert_allgenes)

#Prohozeni nazvu radku a sloupcu v maticich s drahovym mutacnim skorem
pathmutscore_allgenes<-t(pathmutscore_allgenes)
pathmutscore_cancer<-t(pathmutscore_cancer)




##Shlukova analyza Ensemble
##Ensemble clustering

#Tvorba funkci pro vypocet vzdalenosti (binomial,mahalanobis)
rf <- function(x) {
  tmp <- randomForest::randomForest(x = x, y = NULL, ntree = 1000, proximity = TRUE, oob.prox = TRUE)
  stats::as.dist(1 - tmp$proximity)
}
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
cc_cancer<-consensus_cluster(pathmutscore_cancer, nk = 2:4, p.item = 0.9, reps = 5,
                             algorithms = c("hc", "km","pam"),hc.method = "ward.D2",distance = c("binomial","mahalanobis","rf"))

#Vyhodnoceni shlukove analyzy pomoci funkce consensus_evaluate
cc_cancer_evaluate<-consensus_evaluate(pathmutscore_cancer,cc_cancer)

#Vypocet hodnoty PAC
apply(cc_cancer_evaluate$pac[,-1], 1, sum)

#Urèeni cisla shluku pomoci funkce consensus_combine z peti opakovani kazdeho algoritmu
cc_cancer_combine<- consensus_combine(cc_cancer, element = "class")

#Urceni vysledneho cisla shluku pomoci funkce k_modes pro shlukovani do dvou a do ctyr shluku
cancer_ensemble2<-k_modes(cc_cancer_combine$`2`)
cancer_ensemble4<-k_modes(cc_cancer_combine$`4`)

#Porovnani casu do lecby pri shlukove analyze Ensemble do dvou shluku s geny asociovanymi s rakovinou 
patients<-row.names(pathmutscore_cancer)
names(cancer_ensemble2)<-patients
data_cancer<-data[data$id %in% patients,]
data_cancer$cluster<-cancer_ensemble2
treatment <- rep(1,length(data_cancer$donor_relapse_interval))
treatment[is.na(data_cancer$donor_relapse_interval)] <- 0
treatment_interval <- data_cancer$donor_relapse_interval
treatment_interval[is.na(treatment_interval)] <- data_cancer[is.na(treatment_interval),"donor_interval_of_last_followup"]
treatment <- data.frame(is_treatment = treatment, treatment_interval = treatment_interval)
treatment$cluster <- data_cancer$cluster
treatment_obj<-Surv(time=treatment$treatment_interval/365, event=treatment$is_treatment)
treatment_fit<-survfit(treatment_obj~cluster,data=treatment)
ggsurvplot(treatment_fit,data=treatment,legend.title="Legenda",fun="event",pval=TRUE,xlab="Roky",ylab="Pravdìpodobnost léèby")

#Porovnani casu do lecby pri shlukove analyze Ensemble do ctyr shluku s geny asociovanymi s rakovinou 
patients<-row.names(pathmutscore_cancer)
names(cancer_ensemble4)<-patients
data_cancer<-data[data$id %in% patients,]
data_cancer$cluster<-cancer_ensemble4
treatment <- rep(1,length(data_cancer$donor_relapse_interval))
treatment[is.na(data_cancer$donor_relapse_interval)] <- 0
treatment_interval <- data_cancer$donor_relapse_interval
treatment_interval[is.na(treatment_interval)] <- data_cancer[is.na(treatment_interval),"donor_interval_of_last_followup"]
treatment <- data.frame(is_treatment = treatment, treatment_interval = treatment_interval)
treatment$cluster <- data_cancer$cluster
treatment_obj<-Surv(time=treatment$treatment_interval/365, event=treatment$is_treatment)
treatment_fit<-survfit(treatment_obj~cluster,data=treatment)
ggsurvplot(treatment_fit,data=treatment,legend.title="Legenda",fun="event",pval=TRUE,xlab="Roky",ylab="Pravdìpodobnost léèby")


#Shlukova analyza Ensemble pouzita na matici s drahovym mutacnim skorem vsech genu
cc_allgenes<-consensus_cluster(pathmutscore_allgenes, nk = 2:4, p.item = 0.9, reps = 5,
                               algorithms = c("hc","km","pam"),hc.method = "ward.D2",distance = c("binomial","mahalanobis","rf"))

#Vyhodnoceni shlukove analyzy pomoci funkce consensus_evaluate
cc_allgenes_evaluate<-consensus_evaluate(pathmutscore_allgenes,cc_allgenes)

#Vypocet hodnoty PAC
apply(cc_allgenes_evaluate$pac[,-1], 1, sum)

#Urceni cisla shluku pomoci funkce consensus_combine z peti opakovani kazdeho algoritmu
cc_allgenes_combine<- consensus_combine(cc_allgenes, element = "class")

#Urceni vysledneho cisla shluku pomoci funkce k_modes pro shlukovani do dvou a do ctyr shluku
allgenes_ensemble2<-k_modes(cc_allgenes_combine$`2`)
allgenes_ensemble4<-k_modes(cc_allgenes_combine$`4`)

#Porovnani casu do lecby pri shlukove analyze Ensemble do dvou shluku se vsemi geny   
patients<-row.names(pathmutscore_allgenes)
names(allgenes_ensemble4)<-patients
data_allgenes<-data[data$id %in% patients,]
data_allgenes$cluster<-allgenes_ensemble4
treatment <- rep(1,length(data_allgenes$donor_relapse_interval))
treatment[is.na(data_allgenes$donor_relapse_interval)] <- 0
treatment_interval <- data_allgenes$donor_relapse_interval
treatment_interval[is.na(treatment_interval)] <- data_allgenes[is.na(treatment_interval),"donor_interval_of_last_followup"]
treatment <- data.frame(is_treatment = treatment, treatment_interval = treatment_interval)
treatment$cluster <- data_allgenes$cluster
treatment_obj<-Surv(time=treatment$treatment_interval/365, event=treatment$is_treatment)
treatment_fit<-survfit(treatment_obj~cluster,data=treatment)
ggsurvplot(treatment_fit,data=treatment,legend.title="Legenda",fun="event",pval=TRUE,xlab="Roky",ylab="Pravdìpodobnost léèby")

#Porovnani casu do lecby pri shlukove analyze Ensemble do ctyr shluku se vsemi geny  
patients<-row.names(pathmutscore_allgenes)
names(allgenes_ensemble4)<-patients
data_allgenes<-data[data$id %in% patients,]
data_allgenes$cluster<-allgenes_ensemble4
treatment <- rep(1,length(data_allgenes$donor_relapse_interval))
treatment[is.na(data_allgenes$donor_relapse_interval)] <- 0
treatment_interval <- data_allgenes$donor_relapse_interval
treatment_interval[is.na(treatment_interval)] <- data_allgenes[is.na(treatment_interval),"donor_interval_of_last_followup"]
treatment <- data.frame(is_treatment = treatment, treatment_interval = treatment_interval)
treatment$cluster <- data_allgenes$cluster
treatment_obj<-Surv(time=treatment$treatment_interval/365, event=treatment$is_treatment)
treatment_fit<-survfit(treatment_obj~cluster,data=treatment)
ggsurvplot(treatment_fit,data=treatment,legend.title="Legenda",fun="event",pval=TRUE,xlab="Roky",ylab="Pravdìpodobnost léèby")









