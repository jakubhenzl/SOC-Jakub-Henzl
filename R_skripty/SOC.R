
#SOC skript-Jakub Henzl




##Naètení klinických dat pacientù

#Nastavení pracovního adresáøe na složku, kde máte uložené soubory donor.tsv a nature.csv (napø.: "F:/Jakub/SOC")
setwd("F:/...")

#Naètení souborù donor.tsv a nature.csv
donor<-read.delim("donor.tsv")
compl_meta <- read.csv2("nature.csv")




##Statistická analýza vstupních dat

#Naètení knihoven potøebných pro statistickou analýzu
library(plyr)
library(dplyr)
library(survival)
library(ggplot2 )
library(survminer)

#spojení tabulek donor a nature (=compl_meta) podle identifikaèního èísla pacientù
donor <- donor[order(donor$submitted_donor_id),]
data<-merge(donor,compl_meta, by.x="submitted_donor_id", by.y="Case")

#Odstranìní pacientù s neznámým stavem IGHV nebo s IGHV jiným než mutované nebo nemutované
data<-filter(data,IGHV.status=="MUT"|IGHV.status=="UNMUT")

#Vytvoøení histogramu vìku pacientù pøi diagnóze
hist(donor$donor_age_at_diagnosis,
     main = "Histogram vìku pacientù pøi diagnóze",
     xlab = "Vìk",
     ylab = "Poèet pacientù",
     xlim = c(20,90),
     ylim = c(0,200),
     col="red"
)

#Pøevedení sloupce životního stavu pacienta na numerické hodnoty (mrtvý=1, živý=0)
data$donor_vital_status <- as.character(data$donor_vital_status)
data$donor_vital_status[data$donor_vital_status == "deceased"] <- "1"
data$donor_vital_status[data$donor_vital_status == "alive"] <- "0"
data$donor_vital_status <- as.numeric(data$donor_vital_status)

#Vytvoøení Kaplan-Meierových køivek pøežití skupiny s mutovaným a nemutovaným IGHV (porovnání celkové doby pøežití)
survtime<-Surv(time=data$donor_survival_time/365, event=data$donor_vital_status)
fit_temp<-survfit(survtime~IGHV.status,data=data)
ggsurvplot(fit_temp,data=data,pval=TRUE,xlab = "Roky",ylab= "Pravdìpodobnost pøežití")

#Porovnání doby do léèby mezi skupinou s mutovaným a nemutovaným IGHV
treatment <- rep(1,length(data$donor_relapse_interval))
treatment[is.na(data$donor_relapse_interval)] <- 0
treatment_interval <- data$donor_relapse_interval
treatment_interval[is.na(treatment_interval)] <- data[is.na(treatment_interval),"donor_interval_of_last_followup"]
treatment <- data.frame(is_treatment = treatment, treatment_interval = treatment_interval)
treatment$IGHV.status <- data$IGHV.status
treatment_obj<-Surv(time=treatment$treatment_interval/365, event=treatment$is_treatment)
treatment_fit<-survfit(treatment_obj~IGHV.status,data=treatment)
ggsurvplot(treatment_fit,data=treatment,fun="event",pval=TRUE,xlab="Roky",ylab="Pravdìpodobnost léèby")

#Pøevedení sloupce se statusem IGHV na faktor
data$IGHV.status<-as.factor(as.character(data$IGHV.status))




##Naètení dat z celoexomového sekvenování našich pacientù

#Nastavení pracovního adresáøe na složku, kde máte uložené soubory všech pacientù z celoexomového sekvenování  (napø.: "F:/Jakub/SOC/newfilter") 
setwd("F:/...")

#Naètení názvù souborù z celoexomového sekvenování do listu files
files = list.files(pattern="*.extracted")

#Naètení souborù z celoexomového sekvenování do listu s názvem myfiles
myfiles = lapply(files, read.delim, header=FALSE,
                 
                 col.names=c("variation","consequences","impact", "symbol"),
                 
                 colClasses=c("factor","factor","factor","factor"))

#Odstranìní z názvu souborù èást _extracted_gnomAD
names<-as.factor(sub("_extracted_gnomAD","",files))
names(myfiles) <- names




##Vytvoøení matice s genovým mutaèním skórem našich pacientù

#Odstranìní identických transkriptù
#Remove identical transcripts
var_symbol<-lapply(myfiles, function(x) {paste(x$variation,x$symbol)})
for (i in 1:length(myfiles)) {myfiles[[i]]$var_symbol<-var_symbol[[i]]}
unique_genes<-lapply(myfiles, function(x) { x[!duplicated(x$var_symbol),] })

#Vytvoøení vektoru s názvy všech genù, u kterých se u našich pacientù vyskytovaly muatce (=vektor all_genes)
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

#Vytvoøení prázdné matice pro dosazení genového mutaèního skóre
sign_mat <- matrix(0, nrow(gen_mat), length(all_genes))
row.names(sign_mat) <- row.names(gen_mat)
colnames(sign_mat) <- all_genes

#Pøevedení této matice do podoby tabulky
#Transform this matrix into data frame
sign_mat<-as.data.frame(sign_mat)

#Výpoèet genového mutaèního skóre
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




##Propojení tabulky s genovým mutaèním skórem(tabulka sign_mat) s tabulkou s klinickými daty pacientù(tabulka data) 

#Vytvoøení sloupce s identifikaèním èíslem v tabulce data, které bude odpovídat názvùm sloupcù v tabulce s genovým mutaèním skórem
library(stringr)
number<-str_pad(data$submitted_donor_id,3,pad="0")
number<-paste(number,"ND",sep="")
data$id<-number

#Upravení poètu pacientù a jejich poøadí v tabulce data, tak aby odpovídali poøadí pacientùm v tabulce sign_mat
data_mod<-data[data$id %in% rownames(sign_mat),] 
data_order<-data_mod[order(data_mod$id),]

#Upravení poètu pacientù a jejich poøadí v tabulce sign_mat, tak aby odpovídali poøadí pacientùm v tabulce data
sign_mat_mod<-sign_mat[rownames(sign_mat) %in% data$id,]  
sign_order<-sign_mat_mod[order(rownames(sign_mat_mod)),]  




##Hierarchická shluková analýza SAMBAR
##Hierarchical clustering SAMBAR

#Nastavení pracovního adresáøe na složku, kde máte uložený soubor gene_symbols  (napø.: "F:/Jakub/SOC/R")
setwd("F:/...")

#Naètení souboru gene_symbols.txt a pøevedení do formy vektoru
gene_symbols<-read.delim("gene_symbols.txt")
gene_symbols<-as.vector(gene_symbols)

#Naètení knihovny SAMBAR
#Loading SAMBAR package
library(SAMBAR)

#Hierarchická shluková analýza pomocí funkce sambar s geny asociovanými s rakovinou
cancerass_sambar<-sambar(sign_order,signatureset="F:/Jakub/SOC/R/pathway.txt",cangenes=SAMBAR::genes,kmin=2,kmax=7)

#Zobrazení výsledkù funkce sambar s geny asociovanými s rakovinou
for (i in 1:6) {print(table(cancerass_sambar[[i]]))}

#Pøevedení výsledného listu do podoby tabulky
cancerass_sambar<-as.data.frame(cancerass_sambar)


#Hierarchická shluková analýza pomocí funkce sambar se všemi geny (vektor gene_symbols) 
hallgenes_sambar<-sambar(sign_order,signatureset="F:/Jakub/SOC/R/pathway.txt",cangenes=gene_symbols$Approved.symbol,kmin=2,kmax=7)

#Zobrazení výsledkù funkce sambar se všemi geny
for (i in 1:6) {print(table(cancerass_sambar[[i]]))}

#Pøevedení výsledného listu do podoby tabulky
allgenes_sambar<-as.data.frame(allgenes_sambar)

#Porovnání èasu do léèby u pacientù zaøazených do dvou velkých shlukù (n1=99 a n2=309) pøi analýze se všemi geny do pìti shlukù
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




##Výpoèet dráhového mutaèního skóre
##Calculation of pathway mutation scores

#Naètení funkcí potøebných pro výpoèet dráhového mutaèního skóre (vždy úplná cesta k místu uložištì funkce ve vašem poèítaèi napø: source('F:/Jakub/SOC/R/path_mut_score.R'))
source('F:/.../path_mut_score.R')
source('F:/.../path_mut_score_allgenes.R')
source('F:/.../convertgmt_mod.R')

#Vypoètení dráhového mutaèního skóre genù asociovaných s rakovinou
convert_cancerass<-convertgmt_mod(signature="pathway.txt",cagenes=SAMBAR::genes)
pathmutscore_cancer<-path_mut_score(sign_order,convert_cancerass)

#Vypoètení dráhového mutaèního skóre všech genù
convert_allgenes<-convertgmt_mod(signature="pathway.txt",cagenes=gene_symbols$Approved.symbol)
pathmutscore_allgenes<-path_mut_score_allgenes(sign_order,convert_allgenes)

#Prohození názvù øádkù a sloupcù v maticích s dráhovým mutaèním skórem
pathmutscore_allgenes<-t(pathmutscore_allgenes)
pathmutscore_cancer<-t(pathmutscore_cancer)




##Shluková analýza Ensemble
##Ensemble clustering

#Tvorba funkcí pro výpoèet vzdáleností (binomial,mahalanobis)
rf <- function(x) {
  tmp <- randomForest::randomForest(x = x, y = NULL, ntree = 1000, proximity = TRUE, oob.prox = TRUE)
  stats::as.dist(1 - tmp$proximity)
}
#custom binomial function
binomial <- function(x) {
  vegan::vegdist(x, method = "binomial")
}
#custom mahalanobis function
mahalanobis <- function(x) {
  vegan::vegdist(x, method = "mahalanobis")
}

#Naètení knihovny diceR urèené pro provedení shlukové analýzy Ensemble
library(diceR)


#Shluková analýza Ensemble použitá na matici s dráhovým mutaèním skórem genù asociovaných s rakovinou
cc_cancer<-consensus_cluster(pathmutscore_cancer, nk = 2:4, p.item = 0.9, reps = 5,
                             algorithms = c("hc", "km","pam"),hc.method = "ward.D2",distance = c("binomial","mahalanobis","rf"))

#Vyhodnocení shlukové analýzy pomocí funkce consensus_evaluate
cc_cancer_evaluate<-consensus_evaluate(pathmutscore_cancer,cc_cancer)

#Výpoèet hodnoty PAC
apply(cc_cancer_evaluate$pac[,-1], 1, sum)

#Urèení èísla shluku pomocí funkce consensus_combine z pìti opakování každého algoritmu
cc_cancer_combine<- consensus_combine(cc_cancer, element = "class")

#Urèení výsledného èísla shluku pomocí funkce k_modes pro shlukování do dvou a do ètyø shlukù
cancer_ensemble2<-k_modes(cc_cancer_combine$`2`)
cancer_ensemble4<-k_modes(cc_cancer_combine$`4`)

#Porocnání èasu do léèby pøi shlukové analýze Ensemble do dvou shlukù s geny asociovanými s rakovinou 
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

#Porocnání èasu do léèby pøi shlukové analýze Ensemble do ètyø shlukù s geny asociovanými s rakovinou 
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


#Shluková analýza Ensemble použitá na matici s dráhovým mutaèním skórem všech genù
cc_allgenes<-consensus_cluster(pathmutscore_allgenes, nk = 2:4, p.item = 0.9, reps = 5,
                               algorithms = c("hc","km","pam"),hc.method = "ward.D2",distance = c("binomial","mahalanobis","rf"))

#Vyhodnocení shlukové analýzy pomocí funkce consensus_evaluate
cc_allgenes_evaluate<-consensus_evaluate(pathmutscore_allgenes,cc_allgenes)

#Výpoèet hodnoty PAC
apply(cc_allgenes_evaluate$pac[,-1], 1, sum)

#Urèení èísla shluku pomocí funkce consensus_combine z pìti opakování každého algoritmu
cc_allgenes_combine<- consensus_combine(cc_allgenes, element = "class")

#Urèení výsledného èísla shluku pomocí funkce k_modes pro shlukování do dvou a do ètyø shlukù
allgenes_ensemble2<-k_modes(cc_allgenes_combine$`2`)
allgenes_ensemble4<-k_modes(cc_allgenes_combine$`4`)

#Porocnání èasu do léèby pøi shlukové analýze Ensemble do dvou shlukù s geny asociovanými s rakovinou 
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

#Porocnání èasu do léèby pøi shlukové analýze Ensemble do ètyø shlukù s geny asociovanými s rakovinou 
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









