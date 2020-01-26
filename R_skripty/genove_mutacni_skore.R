
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
