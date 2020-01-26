
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
