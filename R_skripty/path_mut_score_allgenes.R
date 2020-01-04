
##Funkce pat_mut_score_allgenes
#Slou�� k v�po�tu dr�hov�ho muta�n�ho sk�re v�ech gen�

path_mut_score_allgenes <- function(mutdata=x, signatureset){
  #Nastaven� pracovn�ho adres��e na slo�ku, v n� je ulo�en vektor gene_symbols
  setwd("F:/Jakub/SOC/R")
  #Na�ten� vektoru gene_symbols
  gene_symbols<-read.delim("gene_symbols.txt")
  gene_symbols<-as.vector(gene_symbols)
  edg <- signatureset
  
  mutlength <- corgenelength(x=mutdata, cagenes=gene_symbols$Approved.symbol, exonsize=exon.size)
  
  mutlength <- t(mutlength)
  patmutrate <- apply(mutlength, 2, sum)
  patmut0 <- which(patmutrate == 0)
  if (length(patmut0) > 0) {
    mutlength <- mutlength[, -patmut0, drop = F]
    patmutrate <- patmutrate[-patmut0]
  }
  mutrate <- mutlength
  for (p in 1:ncol(mutlength)) {
    mutrate[, p] <- mutlength[, p]/patmutrate[p]
  }
  mutrate <- mutrate[which(row.names(mutrate) %in% colnames(edg)), 
                     ]
  genefreq <- apply(edg, 2, sum)
  genefreq <- genefreq[which(names(genefreq) %in% row.names(mutrate))]
  mutratecor <- mutrate/genefreq
  signpat <- desparsify(edgx = edg, mutratecorx = mutratecor)
  return(signpat)
}