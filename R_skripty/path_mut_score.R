
##Funkce pat_mut_score
#Slouží k výpoètu dráhového mutaèního skóre genù asociovaných s rakovinou

path_mut_score <- function(mutdata=x, signatureset){
  library(SAMBAR)
  edg <- signatureset
  #correct number of mutations for gene length (returns gene mutation scores)
  mutlength <- corgenelength(x=mutdata, cagenes=SAMBAR::genes, exonsize=exon.size)
  
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
