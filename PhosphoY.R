##########################################################################
## Section: Loading libraries
##########################################################################
require(stringi)
require(stringr)
require(dplyr)
require(readr)
require(limma)
require(UniProt.ws)
require(org.Hs.eg.db)
require(clusterProfiler)
require(ReactomePA)

##########################################################################
## Section: Import and transform data 
##########################################################################
phospho <- read.csv("Phospho.txt", header=TRUE, sep="\t")
phospho <- dplyr::filter(phospho, Reverse !="+")
phospho <- dplyr::filter(phospho, Localization.prob >= 0.75)
phospho <- phospho[!is.na(phospho$Ratio.H.L.normalized.R1_025) &
                           !is.na(phospho$Ratio.H.L.normalized.R2_025),]
#keep only tyrosines
phosphoY <- phospho[phospho$Amino.acid=="Y",]

##########################################################################
## Section: Check phosphoriylated sites
##########################################################################
dbPAF <- read_tsv("phospho_DB/HUMAN.zip")
qphos <- read_tsv("phospho_DB/qphos_all_data.zip")
PSP <- read.table(gzfile("phospho_DB/Phosphorylation_site_dataset.gz"), 
                  skip=3, header=TRUE, sep="\t")
PSP <- PSP %>% dplyr::filter(ORGANISM=="human")
PSP[, 5] <- gsub("(.*)-.*", "\\1", PSP$MOD_RSD)
PSP$Type <- substring(PSP$MOD_RSD, 1, 1)
PSP[, 5] <- as.numeric(substring(PSP$MOD_RSD, 2))
check <- matrix(nrow=nrow(phosY), ncol=1)
for (i in 1:nrow(phosY)) {
        check[i,] <- if(nrow(dbPAF %>% filter(Uniprot==phosY[i, 1], 
                                              Position==phosY[i, 5], 
                                              Type==phosY[i, 4]))==1){
                "Reported"} else if(
                        nrow(PSP %>% filter(ACC_ID==phosY[i, 1], 
                                            MOD_RSD==phosY[i, 5],
                                            Type==phosY[i, 4]))==1){
                        "Reported"} else if(
                                nrow(qphos %>% filter(UniProt==phosY[i, 1], 
                                                      Position==phosY[i, 5]))==1){
                                "Reported"} else "New"}
check <- as.data.frame(check)
phosphoY$phospho_site <- check$V1

##########################################################################
## Section: Differential expression of proteins
##########################################################################
fit <- lmFit(log2(phosphoY[, c(13,15)]))
eb_fit <- eBayes(fit)
pYome <- topTable(eb_fit, sort.by="logFC", number=Inf, confint=TRUE)
pYome <- pYome %>% mutate(significant=abs(pYome$logFC) >= 1 & 
                                  pYome$adj.P.Val <= .05)
sig_pYome <- pYome %>% filter(significant)

##########################################################################
## Section: GO-term and pathway enrichment analyses
##########################################################################
univ_h <- univ_AnnotDbPkg("org.Hs.eg.db")
IDs <- select(org.Hs.eg.db, keys=loc_pYome$name, columns="ENTREZID",
                    keytype="SYMBOL")
#GO terms. BP
bp_pYome <- enrichGO(unique(IDs_pYome$ENTREZID), keyType="ENTREZID", 
                     OrgDb="org.Hs.eg.db", ont="BP", pvalueCutoff=0.05, 
                     qvalueCutoff=0.05, universe=univ_h, pool=TRUE, 
                     readable=T)
#GO terms. MF
ego_mf <- enrichGO(unique(IDs_pYome$ENTREZID), keyType="ENTREZID", 
                     OrgDb="org.Hs.eg.db", ont="MF", pvalueCutoff=0.1, 
                     qvalueCutoff=0.2, universe=univ_h, pool=TRUE, 
                     readable=T)
#GO terms. CC
ego_cc <- enrichGO(unique(IDs_pYome$ENTREZID), keyType="ENTREZID", 
                     OrgDb="org.Hs.eg.db", ont="CC", pvalueCutoff=0.1, 
                     qvalueCutoff=0.1, universe=univ_h, pool=TRUE, 
                     readable=T)
#Reactome
pathways <- enrichPathway(unique(IDs_pYome$ENTREZID), pvalueCutoff=0.05,
                             qvalueCutoff=0.1, readable=T)




