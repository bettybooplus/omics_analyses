##########################################################################
## Section: Loading libraries
##########################################################################
require(oligo)
require(vsn)
require(pd.clariom.s.human)
require(limma)
require(Biobase)
require(affycoretools)
require(dplyr)
require(org.Hs.eg.db)
require(clusterProfiler)
require(ReactomePA)

##########################################################################
## Section: Loading and reading data
##########################################################################
fns <- list.celfiles(path="Transcriptome", full.names=TRUE)
fns <- fns[-3]
affyraw <- read.celfiles(fns)

##########################################################################
## Section: Normalization and filter low expressed genes out
##########################################################################
raw.rma <- oligo::rma(affyraw)
probids <- oligo::probeNames(affyraw)
bgcorrect <- oligo::backgroundCorrect(affyraw, "rma")
pms <- oligo::pm(bgcorrect)
rownames(pms) <- probids
normalized <- justvsn(pms)
expression <- 2^(oligo::summarize(normalized, method="medianpolish", 
                                  probes=probids))
filter <- apply(expression, 1, function(x){sum(x > 4.5) >= 2})
T_filtered <- subset(expression, filter)

##########################################################################
## Section: Differential expression of genes
##########################################################################
design <- cbind(WT_Scr=c(1,1,0,0), SCC=c(0,0,1,1)) 
fit <- lmFit(T_filtered, design)
cont.matrix <- makeContrasts(SCC-WT_Scr, levels=design)
contrast_fit <- contrasts.fit(fit, cont.matrix)
eb_fit <- eBayes(contrast_fit)
T_all <- topTable(eb_fit, sort.by="logFC", number=Inf, confint=TRUE)
fdr_res <- fdrtool::fdrtool(T_all$t, plot=FALSE, verbose=FALSE)
T_all$qval <- fdr_res$qval
T_all$lfdr <- fdr_res$lfdr
#remove missing data
del <- grep("^hsa_", T_all$ID, value=T)
T_all <- T_all[!T_all$ID %in% del, ]
T_all <- T_all[complete.cases(T_all),]
T_all <- trans_all_AvS %>% mutate(significant=abs(T_all$logFC) >= 1 & 
                T_all$qval <= 0.15)
T_sig <- T_all %>%  filter(significant)

##########################################################################
## Section: GO-term and pathway enrichment analyses 
##########################################################################
univ_h <- univ_AnnotDbPkg("org.Hs.eg.db")
IDs <- bitr(T_sig$SYMBOL, fromType="SYMBOL", toType="ENTREZID", 
                OrgDb="org.Hs.eg.db")
#GO terms. BP
ego_bp <- enrichGO(IDs$ENTREZID, keyType="ENTREZID", 
                       OrgDb="org.Hs.eg.db", ont="BP", pvalueCutoff=0.1, 
                       qvalueCutoff=0.2, universe=univ_h, pool=TRUE, 
                       readable=T)
#GO terms. MF
ego_mf <- enrichGO(IDs$ENTREZID, keyType="ENTREZID", 
                       OrgDb="org.Hs.eg.db", ont="MF", pvalueCutoff=0.1, 
                       qvalueCutoff=0.1, universe=univ_h, pool=TRUE, 
                       readable=T)
#GO terms. CC
ego_cc <- enrichGO(IDs$ENTREZID, keyType="ENTREZID", 
                       OrgDb="org.Hs.eg.db", ont="CC", pvalueCutoff=0.1, 
                       qvalueCutoff=0.1, universe=univ_h, pool=TRUE, 
                       readable=T)
#Reactome
pathways <- enrichPathway(gene=IDs$ENTREZID, pvalueCutoff=.1, readable=T)


##########################################################################
## Section: Gene Set Enrichment Analysis (GSEA)
##########################################################################
gsea <- Tc_mRNA %>% arrange(desc(logFC))
geneList <- gsea$logFC
names(geneList) <- gsea$Symbol
hallmark_T <- GSEA(geneList, TERM2GENE=db_hallmark, minGSSize=10, 
                     pvalueCutoff=0.1, verbose=FALSE)





