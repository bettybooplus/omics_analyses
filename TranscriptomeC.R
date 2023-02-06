##########################################################################
## Section: Loading libraries
##########################################################################
require(Rsamtools)
require(GenomicFeatures)
require(GenomicAlignments)
require(BiocParallel)
require(stringr)
require(limma)
require(edgeR)
require(DESeq2)
require(dplyr)
require(UniProt.ws)
require(org.Hs.eg.db)
require(clusterProfiler)
require(ReactomePA)

##########################################################################
## Section: Import count matrix and annotated genes 
##########################################################################
counts <- read.table("counts.txt", header=T)
genes <- read.delim("genome_GRCh38/GENCODE/gencode_v41/geneInfo.tab", 
                    header=FALSE, skip=1, sep="\t")
colnames(genes) <- c("Ensembl", "Symbol", "Biotype")

##########################################################################
## Section: Differential expression of genes
##########################################################################
count <- counts[, c(7:18)]
group <- rep(c(1,2), each=6)
dge <- DGEList(counts=count, group=group, genes=genes)
keep <- filterByExpr(dge)
dge <- dge[keep, keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)
design <- cbind(WT=1, KOvsWT=rep(c(0,1), each=6))
set.seed(333)
vwts <- voomWithQualityWeights(dge, design, normalize.method="none", 
                               plot=TRUE, col=group)
fit_v2 <- lmFit(vwts)
fit_v2 <- eBayes(fit_v2)
fit_v2$genes <- dge$genes
Tc_mRNA <- topTable(fit_v2, coef=2, number=Inf)
fdr_res <- fdrtool::fdrtool(Tc_mRNA$t, plot=FALSE, verbose=FALSE)
T_mRNA <- T_mRNA %>% mutate(significant=abs(Tc_mRNA$logFC) >= 1.2 & 
                                      Tc_mRNA$qval <= .05)
Tc_sig <- T_mRNA %>% filter(significant & Biotype=="protein_coding")

##########################################################################
## Section: GO-term and pathway enrichment analyses
##########################################################################
univ_h <- univ_AnnotDbPkg("org.Hs.eg.db")
IDs <- AnnotationDbi::select(org.Hs.eg.db, keys=Tc_sig$Symbol, 
                             columns="ENTREZID", keytype="SYMBOL")
#GO terms. BP
ego_bp <- enrichGO(IDs$ENTREZID, keyType="ENTREZID", 
                   OrgDb="org.Hs.eg.db", ont="BP", pvalueCutoff=0.05, 
                   qvalueCutoff=0.1, universe=univ_h, pool=TRUE, 
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
#pathways
react_rna <- enrichPathway(IDs$ENTREZID, pvalueCutoff=0.05,
                           qvalueCutoff=0.1, readable=T)

##########################################################################
## Section: Gene Set Enrichment Analysis (GSEA)
##########################################################################
gsea <- Tc_mRNA %>% arrange(desc(logFC))
geneList <- gsea$logFC
names(geneList) <- gsea$Symbol
hallmark_Tc <- GSEA(geneList, TERM2GENE=db_hallmark, minGSSize=10, 
                     pvalueCutoff=0.1, verbose=FALSE)





