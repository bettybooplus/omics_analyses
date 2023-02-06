##########################################################################
## Section: Loading required libraries
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
require(org.Mm.eg.db)
require(clusterProfiler)

##########################################################################
## Section: Import count matrix and annotated genes 
##########################################################################
counts <- read.table("analysis/counts.txt", header=T)
genes <- read.delim("genome_GRCm39/GENCODE/gencode_v30/geneInfo.tab", 
                    header=FALSE, skip=1, sep="\t")
colnames(genes) <- c("Ensembl", "Symbol", "Biotype")

##########################################################################
## Section: Differential expression of genes
##########################################################################
count <- count_mat[, c(7:ncol(counts))]
group <- rep(c(1,2), each=3)
dge <- DGEList(counts=count, group=group, genes=genes)
keep <- filterByExpr(dge)
dge <- dge[keep, keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)
design <- cbind(WT=1, KOvsWT=rep(c(0,1), each=3))
set.seed(333)
vwts <- voomWithQualityWeights(dge, design, normalize.method="none", 
                               plot=TRUE, col=group)
fit_v2 <- lmFit(vwts)
fit_v2 <- eBayes(fit_v2)
fit_v2$genes <- dge$genes
murine_mRNA <- topTable(fit_v2, coef=2, number=Inf)
fdr_res <- fdrtool::fdrtool(murine_mRNA$t, plot=FALSE, verbose=FALSE)
murine_mRNA <- murine_mRNA %>% mutate(significant=abs(murine_mRNA$logFC) >= 1.2 & 
                                      murine_mRNA$qval <= .05)
murine_sig <- murine_mRNA %>% filter(significant & Biotype=="protein_coding")  

##########################################################################
## Section: GO enrichment analysis
##########################################################################
univ_m <- univ_AnnotDbPkg("org.Mm.eg.db")
IDs <- AnnotationDbi::select(org.Mm.eg.db, keys=murine_sig$Symbol, 
                             columns="ENTREZID",keytype="SYMBOL")
#GO terms. BP
bp_rna <- enrichGO(IDs$ENTREZID, keyType="ENTREZID", 
                   OrgDb="org.Mm.eg.db", ont="BP", pvalueCutoff=0.01, 
                   qvalueCutoff=0.01, universe=univ_m, readable=T)
#GO terms. MF
ego_mf <- enrichGO(IDs$ENTREZID, keyType="ENTREZID", 
                   OrgDb="org.Mm.eg.db", ont="MF", pvalueCutoff=0.01, 
                   qvalueCutoff=0.01, universe=univ_m, readable=T)
#GO terms. CC
ego_cc <- enrichGO(IDs$ENTREZID, keyType="ENTREZID", 
                   OrgDb="org.Mm.eg.db", ont="CC", pvalueCutoff=0.01, 
                   qvalueCutoff=0.01, universe=univ_m, pool=TRUE, 
                   readable=T)




