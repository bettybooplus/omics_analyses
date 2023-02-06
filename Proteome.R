##########################################################################
## Section: Loading libraries
##########################################################################
require(openxlsx)
require(dplyr)
require(DEP)
require(UniProt.ws)
require(org.Hs.eg.db)
require(clusterProfiler)
require(ReactomePA)

##########################################################################
## Section: Import and transform data
##########################################################################
raw_data <- read.xlsx("010817_LFQ.xlsx", sheet=3)
raw_data<- dplyr::filter(raw_data, is.na(Reverse), 
                              is.na(Potential.contaminant),
                              is.na(Only.identified.by.site))
data_unique <- make_unique(raw_hum_alex, "Gene.names", "Protein.IDs", 
                               delim=";")

##########################################################################
## Section: Generate a SummarizedExperiment object 
##########################################################################
label <- colnames(raw_data[c(1:3, 10:15)])
condition <- c(rep("WT", 3), rep("SCC9", 3), rep("SCC14", 3))
replicate <- rep(seq(1:3), 3)
LFQ_columns <- 1:15
exp_design <- data.frame(label, condition, replicate)
se_data <- make_se(data_unique, LFQ_columns, exp_design)

##########################################################################
## Section: Filter proteins with many missing values and normalize data
##########################################################################
filt_data <- filter_missval(se_data, thr=0)
norm_data <- normalize_vsn(filt_data)

##########################################################################
## Section: Imputation of missing values
##########################################################################
prot_MNAR <- get_df_long(norm_data) %>%
        group_by(name, condition) %>%
        dplyr::summarize(NAs=all(is.na(intensity))) %>% 
        filter(NAs) %>% pull(name) %>% unique()
MNAR_names <- names(norm_data) %in% prot_MNAR
mixed_imput <- impute(norm_data, fun="mixed", randna=!MNAR_names,
                               mar="knn", mnar="MinDet")

##########################################################################
## Section: Differential expression of proteins
##########################################################################
diff_data  <- test_diff(mixed_imput, type="control", control="WT")
dep <- add_rejections(diff_data)
dep_P <- get_results(dep)
P_sig <- dep_P %>% filter(significant)

##########################################################################
## Section: GO enrichment analysis 
##########################################################################
univ_h <- univ_AnnotDbPkg("org.Hs.eg.db")
IDs <- select(org.Hs.eg.db, keys=P_sig$name, columns="ENTREZID", keytype="SYMBOL")
#GO terms. BP
ego_bp <- enrichGO(IDs_prot$ENTREZID, keyType="ENTREZID", 
                   OrgDb="org.Hs.eg.db", ont="BP", pvalueCutoff=0.01, 
                   qvalueCutoff=0.2, universe=univ_h, pool=TRUE, 
                   readable=T)
#GO terms. MF
ego_mf <- enrichGO(IDs_prot$ENTREZID, keyType="ENTREZID", 
                   OrgDb="org.Hs.eg.db", ont="MF", pvalueCutoff=0.05, 
                   qvalueCutoff=0.2, universe=univ_h, pool=TRUE, 
                   readable=T)
#GO terms. CC
ego_cc <- enrichGO(IDs_prot$ENTREZID, keyType="ENTREZID", 
                   OrgDb="org.Hs.eg.db", ont="CC", pvalueCutoff=0.05, 
                   qvalueCutoff=0.2, universe=univ_h, pool=TRUE, 
                   readable=T)
#Reactome
pathways <- enrichPathway(gene=IDs$ENTREZID, pvalueCutoff=.1, readable=T)




