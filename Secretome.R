##########################################################################
## Section: Loading libraries
##########################################################################
require(stringi)
require(stringr)
require(dplyr)
require(DEP)
require(UniProt.ws)
require(org.Hs.eg.db)
require(clusterProfiler)
require(ReactomePA)

##########################################################################
## Section: Import and transform data
##########################################################################
raw_data <- read.csv("proteinGroups.txt", header=TRUE, sep="\t")
raw_data <- dplyr::filter(raw_data, Reverse !="+" ,
                          Potential.contaminant != "+" , 
                          Only.identified.by.site != "+")
colnames(raw_data) <- stri_replace_all(colnames(raw_data), "", 
                                       fixed="LFQ.intensity.HS.5.")
colnames(raw_data) <- stri_replace_all(colnames(raw_data), "", fixed="T:.")
data_unique <- make_unique(raw_data, "Gene.names", "Majority.protein.IDs", 
                           delim=";")

##########################################################################
## Section: Generate a SummarizedExperiment object
##########################################################################
label <- colnames(raw_data[c(15:20)])
condition <- c(rep("KO", 3), rep("WT", 3))
replicate <- c(rep(seq(1:3), 2))
LFQ_col <- 15:20
exp_design <- data.frame(label, condition, replicate)
se_data <- make_se(data_unique, LFQ_col, exp_design)

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
mixed_imput <- impute(norm_data, fun="mixed",
                      randna=!MNAR_names, 
                      mar="knn", mnar="MinDet") 

##########################################################################
## Section: Differential expression of proteins
##########################################################################
diff_data <- test_diff(mixed_imput, type="control", control="WT")
dep <- add_rejections(diff_data, lfc=1)
dep_S <- get_results(dep)
S_sig <-  dep_S %>% filter(significant)

##########################################################################
## Section: Filter proteins by cellular location
##########################################################################
human <- UniProt.ws(taxId=9606)
annot_func <- UniProt.ws::select(human, keys=S_sig$ID[-c(58,59)],
                                 columns=c("KEYWORDS", "GO", "ENSEMBL", 
                                           "SUBCELLULAR-LOCATIONS"),
                                 keytypes="UNIPROTKB")
names(annot_func)[c(1,5)] <- c("ID","Location")
annot_func <- annot_func %>% distinct(ID, .keep_all = T)
only_S <- sig_withloc[annot_func$Location=="Secreted", ]

##########################################################################
## Section: GO-term and pathway enrichment analyses 
##########################################################################
univ_h <- univ_AnnotDbPkg("org.Hs.eg.db")
IDs <- select(org.Hs.eg.db, keys=only_S$name, columns="ENTREZID", 
                   keytype="SYMBOL")
#GO terms. BP
ego_bp <- enrichGO(IDs$ENTREZID, keyType="ENTREZID", 
                   OrgDb="org.Hs.eg.db", ont="BP", pvalueCutoff=0.01, 
                   qvalueCutoff=0.2, universe=univ_h, pool=TRUE, 
                   readable=T)
#GO terms. MF
ego_mf <- enrichGO(IDs$ENTREZID, keyType="ENTREZID", 
                   OrgDb="org.Hs.eg.db", ont="MF", pvalueCutoff=0.01, 
                   qvalueCutoff=0.2, universe=univ_h, pool=TRUE, 
                   readable=T)
#GO terms. CC
ego_cc <- enrichGO(IDs$ENTREZID, keyType="ENTREZID", 
                   OrgDb="org.Hs.eg.db", ont="CC", pvalueCutoff=0.01, 
                   qvalueCutoff=0.2, universe=univ_h, pool=TRUE, 
                   readable=T)
#Reactome
pathways <- enrichPathway(gene=IDs$ENTREZID, pvalueCutoff=0.05, readable=T)







