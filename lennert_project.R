# Volcano plots
# Volcano plots visualizing DE genes between primary and met in IMmotion151 and Braunetal
results_dir =  "/work/reznik/xiea1/MIRTH/revision/results_lennert_project/braunetal/primary_met/"

df <- read.csv(paste0(results_dir,"differential_genes.csv"), row.names = 1) %>% filter(log2fc > -100) # filter out values with -Inf log2fc
plot <- EnhancedVolcano(df,
                        lab = rownames(df),
                        x = 'log2fc',
                        y = "p_adj",
                        title = paste0('Genes overexpressed in Primary vs Metastasis, braunetal'),
                        pCutoff = 0.1,
                        FCcutoff = 0,
                        selectLab = c('VHL', 'PBRM1', 'SETD2', 'BAP1', 'MTOR', 'KDM5C', 'PIK3CA', 
                                      'PIK3R1', 'PTEN', 'TP53', 'TSC1', 'TSC2', 'TCEB1', 'SDHB', 'FH'),
                        legendPosition = 'none',
                        pointSize = 3.0,
                        labSize = 4.0,
                        labCol = 'black',
                        boxedLabels = FALSE,
                        colAlpha = 4/5,
                        drawConnectors = TRUE,
                        widthConnectors = 0.75,
                        colConnectors = 'black',
                        typeConnectors = 'closed',
                        endsConnectors = 'last',
                        max.overlaps = Inf
)  
  
  plot <- plot + 
    ggplot2::ylab("-log10(p_adj) ") +
    ggplot2::xlab("log2fc(metastasis/primary)") +
    ggplot2::scale_x_continuous(limits = c(min(df$log2fc)-0.1, max(df$log2fc)+0.1))+
    ggplot2::scale_y_continuous(limits = c(0,max(-log10(df$p_adj))+0.5)) 
  ggsave(file.path(paste0(results_dir,"volcano_DE_genes_braunetal.pdf")),width = 10,height = 10) 


# Scatter plot comparing the results of DE genes in IMmotion151 and Braunetal
results_dir =  "/work/reznik/xiea1/MIRTH/revision/results_lennert_project/immotion151/primary_met/"
tcga_dir = "/work/reznik/xiea1/MIRTH/revision/results_lennert_project/immotion151/primary_met/"
cptac_dir = "/work/reznik/xiea1/MIRTH/revision/results_lennert_project/braunetal/primary_met/"

cptac <- read.csv(paste0(cptac_dir,"differential_genes.csv"), row.names = 1)
tcga <- read.csv(paste0(tcga_dir,"differential_genes.csv"), row.names = 1)
tcga <- dplyr::rename(tcga, log2fc_tcga = log2fc)
tcga <- dplyr::rename(tcga, p_tcga = p)
tcga <- dplyr::rename(tcga, p_adj_tcga = p_adj)
df <- merge(cptac, tcga, by = "row.names")
colnames(df)[1] <- "gene"
df <- df %>% 
  mutate(sig.status = case_when(
    p_adj_tcga >= 0.1 & p_adj >= 0.1 ~ "not significant",
    p_adj_tcga < 0.1 & p_adj >= 0.1   ~ "significant only in IMmotion151",
    p_adj_tcga >= 0.1 & p_adj < 0.1 ~ "significant only in Braunetal",
    log2fc > 0 & log2fc_tcga > 0 ~ "significant in both, consistent",
    log2fc < 0 & log2fc_tcga < 0 ~ "significant in both, consistent",
    T ~ "significant in both, inconsistent")) 

count <- df %>% count(sig.status)
df <- df %>% filter(sig.status!="not significant")

df %>% 
  ggplot(aes(x=log2fc,
             y=log2fc_tcga,
             color=sig.status))+
  geom_point(alpha=0.8)+
  scale_color_manual(values = c( "significant only in IMmotion151" = "#DDCC77",
                                "significant only in Braunetal" = "#88CCEE",
                                "significant in both, consistent" = "#332288",
                                "significant in both, inconsistent" = "red"))+
  #scale_x_continuous(breaks = c(-0.5,0,0.5),limits = c(-0.5,0.5))+
  #scale_y_continuous(breaks = c(-0.5,0,0.5),limits = c(-0.5,0.5))+
  #coord_equal()+
  theme_Publication()+
  labs(x="Braunetal Met/Primary log2(fold change)",
       y="IMmotion151 Met/Primary log2(fold change)",
       title = "Braunetal vs IMmotion151")
ggsave(file.path(paste0(results_dir,"scatterplot_151_br.pdf")),width = 10,height = 10) 
cor.test(df$log2fc, df$log2fc_tcga, method = "spearman")

# Scatter plot comparing the results of DE genes between Responder and Non-responder of Primary and Met 
results_dir =  "/work/reznik/xiea1/MIRTH/revision/results_lennert_project/immotion151/primary_met/"
tcga_dir = "/work/reznik/xiea1/MIRTH/revision/results_lennert_project/immotion151/primary_met/"
cptac_dir = "/work/reznik/xiea1/MIRTH/revision/results_lennert_project/immotion151/primary_met/"

cptac <- read.csv(paste0(cptac_dir,"differential_genes_response_C_P.csv"), row.names = 1)
tcga <- read.csv(paste0(tcga_dir,"differential_genes_response_C_M.csv"), row.names = 1)
tcga <- dplyr::rename(tcga, log2fc_tcga = log2fc)
tcga <- dplyr::rename(tcga, p_tcga = p)
tcga <- dplyr::rename(tcga, p_adj_tcga = p_adj)
df <- merge(cptac, tcga, by = "row.names")
colnames(df)[1] <- "metabolite"
df <- df %>% 
  mutate(sig.status = case_when(
    p_adj_tcga >= 0.1 & p_adj >= 0.1 ~ "not significant",
    p_adj_tcga < 0.1 & p_adj >= 0.1   ~ "significant only in Metastatic samples",
    p_adj_tcga >= 0.1 & p_adj < 0.1 ~ "significant only in Primary samples",
    log2fc > 0 & log2fc_tcga > 0 ~ "significant in both, consistent",
    log2fc < 0 & log2fc_tcga < 0 ~ "significant in both, consistent",
    T ~ "significant in both, inconsistent")) 

count <- df %>% count(sig.status)
#df <- df %>% filter(sig.status!="not significant")
cor = cor.test(df$log2fc, df$log2fc_tcga, method = "spearman")
cor_p = format(cor[["p.value"]], scientific = TRUE)
cor_rho = cor[["estimate"]][["rho"]]

df %>% 
  ggplot(aes(x=log2fc,
             y=log2fc_tcga,
             color=sig.status))+
  geom_point(alpha=0.8)+
  scale_color_manual(values = c( "significant only in Metastatic samples" = "#DDCC77",
                                 "significant only in Primary samples" = "#88CCEE",
                                 "significant in both, consistent" = "#332288",
                                 "significant in both, inconsistent" = "red",
                                 "not significant" = 'grey'))+
  #scale_x_continuous(breaks = c(-0.5,0,0.5),limits = c(-0.5,0.5))+
  #scale_y_continuous(breaks = c(-0.5,0,0.5),limits = c(-0.5,0.5))+
  #coord_equal()+
  theme_Publication()+
  labs(x="Primary Responder/Non-Responder log2(fold change)",
       y="Metastasis Responder/Non-Responder log2(fold change)",
       title = "Braunetal Metastasis vs Primary in Sunitinib arms")+
  geom_text(aes(x = 5, y = -5, label = sprintf("Spearman's rho = %.*f,\n p = %s", 3, cor_rho, cor_p )))
ggsave(file.path(paste0(results_dir,"scatterplot_immotion151_C.pdf")),width = 10,height = 10) 

cor.test(df$log2fc, df$log2fc_tcga, method = "spearman")

####### 
# GSEA Hallmark RNA signature analysis of consistently overexpressed genes in primary/mets
results_dir =  "/work/reznik/xiea1/MIRTH/revision/results_lennert_project/immotion151/primary_met/"
tcga_dir = "/work/reznik/xiea1/MIRTH/revision/results_lennert_project/immotion151/primary_met/"
cptac_dir = "/work/reznik/xiea1/MIRTH/revision/results_lennert_project/braunetal/primary_met/"

cptac <- read.csv(paste0(cptac_dir,"differential_genes_statistic.csv"), row.names = 1)
tcga <- read.csv(paste0(tcga_dir,"differential_genes_statistic.csv"), row.names = 1)
tcga <- dplyr::rename(tcga, log2fc_tcga = log2fc)
tcga <- dplyr::rename(tcga, statistic_tcga = statistic)
tcga <- dplyr::rename(tcga, p_tcga = p)
tcga <- dplyr::rename(tcga, p_adj_tcga = p_adj)
df <- merge(cptac, tcga, by = "row.names")
colnames(df)[1] <- "gene"
df <- df %>% 
  mutate(sig.status = case_when(
    p_adj_tcga >= 0.1 & p_adj >= 0.1 ~ "not significant",
    p_adj_tcga < 0.1 & p_adj >= 0.1   ~ "significant only in IMmotion151",
    p_adj_tcga >= 0.1 & p_adj < 0.1 ~ "significant only in Braunetal",
    log2fc > 0 & log2fc_tcga > 0 ~ "significant in both, consistent",
    log2fc < 0 & log2fc_tcga < 0 ~ "significant in both, consistent",
    T ~ "significant in both, inconsistent")) 
count <- df %>% count(sig.status)
df <- df %>% filter(sig.status=="significant in both, consistent")
df$statistic_ave <- rowMeans(df[c("statistic", "statistic_tcga")])

# GSEA
library(msigdbr)
library(fgsea)

# (1)Load hallmark gene set names
# get a list of gene sets to check
hallmark<- msigdbr(species = "Homo sapiens", category = "H")
hallmark_gene_sets <- split(hallmark$gene_symbol, hallmark$gs_name)
# get a Named vector of gene-level stats
rna_stat<- df$statistic_ave
names(rna_stat) <- df$gene
# Run GSEA
fgseaRes_hallmark_cont <- fgseaMultilevel(pathways = hallmark_gene_sets, 
                                          stats    = rna_stat)
fgseaRes_hallmark_cont$sig <- ifelse(fgseaRes_hallmark_cont$padj < 0.1, "sig","not sig")
fgseaRes_hallmark_cont$pathway <- gsub("HALLMARK_","",fgseaRes_hallmark_cont$pathway)
fgseaRes_hallmark_cont$pathway <- gsub("_"," ",fgseaRes_hallmark_cont$pathway)
pathway_order <- fgseaRes_hallmark_cont[order(fgseaRes_hallmark_cont$NES),1]
fgseaRes_hallmark_cont$pathway <- factor(fgseaRes_hallmark_cont$pathway, levels = pathway_order$pathway)

fgseaRes_hallmark_cont_sig <- fgseaRes_hallmark_cont[fgseaRes_hallmark_cont$padj < 0.1,]
ggplot(fgseaRes_hallmark_cont_sig,aes(x=NES,y=pathway)) + geom_bar(stat="identity") + theme_Publication()
ggsave(file.path(cptac_dir,"hallmark_gsea_ave.pdf"),width = 4,height = 4)

fgseaRes_hallmark_cont$leadingEdge <- NULL  # Remove the list column
write.csv(fgseaRes_hallmark_cont, file.path(cptac_dir,"hallmark_gsea_ave.csv"))


# (2) Load Reactome gene set names
# get a list of gene sets to check
reactome<- msigdbr(species = "Homo sapiens", category = "C2", subcategory = 'CP:REACTOME')
reactome_gene_sets <- split(reactome$gene_symbol, reactome$gs_name)
# get a Named vector of gene-level stats
rna_stat<- df$statistic_ave
names(rna_stat) <- df$gene
# Run GSEA
fgseaRes_hallmark_cont <- fgseaMultilevel(pathways = reactome_gene_sets, 
                                          stats    = rna_stat)
fgseaRes_hallmark_cont$sig <- ifelse(fgseaRes_hallmark_cont$padj < 0.1, "sig","not sig")
fgseaRes_hallmark_cont$pathway <- gsub("REACTOME_","",fgseaRes_hallmark_cont$pathway)
fgseaRes_hallmark_cont$pathway <- gsub("_"," ",fgseaRes_hallmark_cont$pathway)
pathway_order <- fgseaRes_hallmark_cont[order(fgseaRes_hallmark_cont$NES),1]
fgseaRes_hallmark_cont$pathway <- factor(fgseaRes_hallmark_cont$pathway, levels = pathway_order$pathway)

fgseaRes_hallmark_cont_sig <- fgseaRes_hallmark_cont[fgseaRes_hallmark_cont$padj < 0.1,]
ggplot(fgseaRes_hallmark_cont_sig,aes(x=NES,y=pathway)) + geom_bar(stat="identity") + theme_Publication()
ggsave(file.path(cptac_dir,"reactome_gsea_ave.pdf"),width = 10,height = 4)

fgseaRes_hallmark_cont$leadingEdge <- NULL  # Remove the list column
write.csv(fgseaRes_hallmark_cont, file.path(cptac_dir,"reactome_gsea_ave.csv"))



# previous records
# Match gene names with gene set names
df$pathway <- NA
for (i in 1: nrow(df)){
  if (df[i,]$gene %in% hallmark$gene_symbol) {
    df[i,]$pathway <- hallmark[hallmark$gene_symbol==df[i,]$gene,'gs_name']}
  else {next}
}
df <- df[!is.na(df$pathway),] # 152 genes
rna_stat<- df$statistic_tcga
names(rna_stat) <- df$pathway
hallmark_gene_sets <- as.list(unique(df$pathway))
names(hallmark_gene_sets) <- hallmark_gene_sets

