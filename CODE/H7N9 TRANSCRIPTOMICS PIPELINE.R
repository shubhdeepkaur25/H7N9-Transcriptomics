################################################################################
# COMBINED H7N9 TRANSCRIPTOMICS PIPELINE
# Full Analysis & Plotting Execution
################################################################################

## =====================================================
## 0. SETUP & INSTALL MISSING PACKAGES
## =====================================================
setwd("~/Documents/data analysis in immunology/data")

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

packages <- c(
  "ggplot2", "ggVennDiagram", "enrichR", "dplyr", "stringr", "tidyr", 
  "pheatmap", "VennDiagram", "grid", "limma", "ggrepel", 
  "FactoMineR", "factoextra", "ggvenn", "clusterProfiler", "org.Hs.eg.db"
)

for(pkg in packages) {
  if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
    if(pkg %in% c("clusterProfiler", "org.Hs.eg.db")) {
      BiocManager::install(pkg, update = FALSE)
    } else {
      install.packages(pkg)
    }
  }
}

library(ggplot2)
library(ggVennDiagram)
library(enrichR)
library(dplyr)
library(stringr)
library(tidyr)
library(pheatmap)
library(VennDiagram)
library(grid)
library(limma)
library(ggrepel)
library(FactoMineR)
library(factoextra)
library(ggvenn)
library(clusterProfiler)
library(org.Hs.eg.db)

theme_paper <- function() {
  theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      panel.grid.minor = element_blank(),
      legend.title = element_blank()
    )
}

## =====================================================
## 1. READ EXPRESSION MATRIX
## =====================================================
dat_raw <- read.delim(
  "Project-H7N9Human.txt",
  row.names = 1,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
expr_raw <- dat_raw
expr_raw[] <- lapply(expr_raw, function(x) as.numeric(as.character(x)))
rownames(expr_raw) <- make.unique(rownames(expr_raw))

# log2 matrix for t-tests
expr_log2 <- log2(expr_raw + 1)
rownames(expr_log2) <- rownames(expr_raw)

cat("Data preview (raw):\n")
print(head(expr_raw[, 1:4]))
cat("\nData preview (log2):\n")
print(head(expr_log2[, 1:4]))
cat("\n--------------------------------------------------\n")

## =====================================================
## 2. SAMPLE NAMING RULES & DEG THRESHOLDS
## =====================================================
virus_pattern <- c(
  H7N9 = "H7N9",
  H7N7 = "H7N7",
  H3N2 = "H3N2",
  H5N1 = "VN1203",
  mock = "mock"
)

logfc_cutoff <- 1
padj_cutoff <- 0.05

## =====================================================
## 3. VIRUS VS MOCK DEGs (ALL TIMEPOINTS)
## =====================================================
get_deg_vs_mock_all <- function(expr_raw, expr_log2, virus_key, timepoint,
                                logfc_cutoff = 1, padj_cutoff = 0.05) {
  virus_prefix <- virus_pattern[virus_key]
  mock_prefix <- virus_pattern["mock"]
  
  pat_virus <- paste0("^", virus_prefix, "_", timepoint)
  pat_mock <- paste0("^", mock_prefix, "_", timepoint)
  
  idx_virus <- grep(pat_virus, colnames(expr_raw))
  idx_mock <- grep(pat_mock, colnames(expr_raw))
  
  comp_name <- paste0(virus_key, " vs mock (", timepoint, ")")
  
  if (length(idx_virus) == 0 || length(idx_mock) == 0) {
    warning("Skipping: ", comp_name, " - samples not found.")
    return(NULL)
  }
  
  sub_raw <- expr_raw[, c(idx_virus, idx_mock), drop = FALSE]
  sub_log2 <- expr_log2[, c(idx_virus, idx_mock), drop = FALSE]
  group <- c(rep(virus_key, length(idx_virus)), rep("mock", length(idx_mock)))
  
  mean_virus_raw <- rowMeans(sub_raw[, group == virus_key, drop = FALSE])
  mean_mock_raw <- rowMeans(sub_raw[, group == "mock", drop = FALSE])
  log2FC <- log2(mean_virus_raw + 1) - log2(mean_mock_raw + 1)
  
  pvals <- apply(sub_log2, 1, function(x) {
    x1 <- x[group == virus_key]
    x2 <- x[group == "mock"]
    if (length(x1) < 2 || length(x2) < 2 || sd(x1) == 0 || sd(x2) == 0) {
      return(1.0)
    } else {
      out <- try(t.test(x1, x2)$p.value, silent = TRUE)
      if (inherits(out, "try-error") || is.na(out)) return(1.0)
      return(out)
    }
  })
  adj_p <- p.adjust(pvals, method = "BH")
  
  deg_all <- rownames(sub_raw)[adj_p < padj_cutoff & abs(log2FC) > logfc_cutoff]
  cat(comp_name, "| All DEGs:", length(deg_all), "\n")
  return(deg_all)
}

timepoints <- c("24h", "07h", "12h", "03h")
viruses_eval <- c("H7N9", "H3N2", "H7N7", "H5N1")
all_venn_lists <- list()

for (tp in timepoints) {
  cat("\n--- Analyzing DEGs for", tp, "---\n")
  venn_all <- list()
  for (v in viruses_eval) {
    res <- get_deg_vs_mock_all(
      expr_raw = expr_raw, expr_log2 = expr_log2,
      virus_key = v, timepoint = tp,
      logfc_cutoff = logfc_cutoff, padj_cutoff = padj_cutoff
    )
    if (!is.null(res)) {
      venn_all[[v]] <- res
    }
  }
  
  all_venn_lists[[tp]] <- venn_all
  cat("\nNumber of DEGs per virus at", tp, ":\n")
  print(sapply(venn_all, length))
  
  p_all <- ggVennDiagram(venn_all, label_alpha = 0, label = "count") +
    scale_fill_gradient(low = "white", high = "#7CAE00") +
    ggtitle(paste0("All DEGs at ", tp, " (Virus vs mock)"))
  print(p_all)
  
  out_png <- paste0("Venn_", tp, "_All_DEGs_H7N9_H3N2_H7N7_H5N1.png")
  ggsave(out_png, plot = p_all, width = 6, height = 6, dpi = 300)
  cat("Venn diagram saved as:", out_png, "\n")
}

## =====================================================
## 4. SPECIFIC GENES FOR EACH VIRUS (24H)
## =====================================================
venn_24h <- all_venn_lists[["24h"]]

specific_all <- lapply(names(venn_24h), function(v) {
  others <- unlist(venn_24h[names(venn_24h) != v], use.names = FALSE)
  setdiff(venn_24h[[v]], others)
})
names(specific_all) <- names(venn_24h)

cat("\nSpecific ALL-DEG counts at 24h:\n")
print(sapply(specific_all, length))

get_log2fc_vs_mock <- function(expr_raw, virus_key, timepoint) {
  virus_prefix <- virus_pattern[virus_key]
  mock_prefix <- virus_pattern["mock"]
  
  idx_virus <- grep(paste0("^", virus_prefix, "_", timepoint), colnames(expr_raw))
  idx_mock <- grep(paste0("^", mock_prefix, "_", timepoint), colnames(expr_raw))
  
  mean_virus <- rowMeans(expr_raw[, idx_virus, drop = FALSE])
  mean_mock <- rowMeans(expr_raw[, idx_mock, drop = FALSE])
  
  log2(mean_virus + 1) - log2(mean_mock + 1)
}

specific_up <- list()
specific_down <- list()

for (v in names(specific_all)) {
  lfc <- get_log2fc_vs_mock(expr_raw, virus_key = v, timepoint = "24h")
  genes <- specific_all[[v]]
  lfc_g <- lfc[genes]
  
  specific_up[[v]] <- names(lfc_g)[lfc_g > logfc_cutoff]
  specific_down[[v]] <- names(lfc_g)[lfc_g < -logfc_cutoff]
}

cat("\nSpecific UP counts:\n")
print(sapply(specific_up, length))
cat("\nSpecific DOWN counts:\n")
print(sapply(specific_down, length))

# Export H5N1 Specific UP genes for EnrichR
write.table(specific_up[["H5N1"]], "H5N1_SPECIFIC_UP_24h.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)

## =====================================================
## 5. ENRICHR ANALYSIS (PATHWAYS)
## =====================================================
filename <- "H5N1_SPECIFIC_UP_24h.txt"
if(file.exists(filename)){
  gene_list <- read.table(filename, header = FALSE, stringsAsFactors = FALSE)$V1
  gene_list <- unique(gene_list[gene_list != ""])
  cat(paste0("Successfully read ", length(gene_list), " genes.\n"))
  
  setEnrichrSite("Enrichr")
  target_dbs <- c("KEGG_2021_Human", "Reactome_2022")
  cat("Connecting to Enrichr for analysis...\n")
  enriched_results <- enrichr(gene_list, target_dbs)
  
  plot_bubble <- function(enrich_df, db_name, title_color = "red", top_n = 15) {
    data_sig <- enrich_df %>%
      filter(P.value < 0.05) %>%
      arrange(P.value) %>%
      head(top_n)
    
    if (nrow(data_sig) == 0) {
      message(paste("No significantly enriched pathways found in", db_name))
      return(NULL)
    }
    
    data_sig$Count <- as.numeric(sapply(strsplit(data_sig$Overlap, "/"), "[[", 1))
    data_sig$Term <- factor(data_sig$Term, levels = rev(data_sig$Term))
    data_sig$Term_short <- str_trunc(as.character(data_sig$Term), 60, "right")
    
    p <- ggplot(data_sig, aes(x = Odds.Ratio, y = Term, size = Count, color = -log10(P.value))) +
      geom_point(alpha = 0.8) +
      scale_size_continuous(range = c(3, 8)) +
      scale_color_gradient(low = "blue", high = title_color) +
      labs(
        title = paste(db_name, "- Top", top_n, "Pathways"),
        x = "Odds Ratio",
        y = NULL,
        color = "-log10(P-value)",
        size = "Gene Count"
      ) +
      theme_bw() +
      theme(
        axis.text.y = element_text(size = 10, color = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
      )
    return(p)
  }
  
  # 1. KEGG 2021 Human
  p_kegg <- plot_bubble(enriched_results[["KEGG_2021_Human"]], "KEGG 2021 Human", title_color = "red")
  if (!is.null(p_kegg)) {
    print(p_kegg)
    ggsave("KEGG_2021_Bubble.pdf", p_kegg, width = 8, height = 6)
    cat("KEGG image saved as KEGG_2021_Bubble.pdf\n")
  }
  
  # 2. Reactome (2022)
  p_reactome <- plot_bubble(enriched_results[["Reactome_2022"]], "Reactome", title_color = "purple")
  if (!is.null(p_reactome)) {
    print(p_reactome)
    ggsave("Reactome_Bubble.pdf", p_reactome, width = 8, height = 6)
    cat("Reactome image saved as Reactome_Bubble.pdf\n")
  }
}

## =====================================================
## 6. MEAN HEATMAPS: H7N9 VS OTHERS
## =====================================================
mean_mat <- function(mat) rowMeans(mat, na.rm = TRUE)

get_virus_means <- function(v_prefix) {
  cbind(
    `03h` = mean_mat(expr_raw[, grep(paste0(v_prefix, "_03h"), colnames(expr_raw)), drop=FALSE]),
    `07h` = mean_mat(expr_raw[, grep(paste0(v_prefix, "_07h"), colnames(expr_raw)), drop=FALSE]),
    `12h` = mean_mat(expr_raw[, grep(paste0(v_prefix, "_12h"), colnames(expr_raw)), drop=FALSE]),
    `24h` = mean_mat(expr_raw[, grep(paste0(v_prefix, "_24h"), colnames(expr_raw)), drop=FALSE])
  )
}

H3N2_mean <- get_virus_means("H3N2")
H7N9_mean <- get_virus_means("H7N9")
H7N7_mean <- get_virus_means("H7N7")
VN1203_mean <- get_virus_means("VN1203")
mock_mean <- get_virus_means("mock")

plot_comparison_heatmap <- function(mean_data_1, mean_data_2, v_name_1, v_name_2) {
  expr_compare <- cbind(mean_data_1, mean_data_2)
  expr_log <- log2(expr_compare + 1)
  expr_z <- t(scale(t(expr_log)))
  
  virus <- c(rep(v_name_1, 4), rep(v_name_2, 4))
  time <- c("03h", "07h", "12h", "24h", "03h", "07h", "12h", "24h")
  
  annotation_col <- data.frame(Virus = virus, Time = time)
  rownames(annotation_col) <- colnames(expr_z)
  
  pheatmap(
    expr_z,
    annotation_col = annotation_col,
    show_rownames = FALSE,
    cluster_cols = FALSE,
    main = paste("Heatmap –", v_name_1, "vs", v_name_2)
  )
}

plot_comparison_heatmap(H7N9_mean, H3N2_mean, "H7N9", "H3N2")
plot_comparison_heatmap(H7N9_mean, H7N7_mean, "H7N9", "H7N7")
plot_comparison_heatmap(H7N9_mean, VN1203_mean, "H7N9", "VN1203")
plot_comparison_heatmap(H7N9_mean, mock_mean, "H7N9", "mock")

## =====================================================
## 7. 24h VENN ONLY (NO MOCK) & EXPORTS
## =====================================================
comparisons_24h <- c("H3N2", "H7N7", "VN1203")
H7N9_cols_24h <- grep("^H7N9_24h", colnames(expr_raw))
fc_cutoff <- 1
p_cutoff <- 0.05

up_24h_list <- list()
down_24h_list <- list()

safe_ttest <- function(x, idx1, idx2) {
  x1 <- x[idx1]; x2 <- x[idx2]
  if (length(unique(c(x1, x2))) < 2) return(1)
  out <- try(t.test(x1, x2)$p.value, silent = TRUE)
  if (inherits(out, "try-error") || is.na(out)) return(1)
  out
}

for (virus in comparisons_24h) {
  virus_cols_24h <- grep(paste0("^", virus, "_24h"), colnames(expr_raw))
  
  if (length(H7N9_cols_24h) < 2 || length(virus_cols_24h) < 2) {
    message("SKIPPED: ", virus, " (not enough 24h samples)")
    next
  }
  
  mean_H7N9_24h <- rowMeans(expr_raw[, H7N9_cols_24h, drop = FALSE])
  mean_virus_24h <- rowMeans(expr_raw[, virus_cols_24h, drop = FALSE])
  log2FC <- log2(mean_H7N9_24h + 1) - log2(mean_virus_24h + 1)
  
  pvals <- apply(expr_raw, 1, safe_ttest, idx1 = H7N9_cols_24h, idx2 = virus_cols_24h)
  
  volcano_24h <- data.frame(gene = rownames(expr_raw), log2FC = log2FC, pvalue = pvals)
  volcano_24h$status <- "Not Sig"
  volcano_24h$status[volcano_24h$log2FC > fc_cutoff & volcano_24h$pvalue < p_cutoff] <- "Up"
  volcano_24h$status[volcano_24h$log2FC < -fc_cutoff & volcano_24h$pvalue < p_cutoff] <- "Down"
  
  up_24h_list[[virus]] <- subset(volcano_24h, status == "Up")
  down_24h_list[[virus]] <- subset(volcano_24h, status == "Down")
  
  message("H7N9_24h vs ", virus, "_24h — Up: ", nrow(up_24h_list[[virus]]), 
          " | Down: ", nrow(down_24h_list[[virus]]))
}

up_sets_24h_no_mock <- list(
  H7N9_vs_H3N2_24h = up_24h_list[["H3N2"]]$gene,
  H7N9_vs_H7N7_24h = up_24h_list[["H7N7"]]$gene,
  H7N9_vs_VN1203_24h = up_24h_list[["VN1203"]]$gene
)

down_sets_24h_no_mock <- list(
  H7N9_vs_H3N2_24h = down_24h_list[["H3N2"]]$gene,
  H7N9_vs_H7N7_24h = down_24h_list[["H7N7"]]$gene,
  H7N9_vs_VN1203_24h = down_24h_list[["VN1203"]]$gene
)

venn_up_24h_no_mock <- venn.diagram(
  x = up_sets_24h_no_mock, filename = NULL,
  main = "UP genes (24h): H7N9 vs influenza (no mock)",
  col = "black", fill = c("red", "blue", "green"),
  alpha = 0.4, cex = 1.2, cat.cex = 1.2
)

grid.newpage()
grid.draw(venn_up_24h_no_mock)
png("Venn_UP_24h_H7N9_vs_influenza_noMock.png", width = 2000, height = 2000, res = 300)
grid.draw(venn_up_24h_no_mock)
dev.off()

venn_down_24h_no_mock <- venn.diagram(
  x = down_sets_24h_no_mock, filename = NULL,
  main = "DOWN genes (24h): H7N9 vs influenza (no mock)",
  col = "black", fill = c("red", "blue", "green"),
  alpha = 0.4, cex = 1.2, cat.cex = 1.2
)

grid.newpage()
grid.draw(venn_down_24h_no_mock)
png("Venn_DOWN_24h_H7N9_vs_influenza_noMock.png", width = 2000, height = 2000, res = 300)
grid.draw(venn_down_24h_no_mock)
dev.off()

# Export each comparison separately
write.table(up_sets_24h_no_mock$H7N9_vs_H3N2_24h, "H7N9_24h_UP_vs_H3N2_noMock.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(up_sets_24h_no_mock$H7N9_vs_H7N7_24h, "H7N9_24h_UP_vs_H7N7_noMock.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(up_sets_24h_no_mock$H7N9_vs_VN1203_24h, "H7N9_24h_UP_vs_VN1203_noMock.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(down_sets_24h_no_mock$H7N9_vs_H3N2_24h, "H7N9_24h_DOWN_vs_H3N2_noMock.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(down_sets_24h_no_mock$H7N9_vs_H7N7_24h, "H7N9_24h_DOWN_vs_H7N7_noMock.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(down_sets_24h_no_mock$H7N9_vs_VN1203_24h, "H7N9_24h_DOWN_vs_VN1203_noMock.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Export the INTERSECTION
common_up_24h_no_mock <- Reduce(intersect, up_sets_24h_no_mock)
common_down_24h_no_mock <- Reduce(intersect, down_sets_24h_no_mock)

write.table(common_up_24h_no_mock, "H7N9_24h_common_UP_noMock.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(common_down_24h_no_mock, "H7N9_24h_common_DOWN_noMock.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)


## =====================================================
## 8. TRANSCRIPTOMICS: QC, PCA, LIMMA DE, VENN
## =====================================================
sample_names <- colnames(expr_log2)
metadata <- data.frame(
  sample = sample_names,
  virus = sapply(strsplit(sample_names, "_"), `[`, 1),
  time = sapply(strsplit(sample_names, "_"), `[`, 2),
  stringsAsFactors = FALSE
)
metadata$group <- paste(metadata$virus, metadata$time, sep = "_")
metadata$time <- factor(metadata$time, levels = c("03h", "07h", "12h", "24h"))
rownames(metadata) <- metadata$sample

virus_colors <- c("H3N2" = "#1F77B4", "H7N7" = "#FF7F0E", "H7N9" = "#D62728", "VN1203" = "#9467BD", "mock" = "#7F7F7F")
time_colors <- c("03h"="#FEE5D9", "07h"="#FCAE91", "12h"="#FB6A4A", "24h"="#CB181D")
h7n9_cols <- rownames(metadata)[metadata$virus == "H7N9"]

# QC Boxplot
pdf("QC_expression_distribution.pdf", width = 12, height = 6)
boxplot(expr_log2, las = 2, main = "Expression Distribution Across Samples", col = virus_colors[metadata$virus], cex.axis = 0.6)
dev.off()

# QC Correlation Heatmap
cor_matrix <- cor(expr_log2, use = "pairwise.complete.obs")
pdf("QC_sample_correlation.pdf", width = 10, height = 10)
pheatmap(
  cor_matrix,
  main = "Sample–Sample Correlation",
  annotation_col = metadata[, c("virus", "time"), drop = FALSE],
  annotation_colors = list(virus = virus_colors, time = time_colors),
  show_colnames = FALSE, show_rownames = FALSE
)
dev.off()

# PCA (H7N9 only)
meta_h7 <- metadata[h7n9_cols, , drop = FALSE]
pca_h7 <- PCA(t(expr_log2[, h7n9_cols, drop = FALSE]), scale.unit = TRUE, graph = FALSE)
pdf("H7N9_PCA_trajectory.pdf", width = 10, height = 8)
p_h7_pca <- fviz_pca_ind(
  pca_h7,
  habillage = meta_h7$time,
  addEllipses = TRUE,
  ellipse.type = "confidence",
  repel = TRUE,
  palette = time_colors,
  title = "H7N9 temporal trajectory (PCA)"
) + theme_paper()
print(p_h7_pca)
dev.off()

# Differential Expression (Limma)
design <- model.matrix(~0 + group, data = metadata)
colnames(design) <- gsub("^group", "", colnames(design))
fit <- lmFit(expr_log2, design)

contrasts_response <- makeContrasts(
  H7N9_vs_Mock_03h = H7N9_03h - mock_03h,
  H7N9_vs_Mock_07h = H7N9_07h - mock_07h,
  H7N9_vs_Mock_12h = H7N9_12h - mock_12h,
  H7N9_vs_Mock_24h = H7N9_24h - mock_24h,
  levels = design
)
fit_response <- contrasts.fit(fit, contrasts_response)
fit_response <- eBayes(fit_response)

response_results <- list(
  "03h" = topTable(fit_response, coef="H7N9_vs_Mock_03h", number=Inf),
  "07h" = topTable(fit_response, coef="H7N9_vs_Mock_07h", number=Inf),
  "12h" = topTable(fit_response, coef="H7N9_vs_Mock_12h", number=Inf),
  "24h" = topTable(fit_response, coef="H7N9_vs_Mock_24h", number=Inf)
)

pval_cutoff <- 0.05
logfc_cutoff <- 1

sig_genes <- list()
upregulated_genes <- list()
downregulated_genes <- list()

for(tp in names(response_results)) {
  tt <- response_results[[tp]]
  sig_genes[[tp]] <- rownames(tt)[tt$adj.P.Val < pval_cutoff & abs(tt$logFC) > logfc_cutoff]
  upregulated_genes[[tp]] <- rownames(tt)[tt$adj.P.Val < pval_cutoff & tt$logFC > logfc_cutoff]
  downregulated_genes[[tp]] <- rownames(tt)[tt$adj.P.Val < pval_cutoff & tt$logFC < -logfc_cutoff]
  
  write.csv(upregulated_genes[[tp]], paste0("H7N9_upregulated_genes_", tp, ".csv"), row.names = FALSE)
  write.csv(downregulated_genes[[tp]], paste0("H7N9_downregulated_genes_", tp, ".csv"), row.names = FALSE)
}

# DEG counts barplot
deg_counts <- sapply(response_results, function(tt) {
  c(
    Upregulated = sum(tt$adj.P.Val < pval_cutoff & tt$logFC > logfc_cutoff),
    Downregulated = sum(tt$adj.P.Val < pval_cutoff & tt$logFC < -logfc_cutoff)
  )
})
deg_df <- data.frame(
  timepoint = rep(colnames(deg_counts), each = 2),
  direction = rep(c("Upregulated", "Downregulated"), ncol(deg_counts)),
  count = as.vector(deg_counts)
)
deg_df$timepoint <- factor(deg_df$timepoint, levels = c("03h", "07h", "12h", "24h"))

pdf("H7N9_DEG_counts_by_timepoint.pdf", width = 10, height = 6)
p_deg <- ggplot(deg_df, aes(x=timepoint, y=count, fill=direction)) +
  geom_bar(stat="identity", position="dodge", width=0.7) +
  geom_text(aes(label=count), position=position_dodge(0.7), vjust=-0.4, size=4, fontface="bold") +
  scale_fill_manual(values=c("Upregulated"="#D62728", "Downregulated"="#1F77B4")) +
  labs(title="H7N9 differential expression over time",
       subtitle="H7N9 vs mock at each timepoint (FDR<0.05, |log2FC|>1)",
       x="Time", y="Number of DEGs") +
  theme_paper()
print(p_deg)
dev.off()

# Venn Overlap
pdf("H7N9_temporal_overlap_venn.pdf", width = 8, height = 8)
p_venn <- ggvenn(sig_genes, fill_color = c("#FEE5D9", "#FCAE91", "#FB6A4A", "#CB181D"), stroke_size = 0.5, set_name_size = 5)
print(p_venn)
dev.off()

## =====================================================
## 9. GO + KEGG ENRICHMENT + FACETED BARPLOTS
## =====================================================
run_go <- function(genes) {
  if(length(genes) < 5) return(NULL)
  enrichGO(gene = genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
           ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
}

run_kegg <- function(genes) {
  if(length(genes) < 5) return(NULL)
  conv <- suppressMessages(bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db))
  entrez <- unique(conv$ENTREZID)
  if(length(entrez) < 5) return(NULL)
  enrichKEGG(gene = entrez, organism = "hsa", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
}

go_results_up <- list(); go_results_down <- list()
kegg_results_up <- list(); kegg_results_down <- list()

for(tp in names(upregulated_genes)) {
  go_results_up[[tp]] <- run_go(upregulated_genes[[tp]])
  go_results_down[[tp]] <- run_go(downregulated_genes[[tp]])
  kegg_results_up[[tp]] <- run_kegg(upregulated_genes[[tp]])
  kegg_results_down[[tp]] <- run_kegg(downregulated_genes[[tp]])
}

# Save Enrichment Tables
for(tp in names(go_results_up)) {
  if(!is.null(go_results_up[[tp]]) && nrow(as.data.frame(go_results_up[[tp]]))>0)
    write.csv(as.data.frame(go_results_up[[tp]]), paste0("H7N9_GO_upregulated_", tp, ".csv"), row.names=FALSE)
  if(!is.null(go_results_down[[tp]]) && nrow(as.data.frame(go_results_down[[tp]]))>0)
    write.csv(as.data.frame(go_results_down[[tp]]), paste0("H7N9_GO_downregulated_", tp, ".csv"), row.names=FALSE)
  
  if(!is.null(kegg_results_up[[tp]]) && nrow(as.data.frame(kegg_results_up[[tp]]))>0)
    write.csv(as.data.frame(kegg_results_up[[tp]]), paste0("H7N9_KEGG_upregulated_", tp, ".csv"), row.names=FALSE)
  if(!is.null(kegg_results_down[[tp]]) && nrow(as.data.frame(kegg_results_down[[tp]]))>0)
    write.csv(as.data.frame(kegg_results_down[[tp]]), paste0("H7N9_KEGG_downregulated_", tp, ".csv"), row.names=FALSE)
}

prep_facet <- function(enrich_list, direction, max_terms=10) {
  out <- data.frame()
  for(tp in names(enrich_list)) {
    obj <- enrich_list[[tp]]
    if(!is.null(obj) && nrow(as.data.frame(obj))>0) {
      df <- as.data.frame(obj)
      df <- df[order(df$p.adjust), ][1:min(max_terms, nrow(df)), ]
      df$timepoint <- tp
      df$direction <- direction
      out <- rbind(out, df)
    }
  }
  out
}

# GO Faceted Plot
go_up_df <- prep_facet(go_results_up, "Upregulated", 10)
go_down_df <- prep_facet(go_results_down, "Downregulated", 10)
go_df <- rbind(go_up_df, go_down_df)

if (nrow(go_df) > 0) {
  go_df$timepoint <- factor(go_df$timepoint, levels=c("03h", "07h", "12h", "24h"))
  pdf("H7N9_GO_enrichment_faceted.pdf", width = 16, height = 12)
  p_go <- ggplot(go_df, aes(x=reorder(Description, -log10(p.adjust)), y=-log10(p.adjust), fill=direction)) +
    geom_bar(stat="identity", alpha=0.85) +
    geom_text(aes(label=Count), hjust=-0.2, size=2.7) +
    coord_flip() +
    facet_wrap(~timepoint + direction, scales="free_y", ncol=2) +
    scale_fill_manual(values=c("Upregulated"="#D62728", "Downregulated"="#1F77B4")) +
    labs(title="GO Biological Process enrichment across timepoints", x="", y="-log10(adjusted P)") +
    theme_paper() +
    theme(axis.text.y = element_text(size=8), strip.text = element_text(face="bold", size=11), legend.position="bottom")
  print(p_go)
  dev.off()
}

# KEGG Faceted Plot
k_up_df <- prep_facet(kegg_results_up, "Upregulated", 10)
k_down_df <- prep_facet(kegg_results_down, "Downregulated", 10)
k_df <- rbind(k_up_df, k_down_df)

if (nrow(k_df) > 0) {
  k_df$timepoint <- factor(k_df$timepoint, levels=c("03h", "07h", "12h", "24h"))
  pdf("H7N9_KEGG_enrichment_faceted.pdf", width = 16, height = 12)
  p_kegg <- ggplot(k_df, aes(x=reorder(Description, -log10(p.adjust)), y=-log10(p.adjust), fill=direction)) +
    geom_bar(stat="identity", alpha=0.85) +
    geom_text(aes(label=Count), hjust=-0.2, size=2.7) +
    coord_flip() +
    facet_wrap(~timepoint + direction, scales="free_y", ncol=2) +
    scale_fill_manual(values=c("Upregulated"="#D62728", "Downregulated"="#1F77B4")) +
    labs(title="KEGG pathway enrichment across timepoints", x="", y="-log10(adjusted P)") +
    theme_paper() +
    theme(axis.text.y = element_text(size=8), strip.text = element_text(face="bold", size=11), legend.position="bottom")
  print(p_kegg)
  dev.off()
}

## =====================================================
## 10. TEMPORAL PROGRESSION: 24h vs 03h
## =====================================================
contrasts_progression <- makeContrasts(
  H7N9_progression_24h_vs_03h = (H7N9_24h - mock_24h) - (H7N9_03h - mock_03h),
  levels = design
)
fit_prog <- contrasts.fit(fit, contrasts_progression)
fit_prog <- eBayes(fit_prog)

deg_prog <- topTable(fit_prog, coef="H7N9_progression_24h_vs_03h", number=Inf)
deg_prog$gene <- rownames(deg_prog)
deg_prog$group <- "Not Sig"
deg_prog$group[deg_prog$adj.P.Val < 0.05 & deg_prog$logFC > 1] <- "Upregulated"
deg_prog$group[deg_prog$adj.P.Val < 0.05 & deg_prog$logFC < -1] <- "Downregulated"
deg_prog$group <- factor(deg_prog$group, levels = c("Downregulated", "Not Sig", "Upregulated"))

label_met <- c("HK2", "PFKP", "G6PD", "PDHA1", "SDHA", "ACACA", "FASN", "DHODH")
label_imm <- c("ISG15", "MX1", "IFIT1", "STAT1", "OAS2", "IL6", "CXCL10", "CCL5")
label_genes <- unique(c(label_met, label_imm))
deg_prog$label <- ifelse(deg_prog$gene %in% label_genes, deg_prog$gene, NA)

p_prog <- ggplot(deg_prog, aes(x=logFC, y=-log10(adj.P.Val), color=group)) +
  geom_point(alpha=0.65, size=1.5) +
  geom_vline(xintercept=c(-1,1), linetype="dashed", color="grey50") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="grey50") +
  ggrepel::geom_label_repel(aes(label=label), size=4, fill="white", label.size=0.2, color="black", max.overlaps=Inf) +
  scale_color_manual(values=c("Upregulated"="#D62728", "Downregulated"="#1F77B4", "Not Sig"="grey80")) +
  labs(title="H7N9 progression: 24h vs 3h (mock-corrected)", x="log2 fold change", y="-log10 FDR") +
  theme_paper()

ggsave("H7N9_progression_24h_vs_03h_volcano.pdf", p_prog, width=10, height=8)
print(p_prog)

## =====================================================
## 11. MODULE SCORES AND HEATMAPS
## =====================================================
gene_sets <- list(
  ER_Glyco_Lys = c("HSPA5", "PDIA3", "CANX", "CALR", "SEC61A1", "SSR1", "STT3A", "STT3B", "LAMP1", "LAMP2", "CTSB", "CTSD", "CTSL", "CTSS", "ATP6V1A", "ATP6V0D1"),
  ISR = c("ATF4", "ATF3", "DDIT3", "PPP1R15A", "ASNS", "TRIB3", "CHAC1", "GADD45A"),
  Metabolism = c("HK1", "HK2", "PFKP", "PFKM", "GPI", "ALDOA", "GAPDH", "PGK1", "ENO1", "PKM", "LDHA", "G6PD", "PGD", "TKT", "TALDO1", "PDHA1", "CS", "IDH3A", "SDHA", "FH", "MDH2", "ACACA", "FASN", "SCD", "ELOVL6", "CAD", "DHODH", "PPAT", "ATIC", "IMPDH1", "IMPDH2"),
  MHC_I = c("HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "TAPBP", "ERAP1", "ERAP2", "PSMB8", "PSMB9", "NLRC5"),
  MHC_II = c("CIITA", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "CD74", "CTSB", "CTSD", "CTSL", "CTSS")
)

gene_sets_present <- lapply(gene_sets, function(gs) gs[gs %in% rownames(expr_log2)])
modules <- names(gene_sets_present)
expr_z <- t(scale(t(expr_log2)))

module_scores_mat <- sapply(modules, function(m) {
  g <- gene_sets_present[[m]]
  if(length(g) < 3) return(rep(NA, ncol(expr_z)))
  colMeans(expr_z[g, , drop = FALSE], na.rm = TRUE)
})

module_scores <- as.data.frame(module_scores_mat)
module_scores$sample <- colnames(expr_z)
module_scores <- module_scores %>% left_join(metadata, by=c("sample"="sample"))

pref_order <- c("Metabolism", "ER_Glyco_Lys", "MHC_I", "MHC_II", "ISR")
pref_order <- pref_order[pref_order %in% modules]

# A) Barplots (24h mean across replicates, viruses only)
df24_mean <- module_scores %>%
  filter(time == "24h", virus != "mock") %>%
  group_by(virus) %>%
  summarise(across(all_of(pref_order), mean, na.rm=TRUE), .groups="drop") %>%
  pivot_longer(cols = all_of(pref_order), names_to="module", values_to="score")

df24_mean$virus <- factor(df24_mean$virus, levels=c("H3N2", "H7N7", "H7N9", "VN1203"))
df24_mean$module <- factor(df24_mean$module, levels=pref_order)

p_bar24 <- ggplot(df24_mean, aes(x=virus, y=score, fill=virus)) +
  geom_bar(stat="identity", width=0.75) +
  facet_wrap(~module, scales="free_y", ncol=3) +
  scale_fill_manual(values=virus_colors) +
  labs(title="24h module scores (mean across replicates)", x="Virus", y="Mean module score") +
  theme_paper() +
  theme(strip.text=element_text(face="bold"), legend.position="right")

ggsave("Fig_ModuleBarplot_24h_MeanAcrossReplicates.pdf", p_bar24, width=12, height=7)

# B) Summary heatmap (Δ vs mock at same time; 24h; scaled)
ms_delta <- module_scores %>%
  group_by(time) %>%
  mutate(across(all_of(pref_order), ~ . - mean(.[virus=="mock"], na.rm=TRUE), .names="{.col}_delta")) %>%
  ungroup()

ms24 <- ms_delta %>% filter(time=="24h", virus!="mock")
mat24 <- ms24 %>%
  group_by(virus) %>%
  summarise(across(ends_with("_delta"), mean, na.rm=TRUE), .groups="drop") %>%
  as.data.frame()

rownames(mat24) <- mat24$virus
mat24$virus <- NULL
mat24 <- mat24[, paste0(pref_order, "_delta"), drop=FALSE]
colnames(mat24) <- gsub("_delta", "", colnames(mat24))

mat24_z <- t(scale(t(as.matrix(mat24))))

pdf("Fig_ModuleHeatmap_24h_DeltaVsMock_Scaled.pdf", width=7.5, height=4.5)
pheatmap(mat24_z, cluster_rows=FALSE, cluster_cols=FALSE, main="24h module pattern (Δ vs mock, scaled)", fontsize=12)
dev.off()

# C) Gene heatmaps per module (24h)
create_module_gene_heatmap <- function(module_name, genes, timepoint="24h") {
  genes <- genes[genes %in% rownames(expr_log2)]
  if(length(genes) < 3) return(NULL)
  
  samples_tp <- metadata %>% filter(time==timepoint) %>% arrange(virus) %>% pull(sample)
  expr_sub <- expr_log2[genes, samples_tp, drop=FALSE]
  
  annot_col <- metadata[samples_tp, c("virus"), drop=FALSE]
  annot_colors <- list(virus = virus_colors)
  
  pheatmap(
    expr_sub, scale="row", cluster_cols=FALSE, cluster_rows=TRUE,
    show_colnames=FALSE, border_color=NA,
    annotation_col=annot_col, annotation_colors=annot_colors,
    fontsize_row=8, main=paste0(module_name, " genes at ", timepoint)
  )
}

pdf("Fig_ModuleGeneHeatmaps_24h.pdf", width=12, height=10)
for(m in pref_order) {
  create_module_gene_heatmap(m, gene_sets_present[[m]], "24h")
}
dev.off()

cat("\n=== DONE: all figures written to working directory ===\n")