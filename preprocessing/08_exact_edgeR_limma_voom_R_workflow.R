# Exact edgeR + limma-voom workflow for the supplied files
# Run in an R environment with edgeR and limma installed.
suppressPackageStartupMessages({library(edgeR); library(limma)})

out_dir <- "limma_voom_preheatmap_output_R_exact"
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
raw <- read.csv("raw_counts_all_35_ordered(1).csv", check.names=FALSE, stringsAsFactors=FALSE)
meta <- read.csv("sample_metadata_all_35_ordered(1).csv", check.names=FALSE, stringsAsFactors=FALSE)
target <- read.csv("target_gene_annotation_167(1).csv", check.names=FALSE, stringsAsFactors=FALSE)

samples <- meta$sample_id
stopifnot(identical(colnames(raw)[-(1:2)], samples))
stopifnot(all(meta$cohort %in% c("Normal","Obese","Diabetic","HFrEF")))
meta$cohort <- factor(meta$cohort, levels=c("Normal","Obese","Diabetic","HFrEF"))
counts <- as.matrix(raw[, samples])
storage.mode(counts) <- "numeric"
rownames(counts) <- raw$gene_id

y <- DGEList(counts=counts, samples=meta, genes=raw[,c("gene_id","gene_symbol")])
design <- model.matrix(~0 + cohort, data=meta)
colnames(design) <- levels(meta$cohort)
keep <- filterByExpr(y, design=design)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y, method="TMM")

# QC and voom
v <- voom(y, design, plot=TRUE)
fit <- lmFit(v, design)
cont <- makeContrasts(
  Obese_vs_Normal=Obese-Normal,
  Diabetic_vs_Normal=Diabetic-Normal,
  HFrEF_vs_Normal=HFrEF-Normal,
  levels=design
)
fit2 <- eBayes(contrasts.fit(fit, cont))
for (cn in colnames(cont)) {
  tt <- topTable(fit2, coef=cn, number=Inf, sort.by="P")
  tt$gene_id <- y$genes$gene_id[match(rownames(tt), rownames(y))]
  tt$gene_symbol <- y$genes$gene_symbol[match(rownames(tt), rownames(y))]
  write.csv(tt, file.path(out_dir, paste0(cn,"_DE_results.csv")), row.names=FALSE)
}

# Target mapping / normal centering
canon <- function(x) toupper(gsub("\\s+", "", trimws(x)))
target$canon_symbol <- canon(target$gene_symbol)
y$genes$canon_symbol <- canon(y$genes$gene_symbol)
idx <- match(target$gene_id, y$genes$gene_id)
idx_sym <- match(target$canon_symbol, y$genes$canon_symbol)
idx[is.na(idx)] <- idx_sym[is.na(idx)]
map <- target
map$voom_row <- idx
map$retention_status <- ifelse(is.na(idx), "missing_or_filtered_low_expression", "retained")
write.csv(map, file.path(out_dir, "target_mapping_retention_audit.csv"), row.names=FALSE)

ret <- which(!is.na(idx))
Etarget <- v$E[idx[ret],,drop=FALSE]
rownames(Etarget) <- target$gene_symbol[ret]
normal_cols <- which(meta$cohort == "Normal")
normal_mean <- rowMeans(Etarget[,normal_cols,drop=FALSE])
centered <- Etarget - normal_mean
write.csv(Etarget, file.path(out_dir, "target_only_voom_E.csv"))
write.csv(data.frame(gene_symbol=rownames(Etarget), Normal_mean_voom_log2CPM=normal_mean), file.path(out_dir, "target_normal_means.csv"), row.names=FALSE)
write.csv(centered, file.path(out_dir, "final_Normal_centered_target_matrix.csv"))
stopifnot("TMEM94" %in% rownames(centered))
