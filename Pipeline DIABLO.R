# =============================================================================
# Pipeline DIABLO para integração multiômica supervisionada
#
# Input:
#   - X: matrizes transcriptômica, proteômica e de metilação
#   - Y: fenótipo (PR status)
#
# Output:
#   - Tuned DIABLO model
#   - Projeção das amostras
#   - Plots de correlação
#   - Importância de variáveis
#   - ROC curves
#   - Performance e validação cruzada
# =============================================================================

library(mixOmics)

# =============================================================================
# 1. MONTAGEM DO X E Y
# =============================================================================

X <- build_X(
  
  RNA = rna_z,
  Prot = prot_z,
  Meth = meth_z
  
)

Y <- factor(meta2[[outcome]])

# =============================================================================
# 2. DIABLO design matrix (Correlação fraca entre os blocos)
# =============================================================================

design <- matrix(
  0.1,
  nrow=3,
  ncol=3
)

diag(design) <- 0

rownames(design) <- colnames(design) <- names(X)

# =============================================================================
# 3. DATA PROPORTION
# =============================================================================

create_keepX <- function(X, proportions = c(0.02, 0.05, 0.10)) {
  
  keepX <- lapply(X, function(block){
    
    p <- ncol(block)
    vals <- round(p * proportions)
    vals <- unique(pmax(vals, 2))
    vals
    
  })
  
  keepX
  
}

test.keepX <- create_keepX(
  X,
  proportions = c(0.02, 0.05, 0.10)
)

# =============================================================================
# 4. TUNING
# =============================================================================

tune.diablo <- tune.block.splsda(
  X = X,
  Y = Y,
  ncomp = 2,
  design = design,
  test.keepX = test.keepX,
  validation = "Mfold",
  folds = 5,
  nrepeat = 3
)

tune.diablo$choice.keepX

# =============================================================================
# 5. MODELO FINAL
# =============================================================================

list.keepX <- tune.diablo$choice.keepX

diablo <- block.splsda(
  X = X,
  Y = Y,
  ncomp = 2,
  keepX = list.keepX,
  design = design
)


# =============================================================================
# DIABLO - PLOT E ANÁLISES
# =============================================================================

# =====================================
# 1. Projeção das amostras
# =====================================

plotIndiv(
  diablo,
  ind.names = FALSE,
  legend = TRUE,
  title = "DIABLO",
  pch = 16,
  cex = 1.5,
  star = FALSE,
  ellipse = FALSE,
  size.title = rel(1),
  size.subtitle = rel(1),
  size.legend = rel(1)
)


# =====================================
# 2. Correlação entre blocos
# =====================================

plotDiablo(diablo)


# =====================================
# 3. plot arrow (setas longas -> conflito entre blocos)
# =====================================

plotArrow(diablo, ind.names = FALSE)


# =====================================
# 4. Correlation Circle Plot
# =====================================

plotVar(
  diablo,
  comp = c(1, 2),   # obrigatório
  var.names = FALSE,
  col = c("#6c8ebf", "#7f9f7f", "#9a7fbf"),
  pch = c(16, 16, 16),  # Um por bloco
  cex = c(2, 2, 2),
  legend = TRUE
)

# =====================================
# 5. Circos plot
# =====================================
circosPlot(
  diablo,
  comp = 1,
  cutoff = 0.5,
  size.variables = 0.5,
  var.names = list(
    RNA  = gene_symbol[colnames(X$RNA)],
    Prot = colnames(X$Prot),
    Meth = colnames(X$Meth)
  )
)


# =====================================
# 6. Loading plot
# =====================================

plotLoadings(
  diablo,
  block = "RNA",
  name.var = gene_symbol_clean,
  comp = 1,
  contrib = "max",
  method = "mean",
  ndisplay = 20,
  size.name = 0.7,
  size.title = 1
)

# =====================================
# 7. ROC + AUC
# =====================================

library(pROC)

pred <- predict(diablo, newdata = X)

blocks <- c("RNA", "Prot", "Meth")
cols <- c("#6c8ebf", "#7f9f7f", "#9a7fbf")

roc_list <- lapply(blocks, function(b){
  prob <- pred$predict[[b]][,,1]
  roc(Y, prob[, "positive"])
})

names(roc_list) <- blocks

# ---- Plot base ----
op <- par(xpd = FALSE)

plot(
  1 - roc_list[[1]]$specificities,
  roc_list[[1]]$sensitivities,
  type = "l",
  col = cols[1],
  lwd = 3,
  xlab = "False Positive Rate",
  ylab = "True Positive Rate",
  main = "ROC Curve - DIABLO",
  xlim = c(0,1),
  ylim = c(0,1),
  bty = "n"
)

par(op)

# adicionar os outros blocos
for(i in 2:length(roc_list)){
  lines(
    1 - roc_list[[i]]$specificities,
    roc_list[[i]]$sensitivities,
    col = cols[i],
    lwd = 3
  )
}

# diagonal
abline(0, 1, lty = 2, col = "gray70", lwd = 2)

# grid leve
grid(col = "gray90")

# legenda com AUC
auc_vals <- sapply(roc_list, auc)

legend(
  x = 0.40, y = 0.3,
  legend = paste0(names(roc_list), " (AUC=", round(auc_vals, 3), ")"),
  col = cols,
  lwd = 3,
  bty = "n",
  cex = 0.8
)


# =====================================
# 8. Performance e validação cruzada
# =====================================

library(dplyr)
library(ggplot2)

perf.diablo <- perf(
  diablo,
  validation = "Mfold",
  folds = 5,
  nrepeat = 5
)

df <- lapply(names(perf.diablo$error.rate), function(block) {
  
  mat <- perf.diablo$error.rate[[block]]
  
  data.frame(
    block = block,
    comp  = 1:ncol(mat),
    error = colMeans(mat),
    sd    = apply(mat, 2, sd)
  )
  
}) %>% bind_rows()

ggplot(df, aes(x = comp, y = error, color = block)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  
  scale_color_manual(values = c(
    RNA = "#6c8ebf",
    Prot = "#7f9f7f",
    Meth = "#9a7fbf"
  )) +
  
  geom_errorbar(aes(ymin = error - sd, ymax = error + sd), width = 0.1) +
  theme_minimal(base_size = 12) +
  labs(
    x = "Number of components",
    y = "Error rate",
    color = "Omics",
    title = "DIABLO performance"
  )
