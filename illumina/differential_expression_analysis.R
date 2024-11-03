

library(edgeR)
library(readr)
library(ggplot2)
library(ggrepel)
library(xtable)
library(dplyr)

# Read the raw count data from the CSV file
file_path <- "~/Downloads/normalized_sorted_Amerge1.csv"
data <- read_delim(file_path, delim = ";", col_names = FALSE, skip = 1)

# Assign column names as specified
colnames(data) <- c("chr", "start", "end", "A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9")

# Create a count matrix by selecting only the count columns
count_data <- data %>%
  select(A1:A9) %>%
  as.data.frame()

# Create a unique circRNA ID
data <- data %>%
  mutate(circRNA_ID = paste(chr, start, end, sep = "_"))

# Assign row names to the count matrix based on circRNA ID
rownames(count_data) <- data$circRNA_ID

# Create a sample metadata table with "condition" and "time"
sample_metadata <- data.frame(
  row.names = colnames(count_data),
  condition = rep(c("control", "cmv", "zika"), times = 3),
  time = rep(c("24h", "48h", "72h"), each = 3)
)

# Convert the condition and time columns to factors
sample_metadata$condition <- as.factor(sample_metadata$condition)
sample_metadata$time <- as.factor(sample_metadata$time)



# Loop per ogni punto temporale
for (time_point in unique_times) {
  # Filtra i campioni per il punto temporale corrente
  time_samples <- sample_metadata %>%
    filter(time == time_point)
  count_data_time <- count_data[, rownames(time_samples)]
  
  # Crea un DGEList per il tempo specifico
  y_time <- DGEList(counts = count_data_time)
  
  # Normalizzazione TMM
  # y_time <- calcNormFactors(y_time)
  
  # Stima della dispersione
  y_time <- estimateDisp(y_time)
  
  # Creazione della matrice di design in base alla condizione
  design_condition_time <- model.matrix(~ condition, data = time_samples)
  
  # Adattamento del modello
  fit_condition_time <- glmFit(y_time, design_condition_time)
  
  # Test del rapporto di verosimiglianza per i confronti tra condizioni
  lrt_cmv_vs_control <- glmLRT(fit_condition_time, coef = 2)  # CMV vs Control
  lrt_zika_vs_control <- glmLRT(fit_condition_time, coef = 3)  # Zika vs Control
  
  # Salva i risultati in variabili dinamiche per il punto temporale corrente
  assign(paste0("results_cmv_vs_control_", time_point), topTags(lrt_cmv_vs_control, n = nrow(y_time))$table)
  assign(paste0("results_zika_vs_control_", time_point), topTags(lrt_zika_vs_control, n = nrow(y_time))$table)
  
  # Generazione dei grafici Volcano
  # CMV vs Control
  volcano_data_cmv_vs_control <- get(paste0("results_cmv_vs_control_", time_point)) %>%
    mutate(circRNA_ID = rownames(get(paste0("results_cmv_vs_control_", time_point)))) %>%
    mutate(Significant = ifelse(FDR < 0.05 & abs(logFC) > 1, "Yes", "No"))
  
  print(
    ggplot(volcano_data_cmv_vs_control, aes(x = logFC, y = -log10(PValue), color = Significant)) +
      geom_point() +
      geom_text_repel(data = filter(volcano_data_cmv_vs_control, Significant == "Yes"), aes(label = circRNA_ID)) +
      labs(title = paste("Volcano Plot: CMV vs Control at", time_point), x = "Log2 Fold Change", y = "-Log10 P-value") +
      scale_color_manual(values = c("No" = "red", "Yes" = "blue")) +
      theme_minimal()
  )
  
  ggsave(paste0("/Users/martina/Documents/Template_tesi_UniPD_DEI/Immagini/Volcano_CMV_vs_Control_", time_point, ".png"),
         width = 14, height = 8, dpi = 300)
  
  # Zika vs Control
  volcano_data_zika_vs_control <- get(paste0("results_zika_vs_control_", time_point)) %>%
    mutate(circRNA_ID = rownames(get(paste0("results_zika_vs_control_", time_point)))) %>%
    mutate(Significant = ifelse(FDR < 0.05 & abs(logFC) > 1, "Yes", "No"))
  
  ggplot(volcano_data_zika_vs_control, aes(x = logFC, y = -log10(PValue), color = Significant)) +
    geom_point() +
    geom_text_repel(data = filter(volcano_data_zika_vs_control, Significant == "Yes"), aes(label = circRNA_ID)) +
    labs(title = paste("Volcano Plot: Zika vs Control at", time_point), x = "Log2 Fold Change", y = "-Log10 P-value") +
    scale_color_manual(values = c("No" = "red", "Yes" = "blue")) +
    theme_minimal()
  
  ggsave(paste0("/Users/martina/Documents/Template_tesi_UniPD_DEI/Immagini/Volcano_Zika_vs_Control_", time_point, ".png"),
         width = 14, height = 8, dpi = 300)
}

# Funzione per esportare i top circRNA in LaTeX
for (time_point in unique_times) {
  # Esporta i risultati per CMV vs Control
  export_top_circRNA_to_LaTeX(get(paste0("results_cmv_vs_control_", time_point)), "CMV_vs_Control", time_point, top_n = 20)
  
  # Esporta i risultati per Zika vs Control
  export_top_circRNA_to_LaTeX(get(paste0("results_zika_vs_control_", time_point)), "Zika_vs_Control", time_point, top_n = 20)
}
























# Filtra i campioni per Control
control_samples <- sample_metadata %>%
  filter(condition == "control")

# Filtra i dati di conteggio per i campioni Control
count_data_control <- count_data[, rownames(control_samples)]

# Crea un oggetto DGEList per i campioni Control
y_control <- DGEList(counts = count_data_control)

# Normalizzazione (TMM)
# y_control <- calcNormFactors(y_control)

# Stima della dispersione
y_control <- estimateDisp(y_control)

# Crea la matrice di design basata solo sul tempo per i campioni Control
design_time_control <- model.matrix(~ time, data = control_samples)

# Fitting del modello
fit_time_control <- glmFit(y_control, design_time_control)

# Confronto tra 48h e 24h nei campioni Control (coef = 2)
lrt_control_48h_vs_24h <- glmLRT(fit_time_control, coef = 2)

# Confronto tra 72h e 24h nei campioni Control (coef = 3)
lrt_control_72h_vs_24h <- glmLRT(fit_time_control, coef = 3)

# Confronto tra 72h e 48h nei campioni Control
lrt_control_72h_vs_48h <- glmLRT(fit_time_control, contrast = c(0, -1, 1))

# Visualizza i risultati
topTags(lrt_control_48h_vs_24h)
topTags(lrt_control_72h_vs_24h)
topTags(lrt_control_72h_vs_48h)





# Risultati per 48h vs 24h
results_control_48h_vs_24h <- topTags(lrt_control_48h_vs_24h, n = nrow(y_control))$table
volcano_data_control_48h_vs_24h <- results_control_48h_vs_24h %>%
  mutate(circRNA_ID = rownames(results_control_48h_vs_24h)) %>%
  mutate(Significant = ifelse(FDR < 0.05 & abs(logFC) > 1, "Yes", "No"))

# Volcano Plot per 48h vs 24h
ggplot(volcano_data_control_48h_vs_24h, aes(x = logFC, y = -log10(PValue), color = Significant)) +
  geom_point() +
  geom_text_repel(data = filter(volcano_data_control_48h_vs_24h, Significant == "Yes"), aes(label = circRNA_ID)) +
  labs(title = "Volcano Plot: Control 48h vs 24h", x = "Log2 Fold Change", y = "-Log10 P-value") +
  scale_color_manual(values = c("No" = "red", "Yes" = "blue")) +
  theme_minimal()

# Salva il plot
ggsave("/Users/martina/Documents/Template_tesi_UniPD_DEI/Immagini/Volcano_Control_48h_vs_24h.png", 
       width = 14, height = 8, dpi = 300)






# Risultati per 72h vs 24h
results_control_72h_vs_24h <- topTags(lrt_control_72h_vs_24h, n = nrow(y_control))$table
volcano_data_control_72h_vs_24h <- results_control_72h_vs_24h %>%
  mutate(circRNA_ID = rownames(results_control_72h_vs_24h)) %>%
  mutate(Significant = ifelse(FDR < 0.05 & abs(logFC) > 1, "Yes", "No"))

# Volcano Plot per 72h vs 24h
ggplot(volcano_data_control_72h_vs_24h, aes(x = logFC, y = -log10(PValue), color = Significant)) +
  geom_point() +
  geom_text_repel(data = filter(volcano_data_control_72h_vs_24h, Significant == "Yes"), aes(label = circRNA_ID)) +
  labs(title = "Volcano Plot: Control 72h vs 24h", x = "Log2 Fold Change", y = "-Log10 P-value") +
  scale_color_manual(values = c("No" = "red", "Yes" = "blue")) +
  theme_minimal()

# Salva il plot
ggsave("/Users/martina/Documents/Template_tesi_UniPD_DEI/Immagini/Volcano_Control_72h_vs_24h.png", 
       width = 14, height = 8, dpi = 300)





# Risultati per 72h vs 48h
results_control_72h_vs_48h <- topTags(lrt_control_72h_vs_48h, n = nrow(y_control))$table
volcano_data_control_72h_vs_48h <- results_control_72h_vs_48h %>%
  mutate(circRNA_ID = rownames(results_control_72h_vs_48h)) %>%
  mutate(Significant = ifelse(FDR < 0.05 & abs(logFC) > 1, "Yes", "No"))

# Volcano Plot per 72h vs 48h
ggplot(volcano_data_control_72h_vs_48h, aes(x = logFC, y = -log10(PValue), color = Significant)) +
  geom_point() +
  geom_text_repel(data = filter(volcano_data_control_72h_vs_48h, Significant == "Yes"), aes(label = circRNA_ID)) +
  labs(title = "Volcano Plot: Control 72h vs 48h", x = "Log2 Fold Change", y = "-Log10 P-value") +
  scale_color_manual(values = c("No" = "red", "Yes" = "blue")) +
  theme_minimal()

# Salva il plot
ggsave("/Users/martina/Documents/Template_tesi_UniPD_DEI/Immagini/Volcano_Control_72h_vs_48h.png", 
       width = 14, height = 8, dpi = 300)
















# Filtra i campioni per CMV
cmv_samples <- sample_metadata %>%
  filter(condition == "cmv")

# Filtra i dati di conteggio per i campioni CMV
count_data_cmv <- count_data[, rownames(cmv_samples)]

# Crea un oggetto DGEList per i campioni CMV
y_cmv <- DGEList(counts = count_data_cmv)

# Normalizzazione (TMM)
# y_cmv <- calcNormFactors(y_cmv)

# Stima della dispersione
y_cmv <- estimateDisp(y_cmv)

# Crea la matrice di design basata solo sul tempo per i campioni CMV
design_time_cmv <- model.matrix(~ time, data = cmv_samples)

# Fitting del modello
fit_time_cmv <- glmFit(y_cmv, design_time_cmv)

# Confronto tra 48h e 24h nei campioni CMV (coef = 2)
lrt_cmv_48h_vs_24h <- glmLRT(fit_time_cmv, coef = 2)

# Confronto tra 72h e 24h nei campioni CMV (coef = 3)
lrt_cmv_72h_vs_24h <- glmLRT(fit_time_cmv, coef = 3)

# Confronto tra 72h e 48h nei campioni CMV (usando contrast)
lrt_cmv_72h_vs_48h <- glmLRT(fit_time_cmv, contrast = c(0, -1, 1))

# Visualizza i risultati
topTags(lrt_cmv_48h_vs_24h)
topTags(lrt_cmv_72h_vs_24h)
topTags(lrt_cmv_72h_vs_48h)




# Risultati per 48h vs 24h
results_cmv_48h_vs_24h <- topTags(lrt_cmv_48h_vs_24h, n = nrow(y_cmv))$table
volcano_data_cmv_48h_vs_24h <- results_cmv_48h_vs_24h %>%
  mutate(circRNA_ID = rownames(results_cmv_48h_vs_24h)) %>%
  mutate(Significant = ifelse(FDR < 0.05 & abs(logFC) > 1, "Yes", "No"))

# Volcano Plot per 48h vs 24h
ggplot(volcano_data_cmv_48h_vs_24h, aes(x = logFC, y = -log10(PValue), color = Significant)) +
  geom_point() +
  geom_text_repel(data = filter(volcano_data_cmv_48h_vs_24h, Significant == "Yes"), aes(label = circRNA_ID)) +
  labs(title = "Volcano Plot: CMV 48h vs 24h", x = "Log2 Fold Change", y = "-Log10 P-value") +
  scale_color_manual(values = c("No" = "red", "Yes" = "blue")) +
  theme_minimal()

# Salva il plot
ggsave("/Users/martina/Documents/Template_tesi_UniPD_DEI/Immagini/Volcano_CMV_48h_vs_24h.png", 
       width = 14, height = 8, dpi = 300)



# Risultati per 72h vs 24h
results_cmv_72h_vs_24h <- topTags(lrt_cmv_72h_vs_24h, n = nrow(y_cmv))$table
volcano_data_cmv_72h_vs_24h <- results_cmv_72h_vs_24h %>%
  mutate(circRNA_ID = rownames(results_cmv_72h_vs_24h)) %>%
  mutate(Significant = ifelse(FDR < 0.05 & abs(logFC) > 1, "Yes", "No"))

# Volcano Plot per 72h vs 24h
ggplot(volcano_data_cmv_72h_vs_24h, aes(x = logFC, y = -log10(PValue), color = Significant)) +
  geom_point() +
  geom_text_repel(data = filter(volcano_data_cmv_72h_vs_24h, Significant == "Yes"), aes(label = circRNA_ID)) +
  labs(title = "Volcano Plot: CMV 72h vs 24h", x = "Log2 Fold Change", y = "-Log10 P-value") +
  scale_color_manual(values = c("No" = "red", "Yes" = "blue")) +
  theme_minimal()

# Salva il plot
ggsave("/Users/martina/Documents/Template_tesi_UniPD_DEI/Immagini/Volcano_CMV_72h_vs_24h.png", 
       width = 14, height = 8, dpi = 300)






# Risultati per 72h vs 48h
results_cmv_72h_vs_48h <- topTags(lrt_cmv_72h_vs_48h, n = nrow(y_cmv))$table
volcano_data_cmv_72h_vs_48h <- results_cmv_72h_vs_48h %>%
  mutate(circRNA_ID = rownames(results_cmv_72h_vs_48h)) %>%
  mutate(Significant = ifelse(FDR < 0.05 & abs(logFC) > 1, "Yes", "No"))

# Volcano Plot per 72h vs 48h
ggplot(volcano_data_cmv_72h_vs_48h, aes(x = logFC, y = -log10(PValue), color = Significant)) +
  geom_point() +
  geom_text_repel(data = filter(volcano_data_cmv_72h_vs_48h, Significant == "Yes"), aes(label = circRNA_ID)) +
  labs(title = "Volcano Plot: CMV 72h vs 48h", x = "Log2 Fold Change", y = "-Log10 P-value") +
  scale_color_manual(values = c("No" = "red", "Yes" = "blue")) +
  theme_minimal()

# Salva il plot
ggsave("/Users/martina/Documents/Template_tesi_UniPD_DEI/Immagini/Volcano_CMV_72h_vs_48h.png", 
       width = 14, height = 8, dpi = 300)
















# Filtra i campioni per Zika
zika_samples <- sample_metadata %>%
  filter(condition == "zika")

# Filtra i dati di conteggio per i campioni Zika
count_data_zika <- count_data[, rownames(zika_samples)]

# Crea un oggetto DGEList per i campioni Zika
y_zika <- DGEList(counts = count_data_zika)

# Normalizzazione (TMM)
# y_zika <- calcNormFactors(y_zika)

# Stima della dispersione
y_zika <- estimateDisp(y_zika)

# Crea la matrice di design basata solo sul tempo per i campioni Zika
design_time_zika <- model.matrix(~ time, data = zika_samples)

# Fitting del modello
fit_time_zika <- glmFit(y_zika, design_time_zika)

# Confronto tra 48h e 24h nei campioni Zika (coef = 2)
lrt_zika_48h_vs_24h <- glmLRT(fit_time_zika, coef = 2)

# Confronto tra 72h e 24h nei campioni Zika (coef = 3)
lrt_zika_72h_vs_24h <- glmLRT(fit_time_zika, coef = 3)

# Confronto tra 72h e 48h nei campioni Zika
lrt_zika_72h_vs_48h <- glmLRT(fit_time_zika, contrast = c(0, -1, 1))

# Visualizza i risultati
topTags(lrt_zika_48h_vs_24h)
topTags(lrt_zika_72h_vs_24h)
topTags(lrt_zika_72h_vs_48h)


# Risultati per 48h vs 24h
results_zika_48h_vs_24h <- topTags(lrt_zika_48h_vs_24h, n = nrow(y_zika))$table
volcano_data_zika_48h_vs_24h <- results_zika_48h_vs_24h %>%
  mutate(circRNA_ID = rownames(results_zika_48h_vs_24h)) %>%
  mutate(Significant = ifelse(FDR < 0.05 & abs(logFC) > 1, "Yes", "No"))

# Volcano Plot per 48h vs 24h
ggplot(volcano_data_zika_48h_vs_24h, aes(x = logFC, y = -log10(PValue), color = Significant)) +
  geom_point() +
  geom_text_repel(data = filter(volcano_data_zika_48h_vs_24h, Significant == "Yes"), aes(label = circRNA_ID)) +
  labs(title = "Volcano Plot: Zika 48h vs 24h", x = "Log2 Fold Change", y = "-Log10 P-value") +
  scale_color_manual(values = c("No" = "red", "Yes" = "blue")) +
  theme_minimal()

# Salva il plot
ggsave("/Users/martina/Documents/Template_tesi_UniPD_DEI/Immagini/Volcano_Zika_48h_vs_24h.png", 
       width = 14, height = 8, dpi = 300)




# Risultati per 72h vs 24h
results_zika_72h_vs_24h <- topTags(lrt_zika_72h_vs_24h, n = nrow(y_zika))$table
volcano_data_zika_72h_vs_24h <- results_zika_72h_vs_24h %>%
  mutate(circRNA_ID = rownames(results_zika_72h_vs_24h)) %>%
  mutate(Significant = ifelse(FDR < 0.05 & abs(logFC) > 1, "Yes", "No"))

# Volcano Plot per 72h vs 24h
ggplot(volcano_data_zika_72h_vs_24h, aes(x = logFC, y = -log10(PValue), color = Significant)) +
  geom_point() +
  geom_text_repel(data = filter(volcano_data_zika_72h_vs_24h, Significant == "Yes"), aes(label = circRNA_ID)) +
  labs(title = "Volcano Plot: Zika 72h vs 24h", x = "Log2 Fold Change", y = "-Log10 P-value") +
  scale_color_manual(values = c("No" = "red", "Yes" = "blue")) +
  theme_minimal()

# Salva il plot
ggsave("/Users/martina/Documents/Template_tesi_UniPD_DEI/Immagini/Volcano_Zika_72h_vs_24h.png", 
       width = 14, height = 8, dpi = 300)



# Risultati per 72h vs 48h
results_zika_72h_vs_48h <- topTags(lrt_zika_72h_vs_48h, n = nrow(y_zika))$table
volcano_data_zika_72h_vs_48h <- results_zika_72h_vs_48h %>%
  mutate(circRNA_ID = rownames(results_zika_72h_vs_48h)) %>%
  mutate(Significant = ifelse(FDR < 0.05 & abs(logFC) > 1, "Yes", "No"))

# Volcano Plot per 72h vs 48h
ggplot(volcano_data_zika_72h_vs_48h, aes(x = logFC, y = -log10(PValue), color = Significant)) +
  geom_point() +
  geom_text_repel(data = filter(volcano_data_zika_72h_vs_48h, Significant == "Yes"), aes(label = circRNA_ID)) +
  labs(title = "Volcano Plot: Zika 72h vs 48h", x = "Log2 Fold Change", y = "-Log10 P-value") +
  scale_color_manual(values = c("No" = "red", "Yes" = "blue")) +
  theme_minimal()

# Salva il plot
ggsave("/Users/martina/Documents/Template_tesi_UniPD_DEI/Immagini/Volcano_Zika_72h_vs_48h.png", 
       width = 14, height = 8, dpi = 300)






