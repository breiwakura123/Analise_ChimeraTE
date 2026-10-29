#install.packages("devtools")
library(devtools)
#install_github("chrisamiller/fishplot")


library(fishplot)
library(dplyr)


# ============================================================
# CONFIGURAÇÃO
# ============================================================

pasta <- "/home/brenda/Documentos/Brenda/An-lise_chimeraTE/Dados/fishplot_dir"  # ajuste aqui

# Nomes das colunas dos seus arquivos
col_names <- c("gene_id", "gene_strand", "gene_pos", "TE_id", "TE_strand", 
               "TE_pos", "chim_reads", "X8")

# Pacientes disponíveis (extraído automaticamente dos arquivos)
arquivos_pre <- list.files(pasta, pattern = "_TE-exonized_embedded_final_pre.tsv")
pacientes    <- gsub("_TE-exonized_embedded_final_pre.tsv", "", arquivos_pre)

cat("Pacientes encontrados:", paste(pacientes, collapse = ", "), "\n")

# ============================================================
# FUNÇÃO: ler e filtrar um VCF
# ============================================================

# Filtro dos genes da lista do OncoKB
# oncokb <- read.delim(
  # "//wsl.localhost/Ubuntu/home/shuzoarthur/mestrado/R_dic/cancerGeneList.tsv",
  # sep = "\t",
  # header = TRUE,
  # stringsAsFactors = FALSE
# )

# colnames(oncokb)
# genes_oncokb <- unique(oncokb$Hugo.Symbol)


ler_vcf <- function(arquivo) {
  read.table(arquivo, sep = "\t", header = FALSE,
             col.names = col_names, stringsAsFactors = FALSE) %>%
    # filter(
      # region %in% c("exonic", "splicing"),          # só exônicas/splicing
      # !func %in% c("synonymous_SNV","synonymous SNV",
                   # "unknown"),                   # remove sinônimas
      # DP     >= 50,                                  # profundidade mínima
      # VAF    > 0.05,                                     # VAF positivo
      # GENE %in% genes_oncokb
    # ) %>%
    mutate(variant_id = paste(gene_id, gene_strand, gene_pos, TE_id, TE_strand, 
               TE_pos, X8, sep = "_"))
}

# ============================================================
# FUNÇÃO: classificar variantes e montar frac.table
# ============================================================

classificar_variantes <- function(vcf_pre, vcf_pos) {
  
  ids_pre <- vcf_pre$variant_id
  ids_pos <- vcf_pos$variant_id
  
  # 1. Compartilhadas
  compartilhadas <- inner_join(
    vcf_pre %>% select(variant_id, gene_id, chim_reads) %>% rename(chim_reads_pre = chim_reads),
    vcf_pos %>% select(variant_id, chim_reads)        %>% rename(chim_reads_pos = chim_reads),
    by = "variant_id"
  ) %>% mutate(categoria = "compartilhada")
  
  # 2. Exclusivas do pré (eliminadas pela quimio)
  exclusivas_pre <- vcf_pre %>%
    filter(!variant_id %in% ids_pos) %>%
    select(variant_id, gene_id, chim_reads) %>%
    rename(chim_reads_pre = chim_reads) %>%
    mutate(chim_reads_pos  = 0,
           categoria = "exclusiva_pre")
  
  # 3. Exclusivas do pós (emergentes/selecionadas)
  exclusivas_pos <- vcf_pos %>%
    filter(!variant_id %in% ids_pre) %>%
    select(variant_id, gene_id, chim_reads) %>%
    rename(chim_reads_pos = chim_reads) %>%
    mutate(chim_reads_pre  = 0,
           categoria = "exclusiva_pos")
  
  bind_rows(compartilhadas, exclusivas_pre, exclusivas_pos)
}

# ============================================================
# FUNÇÃO: criar e plotar fishplot para um paciente (CHIMERATE)
# ============================================================
plotar_fishplot_reads <- function(paciente, variantes) {
  
  # 1. Preparar os dados
  df <- variantes %>%
    mutate(
      chim_reads_pre = as.numeric(chim_reads_pre),
      chim_reads_pos = as.numeric(chim_reads_pos),
      chim_reads_pre = replace_na(chim_reads_pre, 0),
      chim_reads_pos = replace_na(chim_reads_pos, 0),
      label = paste0(gene_id, "\n(", categoria, ")")
    )
  
  n_clones <- nrow(df)
  
  if (n_clones == 0) {
    cat("Sem variantes para paciente", paciente, "\n")
    return(NULL)
  }
  
  # 2. Calcular a porcentagem relativa de reads por evento
  total_pre <- sum(df$chim_reads_pre)
  total_pos <- sum(df$chim_reads_pos)
  
  df <- df %>%
    mutate(
      frac_pre = if (total_pre > 0) {
        chim_reads_pre / total_pre * 100
      } else {
        0
      },
      frac_pos = if (total_pos > 0) {
        chim_reads_pos / total_pos * 100
      } else {
        0
      }
    ) %>%
    arrange(desc(frac_pre))
  
  # 3. Montar a matriz: linhas = eventos, colunas = amostras
  frac_table <- matrix(
    c(df$frac_pre, df$frac_pos),
    nrow = n_clones,
    ncol = 2,
    byrow = FALSE
  )
  
  # Substituir zeros por um mínimo técnico para visualização
  frac_table[frac_table == 0] <- 0.001
  
  # Garantir que a soma por amostra não ultrapasse 100%
  for (tp in 1:2) {
    soma <- sum(frac_table[, tp])
    
    if (soma > 100) {
      frac_table[, tp] <- frac_table[, tp] / soma * 99
    }
  }
  
  # 4. Todos os eventos são tratados como independentes
  parents <- rep(0, n_clones)
  
  # 5. Cores por categoria
  cores <- case_when(
    df$categoria == "compartilhada" ~ "#E74C3C",
    df$categoria == "exclusiva_pre" ~ "#3498DB",
    df$categoria == "exclusiva_pos" ~ "#2ECC71",
    TRUE ~ "#808080"
  )
  
  # 6. Criar o Fishplot
  tryCatch({
    
    fish <- createFishObject(
      frac.table = frac_table,
      parents = parents,
      timepoints = c(0, 1),
      fix.missing.clones = TRUE,
      clone.labels = df$label
    )
    
    fish <- layoutClones(fish)
    fish <- setCol(fish, col = cores)
    
    fishPlot(
      fish,
      shape = "spline",
      title.btm = paste("Paciente", paciente),
      cex.title = 0.5,
      vlines = c(0, 1),
      vlab = c("Pré-tratamento", "Pós-tratamento"),
      pad.left = 0.2
    )
    
    # 7. Resumo
    cat("\n=== Paciente", paciente, "===\n")
    cat("Total de reads no pré:", total_pre, "\n")
    cat("Total de reads no pós:", total_pos, "\n")
    cat("Compartilhadas:",
        sum(df$categoria == "compartilhada"), "\n")
    cat("Exclusivas pré:",
        sum(df$categoria == "exclusiva_pre"), "\n")
    cat("Exclusivas pós:",
        sum(df$categoria == "exclusiva_pos"), "\n")
    
  }, error = function(e) {
    cat("Erro no paciente", paciente, ":",
        conditionMessage(e), "\n")
  })
}
    
#============================================================
 #  LOOP PRINCIPAL: processar todos os pacientes
#============================================================
    
for (pac in pacientes) {
  
  arquivo_pre <- file.path(
    pasta,
    paste0(pac, "_TE-exonized_embedded_final_pre.tsv")
  )
  
  arquivo_pos <- file.path(
    pasta,
    paste0(pac, "_TE-exonized_embedded_final_pos.tsv")
  )
  
#   Verificar se ambos os arquivos existem
  
  if (!file.exists(arquivo_pre) || !file.exists(arquivo_pos)) {
    cat("Arquivos ausentes para paciente", pac, "\n")
    next
  }
  
#  Importar os dados
  
  vcf_pre <- ler_vcf(arquivo_pre)
  vcf_pos <- ler_vcf(arquivo_pos)
  
 # Classificar as variantes/eventos
  
  variantes <- classificar_variantes(vcf_pre, vcf_pos)
  
  #Gerar o Fishplot usando reads quiméricas
  
  plotar_fishplot_reads(pac, variantes)
}
  