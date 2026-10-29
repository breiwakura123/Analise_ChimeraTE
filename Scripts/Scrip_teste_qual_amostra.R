# =============================================================================
# Script para:
# 1. Identificar, entre dois arquivos brutos do ChimeraTE (ex: paciente 28 e 29),
#    qual deles corresponde ao arquivo já filtrado (chim_reads >= 2, arredondado).
# 2. Aplicar o mesmo filtro no arquivo que ainda não foi filtrado, gerando um
#    .tsv no mesmo padrão (com header) usado para os demais pacientes.
#
# Uso (ajuste os caminhos na seção "CONFIGURAÇÃO" abaixo e rode o script inteiro)
# =============================================================================

library (readr)
library (dplyr)
library (here)

# ----------------------------- CONFIGURAÇÃO ---------------------------------

path_raw28    <- "Dados/temp_teste/TE-exonized_embedded-3828.tsv"
path_raw29    <- "Dados/temp_teste/TE-exonized_embedded-3829.tsv"
path_filtrado <- "Dados/Raw/3828_TE-exonized_embedded_final.tsv"
outdir        <- "Dados/Processados/teste"  

MIN_READS  <- 2          # limiar mínimo de reads
ROUND_FUNC <- round      # round, floor ou trunc, conforme o ChimeraTE usar

# Nomes de coluna esperados para os arquivos BRUTOS (sem header no arquivo)
RAW_COLUMNS <- c("gene_id", "gene_strand", "gene_pos",
                 "TE_id", "TE_strand", "TE_pos", "chim_reads")

# Colunas usadas como chave para comparar/casar linhas entre bruto e filtrado
KEY_COLUMNS <- c("gene_id", "gene_pos", "TE_id", "TE_pos")

# ------------------------------------------------------------------------------

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

load_raw <- function(path) {
  df <- read_tsv(path, col_names = RAW_COLUMNS, show_col_types = FALSE)
  df$chim_reads <- as.numeric(df$chim_reads)
  df
}

load_filtered <- function(path) {
  df <- read_tsv(path, col_names = TRUE, show_col_types = FALSE)
  names(df) <- trimws(names(df))
  df
}

apply_filter <- function(df) {
  df %>%
    mutate(chim_reads_rounded = ROUND_FUNC(chim_reads)) %>%
    filter(chim_reads_rounded >= MIN_READS)
}

make_key <- function(df) {
  apply(df[, KEY_COLUMNS], 1, function(row) paste(row, collapse = "|"))
}

match_score <- function(raw_df, filtered_df) {
  filtered_raw <- apply_filter(raw_df)
  
  keys_filtered_raw  <- unique(make_key(filtered_raw))
  keys_filtered_file <- unique(make_key(filtered_df))
  
  intersection <- intersect(keys_filtered_raw, keys_filtered_file)
  uniao        <- union(keys_filtered_raw, keys_filtered_file)
  
  jaccard <- if (length(uniao) > 0) length(intersection) / length(uniao) else 0
  
  list(
    n_filtered_raw  = nrow(filtered_raw),
    n_filtered_file = nrow(filtered_df),
    n_intersection  = length(intersection),
    jaccard         = jaccard
  )
}

# ------------------------------- EXECUÇÃO ------------------------------------

raw28    <- load_raw(path_raw28)

raw29    <- load_raw(path_raw29)

filtrado <- load_filtered(path_filtrado)

score28 <- match_score(raw28, filtrado)
score29 <- match_score(raw29, filtrado)

cat("=== Comparação ===\n")
cat("Arquivo 28:\n"); str(score28)
cat("Arquivo 29:\n"); str(score29)

if (score28$jaccard > score29$jaccard) {
  ja_filtrado_nome <- "28"
  falta_filtrar_df <- raw29
  falta_filtrar_nome <- "29"
} else if (score29$jaccard > score28$jaccard) {
  ja_filtrado_nome <- "29"
  falta_filtrar_df <- raw28
  falta_filtrar_nome <- "28"
} else {
  stop("Não foi possível identificar com confiança qual arquivo já está filtrado ",
       "(scores iguais ou ambos baixos). Verifique manualmente os valores de 'jaccard' acima.")
}

cat(sprintf("\n>> Arquivo %s corresponde ao já filtrado (jaccard mais alto).\n", ja_filtrado_nome))
cat(sprintf(">> Arquivo %s ainda precisa ser filtrado. Aplicando filtro...\n", falta_filtrar_nome))

resultado <- apply_filter(falta_filtrar_df) %>%
  select(-chim_reads_rounded)


out_path <- file.path(outdir, paste0("paciente_", falta_filtrar_nome, "_filtrado.tsv"))
write_tsv(resultado, out_path)

cat(sprintf("\nArquivo gerado: %s\n", out_path))
cat(sprintf("Linhas no bruto: %d | Linhas após filtro (>= %d reads): %d\n",
            nrow(falta_filtrar_df), MIN_READS, nrow(resultado)))





# =============================================================================
# Filtra em lote vários arquivos brutos do ChimeraTE (chim_reads >= MIN_READS)
#
# Espera que todos os arquivos de entrada:
#   - sejam .tsv
#   - NÃO tenham header
#   - tenham as colunas na ordem: gene_id, gene_strand, gene_pos,
#     TE_id, TE_strand, TE_pos, chim_reads
#
# Gera um arquivo filtrado (com header) para cada arquivo de entrada,
# no diretório de saída.
# =============================================================================

library(readr)
library(dplyr)
library(stringr)

# ----------------------------- CONFIGURAÇÃO ---------------------------------

input_dir  <- "Dados/temp_teste"   # pasta com os brutos
output_dir <- "Dados/Processados/teste"             # pasta de saída

# padrão para encontrar os arquivos a filtrar (ajuste conforme nomenclatura real)
file_pattern <- "29.*\\.tsv$"

MIN_READS  <- 2
ROUND_FUNC <- round   # round, floor ou trunc, conforme confirmado para o paciente 29

RAW_COLUMNS <- c("gene_id", "gene_strand", "gene_pos",
                 "TE_id", "TE_strand", "TE_pos", "chim_reads")

# ------------------------------------------------------------------------------

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

files <- list.files(input_dir, pattern = file_pattern, full.names = TRUE)

if (length(files) == 0) {
  stop("Nenhum arquivo encontrado em '", input_dir, "' com o padrão '", file_pattern, "'.")
}

cat(sprintf("Encontrados %d arquivos para filtrar:\n", length(files)))
cat(paste0(" - ", basename(files), collapse = "\n"), "\n\n")

resumo <- list()

for (f in files) {
  nome_base <- tools::file_path_sans_ext(basename(f))
  
  df <- read_tsv(f, col_names = RAW_COLUMNS, show_col_types = FALSE)
  df$chim_reads <- as.numeric(df$chim_reads)
  
  n_antes <- nrow(df)
  
  filtrado <- df %>%
    mutate(chim_reads = ROUND_FUNC(chim_reads)) %>%
    filter(chim_reads >= MIN_READS)
  
  n_depois <- nrow(filtrado)
  
  out_path <- file.path(output_dir, paste0(nome_base, "_filtrado.tsv"))
  write_tsv(filtrado, out_path)
  
  cat(sprintf("%s: %d -> %d linhas (salvo em %s)\n",
              basename(f), n_antes, n_depois, out_path))
  
  resumo[[nome_base]] <- c(n_antes = n_antes, n_depois = n_depois)
}

cat("\n=== Resumo final ===\n")
print(do.call(rbind, resumo))
