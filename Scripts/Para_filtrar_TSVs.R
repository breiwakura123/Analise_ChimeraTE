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

input_dir  <- "Dados/Sem_filtro_pre_tto"   # pasta com os brutos TROCAR

output_dir <- "Dados/Processados/pre_tto_filtrado"             # pasta de saída

# padrão para encontrar os arquivos a filtrar (ajuste conforme nomenclatura real)
file_pattern <- ".*\\.tsv$"

MIN_READS  <- 2
ROUND_FUNC <- round   # round, floor ou trunc, conforme confirmado para o paciente 29

RAW_COLUMNS <- c("gene_id", "gene_strand", "gene_pos",
                 "TE_id", "TE_strand", "TE_pos", "chim_reads")

# ------------------------------------------------------------------------------

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

files <- list.files(input_dir, pattern = file_pattern, full.names = TRUE)
#faz a lista com todos os arquivos

if (length(files) == 0) {
  stop("Nenhum arquivo encontrado em '", input_dir, "' com o padrão '", file_pattern, "'.")
}
#se não achar nada ele vai avisar!

cat(sprintf("Encontrados %d arquivos para filtrar:\n", length(files)))
cat(paste0(" - ", basename(files), collapse = "\n"), "\n\n")
#lista os arquivos que foram para a variante input -> CONFERIR SE ESTÃO TODOS


resumo <- list()

for (f in files) {
  nome_base <- tools::file_path_sans_ext(basename(f))
  #pega o nome do arquivo
  
  df <- read_tsv(f, col_names = RAW_COLUMNS, show_col_types = FALSE)
  #atribui o nome das colunas
  df$chim_reads <- as.numeric(df$chim_reads)
  #transforma a columa chim_reads em numérica se ela já não estiver
  
  n_antes <- nrow(df)
  
  filtrado <- df %>%
    #mutate(chim_reads = ROUND_FUNC(chim_reads)) %>%
    #arredonda os valores
    filter(chim_reads >= MIN_READS)
  #pega só os maiores que 2
  
  n_depois <- nrow(filtrado)
  
  partes    <- str_match(nome_base, "^(.+)-(\\d+)$")
  novo_nome <- paste0(partes[,3], "_", partes[,2], "_final.tsv")
  out_path  <- file.path(output_dir, novo_nome)
  # para reorganizar o nome do arquivo para ficar igual o final do ChimeraTE
  
  write_tsv(filtrado, out_path)
  
  cat(sprintf("%s: %d -> %d linhas (salvo em %s)\n",
              basename(f), n_antes, n_depois, out_path))
  
  resumo[[nome_base]] <- c(n_antes = n_antes, n_depois = n_depois)
}

{cat("\n=== Resumo final ===\n")

print(do.call(rbind, resumo))}

