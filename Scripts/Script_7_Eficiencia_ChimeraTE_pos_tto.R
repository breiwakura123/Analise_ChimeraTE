
##############################################################################
# ANÁLISE DE QUIMERAS - OUTPUT ChimeraTE 
# Consolidação das tabelas exonized/terminated/initiated + resumos + validação
##############################################################################


##############################################################################
# ANÁLISE DE QUIMERAS - OUTPUT ChimeraTE
##############################################################################

library(tidyverse)
library(here)
library(scales)

#### Definindo tipo de amostra ####

tipo_amostra <- "pre_tto_selecionados"

# Opções:
# "pre_tto"
# "pos_tto"
# "pre_tto_selecionados"

ids_selecionados <- c(
  "7609", "7613", "7614", "9025",
  "9039", "9484", "9488", "9504"
)

#### Definindo pastas de entrada ####

if (tipo_amostra == "pre_tto") {
  
  pasta_raw <- here("Dados", "Raw_pre_tto")
  
} else if (tipo_amostra == "pos_tto") {
  
  pasta_raw <- here("Dados", "Raw_pos_tto")
  
} else if (tipo_amostra == "pre_tto_selecionados") {
  
  pasta_raw <- here("Dados", "Raw_pre_tto")
  
} else {
  
  stop("Tipo de amostra inválido!")
  
}

#### Definindo pasta de saída ####

pasta_saida <- here("Resultados", tipo_amostra)

dir.create(
  pasta_saida,
  recursive = TRUE,
  showWarnings = FALSE
)

#### Função para listar e filtrar arquivos ####

listar_arquivos_chimerate <- function(pasta) {
  
  arquivos <- list.files(
    pasta,
    pattern = "_final\\.tsv$",
    full.names = TRUE
  )
  
  if (tipo_amostra == "pre_tto_selecionados") {
    
    arquivos <- arquivos[
      grepl(
        paste0(
          "^(",
          paste(ids_selecionados, collapse = "|"),
          ")_"
        ),
        basename(arquivos)
      )
    ]
    
  }
  
  if (length(arquivos) == 0) {
    stop(
      "Nenhum arquivo encontrado. Confira a pasta de entrada e os IDs."
    )
  }
  
  return(arquivos)
}

#Pasta Saida
pasta_dados <- pasta_raw
arquivos <- listar_arquivos_chimerate(pasta_dados)

# Lista com os 52 genes de interesse (IDs). Cole aqui a lista completa.
genes_interesse <- c( "APC", "ARID1B", "ATM",
                      "BARD1", "BRCA1", "BRCA2",
                      "CACNA1E", "CACNA1S", "CCDC39", "CDH1", "CEP250", "CHEK2",  "CSMD3",
                      "EGFR", "ERBB2", "ERBB4", "ESR1",
                      "GATA3", "GNAQ",
                      "HFE",
                      "MAP2K4", "MAP3K1", "MLH1", "MLL3", "MSH2", "MSH6", "MUTYH",
                      "NF1", "NIPBL", "NOTCH2",
                      "PALB2", "PIK3CA", "POLD1", "POLE", "PRKAR1A", "PMS2", "PTEN", "PTPN13",
                      "RAD50", "RAD51B", "RAD51C", "RAD51D", "RB1", "RET",
                      "SMC6", "SMURF2", "STAG1", "STK11",
                      "TDG", "TP53", "TYW1",
                      "WIZ")

##############################################################################
# 2. LISTAR ARQUIVOS E EXTRAIR METADADOS DO NOME
##############################################################################

arquivos <- list.files(pasta_dados, pattern = "_final\\.tsv$", full.names = TRUE)


# Regex captura: (1) ID da amostra | (2) tipo | (4) subtipo (opcional, só p/ exonized)
padrao <- "^(\\d+)_TE-(exonized|terminated|initiated)(?:_(embedded|overlapped|intronic))?_final\\.tsv$"

metadados <- tibble(
  caminho = arquivos,
  nome_arquivo = basename(arquivos)
) %>%
  mutate(
    match       = str_match(nome_arquivo, padrao),
    sample      = match[, 2],
    transcript_type = match[, 3],
    subtype     = match[, 4]   # NA quando não se aplica (terminated/initiated)
  ) %>%
  select(-match)

# Checagem: algum arquivo não bateu com o padrão? (ajustar regex se sim)
if (any(is.na(metadados$sample))) {
  warning("Alguns arquivos não bateram com o padrão esperado. Verifique:\n",
          paste(metadados$nome_arquivo[is.na(metadados$sample)], collapse = "\n"))
}

##############################################################################
# 3. LER TODOS OS ARQUIVOS E EMPILHAR
##############################################################################

ler_um_arquivo <- function(caminho, sample, transcript_type, subtype) {
  read_tsv(caminho, show_col_types = FALSE) %>%
    mutate(
      sample = sample,
      transcript_type = transcript_type,
      subtype = subtype
    )
}

df_all <- pmap_dfr(
  list(
    caminho = metadados$caminho,
    sample = metadados$sample,
    transcript_type = metadados$transcript_type,
    subtype = metadados$subtype
  ),
  ler_um_arquivo
)

# Reordenar colunas (metadados primeiro)
df_all <- df_all %>%
  relocate(sample, transcript_type, subtype)

cat("Total de linhas (eventos de quimera) importados:", nrow(df_all), "\n")
cat("Total de amostras detectadas:", n_distinct(df_all$sample), "\n")

# Salvar a tabela consolidada bruta (útil como base para tudo o mais)
write_csv(df_all, file.path(pasta_saida, "chimeras_todas_amostras_consolidado.csv"))

##############################################################################
# 4. TABELAS RESUMO - CONTAGEM DE QUIMERAS
##############################################################################

## 4.1 Total geral
resumo_total <- tibble(total_quimeras = nrow(df_all))

## 4.2 Por amostra
resumo_por_amostra <- df_all %>%
  count(sample, name = "n_quimeras") %>%
  arrange(desc(n_quimeras))

## 4.3 Por tipo de transcrito (exonized/terminated/initiated)
resumo_por_tipo <- df_all %>%
  count(transcript_type, name = "n_quimeras") %>%
  arrange(desc(n_quimeras)) %>%
  rename(
    "Tipo de transcrito" = transcript_type,
    "Número de quimeras" = n_quimeras
  )


## 4.3b Por subtipo (dentro de exonized: embedded/overlapped/intronic)
resumo_por_subtipo <- df_all %>%
  filter(!is.na(subtype)) %>%
  count(transcript_type, subtype, name = "n_quimeras")

## 4.4 Cruzado: amostra x tipo (tabela larga, boa para apresentação)
resumo_amostra_x_tipo <- df_all %>%
  count(sample, transcript_type) %>%
  pivot_wider(names_from = transcript_type, values_from = n, values_fill = 0) %>%
  mutate(total = rowSums(across(where(is.numeric)))) %>%
  arrange(desc(total))

# Exportar
write_csv(resumo_por_amostra, file.path(pasta_saida, "resumo_por_amostra.csv"))
write_csv(resumo_por_tipo, file.path(pasta_saida, "resumo_por_tipo.csv"))
write_csv(resumo_por_subtipo, file.path(pasta_saida, "resumo_por_subtipo.csv"))
write_csv(resumo_amostra_x_tipo, file.path(pasta_saida, "resumo_amostra_x_tipo.csv"))

#Exportar como png
#install.packages("gridExtra")  # só na primeira vez

library(gridExtra)
library(ggplot2)

tabela_tipo <- tableGrob(
  resumo_por_tipo,
  rows = NULL)
tabela_tipo


ggsave( here( "Resultados", "resumo_por_tipo.png"), tabela_tipo, width = 6,
        height = 2.5,
        units = "in",
        dpi = 300)


##############################################################################
# 5. QUIMERAS NOS 52 GENES DE INTERESSE
##############################################################################

df_genes_interesse <- df_all %>%
  filter(gene_id %in% genes_interesse)

## 5.1 Contagem total de quimeras por gene (inclui genes com 0, se quiser ver todos)
quimeras_por_gene <- tibble(gene_id = genes_interesse) %>%
  left_join(
    df_genes_interesse %>% count(gene_id, name = "n_quimeras"),
    by = "gene_id"
  ) %>%
  mutate(n_quimeras = replace_na(n_quimeras, 0)) %>%
  arrange(desc(n_quimeras))

## 5.2 Contagem por gene x tipo de transcrito
quimeras_por_gene_tipo <- df_genes_interesse %>%
  count(gene_id, transcript_type) %>%
  pivot_wider(names_from = transcript_type, values_from = n, values_fill = 0)

## 5.3 Em quantas das 31 amostras cada gene aparece (recorrência)
genes_recorrencia <- df_genes_interesse %>%
  distinct(gene_id, sample) %>%
  count(gene_id, name = "n_amostras_presente") %>%
  right_join(tibble(gene_id = genes_interesse), by = "gene_id") %>%
  mutate(n_amostras_presente = replace_na(n_amostras_presente, 0)) %>%
  arrange(desc(n_amostras_presente))

## 5.4 chim_reads médio/mediano por gene (suporte de evidência)
genes_chimreads <- df_genes_interesse %>%
  group_by(gene_id) %>%
  summarise(
    n_eventos = n(),
    media_chim_reads = mean(chim_reads, na.rm = TRUE),
    mediana_chim_reads = median(chim_reads, na.rm = TRUE),
    max_chim_reads = max(chim_reads, na.rm = TRUE)
  ) %>%
  arrange(desc(n_eventos))

write_csv(quimeras_por_gene, file.path(pasta_saida, "quimeras_por_gene_interesse.csv"))
write_csv(quimeras_por_gene_tipo, file.path(pasta_saida, "quimeras_por_gene_interesse_x_tipo.csv"))
write_csv(genes_recorrencia, file.path(pasta_saida, "genes_interesse_recorrencia_amostras.csv"))
write_csv(genes_chimreads, file.path(pasta_saida, "genes_interesse_chim_reads.csv"))

##############################################################################
# 6. ANÁLISES EXTRAS - "PROVA DE FUNCIONAMENTO" DO ChimeraTE
##############################################################################

## 6.1 Distribuição de chim_reads (suporte de leitura) por tipo
resumo_chimreads_tipo <- df_all %>%
  group_by(transcript_type) %>%
  summarise(
    n = n(),
    media = mean(chim_reads, na.rm = TRUE),
    mediana = median(chim_reads, na.rm = TRUE),
    min = min(chim_reads, na.rm = TRUE),
    max = max(chim_reads, na.rm = TRUE),
    sd = sd(chim_reads, na.rm = TRUE)
  )

# Boxplot chim_reads por tipo de transcrito
p_boxplot_reads <- ggplot(df_all, aes(x = transcript_type, y = chim_reads, fill = transcript_type)) +
  geom_boxplot(outlier.alpha = 0.3) +
  scale_y_log10() +
  labs(title = "Suporte de leitura (chim_reads) por tipo de transcrito",
       x = "Tipo de transcrito", y = "chim_reads (log10)") +
  theme_minimal() +
  theme(legend.position = "none")
p_boxplot_reads

ggsave(file.path(pasta_saida, "boxplot_chimreads_por_tipo.png"), p_boxplot_reads,
       width = 6, height = 5, dpi = 300)

## 6.2 Genes únicos e TEs únicos envolvidos (total e por tipo)
genes_unicos_total <- n_distinct(df_all$gene_id)
tes_unicos_total   <- n_distinct(df_all$TE_id)

genes_tes_unicos_por_tipo <- df_all %>%
  group_by(transcript_type) %>%
  summarise(
    genes_unicos = n_distinct(gene_id),
    TEs_unicos = n_distinct(TE_id)
  )

cat("Genes únicos envolvidos em quimeras (total):", genes_unicos_total, "\n")
cat("TEs únicos envolvidos em quimeras (total):", tes_unicos_total, "\n")

write_csv(genes_tes_unicos_por_tipo, file.path(pasta_saida, "genes_tes_unicos_por_tipo.csv"))

## 6.3 Concordância de strand (gene_strand vs TE_strand)
resumo_strand <- df_all %>%
  mutate(concordancia_strand = if_else(gene_strand == TE_strand, "mesma_strand", "strand_oposta")) %>%
  count(transcript_type, concordancia_strand) %>%
  pivot_wider(names_from = concordancia_strand, values_from = n, values_fill = 0)

write_csv(resumo_strand, file.path(pasta_saida, "resumo_concordancia_strand.csv"))

## 6.4 Recorrência entre amostras: pares gene-TE que aparecem em mais de 1 amostra
recorrencia_gene_te <- df_all %>%
  distinct(gene_id, TE_id, sample) %>%
  count(gene_id, TE_id, name = "n_amostras") %>%
  arrange(desc(n_amostras))

cat("Pares gene-TE encontrados em mais de uma amostra:",
    sum(recorrencia_gene_te$n_amostras > 1), "de", nrow(recorrencia_gene_te), "pares totais\n")

write_csv(recorrencia_gene_te, file.path(pasta_saida, "recorrencia_gene_TE_entre_amostras.csv"))

## 6.5 Gráfico de barras empilhadas: amostra x tipo de transcrito
p_barras_empilhadas <- df_all %>%
  count(sample, transcript_type) %>%
  ggplot(aes(x = sample, y = n, fill = transcript_type)) +
  geom_col() +
  labs(title = "Quimeras por amostra e tipo de transcrito",
       x = "Amostra", y = "Número de quimeras", fill = "Tipo") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))
p_barras_empilhadas

ggsave(file.path(pasta_saida, "barras_empilhadas_amostra_tipo.png"), p_barras_empilhadas,
       width = 10, height = 5, dpi = 300)

##############################################################################
# 7. RESUMO FINAL NO CONSOLE
##############################################################################

cat("\n===== RESUMO GERAL =====\n")
print(resumo_total)
cat("\n===== POR TIPO DE TRANSCRITO =====\n")
print(resumo_por_tipo)
cat("\n===== TOP 10 AMOSTRAS COM MAIS QUIMERAS =====\n")
print(head(resumo_por_amostra, 10))
cat("\n===== TOP 10 GENES DE INTERESSE COM MAIS QUIMERAS =====\n")
print(head(quimeras_por_gene, 10))

cat("\nTodas as tabelas e gráficos foram salvos em:", pasta_saida, "\n")





##############################################################################






########################################################################
# ANÁLISE DE QUIMERAS - OUTPUT ChimeraTE (31 amostras)
# Tabelas e gráficos exportáveis para publicação
########################################################################

# ---- 0. PACOTES ------------------------------------------------------
# install.packages(c("tidyverse","stringr","ComplexUpset","scales"))
library(tidyverse)   # dplyr, readr, tidyr, ggplot2, purrr, stringr
library(scales)
# ComplexUpset é opcional, só usado na seção 7 (Upset plot)
# install.packages("ComplexUpset")
library(ComplexUpset)

# ---- 1. CONFIGURAÇÕES DO USUÁRIO --------------------------------------

# Pasta onde estão todos os .tsv do ChimeraTE (todas as 31 amostras juntas)
input_dir  <-  here("Dados", "Raw")

# Pasta de saída para tabelas e figuras
output_dir <-  here("Resultados", "Geral")

# ---- 2. IMPORTAR E PARSEAR TODOS OS ARQUIVOS --------------------------

# Padrão dos nomes observado:
#   3828_TE-exonized_embedded_final.tsv
#   3828_TE-exonized_overlapped_final.tsv
#   3828_TE-exonized_intronic_final.tsv
#   3828_TE-initiated_final.tsv
#   3828_TE-terminated_final.tsv

files <- list.files(input_dir, pattern = "_final\\.tsv$", full.names = TRUE)

if (length(files) == 0) {
  stop("Nenhum arquivo .tsv encontrado em input_dir. Verifique o caminho.")
}

# Regex: captura sample_id (tudo antes de "_TE-") e o "tipo completo"
# (tudo entre "TE-" e "_final.tsv")
parse_filename <- function(fname) {
  base <- basename(fname)
  m <- str_match(base, "^(.+?)_TE-(.+)_final\\.tsv$")
  sample_id  <- m[, 2]
  full_type  <- m[, 3]   # ex: "exonized_embedded", "initiated", "terminated"
  
  transcript_type <- str_extract(full_type, "^(exonized|initiated|terminated)")
  subtype <- if_else(
    transcript_type == "exonized",
    str_remove(full_type, "^exonized_"),
    NA_character_
  )
  
  tibble(sample_id = sample_id,
         transcript_type = transcript_type,
         exonized_subtype = subtype)
}

col_types_spec <- cols(
  gene_id     = col_character(),
  gene_strand = col_character(),
  gene_pos    = col_character(),
  TE_id       = col_character(),
  TE_strand   = col_character(),
  TE_pos      = col_character(),
  chim_reads  = col_double()
)

read_one <- function(fname) {
  meta <- parse_filename(fname)
  df <- read_delim(fname, delim = "\t", col_types = col_types_spec)
  bind_cols(df, meta[rep(1, nrow(df)), ])   # replica metadados p/ todas as linhas
}

all_data <- map_dfr(files, read_one)

# checagem rápida
message("Total de linhas importadas: ", nrow(all_data))
message("Amostras detectadas: ", n_distinct(all_data$sample_id))
message("Tipos detectados: ", paste(unique(all_data$transcript_type), collapse = ", "))

write_csv(all_data, file.path(output_dir, "todas_quimeras_consolidado.csv"))

# ---- 3. TABELA 1: quimeras por amostra x tipo -------------------------

tab_amostra_tipo <- all_data %>%
  count(sample_id, transcript_type, name = "n_chimeras") %>%
  pivot_wider(names_from = transcript_type, values_from = n_chimeras, values_fill = 0) %>%
  mutate(total = rowSums(across(where(is.numeric)))) %>%
  arrange(sample_id)

# linha de total geral
tab_amostra_tipo_totais <- tab_amostra_tipo %>%
  summarise(across(where(is.numeric), sum)) %>%
  mutate(sample_id = "TOTAL")

tab_amostra_tipo_final <- bind_rows(tab_amostra_tipo, tab_amostra_tipo_totais)

write_csv(tab_amostra_tipo_final, file.path(output_dir, "tabela1_amostra_x_tipo.csv"))

# resumo simples só por tipo (todas amostras somadas)
tab_por_tipo <- all_data %>%
  count(transcript_type, name = "n_chimeras") %>%
  mutate(pct = percent(n_chimeras / sum(n_chimeras)))

write_csv(tab_por_tipo, file.path(output_dir, "tabela1b_total_por_tipo.csv"))

# ---- 4. TABELA 2: genes únicos e TEs únicos por tipo -------------------

tab_unicos <- all_data %>%
  group_by(transcript_type) %>%
  summarise(
    n_eventos     = n(),
    genes_unicos  = n_distinct(gene_id),
    TEs_unicos    = n_distinct(TE_id),
    .groups = "drop"
  )

tab_unicos_geral <- all_data %>%
  summarise(
    transcript_type = "TODOS",
    n_eventos    = n(),
    genes_unicos = n_distinct(gene_id),
    TEs_unicos   = n_distinct(TE_id)
  )

tab_unicos_final <- bind_rows(tab_unicos, tab_unicos_geral)
write_csv(tab_unicos_final, file.path(output_dir, "tabela2_genes_TEs_unicos.csv"))

# ---- 5. TABELA 3: genes de interesse (52 genes) ------------------------

dados_genes_alvo <- all_data %>%
  filter(gene_id %in% genes_interesse)

tab_genes_alvo <- dados_genes_alvo %>%
  group_by(gene_id) %>%
  summarise(
    n_chimeras       = n(),
    n_amostras       = n_distinct(sample_id),
    tipos_encontrados = paste(sort(unique(transcript_type)), collapse = "; "),
    media_chim_reads = round(mean(chim_reads), 2),
    soma_chim_reads  = sum(chim_reads),
    .groups = "drop"
  ) %>%
  arrange(desc(n_amostras), desc(n_chimeras))

# genes da lista que NÃO tiveram nenhuma quimera encontrada (importante reportar!)
genes_sem_chimera <- setdiff(genes_interesse, tab_genes_alvo$gene_id)
if (length(genes_sem_chimera) > 0) {
  tab_genes_alvo <- bind_rows(
    tab_genes_alvo,
    tibble(gene_id = genes_sem_chimera, n_chimeras = 0, n_amostras = 0,
           tipos_encontrados = NA, media_chim_reads = NA, soma_chim_reads = 0)
  )
}

write_csv(tab_genes_alvo, file.path(output_dir, "tabela3_genes_alvo.csv"))

# ---- 6. MÉTRICAS EXTRAS: concordância de strand ------------------------

tab_strand <- all_data %>%
  mutate(concordancia = if_else(gene_strand == TE_strand, "senso (mesma fita)", "antisenso (fita oposta)")) %>%
  count(transcript_type, concordancia, name = "n") %>%
  group_by(transcript_type) %>%
  mutate(pct = percent(n / sum(n))) %>%
  ungroup()

write_csv(tab_strand, file.path(output_dir, "tabela_extra_concordancia_strand.csv"))

# ---- 7. RECORRÊNCIA ENTRE AMOSTRAS (reprodutibilidade) -----------------

# define uma "quimera única" como par gene_id + TE_id
recorrencia <- all_data %>%
  distinct(gene_id, TE_id, sample_id) %>%
  count(gene_id, TE_id, name = "n_amostras_compartilhadas") %>%
  arrange(desc(n_amostras_compartilhadas))

write_csv(recorrencia, file.path(output_dir, "tabela_supl_recorrencia_amostras.csv"))

recorrentes_2ou_mais <- recorrencia %>% filter(n_amostras_compartilhadas >= 2)
message(nrow(recorrentes_2ou_mais), " pares gene-TE recorrentes em >= 2 amostras")

########################################################################
# GRÁFICOS (exportados em PDF e PNG, 300 dpi)
########################################################################

#### Paleta de cores ####
#install.packages("paletteer")
library(paletteer)

scale_colour_paletteer_d("tvthemes::Alexandrite")
scale_color_paletteer_d("tvthemes::Alexandrite")
scale_fill_paletteer_d("tvthemes::Alexandrite")
paletteer_d("tvthemes::Alexandrite")

########


theme_pub <- theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    plot.title = element_text(face = "bold")  )


salvar_plot <- function(p, nome, w = 8, h = 5) {
  ggsave(file.path(output_dir, paste0(nome, ".pdf")), p, width = w, height = h)
  ggsave(file.path(output_dir, paste0(nome, ".png")), p, width = w, height = h, dpi = 300)
}

# ---- Fig 1: barras empilhadas amostra x tipo ---------------------------
p1 <- all_data %>%
  count(sample_id, transcript_type) %>%
  ggplot(aes(x = sample_id, y = n, fill = transcript_type)) +
  geom_col() +
  scale_y_continuous(expand = c(0,0))+
  labs(x = "Amostra", y = "Nº de quimeras", fill = "Tipo de transcrito",
       title = "Quimeras detectadas por amostra") +
  theme_pub
p1

salvar_plot(p1, "fig1_quimeras_por_amostra", w = 12, h = 6)

# ---- Fig 2: total por tipo (com subtipos de exonized) -------------------
p2_data <- all_data %>%
  mutate(categoria = if_else(
    transcript_type == "exonized" & !is.na(exonized_subtype),
    paste0("exonized_", exonized_subtype),
    transcript_type
  )) %>%
  count(categoria)

p2 <- ggplot(p2_data, aes(x = fct_reorder(categoria, n), y = n)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  scale_y_continuous(expand = c(0,0))+
  labs(x = NULL, y = "Nº de quimeras",
       title = "Total de quimeras por tipo/subtipo de transcrito") +
  theme_pub + theme(axis.text.x = element_text(angle = 0))

p2

salvar_plot(p2, "fig2_total_por_tipo", w = 7, h = 5)

# ---- Fig 5: heatmap gene-alvo x amostra (presença/ausência) --------------
heatmap_data <- expand_grid(gene_id = genes_interesse,
                            sample_id = unique(all_data$sample_id)) %>%
  left_join(
    all_data %>% filter(gene_id %in% genes_interesse) %>%
      count(gene_id, sample_id, name = "n_chimeras"),
    by = c("gene_id", "sample_id")
  ) %>%
  mutate(n_chimeras = replace_na(n_chimeras, 0))

p5 <- ggplot(heatmap_data, aes(x = sample_id, y = fct_reorder(gene_id, n_chimeras, .fun = sum),
                               fill = n_chimeras)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient(low = "grey95", high = "firebrick", name = "Nº quimeras") +
  labs(x = "Amostra", y = "Gene", title = "Quimeras dos genes-alvo por amostra") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6))
p5

salvar_plot(p5, "fig5_heatmap_genes_alvo", w = 12, h = 10)

# ---- Fig 6: ranking dos genes-alvo por nº de amostras ---------------------
p6 <- tab_genes_alvo %>%
  filter(n_amostras > 0) %>%
  ggplot(aes(x = fct_reorder(gene_id, n_amostras), y = n_amostras)) +
  geom_col(fill = "darkorange") +
  coord_flip() +
  scale_y_continuous(expand = c(0,0))+
  labs(x = "Gene", y = "Nº de amostras (de 24) com quimera detectada",
       title = "Recorrência dos genes-alvo entre amostras") +
  theme_pub + theme(axis.text.x = element_text(angle = 0))

p6

salvar_plot(p6, "fig6_ranking_genes_alvo", w = 7, h = 10)


########################################################################
message("Concluído! Tabelas e figuras salvas em: ", output_dir)
########################################################################
