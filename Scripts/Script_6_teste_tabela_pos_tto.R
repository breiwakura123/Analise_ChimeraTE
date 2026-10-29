#tudo isso ta no chat : Script_tabela

# ============================================================
# RESUMO DA ANÁLISE ChimeraTE
# ============================================================
#install.packages("gdtools")
#install.packages("flextable")

library(dplyr)
library(readr)
library(here)
library(flextable)
library(officer)
library(tidyverse)
library(writexl)
library(gt)

#### Configuração da análise ####

tipo_amostra <- "pre_tto_selecionados"

# Opções:
# "pre_tto"
# "pos_tto"
# "pre_tto_selecionados"

#### Definindo as pastas ####

if (tipo_amostra == "pre_tto") {
  pasta_processados <- here("Dados", "Processados", "pre_tto")
} else if (tipo_amostra == "pos_tto") {
  pasta_processados <- here("Dados", "Processados", "pos_tto")
} else if (tipo_amostra == "pre_tto_selecionados") {
  pasta_processados <- here(
    "Dados", "Processados", "pre_tto_selecionados"
  )
} else {
  stop("Tipo de amostra inválido!")
}

pasta_resultados <- here("Resultados", tipo_amostra)

# Pasta para as tabelas deste grupo
pasta_tabelas <- file.path(pasta_resultados, "Tabelas_Resumo")

dir.create(
  pasta_tabelas,
  recursive = TRUE,
  showWarnings = FALSE
)

#### Carregando os dados ####

chim_class <- readRDS(
  file.path(pasta_processados, "chim_class.rds")
)

# Conferir se o arquivo possui as colunas necessárias
colnames(chim_class)

if (!all(c(
  "id", "gene_id", "TE_id", "family_id",
  "class_id", "Category", "chim_reads"
) %in% colnames(chim_class))) {
  stop("O arquivo chim_class.rds não contém todas as colunas esperadas.")
}


# ============================================================
# TABELA 1 — RESUMO GERAL
# ============================================================

resumo_geral <- chim_class %>%
  summarise(
    amostras = n_distinct(id),
    genes = n_distinct(gene_id),
    TEs = n_distinct(TE_id),
    familias_TE = n_distinct(family_id),
    classes_TE = n_distinct(class_id),
    quimeras = n(),
    categorias = n_distinct(Category),
    reads_quimericas = sum(chim_reads, na.rm = TRUE),
    reads_media_por_quimera = mean(chim_reads, na.rm = TRUE),
    reads_mediana_por_quimera = median(chim_reads, na.rm = TRUE)
  )

resumo_geral

write_csv(
  resumo_geral,
  file.path(pasta_tabelas, "01_resumo_geral.csv")
)

# ============================================================
# TABELA 2 — RESUMO POR AMOSTRA/PACIENTE
# ============================================================

resumo_por_id <- chim_class %>%
  group_by(id) %>%
  summarise(
    genes = n_distinct(gene_id),
    TEs = n_distinct(TE_id),
    familias_TE = n_distinct(family_id),
    classes_TE = n_distinct(class_id),
    quimeras = n(),
    categorias = n_distinct(Category),
    reads_quimericas = sum(chim_reads, na.rm = TRUE),
    reads_media = mean(chim_reads, na.rm = TRUE),
    reads_mediana = median(chim_reads, na.rm = TRUE),
    .groups = "drop") %>%
  arrange(desc(quimeras))

resumo_por_id

write_csv(
  resumo_por_id,
  file.path(pasta_tabelas, "02_resumo_por_id.csv")
)


# ============================================================
# TABELA 3 — CATEGORIAS DE QUIMERAS
# ============================================================

resumo_categorias <- chim_class %>%
  count(Category, name = "quimeras") %>%
  mutate( percentual = quimeras / sum(quimeras) * 100 ) %>%
  arrange(desc(quimeras))

resumo_categorias

write_csv(
  resumo_categorias,
  file.path(pasta_tabelas, "03_resumo_categorias.csv")
)



# ============================================================
# TABELA 4 — CATEGORIAS POR AMOSTRA
# ============================================================

resumo_categorias_id <- chim_class %>%
  count(id, Category, name = "quimeras") %>%
  group_by(id) %>%
  mutate( percentual = quimeras / sum(quimeras) * 100 ) %>%
  ungroup() %>%
  arrange(id, desc(quimeras))

resumo_categorias_id

write_csv(
  resumo_categorias_id,
  file.path(pasta_tabelas, "04_categorias_por_id.csv")
)

# ============================================================
# TABELA 5 — CLASSES DE TE
# ============================================================

resumo_classes <- chim_class %>%
  
  count(class_id, name = "quimeras") %>%
  mutate(percentual = quimeras / sum(quimeras) * 100) %>%
  arrange(desc(quimeras))

resumo_classes

write_csv(
  resumo_classes,
  file.path(pasta_tabelas, "05_resumo_classes_TE.csv")
)


#Tirar as classes com ? e as NAs

# ============================================================
# TABELA 6 — CLASSES DE TE POR AMOSTRA
# ============================================================

resumo_classes_id <- chim_class %>%
  count(id, class_id, name = "quimeras") %>%
  group_by(id) %>%
  mutate( percentual = quimeras / sum(quimeras) * 100 ) %>%
  ungroup() %>%
  arrange(id, desc(quimeras))

resumo_classes_id

write_csv(
  resumo_classes_id,
  file.path(pasta_tabelas, "06_classes_TE_por_id.csv")
)

# ============================================================
# TABELA 7 — FAMÍLIAS DE TE MAIS FREQUENTES
# ============================================================

resumo_familias <- chim_class %>%
  group_by(family_id, class_id) %>%
  summarise(
    quimeras = n(),
    genes = n_distinct(gene_id),
    TEs = n_distinct(TE_id),
    amostras = n_distinct(id),
    reads_quimericas = sum(chim_reads, na.rm = TRUE),
    .groups = "drop") %>%
  mutate( percentual = quimeras / sum(quimeras) * 100) %>%
  arrange(desc(quimeras))

resumo_familias

write_csv(
  resumo_familias,
  file.path(pasta_tabelas, "07_resumo_familias_TE.csv")
)

# ============================================================
# TABELA 8 — FAMÍLIAS DE TE POR AMOSTRA
# ============================================================

resumo_familias_id <- chim_class %>%
  group_by(id, family_id, class_id) %>%
  summarise(
    quimeras = n(),
    genes = n_distinct(gene_id),
    TEs = n_distinct(TE_id),
    reads_quimericas = sum(chim_reads, na.rm = TRUE),
    .groups = "drop") %>%
  arrange(id, desc(quimeras))

resumo_familias_id

write_csv(
  resumo_familias_id,
  file.path(pasta_tabelas, "08_familias_TE_por_id.csv")
)

# ============================================================
# TABELA 9 — GENES COM MAIOR NÚMERO DE QUIMERAS
# ============================================================

resumo_genes <- chim_class %>%
  group_by(gene_id) %>%
  summarise(
    quimeras = n(),
    familias_TE = n_distinct(family_id),
    TEs = n_distinct(TE_id),
    amostras = n_distinct(id),
    classes_TE = n_distinct(class_id),
    reads_quimericas = sum(chim_reads, na.rm = TRUE),
    .groups = "drop") %>%
  arrange(desc(quimeras))

resumo_genes

write_csv(
  resumo_genes,
  file.path(pasta_tabelas, "09_resumo_genes.csv")
)

# ============================================================
# TABELA 10 — GENES POR AMOSTRA
# ============================================================

resumo_genes_id <- chim_class %>%
  group_by(id, gene_id) %>%
  summarise(
    quimeras = n(),
    familias_TE = n_distinct(family_id),
    TEs = n_distinct(TE_id),
    reads_quimericas = sum(chim_reads, na.rm = TRUE),
    .groups = "drop") %>%
  arrange(id, desc(quimeras))

resumo_genes_id

write_csv(
  resumo_genes_id,
  file.path(pasta_tabelas, "10_genes_por_id.csv")
)

# ============================================================
# TABELA 11 — GENES ASSOCIADOS A MÚLTIPLAS FAMÍLIAS DE TE
# ============================================================

genes_familias <- chim_class %>%
  group_by(gene_id) %>%
  summarise(
    familias_TE = n_distinct(family_id),
    classes_TE = n_distinct(class_id),
    TEs = n_distinct(TE_id),
    quimeras = n(),
    amostras = n_distinct(id),
    .groups = "drop") %>%
  arrange(desc(familias_TE), desc(quimeras))

genes_familias

write_csv(
  genes_familias,
  file.path(pasta_tabelas, "11_genes_multiplas_familias.csv")
)

# ============================================================
# TABELA 12 — COMBINAÇÃO GENE × FAMÍLIA DE TE
# ============================================================

gene_familia <- chim_class %>%
  group_by(gene_id, family_id, class_id) %>%
  summarise(
    quimeras = n(),
    TEs = n_distinct(TE_id),
    amostras = n_distinct(id),
    reads_quimericas = sum(chim_reads, na.rm = TRUE),
    .groups = "drop") %>%
  arrange(desc(quimeras))

gene_familia

write_csv(
  gene_familia,
  file.path(pasta_tabelas, "12_gene_familia_TE.csv")
)

# ============================================================
# TABELA 13 — GENE × CATEGORIA
# ============================================================

gene_categoria <- chim_class %>%
  group_by(gene_id, Category) %>%
  summarise(
    quimeras = n(),
    familias_TE = n_distinct(family_id),
    amostras = n_distinct(id),
    .groups = "drop") %>%
  arrange(desc(quimeras))

gene_categoria

write_csv(
  gene_categoria,
  file.path(pasta_tabelas, "13_gene_categoria.csv")
)

# ============================================================
# TABELA 14 — FAMÍLIA × CATEGORIA
# ============================================================

familia_categoria <- chim_class %>%
  group_by(family_id, class_id, Category) %>%
  summarise(
    quimeras = n(),
    genes = n_distinct(gene_id),
    amostras = n_distinct(id),
    .groups = "drop") %>%
  arrange(desc(quimeras))

familia_categoria

write_csv(
  familia_categoria,
  file.path(pasta_tabelas, "14_familia_categoria.csv")
)

# ============================================================
# TABELA 15 — DIVERSIDADE DE TE POR AMOSTRA
# ============================================================

diversidade_id <- chim_class %>%
  group_by(id) %>%
  summarise(
    genes = n_distinct(gene_id),
    TEs = n_distinct(TE_id),
    familias = n_distinct(family_id),
    classes = n_distinct(class_id),
    categorias = n_distinct(Category),
    quimeras = n(),
    .groups = "drop") %>%
  arrange(desc(familias), desc(quimeras))

diversidade_id

write_csv(
  diversidade_id,
  file.path(pasta_tabelas, "15_diversidade_TE_por_id.csv")
)

# ============================================================
# FINAL
# ============================================================

cat("Tabelas de resumo salvas em:\n")
cat(here("Dados", "Processados", "Tabelas_Resumo"), "\n")

# ============================================================
# TESTE PARA CONCATENAR TUDO
# ============================================================

library(tidyverse)
library(writexl)
library(here)


# ============================================================
# CONSOLIDAR AS TABELAS EM UM ARQUIVO EXCEL
# ============================================================

pasta_origem <- pasta_tabelas

arquivos_csv <- list.files(
  path = pasta_origem,
  pattern = "\\.csv$",
  full.names = TRUE
)

nomes_das_abas <- tools::file_path_sans_ext(
  basename(arquivos_csv)
)

lista_tabelas <- purrr::map(
  arquivos_csv,
  readr::read_csv,
  show_col_types = FALSE
)

names(lista_tabelas) <- nomes_das_abas

caminho_saida <- file.path(
  pasta_tabelas,
  "Tabelas_Consolidadas.xlsx"
)

writexl::write_xlsx(lista_tabelas, caminho_saida)

cat("Arquivo Excel salvo em:\n", caminho_saida, "\n")

#############################################################

#TeSTE

tabela_geral <- chim_class %>%
  summarise(
    `Amostras` = n_distinct(id),
    `Eventos de transcrição quimérica` = n(),
    `Genes envolvidos` = n_distinct(gene_id),
    `TEs envolvidos` = n_distinct(TE_id),
    `Famílias de TE` = n_distinct(family_id),
    `Classes de TE` = n_distinct(class_id),
    `Categorias de eventos` = n_distinct(Category)
  )

tabela_resumo <- data.frame(
  Caracteristica = c(
    "Características da análise",
    "Amostras",
    "Eventos de transcrição quimérica",
    "Genes envolvidos",
    "TEs envolvidos",
    "Famílias de TE",
    "Classes de TE",
    "Categorias de eventos"
  ),
  `n` = c(
    NA,
    n_distinct(chim_class$id),
    nrow(chim_class),
    n_distinct(chim_class$gene_id),
    n_distinct(chim_class$TE_id),
    n_distinct(chim_class$family_id),
    n_distinct(chim_class$class_id),
    n_distinct(chim_class$Category)
  )
)

tabela_ft <- flextable(tabela_resumo)

tabela_ft <- tabela_ft %>%
  theme_booktabs() %>%
  bold(
    i = 1,
    j = 1,
    bold = TRUE
  ) %>%
  align(
    j = 2,
    align = "right"
  ) %>%
  autofit()

save_as_docx(
  tabela_ft,
  path = here(
    "Dados",
    "Processados",
    "Tabelas_Resumo",
    "Tabela_resumo_ChimeraTE.docx"
  )
)

#install.packages("gt")

library(gt)

tabela_gt <- tabela_resumo %>%
  gt() %>%
  tab_header(
    title = "Características da análise"
  ) %>%
  cols_label(
    Caracteristica = "",
    n = "n"
  ) %>%
  tab_options(
    table.font.size = 12,
    table.border.top.style = "solid",
    table.border.bottom.style = "solid",
    column_labels.border.bottom.style = "solid",
    row_group.border.top.style = "solid"
  )

# ARRUMAR A PASTA DE SAIDA. TA SOBRESCREVENDO

tabela_gt
gtsave(
  tabela_gt,
  file.path(pasta_resultados, "Tabela_resumo_ChimeraTE.png"))

gtsave(
  file.path(pasta_resultados, "Tabela_resumo_ChimeraTE.png"),
  plot = tabela_gt,
  width = 9,
  height = 13,
  dpi = 300
)

###############################################################
# ============================================================
# RODAR ATÉ AQUI!!!!
# ============================================================

###############################################################

# ============================================================
# DADOS GERAIS
# ============================================================

n_amostras <- n_distinct(chim_class$id)

n_quimeras <- nrow(chim_class)

n_genes <- n_distinct(chim_class$gene_id)

n_TEs <- n_distinct(chim_class$TE_id)

n_familias <- n_distinct(chim_class$family_id)

n_classes <- n_distinct(chim_class$class_id)

n_categorias <- n_distinct(chim_class$Category)


# ============================================================
# CLASSES DE TE
# ============================================================

classes <- chim_class %>%
  filter(!is.na(class_id)) %>%
  count(class_id, name = "n") %>%
  mutate(
    percentual = n / sum(n) * 100
  ) %>%
  arrange(desc(n))


# ============================================================
# CATEGORIAS DE QUIMERAS
# ============================================================

categorias <- chim_class %>%
  filter(!is.na(Category)) %>%
  count(Category, name = "n") %>%
  mutate(
    percentual = n / sum(n) * 100
  ) %>%
  arrange(desc(n))


# ============================================================
# FAMÍLIAS DE TE
# ============================================================

familias <- chim_class %>%
  filter(!is.na(family_id)) %>%
  count(family_id, class_id, name = "n") %>%
  mutate(
    percentual = n / sum(n) * 100
  ) %>%
  arrange(desc(n))

# ============================================================
# TABELA — CARACTERÍSTICAS GERAIS
# ============================================================

tabela <- data.frame(
  secao = c(
    # Geral
    "Características gerais",
    "Amostras",
    "Eventos de transcrição quimérica",
    "Genes envolvidos",
    "TEs envolvidos",
    "Famílias de TE",
    "Classes de TE",
    "Categorias de eventos",
    
    # Classes
    "Classe de TE",
    
    # Categorias
    "Categoria do evento",
    
    # Famílias
    "Famílias de TE"
  ),
  
  variavel = c(
    "", 
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    
    "",
    
    "",
    
    ""
  ),
  
  n = c(
    NA,
    n_amostras,
    n_quimeras,
    n_genes,
    n_TEs,
    n_familias,
    n_classes,
    n_categorias,
    
    NA,
    
    NA,
    
    NA
  ),
  
  percentual = c(
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    
    NA,
    
    NA,
    
    NA
  )
)

tabela_classes <- data.frame(
  secao = "Classe de TE",
  variavel = paste0("    ", classes$class_id),
  n = classes$n,
  percentual = classes$percentual
)

tabela_categorias <- data.frame(
  secao = "Categoria do evento",
  variavel = paste0("    ", categorias$Category),
  n = categorias$n,
  percentual = categorias$percentual
)

top_familias <- familias %>%
  slice_head(n = 10)

tabela_familias <- data.frame(
  secao = "Famílias de TE",
  variavel = paste0("    ", top_familias$family_id),
  n = top_familias$n,
  percentual = top_familias$percentual
)

tabela_final <- bind_rows(
  tabela,
  tabela_classes,
  tabela_categorias,
  tabela_familias
)


tabela_final <- bind_rows(
  tabela,
  tabela_classes,
  tabela_categorias,
  tabela_familias
) %>%
  mutate(
    `n (%)` = case_when(
      is.na(n) ~ "",
      is.na(percentual) ~ as.character(n),
      TRUE ~ paste0(
        n, " (", sprintf("%.1f", percentual), "%)"
      )
    )
  ) %>%
  select(
    secao,
    variavel,
    `n (%)`
  )

names(tabela_final) <- c(
  "Seção",
  "Característica",
  "n (%)"
)

ft <- flextable(tabela_final)

ft <- ft %>%
  theme_booktabs() %>%
  font(fontname = "Arial", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  bold(
    i = which(
      tabela_final$Característica == ""
    ),
    j = 1,
    bold = TRUE
  ) %>%
  align(
    j = 3,
    align = "right",
    part = "all"
  ) %>%
  align(
    j = 1:2,
    align = "left",
    part = "all"
  ) %>%
  autofit()

secoes <- c(
  "Características gerais",
  "Classe de TE",
  "Categoria do evento",
  "Famílias de TE"
)

for (secao_atual in secoes) {
  
  linhas <- which(
    tabela_final$Seção == secao_atual &
      tabela_final$Característica == ""
  )
  
  if (length(linhas) > 0) {
    ft <- ft %>%
      bold(
        i = linhas,
        j = 1:3,
        bold = TRUE
      )
  }
}

tabela_final$Característica[
  tabela_final$Seção == "Características gerais" &
    tabela_final$Característica == ""
] <- "Características gerais"

tabela_final$Característica[
  tabela_final$Seção == "Classe de TE" &
    tabela_final$Característica == ""
] <- "Classe de TE"

tabela_final$Característica[
  tabela_final$Seção == "Categoria do evento" &
    tabela_final$Característica == ""
] <- "Categoria do evento"

tabela_final$Característica[
  tabela_final$Seção == "Famílias de TE" &
    tabela_final$Característica == ""
] <- "Famílias de TE"

ft <- flextable(
  tabela_final %>%
    select(
      Característica,
      `n (%)`
    )
)

linhas_secao <- which(
  tabela_final$Característica %in% c(
    "Características gerais",
    "Classe de TE",
    "Categoria do evento",
    "Famílias de TE"
  )
)

ft <- ft %>%
  bold(
    i = linhas_secao,
    j = 1:2,
    bold = TRUE
  )

save_as_docx(
  ft,
  path = here(
    "Dados",
    "Processados",
    "Tabelas_Resumo",
    "Tabela_Resumo_ChimeraTE.docx"
  )
)
ft




