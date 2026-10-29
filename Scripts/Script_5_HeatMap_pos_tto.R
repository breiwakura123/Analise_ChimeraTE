## =========================================================
## Heatmap gene_id x family_id (ChimeraTE)
## Cor = quantidade de vezes que a combinacao gene/familia
##       aparece na tabela (contagem por linha) OU
##       quantidade de "id" distintos associados a ela
## =========================================================

## -----------------------------------------------------
## 1. Carregar os dados
## -----------------------------------------------------

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(patchwork)

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
  pasta_processados <- here("Dados", "Processados", "pre_tto_selecionados")
} else {
  stop("Tipo de amostra inválido!")
}

pasta_resultados <- here("Resultados", tipo_amostra)

dir.create(
  pasta_resultados,
  recursive = TRUE,
  showWarnings = FALSE
)

#### Carregando os dados ####

chim <- readRDS(
  file.path(pasta_processados, "chim_class.rds")
)

# Conferindo a estrutura
colnames(chim)
head(chim)

## -----------------------------------------------------
## 2. Filtrar apenas os genes de interesse
## -----------------------------------------------------
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

chim_genes <- chim %>%
  filter(gene_id %in% genes_interesse) %>%
  filter(!is.na(family_id))  # remove linhas sem familia de TE (NA)


# conferindo se todos os 52 genes foram encontrados na tabela
setdiff(genes_interesse, unique(chim_genes$gene_id))

## -----------------------------------------------------
## 3. Selecionar o TOP 10 de family_id mais frequentes
##    (baseado no numero de linhas em que aparecem,
##     considerando somente os genes de interesse)
## -----------------------------------------------------
top10_families <- chim_genes %>%
  count(family_id, sort = TRUE) %>%
  slice_head(n = 10) %>%
  pull(family_id)

chim_top10 <- chim_genes %>%
  filter(family_id %in% top10_families)

## =========================================================
## OPCAO A: contagem por LINHA
## (quantas vezes esse par gene_id/family_id aparece na tabela,
##  contando repeticoes de "id")
## =========================================================
heat_linhas <- chim_top10 %>%
  count(gene_id, family_id, name = "n_linhas")

## garantir que todas as combinacoes gene x familia existam
## (mesmo que com valor 0), para o heatmap ficar completo
heat_linhas_full <- heat_linhas %>%
  complete(gene_id = genes_interesse, family_id = top10_families, fill = list(n_linhas = 0)) %>%
  mutate(family_id = factor(family_id, levels = top10_families))  # ordena por frequencia (maior -> menor)

p_linhas <- ggplot(heat_linhas_full,
                   aes(x = family_id, y = gene_id, fill = n_linhas)) +
  geom_tile(color = NA, width = 0.9, height = 0.9) +
  scale_fill_gradient(low = "white", high = "darkred", name = "N de linhas") +
  labs(
    title = "Heatmap: contagem por linha",
    x = "Familia de TE (top 10)",
    y = "Gene"
  ) +
  theme_gray(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13)
  )

print(p_linhas)


ggsave(
  file.path(pasta_resultados, "heatmap_gene_family_por_linha.png"),
  plot = p_linhas,
  width = 9,
  height = 13,
  dpi = 300
)

## =========================================================
## OPCAO B: contagem por ID DISTINTO
## (quantos valores UNICOS de "id" existem para esse par
##  gene_id / family_id)
## =========================================================
heat_id_distinto <- chim_top10 %>%
  group_by(gene_id, family_id) %>%
  summarise(n_id_distintos = n_distinct(id), .groups = "drop")

heat_id_distinto_full <- heat_id_distinto %>%
  complete(gene_id = genes_interesse, family_id = top10_families, fill = list(n_id_distintos = 0)) %>%
  mutate(family_id = factor(family_id, levels = top10_families))  # ordena por frequencia (maior -> menor)

p_id <- ggplot(heat_id_distinto_full,
               aes(x = family_id, y = gene_id, fill = n_id_distintos)) +
  geom_tile(color = NA, width = 0.9, height = 0.9) +
  scale_fill_gradient(low = "white", high = "darkblue", name = "N de pacientes") +
  labs(
    title = "Heatmap: contagem por paciente",
    x = "Familia de TE (top 10)",
    y = "Gene"
  ) +
  theme_gray(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13)
  )

print(p_id)


ggsave(
  file.path(pasta_resultados, "heatmap_gene_family_id_distinto.png"),
  plot = p_id,
  width = 9,
  height = 13,
  dpi = 300
)


## -----------------------------------------------------
## (Opcional) Salvar as tabelas usadas nos heatmaps
## -----------------------------------------------------

write.csv(
  heat_linhas_full,
  file.path(pasta_processados, "tabela_heatmap_por_linha.csv"),
  row.names = FALSE
)

write.csv(
  heat_id_distinto_full,
  file.path(pasta_processados, "tabela_heatmap_id_distinto.csv"),
  row.names = FALSE
)


## -----------------------------------------------------
## Unindo gráficos
## -----------------------------------------------------

library(patchwork)


figura_2_final <- p_linhas + p_id +
  plot_layout(widths = c(2, 3)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(
        size = 20,
        face = "bold"
      )
    )
  )

print(figura_2_final)

ggsave(
  file.path(pasta_resultados, "Figura_2_final.png"),
  plot = figura_2_final,
  width = 16,
  height = 10,
  dpi = 300
)



