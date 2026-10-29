#Análise dos Genes de Interesse

####Pacotes ####
require("tidyverse")
require("readr")
require("data.table")
require("stringr")
require("here")

#### Organizando e Selecionando Dadps ####

#### Configuração da análise ####

tipo_amostra <- "pos_tto"
# Opções:
# "pre_tto"
# "pos_tto"
# "pre_tto_selecionados"


#### Definindo pastas ####

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

dir.create(
  pasta_resultados,
  recursive = TRUE,
  showWarnings = FALSE
)


#### Carregando os dados da análise escolhida ####

chim_class <- readRDS(
  file.path(pasta_processados, "chim_class.rds")
)

#### Seleção dos Genes ####
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
#######################################################################
#### DADOS ####
##### Gene x Familia - número de eventos ##### ####

familias_genes <- chim_class %>%
  filter(gene_id %in% genes_interesse) %>%  #vê se o gene das amostras está na nossa lista de ref de genes
  filter(!is.na(family_id)) %>% #tem que ter uma família
  group_by(gene_id,family_id) %>%
  summarise ( eventos = n (), .groups= "drop") %>% # resumo contando cada linha de 
                                                #cada grupo. O drop "apaga" o grupo criado por group_by.
  arrange(desc(eventos)) #organizar por gene_id?

head(familias_genes) # não está descrescente?

familias_genes

##### Gene x Familia x Paciente - número de eventos #####  

familias_genes_pacientes <- chim_class %>%
  #filter(gene_id %in% genes_interesse) %>% # PORQUE EU TIREI ISSO? CONFERIR
  filter (!is.na(family_id)) %>%
  group_by(gene_id, family_id, id) %>%
  summarise (eventos = n (), .groups = "drop") %>%
  arrange (desc(eventos))

head(familias_genes_pacientes)
  
##### Heatmap#####

heatmap_gene <- chim_class %>%
  filter (gene_id %in% genes_interesse) %>% #vê se o gene das amostras está na nossa lista de ref de genes
  filter(!is.na (family_id)) %>%
  group_by(gene_id,family_id) %>%
  summarise( eventos = n(), .groups = "drop")

#==============================================================================

#### GRÁFICOS ####
##### Famílias de TE por gene  #####

####### 1 - Com o facet_wrap para avaliar entre os genes######

familias_genes_plot2 <- ggplot(familias_genes, aes(x = family_id, y = eventos)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~gene_id, scales = "free_y") +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(
    strip.text = element_text(size = 13, face = "bold")  # ajuste o tamanho aqui
  ) +
  labs(x = "Família TE", y = "Número de eventos", title = "Famílias de TE por gene")

familias_genes_plot2

ggsave(
  file.path(pasta_resultados, "familias_genes_plot2.png"),
  familias_genes_plot2,
  width = 14,
  height = 12,
  dpi = 300
)

###### 2 Sem facet######

familias_agrupadas <- familias_genes %>%
  group_by(family_id) %>%
  summarise(eventos = sum(eventos), .groups = "drop")

familias_genes_plot3 <- ggplot(familias_agrupadas, aes(x = reorder(family_id, eventos), y = eventos)) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 12)  # ajuste o tamanho aqui
  ) +
  labs(x = "Família TE", y = "Número de eventos", title = "Famílias de TE")

familias_genes_plot3

ggsave(
  file.path(pasta_resultados, "familias_genes_plot3.png"),
  familias_genes_plot3,
  width = 16,
  height = 10,
  dpi = 300
)



  
  
  