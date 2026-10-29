####Organizando e chamando pacotes####

require("tidyverse")
require("readr")
require("data.table")
require("stringr")
require("here")

conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")

####Carregando os TEs referência####

gtf <- read.delim(
  here("Dados", "Referências", "hg38_rmsk_TE.gtf"),
  header = FALSE,
  comment.char = "#"
)

colnames(gtf) <- c(
  "chr","source","feature","start","end",
  "score","strand","frame","attributes"
)

head(gtf)

#### Extrair informações####

te_info <- gtf %>% 
  mutate(TE_id = str_extract(attributes, "(?<=gene_id )[^;]+"), 
         family_id = str_extract(attributes, "(?<=family_id )[^;]+"),
         class_id = str_extract(attributes, "(?<=class_id )[^;]+")) %>%
  select(TE_id, family_id, class_id) %>%
  distinct()


#### Salvar ####
saveRDS(te_info,here( "Dados","Processados","te_info.rds" ))

