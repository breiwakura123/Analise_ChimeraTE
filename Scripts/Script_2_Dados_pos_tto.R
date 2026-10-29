
####Importando pacotes####
library(tidyverse)
library(readr)
library(here)

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


#### Definindo pastas ####

if (tipo_amostra == "pre_tto") {
  
  pasta_raw <- here("Dados", "Raw_pre_tto")
  pasta_processados <- here("Dados", "Processados", "pre_tto")
  
} else if (tipo_amostra == "pos_tto") {
  
  pasta_raw <- here("Dados", "Raw_pos_tto")
  pasta_processados <- here("Dados", "Processados", "pos_tto")
  
} else if (tipo_amostra == "pre_tto_selecionados") {
  
  pasta_raw <- here("Dados", "Raw_pre_tto")
  pasta_processados <- here(
    "Dados", "Processados", "pre_tto_selecionados"
  )
  
} else {
  
  stop("Tipo de amostra inválido!")
  
}


#### Definindo output ####

pasta_resultados <- here("Resultados", tipo_amostra)

dir.create(
  pasta_processados,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  pasta_resultados,
  recursive = TRUE,
  showWarnings = FALSE
)


#### Carregar arquivos ####

arquivos <- list.files(
  pasta_raw,
  pattern = "final\\.tsv$",
  full.names = TRUE
)


#### Filtrar IDs, se necessário ####

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


#### Conferir arquivos selecionados ####

if (length(arquivos) == 0) {
  stop("Nenhum arquivo encontrado! Confira a pasta e os IDs.")
}

basename(arquivos)

####Importar os arquivos e classifica-los####

chim_all <- purrr::map_dfr(
  arquivos,  function(arq){  df <- read.delim(arq)
  nome <- basename(arq)
  df$id <- substr(nome, 1, 4)
    
    df$Category <- case_when(
      str_detect(nome,"embedded") ~ "Exonized",
      str_detect(nome,"overlapped") ~ "Exonized",
      str_detect(nome,"exonized_intronic") ~ "Exonized",
      str_detect(nome,"initiated") ~ "Initiated",
      str_detect(nome,"terminated") ~ "Terminated"  )
        df
  }
  )

saveRDS(
  chim_all,
  file.path(pasta_processados, "chim_all.rds")
)

chim_all


####Teste para separar os exonificados ####
#com a classificação toda

chim_all_full <- purrr::map_dfr(
  arquivos,   function(arq){  df <- read.delim(arq)
  nome <- basename(arq)
  df$id <- substr(nome,  1,  4  )
  
  df$Category <- case_when(
    str_detect(nome,"embedded") ~ "Embedded",
    str_detect(nome,"overlapped") ~ "Overlapped",
    str_detect(nome,"exonized_intronic") ~ "Intronic",
    str_detect(nome,"initiated") ~ "Initiated",
    str_detect(nome,"terminated") ~ "Terminated"  )
  df
  }
  )

saveRDS(
  chim_all_full,
  file.path(pasta_processados, "chim_all_full.rds")
)

chim_all_full


#==============================================================================


####Importando os dados processados####

te_info <- readRDS( here( "Dados","Processados","te_info.rds" ))


####Juntando####
chim_class <- left_join(chim_all,  te_info,  by = "TE_id") %>%
  filter(!is.na(Category)) 


chim_class

saveRDS(
  chim_class,
  file.path(pasta_processados, "chim_class.rds")
)

head(chim_class)

##### Revisando Perdas #####
sum(is.na(chim_class$class_id))

####Juntando_ALL####
chim_class_full <- left_join(chim_all_full,  te_info,  by = "TE_id") %>%
  filter(!is.na(Category)) 


chim_class_full

saveRDS(
  chim_all,
  file.path(pasta_processados, "chim_class_full..rds")
)

##### Revisando Perdas #####
sum(is.na(chim_class_full$class_id))

