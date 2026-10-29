
#### Importando pacotes ####

library(tidyverse)
library(readr)
library(here)
library(paletteer)
library(cowplot)

#### Definindo tipo de amostra ####

## TROCAR O TIPO DE AMOSTRA !!!! ##
tipo_amostra <- "pre_tto_selecionados"
# Opções:
# "pre_tto"
# "pos_tto"
# "pre_tto_selecionados"

#### Definindo pastas ####

pasta_processados <- here(
  "Dados",
  "Processados",
  tipo_amostra
)


pasta_resultados <- here(
  "Resultados",
  tipo_amostra
)


####Tabela ####

#Eventos por categoria e classe:
class_counts <- chim_class %>%
  count(  Category,    class_id) %>%
  filter(class_id != "DNA?") %>%
  filter(class_id != "LTR?") %>%
  filter(class_id != "RC?") %>%
  filter(!is.na(Category)) 

class_counts

write.csv(
  class_counts,
  file.path(pasta_resultados, "TE_class_counts.csv"),
  row.names = FALSE
)

#Eventos por categoria e classe FULL:
class_counts_full <- chim_class_full %>%
  count(Category,    class_id) %>%
  filter(class_id != "DNA?") %>%
  filter(class_id != "LTR?") %>%
  filter(class_id != "RC?") %>%
  filter(!is.na(Category)) 

class_counts_full

write.csv(
  class_counts_full,
  file.path(pasta_resultados, "TE_class_counts_full.csv"),
  row.names = FALSE
)

#### Teste primeiro gráfico ####

#Paleta

library(paletteer)

set_theme()
scale_colour_paletteer_d("tvthemes::Alexandrite")
scale_color_paletteer_d("tvthemes::Alexandrite")
scale_fill_paletteer_d("tvthemes::Alexandrite")
paletteer_d("tvthemes::Alexandrite")

# gráfico 1 Número de quimeras por pacientes por classificação

resultado_chimeraTE_table <- chim_class %>%
  count(id, Category)

resultado_chimeraTE <- ggplot( resultado_chimeraTE_table, aes(x= id, y = n, fill= Category))+
  geom_col(position = "dodge", colour = "gray30")+
  scale_fill_paletteer_d ("nationalparkcolors::Arches")+
  scale_y_continuous( expand = c(0,0), limits = c(0,51000))+
  theme_classic()+
  labs(x= "Paciente", y = "Número de eventos", fill = "Categoria")+
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 18),
    
    legend.position = c(0.80, 0.80),
    legend.background = element_rect(
      fill = "white",
      colour = "gray30"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12) )+
  coord_flip()+
  geom_text( aes (label = n), position = position_dodge(width = 0.9), hjust = -0.1)

resultado_chimeraTE

saveRDS(
  resultado_chimeraTE,
  file.path(pasta_processados, "resultado_chimeraTE.rds")
)

ggsave(
  file.path(pasta_resultados, "Resultado_Corrida.png"),
  resultado_chimeraTE,
  width = 8,
  height = 12,
  dpi = 300
)

#gráfico 1 Número de quimeras por pacientes por classificação FULL

resultado_chimeraTE_table_full <- chim_class_full %>%
  count(id, Category)

resultado_chimeraTE_full <- ggplot( resultado_chimeraTE_table_full, aes(x= id, y = n, fill= Category))+
  geom_col(position = "dodge", colour = "gray30")+
  scale_fill_paletteer_d ("nationalparkcolors::Arches")+
  scale_y_continuous( expand = c(0,0), limits = c(0,51000))+
  theme_classic()+
  labs(x= "Paciente", y = "Número de eventos", fill = "Categoria")+
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 18),
    
    legend.position = c(0.80, 0.80),
    legend.background = element_rect(
      fill = "white",
      colour = "gray30"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12) )+
  coord_flip()+
  geom_text( aes (label = n), position = position_dodge(width = 0.9), hjust = -0.1)

resultado_chimeraTE_full

saveRDS(
  resultado_chimeraTE_full,
  file.path(pasta_processados, "resultado_chimeraTE_full.rds")
)

ggsave(
  file.path(pasta_resultados, "Resultado_Corrida_full.png"),
  resultado_chimeraTE,
  width = 8,
  height = 12,
  dpi = 300
)

#### Gráfico 2 ####

set_theme(theme_classic())

###### Classes de TEs pelos tipos de chimeras ######

TE_class_plot <- ggplot(class_counts,  aes(x = reorder(Category, -n), y = n, fill = class_id)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,100000))+ #Arrumar o limits
  geom_col(position = 'dodge', colour = "gray30") +
  theme_classic() +
  labs(
    x = "",
    y = "Número de eventos",
    fill = "TE class")

dev.off()
TE_class_plot


ggsave(
  file.path(pasta_resultados, "TE_class_by_category.png"),
  TE_class_plot,
  width = 8,
  height = 8,
  dpi = 300
)

####Gráfico 3 ####

#Para ver as diferentes famílias de TE, mas muy poluido

#family_counts <- chim_class %>%
#  filter(!is.na(family_id)) %>%
#  group_by(id, Category, family_id) %>%
#  summarise(n = n(), .groups = "drop")

#familias_TE <- ggplot(family_counts, aes(x = id, y =n, fill = family_id)) +
  #  geom_col(position = "dodge", colour = "gray30") +
#  theme_classic()+
#  scale_y_continuous(expand = c(0,0))+
#  labs(title = "Famílias de TE por categoria",
#    x = "",    
#    y = "Número de eventos",    
#    fill = "TE family") +
#  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
#  facet_grid(Category ~ ., scale ='free_y')

#dev.off()
#familias_TE

#ggsave( here(  "Resultados","familias_TE.png" ),familias_TE , width = 16, height = 8, dpi = 300)


####Gráfico 4 ####

top_families <- chim_class %>%
  filter(!is.na(family_id)) %>%
  group_by(Category, family_id) %>%
  summarise(n = n()) %>%
  #group_by(Category) %>%
  #slice_max(n, n = 20, with_ties = FALSE) %>%
  #ungroup()
  slice_max(n, n = 10)


top_familias <- ggplot(top_families, aes(x = reorder(family_id, n), y = n, fill = Category)) +
  geom_col(colour = 'gray40') +
  coord_flip() +
  theme_classic() +
  scale_y_continuous(expand = c(0,0))+
  scale_fill_discrete(palette = "Set2")+
  labs( title = "Top famílias de TE por categoria",   x = "",   y = "Número de eventos" ) +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 18),
    legend.position = c(0.80, 0.20),
    legend.background = element_rect(
      fill = "white",
      colour = "gray30"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12) )
top_familias


ggsave(
  file.path(pasta_resultados, "top_familias.png" ),
  top_familias,
  width = 8,
  height = 6,
  dpi = 300
)

#### Grafico por TE#####

top_TE <- chim_class %>%
  filter(!is.na(family_id)) %>%
  group_by(Category, TE_id) %>%
  summarise(n = n()) %>%
  #group_by(Category) %>%
  #slice_max(n, n = 20, with_ties = FALSE) %>%
  #ungroup()
  slice_max(n, n = 20)


top_TE_g <- ggplot(top_TE, aes(x = reorder(TE_id, n), y = n, fill = Category)) +
  geom_col(colour = 'gray40', linewidth = 0.2) +
  coord_flip() +
  theme_classic() +
  scale_y_continuous(expand = c(0,0))+
  scale_fill_paletteer_d("tvthemes::Alexandrite") +
  labs( title = "Top 20 TEs por categoria",   x = "",   y = "Número de eventos" ) +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 18),
    legend.position = c(0.80, 0.20),
    legend.background = element_rect(
      fill = "white",
      colour = "gray30"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12) )+
  facet_grid(Category ~ ., scale ='free_y')
top_TE_g


ggsave(
  file.path(pasta_resultados, "top_TE_g.png" ),
  top_TE_g,
  width = 8,
  height = 6,
  dpi = 300
)


####Gráfico 5 ####
#problema com cores que precisa ser resolvido

#classes_TE <- chim_class %>%
#  filter(!is.na(class_id)) %>%
#  filter(class_id != "DNA?") %>%
#  filter(class_id != "LTR?") %>%
#  filter(class_id != "RC?") %>%
#  group_by(id, class_id) %>%
#  summarise(n = n())

#TE_class_plot <- ggplot(classes_TE, aes(x = reorder(class_id, -n), y = n,  fill = id)) +
#  geom_col(position = "dodge") +
#  theme_classic() +
#  theme(
#    axis.text.x = element_text(angle = 45, hjust = 1),
#    axis.title.x = element_blank()) +
  #scale_fill_brewer(palette = "Set1") +
#  labs(
#    title = "Classes de TE por paciente",
#    y = "Contagem de eventos",
#    fill = "Paciente")+ 
#  scale_y_continuous(expand = c(0,0))

#TE_class_plot

#aarumar o ggsave se for usar
#ggsave(here( "Resultados", "TE_class_plot.png"),  TE_class_plot,  width = 8,  height = 6,  dpi = 300)

#### ==========================================================================


#install.packages("cowplot")
require(cowplot)

plot_grid(TE_class_plot, TE_class_plot_full, top_familias, align = 'hv', ncol = 1, labels = c("A","B","C"))


####Gráfico BOXPLOT Exonized ####

top_families_exon <- chim_class %>%
  filter(!is.na(family_id)) %>%
  group_by(Category, family_id, id) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n>200) %>%
  filter(Category == "Exonized")

top_exonized<- ggplot(top_families_exon, aes(x = reorder(family_id, n), y = n, fill = Category)) +
  geom_boxplot(colour = 'gray30', fill = 'tomato3') +
 # coord_flip() +
  theme_classic() +
  scale_y_continuous(expand = c(0,0), breaks = c(0,2500,5000,7500,10000, 12500, 15000), limits = c(0,15000))+
  labs( title = "Top famílias de TE Exonized",   x = "",   y = "Número de eventos" ) +
  theme(
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 18),
    
    legend.position = c(0.20, 0.80),
    legend.background = element_rect(
      fill = "white",
      colour = "gray30"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12) )
top_exonized


ggsave(
  file.path(pasta_resultados, "top_familias_2.png"),
  top_families_exon,
  width = 8,
  height = 6, #Conferir tamanho
  dpi = 300
)

####Gráfico BOXPLOT Terminated ####

top_families_terminated <- chim_class %>%
  filter(!is.na(family_id)) %>%
  group_by(Category, family_id, id) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n>50) %>%
  filter(Category == "Terminated")

top_terminated <- ggplot(top_families_terminated, aes(x = reorder(family_id, n), y = n, fill = Category)) +
  geom_boxplot(colour = 'gray30', fill = 'cornflowerblue') +
  # coord_flip() +
  theme_classic() +
  scale_y_continuous(expand = c(0,0), breaks = c(0,250,500,750,1000, 1250, 1500), limits = c(0,1500))+
  labs( title = "Top famílias de TE - Terminated",   x = "",   y = "Número de eventos" ) +
  theme(
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 18),
    
    legend.position = c(0.20, 0.80),
    legend.background = element_rect(
      fill = "white",
      colour = "gray30"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12) )

top_terminated

#conferir
ggsave(
  file.path(pasta_resultados, "top_terminated.png"),
  top_terminated,
  width = 15,
  height = 12,
  dpi = 300
)

####Gráfico BOXPLOT Terminated ####

top_families_initiated <- chim_class %>%
  filter(!is.na(family_id)) %>%
  group_by(Category, family_id, id) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n>3) %>%
  filter(Category == "Initiated")

top_initiated <- ggplot(top_families_initiated, aes(x = reorder(family_id, n), y = n, fill = Category)) +
  geom_boxplot(colour = 'gray30', fill = 'plum3') +
  # coord_flip() +
  theme_classic() +
  scale_y_continuous(expand = c(0,0), breaks = c(0,20,40,60,80, 100), limits = c(0,100))+
  labs( title = "Top famílias de TE - Initiated",   x = "",   y = "Número de eventos" ) +
  theme(
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 18),
    
    legend.position = c(0.20, 0.80),
    legend.background = element_rect(
      fill = "white",
      colour = "gray30"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12) )

top_initiated

#conferir
ggsave(
  file.path(pasta_resultados, "top_initiated.png"),
  top_initiated,
  width = 15,
  height = 12,
  dpi = 300
)

library(cowplot)
grafs <- plot_grid(top_exonized, top_terminated, top_initiated, align = 'hv', ncol = 3, labels = c("A","B","C"))

png("TOP_families_boxplot2.png", res = 300, width = 3600, height = 1520)
grafs
dev.off()
