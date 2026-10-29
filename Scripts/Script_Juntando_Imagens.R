#===========================================================================
### Juntando gráficos ####
#===========================================================================
library(patchwork)
library(gridExtra)
library(grid)

####Número de eventos/ Eficiência ChimeraTE####
resumo_por_tipo <- read.csv (here("Resultados", "Geral", "resumo_por_tipo.csv"))
resultado_chimeraTE_ <- readRDS(here("Dados", "Processados", "resultado_chimeraTE.rds"))

# Transformar a tabela em uma tabela gráfica
tabela_tipo <- tableGrob(
  resumo_por_tipo,
  rows = NULL
)

# Juntar tabela + gráfico
figura_1_final <- wrap_elements(tabela_tipo) +
  resultado_chimeraTE_ +
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

figura_1_final

#Salvando
ggsave( here("Resultados", "Figura_1_final_2.png"),
  figura_1_final,
  width = 16,
  height = 10,
  dpi = 300)
