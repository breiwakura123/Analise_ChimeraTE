## =========================================================
## Painel A - "Manual curation" (estilo IGV) para 1 amostra
## Tracks: anotacao genica (exons/UTRs) + cobertura RNA-seq
##         + anotacao do TE (barra vermelha)
## =========================================================

#TESTE

## ---------------------------------------------------------
## 0) INSTALACAO (rodar uma vez)
## ---------------------------------------------------------
# if (!require("BiocManager")) install.packages("BiocManager")
# BiocManager::install(c("Gviz", "GenomicFeatures", "rtracklayer",
#                         "Rsamtools", "GenomicAlignments"))

library(Gviz)
library(GenomicFeatures)
library(rtracklayer)
library(Rsamtools)
library(GenomicAlignments)
library(GenomeInfoDb)
library(paletteer)
library(here)

## ---------------------------------------------------------
## 0b) PALETA DE CORES (paletteer)
## ---------------------------------------------------------
## Extrai as cores da paleta discreta "tvthemes::Alexandrite"
## e atribui manualmente qual cor vai para cada track.
## Rode paletteer_d("tvthemes::Alexandrite") para ver quantas
## cores tem e escolher os indices que preferir.

alexandrite <- paletteer_d("tvthemes::Alexandrite")
print(alexandrite)  # confira quantas cores existem e a ordem

gene_color <- as.character(alexandrite[1])   # cor do modelo genico
cov_color  <- as.character(alexandrite[2])   # cor da cobertura RNA-seq
te_color   <- as.character(alexandrite[3])   # cor do TE

## ---------------------------------------------------------
## 1) ARQUIVOS NECESSARIOS (o que voce precisa ter em maos)
## ---------------------------------------------------------
## a) BAM alinhado e ORDENADO/INDEXADO da amostra
##    (ex: saida do STAR do proprio ChimeraTE, ou um BAM seu)
##    -> precisa existir tambem o arquivo .bai na MESMA PASTA,
##       com o mesmo nome + ".bai" (gerado via `samtools index`)
bam_file <- here("Dados", "bam_bai", "7877_T_dup_bqsr.bam")

## b) Anotacao de genes em GTF/GFF3 (a mesma usada no ChimeraTE)
gtf_file <- here("Dados", "Referências", "GCF_000001405.40_GRCh38.p14_genomic.gtf")

## c) Anotacao de TEs em formato GTF
##    (saida do RepeatMasker convertida p/ GTF, ou o TE.gtf
##    que voce ja usa como input do ChimeraTE)
te_gtf_file <- here("Dados", "Referências", "hg38_rmsk_TE.gtf")

## d) O locus de interesse: cromossomo + coordenadas da regiao
##    a ser plotada. Exemplo abaixo: gene TP53 (GRCh38/hg38),
##    com uma margem de ~2 kb para cada lado pra dar contexto.
##    Coordenadas do gene em si: chr17:7,668,421-7,687,490 (fita -)
##
##    IMPORTANTE: seu BAM usa nomenclatura UCSC (chr17), enquanto
##    seus GTFs usam nomenclatura RefSeq (NC_000017.11). Vamos
##    padronizar tudo para UCSC (chr17) mais abaixo, entao aqui
##    ja usamos o nome que bate com o BAM:
chr_of_interest    <- "chr17"
region_start       <- 7666000
region_end         <- 7689500
gene_id_of_interest <- "TP53"   # so para o titulo/legenda

## ---------------------------------------------------------
## 2) CHECAGENS BASICAS
## ---------------------------------------------------------
stopifnot(file.exists(bam_file))
stopifnot(file.exists(paste0(bam_file, ".bai")) |
            file.exists(sub("\\.bam$", ".bai", bam_file)))
stopifnot(file.exists(gtf_file))
stopifnot(file.exists(te_gtf_file))

## ---------------------------------------------------------
## 3) TRACK DE ANOTACAO GENICA (exons/UTRs, setas de sentido)
## ---------------------------------------------------------
txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")

## Renomeia os cromossomos de RefSeq (NC_000017.11) para UCSC (chr17),
## para bater com a nomenclatura do BAM
seqlevelsStyle(txdb) <- "UCSC"

geneTrack <- GeneRegionTrack(
  txdb,
  chromosome = chr_of_interest,
  start = region_start,
  end   = region_end,
  name  = "Gene model",
  fill  = gene_color,
  col   = gene_color,
  fontcolor.group = "black",
  transcriptAnnotation = "gene"
)

## ---------------------------------------------------------
## 4) TRACK DE COBERTURA DE RNA-seq (histograma preto)
## ---------------------------------------------------------
covTrack <- DataTrack(
  range = bam_file,
  type = "histogram",
  chromosome = chr_of_interest,
  start = region_start,
  end   = region_end,
  name  = "RNA-seq coverage",
  col.histogram = cov_color,
  fill.histogram = cov_color,
  window = -1              # resolucao maxima (por base)
)

## ---------------------------------------------------------
## 5) TRACK DO TE (barra vermelha)
## ---------------------------------------------------------
te_gr_all <- import(te_gtf_file, format = "GTF")

## Renomeia os cromossomos de RefSeq (NC_000017.11) para UCSC (chr17),
## mesma logica aplicada ao txdb acima
seqlevelsStyle(te_gr_all) <- "UCSC"

## O GTF do TE costuma trazer o nome/familia em algum atributo
## (gene_id, transcript_id, repeat_name, etc). Confira com:
##   colnames(mcols(te_gr_all))
##   head(mcols(te_gr_all))
## e ajuste a coluna abaixo (usada so como rotulo na track):
te_name_col <- "gene_id"   # <-- troque se o nome do atributo for outro

te_gr <- subsetByOverlaps(
  te_gr_all,
  GRanges(chr_of_interest, IRanges(region_start, region_end))
)

if (te_name_col %in% colnames(mcols(te_gr))) {
  names(te_gr) <- mcols(te_gr)[[te_name_col]]
}

teTrack <- AnnotationTrack(
  te_gr,
  chromosome = chr_of_interest,
  name  = "TE",
  fill  = te_color,
  col   = te_color,
  stacking = "dense"
)

## ---------------------------------------------------------
## 6) EIXO DE COORDENADAS
## ---------------------------------------------------------
axisTrack <- GenomeAxisTrack()

## ---------------------------------------------------------
## 7) PLOTAR TUDO JUNTO (ordem de cima p/ baixo = igual figura)
## ---------------------------------------------------------
png(paste0("panelA_", gene_id_of_interest, ".png"),
    width = 1600, height = 500, res = 150)

plotTracks(
  list(axisTrack, geneTrack, covTrack, teTrack),
  chromosome = chr_of_interest,
  from = region_start,
  to   = region_end,
  sizes = c(1, 1.5, 2.5, 0.6),   # proporcao de altura de cada track
  background.title = "white",
  col.title = "black",
  fontsize = 12,
  main = gene_id_of_interest
)

dev.off()

## Tambem plota direto na tela (RStudio / interativo)
plotTracks(
  list(axisTrack, geneTrack, covTrack, teTrack),
  chromosome = chr_of_interest,
  from = region_start,
  to   = region_end,
  sizes = c(1, 1.5, 2.5, 0.6),
  background.title = "white",
  col.title = "black",
  main = gene_id_of_interest
)