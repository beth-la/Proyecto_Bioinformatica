# Project
#   Analisis de expresion diferencial:Transcriptional and Chromatin Dynamics of Muscle Regeneration
#   after Severe Trauma

# Author
#   Lopez Angeles B. Elizabeth

# Description
#   Descargando datos del proyecto

# Path

#

## Load recount3 R package
library("recount3")

## Cambiar el URL de recount3 a Amazon (AWS)

getOption(
  "recount3_url",
  "http://duffel.rail.bio/recount3"
)

options(recount3_url = "https://recount-opendata.s3.amazonaws.com/recount3/release")

## Confirmando que se cambió el URL
getOption(
  "recount3_url",
  "http://duffel.rail.bio/recount3"
)

## Vemos los proyectos disponibles de raton
mouse_projects <- available_projects(organism = "mouse")

#Descargando datos del proyecto SRP118914
project_info <- subset(
  mouse_projects,
  project == "SRP118914"
)
rse_gene_SRP118914 <- create_rse(project_info)
assay(rse_gene_SRP118914, "counts") <- compute_read_counts(rse_gene_SRP118914)

# Explorando objeto RSE
rse_gene_SRP118914

# Convirtiendo las cuentas por nucleotidos a cuentas por lectura
assay(rse_gene_SRP118914, "counts") <- compute_read_counts(rse_gene_SRP118914)

# Formateando información del experimento.
rse_gene_SRP118914 <- expand_sra_attributes(rse_gene_SRP118914)
colData(rse_gene_SRP118914)[
  ,
  grepl("^sra_attribute", colnames(colData(rse_gene_SRP118914)))
]

# Pasar de caracter a factor
rse_gene_SRP118914$sra_attribute.cell_type <- factor(rse_gene_SRP118914$sra_attribute.cell_type)
rse_gene_SRP118914$sra_attribute.strain <- factor(rse_gene_SRP118914$sra_attribute.strain)
rse_gene_SRP118914$sra_attribute.source_name <- factor(rse_gene_SRP118914$sra_attribute.source_name)
rse_gene_SRP118914$sra_attribute.tissue <- factor(rse_gene_SRP118914$sra_attribute.tissue)

# Eliminando "horas" de timepoint
# Primero extraemos solo los numeros con expresiones regulares.
just_num <- regmatches(rse_gene_SRP118914$sra_attribute.timepoint, gregexpr("\\d+", rse_gene_SRP118914$sra_attribute.timepoint))

# Los asignamos como numeric.
rse_gene_SRP118914$sra_attribute.timepoint <- as.numeric(just_num)

# Generamos variables para nuestro analisis:
# Muestras early y middle
# Early (3- 24hr luego del estimulo)
# Middle (48 - 168 hr lueago del estimulo)

rse_gene_SRP118914$stage <- factor(ifelse(rse_gene_SRP118914$sra_attribute.timepoint <= 24, "early", "middle"))
table(rse_gene_SRP118914$stage)

# Limpieza de muestras de baja calidad (si es que existen)
# Guardamos nuestros datos
rse_gene_SRP118914_unfiltered <- rse_gene_SRP118914
rse_gene_SRP118914$assigned_gene_prop <- rse_gene_SRP118914$recount_qc.gene_fc_count_all.assigned / rse_gene_SRP118914$recount_qc.gene_fc_count_all.total
summary(rse_gene_SRP118914$assigned_gene_prop)

with(colData(rse_gene_SRP118914), tapply(assigned_gene_prop, stage, summary))

# Revisando calidad de las muestras
hist(rse_gene_SRP118914$assigned_gene_prop)
table(rse_gene_SRP118914$assigned_gene_prop < 0.3)

# Calculando niveles medios de la expresion de los genes en las muestras:
gene_means <- rowMeans(assay(rse_gene_SRP118914, "counts"))
summary(gene_means)

# Utilizamos solo awurellos genes cuyo valor medio de expresion sea mayor a 0.1
rse_gene_SRP118914 <- rse_gene_SRP118914[gene_means > 0.1, ]

## Porcentaje de genes que se retuvieron
round(nrow(rse_gene_SRP118914) / nrow(rse_gene_SRP118914_unfiltered) * 100, 2)
