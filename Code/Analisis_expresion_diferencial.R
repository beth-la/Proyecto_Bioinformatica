# Project
#   Analisis de expresion diferencial:Transcriptional and Chromatin Dynamics of Muscle Regeneration
#   after Severe Trauma

# Author
#   Lopez Angeles B. Elizabeth

# Description
#   Codigo de analisis de expresion diferencial.

# Añadiendo librerias:
library("recount3")
library("limma")
library("edgeR") # BiocManager::install("edgeR", update = FALSE)
library("ggplot2")
library("pheatmap")
library("RColorBrewer")

# Cambienado el URL de recount3 a Amazon (AWS)

getOption(
  "recount3_url",
  "http://duffel.rail.bio/recount3"
)

options(recount3_url = "https://recount-opendata.s3.amazonaws.com/recount3/release")

getOption(
  "recount3_url",
  "http://duffel.rail.bio/recount3"
)

# Proyectos de ratón disponibles en recount3
mouse_projects <- available_projects(organism = "mouse")

# Descargando datos del proyecto SRP118914
project_info <- subset(
  mouse_projects,
  project == "SRP118914"
)

# Creando objeto
rse_gene_SRP118914 <- create_rse(project_info)

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

# Eliminando "horas" de timepoint para tomarlo como un numeric
# Primero extraemos solo los numeros con expresiones regulares.
just_num <- regmatches(rse_gene_SRP118914$sra_attribute.timepoint, gregexpr("\\d+", rse_gene_SRP118914$sra_attribute.timepoint))

# Los asignamos como numeric.
rse_gene_SRP118914$sra_attribute.timepoint <- as.numeric(just_num)

# Generamos variables para nuestro analisis:
# Muestras early y middle
# Early (3- 24hr luego del estimulo)
# Middle (48 - 168 hr luego del estimulo)

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

# Utilizamos solo aquellos genes cuyo valor medio de expresion sea mayor a 0.1
rse_gene_SRP118914 <- rse_gene_SRP118914[gene_means > 0.1, ]

# Porcentaje de genes que se retuvieron
round(nrow(rse_gene_SRP118914) / nrow(rse_gene_SRP118914_unfiltered) * 100, 2)

# Usaremos ExploreModelMatrix para el analisis de nuestro modelo
# estadisitico.
# Se quiere observar la contribucion del tiempo segun el estado (early o middle)
# para el cual se este midiendo la expresion genica.
ExpModelMatrix_SRP118914 <- ExploreModelMatrix::VisualizeDesign(
  sampleData = colData(rse_gene_SRP118914),
  designFormula = stage ~ sra_attribute.timepoint,
  textSizeFitted = 4
)

# Visualizar las imagenes:
cowplot::plot_grid(plotlist = ExpModelMatrix_SRP118914$plotlist)

# Modelo estadistico: Decimos que la variable stage se encuentra explicada
# por el tiempo (timepoint) y la proporcion de genes asignada.
model <- model.matrix( ~ sra_attribute.timepoint + assigned_gene_prop,
                       data = colData(rse_gene_SRP118914)
)
colnames(model)

# Imagenes:
cowplot::plot_grid(plotlist = ExpModelMatrix_SRP118914$plotlist)

# Normalizacion de los datos:

dge <- DGEList(
  counts = assay(rse_gene_SRP118914, "counts"),
  genes = rowData(rse_gene_SRP118914)
)
dge <- calcNormFactors(dge)

# Analisis de expresion diferencial
# Genrando boxplot para analisis de expresion de genes entre los estados tempranos y medios.

ggplot(as.data.frame(colData(rse_gene_SRP118914)), aes(y = assigned_gene_prop, x = stage)) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  ylab("Assigned Gene Prop") +
  xlab("Stage")

# Graficando: Mean-variance trend

vGene <- voom(dge, model, plot = TRUE)

eb_results <- eBayes(lmFit(vGene))

de_results <- topTable(
  eb_results,
  coef = 2,
  number = nrow(rse_gene_SRP118914),
  sort.by = "none"
)
dim(de_results)

# Genes diferencialmente expresados entre early y middle con FDR < 5%
table(de_results$adj.P.Val < 0.05)

# Visualizar datos
plotMA(eb_results, coef = 2)
volcanoplot(eb_results, coef = 2, highlight = 3, names = de_results$gene_name)

de_results[de_results$gene_name %in% c("Postn", "Atp5b", "Gm4735"), ]

# HEATMAP: Los 50 genes mas significativos.
# Extraer valores de los genes de interés utilizando rank
exprs_heatmap <- vGene$E[rank(de_results$adj.P.Val) <= 50, ]

# Creando un dtaframe con información de las muestrasy con nombres de columnas de interes.
df <- as.data.frame(colData(rse_gene_SRP118914)[, c("stage", "sra_attribute.timepoint")])
colnames(df) <- c("STAGE", "HOURS")

# Generamos una funcion con lappy que nos permita buscar el nomre del gene segun el ID que se encuentra asociado a este.
gene_names <- unlist(lapply(row.names(exprs_heatmap), function(id){
  index <- match(id,de_results$gene_id)
  return(de_results$gene_name[index])
}))

# Asegurarnos de que se obtuvieron los nombres correctamente:
length(gene_names) == length(row.names(exprs_heatmap))

# Asignamos los nombres de los genes.
row.names(exprs_heatmap) <- gene_names

pheatmap(
  exprs_heatmap,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = df,
  color = hcl.colors(50, "RdPu")
)

# PERFILES DE EXPRESION

col.group <- df$STAGE
levels(col.group) <- brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)

plotMDS(vGene$E, labels = df$HOURS, col = col.group)

# Obtendremos los counts de los 10 genes mas significativos para relizar un barplot
# Primero obtenemos los genes que se encuentran en las primeras 10 posiciones
# Obtenemos los IDs y generamos un dataframe con esta informacion
gene_ids <- vGene$E[rank(de_results$adj.P.Val) <= 10, ]
gene_ids <- row.names(gene_ids)
counts_first10genes <- data.frame(dge$counts[gene_ids,])

# Renombrando rows de dataframe con los nombres de genes.
gene_10names <- unlist(lapply(gene_ids, function(id){
  index <- match(id,de_results$gene_id)
  return(de_results$gene_name[index])
}))
row.names(counts_first10genes) <- gene_10names

# Renombramos con el tiempo de cada muestra
colnames(counts_first10genes) <- rse_gene_SRP118914$sra_attribute.timepoint

# Redimensionamos el dataframe para una mejor visualizacion de la grafica
counts_first10genes <- apply(counts_first10genes, MARGIN = c(2,1), mean)

# Generamos el barplot
barplot(t(as.matrix(counts_first10genes)),beside=TRUE, legend.text = gene_10names, col = rainbow(11), main = "Expresion diferencial: Top 10")
