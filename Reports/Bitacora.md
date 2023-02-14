# Bitácora: Proyecto 

López Ángeles Brenda Elizabeth

Línea de código para establecer directorio de trabajo 

```R
usethis::proj_set("C:/Users/Brenda Elizabeth L/Documents/Proyecto_Bioinformatica/")
```

***Miércoles 8 Feb 2023***

- Con la intención de elegir un proyecto disponible en Reacount3 para realizar un análisis de expresión diferencial, se consultó la pagina de [recount3](https://jhubiostatistics.shinyapps.io/recount3-study-explorer/). 
- Se eligió el proyecto SRP118914 "Transcriptional and Chromatin Dynamics of Muscle Regeneration after Severe Trauma".
- Se descargaron los datos desde Recount3.

## Limpieza y normalización de datos

***Jueves 8 Feb 2023***

Una vez habiendo descargado los datos, podemos explorarlos obteniendo que:

```R
rse_gene_SRP118914
'''
class: RangedSummarizedExperiment 
dim: 55421 10 
metadata(8): time_created recount3_version ... annotation recount3_url
assays(2): raw_counts counts
rownames(55421): ENSMUSG00000079800.2 ENSMUSG00000095092.1 ... ENSMUSG00000096850.1
  ENSMUSG00000099871.1
rowData names(11): source type ... havana_gene tag
colnames(10): SRR7184296 SRR7184297 ... SRR7184304 SRR7184305
colData names(177): rail_id external_id ... recount_pred.curated.cell_line BigWigURL
'''
> rse_gene_SRP118914$sra.sample_attributes
'''
 [1] "cell type;;muscle satellite cells (MuSCs)|source_name;;168hr_MuSCs|strain;;C57BL/6J|timepoint;;168 hours|tissue;;tibialis anterior muscle"
 [2] "cell type;;muscle satellite cells (MuSCs)|source_name;;168hr_MuSCs|strain;;C57BL/6J|timepoint;;168 hours|tissue;;tibialis anterior muscle"
 [3] "cell type;;muscle satellite cells (MuSCs)|source_name;;24hr_MuSCs|strain;;C57BL/6J|timepoint;;24 hours|tissue;;tibialis anterior muscle"  
 [4] "cell type;;muscle satellite cells (MuSCs)|source_name;;24hr_MuSCs|strain;;C57BL/6J|timepoint;;24 hours|tissue;;tibialis anterior muscle"  
 [5] "cell type;;muscle satellite cells (MuSCs)|source_name;;3hr_MuSCs|strain;;C57BL/6J|timepoint;;3 hours|tissue;;tibialis anterior muscle"    
 [6] "cell type;;muscle satellite cells (MuSCs)|source_name;;3hr_MuSCs|strain;;C57BL/6J|timepoint;;3 hours|tissue;;tibialis anterior muscle"    
 [7] "cell type;;muscle satellite cells (MuSCs)|source_name;;48hr_MuSCs|strain;;C57BL/6J|timepoint;;48 hours|tissue;;tibialis anterior muscle"  
 [8] "cell type;;muscle satellite cells (MuSCs)|source_name;;48hr_MuSCs|strain;;C57BL/6J|timepoint;;48 hours|tissue;;tibialis anterior muscle"  
 [9] "cell type;;muscle satellite cells (MuSCs)|source_name;;72hr_MuSCs|strain;;C57BL/6J|timepoint;;72 hours|tissue;;tibialis anterior muscle"  
[10] "cell type;;muscle satellite cells (MuSCs)|source_name;;72hr_MuSCs|strain;;C57BL/6J|timepoint;;72 hours|tissue;;tibialis anterior muscle"  
'''
```

- Se estará trabajando con 55421 genes y 10 muestras. 
- No se necesita formatear la información en las muestras ya que todas tienen la misma cantidad de atributos.

```R
# Al correr este bloque de codigo nos damos cuenta de que tenemos un problema, debemos pasar de caracteres a numeric o factor.
colData(rse_gene_SRP118914)[
  ,
  grepl("^sra_attribute", colnames(colData(rse_gene_SRP118914)))
]
'''
DataFrame with 10 rows and 5 columns
           sra_attribute.cell_type sra_attribute.source_name sra_attribute.strain
                       <character>               <character>          <character>
SRR7184296  muscle satellite cel..               168hr_MuSCs             C57BL/6J
SRR7184297  muscle satellite cel..               168hr_MuSCs             C57BL/6J
SRR7184298  muscle satellite cel..                24hr_MuSCs             C57BL/6J
SRR7184299  muscle satellite cel..                24hr_MuSCs             C57BL/6J
SRR7184300  muscle satellite cel..                 3hr_MuSCs             C57BL/6J
SRR7184301  muscle satellite cel..                 3hr_MuSCs             C57BL/6J
SRR7184302  muscle satellite cel..                48hr_MuSCs             C57BL/6J
SRR7184303  muscle satellite cel..                48hr_MuSCs             C57BL/6J
SRR7184304  muscle satellite cel..                72hr_MuSCs             C57BL/6J
SRR7184305  muscle satellite cel..                72hr_MuSCs             C57BL/6J
           sra_attribute.timepoint   sra_attribute.tissue
                       <character>            <character>
SRR7184296               168 hours tibialis anterior mu..
SRR7184297               168 hours tibialis anterior mu..
SRR7184298                24 hours tibialis anterior mu..
SRR7184299                24 hours tibialis anterior mu..
SRR7184300                 3 hours tibialis anterior mu..
SRR7184301                 3 hours tibialis anterior mu..
SRR7184302                48 hours tibialis anterior mu..
SRR7184303                48 hours tibialis anterior mu..
SRR7184304                72 hours tibialis anterior mu..
SRR7184305                72 hours tibialis anterior mu..
'''

# Explorando cada una de las columnas tambien notamos que el tiempo en todas las muestras fue medido en horas, por lo que se decide, para fines mas practicos y anticipando cualquier problema posterior, eliminar la palabra "horas" y tomar a sra_attribute.timepoint como tipo numeric. 

# Hacemos lo anterior con expresiones regulares.
just_num <- regmatches(rse_gene_SRP118914$sra_attribute.timepoint, gregexpr("\\d+", rse_gene_SRP118914$sra_attribute.timepoint))

# Los asignamos como numeric.
rse_gene_SRP118914$sra_attribute.timepoint <- as.numeric(just_num)
```

En los datos de las muestras que se tiene, únicamente hay de 3 a 72 horas después de causar el estimulo en los ratones, estos, según el articulo del cual obtuvimos los datos, disponible en el [repositorio/Reports](https://github.com/beth-la/Proyecto_Bioinformatica/blob/master/Reports/Associated_paper.pdf) corresponden a las etapas tempranas y medias, siendo estas aquellas con las que trabajaremos en este estudio de análisis diferencial. En el estudio citado, se trabaja de 3 a 672 horas. 

```R
> table(rse_gene_SRP118914$stage)
 early middle 
     4      6 
```

***Para eliminar muestras de mala calidad:***

```R
rse_gene_SRP118914$assigned_gene_prop <- rse_gene_SRP118914$recount_qc.gene_fc_count_all.assigned / rse_gene_SRP118914$recount_qc.gene_fc_count_all.total
summary(rse_gene_SRP118914$assigned_gene_prop)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.6786  0.6802  0.6916  0.7251  0.7809  0.7939 

with(colData(rse_gene_SRP118914), tapply(assigned_gene_prop, stage, summary))
$early
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.6800  0.6805  0.6859  0.6860  0.6913  0.6921 

$middle
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.6786  0.7045  0.7809  0.7511  0.7904  0.7939

# Revisando calidad de las muestras
> rse_gene_SRP118914$assigned_gene_prop > 0.3
 [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
hist(rse_gene_SRP118914$assigned_gene_prop)
```

![image-20230209210130600](C:\Users\Brenda Elizabeth L\AppData\Roaming\Typora\typora-user-images\image-20230209210130600.png)

***Filtrado de genes:***

```R
# Una vez eliminados aquellos genes cuyo valor medio de expresion era menor a 0.1, observamos la cantidad de genes retenidos.
# Eliminamos genes:
gene_means <- rowMeans(assay(rse_gene_SRP118914, "counts"))
rse_gene_SRP118914 <- rse_gene_SRP118914[gene_means > 0.1, ]
round(nrow(rse_gene_SRP118914) / nrow(rse_gene_SRP118914_unfiltered) * 100, 2)
[1] 48.46

# Entonces:
> dim(rse_gene_SRP118914_unfiltered)
[1] 55421    10
> dim(rse_gene_SRP118914)
[1] 26857    10
# De comenzar con 55421 genes, terminaremos trabajando unicamente con 26857
```

***Modelo estadístico***

El set de datos que se esta analizando cuenta con diferentes características para cada una de las muestras a analizar, debemos tener en claro que luego del filtrado, nos encontraremos trabajando únicamente con 26857 genes y con 10 muestras, cada muestra cuenta con:

- cell type: muscle satellite cells (MuSCs)
- source_name: 168hr_MuSCs
- strain: C57BL/6
- timepoint: (3, 24, 48, 72, 168) hours
- tissue: tibialis anterior muscle

Nuestro modelo estadisitco, por lo tanto, se encuentra asociado al timepoint (la cual es la unica variable que cambia para cada muestra) y que será dividida en dos clasificaciones distintas: 

- Early: 3 - 24 hrs luego del estimulo
- Middle: 48 - 168 luego del estimulo 

```R
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
# por el (tiempo timepoint) y la proporsion de genes asignada.
model <- model.matrix( ~ sra_attribute.timepoint + assigned_gene_prop,
                      data = colData(rse_gene_SRP118914)
)

# Normalización de datos
library("limma")
library("edgeR") # BiocManager::install("edgeR", update = FALSE)
dge <- DGEList(
  counts = assay(rse_gene_SRP118914, "counts"),
  genes = rowData(rse_gene_SRP118914)
)
dge <- calcNormFactors(dge)
```

***Lista de tareas: Viernes 9 Feb 2023***

- [x] Actualizar y mejorar README del proyecto, agregar link hacia recount3.
- [x] Añadir licencia al repositorio.
- [x] Realizar pasos iniciales del análisis de expresión.

## Análisis de expresión diferencial

```R
# Comenzamos generando un boxplot
library("ggplot2")
ggplot(as.data.frame(colData(rse_gene_SRP118914)), aes(y = assigned_gene_prop, x = stage)) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  ylab("Assigned Gene Prop") +
  xlab("Stage")
```

![image-20230211200304176](C:\Users\Brenda Elizabeth L\AppData\Roaming\Typora\typora-user-images\image-20230211200304176.png)

```R
#Se decidió tomar las 10 primeras posiciones de genes significativos para realizar un barplot
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

```

# Cambios pedientes

- [x] Hacer barplots con los genes
- [x] Agregar  citas hacia códigos consultados y fuentes de información como recount3.
- [x] Tomar genes mas significativos y realizar análisis 
