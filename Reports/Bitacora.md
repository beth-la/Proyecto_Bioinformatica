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

## Análisis de expresión diferencial 

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
```

