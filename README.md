# Análisis de expresión diferencial: "Transcriptional and Chromatin Dynamics of Muscle Regeneration after Severe Trauma"

Análisis de expresión diferencial utilizando datos del estudio: *Transcriptional and Chromatin Dynamics of Muscle Regeneration after Severe Trauma*, obtenidos desde [recount3](https://jhubiostatistics.shinyapps.io/recount3-study-explorer/). 

## Contents

- Code:
   - Analisis_expresion_diferencial.R: Código con el cual se llevó a cabo el analisis de expresión diferencial.
   - Explorando_proyectosRecountR: Código en el que se exploraron los proyectos disponibles en recount3
   
- Reports 
   - Associated_paper: *Transcriptional and Chromatin Dynamics of Muscle Regeneration after Severe Trauma*
   - Bitacora: Bitacora asociada al codigo que se encuentra en Proyecto_bioinformatica/Code/Analisis_expresion_diferencial.R, unicamente con el proposito de organizar y tener mas control sobre los pasos del proyecto. 
   - Report_results.md: Reporte de análsisi de los resultados en formato .md
   - Experimental-data_illustration: Ilustación de los datos para el análisis de expresión diferencial.
   
- Plots 
    - Heatmap.png: heatmap de los 50 genes mas significativos diferencialmente expresados.
    - ExploreModelMatrix.png: Output resultante de Explore Model matrix para el modelo estadistico.
    - Histogram_rse_gene_SRP118914$assigned_gene_prop: Histograma de expresion de genes en estados early y middle
    - Boxplot_early-middle-stage_gene-prop.png: Boxplot de expresión para condiciones early y middle.
    - MDS_plot: Gráfica MDS.
    - Top_10genes_couts: Grafica de barras que muestra la expresión del top 10 de genes mas significativos en contraste con el tiempo en el que fueron tomadas las muestras. 
    - Volcano_plot: Volcano plot de los genes diferencialmente expresados.
    - Mean_variance
## Experimental data 
![experimental-data_ilustration](https://user-images.githubusercontent.com/100377667/218641062-63b562a5-5ad5-45a8-acf7-86fdef1c5806.png)
>>>>>>> 30b9a3e11fcb2bda2cb4f38ba2c0d02ef94e9fa6

## Sources

Este proyecto se realizó siguiendo el protocolo del Dr. Leonardo Collado Torres disponible en [github](https://github.com/lcolladotor/rnaseq_LCG-UNAM_2023).

- Collado-Torres L (2023). *Explore and download data from the recount3 project*. doi: [10.18129/B9.bioc.recount3](https://doi.org/10.18129/B9.bioc.recount3), https://github.com/LieberInstitute/recount3 - R package version 1.8.0, http://www.bioconductor.org/packages/recount3.
- Wilks C, Zheng SC, Chen FY, Charles R, Solomon B, Ling JP, Imada EL, Zhang D, Joseph L, Leek JT, Jaffe AE, Nellore A, Collado-Torres L, Hansen KD, Langmead B (2021). “recount3: summaries and queries for large-scale RNA-seq expression and splicing.” *Genome Biol*. doi: [10.1186/s13059-021-02533-6](https://doi.org/10.1186/s13059-021-02533-6), https://doi.org/10.1186/s13059-021-02533-6.
- Aguilar, C. A., Pop, R., Shcherbina, A., Watts, A., Matheny, R. W., Cacchiarelli, D., ... & Meissner, A. (2016). Transcriptional and chromatin dynamics of muscle regeneration after severe trauma. *Stem cell reports*, *7*(5), 983-997.

*Los datos utilizados para el analisis de expresión diferencial en el presente proyecto no me pertenecen, fueron obtenidos desde recount3.*

## LICENSE 

MIT License
