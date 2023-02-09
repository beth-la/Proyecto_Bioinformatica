# Buscando Proyecto
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

#Probando proyecto 147 SRP118914
project_info <- subset(
  mouse_projects,
  project == "SRP118914"
)
rse_gene_SRP118914 <- create_rse(project_info)
assay(rse_gene_SRP118914, "counts") <- compute_read_counts(rse_gene_SRP118914)

# Se decidio usar el proyecto anterior
