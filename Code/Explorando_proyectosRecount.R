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
#PROYECTO 147
