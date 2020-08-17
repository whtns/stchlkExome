library(httr)
library(jsonlite)
library(xml2)

server <- "https://grch37.rest.ensembl.org"
ext <- "/vep/homo_sapiens/region"
r <- POST(paste(server, ext, sep = ""), content_type("application/json"), accept("application/json"), body = '{ "variants" : ["21  26960070  . G A . . .", "21  26965146  . GG A . . ." ] }')

stop_for_status(r)

# use this if you get a simple nested list back, otherwise inspect its structure

# head(data.frame(t(sapply(content(r),c))))
test0 <- head(fromJSON(toJSON(content(r))))
