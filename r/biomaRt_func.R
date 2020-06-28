library(biomaRt)
library(dplyr)


ensembl_hsapiens   = useEnsembl(biomart="ensembl", dataset=c("hsapiens_gene_ensembl"), host = "www.ensembl.org")
ensembl_oaries     = useEnsembl(biomart="ensembl", dataset=c("oaries_gene_ensembl"), host = "www.ensembl.org")
ensembl_btaurus    = useEnsembl(biomart="ensembl", dataset=c("btaurus_gene_ensembl"), host = "www.ensembl.org")
ensembl_mmusculus  = useEnsembl(biomart="ensembl", dataset=c("mmusculus_gene_ensembl"), host = "www.ensembl.org")
ensembl_rnorvegicus  = useEnsembl(biomart="ensembl", dataset=c("rnorvegicus_gene_ensembl"), host = "www.ensembl.org")


datasets <- listDatasets(useMart("ensembl"))
head(datasets)

getGeneGOs <- function(gene_names){
  
  temptab <- list()
  
  gene_names <- na.omit(gene_names)
  
  try(temptab[[1]] <- cbind(Species = "hsapiens",
                            getBM(attributes=c('ensembl_gene_id', 'external_gene_name','description', 'phenotype_description', 'go_id', 'name_1006', 'definition_1006'),
                                  filters = 'external_gene_name',
                                  values = gene_names, mart = ensembl_hsapiens)))
  try(temptab[[2]] <- cbind(Species = "oaries",
                            getBM(attributes=c('external_gene_name','description', 'phenotype_description', 'go_id', 'name_1006', 'definition_1006'),
                                  filters = 'external_gene_name',
                                  values = gene_names, mart = ensembl_oaries)))
  try(temptab[[3]] <- cbind(Species = "btaurus",
                            getBM(attributes=c('external_gene_name','description', 'phenotype_description', 'go_id', 'name_1006', 'definition_1006'),
                                  filters = 'external_gene_name',
                                  values = gene_names, mart = ensembl_btaurus)))
  try(temptab[[4]] <- cbind(Species = "mmusculus",
                            getBM(attributes=c('external_gene_name','description', 'phenotype_description', 'go_id', 'name_1006', 'definition_1006'),
                                  filters = 'external_gene_name',
                                  values = gene_names, mart = ensembl_mmusculus)))
  
  temptab <- bind_rows(temptab)
  
}

getOrthoIDs <- function(gene_ids){
  
  gene.orthos <- NULL
  
  # mice
  x <- getLDS(attributes=c("ensembl_gene_id"),
              filters="ensembl_gene_id", values=gene_ids, mart=ensembl_oaries,
              attributesL=c("ensembl_gene_id"), martL=ensembl_mmusculus)
  
  names(x) <- c("sheep_gene_id", "ensembl_gene_id")
  x$Species <- "mmmusculus"
  
  x <- join(x, getBM(attributes=c('ensembl_gene_id', 'external_gene_name'),
                     filters = 'ensembl_gene_id',
                     values = x$ensembl_gene_id, mart = ensembl_mmusculus))
  
  gene.orthos <- rbind(gene.orthos, x)
  rm(x)
  
  # humans
  x <- getLDS(attributes=c("ensembl_gene_id"),
              filters="ensembl_gene_id", values=gene_ids, mart=ensembl_oaries,
              attributesL=c("ensembl_gene_id"), martL=ensembl_hsapiens)
  
  names(x) <- c("sheep_gene_id", "ensembl_gene_id")
  x$Species <- "hsapiens"
  
  x <- join(x, getBM(attributes=c('ensembl_gene_id', 'external_gene_name'),
                     filters = 'ensembl_gene_id',
                     values = x$ensembl_gene_id, mart = ensembl_hsapiens))
  
  gene.orthos <- rbind(gene.orthos, x)
  rm(x)
  
  # cattle
  x <- getLDS(attributes=c("ensembl_gene_id"),
              filters="ensembl_gene_id", values=gene_ids, mart=ensembl_oaries,
              attributesL=c("ensembl_gene_id"), martL=ensembl_btaurus)
  
  names(x) <- c("sheep_gene_id", "ensembl_gene_id")
  x$Species <- "btaurus"
  
  x <- join(x, getBM(attributes=c('ensembl_gene_id', 'external_gene_name'),
                     filters = 'ensembl_gene_id',
                     values = x$ensembl_gene_id, mart = ensembl_btaurus))
  
  gene.orthos <- rbind(gene.orthos, x)
  rm(x)
  
  #~~ Format the table:
  
  x <- gene.orthos[,c("sheep_gene_id", "external_gene_name")]
  x$external_gene_name <- toupper(x$external_gene_name)
  x <- unique(x)
  x <- subset(x, external_gene_name != "")
  x
  
}

getSheepGOs <- function(chrid, chrstart, chrstop){
  
  message("Querying sheep...")
  x <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'description',
                            'phenotype_description', 'go_id', 'name_1006', 'definition_1006',
                            'chromosome_name', 'start_position', 'end_position', 'strand'), 
             filters = c('chromosome_name','start','end'),
             values = list(chrid, chrstart, chrstop), 
             mart = ensembl_oaries, quote = "")
  x$Species <- "oaries"
  
  #~~ Mice
  message("Querying mouse...")
  
  x.mice <- getLDS(attributes=c("ensembl_gene_id"),
                   filters="ensembl_gene_id", values=unique(x$ensembl_gene_id), mart=ensembl_oaries,
                   attributesL=c("ensembl_gene_id"), martL=ensembl_mmusculus)
  
  names(x.mice) <- c("sheep_gene_id", "ensembl_gene_id")
  
  x.mice <- join(x.mice, getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description',
                                            'phenotype_description', 'go_id', 'name_1006', 'definition_1006',
                                            'chromosome_name', 'start_position', 'end_position', 'strand'),
                               filters = 'ensembl_gene_id',
                               values = x.mice$ensembl_gene_id, mart = ensembl_mmusculus))
  
  #~~ Humans
  message("Querying humans...")
  
  x.humans <- getLDS(attributes=c("ensembl_gene_id"),
                     filters="ensembl_gene_id", values=unique(x$ensembl_gene_id), mart=ensembl_oaries,
                     attributesL=c("ensembl_gene_id"), martL=ensembl_hsapiens)
  
  names(x.humans) <- c("sheep_gene_id", "ensembl_gene_id")
  
  x.humans <- join(x.humans, getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description',
                                                'phenotype_description', 'go_id', 'name_1006', 'definition_1006',
                                                'chromosome_name', 'start_position', 'end_position', 'strand'),
                                   filters = 'ensembl_gene_id',
                                   values = x.humans$ensembl_gene_id, mart = ensembl_hsapiens))
  
  #~~ Cattle
  message("Querying cattle...")
  
  x.cow <- getLDS(attributes=c("ensembl_gene_id"),
                  filters="ensembl_gene_id", values=unique(x$ensembl_gene_id), mart=ensembl_oaries,
                  attributesL=c("ensembl_gene_id"), martL=ensembl_btaurus)
  
  names(x.cow) <- c("sheep_gene_id", "ensembl_gene_id")
  
  x.cow <- join(x.cow, getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description',
                                          'phenotype_description', 'go_id', 'name_1006', 'definition_1006',
                                          'chromosome_name', 'start_position', 'end_position', 'strand'),
                             filters = 'ensembl_gene_id',
                             values = x.cow$ensembl_gene_id, mart = ensembl_btaurus))
  
  #~~ Get sheep positions
  
  x$sheep_gene_id <- x$ensembl_gene_id
  
  x.for.joining <- subset(x, select = c(sheep_gene_id, chromosome_name, start_position, end_position))
  
  names(x.for.joining)[2:4] <- c("sheep_chromosome_name", "sheep_start_position", "sheep_stop_position")
  
  x <- rbind(x,
             cbind(x.mice, Species = "mmusculus"),
             cbind(x.humans, Species = "hsapiens"),
             cbind(x.cow, Species = "btaurus"))
  
  x <- join(x, x.for.joining)
  
  x <- unique(x)
  message("...done.")
  
  x
  
}



getJustSheepGOs <- function(chrid, chrstart, chrstop){
  
  x <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'description',
                            'phenotype_description', 'go_id', 'name_1006', #'definition_1006',
                            'chromosome_name', 'start_position', 'end_position', 'strand'), 
             filters = c('chromosome_name','start','end'),
             values = list(chrid, chrstart, chrstop), 
             mart = ensembl_oaries, quote = "")

  x
  
}
