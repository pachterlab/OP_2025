library(readr)
library(Seurat)
library(Matrix)
library(tidyverse)
library(tximport)
library(HDF5Array)
library(tidyr)


triple_count_mat <- function(cells, 
                             t2g, 
                             tx=TRUE, 
                             filter=10, 
                             tpm='lengthScaledTPM'
                             ) {
    txi <- tximport(cells$amb_path, type = "kallisto", tx2gene = t2g, countsFromAbundance=tpm, ignoreTxVersion = TRUE, txOut = tx)
    amb_counts <- round(txi$counts)
    rownames(amb_counts) <- as.character(map(strsplit(rownames(amb_counts), split = ".", fixed=TRUE), 1))
    
    txi <- tximport(cells$nac_path, type = "kallisto", tx2gene = t2g, countsFromAbundance=tpm, ignoreTxVersion = TRUE, txOut = tx)
    nac_counts <- round(txi$counts)
    rownames(nac_counts) <- as.character(map(strsplit(rownames(nac_counts), split = ".", fixed=TRUE), 1))
    
    txi <- tximport(cells$mat_path, type = "kallisto", tx2gene = t2g, countsFromAbundance=tpm, ignoreTxVersion = TRUE, txOut = tx)
    mat_counts <- round(txi$counts)
    rownames(mat_counts) <- as.character(map(strsplit(rownames(mat_counts), split = ".", fixed=TRUE), 1))
    
    
    
    SE <- SummarizedExperiment(assays=list(total = mat_counts+nac_counts+amb_counts,
                                           spliced = mat_counts, 
                                           unspliced = nac_counts,
                                           ambiguous = amb_counts),
                               rowData = data.frame(gene_id = rownames(amb_counts),
                                                    has_U_tr = rowSums(nac_counts)>0,
                                                    eff_len = txi$length[,1]),
                               colData = cells)
    
    ensembl <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "rnorvegicus_gene_ensembl")  #* If service is down, try mirror = "useast"
    mt_genes <- biomaRt::getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), filters = 'chromosome_name', values = 'MT', mart = ensembl)
    
    rowData(SE)$mt <- rowData(SE)$gene_id %in% mt_genes['ensembl_gene_id']
    
    if (tx){
        rowData(SE) <- merge(rowData(SE), t2g, 
                         by.x = 'gene_id',
                         by.y='ensembl_transcript_id', all.x=TRUE)
    }else{
        t2g <- t2g[t2g$ensembl_gene_id==t2g$ensembl_transcript_id,]
        rowData(SE) <- merge(rowData(SE), t2g, 
                             by.x = 'gene_id',
                             by.y='ensembl_gene_id', all.x=TRUE)
        
    }
    
    
    sum_data <- data.frame(rowData(SE)) %>% group_by(gene_id) %>% reframe(max(has_U_tr), across())
    rowData(SE)$has_U_tr <- as.logical(unlist(data.frame(sum_data['max(has_U_tr)'])))
    sel = rowSums(assays(SE)$total) > filter
    SE = SE[sel,]
    
    return(SE)}

genes_in_original <- function(gene_list, 
                              og_de_path ='output/original_de.csv', 
                              t2go_path = 'metadata_csvs/t2go.csv', 
                              tx = FALSE){
    og_de = read_csv(og_de_path)
    t2go = read_csv(t2go_path)
    if (tx){
        filtered = t2go %>% filter_at(vars('ensembl_transcript_id'), any_vars(. %in% gene_list))
        gene_names = filtered$X3
    } else{
        filtered = t2go %>% filter_at(vars('ensembl_gene_id'), any_vars(. %in% gene_list))
        gene_names = filtered$X3
    }
    return(og_de %>% filter_at(vars('gene_symbol'), any_vars(. %in% gene_names)))
    }