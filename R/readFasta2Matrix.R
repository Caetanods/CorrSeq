get.seq <- function(X, a, b){
    seq.list <- X[a:b]
    seq.string <- do.call(c, seq.list)
    return( seq.string )
}

##' Reads a FASTA file of any type and returns a matrix with the species as rows and the sites as columns.
##'
##' Each site of the FASTA file is transformed into a "character" in the matrix. The function ony works for sequences of the same length (i.e., alignments).
##' @title Read a FASTA file and returns a matrix
##' @param file the FASTA file to be read.
##' @return a matrix
##' @author daniel
##' @export
readFasta2Matrix <- function(file){
    raw.matrix <- read.csv(file = file, as.is = TRUE, header = FALSE)
    parsed <- apply(raw.matrix, 1, function(x) strsplit(x, split="")[[1]])
    spp.name <- sapply(parsed, function(x) x[1] == ">")
    species.list <- parsed[spp.name]
    species <- sapply(species.list, function(x) paste0(x[-1], collapse="") )
    spp.id <- c(which(spp.name), length(spp.name))
    seq.id.table <- cbind(spp.id[1:(length(spp.id)-1)]+1, c(spp.id[2:(length(spp.id)-1)]-1, length(spp.name)))
    seq.data <- t( sapply(1:nrow(seq.id.table), function(x) get.seq(X=parsed, a=seq.id.table[x,1], b=seq.id.table[x,2]) ) )
    rownames( seq.data ) <- species
    colnames( seq.data ) <- paste0("pos", 1:ncol(seq.data))
    return( seq.data )   
}
