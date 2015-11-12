#!/usr/bin/Rscript

#YN version4

argspec <- c("aggregate_dups.r
Usage:
    aggregate_dups.r input.tab output.tab
Example:
	Rscript aggregate_dups.r test_matrix.tab output_matrix.tab
Options:
	input matrix (annotated by row and column names)
	output file
	")

read_matrix <- function(in_file){
	header <- strsplit(readLines(con=in_file, n=1), "\t")[[1]]
	cl.cols<- 1:length(header) > 1
	data_matrix.df <- read.delim(in_file, header=TRUE, row.names=NULL, stringsAsFactors=FALSE, na.strings="NA", check.names=FALSE, quote="")
	data_matrix <- as.matrix(data_matrix.df[,cl.cols])
	colnames(data_matrix) <- colnames(data_matrix.df)[cl.cols]
	rownames(data_matrix) <- data_matrix.df[,1]
	return(data_matrix)
}

write_matrix <- function(data_matrix){
	header <- append(c("samples"), colnames(data_matrix))
	if(is.null(colnames(data_matrix)) || length(colnames(data_matrix))==1)
	{
		header <- c("samples", "corr")
	}
	write.table(t(header), stdout(), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
	write.table(data_matrix, stdout(), quote=FALSE, sep="\t", row.names=TRUE, col.names=FALSE)
}

write_matrix_to_file <- function(data_matrix, output_file){
	header <- append(c("samples"), colnames(data_matrix))
	if(is.null(colnames(data_matrix)) || length(colnames(data_matrix))==1)
	{
		header <- c("samples", "corr")
	}
	write.table(t(header), output_file, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
	write.table(data_matrix, output_file, quote=FALSE, sep="\t", row.names=TRUE, col.names=FALSE, append=TRUE)
}

rowMean <- function(i, m){
	subset <- m[rownames(m) == i,]
	combined <- apply(subset, 2, mean)
	return(combined)
}

main <- function(argv) {  
	ptm <- proc.time()
	
	if(length(argv) == 1){
		if(argv==c('--help')){ 
			write(argspec, stderr());
			q();
		}
	}
		
	if(!(length(argv) == 2)){
		write("ERROR: invalid number of arguments is specified", stderr());
		q();
	}
	
	#store command line arguments in variables:
	input_file <- argv[1]
	output_file <- argv[2]
	
	write("Reading the data", stderr());
	data_matrix <- read_matrix(input_file)
	
	write("Aggregating the data", stderr());
	dup.names <- duplicated(rownames(data_matrix))
	names(dup.names) <- rownames(data_matrix)
	data_matrix.true.tmp <- data_matrix[which(dup.names == TRUE),]
	dup.rows <- unique(rownames(data_matrix.true.tmp))
	data_matrix.true <- data_matrix[rownames(data_matrix) %in% dup.rows,]
	data_matrix.false <- data_matrix[!(rownames(data_matrix) %in% dup.rows),]

	if(nrow(data_matrix.true) == 0){
		result.mtrx <- data_matrix.false
		write("WARNING: no duplicate rows were found", stderr());
	}else{
		d.n <- unique(rownames(data_matrix.true))
		combined.rows <- t(rbind(sapply(d.n, rowMean, m=data_matrix.true)))
		result.mtrx <- rbind(data_matrix.false, combined.rows)
	}
	
	write("Writing the data", stderr());
	write_matrix_to_file(result.mtrx, output_file)
	
	write("Execution time: ", stderr());
	write(proc.time() - ptm, stderr());
}

main(commandArgs(TRUE))
