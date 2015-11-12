#!/usr/bin/Rscript

#YN version1
argspec <- c("extract_matrix_columns.r
Usage:
    extract_matrix_columns.r cmap.tab output.tab columns.tab include_flag
Example:
	Rscript extract_matrix_columns.r cmap.tab output_matrix.tab columns.tab TRUE
Options:
	input matrix
	output file name
	file name that contains the list of row headers
	TRUE or FALSE flag to include or exclude the samples in the list")

read_matrix <- function(in_file){
	header <- strsplit(readLines(con=in_file, n=1), "\t")[[1]]
	cl.cols<- 1:length(header) > 1
	data_matrix.df <- read.delim(in_file, header=TRUE, row.names=NULL, stringsAsFactors=FALSE, na.strings="NA", check.names=FALSE, quote="")
	data_matrix <- as.matrix(data_matrix.df[,cl.cols])
	rownames(data_matrix) <- data_matrix.df[,1]
	return(data_matrix)
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

main <- function(argv) {  
	ptm <- proc.time()
	
	if(length(argv) == 1){
		if(argv==c('--help')){ 
			write(argspec, stderr());
			q();
		}
	}
		
	if(!(length(argv) == 4)){
		write("ERROR: invalid number of arguments is specified", stderr());
		q();
	}
	
	#store command line arguments in variables:
	input_file1 <- argv[1]
	output_file <- argv[2]
	druglist_file <- argv[3]
	include_flag_str <- toupper(argv[4])
	if(!(include_flag_str %in% c("TRUE", "FALSE"))){
		write("ERROR: include flag must be TRUE or FALSE", stderr());
		q();
	}
	include_flag <- FALSE
	if(include_flag_str == "TRUE"){
		include_flag <- TRUE
	}	
	
	#read the input file(s):
	data_matrix <- read_matrix(input_file1)
	
	druglist <- read.delim(druglist_file, sep="\t", header=FALSE, quote = "",na.strings = "NA",fill=TRUE,as.is = TRUE)$V1
	if(include_flag == TRUE){
		data_matrix <- data_matrix[,colnames(data_matrix) %in% druglist]
	}else{
		data_matrix <- data_matrix[,!(colnames(data_matrix) %in% druglist)]
	}
	write_matrix_to_file(data_matrix, output_file)

	write("Execution time: ", stderr());
	write(proc.time() - ptm, stderr());
}

main(commandArgs(TRUE))
