#!/usr/bin/Rscript

#YN 20130123

#usage, options and doc goes here
argspec <- c("filter_rows_on_sd.r - filters rows based on specified sd
Usage:
    filter_rows_on_sd.r input_data_file cutoff_type cutoff filtered_output_file
Example:
	Rscript filter_rows_on_sd.r input_data.tab TOP 6000 filtered_output.tab
	Rscript filter_rows_on_sd.r input_data.tab SD 0.03 filtered_output.tab
Options:
	input matrix (annotated by row and column names)
	type of cutoff (TOP, SD)
	standard diviation cutoff or the top n cutoff
	output matrix file")

read_matrix <- function(in_file){
	header <- strsplit(readLines(con=in_file, n=1), "\t")[[1]]
	cl.cols<- 1:length(header) > 1
	data_matrix.df <- read.delim(in_file, header=TRUE, row.names=NULL, stringsAsFactors=FALSE, na.strings="NA", check.names=FALSE)
	data_matrix <- as.matrix(data_matrix.df[,cl.cols])
	rownames(data_matrix) <- data_matrix.df[,1]
	return(data_matrix)
}

write_matrix <- function(data_matrix, file_name){
	header <- append(c("gene"), colnames(data_matrix))
	write.table(t(header), file_name, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
	write.table(data_matrix, file_name, quote=FALSE, sep="\t", row.names=TRUE, col.names=FALSE, append=TRUE)
}

filter_low_sd <- function(mtrx, cutoff_type, cutoff){
	if(cutoff_type == "SD"){
		mtrx.new <- c()
		rownames.new <- c()
		for(i in 1:nrow(mtrx)){
			row_vals <- mtrx[i,]
			row_vals.sd <- sd(row_vals)
			if(row_vals.sd >= cutoff){
				rownames.new <- append(rownames.new, rownames(mtrx)[i])
				if(i == 1){
					mtrx.new <- row_vals
				}else{
					mtrx.new <- rbind(mtrx.new, row_vals)
				}
			}
		}
		rownames(mtrx.new) <- rownames.new
		colnames(mtrx.new) <- colnames(mtrx)
	}else if(cutoff_type == "TOP"){
		sds <- abs(apply(mtrx, 1, sd))
		sds.sorted <- sds[order(sds,decreasing = TRUE)]
		sds.sorted.selected <- sds.sorted[1:cutoff]
		mtrx.new <- mtrx[rownames(mtrx) %in% names(sds.sorted.selected), ]
	}
	return(mtrx.new)
}

main <- function(argv) {  
	valid_cutoffs <- c("TOP", "SD")	
	
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
	input_file <- argv[1]
	cutoff_type <- toupper(argv[2])
	cutoff <- as.numeric(argv[3])
	output_file <- argv[4]
	
	if(!(cutoff_type %in% valid_cutoffs)){
		write("ERROR: invalid type of cutoff. Valid types are:", stderr());
		write(valid_cutoffs, stderr());
		q();
	}	
	
	#read the input file(s):
	data_matrix <- read_matrix(input_file)
	data_matrix.new <- filter_low_sd(data_matrix, cutoff_type, cutoff)
	write_matrix(data_matrix.new, output_file)
}

main(commandArgs(TRUE))
