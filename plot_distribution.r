#!/usr/bin/Rscript

#usage, options and doc goes here
argspec <- c("plot_distribution.r - plots distribution of the value in the list or the matrix. Assumes the first line and the first column are annotations.
Usage:
    plot_distribution.r input_matrix.tab output_file.pdf
Example:
	Rscript plot_distribution.r input_matrix.tab output_file.pdf
Options:
	input file name
	output file name")

read_matrix <- function(in_file){
	header <- strsplit(readLines(con=in_file, n=1), "\t")[[1]]
	cl.cols<- 1:length(header) > 1
	data_matrix.df <- read.delim(in_file, header=TRUE, row.names=NULL, stringsAsFactors=FALSE, na.strings="NA", check.names=FALSE)
	data_matrix <- as.matrix(data_matrix.df[,cl.cols])
	rownames(data_matrix) <- data_matrix.df[,1]
	return(data_matrix)
}

main <- function(argv) {
	in_file <- argv[1]
	out_file <- argv[2]
	sink('/dev/null') 
	
	input_data <- read_matrix(in_file)
	input_data.df <- as.data.frame(input_data)
	input_data.lst <- as.list(input_data.df)
	input_data.unlst <- unlist(input_data.lst)
	input_data.nona <- input_data.unlst[!is.na(input_data.unlst)]
	
	pdf(out_file, bg="white")
	par(mfrow=c(1,1))
	hist(input_data.nona, col="lightblue", labels=TRUE, main="Histogram", xlab="")
	plot(density(input_data.nona), type="l", col="blue", main="Density")
	dev.off()
}
main(commandArgs(TRUE))	