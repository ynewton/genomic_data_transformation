#!/usr/bin/Rscript

#Yulia Newton, last updated 20150102 v.1

#usage, options and doc goes here
argspec <- c("limma.r - computes sam differential
Usage:
    limma.r input.tab samples.class.1 samples.class.0 output.tab
Example:
	Rscript limma.r test_matrix.tab samples.class.1.tab samples.class.0.tab output.tab
Options:
	input matrix (annotated by row and column names)
	list of samples of class 1
	list of samples of class 0
	output file prefix (.tab and .pdf will be added to the output files)")

read_matrix <- function(in_file){
	header <- strsplit(readLines(con=in_file, n=1), "\t")[[1]]
	cl.cols<- 1:length(header) > 1
	data_matrix.df <- read.delim(in_file, header=TRUE, row.names=NULL, stringsAsFactors=FALSE, na.strings="NA", check.names=FALSE)
	data_matrix <- as.matrix(data_matrix.df[,cl.cols])
	rownames(data_matrix) <- data_matrix.df[,1]
	return(data_matrix)
}

write_matrix <- function(data_matrix){
	header <- append(c("Sigs"), colnames(data_matrix))
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

read_list <- function(in_file){
	return(as.matrix(read.table(in_file, header=FALSE, sep="", as.is = TRUE))[, 1])
	#return(as.matrix(read.table(in_file, header=FALSE, sep="", as.is = TRUE)))
}

main <- function(argv) {
	if(length(argv) == 1){
		if(argv==c('--help')){ 
			write(argspec, stderr());
			q();
		}
	}

	#source("http://bioconductor.org/biocLite.R")
	#biocLite(c("limma","edgeR"))
	library(edgeR)
	library(limma)	
		
	if(!(length(argv) == 4)){
		write("ERROR: invalid number of arguments is specified", stderr());
		q();
	}
	
	#store command line arguments in variables:
	input_file <- argv[1]
	samples.1.file <- argv[2]
	samples.0.file <- argv[3]
	output_file <- argv[4]
		
	#read the input file(s):
	data_matrix <- read_matrix(input_file)
	samples.1 <- read_list(samples.1.file)
	samples.1 <- samples.1[samples.1 %in% colnames(data_matrix)]
	samples.0 <- read_list(samples.0.file)
	samples.0 <- samples.0[samples.0 %in% colnames(data_matrix)]
	
	limma.samples <- append(samples.1,samples.0)
	data_matrix.limma <- data_matrix[,colnames(data_matrix) %in% limma.samples]
	data_matrix.limma <- data_matrix.limma[,limma.samples]
	
	#followed the steps here:
	#http://master.bioconductor.org/help/course-materials/2005/BioC2005/labs/lab01/estrogen/
	
	#design matrix:
	#samples.1.class <- rep(1,length(samples.1))
	#samples.0.class <- rep(2,length(samples.0))
	#samples.class <- append(samples.1.class,samples.0.class)
	#names(samples.class) <- limma.samples
	
	#samples.class <- cbind(limma.samples, samples.class)
	#design.matrix <- model.matrix(~samples.class[,2])
		
	samples.1.class <- rep(1,length(samples.1))
	samples.1.not.class <- rep(0,length(samples.0))
	samples.1.class_col <- append(samples.1.class,samples.1.not.class)
	samples.0.class <- rep(1,length(samples.0))
	samples.0.not.class <- rep(0,length(samples.1))
	samples.0.class_col <- append(samples.0.not.class,samples.0.class)
	samples.class.mtrx <- cbind(samples.1.class_col,samples.0.class_col)
	colnames(samples.class.mtrx) <- c("class1","class0")
	rownames(samples.class.mtrx) <- limma.samples
	
	fit <- lmFit(data_matrix.limma, samples.class.mtrx)
	
	#contrast matrix:
	contrast.matrix = matrix(1,ncol=1,nrow=2)
	contrast.matrix[2,1] <- -1
	
	fit2  <- contrasts.fit(fit, contrast.matrix)
	fit2  <- eBayes(fit2)
	
	tt = topTable(fit2, coef=1, number=length(rownames(data_matrix.limma)))
	
	write_matrix_to_file(tt, output_file)
}

main(commandArgs(TRUE))
