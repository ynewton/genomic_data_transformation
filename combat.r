#!/usr/bin/Rscript

#Yulia Newton

#usage, options and doc goes here
argspec <- c("combat.r - Takes any matrix and returns a matrix with major batch effect remove. Requires a single header line and a single cloumn of annotation in the input matrix. Output will be a matrix of the same dimension.
Usage:
    combat.r input.tab source.tab output.tab
Example:
	Rscript combat.r input_table.tab source.tab output.tab
	input matrix (annotated by row and column names)
	annotation matrix (first column is the batch effect to remove, other columns are covariates for true biological structure
	output matrix
	")

read_dataframe <- function(in_file){
    df  <- read.delim(in_file,header=TRUE,row.names=1,stringsAsFactors=F,check.names=F)                        
    return (as.data.frame(df))
}
read_matrix <- function(in_file){
	header <- strsplit(readLines(con=in_file, n=1), "\t")[[1]]
	cl.cols<- 1:length(header) > 1
	data_matrix.df <- read.delim(in_file, header=TRUE, row.names=NULL, stringsAsFactors=FALSE, na.strings="NA", check.names=FALSE)
	data_matrix <- as.matrix(data_matrix.df[,cl.cols])
	rownames(data_matrix) <- data_matrix.df[,1]
	return(data_matrix)
}

write_dataframe <- function (data, output_file) {
    write.table(data,file=output_file,sep="\t")
}

write_matrix <- function(data_matrix,output_file){
	header <- append(c(""), colnames(data_matrix))
	write.table(t(header), output_file, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
	write.table(data_matrix, output_file, quote=FALSE, sep="\t", row.names=TRUE, col.names=FALSE)
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

read_col_names <- function(in_file){
	return(read.table(in_file, header=FALSE, sep="", as.is = TRUE)$V1)
}

remove_batch_effect <- function(pre.COMBAT, batches){
#cool info: http://genomicsclass.github.io/book/pages/expressing_design_formula.html
#https://support.bioconductor.org/p/62058/
	if(ncol(batches) == 1){
    	post.COMBAT <- sva::ComBat(dat=pre.COMBAT,batch=batches[,1],mod=NULL)
	}else{
		#mod <- model.matrix(~1,data=batches)
		#mod <- model.matrix(~as.factor(batches[,2]) + ~as.factor(batches[,3]) + 0,data=batches)
		
		ncovars <- ncol(batches) - 1
		for(c in 1:ncovars){
			mod_c <- model.matrix(~as.factor(batches[,(c+1)]),data=batches)
			if(c == 1){
				mod <- mod_c
			} else{
				mod <- cbind(mod,mod_c[,-1])
			}
		}
		
    	post.COMBAT <- sva::ComBat(dat=pre.COMBAT,batch=batches[,1],mod=mod,numCovs=NULL)
    }
    return (post.COMBAT)
}

main <- function(argv) {  
    #library(sva)
    source("http://bioconductor.org/biocLite.R")
    biocLite("sva")
    library("sva")
    biocLite("mgcv")
    library("mgcv")
    
	if(length(argv) != 3){
		if(argv==c('--help')){ 
			write(argspec, stderr());
			q();
		}
	}
		
	if(!(length(argv) == 3)){
		write("ERROR: invalid number of arguments is specified", stderr());
		q();
	}
	
	#store command line arguments in variables:
	input_file <- argv[1]
	source_file <- argv[2]
	output_file <- argv[3]
	
	
	#read the input file(s):
	data_matrix <- read_dataframe(input_file)
	gene.sd <- unlist(apply(data_matrix,1,sd))
	rows.for.removal <- names(which(gene.sd == 0.0))
	data_matrix <- data_matrix[!(rownames(data_matrix) %in% rows.for.removal),]
	
	#if (min(data_matrix) < 0) {
	#	write(c("matrix cannot have negative numbers, min value", min(data_matrix)), stderr());
	#	q();
	#}
	
	batches <- read_dataframe(source_file)
	common_columns <- intersect(rownames(batches), colnames(data_matrix))
	data_matrix <- data_matrix[,colnames(data_matrix) %in% common_columns]
	batches.filtered <- as.data.frame(batches[rownames(batches) %in% common_columns,])
	#rownames(batches.filtered) <- rownames(batches)
	data_matrix <- data_matrix[,order(colnames(data_matrix))]
	batches <- as.data.frame(batches.filtered[order(rownames(batches.filtered)),])
	
	result <- remove_batch_effect((data_matrix+1), batches)
	
	#write.table(post.COMBAT,file="output.tab",sep="\t")
    #write_dataframe(result, output_file)
    write_matrix_to_file(result, output_file)
}

#batches <- read.delim("/data/import/UCSCprostate/data/wrangled-current/joint_source.tab",header=TRUE,row.names=1,stringsAsFactors=F,check.names=F)                      
#need a model - check the sva vignette
#since we are only using one adjustment variable (the batch one)
#and the batch variable is added in later on in COMBAT, the model only has  an intercept
main(commandArgs(TRUE))
