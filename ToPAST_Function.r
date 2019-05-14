#Read and write PAST compatible files?

#Author: Manuel Weinkauf (Manuel.Weinkauf@unige.ch)
#Version: 1.0
#Date: 14 May 2019

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.#
#To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/.                   #
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#Creation of test dataset
#Dat<-matrix(runif(30*30), 30, 30)

#**************************************************************************************
#Setting working directory
#setwd("C:/R_TestData/PAST")

#########################################################################
# Read PAST data files                                                  #
# Necessary input variables:                                            #
#    File: PAST file to be read.                                        #
#          *character*                                                  #
#    version: Which version of PAST was used to export the file...      #
#             Either 2 or 3.                                            #
#             *numeric (integer)*                                       #
#             default=2                                                 #
# Output data: Matrix object.                                           #
# Input dataset: Data file from PAST (.dat-format).                     #
#########################################################################

Read.PAST.Data<-function (File, version=2) {
	#Test input consistency
	if (version!=2 & version!=3) {stop("Version must be either '2' or '3'!")}

	#Read and prepare file
	{if (version==2) {
		PAST<-read.table(File, header=FALSE, sep="\t", colClasses="character")
		if (all(PAST[,ncol(PAST)]=="")) {PAST<-PAST[,1:(ncol(PAST)-1)]}
		colnames(PAST)<-PAST[1,]; PAST<-PAST[2:nrow(PAST),]
		rownames(PAST)<-PAST[,1]; PAST[,1]<-NULL
		PAST<-apply(PAST, c(1, 2), as.numeric)
	}
	else {
		PAST<-read.table(File, header=FALSE, sep="\t", colClasses="character")
		Metadata<-list()
		Metadata$Colour<-PAST[3:nrow(PAST),1]
		Metadata$Symbol<-PAST[3:nrow(PAST),2]
		colnames(PAST)<-PAST[2,]; PAST<-PAST[3:nrow(PAST),]
		rownames(PAST)<-PAST[,3]; PAST<-PAST[,4:ncol(PAST)]
		PAST<-apply(PAST, c(1, 2), as.numeric)
		if (all(is.na(PAST[,ncol(PAST)]))) {PAST<-PAST[,1:(ncol(PAST)-1)]}
	}
	}
	
	{if (version==2) {return(PAST)}
	else {return(list("Data"=PAST, "Metadata"=Metadata))}}
}

#########################################################################
# Write PAST data files                                                 #
# Necessary input variables:                                            #
#    Input: Object to be exported.                                      #
#           *matrix*                                                    #
#    Row.labels: Row labels which may be printed in the	 file if...     #
#                desired.                                               #
#                *vector*                                               #
#                default=1:nrow(Input)                                  #
#    Col.labels: Column labels which may be printed in the file if...   #
#                desired.                                               #
#                *vector*                                               #
#                default=NULL                                           #
#    Output: Name of output file (including extension).                 #
#            *string*                                                   #
#    version: Which version of PAST the export file should be...        #
#             compatible with? Either 2 or 3.                           #
#             *numeric (integer)*                                       #
#             default=2                                                 #
#    Col.class: Classes of the data columns. Either of "-"...           #
#               (numerical), "Group", "Ordinal", "Nominal", or "Binary".#
#               *character*                                             #
#               default: "-" for all columns                            #
#    Col: Colour information for PAST v. 3.x (ignored if version==2).   #
#         *character*                                                   #
#         default="Black"                                               #
#    Sym: Symbol information for PAST v. 3.x (ignored if version==2).   #
#         *character*                                                   #
#         default="Dot"                                                 #
# Output data: Data file in PAST's .dat-format.                         #
# Input dataset: Matrix object.                                         #
#########################################################################

Write.PAST.Data<-function (Input, Row.labels=as.character(1:nrow(Input)), Col.labels=NULL, Output, version=2, Col.class=rep("-", ncol(Input)), Col=rep("Black", nrow(Input)), Sym=rep("Dot", nrow(Input))) {
	#Check for consistency
	if (version!=2 & version!=3) {stop("Version must be either '2' or '3'!")}
	if (!is.null(Row.labels) && length(Row.labels)!=nrow(Input)) {stop("Row.labels must correspond in length to number of specimens in Input!")}
	if (!is.null(Col.labels) && length(Col.labels)!=(ncol(Input)[1])) {stop("Col.labels must correspond in length to number of coordinates per specimen in Input!")}
	if (length(Col.class)!=ncol(Input)) {stop("One value of Col.class per parameter needed.")}
	if (!all(Col.class %in% c("-", "Group", "Ordinal", "Nominal", "Binary"))) {stop("Col.class must be either of '-', 'Group', 'Ordinal', 'Nominal', or 'Binary'!")}
	if (length(Col)!=nrow(Input)) {stop("One value of Col per specimen needed.")}
	if (length(Sym)!=nrow(Input)) {stop("One value of Sym per specimen needed.")}
	
	#Label data
	{if (!is.null(Row.labels)) {rownames(Input)<-Row.labels} else {rownames(Input)<-as.character(1:nrow(Input))}}
	ml<-vector(mode="list", length=3)
	ml[[1]]<-LETTERS
	for (i in 2:3) {
		ml[[i]]<-c(sapply(ml[[i-1L]], function(y) paste0(y, LETTERS)))
	}
	ml<-unlist(ml)
	{if (!is.null(Col.labels)) {colnames(Input)<-Col.labels} else {colnames(Input)<-ml[1:ncol(Input)]}}
	
	#Export data
	{if (version==2) {
		write.table(Input, Output, quote=FALSE, sep="\t", col.names=c(paste(".", colnames(Input)[1], sep="\t"), colnames(Input)[2:ncol(Input)]))
	}
	else {
		Input<-cbind(Col, Sym, rownames(Input), Input)
		Input<-rbind(c("", "", "", colnames(Input)[4:ncol(Input)]), Input)
		Input<-rbind(c(":", "", "", Col.class), Input)
		write.table(Input, Output, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
	}
	}
}

#--------------------------------------------

#Examples
#setwd("C:/R_TestData/PAST")
#Write.PAST.Data(Dat, Output="Past2_Export.dat")
#Write.PAST.Data(Dat, Output="Past3_Export.dat", version=3)
#Read.PAST.Data("Past2_Export.dat")
#Read.PAST.Data("Past3_Export.dat", version=3)

#--------------------------------------------
#--------------------------------------------
#Version History
#1.0	Finished Program
#--------------------------------------------
#--------------------------------------------
