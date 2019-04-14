#Reading and writing morphometric standard files
#Input data set:
#	-For write functions, morphometric data as shapes object
#Further Reading:
#Hammer, O. (1999-2012) "PAST: PAlaeontological STatistics Version 2.17 Reference...
#	Manual". (University of Oslo: Oslo).
#Rohlf, F. J. (2004) "NTSYSpc: Numerical Taxonomy and Multivariate Analysis...
#	System v. 2.1". (Exeter Software: Setauket).
#	http://www.exetersoftware.com/downloads/ntsysguide21.pdf
#Sheets, H. D. (2004) "IMP6-Integrated Morphometrics Package". (Canisius College:...
#	Buffalo)
#	http://www.canisius.edu/~sheets/morphsoft.html

#Author: Manuel Weinkauf  (Manuel.Weinkauf@unige.ch)
#Version: 1.2
#Date: 2 November 2017

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.#
#To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/.                   #
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#**************************************************************************************
#Setting working directory
#setwd("C:/R_TestData/Outline")

#########################################################################
# Write NTS file from shapes object                                     #
# Necessary input variables:                                            #
#    Input: Shapes object to be exported.                               #
#           *array*                                                     #
#    Matrix.type: Numeric code for matrix type.                         #
#                 1: rectangular data matrix                            #
#                 2: symmetric dissimilarity matrix                     #
#                 3: symmetric similarity matrix                        #
#                 4: diagonal matrix                                    #
#                 5: tree matrix for dissimilarity data                 #
#                 6: tree matrix for similarity data                    #
#                 7: graph matrix for dissimilarity data                #
#                 8: graph matrix for similarity data                   #
#                 *numeric (integer)*                                   #
#                 default=1                                             #
#    Row.labels: Row labels which may be printed in the file if desired.#
#                *vector*                                               #
#                default=1:dim(Input)[3]                                #
#    Col.labels: Column labels which may be printed in the file if...   #
#                desired.                                               #
#                *vector*                                               #
#                default=NULL                                           #
#    Missing.values: Do the data contain missing values.                #
#                    *logical*                                          #
#                    default=FALSE                                      #
#    Missing: Code for missing values, if present.                      #
#             *numeric*                                                 #
#             default=-999                                              #
#    Dimensions: Dimensions of dataset (2D or 3D).                      #
#                *numeric (integer)*,  either "2" or "3"                #
#                default=2                                              #
#    Output: Name of output file (including extension).                 #
#            *character*                                                #
# Output data: Morphometric data file in NTS format.                    #
# Input dataset: Morphometric object suitable for shapes-package.       #
#########################################################################

Write.NTS<-function (Input, Matrix.type=1, Row.labels=1:dim(Input)[3], Col.labels=NULL, Missing.values=FALSE, Missing=-999, Dimensions=2, Output) {
	#Check for consistency
	if (Matrix.type%in%1:8==FALSE) {stop("Matrix.type must be a number between 1 and 8!")}
	if (!is.null(Row.labels) && length(Row.labels)!=dim(Input)[3]) {stop("Row.labels must correspond in length to number of specimens in Input!")}
	if (!is.null(Col.labels) && length(Col.labels)!=(dim(Input)[1]*dim(Input)[2])) {stop("Col.labels must correspond in length to number of coordinates per specimen in Input!")}
	if (Dimensions%in%2:3==FALSE) {stop("Dimensions of dataset must be either 2D or 3D!")}
	
	#Prepare metadata
	{if (!is.null(Row.labels)) {Row.info<-paste(dim(Input)[3], "L", sep="")}
	else {Row.info<-as.character(dim(Input)[3])}}
	{if (!is.null(Col.labels)) {Col.info<-paste(dim(Input)[1]*dim(Input)[2], "L", sep="")}
	else {Col.info<-as.character(dim(Input)[1]*dim(Input)[2])}}
	{if (Missing.values==FALSE) {Miss.info<-as.character(0)}
	else {Miss.info<-paste(1, Missing, sep=" ")}}
	firstl<-paste(Matrix.type, Row.info, Col.info, Miss.info, paste("dim=", Dimensions, sep=""), sep=" ")
	##Create row and column label lines
	if (!is.null(Row.labels)) {Row.line<-paste(Row.labels, sep="", collapse=" ")}
	if (!is.null(Col.labels)) {Col.line<-paste(Col.labels, sep="", collapse=" ")}
	
	#Write output file
	cat(firstl, file=Output, sep="\n", append=FALSE)
	if (!is.null(Row.labels)) {cat(Row.line, file=Output, sep="\n", append=TRUE)}
	if (!is.null(Col.labels)) {cat(Col.line, file=Output, sep="\n", append=TRUE)}
	#Create data body
	for (j in 1:(dim(Input)[3])) {
		B<-NULL
		for (i in 1:(dim(Input)[1])) {
			B<-append(B, paste(Input[i,,j], sep="", collapse=" "))
		}
		B<-paste(B, sep="", collapse=" ")
		cat(B, file=Output, sep="\n", append=TRUE)
	}
}

#########################################################################
# Read NTS file	                                                        #
# Necessary packages: stringr                                           #
# Necessary input variables:                                            #
#    File: NTS file to be read.                                         #
#          *string*                                                     #
#    na.remove: Shall specimens that contain any NA values be removed?  #
#               *logical*                                               #
#               TRUE: remove those specimens                            #
#               FALSE: leave data as is                                 #
#               default=TRUE                                            #
# Output data: Morphometric object suitable for shapes-package.         #
# Input dataset: Morphometric data file in NTS format.                  #
#########################################################################

#Load packages
require(stringr)

Read.NTS<-function (File, na.remove=TRUE) {
	#Read the file as is into vector
	NTS<-scan(File, what="char", quote="", sep="\n", comment.char="\"")
	
	#Read and interpret header
	NTS.head<-NTS[1]
	NTS.meta<-list()
	NTS.meta$Dimension<-vector(mode="numeric", length=1)
	NTS.meta$Specimens<-vector(mode="numeric", length=2)
	NTS.meta$Landmarks<-vector(mode="numeric", length=2)
	names(NTS.meta$Specimens)<-names(NTS.meta$Landmarks)<-c("Number", "Labels")
	NTS.meta$Missing<-vector(mode="numeric", length=2)
	names(NTS.meta$Missing)<-c("Logical", "Value")
	##Read dimensions of outline data
	DIM<-unlist(strsplit(NTS.head, "="))
	{if (length(DIM)==1) {NTS.meta$Dimension[1]<-2}
	else {NTS.meta$Dimension[1]<-as.numeric(DIM[2])}}
	##Read dimensions of data file
	Meta<-unlist(strsplit(NTS.head, " "))
	NTS.meta$Specimens[1]<-as.numeric(str_extract(Meta[2],"[[:digit:]]+"))
	{if (is.na(str_extract(Meta[2],"[[:upper:]]"))) {NTS.meta$Specimens[2]<-0}
	else {NTS.meta$Specimens[2]<-1}}
	NTS.meta$Landmarks[1]<-as.numeric(str_extract(Meta[3],"[[:digit:]]+"))/NTS.meta$Dimension
	{if (is.na(str_extract(Meta[3],"[[:upper:]]"))) {NTS.meta$Landmarks[2]<-0}
	else {NTS.meta$Landmarks[2]<-1}}
	NTS.meta$Missing[1]<-as.numeric(Meta[4])
	if (NTS.meta$Missing[1]==1) {
		{if (is.na(as.numeric(Meta[5]))) {stop("NAs indicated, but no corresponding numerical value is given!")}
		else {NTS.meta$Missing[2]<-as.numeric(Meta[5])}}
	}
	
	#Interpret data file
	##Eliminate coding lines
	{if (NTS.meta$Specimens["Labels"]==1 && NTS.meta$Landmarks["Labels"]==1) {
		Specimens<-strsplit(NTS[2], " ")
		Landmarks<-strsplit(NTS[3], " ")
		NTS<-NTS[-(1:3)]
	}
	else if (NTS.meta$Specimens["Labels"]==1) {
		Specimens<-strsplit(NTS[2], " ")
		NTS<-NTS[-(1:2)]
	}
	else if (NTS.meta$Landmarks["Labels"]==1) {
		Landmarks<-strsplit(NTS[2], " ")
		NTS<-NTS[-(1:2)]
	}
	else {NTS<-NTS[-1]}}
	
	##Set up shapes object
	LMData<-array(NA, dim=c(NTS.meta$Landmarks["Number"], NTS.meta$Dimension, NTS.meta$Specimens["Number"]))
	{if (exists("Specimens") && exists ("Landmarks")) {dimnames(LMData)<-list(NULL, letters[24:(24+NTS.meta$Dimension-1)], Specimens[[1]])}
	else if (exists("Specimens")) {dimnames(LMData)<-list(NULL, letters[24:(24+NTS.meta$Dimension-1)], Specimens[[1]])}
	else if (exists ("Landmarks")) {dimnames(LMData)<-list(Landmarks[[1]], letters[24:(24+NTS.meta$Dimension-1)], NULL)}
	else {dimnames(LMData)<-list(NULL, letters[24:(24+NTS.meta$Dimension-1)], NULL)}
	}
	
	#Write data into shapes object
	{if (NTS.meta$Dimension==2) {
		Seq.x<-seq.int(from=1, to=NTS.meta$Landmarks["Number"]*NTS.meta$Dimension-1, by=2)
		Seq.y<-seq.int(from=2, to=NTS.meta$Landmarks["Number"]*NTS.meta$Dimension, by=2)
		for (i in 1:length(NTS)) {
			Temp<-strsplit(NTS[i], " ")
			LMData[,1,i]<-as.numeric(Temp[[1]][Seq.x])
			LMData[,2,i]<-as.numeric(Temp[[1]][Seq.y])
		}
	}
	else if (NTS.meta$Dimension==3) {
		Seq.x<-seq.int(from=1, to=NTS.meta$Landmarks["Number"]*NTS.meta$Dimension-2, by=3)
		Seq.y<-seq.int(from=2, to=NTS.meta$Landmarks["Number"]*NTS.meta$Dimension-1, by=3)
		Seq.z<-seq.int(from=3, to=NTS.meta$Landmarks["Number"]*NTS.meta$Dimension, by=3)
		for (i in 1:length(NTS)) {
			Temp<-strsplit(NTS[i], " ")
			LMData[,1,i]<-as.numeric(Temp[[1]][Seq.x])
			LMData[,2,i]<-as.numeric(Temp[[1]][Seq.y])
			LMData[,3,i]<-as.numeric(Temp[[1]][Seq.z])
		}
	}
	else {stop("Dimensionality must be either 2 or 3!")}
	}
	
	#Eliminate NA's
	if (NTS.meta$Missing[1]==1) {LMData[which(LMData==NTS.meta$Missing[2])]<-NA}
	ConTest<-vector(mode="logical", length=1)
	{if ((!anyNA(LMData) && NTS.meta$Missing[1]==0) | (anyNA(LMData) && NTS.meta$Missing[1]==1)) {ConTest<-TRUE}
	else {ConTest<-FALSE}}
	{if (na.remove==TRUE && ConTest==TRUE) {
		NA.entries<-vector(mode="logical", length=dim(LMData)[3])
		for (i in 1:(dim(LMData)[3])) {
			if (anyNA(LMData[,,i])) {NA.entries[i]<-TRUE}
		}
		if (any(NA.entries==TRUE)) {LMData<-LMData[,,-which(NA.entries==TRUE)]}
	}
	else if (na.remove==TRUE && ConTest==FALSE) {warning("Missing data detected, but not indicated in header! No data removed, please check data encoding.")}}
	
	return(LMData)
}

#########################################################################
# Write PAST morphology file from shapes object                         #
# Necessary input variables:                                            #
#    Input: Shapes object to be exported.                               #
#           *array*                                                     #
#    Row.labels: Row labels which may be printed in the	 file if...     #
#                desired.                                               #
#                *vector*                                               #
#                default=1:dim(Input)[3]                                #
#    Col.labels: Column labels which may be printed in the file if...   #
#                desired.                                               #
#                *vector*                                               #
#                default=paste(rep(letters[seq(from=24, to=25)],...     #
#                    dim(Input)[1]), 1:dim(Input)[1], sep="")           #
#    Output: Name of output file (including extension).                 #
#            *string*                                                   #
# Output data: Morphometric data file in PAST moprhology file format.   #
# Input dataset: Morphometric object suitable for shapes-package.       #
#########################################################################

Write.PAST<-function (Input, Row.labels=1:dim(Input)[3], Col.labels=paste(rep(letters[seq(from=24, to=25)], dim(Input)[1]), rep(1:dim(Input)[1], each=2), sep=""), Output) {
	#Check for consistency
	if (!is.null(Row.labels) && length(Row.labels)!=dim(Input)[3]) {stop("Row.labels must correspond in length to number of specimens in Input!")}
	if (!is.null(Col.labels) && length(Col.labels)!=(dim(Input)[1]*dim(Input)[2])) {stop("Col.labels must correspond in length to number of coordinates per specimen in Input!")}

	#Write output file
	##Create data body
	Res<-matrix(NA, dim(Input)[3], dim(Input)[1]*dim(Input)[2])
	for (j in 1:(dim(Input)[3])) {
		for (i in 1:(dim(Input)[1])) {
			Pos1<-i+(i-1)
			Pos2<-Pos1+1
			Res[j,Pos1]<-Input[i,1,j]
			Res[j,Pos2]<-Input[i,2,j]
		}
	}
	##Label and export data
	if (!is.null(Row.labels)) {rownames(Res)<-Row.labels}
	if (!is.null(Col.labels)) {colnames(Res)<-Col.labels}
	write.table(Res, Output, sep="\t", col.names=NA)
}

#########################################################################
# Read PAST morphology file                                             #
# Necessary input variables:                                            #
#    File: PAST file to be read.                                        #
#          *character*                                                  #
# Output data: Morphometric object suitable for shapes-package.         #
# Input dataset: Morphometric data file in PAST morphology file format. #
#########################################################################

Read.PAST<-function (File) {
	#Read and prepare file
	PAST<-read.table(File, header=FALSE, sep="\t", skip=1)
	PAST[,1]<-NULL
	colnames(PAST)<-paste(rep(letters[seq(from=24, to=25)], ncol(PAST)/2), rep(1:(ncol(PAST)/2), each=2), sep="")
	
	#Coerce data
	LMData<-array(NA, dim=c(ncol(PAST)/2, 2, nrow(PAST)), dimnames=list(NULL, c("x", "y"), NULL))
	Seq.x<-seq.int(from=1, to=ncol(PAST)-1, by=2)
	Seq.y<-seq.int(from=2, to=ncol(PAST), by=2)
	for (j in 1:nrow(PAST)) {
		for (i in 1:length(Seq.x)) {
			LMData[i,1,j]<-PAST[j,Seq.x[i]]
			LMData[i,2,j]<-PAST[j,Seq.y[i]]
		}
	}
	
	return(LMData)
}

#########################################################################
# Read IMP morphology file                                              #
# Necessary input variables:                                            #
#    File: IMP file to be read.                                         #
#          *character*                                                  #
# Output data: Morphometric object suitable for shapes-package.         #
# Input dataset: Morphometric data file in IMP morphology file format.  #
#########################################################################

Read.IMP<-function (File) {
	#Read and prepare file
	IMP<-read.table(File, header=FALSE, sep=" ", row.names=NULL)
	colnames(IMP)<-c(paste(rep(letters[seq(from=24, to=25)], (ncol(IMP)-1)/2), rep(1:((ncol(IMP)-1)/2), each=2), sep=""), "CS")
	
	#Coerce data
	Output<-list()
	Output$Centroid.Size<-IMP[,ncol(IMP)]
	IMP[,ncol(IMP)]<-NULL
	Output$LMData<-array(NA, dim=c(ncol(IMP)/2, 2, nrow(IMP)), dimnames=list(NULL, c("x", "y"), NULL))
	Seq.x<-seq.int(from=1, to=ncol(IMP)-1, by=2)
	Seq.y<-seq.int(from=2, to=ncol(IMP), by=2)
	for (j in 1:nrow(IMP)) {
		for (i in 1:length(Seq.x)) {
			Output$LMData[i,1,j]<-IMP[j,Seq.x[i]]
			Output$LMData[i,2,j]<-IMP[j,Seq.y[i]]
		}
	}
	
	return(Output)
}

#########################################################################
# Write TPS file from shapes object                                     #
# loosely based on Zelditch et al. (2012), write.tps()                  #
# Necessary input variables:                                            #
#    Input: Shapes object to be exported.                               #
#           *array*                                                     #
#    Centroids: List of centroid sizes corresponding to objects in...   #
#               Input.                                                  #
#               *vector*                                                #
#    Filenames: List of filenames corresponding to objects in Input.    #
#               *vector*                                                #
#    Scaling: Should the landmark data be rescaled back before export.  #
#             NOTE: This is only meaningful when the data have been...  #
#                   scaled e.g. on import using the Scaling=TRUE...     #
#                   option with Read.TPS. It is possible to have a...   #
#                   SCALE parameter present while at the same time...   #
#                   disabling the scaling of the landmarks upon...      #
#                   import, in which case Scaling=TRUE in this...       #
#                   function will produce wrong values!                 #
#                   Depending on how the data were treated before...    #
#                   there is still much room for errors, and the user...#
#                   is strictly advised to thoroughly think through...  #
#                   which option is to use here!                        #
#             *logical*                                                 #
#             TRUE: Landmarks will be scaled back to pixel values       #
#             FALSE: Landmarks will be left as is                       #
#             default=TRUE                                              #
#    Scale: List of scaling factors corresponding to objects in Input.  #
#           *vector*                                                    #
#    Output: Name of output file (including extension).                 #
#            *character*                                                #
# Output data: Morphometric data file in TPS format.                    #
# Input dataset: Morphometric object suitable for shapes-package.       #
#########################################################################

Write.TPS<-function (Input, Centroids=NULL, Filenames=NULL, Scaling=TRUE, Scale=NULL, Output) {
	#Check for existing file
	if (file.exists(Output)) {file.remove(Output)}
	
	#Check for consistency
	if (Scaling==TRUE & any(is.null(Scale), length(Scale)==0)) {stop("Scaling requested but no Scale provided!")}
	if (!is.null(Centroids)) {if (length(Centroids)!=(dim(Input)[3])) {stop("Centroids vector of wrong length!")}}
	if (!is.null(Filenames)) {if (length(Filenames)!=(dim(Input)[3])) {stop("Filenames vector of wrong length!")}}
	if (!is.null(Scale)) {if (length(Scale)!=(dim(Input)[3])) {stop("Scale vector of wrong length!")}}
	if (!is.null(Scale) & Scaling==FALSE) {warning("Scale present but no scaling requested! Is this correct?")}
	
	#Prepare metadata
	if (is.null(Centroids)) {Centroids<-rep(NA, dim(Input)[3])}
	if (is.null(Filenames)) {Filenames<-rep(NA, dim(Input)[3])}
	
	#Scale data
	if (Scaling==TRUE) {
		for (i in 1:(dim(Input)[3])) {
			Input[,,i]<-Input[,,i]/Scale[i]
		}
	}
	
	#Write output file
	firstl<-paste("LM=", dim(Input)[1], sep="")
	##Create file
	for (j in 1:(dim(Input)[3])) {
		cat(firstl, file=Output, sep="\n", append=TRUE)
		for (i in 1:(dim(Input)[1])) {
			B<-paste(Input[i,,j], sep="", collapse=" ")
			cat(B, file=Output, sep="\n", append=TRUE)
		}
		cat(paste("IMAGE=", Filenames[j], sep=""), file=Output, sep="\n", append=TRUE)
		cat(paste("ID=", dimnames(Input)[[3]][j], sep=""), file=Output, sep="\n", append=TRUE)
		if (!is.null(Scale)) {cat(paste("SCALE=", Scale[j], sep=""), file=Output, sep="\n", append=TRUE)}
	}
}

#########################################################################
# Read TPS file                                                         #
# based on Zelditch et al. (2012), read.tps2()                          #
# Necessary input variables:                                            #
#    File: TPS file to be read.                                         #
#          *string*                                                     #
#    Scale: Shall landmarks be scaled according to provided scale...    #
#           parameter?                                                  #
#           *logical*                                                   #
#           TRUE: Scale landmark coordinates                            #
#           FALSE: Leave landmark coordinates as they are               #
#           default=TRUE                                                #
#    na.remove: Shall specimens that contain any NA values be removed?  #
#               *logical*                                               #
#               TRUE: Remove those specimens                            #
#               FALSE: Leave data as is                                 #
#               default=TRUE                                            #
# Output data: Morphometric object suitable for shapes-package.         #
# Input dataset: Morphometric data file in TPS format.                  #
#########################################################################

Read.TPS<-function (File, Scale=TRUE, na.remove=TRUE) {
	#Read the file as is
	TPS<-readLines(File)
	
	#Gather meta data
	TPS.meta<-list()
	TPS.meta$Dimension<-vector(mode="numeric", length=1)
	TPS.meta$Specimens<-vector(mode="numeric", length=1)
	TPS.meta$Landmarks<-vector(mode="numeric", length=1)
	TPS.meta$Missing<-vector(mode="numeric", length=1)
	TPS.meta$ID<-vector(mode="character", length=0)
	TPS.meta$Image<-vector(mode="character", length=0)
	TPS.meta$Scale<-vector(mode="numeric", length=0)
	##Read dimensions of data file
	TPS.meta$Dimension<-length(strsplit(TPS[2], " ")[[1]])
	TPS.meta$Landmarks[1]<-as.numeric(strsplit(TPS[1], "=")[[1]][2])
	for(i in 1:length(TPS)) {
		if(substr(TPS[i], 1, 3)=="ID=") {TPS.meta$ID<-append(TPS.meta$ID, strsplit(TPS[i], "=")[[1]][2])}
		if(substr(TPS[i], 1, 6)=="IMAGE=") {TPS.meta$Image<-append(TPS.meta$Image, strsplit(TPS[i], "=")[[1]][2])}
		if(substr(TPS[i], 1, 6)=="SCALE=") {TPS.meta$Scale<-append(TPS.meta$Scale, as.numeric(strsplit(TPS[i], "=")[[1]][2]))}
		if(substr(TPS[i], 1, 3)=="LM=") {TPS.meta$Specimens<-TPS.meta$Specimens+1}
	}
	##Subset file to data-content only
	lenRecord<-length(TPS)/TPS.meta$Specimens[1]
	posLands<-rep(c(FALSE, rep(TRUE, TPS.meta$Landmarks[1]), rep(FALSE, lenRecord-TPS.meta$Landmarks[1]-1)), TPS.meta$Specimens)
	TPS.data<-TPS[posLands]
	
	#Set up output object
	Output<-list()
	Output$Filenames<-TPS.meta$Image
	Output$Scale<-TPS.meta$Scale
	Output$LMData<-array(NA, dim=c(TPS.meta$Landmarks, TPS.meta$Dimension, TPS.meta$Specimens), dimnames=list(NULL, letters[24:(24+TPS.meta$Dimension-1)], TPS.meta$ID))
	
	#Write data into shapes object
	LineCount<-1
	SpecCount<-1
	for (i in 1:(length(TPS.data))) {
		Output$LMData[LineCount,,SpecCount]<-as.numeric(strsplit(TPS.data[i], " ")[[1]])
		{if (LineCount<TPS.meta$Landmarks) {LineCount<-LineCount+1}
		else LineCount<-1}
		{if (LineCount==1) {SpecCount<-SpecCount+1}}
	}

	#Scale landmarks
	{if (Scale==TRUE & length(TPS.meta$Scale)==0) {warning("Scaling requested but no SCALE parameter provided, data were left as they are!")}
	else if (Scale==TRUE) {
		for (i in 1:(dim(Output$LMData)[3])) {
			Output$LMData[,,i]<-Output$LMData[,,i]*TPS.meta$Scale[i]
		}
	}
	}

	#Eliminate NA's
	if (na.remove==TRUE) {
		na.pos<-vector(mode="logical", length=dim(Output$LMData)[3])
		for (i in 1:(dim(Output$LMData)[3])) {
			{if (anyNA(Output$LMData[,,i])) {na.pos[i]<-TRUE}
			else {na.pos[i]<-FALSE}}
		}
		Output$Filenames<-Output$Filenames[which(na.pos==FALSE)]
		Output$Scale<-Output$Scale[which(na.pos==FALSE)]
		Output$LMData<-Output$LMData[,,which(na.pos==FALSE)]
	}

	return(Output)
}

#########################################################################
# Read Spiral file (as exported by function SpiralExtraction)           #
# Necessary input variables:                                            #
#    File: SPIRAL file to be read.                                      #
#          *string*                                                     #
# Output data: List object containing coordinates and polar coordi-...  #
#              nates of all specimens (one list entry per specimen).    #
# Input dataset: Morphometric data file in SPIRAL format.               #
#########################################################################

Read.Spiral<-function (File) {
	#Read file
	Dat<-read.table(File, header=FALSE, sep=",", row.names=1)
	
	#Prepare results object
	Output<-list()
	
	#Converse data into Output
	Pos<-1
	for (i in 1:nrow(Dat)) {
		{if (i%%4==1) {
			Temp<-Dat[i:(i+3),]
			Length<-max(which(!is.na(Temp[1,])))
			Output[[Pos]]<-matrix(NA, Length, 4)
			colnames(Output[[Pos]])<-c("x", "y", "t", "theta")
			Output[[Pos]][,"x"]<-as.vector(as.matrix(Temp[1,1:Length]))
			Output[[Pos]][,"y"]<-as.vector(as.matrix(Temp[2,1:Length]))
			Output[[Pos]][,"t"]<-as.vector(as.matrix(Temp[3,1:Length]))
			Output[[Pos]][,"theta"]<-as.vector(as.matrix(Temp[4,1:Length]))
			Pos<-Pos+1
		}}
	}
	
	return(Output)
}

#--------------------------------------------

#Examples
#setwd("C:/R_TestData/GeometricMorphometrics")
#Test1<-array(runif(100), dim=c(5, 2, 10), dimnames=list(NULL, c("x", "y"), 1:10))
#Test2<-matrix(NA, 10, 11)
#Test2[,1:10]<-runif(100)
#Test2[,11]<-sample(30, 10)
#write.table(Test2, "IMPExample_File.txt", sep=" ", row.names=FALSE, col.names=FALSE)
#Test3<-matrix(NA, 40, 20)
#rownames(Test3)<-paste(c("x", "y", "t", "theta"), rep(1:10, each=4), sep=".")
#for (i in 1:(nrow(Test3)-4)) {
#	{if (i%%4==1) {
#		VecLength<-sample(4:ncol(Test3), 1)
#		Test3[i,1:VecLength]<-sample(100, VecLength, replace=TRUE)
#	}
#	else if (i%%4==2) {Test3[i,1:VecLength]<-sample(100, VecLength, replace=TRUE)}
#	else if (i%%4==3) {Test3[i,1:VecLength]<-sort(runif(VecLength, min=0, max=1))} 
#	else {Test3[i,1:VecLength]<-runif(VecLength, min=0, max=1)}
#	}
#}
#Test3[(nrow(Test3)-3),1:ncol(Test3)]<-sample(100, ncol(Test3), replace=TRUE)
#Test3[(nrow(Test3)-2),1:ncol(Test3)]<-sample(100, ncol(Test3), replace=TRUE)
#Test3[(nrow(Test3)-1),1:ncol(Test3)]<-sort(runif(ncol(Test3), min=0, max=1))
#Test3[(nrow(Test3)),1:ncol(Test3)]<-runif(ncol(Test3), min=0, max=1)
#write.table(Test3, "SPIRALExample_File.spiral", sep=",", col.names=FALSE)

#NTS file format
#Write.NTS(Test1, Output="NTSExample_File.nts")
#Data<-Read.NTS("NTSExample_File.nts")

#File format for PAST 2.x and 3.x
#Write.PAST(Test1, Output="PASTExample_File.txt")
#Data<-Read.PAST("PASTExample_File.txt")

#IMP file format
#Data<-Read.IMP("IMPExample_File.txt")

#TPS file format
#Write.TPS(Test1, Filenames=paste("Image", 1:10, ".jpg", sep=""), Scaling=FALSE, Output="TPSExample_File.tps")
#Data<-Read.TPS("TPSExample_File.tps", Scale=FALSE, na.remove=FALSE)

#SPIRAL file format
#Data<-Read.Spiral("SPIRALExample_File.spiral")

#--------------------------------------------
#--------------------------------------------
#Version History
#1.0	Finished Program
#1.1	Added Write.PAST, Read.PAST, and Read.IMP functions
#	added colnames check in Write.NTS
#	added file cleanup for Write.TPS
#1.2	added Read.Spiral
#--------------------------------------------
#--------------------------------------------

