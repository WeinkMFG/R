#Convert data files into GMT compatible files?
#Input data set:
#	-World Ocean Atlas .csv file
#	-Ocean Data View spreadsheet .txt file

#Author: Manuel Weinkauf (Manuel.Weinkauf@unige.ch)
#Version: 1.2
#Date: 20 May 2019

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.#
#To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/.                   #
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#Creation of test dataset
#Lon<-sort(runif(20, min=10, max=50))
#Lat<-sort(runif(20, min=20, max=70))
#Raw.Path<-cbind(Lon, Lat)

#**************************************************************************************
#Setting working directory
#setwd("C:/R_TestData/GMT")

#########################################################################
# Function to summarize information from WOA .csv file                  #
# Necessary input variables:                                            #
#   File: Name of the input file (.csv)                                 #
#         *character*                                                   #
# Output data: Summary of the depth present in the dataset. This is...  #
#              deemed helpful for the determination of Depth.interval...#
#              in function WOA.to.netCDF.                               #
# Input dataset: Comma separated file, as downloaded from World Ocean...#
#                Atlas.                                                 #
#########################################################################

WOA.Info<-function (File) {
	#Read file
	Depths<-scan(File, sep=",", skip=1, nlines=1, what="raw", quiet=TRUE)
	Depths<-c(as.numeric(strsplit(Depths[3], split=":", fixed=TRUE)[[1]][2]), as.numeric(Depths[4:length(Depths)]))
	Input<-read.csv(File, header=FALSE, sep=",", comment.char="#", col.names=c("Lat", "Lon", as.character(Depths)), fill=TRUE)
	
	#Summarize data
	print(paste("Number of data entries:", nrow(Input), sep=" "))
	Depth.Info<-list()
	Depth.Info[[1]]<-Depths
	names(Depth.Info)<-"Available depths:"
	print(Depth.Info)
}

#########################################################################
# Function to read WOA .csv file                                        #
# Necessary input variables:                                            #
#   File: Name of the input file (.csv)                                 #
#         *character*                                                   #
# Output data: Matrix with the data in the WOA file.                    #
# Input dataset: Comma separated file, as downloaded from World Ocean...#
#                Atlas.                                                 #
#########################################################################

WOA.Read<-function (File) {
	#Read file
	Depths<-scan(File, sep=",", skip=1, nlines=1, what="raw", quiet=TRUE)
	Depths<-c(as.numeric(strsplit(Depths[3], split=":", fixed=TRUE)[[1]][2]), as.numeric(Depths[4:length(Depths)]))
	Input<-read.csv(File, header=FALSE, sep=",", comment.char="#", col.names=c("Lat", "Lon", paste(as.character(Depths), "m", sep=".")), fill=TRUE)
	
	#Return results
	return(Input)
}

#########################################################################
# Function to convert WOA .csv file into 1D netCDF file                 #
# Necessary packages: raster, ncdf4                                     #
# Necessary input variables:                                            #
#   File: Name of the input file (.csv)                                 #
#         *character*                                                   #
#   Output.File: Name of the output file.                               #
#                *character*                                            #
#   Depth.interval: Range of the depth interval for which data should...#
#                   be extracted. Note that data are averaged over...   #
#                   the given interval.                                 #
#                   *numeric (integer)* of length 2                     #
#   Resolution: Resolution of the WOA file version in degrees...        #
#               (normally 5, 1, or 1/4 degree).                         #
#               *numeric (real)*                                        #
#               default=1                                               #
#   Average: Method to average data across layers.                      #
#            one of min, max, mean                                      #
#            default=mean                                               #
# Output data: 1D netCDF file compatible with GMT.                      #
# Input dataset: Comma separated file, as downloaded from World Ocean...#
#                Atlas.                                                 #
#########################################################################

#Loading packages
require("raster")
require("ncdf4")

WOA.to.GMT<-function (File, Output.File, Depth.interval, Resolution=1, Average=mean) {
	#Read file
	Depths<-scan(File, sep=",", skip=1, nlines=1, what="raw", quiet=TRUE)
	Depths<-c(as.numeric(strsplit(Depths[3], split=":", fixed=TRUE)[[1]][2]), as.numeric(Depths[4:length(Depths)]))
	Input<-read.csv(File, header=FALSE, sep=",", comment.char="#", col.names=c("Lat", "Lon", as.character(Depths)), fill=TRUE)
	
	#Test data consistency
	if (length(Depth.interval)!=2) {stop("Depth.interval must have length two!")}
	if (!all(Depth.interval%in%Depths)) {stop("Depth.interval may only contain available depths. Consult WOA.Info for help.")}
	if (Depth.interval[1] > Depth.interval[2]) {
		Temp<-Depth.interval[1]
		Depth.interval[1]<-Depth.interval[2]
		Depth.interval[2]<-Temp
		warning("Lower Depth.interval was larger than upper Depth.interval. Values have been swapped!")
	}
	
	#Reduce input data to intended depth interval
	Coord<-Input[,1:2]
	Input[,1:2]<-NULL
	Input<-Input[,which(Depths==Depth.interval[1]):which(Depths==Depth.interval[2])]
	
	#Calculate mean over depths
	Means<-apply(Input, 1, mean, na.rm=TRUE)
	Output<-as.matrix(cbind(Coord, Means))
	colnames(Output)<-c("Lat", "Lon", "Mean")
	Output[which(is.nan(Output[,"Mean"])),"Mean"]<-NA
	
	#Export as netCDF file
	r<-raster(resolution=c(Resolution, Resolution))
	Output.r<-rasterize(Output[,2:1], r, field=Output[,3], fun=Average, background=NA)
	writeRaster(Output.r, Output.File, format="CDF")
}

#########################################################################
# Function to summarize information ODV spreadsheet .txt file           #
# Necessary input variables:                                            #
#   File: Name of the input file (.txt)                                 #
#         *character*                                                   #
# Output data: Summary of the variables present in the dataset. This... #
#              is deemed helpful for the determination of Lon, Lat,...  #
#              and Variable in function ODV.to.GMT.                     #
# Input dataset: Tabstop delimited .txt file, as exported by ODV via... #
#                "Export as ODV spreadsheet".                           #
#########################################################################

ODV.Info<-function (File) {
	#Read file
	Input<-read.table(File, header=TRUE, sep="\t", comment.char="/")
	
	#Summarize data
	print(paste("Number of data entries:", nrow(Input), sep=" "))
	Var.Info<-list()
	Var.Info[[1]]<-colnames(Input)
	names(Var.Info)<-"Available variables:"
	print(Var.Info)
}

#########################################################################
# Function to read ODV spreadsheet .txt file                            #
# Necessary input variables:                                            #
#   File: Name of the input file (.txt)                                 #
#         *character*                                                   #
# Output data: Matrix with the data in the ODV file.                    #
# Input dataset: Tabstop delimited .txt file, as exported by ODV via... #
#                "Export as ODV spreadsheet".                           #
#########################################################################

ODV.Read<-function (File) {
	#Read file
	Input<-read.table(File, header=TRUE, sep="\t", comment.char="/")
	
	#Return data
	return(Input)
}

#########################################################################
# Function to convert ODV spreadsheet .txt file into xyz txt files...   #
#    suitable for interpolation in GMT (surface, triangulate)           #
# Necessary input variables:                                            #
#   File: Name of the input file (.csv)                                 #
#         *character*                                                   #
#   Output.File: Name of the output file.                               #
#                *character*                                            #
#   Variable: Index of the variable of interest.                        #
#             *numeric (integer)*                                       #
#   CoordType: Are coordinates given as -180/180/-90/90 degree...       #
#              system or as 0/360/-90/90 degree system?                 #
#              *character*, one of "180" or "360"                       #
#   Lon: Index of the column containing longitude data.                 #
#        *numeric (integer)*                                            #
#        default=5                                                      #
#   Lat: Index of the column containing latitude data.                  #
#        *numeric (integer)*                                            #
#        default=6                                                      #
# Output data: Tabstop delimited xyz file for further use in GMT.       #
# Input dataset: Tabstop delimited .txt file, as exported by ODV via... #
#                "Export as ODV spreadsheet".                           #
#########################################################################

ODV.to.GMT<-function (File, Output.File, Variable, CoordType, Lon=5, Lat=6) {
	#Read file
	Input<-read.table(File, header=TRUE, sep="\t", comment.char="/")
	
	#Test data consistency
	if (CoordType!="180" & CoordType!="360") {stop("CoordType must be either '180' or '360'!")}
	
	#Prepare export file
	Output<-matrix(NA, nrow(Input), 3)
	colnames(Output)<-c("#Lon", "Lat", "Value")
	
	#Convert longitude data
	{if (CoordType=="360") {
		Lon.Dat<-Input[,Lon]
		Lon.Cor<-vector(mode="numeric", length=length(Lon.Dat))
		East<-which(Lon.Dat>=0 & Lon.Dat<=180)
		West<-which(Lon.Dat>180)
		Other<-which(is.na(Lon.Dat))
		Lon.Cor[East]<-Lon.Dat[East]
		Diff<-Lon.Dat[West]-180
		Lon.Cor[West]<--180+Diff
		Lon.Cor[Other]<-NA
		Output[,"#Lon"]<-Lon.Cor
	}
	else {Output[,"#Lon"]<-Input[,Lon]}}
	Output[,"Lat"]<-Input[,Lat]
	Output[,"Value"]<-Input[,Variable]
	
	#Export data
	write.table(Output, Output.File, sep="\t", quote=FALSE, row.names=FALSE)
}

#########################################################################
# Function to smooth digitized paths for plotting in GMT                #
# Necessary input variables:                                            #
#   Coord: Dataset containing raw coordinates.                          #
#      *matrix* with two columns                                        #
#   iter: number of smoothing iterations.                               #
#         *numeric (integer)*                                           #
# Output data: Smoothed coordinates.                                    #
# Input dataset: Raw coordinates to smooth.                             #
#########################################################################

Path.smoothout<-function(Coord, iter){
	p<-dim(Coord)[1]
	a<-0
	while(a<=iter){
		a<-a+1
		Ms<-rbind(Coord[1,], Coord[-p,])
		Mi<-rbind(Coord[-1,], Coord[p,])
		Coord<-Coord/2+Ms/4+Mi/4
	}
	Coord
}


#--------------------------------------------

#Example
##World Ocean Atlas
#WOA.Info("woa13_decav_t13mn01v2.csv")
#WOA.to.GMT("woa13_decav_t13mn01v2.csv", "Test.nc", c(0, 50))

##Ocean Data View
#ODV.Info("data_from_MARGO_ModernPF_Synthesis.txt")
#ODV.to.GMT("data_from_MARGO_ModernPF_Synthesis.txt", "MARGO.txt", Variable=15, CoordType="360")

##Path smoothing
#Smooth.Path<-Path.smoothout(Raw.Path, 3)
#plot(Raw.Path, type="l", col="grey", lwd=2)
#lines(Smooth.Path, col="blue", lwd=2)

#--------------------------------------------
#--------------------------------------------
#Version History
#1.0	Finished Program
#1.1	Added function Path.smoothout
#1.2	Added WOA.Read and ODV.Read functions
#--------------------------------------------
#--------------------------------------------
