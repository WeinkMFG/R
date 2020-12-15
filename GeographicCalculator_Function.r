#Set of functions for geographical calculations
#Further reading: http://www.movable-type.co.uk/scripts/latlong.html	

#Author: Manuel Weinkauf (weinkaufmanuel@gmail.com)
#Version: 2.0
#Date: 15 December 2020

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.#
#To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/.                   #
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#**************************************************************************************
#Setting working directory
setwd("C:/R_TestData")

#########################################################################
# Function to translate geographical coordinates in another system      #
# Necessary input variables:                                            #
#   OutputType: Type of desired output coordinate format.               #
#               *character*                                             #
#               possible values: DMS=Degree--Minute--(Second)           #
#                                DMM=Degree--Minute.Minute-Fraction     #
#                                DD=Decimal Degrees                     #
#   Input: Name of input matrix, header assumed, missing values coded...#
#          as NA.                                                       #
#          *character*                                                  #
#   Output: Name of output matrix.                                      #
#           *character*                                                 #
#   FirstMeta: Does a first column with metadata exist?                 #
#              *logical*                                                #
#              TRUE: First column contains metadata.                    #
#              FALSE: First column contains coordinate type.            #
#              default=TRUE                                             #
#   NSEW: Is the input data coding with minus (south and west) and...   #
#         plus (north and east), or with "N, S, E, W"? In the latter... #
#         case, two additonal columns are expected in the input data,...#
#         which code the earth sector for the respective coordinate.... #
#         Note that if set to TRUE, NOT ALL data have to be coded...    #
#         that way; NA values then mark coordinates in plus/minus...    #
#         system.                                                       #
#         *logical*                                                     #
#         TRUE: At least some coordinates are coded as "N, S, E, W".    #
#         FALSE: ALL coordinates are coded with plus/minus.             #
#         default=FALSE                                                 #
#   Collapse: Should DMS and DMM be collapsed to text string format?    #
#             *logical*                                                 #
#             TRUE: Collapse e.g. 35 12 34 (three columns, numbers)...  #
#                   to 35°12'34'' (one column, text).                   #
#             FALSE: Leave data as pure numeral values, one column...   #
#                    per element.                                       #
#             default=FALSE                                             #
# Output data: Matrix file with converted coordinates.                  #
# Input dataset: Coordinate data in matrix. First column can contain... #
#                metadata (like core number), first or second column... #
#                must contain input format of coordinates in that row...#
#                (DMS, DMM, or DD; see OutputType)! All consecutive...  #
#                columns are expected to contain coordinates...         #
#                (latitude and longitude, in variable order) for...     #
#                that point.                                            #
#                Tow additional columns for decimal degrees data,...    #
#                four additional columns for degree--minute....         #
#                minute-fraction data (first degrees, second minute.... #
#                minute-fraction), six additional columns for...        #
#                degree--minute--second data (first degrees, second...  #
#                minutes, third seconds [optional]).                    #
#                If some coordinates are coded as "N, S, E, W", two...  #
#                additional columns are expected, one after each last...#
#                element of the respective coordinate, containing...    #
#                either of "N", "S", "E", or "W"; or NA for such...     #
#                coordinates which are coded by plus/minus.             #
#                All missing values must be coded as NA. See examples...#
#                for more information.                                  #
#########################################################################

#Input data set: A file containing coordinates to be translated. The whole code was reworked from v. 1.0, so...
#                data handling is completely different now. Coordinates of various types can be mixed in the same...
#                file and converted at once into ONE chosen output type. To enable that, an additional column was...
#                introduced, that codes the coordinate format. Furthermore, the coding of coordinates as "N", "S",...
#                "E", and "W" is now supported. The output, however, will in any case code northern latitudes and...
#                eastern longitudes positive, and southern latitudes and western longitudes negative.
#                Each logical element of a coordinate is contained in one column.
#                The file must naturally always contain as many columns as the longest data type necessitates; unused...
#                columns are then coded NA in other data entries. This means, the total number of columns must be somewhere...
#                between 3 (all decimal degrees without meta data) and 10 (at least one degree--minute--second coding...
#                with "N, S, E, W" coding for hemisphere, and metadata). Please refer to the following descriptions and...
#                the examples for more information about data formatting.
#                   -Metadata (e.g. core numbers) *optional: If present, they go in column 1 and FirstMeta is set to TRUE.
#                   -Input coding (DMS, DMM, DD): They go on a per row basis in column two, if metadata are present, else...
#                    in column one.
#                For the rest of the description, the term coordinate block will be used, to describe all data that belong...
#                to one part (either latitude or longitude) of one coordinate. The coordinate block consists of the...
#                obligatory C-Part, that contains all numerical values, and the optional H-Part, that contains the...
#                "N, S, E, W" coding for the hemisphere.
#                   -The H-Part *optional: The H-Part always goes in the last column belonging to that coordinate block, after...
#                    the respective C-Part. If present in at least one column set, it thus enlarges the number of columns...
#                    by two. If used in at least one row, NSEW must be set to TRUE. This does not necessitate to introduce...
#                    such a data coding in all rows, just using NA in that column for entries where the hemispheres are coded...
#                    by plus and minus is perfectly fine.
#                   -Degree--minute--second data (DMS): Such data necessitate the longest C-Part per se, which are three...
#                    columns per C-Part: Column 1--degree, Column 2--minutes, Column 3--seconds.
#                    If used together with "N, S, E, W" coding they produce coordinate blocks four columns wide,...
#                    thus the whole coordinate then requires eight columns.
#                    Example: -12°35'28'' or 12°35'28'' S
#                   -Degree--minute.minute-fraction data (DMM): Such data require two columns per C-Part, thus at least four...
#                    columns per coordinate: Column 1: degree, Column 2: minutes.
#                    Including an H-Part they need six columns per coordinate.
#                    Note that DMS data which were only measured to the minutes level (but without the precision of giving...
#                    the minutes as a non-integer value) also fall in that category.
#                    Example: -12° 35.48' or 12° 35.48' S, but also -12° 35'
#                   -Decimal degree data (DD): Those give coordinates as real value, thus only needing one column per C-Part:...
#                    Column 1--degree.
#                    If they come with an H-Part, which is rarely the case in those data, they need four columns per coordinate.
#                    Example: -12.53° or (seldom) 12.53° S

Trans<-function(OutputType, Input, Output=NULL, FirstMeta=TRUE, NSEW=FALSE, Collapse=FALSE) {
	#Read data
	Input<-read.table(Input, header=TRUE)
	{if (FirstMeta==TRUE) {Type<-2} else {Type<-1}}
	
	#Check data for consistency
	#Test if OutputType is meaningful
	if (!any((OutputType=="DMS"), (OutputType=="DMM"), (OutputType=="DD"))) {stop("OutputType must be either DMS, DMM, or DD!")}
	#Test if input types are meaningful
	T1<-Input[,Type]=="DMS"
	T2<-Input[,Type]=="DMM"
	T3<-Input[,Type]=="DD"
	T<-(T1 | T2 | T3)
	if (any(T==FALSE)==TRUE) {stop("Input type must be either DMS, DMM, or DD!")}
	
	#Preparing calculations
	#Initialising objects
	{if (OutputType=="DMS") {Res<-matrix(NA, nrow(Input), 7)}
		else if (OutputType=="DMM") {Res<-matrix(NA, nrow(Input), 5)}
		else {Res<-matrix(NA, nrow(Input), 3)}
	}
	Res[,1]<-(1:(dim(Input)[1]))
	
	#Checking data distribution
	{if (FirstMeta==TRUE) {COL<-ncol(Input)-2} else {COL<-ncol(Input)-1}}
	if (COL%%2!=0) {stop("Number of columns containing data must be even!")}
		##Determing range of first and second coordinate data
		{if (NSEW==FALSE) {
			Dx<-(1:(COL/2))
			Dy<-(((COL/2)+1):COL)
		}
		else {
			Dx<-(1:((COL/2)-1))
			Dy<-(((COL/2)+1):(COL-1))
		}
		}
		{if (FirstMeta==TRUE) {Dx<-Dx+2; Dy<-Dy+2} else {Dx<-Dx+1; Dy<-Dy+1}}
	
	#Test for meaningfulness of data
	T1<-which(Input[,Dx[1]] > 180 | Input[,Dx[1]] < -180 | Input[,Dy[1]] > 180 | Input[,Dy[1]] < -180)
	if (length(Dx)==2) {
		T1<-append(T1, (which((Input[,Type]=="DMS" & Input[,Dx[2]]>=60) | (Input[,Type]=="DMS" & Input[,Dy[2]]>=60))))
		T1<-append(T1, (which((Input[,Type]=="DMM" & Input[,Dx[2]]>=60) | (Input[,Type]=="DMM" & Input[,Dy[2]]>=60))))
	}
	if (length(Dx)==3) {
		T1<-append(T1, (which((Input[,Type]=="DMS" & Input[,Dx[3]]>=60) | (Input[,Type]=="DMS" & Input[,Dy[3]]>=60))))
	}
	if (length(T1>0)) {stop("There are degree values > 180, or minutes or seconds values > 60!")}
	
	#Converting coordinate data
	for (i in 1:nrow(Input)) {
		Temp<-vector()
		#Code for input in degree--minute--second...
		{if (Input[i,Type]=="DMS") {
			#...and output in degree--minute.minute-fraction
			{if (OutputType == "DMM") {
				Temp[1]<-Input[i,Dx[1]]
				Temp[2]<-Input[i,Dx[2]]
				{if (is.na(Input[i,Dx[3]])==TRUE) {Temp[3]<-0} else {Temp[3]<-round(((Input[i,Dx[3]])/60)*10)}}
				Res[i,2]<-Temp[1]
				Res[i,3]<-as.numeric(paste(Temp[2], ".", Temp[3], sep=""))
				Temp[1]<-Input[i,Dy[1]]
				Temp[2]<-Input[i,Dy[2]]
				{if (is.na(Input[i,Dy[3]])==TRUE) {Temp[3]<-0} else {Temp[3]<-round(((Input[i,Dy[3]])/60)*10)}}
				Res[i,4]<-Temp[1]
				Res[i,5]<-as.numeric(paste(Temp[2], ".", Temp[3], sep=""))
			}
			#...and output in decimal degrees
			else if (OutputType=="DD") {
				Temp[1]<-Input[i,Dx[1]]
				Temp[2]<-(Input[i,Dx[2]])/(60/100)
				{if (is.na(Input[i,Dx[3]])==TRUE) {Temp[3]<-0} else {Temp[3]<-(Input[i,Dx[3]])/((60*60)/100)}}
				Temp[4]<-round((Temp[2]+Temp[3])*100)
				{if (Input[i,Dx[2]]>=6) {Res[i,2]<-as.numeric(paste(Temp[1], ".", Temp[4], sep=""))} else {Res[i,2]<-as.numeric(paste(Temp[1], ".0", Temp[4], sep=""))}}
				Temp[1]<-Input[i,Dy[1]]
				Temp[2]<-(Input[i,Dy[2]])/(60/100)
				{if (is.na(Input[i,Dy[3]])==TRUE) {Temp[3]<-0} else {Temp[3]<-(Input[i,Dy[3]])/((60*60)/100)}}
				Temp[4]<-round((Temp[2]+Temp[3])*100)
				{if (Input[i,Dy[2]]>=6) {Res[i,3]<-as.numeric(paste(Temp[1], ".", Temp[4], sep=""))} else {Res[i,3]<-as.numeric(paste(Temp[1], ".0", Temp[4], sep=""))}}
			}
			else {
				Res[i,2]<-Input[i,Dx[1]]
				Res[i,3]<-Input[i,Dx[2]]
				{if (is.na(Input[i,Dx[3]])==TRUE) {Res[i,4]<-0} else {Res[i,4]<-Input[i,Dx[3]]}}
				Res[i,5]<-Input[i,Dy[1]]
				Res[i,6]<-Input[i,Dy[2]]
				{if (is.na(Input[i,Dy[3]])==TRUE) {Res[i,7]<-0} else {Res[i,7]<-Input[i,Dy[3]]}}
			}
			}
			}
		##Code for input in degree--minute.minute-fraction...
		else if (Input[i,Type]=="DMM") {
			#...and output in degree--minute--second
			{if (OutputType=="DMS") {
				Res[i,2]<-Input[i,Dx[1]]
				Res[i,3]<-trunc(Input[i,Dx[2]])
				Res[i,4]<-(Input[i,Dx[2]]-trunc(Input[i,Dx[2]]))*60
				Res[i,5]<-Input[i,Dy[1]]
				Res[i,6]<-trunc(Input[i,Dy[2]])
				Res[i,7]<-(Input[i,Dy[2]]-trunc(Input[i,Dy[2]]))*60
			}
			#...and output in decimal degrees
			else if (OutputType=="DD") {
				Temp[1]<-Input[i,Dx[1]]
				Temp[2]<-trunc(Input[i,Dx[2]])
				Temp[3]<-round((Input[i,Dx[2]]-trunc(Input[i,Dx[2]]))*60)
				Temp[4]<-(Temp[2])/(60/100)
				Temp[5]<-(Temp[3])/((60*60)/100)
				Temp[6]<-round((Temp[4]+Temp[5])*100)
				{if (Input[i,Dx[2]]>=6) {Res[i,2]<-as.numeric(paste(Temp[1], ".", Temp[6], sep=""))} else {Res[i,2]<-as.numeric(paste(Temp[1], ".0", Temp[6], sep=""))}}
				Temp[1]<-Input[i,Dy[1]]
				Temp[2]<-trunc(Input[i,Dy[2]])
				Temp[3]<-round((Input[i,Dy[2]]-trunc(Input[i,Dy[2]]))*60)
				Temp[4]<-(Temp[2])/(60/100)
				Temp[5]<-(Temp[3])/((60*60)/100)
				Temp[6]<-round((Temp[4]+Temp[5])*100)
				{if (Input[i,Dy[2]]>=6) {Res[i,3]<-as.numeric(paste(Temp[1], ".", Temp[6], sep=""))} else {Res[i,3]<-as.numeric(paste(Temp[1], ".0", Temp[6], sep=""))}}
			}
			else {
				Res[i,2]<-Input[i,Dx[1]]
				Res[i,3]<-Input[i,Dx[2]]
				Res[i,4]<-Input[i,Dy[1]]
				Res[i,5]<-Input[i,Dy[2]]
			}
			}
			}
		
		#Code for input in decimal degrees...
		else if (Input[i,Type]=="DD") {
			#...and output in degree--minute--second
			{if (OutputType=="DMS") {
				Res[i,2]<-trunc(Input[i,Dx[1]])
				Res[i,3]<-trunc(abs(Input[i,Dx[1]])*60)%%60
				Res[i,4]<-(abs(Input[i,Dx[1]])*(60*60))%%60
				Res[i,5]<-trunc(Input[i,Dy[1]])
				Res[i,6]<-trunc(abs(Input[i,Dy[1]])*60)%%60
				Res[i,7]<-(abs(Input[i,Dy[1]])*(60*60))%%60
			}
			#...and output in degree--minute.minute-fraction
			else if (OutputType=="DMM") {
				Temp[1]<-trunc(Input[i,Dx[1]])
				Temp[2]<-trunc(abs(Input[i,Dx[1]])*60)%%60
				Temp[3]<-(abs(Input[i,Dx[1]])*(60*60))%%60
				Temp[4]<-round(((Temp[3])/60)*10)
				Res[i,2]<-Temp[1]
				Res[i,3]<-as.numeric(paste(Temp[2], ".", Temp[4], sep=""))
				Temp[1]<-trunc(Input[i,Dy[1]])
				Temp[2]<-trunc(abs(Input[i,Dy[1]])*60)%%60
				Temp[3]<-(abs(Input[i,Dy[1]])*(60*60))%%60
				Temp[4]<-round(((Temp[3])/60)*10)
				Res[i,4]<-Temp[1]
				Res[i,5]<-as.numeric(paste(Temp[2], ".", Temp[4], sep=""))
			}
			else {
				Res[i,2]<-Input[i,Dx[1]]
				Res[i,3]<-Input[i,Dy[1]]
			}
			}
			}
		}
	}
	
	#If some data are given as NSEW, perform transformation of all those data
	if (NSEW==TRUE) {
		#Determine columns where NSEW info is given
		First<-COL/2
		Second<-COL
		{if (FirstMeta==TRUE) {First<-First+2; Second<-Second+2} else {First<-First+1; Second<-Second+1}}
		for (i in 1:nrow(Res)) {
			#For output type DMS
			{if (OutputType=="DMS") {
				{if (is.na(Input[i,First])==TRUE | Input[i,First]=="E" | Input[i,First]=="N") {}
				else if (Input[i,First]=="S" | Input[i,First]=="W") {Res[i,2]<-Res[i,2]*(-1)}
				else {stop("Hemisphere coding must be N, S, E, W, or NA!")}
				}
				{if (is.na(Input[i,Second])==TRUE | Input[i,Second]=="E" | Input[i,Second]=="N") {}
				else if (Input[i,Second]=="S" | Input[i,Second]=="W") {Res[i,5]<-Res[i,5]*(-1)}
				else {stop("Hemisphere coding must be N, S, E, W, or NA!")}
				}
			}
			#For output type DMM
			else if (OutputType=="DMM") {
				{if (is.na(Input[i,First])==TRUE | Input[i,First]=="E" | Input[i,First]=="N") {}
				else if (Input[i,First]=="S" | Input[i,First]=="W") {Res[i,2]<-Res[i,2]*(-1)}
				else {stop("Hemisphere coding must be N, S, E, W, or NA!")}
				}
				{if (is.na(Input[i,Second])==TRUE | Input[i,Second]=="E" | Input[i,Second]=="N") {}
				else if (Input[i,Second]=="S" | Input[i,Second]=="W") {Res[i,4]<-Res[i,4]*(-1)}
				else {stop("Hemisphere coding must be N, S, E, W, or NA!")}
				}
			}
			#For output type DD
			else if (OutputType=="DD") {
				{if (is.na(Input[i,First])==TRUE | Input[i,First]=="E" | Input[i,First]=="N") {}
				else if (Input[i,First]=="S" | Input[i,First]=="W") {Res[i,2]<-Res[i,2]*(-1)}
				else {stop("Hemisphere coding must be N, S, E, W, or NA!")}
				}
				{if (is.na(Input[i,Second])==TRUE | Input[i,Second]=="E" | Input[i,Second]=="N") {}
				else if (Input[i,Second]=="S" | Input[i,Second]=="W") {Res[i,3]<-Res[i,3]*(-1)}
				else {stop("Hemisphere coding must be N, S, E, W, or NA!")}
				}
			}
			}
		}
	}
	
	#Exporting results
	Results<-Res[,-1]
	{if (FirstMeta==TRUE) {rownames(Results)<-Input[,1]} else {rownames(Results)<-paste("Coord. Set", 1:nrow(Input), sep=" ")}}
	{if (Collapse==FALSE) {
		{if (OutputType=="DMS") {colnames(Results)<-c("Degrees", "Minutes", "Seconds", "Degrees", "Minutes", "Seconds")}
		else if (OutputType=="DMM") {colnames(Results)<-c("Degrees", "Minutes", "Degrees", "Minutes")}
		else {colnames(Results)<-c("Degrees", "Degrees")}
		}
	}
	else {
		ResT<-matrix(NA,dim(Results)[1],2)
		{if (OutputType=="DMS") {
			ResT[,1]<-paste(Results[,1], "° ", Results[,2], "' ", Results[,3], "''", sep="")
			ResT[,2]<-paste(Results[,4], "° ", Results[,5], "' ", Results[,6], "''", sep="")
			Results<-ResT
			colnames(Results)<-c("Coordinate 1", "Coordinate 2")
		}
		else if (OutputType=="DMM") {
			ResT[,1]<-paste(Results[,1], "° ", Results[,2], "'", sep="")
			ResT[,2]<-paste(Results[,3], "° ", Results[,4], "'", sep="")
			Results<-ResT
			colnames(Results)<-c("Coordinate 1", "Coordinate 2")
		}
		else {stop("Collapsing is only meaningfull for conversion into DMS and DMM systems!")}
		}
	}
	}
	write.table(Results, Output, sep="\t")
}

#########################################################################
# Function to calculate corner coordinates of rectangle around given... #
#    center point                                                       #
# Necessary packages: maptools, raster                                  #
# Necessary input variables:                                            #
#    Start: Latitude and Longitude of starting point in decimal...      #
#           degrees (-90:90°, -180:180°).                               #
#          *numeric (real)*, length=2                                   #
#    Dist: Distance from starting point into directions N, E, S, W...   #
#          (in this order) in km.                                       #
#          *numeric (real)*, length=4                                   #
# Output data: Map and coordinates of rectangle borders towards N, E,...#
#              S, and W around starting point.                          #
# Input data: Information about position and size of area.              #
#########################################################################

#Loading packages
require(maptools)
data(wrld_simpl)
require(raster)

GeoRect<-function (Start, Dist) {
	#Check data for consistency
	if (length(Start)!=2) {stop("Start coordinates (Lat, Lon) must be provided as vector of length 2")}
	if (length(Dist)!=4) {stop("Distances to N, E, S, and W must be provided as vector of length 4")}
	if (Start[1]<(-90) | Start[1]>90) {stop("Start latitude must lie between -90 and 90")}
	if (Start[2]<(-180) | Start[2]>180) {stop("Start longitude must lie between -180 and 180")}
	
	#Convert Start from degrees to radians
	Start<-Start*(pi/180)
	
	#Calculate parameters of rectangle
	Params<-matrix(NA, 2, 4)
	##Bearings
	Params[1,]<-c(0*(pi/180), 90*(pi/180), 180*(pi/180), 270*(pi/180))
	##Distances
	Params[2,]<-Dist

	#Calculate latitude and longitude of rectangle borders
	Borders<-vector(mode="numeric", length=4)
	names(Borders)<-c("Lat.N", "Lon.E", "Lat.S", "Lon.W")
	##To North
	Borders[1]<-asin((sin(Start[1])*cos(Params[2,1]/6371))+(cos(Start[1])*sin(Params[2,1]/6371)*cos(Params[1,1])))*(180/pi)

	##To East
	Aux<-asin((sin(Start[1])*cos(Params[2,2]/6371))+(cos(Start[1])*sin(Params[2,2]/6371)*cos(Params[1,2])))*(180/pi)
	Borders[2]<-(Start[2]+atan2(sin(Params[1,2])*sin(Params[2,2]/6371)*cos(Start[1]), cos(Params[2,2]/6371)-(sin(Start[1])*sin(Aux))))*(180/pi)

	##To South
	Borders[3]<-asin((sin(Start[1])*cos(Params[2,3]/6371))+(cos(Start[1])*sin(Params[2,3]/6371)*cos(Params[1,3])))*(180/pi)

	##To West
	Aux<-asin((sin(Start[1])*cos(Params[2,4]/6371))+(cos(Start[1])*sin(Params[2,4]/6371)*cos(Params[1,4])))*(180/pi)
	Borders[4]<-(Start[2]+atan2(sin(Params[1,4])*sin(Params[2,4]/6371)*cos(Start[1]), cos(Params[2,4]/6371)-(sin(Start[1])*sin(Aux))))*(180/pi)

	##Plot map
	plot(wrld_simpl, xlim=c(Borders[4]-20, Borders[2]+20), ylim=c(Borders[3]-20, Borders[1]+20), axes=TRUE, col="light yellow")
	box()
	Area<-extent(Borders[4], Borders[2], Borders[3], Borders[1])
	plot(Area, add=TRUE, col="red", lwd=2)
	points(Start[2]*(180/pi), Start[1]*(180/pi), pch=23, col="darkgreen", bg="yellow", cex=2)
	
	##Return coordinates
	return(Borders)
}

#--------------------------------------------

#Example
#Trans("DMS",Input="NESWTrue.txt",Output="Test1.txt",NSEW=TRUE)
#Trans("DMM",Input="NESWTrue.txt",Output="Test2.txt",NSEW=TRUE)
#Trans("DD",Input="NESWTrue.txt",Output="Test3.txt",NSEW=TRUE)
#Trans("DMS",Input="NESWTrue.txt",Output="Test1-1.txt",NSEW=TRUE,Collapse=TRUE)
#Trans("DMM",Input="NESWTrue.txt",Output="Test2-1.txt",NSEW=TRUE,Collapse=TRUE)
#Trans("DD",Input="NESWTrue.txt",Output="Test3-1.txt",NSEW=TRUE,Collapse=TRUE)#MUST PRODUCE ERROR


#Trans("DMS",Input="NESWFalse.txt",Output="Test4.txt")
#Trans("DMM",Input="NESWFalse.txt",Output="Test5.txt")
#Trans("DD",Input="NESWFalse.txt",Output="Test6.txt")

#GeoRect(c(33, -22), c(300, 150, 150, 300))

#--------------------------------------------
#--------------------------------------------
#Version History
#1.0	Finished Program
#2.0	Incorporated code from the old CoordinateTranslation_Function.r
#--------------------------------------------
#--------------------------------------------
