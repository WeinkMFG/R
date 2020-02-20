#Extracting data points coordinates from a diagram (image)
#Input data set: A digitised image of a diagram (or any other image type from which some kind of coordinates...
#                should be retrieved)

#Limitations: The code to work with maps requires to digitize every latitude and longitude tick present in the map
#             as support point for maximum accuracy. Between those support points, coordinates are interpolated
#             linearily. For map projections with equal distances between grid lines, this is no problem (in fact,
#             digitizing just two points as in diagrams would theoretically suffice here). For map projections
#             where grid lines are not equidistant the resulting values are still not exact in the current
#             implementation, but their accuracy can be increased by digitizing more support points.
#             The code in its current form only works with map projections in which all grid lines are parallel/
#             perpendicular to each other (e.g. Mercator).

#Author: Manuel Weinkauf (Manuel.Weinkauf@unige.ch)
#Version: 1.4
#Date: 20 February 2020

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.#
#To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/.                   #
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#**************************************************************************************
#Setting working directory
#setwd("C:/R_TestData/Diagrams")

#########################################################################
# Image picking                                                         #
# Necessary input variables:                                            #
#    Dir: Directory path.                                               #
#         *character*                                                   #
#    pattern: Which file types to look for (i.e. the types of your)...  #
#             image files. By default scans for the common file types...#
#             "jpg", "png", "tif", "bmp", "gif", and the common type... #
#             used in moprhometrics "ppm".                              #
#             *character*                                               #
# Output data: Vector of all image files to use.                        #
# Input dataset: Computer directory.                                    #
#########################################################################

ImagePicking<-function (Dir, pattern=c("bmp", "gif", "jpg", "jpeg", "png", "tif", "tiff", "ppm")) {
	#Set up pattern parameter
	pattern<-paste("\\.", pattern, "$", sep="", collapse="|")

	#List all files
	File.list<-list.files(Dir, pattern=pattern, ignore.case=TRUE)
	names(File.list)<-1:length(File.list)
	print(File.list)
	
	#Provide the choice of files by number
	writeLines("Please input the numbers of the files you wish to use. \nSeveral numbers can be input separated by commas. \nDO NOT ENTER SPACES!")
	Image.No<-readline(prompt="Enter the image numbers: ")
	Image.No<-scan(text=Image.No, quiet=TRUE, sep=",")
	Images<-File.list[Image.No]
	
	#Output image names
	Images<-as.data.frame(Images)
	Images[]<-lapply(Images, as.character)
	return(Images)
}

#########################################################################
# Image conversion                                                      #
# Necessary programs: Image Magic                                       #
# Necessary input variables:                                            #
#    Images: A list of all the images to be converted.                  #
#            *data frame* with image names in first column.             #
#    Scaling: Relative scaling of the image during conversion.          #
#             *numeric (integer)*                                       #
#             default=100                                               #
#    OutputType: Image type for output as "xyz".                        #
#                *character*                                            #
# Output data: Images in OutputType format.                             #
# Input dataset: Images to convert.                                     #
#########################################################################

ImageConversion<-function(Images, Scaling=100, OutputType) {
	#Test data consistency
	if (!is.data.frame(Images)) {stop("'Images' must be a data frame!")}
	if (ncol(Images)!=1) {warning("'Images' has more than one column. Only first column will be used!")}
	
	#Read image list into vector
	Images<-as.vector(Images[,1])
	
	#Convert images
	for (i in 1:length(Images)){
		Name<-strsplit(Images[i], split=".", fixed=TRUE)[[1]][1]
		com<-paste("convert ", Images[i], " ", "-resize ", Scaling, "%", " ", Name, ".", OutputType, sep="")
		shell(com)
	}
}

#########################################################################
# Function to extract point coordinates from a given diagram            #
# Necessary packages: pixmap, rtiff, jpeg, png, sp                      #
# Necessary input variables:                                            #
#    Image: Name of input image. It can be in one the following...      #
#           formats: .ppm, .tif, .jpg, .png                             #
#           *character*                                                 #
#    Output: Name of results object.                                    #
#            *character* if Export==TRUE                                #
#    Export: Shall results be exported as file.                         #
#            *logical*                                                  #
#            TRUE: Store results in file                                #
#            FALSE: Store results in variable                           #
#            default=TRUE                                               #
#    LogX: Is the x-axis logarithmically scaled?                        #
#          *logical*                                                    #
#          default=FALSE                                                #
#    LogY: Is the y-axis logarithmically scaled?                        #
#          *logical*                                                    #
#          default=FALSE                                                #
#    Equidist: Should a set of digitized points be converted to a...    #
#              line with equidistant points? This implicitly assumes... #
#              that a curve (rather than distinct points) has been...   #
#              digitized, and that the difference in sampling density...#
#              is an artifact of the manual tracing and not...          #
#              appreciated. If wanted, Equidist must give the number... #
#              of points in the interpolated curve.                     #
#              *numeric (integer)*                                      #
#              default=NULL                                             #
#    Map: Is the plot a map (TRUE) or a diagram (FALSE)?                #
#         *logical*                                                     #
#         default=FALSE                                                 #
# Output data: Matrix with x and y coordinates of the digitized...      #
#              points in diagram or map.                                #
# Input dataset: Digitised image of diagram or map for which data...    #
#                point coordinates should be extracted.                 #
#########################################################################

#Loading packages
require(pixmap)
require(rtiff)
require(jpeg)
require(png)
require(sp)

CoordExt<-function (Image, Output=NULL, Export=TRUE, LogX=FALSE, LogY=FALSE, Equidist=NULL, Map=FALSE) {
	#Find image type
	Type<-unlist(strsplit(Image, "[.]"))
	Type<-Type[length(Type)]
	if (!Type%in%c("ppm", "jpg", "jpeg", "tif", "tiff", "png")) {stop("Image type must be either of .ppm, .jpg, .png, or .tif!")}
	
	#Test for consistency
	if (Export==TRUE & is.null(Output)) {stop("Output must be specified if Export==TRUE!")}
	if (!is.null(Equidist) & !is.numeric(Equidist)) {stop("Equidist must be a number if provided!")}
	
	#Reading and plotting image
	{if (Type=="ppm") {
		y<-read.pnm(Image)
		par(mar=c(1, 1, 1, 1))
		plot(y)
	}
	else if (Type=="tif" | Type=="tiff") {
		y<-readTiff(Image)
		par(mar=c(1, 1, 1, 1))
		plot(y)
	}
	else if (Type=="jpg" | Type=="jpeg") {
		y<-readJPEG(Image)
		par(mar=c(1, 1, 1, 1))
		plot(1, 1, type="n", xlim=c(0, dim(y)[2]), ylim=c(0, dim(y)[1]), asp=1, axes=FALSE)
		rasterImage(y, xleft=0, ybottom=0, xright=dim(y)[2], ytop=dim(y)[1])
	}
	else if (Type=="png") {
		y<-readPNG(Image)
		transparent<-y[,,4]==0
		y<-as.raster(y[,,1:3])
		y[transparent]<-NA
		par(mar=c(1, 1, 1, 1))
		plot(1, 1, type="n", xlim=c(0, dim(y)[2]), ylim=c(0, dim(y)[1]), asp=1, axes=FALSE)
		rasterImage(y, xleft=0, ybottom=0, xright=dim(y)[2], ytop=dim(y)[1], interpolate=FALSE)
	}
	}
	
	#Setting x-scale
	{if (Map==FALSE) {
		XMat<-matrix(NA, 2, 2)
		XScale<-list(m=NA, b=NA)
		writeLines("Please choose the x-axis scale. \nClick on two arbitrary points along the x-axis.")
		flush.console()
		X<-locator(2)
		XMat[1,1]<-X$x[1]
		XMat[2,1]<-X$x[2]
		bringToTop(-1)
		XMat[1,2]<-as.numeric(readline("Please give now the x-value for the first point:"))
		XMat[2,2]<-as.numeric(readline("Please give now the x-value for the second point:"))
		#Calculate scale
		{if (LogX==FALSE) {
			XScale$m<-(XMat[2,2]-XMat[1,2])/(XMat[2,1]-XMat[1,1])
			XScale$b<-XMat[1,2]-XScale$m*XMat[1,1]
		}
		else {
			XScale$m<-(log(XMat[2,2], base=10)-log(XMat[1,2], base=10))/(XMat[2,1]-XMat[1,1])
			XScale$b<-log(XMat[1,2], base=10)-XScale$m*XMat[1,1]
		}
		}
	}
	else {
		bringToTop(-1)
		XTicks<-as.numeric(readline("Please give the number of ticks along the longitude axis:"))
		XMat<-matrix(NA, XTicks, 2)
		XScale<-list()
		writeLines(paste("Please choose the longitude scale. \nClick on ", XTicks, " points from left to right along the x-axis.", sep=""))
		flush.console()
		X<-locator(XTicks)
		XMat[,1]<-X$x
		bringToTop(-1)
		for (i in 1:XTicks) {
			XMat[i,2]<-as.numeric(readline(paste("Please give now the lon-value for the ", i, " point:", sep="")))
		}
		#Calculate scale
		XScale$m<-XScale$b<-vector(length=(nrow(XMat)-1), mode="numeric")
		for (i in 2:nrow(XMat)) {
			XScale$m[i-1]<-(XMat[i,2]-XMat[i-1,2])/(XMat[i,1]-XMat[i-1,1])
			XScale$b[i-1]<-XMat[i-1,2]-XScale$m[i-1]*XMat[i-1,1]
		}
		
		#Check for correct order
		if (any(diff(XMat[,1])<0)) {stop("Longitude points must be digitized strictly from left to right!")}
	}
	}
	
	#Setting y-scale
	{if (Map==FALSE) {
		YMat<-matrix(NA, 2, 2)
		YScale<-list(m=NA, b=NA)
		writeLines("Please choose the y-axis scale. \nClick on two arbitary points along the y-axis.")
		flush.console()
		Y<-locator(2)
		YMat[1,1]<-Y$y[1]
		YMat[2,1]<-Y$y[2]
		bringToTop(-1)
		YMat[1,2]<-as.numeric(readline("Please give now the y-value for the first point:"))
		YMat[2,2]<-as.numeric(readline("Please give now the y-value for the second point:"))
		#Calculate scale
		{if (LogY==FALSE) {
			YScale$m<-(YMat[2,2]-YMat[1,2])/(YMat[2,1]-YMat[1,1])
			YScale$b<-YMat[1,2]-YScale$m*YMat[1,1]
		}
		else {
			YScale$m<-(log(YMat[2,2], base=10)-log(YMat[1,2], base=10))/(YMat[2,1]-YMat[1,1])
			YScale$b<-log(YMat[1,2], base=10)-YScale$m*YMat[1,1]
		}
		}
	}
	else {
		bringToTop(-1)
		YTicks<-as.numeric(readline("Please give the number of ticks along the latitude axis:"))
		YMat<-matrix(NA, YTicks, 2)
		YScale<-list()
		writeLines(paste("Please choose the latitude scale. \nClick on ", YTicks, " points from bottom to top along the y-axis.", sep=""))
		flush.console()
		Y<-locator(YTicks)
		YMat[,1]<-Y$y
		bringToTop(-1)
		for (i in 1:YTicks) {
			YMat[i,2]<-as.numeric(readline(paste("Please give now the lat-value for the ", i, " point:", sep="")))
		}
		#Calculate scale
		YScale$m<-YScale$b<-vector(length=(nrow(YMat)-1), mode="numeric")
		for (i in 2:nrow(YMat)) {
			YScale$m[i-1]<-(YMat[i,2]-YMat[i-1,2])/(YMat[i,1]-YMat[i-1,1])
			YScale$b[i-1]<-YMat[i-1,2]-YScale$m[i-1]*YMat[i-1,1]
		}
		
		#Check for correct order
		if (any(diff(YMat[,1])<0)) {stop("Latitude points must be digitized strictly from bottom to top!")}
	}
	}
	
	#Setting up the results table and extracting data point coordinates
	par(cex=2)
	writeLines("Start Digitizing. \nPlease click on all points to be digitized, one after another. \nWhen you are finished, right-click and chose stop.")
	flush.console()
	C<-locator(type="p", pch=4, col="red", lwd=2)
	Res.Raw<-matrix(unlist(C), length(C$x), 2)
	#Interpolate curve if required
	if (!is.null(Equidist)) {
		Equidist<-round(Equidist, digits=0)
		Res.Line<-Line(Res.Raw)
		Res.Raw<-rbind(spsample(Res.Line, n=Equidist, type="regular", offset=0)@coords, Res.Raw[nrow(Res.Raw),])
	}
	
	Res<-matrix(NA, nrow(Res.Raw), 2)
	{if (Map==FALSE) {colnames(Res)<-c("X", "Y")}
	else {colnames(Res)<-c("Lon", "Lat")}}
	
	#Scale data points to axes scales
	{if (Map==FALSE) {
		for (i in 1:nrow(Res.Raw)) {
			Res[i,1]<-Res.Raw[i,1]*XScale$m+XScale$b
			Res[i,2]<-Res.Raw[i,2]*YScale$m+YScale$b
		}
	}
	else {
		for (i in 1:nrow(Res.Raw)) {
			Int.x<-findInterval(Res.Raw[i,1], XMat[,1])
			if (Int.x==0) {Int.x<-1}
			if (Int.x==nrow(XMat)) {Int.x<-nrow(XMat)-1}
			Int.y<-findInterval(Res.Raw[i,2], YMat[,1])
			if (Int.y==0) {Int.y<-1}
			if (Int.y==nrow(YMat)) {Int.y<-nrow(YMat)-1}
			Res[i,1]<-Res.Raw[i,1]*XScale$m[Int.x]+XScale$b[Int.x]
			Res[i,2]<-Res.Raw[i,2]*YScale$m[Int.y]+YScale$b[Int.y]
		}
	}
	}
	
	#Transform logarithmic data back to cartesian
	if (LogX==TRUE) {Res[,1]<-10^Res[,1]}
	if (LogY==TRUE) {Res[,2]<-10^Res[,2]}
	
	#Exporting results
	{if (Export==TRUE) {write.table(Res, Output, sep="\t")}
	else {return(Res)}}
}

#--------------------------------------------

#Example

#Image conversion
#ImageConversion("Diagram", 1, 1, ".jpg", ".ppm")

#Coordinate extraction from diagrams
#Expected values
#	x	y
#	1	291.6384
#	2	1721.8860
#	3	834.5194
#	4	292.4391
#	5	344.1404
#	6	809.0723
#	7	219.9384
#	8	1553.3091
#	9	1959.9802
#	10	955.9337

#TT1<-CoordExt("Plot_Linear.jpg", Export=FALSE)
#TT2<-CoordExt("Plot_Logarithmic.tif", Export=FALSE, LogY=TRUE)

#Coordinate extraction from maps
#TT3<-CoordExt("GeoMap.png", Export=FALSE, Map=TRUE)#For the core positions
#TT4<-CoordExt("GeoMap.png", Export=FALSE, Map=TRUE, Equidist=100)#For the oceanic front

#--------------------------------------------
#--------------------------------------------
#Version History
#1.0	Finished Program
#1.1	Added the possibility not to export the results
#1.2	Added the functionality to also work with logarithmic axes
#1.3	Changed handling of log axes, added functionality for more image types,
#	implemented points extraction for geographical maps
#1.4	Updated image converter corresponding to morphometric functions
#--------------------------------------------
#--------------------------------------------
