#Extracting morphometric data from images in R
#Input data set: Images, possibly binarized

#Author: Manuel Weinkauf  (Manuel.Weinkauf@unige.ch)
#Version: 1.2.1
#Date: 15 May 2019

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.#
#To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/.                   #
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#**************************************************************************************
#Setting working dierctory
#setwd("C:/R_TestData/GeometricMorphometrics")

#########################################################################
# Image conversion                                                      #
# Necessary programs: Image Magic                                       #
# Necessary input variables:                                            #
#    ImageName: Name part of images (same for all immages).             #
#               *character*                                             #
#    StartNum: Smallest of row of continuous numbers used to number...  #
#              images.                                                  #
#              *numeric (integer)*                                      #
#    StopNum: Largest of row of continuous numbers used to number...    #
#             images.                                                   #
#             *integer*                                                 #
#    Scaling: Relative scaling of the image during conversion.          #
#             *numeric (integer)*                                       #
#             default=100                                               #
#    ImageType: Image type as ".xyz".                                   #
#               *character*                                             #
#    OutputType: Image type for output as ".xyz".                       #
#                *character*                                            #
# Output data: Images in OutputType format.                             #
# Input dataset: Images in ImageType format.                            #
#########################################################################

ImageConversion<-function(ImageName, StartNum, StopNum, Scaling=100, ImageType, OutputType) {
	for (i in (StartNum:StopNum)){
		com<-paste("convert ", ImageName, i, ImageType, " ", "-resize ", Scaling, "%", " ", ImageName, i, OutputType, sep="")
		shell(com)
	}
}

#########################################################################
# Outline rastering                                                     #
# based on Claude (2008), p. 47                                         #
# Necessary input variables:                                            #
#    x: Manually chosen starting point                                  #
#       *integer*                                                       #
#    imagematrix: Name of image to use                                  #
#                 *string*                                              #
# Output data: List of x-y coordinates of all pixel along outline.      #
# Input dataset: Image in ppm format.                                   #
#########################################################################

Conte<-function (x, imagematrix) {
	#Reading image as image matrix and calculate dimensions
	I<-imagematrix
	x<-rev(x)
	x[1]<-dim(I)[1]-x[1]
	
	#Going from the starting point to the left, until reaching a pixel that significantly differs in grey value from the pixel next to it
	while (abs(I[x[1], x[2]]-I[x[1], (x[2]-1)])<0.1) {x[2]<-x[2]-1}
	#Defining that pixel as starting point of the outline
	a<-1
	
	#Writing indices of eight pixels surrounding starting position (0,0) into matrix
	M<-matrix(c(0, -1, -1, -1, 0, 1, 1, 1, 1, 1, 0, -1, -1, -1, 0, 1), 2, 8, byrow=TRUE)
	#Attaching first and last column of matrix as ninth and tenth column on matrix
	M<-cbind(M[,8], M, M[,1])
	
	#Initialising parameters
	X<-0
	Y<-0
	x1<-x[1]
	x2<-x[2]
	SS<-NA
	S<-6
	
	#Allocating all pixels that belong to the outline (running clockwise)
	while ((any(c(X[a], Y[a])!=c(x1, x2)) | length(X)<3)) {
		if(abs(I[x[1]+M[1,S+1], x[2]+M[2,S+1]]-I[x[1], x[2]])<0.1) {a<-a+1; X[a]<-x[1]; Y[a]<-x[2]; x<-x+M[,S+1]; SS[a]<-S+1; S<-(S+7)%%8}
		else if (abs(I[x[1]+M[1,S+2], x[2]+M[2,S+2]]-I[x[1], x[2]])<0.1) {a<-a+1; X[a]<-x[1]; Y[a]<-x[2]; x<-x+M[,S+2]; SS[a]<-S+2; S<-(S+7)%%8}
		else if (abs(I[x[1]+M[1,(S+3)], x[2]+M[2,(S+3)]]-I[x[1], x[2]])<0.1) {a<-a+1; X[a]<-x[1]; Y[a]<-x[2]; x<-x+M[,(S+3)]; SS[a]<-S+3; S<-(S+7)%%8}
		else S<-(S+1)%%8
	}
	list(X = (Y[-1]), Y = ((dim(I)[1]-X))[-1])
}

#########################################################################
# Curvelinear fitting                                                   #
# based on Claude (2008), p. 52                                         #
# Necessary input variables:                                            #
#    Input: Raw pixel coordinates of outline as extracted with...       #
#           function Conte.                                             #
#           *vector*                                                    #
#    n: Desired number of curvelinearily equally spaced points along... #
#       the  outline.                                                   #
#       *integer*                                                       #
# Output data: Matrix with equally spaced coordinates of outline.       #
# Input dataset: Vector of raw x and y coordinates.                     #
#########################################################################

EquiDist<-function (Input, n) {
	RcEqX<-(Input$X[seq(1, length(Input$X), length=(n+1))])[-1]
	RcEqY<-(Input$Y[seq(1, length(Input$Y), length=(n+1))])[-1]
	X_Temp<-as.matrix(RcEqX)
	Y_Temp<-as.matrix(RcEqY)
	M<-cbind(X_Temp, Y_Temp)
	colnames(M)<-c("X", "Y")
	rownames(M)<-paste("Coord", 1:n, sep="")
	M
}

#########################################################################
# Smoothing outline                                                     #
# based on Claude (2008), p. 55                                         #
# Necessary input variables:                                            #
#    M: Matrix of equally spaced outline points.                        #
#       *matrix*                                                        #
#    n: Number of iterations.                                           #
#       *integer*                                                       #
# Output data: List of smoothed x-y coordinates of outline.             #
# Input dataset: Curvelinear equally spaced x-y coordinates of outline. #
#########################################################################

smoothout<-function(M, n){
	p<-dim(M)[1]
	a<-0
	while(a<=n){
		a<-a+1
		Ms<-rbind(M[p,], M[-p,])
		Mi<-rbind(M[-1,], M[1,])
		M<-M/2+Ms/4+Mi/4
	}
	M
}

#########################################################################
# Outline extraction                                                    #
# Necessary packages: pixmap, rtiff, splancs                            #
# Necessary functions: Conte, EquiDist, smoothout                       #
# Necessary input variables:                                            #
#    ImageName: Name part of images (same for all images).              #
#               *string*                                                #
#    InputType: Type of input images. It can be in one the following... #
#               formats: "ppm", "tif"/"tiff"                            #
#               *character*                                             #
#               default="ppm"                                           #
#    StartNum: Smallest of row of continuous numbers used to number...  #
#              images.                                                  #
#              *integer*                                                #
#    StopNum: Largest of row of continuous numbers used to number...    #
#             images.                                                   #
#             *integer*                                                 #
#    Output: Name of the output file (excluding extension).             #
#            *string*                                                   #
#    Specimen.Labels: Names for the specimens. If NULL, numbers will... #
#                     be used.                                          #
#                     default=NULL                                      #
#    OutlinePoints: Desired number of equidistantly spaced points...    #
#                   along outline (see function EquiDist).              #
#                   *integer*                                           #
#                   default=100                                         #
#    Smoothing: Desired number of iterations for outline smoothing...   #
#               (see function smoothout).                               #
#               *integer*                                               #
#               default=1                                               #
#    Baseline: If True, you have to provide two points which can...     #
#              serve as baseline for uniform rotation during the NEF... #
#              (with option Rotation="Baseline"). Note that the two...  #
#              points have to be digitised in the same order in each... #
#              image.                                                   #
#              *logical*                                                #
#              TRUE=Extract baseline coordinates                        #
#              FALSE=Do not extract baseline coordinates (uniform...    #
#                    rotation will be achieved on the basis of the...   #
#                    longest axis of the first harmonic during NEF).    #
#                    deafult=FALSE                                      #
#    Scale: Do the images contain a scale bar that should should be...  #
#           used to extract the size of the object as cross-sectional...#
#           area?                                                       #
#           *logical*                                                   #
#           TRUE=Include step to extract size information               #
#           FALSE=Extract outline coordinates only                      #
#           default=FALSE                                               #
#    RawVersion: Defines whether or not the raw (i.e. un smoothed...    #
#                version of the outline coordinates should be ex-...    #
#                ported as well.                                        #
#                *logical*                                              #
#                TRUE=Export raw version                                #
#                FALSE=Do not export raw version                        #
#                default=FALSE                                          #
# Output data: Matrix containing x and y coordinates of specified...    #
#              number of equidistant points along outline.              #
# Input dataset: Several logically named image files in .ppm format.    #
#########################################################################

#Loading packages
require(pixmap)
require(rtiff)
require(splancs)

OutlineExtraction<-function (ImageName, InputType="ppm", StartNum, StopNum, Output, Specimen.Labels=NULL, OutlinePoints=100, Smoothing=1, Baseline=FALSE, Scale=FALSE, RawVersion=FALSE) {
	#Test data for consistency
	if (!InputType%in%c("ppm", "tif", "tiff")) {stop("Image type must be either of .ppm or .tif!")}
	if (!is.null(Specimen.Labels) & length(Specimen.Labels)!=StopNum-StartNum+1) {stop("Specimen.Labels must have same length as number of images!")}
	
	#Setting up control matrix
	ExtFail<-matrix(NA, (StopNum-StartNum)+1, 1)
	colnames(ExtFail)<-c("Success")
	rownames(ExtFail)<-paste(ImageName, c(StartNum:StopNum), sep="_")
	
	#Setting up sizes matrix
	if (Scale==TRUE) {
		Sizes<-matrix(NA, (StopNum-StartNum)+1, 1)
		colnames(Sizes)<-("Area")
		rownames(Sizes)<-paste(ImageName, c(StartNum:StopNum), sep="_")
	}
	
	#Setting up baseline matrix
	if (Baseline==TRUE) {
		Base<-matrix(NA, (StopNum-StartNum)+1, 4)
		colnames(Base)<-c("x1", "y1", "x2", "y2")
		rownames(Base)<-paste(ImageName, c(StartNum:StopNum), sep="_")
	}
	
	#Setting up results matrix
	if (RawVersion==TRUE) {CoordRes.Raw<-matrix(NA, (StopNum-StartNum+1), OutlinePoints*2)}
	CoordRes<-matrix(NA, (StopNum-StartNum+1), OutlinePoints*2)
	xc<-seq.int(from=1, to=OutlinePoints*2, by=2)
	yc<-seq.int(from=2, to=OutlinePoints*2, by=2)
	
	#Start data acquisition
	MatPos<-1
	for (i in StartNum:StopNum) {
		print(i)
		Image<-paste(ImageName, i, ".", InputType, sep="")
		#Reading image
		{if (InputType=="ppm") {
			y<-read.pnm(Image)
		}
		else if (InputType=="tif" | InputType=="tiff") {
			y<-readTiff(Image)
		}}
		#Converting image to grey scale
		y<-as(y, "pixmapGrey")
		#Converting greyscale to grey on white
		y@grey[which(y@grey>=0.9)]<-1
		y@grey[which(y@grey<0.9)]<-0.7
		
		#Setting margins and plotting image
		par(mar=c(1, 1, 1, 1))
		plot(y)
		
		#Manually chosing starting point and decide whether or not outline is well fitted
		cont<-NA
		while(any(is.na(cont), cont=="n", (cont!="y" && cont!="c"))) {
			#Set scale
			if (Scale==TRUE) {
				writeLines("Please provide scale of the image. \nClick on two points on the scalebar.")
				flush.console()
				a<-locator(2, type="o", pch=8, lwd=2, col="grey60", lty=11)
				scale.px<-sqrt(sum(diff(a$x)**2+diff(a$y)**2))
				bringToTop(-1)
				scale.length<-as.numeric(readline("How long is the implied line? "))
			}
		
			#Extract baseline
			if (Baseline==TRUE) {
				writeLines("Please provide baseline. \nClick on two landmarks of the object (order sensitive).")
				flush.console()
				a<-locator(2, type="o", pch=8, lwd=2, col="grey60", lty=11)
				Base[MatPos,1]<-a$x[1]
				Base[MatPos,2]<-a$y[1]
				Base[MatPos,3]<-a$x[2]
				Base[MatPos,4]<-a$y[2]
			}
		
			#Extract outline
			print("Click within the object to the right of the starting point!")
			flush.console()
			start<-locator(1)
			Rc<-Conte(c(round(start$x), round(start$y)), y@grey)
			lines(Rc$X, Rc$Y, lwd=4)
			arrows(0, Rc$Y[1], Rc$X[1], Rc$Y[1], length=0.1)
			bringToTop(-1)
			cont<-readline("Is the outline correct? (y=proceed,n=try again,c=cancel and proceed)")
			if (cont!="y" && cont!="c") {par(mar=c(1, 1, 1, 1)); plot(y)}
		}
		{if (cont=="y") {ExtFail[MatPos,1]<-1}
			else {ExtFail[MatPos,1]<-0}
		}

		#Normalising curve for curvelinear equally spaced number of points
		RcEqual<-EquiDist(Rc, OutlinePoints)
		
		#Smoothing outline
		RcSmooth<-smoothout(RcEqual, Smoothing)
		
		#Calculating sizes
		if (Scale==TRUE) {
			{if (cont=="y") {
				RcScaled<-cbind(Rc$X,Rc$Y)
				for (j in 2:nrow(RcScaled)) {
					XD<-RcScaled[j,1]-RcScaled[1,1]
					YD<-RcScaled[j,2]-RcScaled[1,2]
					RcScaled[j,1]<-RcScaled[1,1]+(XD/scale.px)*scale.length
					RcScaled[j,2]<-RcScaled[1,2]+(YD/scale.px)*scale.length
				}
				RcScaled<-as.matrix(rbind(RcScaled, RcScaled[1,]))
				Sizes[MatPos,1]<-areapl(RcScaled)
			}
			else {Sizes[MatPos,1]<-NA}
			}
		}
		
		#Writing x-y coordinates of outline into variable
		{if (cont=="y") {
			if (RawVersion==TRUE) {
				CoordRes.Raw[i,xc]<-RcEqual[,"X"]
				CoordRes.Raw[i,yc]<-RcEqual[,"Y"]
			}
			CoordRes[i,xc]<-RcSmooth[,"X"]
			CoordRes[i,yc]<-RcSmooth[,"Y"]
		}
		else if (cont=="c") {
			CoordRes[i,]<--999
		}
		}
		
		MatPos<-MatPos+1
	}
	
	#Save outline coordinates as NTS file
	if (RawVersion==TRUE) {
		FileName<-paste(Output, "_Raw.nts", sep="")
		{if (any(ExtFail[,1]==0)) {firstl<-paste(1, paste(StopNum-StartNum+1, "L", sep=""), OutlinePoints*2, 1, -999, "dim=2", sep=" ")}
		else {firstl<-paste(1, paste(StopNum-StartNum+1, "L", sep=""), OutlinePoints*2, 0, "dim=2", sep=" ")}}
		{if (is.null(Specimen.Labels)) {L<-StartNum:StopNum}
		else {L<-Specimen.Labels}}
		secondl<-paste(L, sep="", collapse=" ")
		##Create file and write header
		cat(firstl, secondl, file=FileName, sep="\n", append=FALSE)
		##Create data body
		for (i in 1:(nrow(CoordRes.Raw))) {
			B<-paste(CoordRes.Raw[i,], sep="", collapse=" ")
			cat(B, file=FileName, sep="\n", append=TRUE)
		}
	}
	FileName<-paste(Output, ".nts", sep="")
	{if (any(ExtFail[,1]==0)) {firstl<-paste(1, paste(StopNum-StartNum+1, "L", sep=""), OutlinePoints*2, 1, -999, "dim=2", sep=" ")}
	else {firstl<-paste(1, paste(StopNum-StartNum+1, "L", sep=""), OutlinePoints*2, 0, "dim=2", sep=" ")}}
	{if (is.null(Specimen.Labels)) {L<-StartNum:StopNum}
	else {L<-Specimen.Labels}}
	secondl<-paste(L, sep="", collapse=" ")
	##Create file and write header
	cat(firstl, secondl, file=FileName, sep="\n", append=FALSE)
	##Create data body
	for (i in 1:(nrow(CoordRes))) {
		B<-paste(CoordRes[i,], sep="", collapse=" ")
		cat(B, file=FileName, sep="\n", append=TRUE)
	}
	
	#Save list of successes, i.e. in which specimens did outline extraction fail, size information, and baseline coordinates
	FileName<-paste(Output, "_Success.txt", sep="")
	write.table(ExtFail, FileName, sep="\t")
	if (Scale==TRUE) {FileName<-paste(Output, "_Area.txt", sep=""); write.table(Sizes, FileName, sep="\t")}
	if (Baseline==TRUE) {FileName<-paste(Output, "_Baseline.txt", sep=""); write.table(Base, FileName, sep="\t")}
}

#########################################################################
# Landmark extraction                                                   #
# Necessary packages: pixmap, rtiff                                     #
# Necessary input variables:                                            #
#    Image: Name-part of the images (same for all images).              #
#           *character*                                                 #
#    InputType: Type of input images. It can be in one the following... #
#               formats: "ppm", "tif"/"tiff"                            #
#               *character*                                             #
#               default="ppm"                                           #
#    StartNum: Smallest of row of continuous numbers used to number...  #
#              images.                                                  #
#              *numeric (integer)*                                      #
#    StopNum: Largest of row of continuous numbers used to number...    #
#             images.                                                   #
#             *numeric (integer)*                                       #
#    Output: Name of the output file (excluding extension).             #
#            *character*                                                #
#    Specimen.Labels: Names for the specimens. If NULL, numbers will... #
#                     be used.                                          #
#                     default=NULL                                      #
#    Scale: Does the image contain a scale for which the the point...   #
#           coordinates should be normalized?                           #
#           NOTE: if Scale==TRUE and Export=="NTS" the coordinates...   #
#           will be exported in scaled format. However, if...           #
#           Export=="TPS" the coordinates will be exported as pixel...  #
#           coordinates and the SCALE= parameter will be included for...#
#           later conversion.                                           #
#           *logical*                                                   #
#           TRUE: Scale coordinates                                     #
#           FALSE: Leave coordinates as they are                        #
#           default=TRUE                                                #
#    ScaleParam: Vector containing the known 1 px=y units conversion... #
#                for each image.                                        #
#                *numeric (real)*, length=number of images              #
#                default=NULL                                           #
#    N: Number of landmark points to extract.                           #
#       *numeric (integer)*                                             #
#    Export: In which format should results be exported? It can be...   #
#            one of the following: "NTS", "TPS"                         #
#            *character*                                                #
#            default="TPS"                                              #
# Output data: List of x-y coordinates of landmarks in .nts or .tps...  #
#              file.                                                    #
# Input dataset: Images in InputType format.                            #
#########################################################################

#Load packages
require(pixmap)
require(rtiff)

LMExtract<-function(Image, InputType="ppm", StartNum, StopNum, Output, Specimen.Labels=NULL, Scale=TRUE, ScaleParam=NULL, N, Export="TPS") {
	#Test data for consistency
	if (!InputType%in%c("ppm", "tif", "tiff")) {stop("Image type must be either of .ppm or .tif!")}
	if (Export!="NTS" & Export!="TPS") {stop("Export format must be either NTS or TPS!")}
	if (Scale==TRUE & !is.null(ScaleParam)) {Scale<-FALSE; warning("Cannot have scalebar and scale parameter at the same time, Scale has been set to FALSE")}
	if (!is.null(ScaleParam) & length(ScaleParam)!=StopNum-StartNum+1) {stop("Must supply one value of ScaleParam per image")}
	if (!is.null(Specimen.Labels) & length(Specimen.Labels)!=StopNum-StartNum+1) {stop("Specimen.Labels must have same length as number of images!")}
	
	#Set up results matrices
	{if (Export=="NTS") {
		CoordRes<-matrix(NA, (StopNum-StartNum+1), N*2)
		xc<-seq.int(from=1, to=N*2, by=2)
		yc<-seq.int(from=2, to=N*2, by=2)
	}
	else {
		CoordRes<-array(NA, dim=c(N, 2, StopNum-StartNum+1))
		Meta<-list()
		Meta$Files<-vector(mode="character", length=StopNum-StartNum+1)
		Meta$ID<-vector(mode="numeric", length=StopNum-StartNum+1)
		Meta$Scale<-vector(mode="numeric", length=StopNum-StartNum+1)
	}}
	
	#Setup success report matrix
	ExtFail<-matrix(NA, (StopNum-StartNum)+1, 1)
	rownames(ExtFail)<-paste(Image, StartNum:StopNum, sep=".")
	
	for (i in StartNum:StopNum) {
		print(i)
		FileName<-paste(Image, i, ".", InputType, sep="")

		#Read and plot image
		{if (InputType=="ppm") {
			y<-read.pnm(FileName)
		}
		else if (InputType=="tif" | InputType=="tiff") {
			y<-readTiff(FileName)
		}}
		cont<-NA
		while(any(is.na(cont), cont=="n", (cont!="y" && cont!="c"))){
			par(mar=(c(1, 1, 1, 1)))
			plot(y)
	
			#Set scale
			if (Scale==TRUE) {
				writeLines("Please provide scale of the image. \nClick on two points on the scalebar.")
				flush.console()
				a<-locator(2, type="o", pch=8, lwd=2, col="grey60", lty=11)
				scale.px<-sqrt(sum(diff(a$x)**2+diff(a$y)**2))
				bringToTop(-1)
				scale.length<-as.numeric(readline("How long is the implied line? "))
			}
	
			#Extract and label landmarks
			print(paste("Digitize ", N, " landmarks by clicking into the image!"))
			flush.console()
			LM<-locator(n=N, type="p", pch=3, col="red", lwd=2)
			text(LM, pos=2, labels=1:N, col="red", font=2)
			
			bringToTop(-1)
			cont<-readline("Are the landmarks correct? (y=proceed, n=try again, c=cancel and proceed)")
		}
		
		#Save copy of image with points for later comparison
		dev.copy(png, filename=paste(Image, i, "_Points.png", sep=""));
		dev.off();
		
		#Write success report
		{if (cont=="y") {ExtFail[i,1]<-1}
		else {ExtFail[i,1]<-0}}
	
		#Write landmark coordinates into table
		{if (Scale==TRUE & Export=="NTS") {LM$x<-LM$x*(scale.length/scale.px); LM$y<-LM$y*(scale.length/scale.px)}
		else if (Scale==FALSE & !is.null(ScaleParam) & Export=="NTS") {LM$x<-LM$x*ScaleParam[i]; LM$y<-LM$y*ScaleParam[i]}
		else if (Scale==TRUE & Export=="TPS") {ScaleTPS=scale.length/scale.px}
		else if (Scale==FALSE & !is.null(ScaleParam) && Export=="TPS") {ScaleTPS=ScaleParam[i]}}
	
		{if (Export=="NTS") {
			{if (cont=="y") {
				CoordRes[i,xc]<-LM$x
				CoordRes[i,yc]<-LM$y
			}
			else if (cont=="c") {
				CoordRes[i,]<--999
			}
			}
		}
		else {
			if (cont=="y") {
				CoordRes[,,i]<-cbind(LM$x, LM$y)
				if (Scale==TRUE | !is.null(ScaleParam)) {Meta$Scale[i]<-ScaleTPS}
			}
			Meta$Files[i]<-FileName
			{if (is.null(Specimen.Labels)) {Meta$ID[i]<-i}
			else {Meta$ID[i]<-Specimen.Labels[i]}}
		}}
	}
	
	#Export coordinates
	{if (Export=="NTS") {
		FileName<-paste(Output, ".nts", sep="")
		{if (any(ExtFail[,1]==0)) {firstl<-paste(1, paste(StopNum-StartNum+1, "L", sep=""), N*2, 1, -999, "dim=2", sep=" ")}
		else {firstl<-paste(1, paste(StopNum-StartNum+1, "L", sep=""), N*2, 0, "dim=2", sep=" ")}}
		{if (is.null(Specimen.Labels)) {L<-StartNum:StopNum}
		else {L<-Specimen.Labels}}
		secondl<-paste(L, sep="", collapse=" ")
		##Create file and write header
		cat(firstl, secondl, file=FileName, sep="\n", append=FALSE)
		##Create data body
		for (i in 1:(nrow(CoordRes))) {
			B<-paste(CoordRes[i,], sep="", collapse=" ")
			cat(B, file=FileName, sep="\n", append=TRUE)
		}
	}
	else {
		FileName<-paste(Output, ".tps", sep="")
		firstl<-paste("LM=", N, sep="")
		##Create file
		for (j in 1:(dim(CoordRes)[3])) {
			cat(firstl, file=FileName, sep="\n", append=TRUE)
			for (i in 1:(dim(CoordRes)[1])) {
				B<-paste(CoordRes[i,,j], sep="", collapse=" ")
				cat(B, file=FileName, sep="\n", append=TRUE)
			}
			cat(paste("IMAGE=", Meta$Files[j], sep=""), file=FileName, sep="\n", append=TRUE)
			cat(paste("ID=", Meta$ID[j], sep=""), file=FileName, sep="\n", append=TRUE)
			if (Scale==TRUE | !is.null(ScaleParam)) {cat(paste("SCALE=", Meta$Scale[j], sep=""), file=FileName, sep="\n", append=TRUE)}
		}
	}
	}
	
	#Export failure report
	write.table(ExtFail, paste(Output, "_SuccessReport.txt", sep=""), sep="\t")
}

#########################################################################
# Plotting image, digitizing points along spiral                        #
# Necessary packages: pixmap, rtiff                                     #
# Necessary input variables:                                            #
#    ImageName: Name part of images (same for all images).              #
#               *character*                                             #
#    InputType: Type of input images. It can be in one the following... #
#               formats: "ppm", "tif"/"tiff"                            #
#               *character*                                             #
#               default="ppm"                                           #
#    StartNum: Smallest of row of continuous numbers used to number...  #
#              images.                                                  #
#              *numeric (integer)*                                      #
#    StopNum: Largest of row of continuous numbers used to number...    #
#             images.                                                   #
#             *numeric (integer)*                                       #
#    Output: Name of the output file (excluding extension).             #
#            *character*                                                #
#    Guidelines: Shall guiders be drawn to enable to take take points...#
#                at equidistant spaces? Not recommended if measuring... #
#                structures that have natural guides (like chambers).   #
#                *logical*                                              #
#                TRUE=Guidelines will be drawn                          #
#                FALSE=No guidelines will be drawn                      #
#                default=TRUE                                           #
#    Density: Sampling density (degrees). Only meaningful if...         #
#             Guidelines=TRUE.                                          #
#             *numeric (integer)*                                       #
#             default=45                                                #
#    Equidistant: Shall points along spiral be equidistant?             #
#                 *logical*                                             #
#                 TRUE=points must be equidistant                       #
#                 FALSE=points do not need to be equidistant            #
#                 default=TRUE                                          #
#                 IMPORTANT: If FALSE is chosen points along the...     #
#                            spiral outline are allowed to be not...    #
#                            equidistant (whether they actually are...  #
#                            depends on the value of Density). This...  #
#                            could likely bias the results of the...    #
#                            further analysis!                          #
#    Normalize: Should the spirals be normalized for unit size before...#
#               coordinate export?                                      #
#               *logical*                                               #
#               TRUE: Normalize size to radius 1                        #
#               FALSE: Leave size as is                                 #
#               default=TRUE                                            #
#    Double: Should a double spiral be extracted, for instance on...    #
#            the external and internal side of the same shell?          #
#            *logical*                                                  #
#            default=FALSE                                              #
# Output data: File of type .spiral (effectively a .csv file without... #
#              row names) containing x- and y-coordinates and polar...  #
#              coordinates (distance from center t, angle theta in...   #
#              radians) for digitized points along spiral.              #
# Input dataset: Images in InputType format.                            #
#########################################################################

#Loading packages
require(pixmap)
require(rtiff)

SpiralExtraction<-function(ImageName, InputType="ppm", StartNum, StopNum, Output, Guidelines=TRUE, Density=45, Equidistant=TRUE, Normalize=TRUE, Double=FALSE){
	#Test data for consistency
	if (!InputType%in%c("ppm", "tif", "tiff")) {stop("Image type must be either of .ppm or .tif!")}
	if(Equidistant==TRUE){
		if(180%%Density!=0){stop(paste("180 cannot be divided by ", Density, " without remainder. Choose another value to get equidistant points along the spiral!", sep=""))}
	}
	
	#Set up temporary results objects
	Res<-list()
	if (Double==TRUE) {Res2<-list()}
	
	#Setup success report matrix
	ExtFail<-matrix(NA, (StopNum-StartNum)+1, 1)
	rownames(ExtFail)<-paste(ImageName, StartNum:StopNum, sep=".")
	
	#Extract spiral
	Ind<-1
	for (k in StartNum:StopNum) {
		print(k)
		Image<-paste(ImageName, k, ".", InputType, sep="")
		
		#Read image
		{if (InputType=="ppm") {
			y<-read.pnm(Image)
		}
		else if (InputType=="tif" | InputType=="tiff") {
			y<-readTiff(Image)
		}}
		cont<-NA
		while(any(is.na(cont), cont=="n", (cont!="y" && cont!="c"))){
			#Plot image
			par(mar=c(1, 1, 1, 1))
			plot(y)
			
			#Set midpoint for spiral
			writeLines("Please choose the center of the spiral. \nClick within the image.")
			flush.console()
			start<-locator(1)
			
			#Plot reference lines
			{if (Guidelines==TRUE) {
				abline(h=start$y, col="pink", lwd=2)
				if(90%%Density==0){abline(v=start$x, col="pink", lwd=2)}
				Ang<-0
				while(Ang<(180-Density)){
					Ang<-Ang+Density
					{if(Ang!=90){abline(a=start$y-start$x*tan(Ang*(pi/180)), b=tan(Ang*(pi/180)), col="pink", lwd=2)}
					else{}}
				}
			}
			else {
				points(start$x, start$y, pch=16, col="black")
			}}
			
			#Digitize spiral outline
			writeLines("Please digitize the points along the spiral. \nWhen you are finished right-click and choose stop.")
			flush.console()
			Coord<-locator(n=1000, type="o", lwd=2, pch=3, col="red")
			if (Double==TRUE) {
				writeLines("Please digitize the points along the other side of the spiral. \nWhen you are finished right-click and choose stop.")
				flush.console()
				Coord2<-locator(n=1000, type="o", lwd=2, pch=3, col="yellow")
			}
			
			bringToTop(-1)
			cont<-readline("Are the landmarks correct? (y=proceed, n=try again, c=cancel and proceed)")
		}
		
		#Save image with spiral on for later comparison
		dev.copy(jpeg, filename=paste(ImageName, k, "_Spiral.jpg", sep=""));
		dev.off();
		
		#Write success report
		{if (cont=="y") {ExtFail[Ind,1]<-1}
		else {ExtFail[Ind,1]<-0}}
		
		#Calculate lengths and angles of data
		Res[[Ind]]<-matrix(NA, length(Coord$x), 4)
		colnames(Res[[Ind]])<-c("x", "y", "t", "theta")
		Res[[Ind]][,"x"]<-Coord$x
		Res[[Ind]][,"y"]<-Coord$y
		for (i in 1:(length(Coord$x))) {
			Res[[Ind]][i,"t"]<-sqrt((Coord$x[i]-start$x)^2+(Coord$y[i]-start$y)^2)
			{if((Coord$x[i]-start$x)>=0 && (Coord$y[i]-start$y)>=0){Res[[Ind]][i,"theta"]<-atan(abs((Coord$y[i]-start$y))/abs((Coord$x[i]-start$x)))}
			else if((Coord$x[i]-start$x)<0 && (Coord$y[i]-start$y)>=0){Res[[Ind]][i,"theta"]<-pi-atan(abs((Coord$y[i]-start$y))/abs((Coord$x[i]-start$x)))}
			else if((Coord$x[i]-start$x)<0 && (Coord$y[i]-start$y)<0){Res[[Ind]][i,"theta"]<-pi+atan(abs((Coord$y[i]-start$y))/abs((Coord$x[i]-start$x)))}
			else {Res[[Ind]][i,"theta"]<-(2*pi)-atan(abs((Coord$y[i]-start$y))/abs((Coord$x[i]-start$x)))}
			}
		}
		if (Double==TRUE) {
			Res2[[Ind]]<-matrix(NA, length(Coord2$x), 4)
			colnames(Res2[[Ind]])<-c("x", "y", "t", "theta")
			Res2[[Ind]][,"x"]<-Coord2$x
			Res2[[Ind]][,"y"]<-Coord2$y
			for (i in 1:(length(Coord2$x))) {
				Res2[[Ind]][i,"t"]<-sqrt((Coord2$x[i]-start$x)^2+(Coord2$y[i]-start$y)^2)
				{if((Coord2$x[i]-start$x)>=0 && (Coord2$y[i]-start$y)>=0){Res2[[Ind]][i,"theta"]<-atan(abs((Coord2$y[i]-start$y))/abs((Coord2$x[i]-start$x)))}
				else if((Coord2$x[i]-start$x)<0 && (Coord2$y[i]-start$y)>=0){Res2[[Ind]][i,"theta"]<-pi-atan(abs((Coord2$y[i]-start$y))/abs((Coord2$x[i]-start$x)))}
				else if((Coord2$x[i]-start$x)<0 && (Coord2$y[i]-start$y)<0){Res2[[Ind]][i,"theta"]<-pi+atan(abs((Coord2$y[i]-start$y))/abs((Coord2$x[i]-start$x)))}
				else {Res2[[Ind]][i,"theta"]<-(2*pi)-atan(abs((Coord2$y[i]-start$y))/abs((Coord2$x[i]-start$x)))}
				}
			}
		}
		
		#Normalize spiral outline
		##Normalize size for radius=1
		if (Normalize==TRUE) {
			Size<-max(c(Res[[Ind]][,"t"], Res2[[Ind]][,"t"]))
			Res[[Ind]][,"t"]<-Res[[Ind]][,"t"]/Size
			if (Double==TRUE) {
				Res2[[Ind]][,"t"]<-Res2[[Ind]][,"t"]/Size
			}
		}
		##Normalize rotation for start-point at radian=0
		Rotation<-Res[[Ind]][1,"theta"]
		for (i in 1:nrow(Res[[Ind]])) {
			{if (Res[[Ind]][i,"theta"]>=Rotation) {Res[[Ind]][i,"theta"]<-Res[[Ind]][i,"theta"]-Rotation}
			else {Res[[Ind]][i,"theta"]<-(2*pi)+(Res[[Ind]][i,"theta"]-Rotation)}}
		}
		if (Double==TRUE) {
			for (i in 1:nrow(Res2[[Ind]])) {
				{if (Res2[[Ind]][i,"theta"]>=Rotation) {Res2[[Ind]][i,"theta"]<-Res2[[Ind]][i,"theta"]-Rotation}
				else {Res2[[Ind]][i,"theta"]<-(2*pi)+(Res2[[Ind]][i,"theta"]-Rotation)}}
			}
		}
		
		#Increase row counter
		Ind<-Ind+1
	}
	
	#Prepare output file
	Res.Final<-matrix(NA, length(Res)*4, max(unlist(lapply(lapply(Res, dim), '[[', 1))))
	rownames(Res.Final)<-paste(c("x", "y", "t", "theta"), rep(StartNum:StopNum, each=4), sep=".")
	for (i in 1:length(Res)) {
		Start.Line<-i+((i-1)*3)
		L<-nrow(Res[[i]])
		Res.Final[Start.Line,1:L]<-Res[[i]][,"x"]
		Res.Final[Start.Line+1,1:L]<-Res[[i]][,"y"]
		Res.Final[Start.Line+2,1:L]<-Res[[i]][,"t"]
		Res.Final[Start.Line+3,1:L]<-Res[[i]][,"theta"]
	}
	if (Double==TRUE) {
		Res.Final2<-matrix(NA, length(Res2)*4, max(unlist(lapply(lapply(Res2, dim), '[[', 1))))
		rownames(Res.Final2)<-paste(c("x", "y", "t", "theta"), rep(StartNum:StopNum, each=4), sep=".")
		for (i in 1:length(Res2)) {
			Start.Line<-i+((i-1)*3)
			L<-nrow(Res2[[i]])
			Res.Final2[Start.Line,1:L]<-Res2[[i]][,"x"]
			Res.Final2[Start.Line+1,1:L]<-Res2[[i]][,"y"]
			Res.Final2[Start.Line+2,1:L]<-Res2[[i]][,"t"]
			Res.Final2[Start.Line+3,1:L]<-Res2[[i]][,"theta"]
		}
	}
	
	#Export results
	{if (Double==FALSE) {write.table(Res.Final, paste(Output, ".spiral", sep=""), sep=",", col.names=FALSE)}
	else {
		write.table(Res.Final, paste(Output, "_SideRed.spiral", sep=""), sep=",", col.names=FALSE)
		write.table(Res.Final2, paste(Output, "_SideYellow.spiral", sep=""), sep=",", col.names=FALSE)
	}
	}
	write.table(ExtFail, paste(Output, "_SuccessReport.txt", sep=""), sep="\t")
}

#--------------------------------------------

#Examples

#Converting images
#ImageConversion("Sp", 1, 3, ImageType=".png", OutputType=".ppm")
#ImageConversion("Sp", 1, 3, ImageType=".png", OutputType=".tif")
#ImageConversion("Spiral", 1, 3, ImageType=".png", OutputType=".ppm")
#ImageConversion("Spiral", 1, 3, ImageType=".png", OutputType=".tif")

#Extract outlines
#setwd("C:/R_TestData/GeometricMorphometrics/Outlines")
#OutlineExtraction("Sp", StartNum=1, StopNum=3, Output="Stars_Rep1", Specimen.Labels=paste("Spec", 1:3, sep="."), OutlinePoints=60, Smoothing=1, Scale=TRUE, RawVersion=TRUE)
#OutlineExtraction("Sp", InputType="tif", StartNum=1, StopNum=3, Output="Stars_Rep2", Specimen.Labels=paste("Spec", 1:3, sep="."), OutlinePoints=70, Smoothing=1, Scale=TRUE, RawVersion=TRUE)

#Extract landmarks
#setwd("C:/R_TestData/GeometricMorphometrics/Landmarks")
#LMExtract(Image="Sp", StartNum=1, StopNum=3, Output="Landmarks1", N=10, Export="NTS")
#LMExtract(Image="Sp", StartNum=1, StopNum=3, ScaleParam=rep(5.84818294686911, 3), Output="Landmarks2", N=10, Export="NTS")
#LMExtract(Image="Sp", InputType="tif", StartNum=1, StopNum=3, Output="Landmarks3_Rep1", N=10, Export="TPS")
#LMExtract(Image="Sp", InputType="tif", StartNum=1, StopNum=3, Output="Landmarks3_Rep2", N=10, Export="TPS")

#Extract spiral morphology
#setwd("C:/R_TestData/GeometricMorphometrics/Spirals")
#SpiralExtraction("Spiral", StartNum=1, StopNum=3, Output="SpiralForm", Guidelines=TRUE, Density=23, Equidistant=TRUE, Normalize=TRUE)
#SpiralExtraction("Spiral", StartNum=1, StopNum=3, Output="SpiralForm", Guidelines=TRUE, Density=45, Equidistant=TRUE, Normalize=TRUE)
#SpiralExtraction("Spiral", InputType="tif", StartNum=1, StopNum=3, Output="SpiralForm", Guidelines=TRUE, Density=45, Equidistant=TRUE, Normalize=TRUE)

#--------------------------------------------
#--------------------------------------------
#Version History
#1.0	Finished Program
#1.1	Added function SpiralExtraction
#1.1.1	Added possibility to provide specimen labels manually in OutlineExtraction and LMExtract
#1.1.2	Numbering of specimens in Spiral.Extraction now based on start and stop number
#1.2	Added functionality to Spiral.Extraction to extract two parallel spirals and enhanced visuals
#1.2.1	Fixed an error where Spiral.Extraction would fail if StartNum was different from 1
#--------------------------------------------
#--------------------------------------------






































