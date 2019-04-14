#Apply traditional morphometrics
#Input data set:
#	-Landmarks stored in a TPS or NTS File
#	-Coordinates of landmarks in a matrix with two columns (x, y) or array
#	-Length measurements
#Further reading:
#	Claude, J. (2008) "Morphometrics with R". Gentleman, R., Hornik, K., and...
#		Parmigiani, G. (eds) "Use R!", vol. ii, 316 pp. (Springer).
#	Jolicoeur, P. (1963) "The degree of generality of robustness in Martes americana"...
#		Growth 27: 1--28.
#	Legendre, P. (sine anno) "Model II Regression User's Guide, R Edition"
#		http://cran.r-project.org/web/packages/lmodel2/vignettes/mod2user.pdf

#Author: Manuel Weinkauf (Manuel.Weinkauf@unige.ch)
#Version: 1.3
#Date: 9 November 2014

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.#
#To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/.                   #
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#Creation of test dataset
#LM<-matrix(c(rep(1, 10), rep(12, 10), rep(8, 10), rep(7, 10)), 10, 4)
#LM<-LM-runif(40)
#colnames(LM)<-c("x1", "y1", "x2", "y2")
#rownames(LM)<-paste("Specimen", 1:10, sep=".")
#plot(LM[,1], LM[,2], col="blue", xlim=c(0, 10), ylim=c(5, 15), asp=TRUE)
#points(LM[,3], LM[,4], col="red")
#LM2<-cbind(LM, c(0.12, 0.12, 0.16, 0.12, 0.22, 0.16, 0.10, 0.22, 0.22, 0.12))
#colnames(LM2)<-c("x1", "y1", "x2", "y2", "Scale")
#LM3<-cbind(LM[1:5,], LM[6:10,])
#rownames(LM3)<-paste("Angle", 1:5, sep=".")

#Bivar<-matrix(c(rep(1:2, 10), 2:21, seq(from=0.2, to=0.77, by=0.03)), 20, 3)
#Bivar[seq(from=1, to=19, by=2),2]<-Bivar[seq(from=1, to=19, by=2),2]+runif(10, min=1, max=1.2)
#Bivar[seq(from=2, to=20, by=2),2]<-Bivar[seq(from=2, to=20, by=2),2]+runif(10, min=0.4, max=3.1)
#Bivar[seq(from=1, to=19, by=2),3]<-Bivar[seq(from=1, to=19, by=2),3]+runif(10, min=0, max=0.01)
#Bivar[seq(from=2, to=20, by=2),3]<-Bivar[seq(from=2, to=20, by=2),3]+runif(10, min=0.02, max=0.06)
#colnames(Bivar)<-c("Species", "Body.Length", "Brain.Size")
#rownames(Bivar)<-paste("Specimen", 1:20, sep=".")

#Multivar<-matrix(c(rep(1:2, 10), 2:21, seq(from=0.2, to=0.77, by=0.03), seq(from=0.03, to=0.0585, by=0.0015)), 20, 4)
#Multivar[seq(from=1, to=19, by=2),2]<-Multivar[seq(from=1, to=19, by=2),2]+runif(10, min=1, max=1.2)
#Multivar[seq(from=2, to=20, by=2),2]<-Multivar[seq(from=2, to=20, by=2),2]+runif(10, min=0.4, max=3.1)
#Multivar[seq(from=1, to=19, by=2),3]<-Multivar[seq(from=1, to=19, by=2),3]+runif(10, min=0, max=0.01)
#Multivar[seq(from=2, to=20, by=2),3]<-Multivar[seq(from=2, to=20, by=2),3]+runif(10, min=0.02, max=0.06)
#Multivar[seq(from=1, to=19, by=2),4]<-Multivar[seq(from=1, to=19, by=2),3]+runif(10, min=0, max=0.001)
#Multivar[seq(from=2, to=20, by=2),4]<-Multivar[seq(from=2, to=20, by=2),3]+runif(10, min=0.002, max=0.006)
#colnames(Multivar)<-c("Species", "Body.Length", "Brain.Size", "Eye.Size")
#rownames(Multivar)<-paste("Specimen", 1:20, sep=".")

#**************************************************************************************
#Setting working directory
#setwd("C:/R_TestData")

#Loading required functions
#source("GeometricMorphometrics_Functions.r")#NOTE: Change directory accordingly, only needed if data are in TPS/NTS format

#########################################################################
# Function to calculate distances between landmarks                     #
# based on Claude (2008), p. 49                                         #
# Necessary input variables:                                            #
#   Data: X and y coordinates of landmarks between which distances...   #
#         should be calculated.                                         #
#         *matrix*                                                      #
#   Scale: Either FALSE, if the input coordinates already correspond... #
#          to the measurement unit, or a real value for conversion...   #
#          (via multiplication), if the coordinates are in pixel with...#
#          a fixed scaling factor to the anticipated unit of...         #
#          measurement. Write either as equation or its result of the...#
#          form unitlength/pixelnumber (e.g. "1/10" or "0.1" if 10px... #
#          equal 1mm).                                                  #
#          *logical* or *numeric (real)*                                #
#          default=FALSE                                                #
#   Scale.Column: Set to true if scaling is different for each...       #
#                 observation, in which case this is provided in the... #
#                 fifth column of Data.                                 #
#                 NOTE: This is a fallback solution that is not...      #
#                       recommended. Since the measurement error is...  #
#                       also a function of magnification it is highly...#
#                       recommended to use the same magnification...    #
#                       (i.e. scaling factor) for all observations.     #
#                 *logical*                                             #
#                 default=TRUE                                          #
# Output data: Data matrix containing the calculated distances for...   #
#              each observation (if scaling factor was provided then... #
#              the raw pixel distances are also returned).              #
# Input dataset: Matrix with four or five columns. Each row...          #
#                corresponds to one observation. Columns 1 and 2 hold...#
#                x and y-coordinates for landmark 1, columns 3 and 4... #
#                hold x and y coordinates of landmark 2. If...          #
#                Scale.Column==TRUE a fifth column providing the...     #
#                scaling factor for each observation separately is...   #
#                expected.                                              #
#########################################################################

#Loading packages

LandmarkDist<-function (Data, Scale=FALSE, Scale.Column=FALSE) {
	#Test data for consistency
	if (Scale.Column==FALSE & ncol(Data)!=4) {stop("Exactly four columns needed in input data")}
	if (Scale.Column==TRUE & ncol(Data)!=5) {stop("Exactly five columns needed in input data")}
	if (is.numeric(Scale) & Scale.Column==TRUE) {stop("Cannot provide global and per-observation scaling-factor at the same time")}
	if (Scale==TRUE) {stop("If Scale!=FALSE it must be provided as numeric value")}
	
	#Setup output matrix
	{if (Scale==FALSE & Scale.Column==FALSE) {Dist<-matrix(NA, nrow(Data), 1); colnames(Dist)<-"Distance.unit"}
	else {Dist<-matrix(NA, nrow(Data), 2); colnames(Dist)<-c("Distance.px", "Distance.unit")}}
	if (!is.null(rownames(Data))) {rownames(Dist)<-rownames(Data)}
	
	#Calculate distances
	for (i in 1:nrow(Data)) {
		Dist[i,1]<-sqrt(sum((c(Data[i,1], Data[i,2])-c(Data[i,3], Data[i,4]))^2))
	}
	##Recalculate for scaled landmarks
	{if (is.numeric(Scale) & Scale.Column==FALSE) {
		Data<-Data*Scale
		for (i in 1:nrow(Data)) {
			Dist[i,2]<-sqrt(sum((c(Data[i,1], Data[i,2])-c(Data[i,3], Data[i,4]))^2))
		}
	}
	else if (Scale==FALSE & Scale.Column==TRUE) {
		Data<-Data[,1:4]*Data[,5]
		for (i in 1:nrow(Data)) {
			Dist[i,2]<-sqrt(sum((c(Data[i,1], Data[i,2])-c(Data[i,3], Data[i,4]))^2))
		}
	}
	}
	
	#Return result
	return(Dist)
}

#########################################################################
# Function to calculate angles between landmarks                        #
# based on Claude (2008), pp. 50f                                       #
# Necessary input variables:                                            #
#   Data: X and y coordinates of landmarks between which distances...   #
#         should be calculated.                                         #
#         *matrix*                                                      #
# Output data: Data matrix containing the calculated angle for each...  #
#              observation.                                             #
# Input dataset: Matrix with eight or nine columns. Each row...         #
#                corresponds to one observation. Columns 1 to 4 hold... #
#                x and y-coordinates for landmarks 1 and 2 of vector... #
#                1, columns 5 to 8 hold x and y coordinates of...       #
#                landmarks 1 and 2 of vector 2. If...                   #
#                Scale.Column==TRUE a fifth column providing the...     #
#                scaling factor for each observation separately is...   #
#                expected.                                              #
#########################################################################

#Loading packages

LandmarkAngle<-function (Data) {
	#Test data for consistency
	if (ncol(Data)!=8) {stop("Exactly eight columns needed in input data")}
	
	#Transform coordinates to vectors
	Vectors<-array(NA, dim=c(nrow(Data), 2, 2), dimnames=list(rownames(Data), c("x", "y"), c("Vector1", "Vector2")))
	Vectors[,1,1]<-Data[,3]-Data[,1]
	Vectors[,2,1]<-Data[,4]-Data[,2]
	Vectors[,1,2]<-Data[,7]-Data[,5]
	Vectors[,2,2]<-Data[,8]-Data[,6]
	
	#Calculate angles
	Angles<-matrix(NA, nrow(Data), 1)
	colnames(Angles)<-"Angle.deg"
	if (!is.null(rownames(Data))) {rownames(Angles)<-rownames(Data)}
	for (i in 1:(dim(Vectors)[1])) {
		Angles[i,1]<-(acos(sum(Vectors[i,,1]*Vectors[i,,2])/(sqrt(sum(Vectors[i,,1]^2))*sqrt(sum(Vectors[i,,2]^2)))))*(180/pi)
	}
	
	return(Angles)
}

#########################################################################
# Function to calculate confidence ellipse for bivariate data           #
# based on Claude (2008), p. 85                                         #
# Necessary input variables:                                            #
#   Data: Bivariate morphometric measurements for which the bivariate...#
#         confidence ellipse should be calculated.                      #
#         *matrix*                                                      #
#   conf: Confidence level to calculate ellipse for.                    #
#         *numeric (real)*                                              #
#         default=0.95                                                  #
#   NPerim: Number of perimeter points to calculate.                    #
#           *numeric (integer)*                                         #
#           default=50                                                  #
# Output data: Data matrix containing the coordinates of points along...#
#              the confidence ellipse for plotting.                     #
# Input dataset: Matrix with two columns, containing x and y for a...   #
#                bivariate dataset.                                     #
#########################################################################

#Loading packages

BivarEllipse<-function (Data, conf=0.95, NPerim=50) {
	#Test data for consistency
	if (ncol(Data)!=2) {stop("Exactly two columns needed in input data")}
	if (conf>=1 | conf<=0) {stop("Confidence level must vary between 0 and 1")}
	NPerim<-round(NPerim, digits=0)
	
	#Calculate data parameters
	centroid<-apply(Data, 2, mean)
	ang<-seq(0, 2*pi, length=NPerim)
	z<-cbind(cos(ang), sin(ang))
	
	#Calculate confidence ellipse
	radiuscoef<-qnorm((1-conf)/2, lower.tail=F)
	r.xy<-cor(Data[,1], Data[,2])
	M1<-matrix(c(1, 1, -1, 1), 2, 2)
	M2<-matrix(c(var(Data[,1]), var(Data[,2])), 2, 2)
	M3<-matrix(c(1+r.xy, 1-r.xy), 2, 2, byrow=TRUE)
	ellpar<-M1*sqrt(M2*M3/2)
	return(t(centroid + radiuscoef * ellpar %*% t(z)))
}

#########################################################################
# Function to test bivariate data for allometry                         #
# based on Claude (2008), p. 97                                         #
# Necessary packages: smatr                                             #
# Necessary input variables:                                            #
#   Data: Bivariate morphometric measurements which should be tested... #
#         for allometry.                                                #
#         *matrix*                                                      #
#   Log: Should data be log-transformed?                                #
#        NOTE: Log-transformation is necessary for allometry test....   #
#              Only set to FALSE if input are already log-transformed.  #
#        *logical*                                                      #
#        default=TRUE                                                   #
#   Method: Which method (of "MA", "SMA", and "OLS", see the smatr...   #
#           package documentation for details) should be used for...    #
#           regression.                                                 #
#           NOTE: Legendre (sine anno) gives an overview which...       #
#                 methods are applicable under certain circumstances.   #
#           *string*                                                    #
#           default="MA"                                                #
# Output data: List with used method, regression model parameters,...   #
#              slope test results, and result of allometry test.        #
# Input dataset: Matrix with bivariate data (two columns)....           #
#                Observations in rows, x and y in columns.              #
#########################################################################

#Loading packages
require(smatr)

Allometry<-function (Data, Log=TRUE, Method="MA") {
	#Check data for consistency
	if (Method!="MA" & Method!="SMA" & Method!="OLS") {stop("Method must be one of 'MA', 'SMA', or 'OLS'")}
	if (ncol(Data)!=2) {stop("Data must have exactly two columns")}
	
	#Transform data
	if (Log==TRUE) {Data<-log(Data)}
	
	#Calculate regression
	##Calculate model
	RegMod<-line.cis(Data[,2], Data[,1], method=Method)
		
	##Test for isometry
	IsoTest<-unlist(slope.test(Data[,2], Data[,1], test.value=1, method=Method))
	{if (IsoTest["p"]<0.05) {IsoRes<-"Allometry assumed"}
	else {"Isometry assumed"}}
		
	##Plot model
	{if (!is.null(colnames(Data))) {AxesNames<-paste("Log", colnames(Data), sep=" ")}
	else {AxesNames<-c("Log var. 1", "Log var. 2")}}
	plot(Data[,1], Data[,2], pch=16, col="grey50", xlab=AxesNames[1], ylab=AxesNames[2])
	curve(RegMod[2,1]*x+RegMod[1,1], add=TRUE, col="black", lwd=2, lty=1)
	curve(RegMod[2,2]*x+RegMod[1,3], add=TRUE, col="grey50", lwd=1, lty=2)
	curve(RegMod[2,3]*x+RegMod[1,2], add=TRUE, col="grey50", lwd=1, lty=2)
		
	##Output results
	return(list(Regression.Type=Method, Regression.Model=RegMod, Isometry.Test=IsoTest, Isometry.Result=IsoRes))
}

#########################################################################
# Function to compare allometry in bivariate data between groups        #
# based on Claude (2008), pp. 97f                                       #
# Necessary packages: smatr                                             #
# Necessary input variables:                                            #
#   Data: Bivariate morphometric measurements which should be...        #
#         compared concerning their allometry.                          #
#         *matrix*                                                      #
#   Groups: Vector coding the different groups.                         #
#           *factor*                                                    #
#   Log: Should data be log-transformed?                                #
#        NOTE: Log-transformation is necessary for allometry test....   #
#              Only set to FALSE if input are already log-transformed.  #
#        *logical*                                                      #
#        default=TRUE                                                   #
#   Method: Which method (of "MA" and "SMA", see the smatr package...   #
#           documentation for details) should be used for regression.   #
#           NOTE: Legendre (sine anno) gives an overview which...       #
#                 methods are applicable under certain circumstances.   #
#           *string*                                                    #
#           default="MA"                                                #
# Output data: Data matrix containing the calculated angle for each...  #
#              observation.                                             #
# Input dataset: Matrix with bivariate data (two columns)....           #
#                Observations in rows, x and y in columns.              #
#########################################################################

#Loading packages
require(smatr)

Allometry.Comparison<-function (Data, Groups, Log=TRUE, Method="MA") {
	#Check data for consistency
	if (Method!="MA" & Method!="SMA") {stop("Method must be one of 'MA' or 'SMA'")}
	if (ncol(Data)!=2) {stop("Data must have exactly two columns")}
	
	#Transform data
	if (Log==TRUE) {Data<-log(Data)}
	
	#Compare allometry between groups
	AlloComp<-slope.com(Data[,2], Data[,1], groups=Groups)
	
	##Plot model
	{if (!is.null(colnames(Data))) {AxesNames<-paste("Log", colnames(Data), sep=" ")}
	else {AxesNames<-c("Log var. 1", "Log var. 2")}}
	Col<-rainbow(length(unique(Groups)))
	XLIM<-c(max(Data[,1]), min(Data[,1]))
	YLIM<-c(max(Data[,2]), min(Data[,2]))
	for (i in 1:length(unique(Groups))) {
		PlotData<-Data[which(Groups==unique(Groups)[i]),]
		{if (i==1) {plot(PlotData[,1], PlotData[,2], pch=16, col=Col[i], xlab=AxesNames[1], ylab=AxesNames[2], xlim=XLIM, ylim=YLIM)}
		else {points(PlotData[,1], PlotData[,2], pch=16, col=Col[i])}}
		RegMod<-line.cis(PlotData[,2], PlotData[,1], method=Method)
		curve(RegMod[2,1]*x+RegMod[1,1], add=TRUE, col=Col[i], lwd=2, lty=1)
	}
		
	##Output results
	return(AlloComp)
}

#########################################################################
# Function to test multivariate data for allometry                      #
# based on Claude (2008), p. 110                                        #
# Necessary input variables:                                            #
#   Data: Morphometric dataset with several parameters.                 #
#         *matrix*                                                      #
#   Log: Should data be log-transformed?                                #
#        NOTE: The test requires log-transformed data. Set to FALSE...  #
#              only if Data already contains log-transformed data.      #
#        *logical*                                                      #
#        default=TRUE                                                   #
#   Method: Which method should be used to test for allometry. Must...  #
#           be either "Jolicoeur" for a chi-square approach Jolicoeur...#
#           (1963) or "Angle" for a permutation angle comparison.       #
#           *string*, either of "Jolicoeur" or "Angle"                  #
#           default="Jolicoeur"                                         #
# Output data: Results of a chi-square test or permutation test for...  #
#              allometry.                                               #
# Input dataset: Matrix with morphometric data. Observations in rows,...#
#                parameters in columns.                                 #
#########################################################################

#Loading packages

Multivar.Allometry<-function (Data, Log=TRUE, Method="Jolicoeur") {
	#Test data consistency
	if (ncol(Data)>nrow(Data)) {stop("Data must contain more observations than variables")}
	if (Method!="Jolicoeur" & Method!="Angle") {stop("Method must be one of 'Jolicoeur' or 'Angle'")}
	
	#Retrieve data dimension
	n<-dim(Data)[1]
	p<-dim(Data)[2]
	
	#Transform data
	if (Log==TRUE) {Data<-log(Data)}
	
	#Define function for angle comparison
	coli<-function (ev1, ev2, nperm=10000) {
		Dist<-numeric(nperm)
		n<-length(ev1)
		Angle<-function (v1, v2) {90-abs((180*(acos(sum(v1*v2)/(sqrt(sum(v1^2))*sqrt(sum(v2^2)))))/pi)-90)}
		for (i in 1:nperm) {
			X1<-runif(n, -1, 1)
			X2<-runif(n, -1, 1)
			Dist[i]<-Angle(X1, X2)
		}
		zobs<-Angle(ev1, ev2)
		pv<-length(Dist[Dist<zobs])/nperm
		list(angle=zobs, p=pv)
	}
	
	{if (Method=="Jolicoeur") {
		#Compute variance-covariance-matrix
		S<-var(Data)
		V1<-rep(sqrt(1/p), p)
	
		#Calculate first singular value
		L1<-svd(S)$d[1]
		chiobs<-(n-1)*(L1*t(V1)%*%ginv(S)%*%V1+(1/L1)*t(V1)%*%S%*%V1-2)
	
		#Return results
		return(unlist(list(Chisq=chiobs, p=pchisq(chiobs,p-1,lower.tail=F))))
	}
	else {
		return(unlist(coli(prcomp(Data)[[2]][,1], rep(sqrt(1/p), p))))
	}
	}
}

#########################################################################
# Plots the whole dataset on the PCA axes defined by one group in...    #
#   the dataset                                                         #
# based on Claude (2008), p. 110                                        #
# Required functions: BivarEllipse                                      #
# Necessary input variables:                                            #
#   Data: Morphometric dataset with several parameters.                 #
#         *matrix*                                                      #
#   Groups: Vector, encoding groups of Data.                            #
#           *factor* of length nrow(Data)                               #
#   Proj.group: Name (according to Groups) of the group into which...   #
#               PCA space the data should be projected.                 #
#               *string*                                                #
# Output data: Partial PCA plot                                         #
# Input dataset: Matrix with morphometric data in more than one group...#
#                Observations in rows, parameters in columns.           #
#########################################################################

#Loading packages

Partial.PCA<-function (Data, Groups, Proj.group) {
	#Test data consistency
	if (length(Groups)!=nrow(Data)) {stop("Groups vector must be of same length as Data")}
	if (!any(Groups==Proj.group)) {stop("Proj.group must be contained in Groups")}
	
	#Perform PCA
	Groups<-as.factor(Groups)
	Proj.data<-Data[Groups==Proj.group,]
	if (nrow(Proj.data)<=ncol(Proj.data)) {stop("Proj.group must have at least as many observations as variables")}
	PCA<-prcomp(Proj.data)
	Proj<-Data%*%PCA[[2]]
	
	#Plot PCA
	ColSet<-rainbow(length(levels(Groups)))
	Cols<-ColSet[match(Groups, levels(Groups))]
	plot(Proj[,2:3], col=Cols, pch=16, xlab="PC 2", ylab="PC 3", main=paste("Partial PCA on Group ", Proj.group, " space", sep=""))
	for (i in 1:(length(levels(Groups)))) {
		lines(BivarEllipse(Proj[which(Groups==levels(Groups)[i]),2:3]), col=ColSet[i])
	}
	legend("topleft", pch=16, col=ColSet, legend=paste("Group ", levels(Groups), sep=""), cex=0.7)
}

#########################################################################
# Function to calculate Mosimann shape vectors                          #
# based on Claude (2008), pp. 99f                                       #
# Necessary input variables:                                            #
#   Data: Morphometric dataset with several parameters.                 #
#         *matrix*                                                      #
#   Type: Type of combining measurements:                               #
#         "Arithmetic": arithmetic mean                                 #
#         "Geometric": geometric mean                                   #
#         "SumSquares": square-root of sum of squared measurements      #
#         *string*                                                      #
#         default="Geometric"                                           #
#   Groups: Optional vector to code groups in data.                     #
#           *factor*                                                    #
#           default=NULL                                                #
# Output data: Data matrix containing scaled parameters.                #
# Input dataset: Matrix with morphometric data. Observations in rows,...#
#                parameters in columns.                                 #
#########################################################################

#Loading packages

Mosimann<-function (Data, Type="Geometric", Groups=NULL) {
	#Check data for consistency
	if (Type!="Arithmetic" & Type!="Geometric" & Type!="SumSquares") {stop("Type must be one of 'Arithmetic', 'Geometric', or 'SumSquares'")}
	if (!is.null(Groups)) {
		Levels<-unique(Groups)
		Sym<-match(Groups, Levels)
	}
	
	#Calculate mean
	{if (Type=="Arithmetic") {Scale<-apply(Data, 1, mean)}
	else if (Type=="Geometric") {Scale<-apply(Data, 1, prod)^(1/ncol(Data))}
	else {Scale<-sqrt(apply(Data, 1, sum)^2)}
	}
	
	#Calculate Mosimann vector
	Shape<-as.matrix(Data/Scale)
	
	#Plot Mosiman vectors
	{if (is.null(Groups)) {pairs(log(Shape), pch=16)}
	else {
		pairs(log(Shape), pch=Sym)
		par(xpd=TRUE)
		legend(0, 1, pch=unique(Sym), legend=Levels, bg="white", cex=0.7)
	}}
	
	#Return scaled size parameters
	return(Shape)
}

#########################################################################
# Function to calculate residual-variance shape vectors                 #
# based on Claude (2008), pp. 101f                                      #
# Necessary input variables:                                            #
#   Data: Morphometric dataset with several parameters.                 #
#         *matrix*                                                      #
#   Type: Type of combining measurements:                               #
#         "Arithmetic": arithmetic mean                                 #
#         "Geometric": geometric mean                                   #
#         "SumSquares": square-root of sum of squared measurements      #
#         *string*                                                      #
#         default="Geometric"                                           #
#   Allometric: Should the model be allowed to allocate some...         #
#               allometric variation into size?                         #
#               *logical*                                               #
#               default=FALSE                                           #
#   Groups: Optional vector to code groups in data.                     #
#           *factor*                                                    #
#           default=NULL                                                #
# Output data: Data matrix containing scaled parameters and fraction... #
#              of variance explained by size. If Allometry==TRUE...     #
#              additionally the data are checked for contained...       #
#              allometric growth.                                       #
# Input dataset: Matrix with morphometric data. Observations in rows,...#
#                parameters in columns.                                 #
#########################################################################

#Loading packages

ResVar.Shape<-function (Data, Type="Geometric", Allometric=FALSE, Groups=NULL) {
	#Check data for consistency
	if (Type!="Arithmetic" & Type!="Geometric" & Type!="SumSquares") {stop("Type must be one of 'Arithmetic', 'Geometric', or 'SumSquares'")}
	if (!is.null(Groups)) {
		Levels<-unique(Groups)
		Sym<-match(Groups, Levels)
	}
	
	#Calculate mean
	{if (Type=="Arithmetic") {Scale<-apply(Data, 1, mean)}
	else if (Type=="Geometric") {Scale<-apply(Data, 1, prod)^(1/ncol(Data))}
	else {Scale<-sqrt(apply(Data, 1, sum)^2)}
	}
	
	#Calculate residual-variance shape vector
	{if (Allometric==TRUE) {Shape<-lm(as.matrix(Data)~Scale)}
	else {Shape<-lm(as.matrix(Data)~Scale-1)}}
	
	#Plot shape vectors
	{if (is.null(Groups)) {pairs(Shape$residuals, pch=16)}
	else {
		pairs(Shape$residuals, pch=Sym)
		par(xpd=TRUE)
		legend(0, 1, pch=unique(Sym), legend=Levels, bg="white", cex=0.7)
	}}
	
	#Test for allometry
	##Define function
	if (Allometric==TRUE) {win.graph(10, 10, 12)}
	coli<-function (ev1, ev2, nperm=10000) {
		Dist<-numeric(nperm)
		n<-length(ev1)
		Angle<-function (v1, v2) {90-abs((180*(acos(sum(v1*v2)/(sqrt(sum(v1^2))*sqrt(sum(v2^2)))))/pi)-90)}
		for (i in 1:nperm) {
			X1<-runif(n, -1, 1)
			X2<-runif(n, -1, 1)
			Dist[i]<-Angle(X1, X2)
		}
		zobs<-Angle(ev1, ev2)
		pv<-length(Dist[Dist<zobs])/nperm
		hist(Dist, breaks=50, main="Distribution of angles between 2 random vectors", xlab="Z statistic", ylab="# of random vect.", sub=paste("Actual z-obs = ",round(zobs, 5),": p < ", round(pv, digits=4), sep=""), col="grey50")
		abline(v=zobs, col="blue", lwd=2)
		list(angle=zobs, p=pv)
	}
	##Perform test
	if (Allometric==TRUE) {Colinear<-unlist(coli(Shape$coefficients[2,], rep(1, ncol(Shape$coefficients))))}
	
	#Calculate percentage of variance explained
	Variance<-sum(diag(var(Shape$fitted.values)))/sum(diag(var(as.matrix(Data))))
	
	#Return scaled size parameters
	{if (Allometric==TRUE) {return(list(Shape.Vector=Shape$fitted.values, Allometry=Colinear, Variation.Explained=Variance))}
	else {return(list(Shape.Vector=Shape$fitted.values, Variation.Explained=Variance))}}
}

#########################################################################
# Form and shape LDA/CVA                                                #
# based on Claude (2008), p. 112ff                                      #
# Required function: Mosimann                                           #
# Required packages: MASS                                               #
# Necessary input variables:                                            #
#   Data: Morphometric dataset with several parameters.                 #
#         *matrix*                                                      #
#   Groups: Vector, encoding groups of Data.                            #
#           *vector* of length nrow(Data)                               #
#   Type: Type of CVA, either "Form" for raw values or "Shape" for...   #
#         Mosimann shape vectors.                                       #
#         *string*                                                      #
#         default="Form"                                                #
# Output data: CVA plot and statistics for significance of CVA.         #
# Input dataset: Matrix with morphometric data in more than one group...#
#                Observations in rows, parameters in columns.           #
#########################################################################

#Loading packages
require(MASS)

Morpho.CVA<-function (Data, Groups, Type="Form") {
	#Test data consistency
	if (Type!="Form" & Type!="Shape") {stop("Type must be either of 'Form' or 'Shape'")}
	
	#Perform CVA
	Groups<-as.factor(Groups)
	if (Type=="Shape") {Data<-Mosimann(Data)}
	CVA<-lda(Data, Groups)
	MANOVA<-manova(Data~Groups)
	Stat<-(summary(MANOVA, test="Hotelling"))
	Proj<-Data%*%CVA$scaling
	{if (length(levels(Groups))==2) {
		YLIM<-c(0, max(c(max(hist(Proj[Groups==levels(Groups)[1],1], plot=FALSE)$density), max(hist(Proj[Groups==levels(Groups)[2],1], plot=FALSE)$density))))
		hist(Proj[Groups==levels(Groups)[1],1], xlim=c(min(Proj[,1]), max(Proj[,1])), ylim=YLIM, xlab="LD 1", main="Discriminant scores", freq=FALSE, col="red")
		hist(Proj[Groups==levels(Groups)[2],1], xlim=c(min(Proj[,1]), max(Proj[,1])), freq=FALSE, col=rgb(0, 1, 0, 0.5), add=TRUE)
		legend("topright", legend=c(paste("Group", levels(Groups)[1], sep=" "), paste("Group", levels(Groups)[2], sep=" ")), fill=c("red", rgb(0, 1, 0, 0.5)), bg="white", cex=0.7)
	}
	else {
		ColSet<-rainbow(length(levels(Groups)))
		Cols<-ColSet[match(Groups, levels(Groups))]
		{if (Type=="Form") {plot(Proj[,1:2], col=Cols, pch=16, xlab="FD 1", ylab="FD 2", main="Form FDA")}
		else {plot(Proj[,1:2], col=Cols, pch=16, xlab="FD 1", ylab="FD 2", main="Shape FDA")}
		}
		legend("topleft", pch=16, col=ColSet, legend=paste("Group ", levels(Groups), sep=""), cex=0.7)
	}
	}
	return(Stat)
}

#########################################################################
# Burnaby's nuisance correction                                         #
# based on Claude (2008), p. 118f                                       #
# Necessary input variables:                                            #
#   Data: Morphometric dataset with several parameters.                 #
#         *matrix*                                                      #
#   Groups: Vector, encoding groups of Data.                            #
#           *factor* of length nrow(Data)                               #
#   Nuisance: Discriminant axis that corresponds with nuisance factor...#
#             (e.g. size, environment).                                 #
#             *numeric (integer)*                                       #
#             default=1                                                 #
#   Log: Should data be log-transformed.                                #
#        NOTE: The method requires log-transformed data. Set to FALSE...#
#              only if Data already contains log-transformed data.      #
#        *logical*                                                      #
#        default=TRUE                                                   #
# Output data: Test statistic for difference between groups other...    #
#              than caused by the nuisance factor.                      #
# Input dataset: Matrix with morphometric data in more than one group...#
#                Observations in rows, parameters in columns.           #
#########################################################################

#Loading packages

Burnaby<-function (Data, Groups, Nuisance=1, Log=TRUE) {
	#Test data consistency
	if (length(Groups)!=nrow(Data)) {stop("Groups vector must be of same length as Data")}
	if (Nuisance>ncol(Data)) {stop(paste("With ", ncol(Data), " variables there cannot be ", Nuisance, " principal axes", sep=""))}
	
	#Transform data
	Groups<-as.factor(Groups)
	if (Log==TRUE) {Data=log(Data)}
	if ((ncol(Data)-length(levels(Groups)))<1) {stop(paste("With ", ncol(Data), " and ", length(levels(Groups)), " the approach would yield less than one in-fine dimension", sep=""))}
	
	#Prepare data
	G<-prcomp(Data[which(Groups==levels(Groups)[1]),])[[2]][,Nuisance]
	for (i in 2:(length(levels(Groups)))) {
		G<-cbind(G, prcomp(Data[which(Groups==levels(Groups)[i]),])[[2]][,Nuisance])
	}
	I<-diag(1, ncol(Data))
	ortho<-I-G%*%ginv(t(G)%*%G)%*%t(G)
	
	#Calculate results
	newdata<-Data%*%prcomp(Data%*%ortho)[[2]][,1:(ncol(Data)-ncol(G))]
	{if (ncol(newdata)==1) {Stat<-summary(aov(newdata~Groups))}
	else {Stat<-summary(manova(newdata~Groups, test="Hotelling"))}
	}
	return(Stat)
}

#########################################################################
# Clustering of morphometric data                                       #
# Required packages: MASS, shapes, ape, cluster                         #
# based on Claude (2008), pp. 122ff                                     #
# Necessary input variables:                                            #
#   Data: Bi- or multivariate morphometric data.                        #
#         *matrix*                                                      #
#   Groups: Coding of groups in Data.                                   #
#           *factor*                                                    #
#   K: Number of clusters for Kmeans clustering.                        #
#      *numeric (integer)*                                              #
#      default=3                                                        #
# Output data: Plots showing UPGMA and complete clustering, elbow-...   #
#              plot, and partitional clsutering plot.                   #
# Input dataset: Matrix with morphometric data in more than one group...#
#                Observations in rows, parameters in columns.           #
#########################################################################

#Load packages
require(MASS)
require(shapes)
require(ape)
require(cluster)

Clustering<-function (Data, Groups, K=3) {
	#Test data consistency
	K<-round(K, digits=0)
	if (length(Groups)!=dim(Data)[1]){stop("Vector Groups must be of same length as Data")}
	
	#Transform groups
	Groups<-as.factor(Groups)
	
	#Plot results
	win.graph(16, 14, 10)
	par(mar=c(4, 4, 3, 0), oma=c(0, 0, 0, 0))
	layout(matrix(1:4, 2, 2))
	##Define colour-coding
	ColSet<-rainbow(length(levels(Groups)))
	Cols<-ColSet[match(Groups, levels(Groups))]
	##Prepare and plot clusters
	Dist.Mat<-dist(Data, method="euclidean")
	Clust<-hclust(Dist.Mat, method="average")
	Clust<-as.phylo(Clust)
	Clust$tip.label<-paste(as.character(Groups), Clust$tip.label, sep="-")
	plot(Clust, type="fan", tip.color=Cols, main="UPGMA")
	Clust<-hclust(Dist.Mat, method="complete")
	Clust<-as.phylo(Clust)
	Clust$tip.label<-paste(as.character(Groups), Clust$tip.label, sep="-")
	plot(Clust, type="fan", tip.color=Cols, main="Complete")
	
	#Perform partitional clustering
	d.f<-nrow(Data)-1
	SStot<-sum(diag(var(Data)))*d.f
	SSratio<-0
	for (i in 2:(round(dim(Data)[1]/3, digits=0))) {
		mod<-kmeans(Data, i)
		SSratio[i]<-(SStot-sum(mod$withinss))/SStot
	}
	plot(1:(round(dim(Data)[1]/3, digits=0)), SSratio, xlab="Number of clusters", ylab="Explained variance (fraction)", pch=16, type="b")
	KMeans<-pam(dist(Data), k=K)
	Class<-length(which(KMeans$clustering==as.numeric(Groups)))/length(Groups)
	plot(KMeans, which.plot=1, shade=TRUE, col.p=Cols, col.txt=Cols, col.clus="grey90", main=paste("Correct classification: ", round(Class*100, digits=2), "%", sep=""))
}

#########################################################################
# Two-block partial-least-squares approach for relationship between...  #
#   two sets of variables                                               #
# based on Claude (2008), pp. 124ff                                     #
# Necessary input variables:                                            #
#   Shape.Data: Bi- or multivariate morphometric data.                  #
#               *matrix*                                                #
#   Parameter.Matrix: Matrix containing any predictor variable. Can,... #
#                     but does not have to be, morphometric data. The...#
#                     first column is expected to contain the...        #
#                     sample-coding.                                    #
#                     *matrix*                                          #
# Output data: Test statistics and model for partial least squares...   #
#              regression between the two datasets.                     #
# Input dataset: Matrix with morphometric data. Observations in rows,...#
#                parameters in columns.                                 #
#########################################################################

#Load packages
require(pls)

TBPLS<-function (Shape.Data, Parameter.Matrix) {
	#Prepare shape matrix, parameter matrix, and sample coding
	Shapes<-colnames(Shape.Data)
	Shape.Data<-as.matrix(Shape.Data)
	NA.rows<-complete.cases(Parameter.Matrix[,-1])
	Samples<-Parameter.Matrix[NA.rows,1]
	{if (is.null(colnames(Parameter.Matrix))) {Params<-paste("Param.", 1:(ncol(Parameter.Matrix)-1))}
	else {Params<-colnames(Parameter.Matrix)[-1]}}
	Parameter.Matrix<-as.matrix(Parameter.Matrix[,-1])
	colnames(Parameter.Matrix)<-Params
	Param.N<-ncol(Parameter.Matrix)
	
	#Test data for consistency
	if (nrow(Shape.Data)!=nrow(Parameter.Matrix)) {stop("Shape and parameter matrix do not have the same length!")}
	
	#Define function to calculate trace of square matrix
	Trace<-function (Mat) {
		Vals<-vector(mode="numeric", length=ncol(Mat))
		for (i in 1:(ncol(Mat))) {
			Vals[i]<-Mat[i,i]
		}
		Tr<-sum(Vals)
		return(Tr)
	}
	
	#Combine data into data frame
	PLS.Data<-data.frame(matrix(integer(0), nrow=nrow(Shape.Data)))
	PLS.Data$Shape<-Shape.Data
	PLS.Data<-cbind(PLS.Data, Parameter.Matrix)
	
	#Calculate significance of correlation between blocks
	R1<-var(Shape.Data, na.rm=TRUE)
	R2<-var(Parameter.Matrix, na.rm=TRUE)
	R12<-cov(Shape.Data, Parameter.Matrix, use="complete.obs")
	##Calculate Escoufiers coefficient
	RV<-Trace(R12%*%t(R12))/(sqrt(Trace(R1%*%t(R1))*Trace(R2%*%t(R2))))
	##Perform randomization test
	RV.Rand<-vector(mode="numeric", length=1000)
	for (i in 1:1000) {
		Rand.Vec<-sample(nrow(Shape.Data), replace=FALSE)
		Rand.Shape<-Shape.Data[Rand.Vec,]
		R1<-var(Rand.Shape, na.rm=TRUE)
		R12<-cov(Rand.Shape, Parameter.Matrix, use="complete.obs")
		RV.Rand[i]<-Trace(R12%*%t(R12))/(sqrt(Trace(R1%*%t(R1))*Trace(R2%*%t(R2))))
	}
	Large<-length(which(abs(RV.Rand)>=abs(RV)))
	Small<-length(which(abs(RV.Rand)<abs(RV)))
	P<-Large/(Large+Small)
	
	#Calculate partial least squares
	Indep<-paste(colnames(Parameter.Matrix), collapse="+")
	PLS<-plsr(as.formula(paste("Shape ~ ", Indep)), data=PLS.Data, na.action=na.exclude, validation="LOO", jackknife=TRUE)
	
	#Return results
	return(list("Significance"=P, "Morphology.Values"=Shapes, "Number.Parameters"=Param.N, "Parameter.Names"=Params, "Sample.Names"=Samples, "Model"=PLS))
}

#########################################################################
# Two-block partial least squares visualization                         #
# Required packages: pls                                                #
# Necessary input variables:                                            #
#   Model: Output as produced by TBPLS.                                 #
#   PC1: First principal component to use.                              #
#        *numeric (integer)*                                            #
#        default=1                                                      #
#   PC2: Second principal component to use.                             #
#        *numeric (integer)*                                            #
#        default=2                                                      #
#   Type: Which type of diagnostic plots should be produced. Either...  #
#         "Scoreplot" (for scores- and loadings- plots), "Biplot"...    #
#         (for residual variations and biplot), "T.U.Plot" (for...      #
#         correlation plot between T and U scores).                     #
#         *string*                                                      #
#         default=NULL                                                  #
#   Multiplic: Optional multiplicator for loadings vectors. Only...     #
#              meaningful for Type=="Biplot".                           #
#              *numeric (real)*                                         #
#              default=NULL                                             #
# Output data: A series of diagnostic plots for the PLS.                #
# Input dataset: PLS Model list as produced by PLS.Shape.               #
#########################################################################

#Load packages
require(pls)

TBPLS.Plot<-function (Model, PC1=1, PC2=2, Type=NULL, Multiplic=NULL) {
	#Test data consistency
	PC1<-round(PC1, digits=0)
	PC2<-round(PC2, digits=0)
	if (is.null(Type)) {stop("'Type' must be provided!")}
	if (Type!="Scoreplot" & Type!="Biplot" & Type!="T.U.Plot") {stop("Type must be either of 'Scoreplot', 'Biplot', or 'T.U.Plot'!")}

	#Prepare data
	Mod<-Model$Model
	
	#Gather data
	PLS.scoresX<-scores(Mod)
	PLS.scoresY<-Yscores(Mod)
	PLS.loadingsX<-loadings(Mod)
	PLS.loadingsY<-Yloadings(Mod)
	PLS.residuals<-explvar(Mod)
	
	#PLS1 visualization
	if (Model$Number.Parameters==1) {
		if (PC1!=1) {PC1<-1; warning("PLS1 model, 'PC1' set to 1.")}
		#Scoreplot/loadingsplot
		if (Type=="Scoreplot") {
			win.graph(19, 10, 10)
			layout(matrix(c(1, 2), 1, 2))
			plot(max(PLS.scoresX[,PC1])*1000, 1000, xlim=range(PLS.scoresX[,PC1]), xlab=paste("PC ", PC1, sep=""), ylim=c(-1, 1), ylab="Meaningless", main="Scores")
			text(PLS.scoresX[,PC1], 0, labels=Model$Sample.Names)
			plot(max(c(PLS.loadingsX[,PC1], PLS.loadingsY[,PC1]))*1000, 1000, xlim=range(c(PLS.loadingsX[,PC1], PLS.loadingsY[,PC1])), xlab=paste("PC ", PC1, sep=""), ylim=c(1, 2), ylab="Meaningless", main="Loadings")
			text(PLS.loadingsX[,PC1], 1, labels=rownames(PLS.loadingsX), col="red")
			text(PLS.loadingsY[,PC1], 2, labels=Model$Morphology.Values, col="blue")
		}
	
		#Residual variation/biplot
		if (Type=="Biplot") {
			hist(PLS.scoresX[,PC1], xlab=paste("PC ", PC1, " (", round(PLS.residuals[PC1], digits=1), "%)", sep=""), col="grey50", main="Biplot", sub=paste("Positive with ", rownames(PLS.loadingsX), sep=""))
		}
	
		#T-U-scoreplots
		if (Type=="T.U.Plot") {
			TU1<-lm(PLS.scoresY[,PC1]~PLS.scoresX[,PC1])
			R<-summary(TU1)$adj.r.squared
			P<-summary(TU1)$coefficients[2,"Pr(>|t|)"]
			plot(PLS.scoresX[,PC1], PLS.scoresY[,PC1], xlab="T scores", ylab="U scores", sub=paste("R2 = ", R, ", p = ", P, sep=""), main=paste("PC ", PC1, sep=""))
			curve(TU1$coefficients[2]*x+TU1$coefficients[1], add=TRUE, lwd=2)
		}
	}
	
	#PLS2 visualization
	if (Model$Number.Parameters>1) {
		#Scoreplot/loadingsplot
		if (Type=="Scoreplot") {
			win.graph(19, 10, 10)
			layout(matrix(c(1, 2), 1, 2))
			plot(max(PLS.scoresX[,PC1])*1000, max(PLS.scoresX[,PC2])*1000, xlim=range(PLS.scoresX[,PC1]), xlab=paste("PC ", PC1, sep=""), ylim=range(PLS.scoresX[,PC2]), ylab=paste("PC ", PC2, sep=""), main="Scores")
			text(PLS.scoresX[,PC1], PLS.scoresX[,PC2], labels=Model$Sample.Names)
			plot(max(c(PLS.loadingsX[,PC1], PLS.loadingsY[,PC1]))*1000, max(c(PLS.loadingsX[,PC2], PLS.loadingsY[,PC2]))*1000, xlim=range(c(PLS.loadingsX[,PC1], PLS.loadingsY[,PC1])), xlab=paste("PC ", PC1, sep=""), ylim=range(c(PLS.loadingsX[,PC2], PLS.loadingsY[,PC2])), ylab=paste("PC ", PC2, sep=""), main="Loadings")
			text(PLS.loadingsX[,PC1], PLS.loadingsX[,PC2], labels=rownames(PLS.loadingsX), col="red")
			text(PLS.loadingsY[,PC1], PLS.loadingsY[,PC2], labels=Model$Morphology.Values, col="blue")
		}
	
		#Residual variation/biplot
		if (Type=="Biplot") {
			win.graph(19, 10, 10)
			layout(matrix(c(1, 2), 1, 2))
			barplot(PLS.residuals, col="blue", ylim=c(0, 100), ylab="Variation explained (%)", main="Variance explained")
			plot(PLS.scoresX[,PC1], PLS.scoresX[,PC2], pch=16, xlab=paste("PC ", PC1, " (", round(PLS.residuals[PC1], digits=1), "%)", sep=""), ylab=paste("PC ", PC2, " (", round(PLS.residuals[PC2], digits=1), "%)", sep=""), main="Biplot")
			{if (is.null(Multiplic)) {text(PLS.loadingsX[,PC1], PLS.loadingsX[,PC2], labels=rownames(PLS.loadingsX), col="red")}
			else {text(PLS.loadingsX[,PC1]*Multiplic, PLS.loadingsX[,PC2]*Multiplic, labels=rownames(PLS.loadingsX), col="red")}
			}
		}
	
		#T-U-scoreplots
		if (Type=="T.U.Plot") {
			win.graph(19, 10, 10)
			layout(matrix(c(1, 2), 1, 2))
			TU1<-lm(PLS.scoresY[,PC1]~PLS.scoresX[,PC1])
			R<-summary(TU1)$adj.r.squared
			P<-summary(TU1)$coefficients[2,"Pr(>|t|)"]
			plot(PLS.scoresX[,PC1], PLS.scoresY[,PC1], xlab="T scores", ylab="U scores", sub=paste("R2 = ", R, ", p = ", P, sep=""), main=paste("PC ", PC1, sep=""))
			curve(TU1$coefficients[2]*x+TU1$coefficients[1], add=TRUE, lwd=2)
			TU2<-lm(PLS.scoresY[,PC2]~PLS.scoresX[,PC2])
			R<-summary(TU2)$adj.r.squared
			P<-summary(TU2)$coefficients[2,"Pr(>|t|)"]
			plot(PLS.scoresX[,PC2], PLS.scoresY[,PC2], xlab="T scores", ylab="U scores", sub=paste("R2 = ", R, ", p = ", P, sep=""), main=paste("PC ", PC2, sep=""))
			curve(TU2$coefficients[2]*x+TU2$coefficients[1], add=TRUE, lwd=2)
		}
	}
}

#--------------------------------------------

#Examples
##Landmark distances
#LandmarkDist(LM)
#LandmarkDist(LM, Scale=0.1)
#LandmarkDist(LM2, Scale.Column=TRUE)

##Landmark angles
#LandmarkAngle(LM3)

##Bivariate confidence ellipse
#ELLI1<-BivarEllipse(Bivar[which(Bivar[,1]==1),2:3], conf=0.66)
#ELLI2<-BivarEllipse(Bivar[which(Bivar[,1]==2),2:3], conf=0.66)
#plot(Bivar[which(Bivar[,1]==1),2], Bivar[which(Bivar[,1]==1),3], col="blue", pch=16)
#points(Bivar[which(Bivar[,1]==2),2], Bivar[which(Bivar[,1]==2),3], col="red", pch=16)
#lines(ELLI1, col="blue")
#lines(ELLI2, col="red")

##Allometry test
#D1<-Bivar[which(Bivar[,"Species"]==1),-1]
#Allometry(D1, Method="OLS")
#Allometry.Comparison(Data=Bivar[,-1], Groups=Bivar[,1])
#Multivar.Allometry(Multivar[which(Multivar[,1]==1),-1])
#Multivar.Allometry(Multivar[which(Multivar[,1]==1),-1], Method="Angle")

##Shape vectors
#Mosimann(Bivar[which(Bivar[,1]==1),-1])
#Mosimann(Bivar[which(Bivar[,1]==1),-1], Groups=Bivar[,1])
#ResVar.Shape(Bivar[which(Bivar[,1]==1),-1])
#ResVar.Shape(Bivar[which(Bivar[,1]==1),-1], Allometric=TRUE)

##Shape comparisons
#Partial.PCA(Multivar[,-1], Groups=Multivar[,1], Proj.group=1)
#Partial.PCA(Multivar[,-1], Groups=Multivar[,1], Proj.group=2)
#Morpho.CVA(Multivar[,-1], Groups=Multivar[,1])
#Morpho.CVA(Multivar[,-1], Groups=Multivar[,1], Type="Shape")
#Morpho.CVA(rbind(Multivar[,-1], cbind(Multivar[1:10,2]+runif(10, min=-0.5, max=2), Multivar[1:10,3]+runif(10, min=-0.06, max=0.12), Multivar[1:10,4]+runif(10, min=-0.003, max=0.010))), Groups=c(Multivar[,1], rep(3, 10)))
#Burnaby(Multivar[,-1], Groups=Multivar[,1])
#Clustering(Multivar[,-1], Groups=Multivar[,1])
#TBPLS(cbind(Multivar[,"Body.Length"], Multivar[,"Eye.Size"]), Predictor=as.matrix(Multivar[,"Brain.Size"]))

#--------------------------------------------
#--------------------------------------------
#Version History
#1.0	Finished Program
#1.1	Included test for contained allometry in ResVar.Shape
#1.2	Several new functions included
#1.2.1	Coding of groups in Mosimann and ResVar.Shape supported
#1.2.2	Test statistic changed to Hotelling for CVA
#1.3	Rewrote function TBPLS analogue to geometric morphometrics version, added function TBPLS.Plot
#--------------------------------------------
#--------------------------------------------
