#Spiral Growth Analysis
#Measuring whether or not a spiral growth process follows the logarithmic spiral.
#Input data set:
#	-Images in the ppm format

#Author: Manuel Weinkauf (Manuel.Weinkauf@unige.ch)
#Version: 1.4.2
#Date: 19 March 2020

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.#
#To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/.                   #
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#**************************************************************************************
#Setting working directory
#setwd("C:/R_TestData/SGA")

#########################################################################
# Analysing data of spiral outlines                                     #
# Necessary packages: lmodel2                                           #
# Necessary input variables:                                            #
#   Input: Dataset to be analysed.                                      #
#          *list*                                                       #
#   T: Name of the column containing the radial distance value.         #
#      *string*                                                         #
#      default="t"                                                      #
#   Theta: Name of the column containing the angle value.               #
#          *string*                                                     #
#          default="theta"                                              #
#   Raw: Should the raw data for the regression be exported as well?    #
#        *logical*                                                      #
#        default=TRUE                                                   #
# Input dataset: A set of data points along the spiral, given as...     #
#                radial distance t and  angle theta in radians....      #
#                Values are supposed to be normalized for radius one... #
#                and rotated such that theta[1] is zero! The object...  #
#                "Input" is supposed to be a list as produced by...     #
#                Read.Spiral.                                           #
# Output dataset: Matrix with dip, intercept, and R2 for model II...    #
#                 linear regression (ranged major axis) of log-...      #
#                 transformed distances t dependent on theta.           #
#########################################################################

#Loading packages
require(lmodel2)

SpiralAnalysis<-function(Input, T="t", Theta="theta", Raw=TRUE) {
	#Test data consistency
	Input.Colnames<-lapply(Input, colnames)
	Col.checker<-list()
	Col.checker$T<-Col.checker$Theta<-vector(mode="logical", length=length(Input))
	for (i in 1:length(Input)) {
		if (T%in%Input.Colnames[[i]]) {Col.checker$T[i]<-TRUE}
		if (Theta%in%Input.Colnames[[i]]) {Col.checker$Theta[i]<-TRUE}
	}
	if (!all(Col.checker$T==TRUE)) {stop("No data for radial distance in all datasets!")}
	if (!all(Col.checker$Theta==TRUE)) {stop("No data for angle in all datasets!")}

	#Setting up results matrix
	Res<-matrix(NA, length(Input), 3)
	colnames(Res)<-c("Dip", "Intercept", "Adjusted R squared")
	rownames(Res)<-names(Input)
	for (k in 1:length(Input)) {
		#Choosing data subset
		Data<-Input[[k]]
		Data<-cbind(Data[,T], Data[,Theta])
	
		#Determining whether spiral is left or right winding
		{if (Data[3,2]>Data[2,2]) {Left<-FALSE}
		else {Left<-TRUE}}
		
		#Chaining radial values to continuous chain, i.e. eliminating effect of crossing the 2pi boundary
		{if (Left==FALSE) {
			tlist<-vector()
			for (i in 2:(dim(Data)[1])) {
				if (Data[i,2]<Data[(i-1),2]) {tlist<-append(tlist, i)}
			}
			tlist<-append(tlist, (dim(Data)[1]+1))
			tc<-1
			for (i in 1:((length(tlist))-1)) {
				for (j in (tlist[i]):((tlist[i+1])-1)) {
					Data[j,2]<-Data[j,2]+(tc*(pi*2))
				}
				tc<-tc+1
			}
		}
		else {
			tlist<-2
			for (i in 3:(dim(Data)[1])) {
				if (Data[i,2]>Data[(i-1),2]) {tlist<-append(tlist, i)}
			}
			tlist<-append(tlist,dim(Data)[1]+1)
			tc<-1
			for (i in 1:((length(tlist))-1)) {
				for (j in (tlist[i]):((tlist[i+1])-1)) {
					Data[j,2]<-(tc*(pi*2))-Data[j,2]
				}
				tc<-tc+1
			}
		}
		}
		
		#Log transforming t for applicability of linear regression
		Data[,1]<-log(Data[,1])
		#Storing raw data
		if (Raw==TRUE) {
			Res.Temp<-as.data.frame(matrix(NA, nrow(Data), 3))
			Res.Temp[,1]<-rep(names(Input)[k], nrow(Data))
			Res.Temp[,2:3]<-Data
			{if (k==1) {Res2<-Res.Temp}
			else {Res2<-rbind(Res2, Res.Temp)}}
		}
		#Calculating model II linear regression
		Fit<-lmodel2(Data[,1]~Data[,2], range.x="relative", range.y="interval")
		Res[k,1]<-Fit$regression.results[which(Fit$regression.result[,"Method"]=="RMA"),"Slope"]
		Res[k,2]<-Fit$regression.results[which(Fit$regression.result[,"Method"]=="RMA"),"Intercept"]
		##Calculating R2 of RMA
		y.mean<-mean(Data[,1])
		SS.tot<-sum((Data[,1]-y.mean)^2)
		SS.res<-sum((Data[,1]-(Data[,2]*Res[k,1]+Res[k,2]))^2)
		Res[k,3]<-1-(SS.res/SS.tot)
	}
	
	#Return results
	if (Raw==TRUE) {
		colnames(Res2)<-c("ID", "log(T)", "Theta")
		Res<-list(Regression=Res, Raw.values=Res2)
	}
	return(Res)
}

#--------------------------------------------

#Example
#ImageConversion("Scit",1,2,".tif",".ppm")
#SpiralExtraction("Scit",1,2)
#SpiralAnalysis("Scit",1,2)

#--------------------------------------------
#--------------------------------------------
#Version History
#1.0	Finished Program
#1.1	Adds option in SpiralExtraction to skip size normalization
#1.2	Changes OLS to RMA in SpiralAnalysis
#1.3	Removed ImageConversion and SpiralExtraction functions to include in separate functions file MorphometricExtraction_Functions.r
#1.4	Modified SpiralAnalysis to be integrated with the new data file structure of Spiral files
#1.4.1	Added option in to also export the raw values for the regression
#1.4.2	Modified SpiralAnalysis, so that actually T is the dependent variable to Theta, not the other way around. This woild not be wrong, 
#	but make results less intiutive to interpret.
#--------------------------------------------
#--------------------------------------------


