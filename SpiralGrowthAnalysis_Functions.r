#Spiral Growth Analysis
#Measuring whether or not a spiral growth process follows the logarithmic spiral.
#Input data set:
#	-Images in the ppm format

#Author: Manuel Weinkauf (Manuel.Weinkauf@unige.ch)
#Version: 1.3
#Date: 1 Novenber 2014

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
#   DataName: Name part of datasets (same for all files).               #
#             *string*                                                  #
#   StartNum: Smallest of row of continuous numbers used to number...   #
#             files.                                                    #
#             *numeric (integer)*                                       #
#    StopNum: Largest of row of continuous numbers used to number files.#
#             *numeric (integer)*                                       #
# Input dataset: A set of data points along the spiral, given as...     #
#                radial distance t (first column) and  angle theta in...#
#                radians (second column). Values are supposed to be...  #
#                normalized for radius one and rotated such that...     #
#                theta[1] is zero!                                      #
# Output dataset: Matrix with dip, intercept, and R2 for model II...    #
#                 linear regression (ranged major axis) of log-...      #
#                 transformed distances t dependent on theta.           #
#########################################################################

#Loading packages
require(lmodel2)

SpiralAnalysis<-function(DataName, StartNum, StopNum) {
	#Setting up results matrix
	Res<-matrix(NA, (StopNum-StartNum+1), 3)
	colnames(Res)<-c("Dip", "Intercept", "Adjusted R squared")
	for (k in StartNum:StopNum) {
		#Reading data
		DName<-paste(DataName, k, ".txt", sep="")
		Data<-read.table(DName, header=TRUE)
		
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
		#Calculating model II linear regression
		Fit<-lmodel2(Data[,2]~Data[,1], range.x="interval", range.y="relative")
		Res[k,1]<-Fit$regression.results[which(Fit$regression.result[,"Method"]=="RMA"),"Slope"]
		Res[k,2]<-Fit$regression.results[which(Fit$regression.result[,"Method"]=="RMA"),"Intercept"]
		##Calculating R2 of RMA
		y.mean<-mean(Data[,2])
		SS.tot<-sum((Data[,2]-y.mean)^2)
		SS.res<-sum((Data[,2]-(Data[,1]*Res[k,1]+Res[k,2]))^2)
		Res[k,3]<-1-(SS.res/SS.tot)
	}
	
	#Export results
	write.table(Res,"SpiralFitting.txt",sep="\t")
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
#--------------------------------------------
#--------------------------------------------


