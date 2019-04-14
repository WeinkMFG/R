#Bootstrapping of confidence intervals for means and calculating standard deviations and associated...
#	confidence intervals of large quantities of data
#Input data set:
#	-A number of data values, obtained from measurement of several species per sample
#	-One value per variable per specimen, at least three values per sample
#	The whole function structure has been reworked since version 1. It is no longer...
#	the case, that samples go in columns and values in rows. This old input file structure...
#	had two significant disadvantages: (1) It allowed only one measured parameter per file...
#	to be treated, and (2) it necessitated to arrange all data in a block matrix with the...
#	number of rows according to the sample with most data in it, and coding the last rows...
#	in all samples with fewer data as NA.
#	With the new file structure both disadvantages has been eliminated. The first column now...
#	codes the sample in an unambiguous way (e.g. numerical or alphanumerical), while as many...
#	measured parameters as desired can be included in the 2nd-nth column. Since NA values are...
#	still supported, not all data values must be present in all samples.
#Further reading:
#	Dixon, Ph. M. (2002) "Bootstrap resampling". IN "Encyclopedia of Environmetrics" (John...
#		Wiley and Sons: Chichester), pp. 212-20.
#	Tabachnik, B. G. and Fidell, L. S. (1996) "Using Multivariate Statistics", 3rd ed. (Harper...
#		Collins: New York).
#	Tabor, J. (2010) "Investigating the Investigative Task: Testing for Skewness...
#		An Investigation of Different Test Statistics and their Power to Detect Skewness"...
#		Journal of Statistics Education 18/2: 1-13.
#	Sheskin, D. J. (2011) "Handbook of parametric and nonparametric statistical procedures",...
#		 1886 pp. (Boca Raton, London, New York: Chapman & Hall/CRC Press).

#Author: Manuel Weinkauf (Manuel.Weinkauf@unige.ch)
#Version: 3.0
#Date: 26 March 2014

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.#
#To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/.                   #
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#Creation of test dataset
#setwd("C:/R_TestData/BootCI")
#Test1<-matrix(sample(1:100),100,1)
#Test2<-matrix(sample(1:8,replace=TRUE),100,1)
#Sam<-c(rep.int(1,50),rep.int(4,50))
#Test<-cbind(Sam,Test1,Test2)
#colnames(Test)<-c("Sample","Value 1", "Value 2")
#write.table(Test,"Test.txt",sep="\t")
#file.show("Test.txt")

#**************************************************************************************
#Setting working directory
#setwd("C:/R_TestData/BootCI")

#########################################################################
# Function to apply bootstrapping CI on data                            #
# Necessary packages: boot, moments, pastecs                            #
# Necessary input variables:                                            #
#   Input: Name of input matrix (in case of file header assumed).       #
#          *character* if Import==TRUE                                  #
#          *object* if Import==FALSE                                    #
#   Output: Namepart of output file (without suffix), also used as...   #
#           name for output variables.                                  #
#           *character*                                                 #
#   ConfLevel: Desired confidence interval.                             #
#              *numric (real)*, 0<ConfLevel<1                           #
#              default=0.95                                             #
#   Replicates: Number of desired replications                          #
#               *numeric (integer)*                                     #
#               default=10000                                           #
#   Adapt: Shall bootstrapping adapt to data?                           #
#          *logical*                                                    #
#          TRUE: Use BCA when possible, else basic bootstrapping        #
#          FALSE: Always use basic bootstrapping                        #
#          default=TRUE                                                 #
#   NAV: Coding for missing values in data.                             #
#        default=NA                                                     #
#   Import: Should the data be imported from a txt file?                #
#           *logical*                                                   #
#           TRUE: Read data from file                                   #
#           FALSE: read data from R-object                              #
#           default=FALSE                                               #
#   Export: Should results be exported as .txt file?                    #
#           *logical*                                                   #
#           TRUE: Export results as file                                #
#           FALSE: Stores results in R-object                           #
#           default=TRUE                                                #
# Output data: Mean and 95% confidence intervals; one result object...  #
#              will be produced for each variable, containing one...    #
#              line each per sample.                                    #
# Input dataset: Matrix object, values in rows, variables in columns... #
#                first column contains sample encoding, all other...    #
#                columns contain measurements for which confidence...   #
#                intervals shall be calculated.                         #
#########################################################################

#Loading packages
require(boot)
require(moments)
require(pastecs)

ConfInt<-function(Input, Output, ConfLevel=0.95, Replicates=10000, Adapt=TRUE, NAV=NA, Import=FALSE, Export=TRUE){
	#Creation of function for 'statistic'
	my_mean<-function(data, indices){
		d<-data[indices]
		mean(d, na.rm=TRUE)
	}

	#Reading input data
	{if (Import==TRUE) {Dat<-read.table(Input, header=TRUE, na.strings=NAV)}
	else {Dat<-Input}}
	NumSum<-unique(Dat[,1])
	FileNames<-paste(Output,"_", colnames(Dat)[-1], ".txt", sep="")
	VariableNames<-paste(Output, colnames(Dat)[-1], sep=".")
	VariableNames<-gsub(" ", ".", VariableNames, fixed=TRUE)#Eliminates white spaces in variable names

	for (x in 2:(dim(Dat)[2])) {
		#Creating output objects
		Sam<-list()
		Skewness<-TRUE
		Results<-as.data.frame(matrix(NA, length(NumSum), 8))
		colnames(Results)<-c("Mean", "Mean.LCI", "Mean.UCI", "Boot.Method", "q.level", "Stand.Dev", "Stand.Dev.LCI", "Stand.Dev.UCI")
		rownames(Results)<-NumSum

		#Performing analyses
		for (i in 1:length(NumSum)) {
			#Writing datasubset into Sam
			Sam$Vals<-Dat[which(Dat[,1]==NumSum[i]), x]
			#Testing suitability of input data
			NVal<-stat.desc(Sam$Vals)[[1]]
			if (Adapt==TRUE && NVal>=Replicates) {stop("If you want to use accelerated bootstrap Replicates must be larger than the largest sample size!")}
			
			{if (NVal>=3) {
				#Performing bootstrap
				Boot_Res<-boot(data=Sam$Vals, statistic=my_mean, R=Replicates, sim="ordinary")
				plot(Boot_Res)

				#Checking for skewness of data
				xPrime<-mean(Sam$Vals, na.rm=TRUE)
				Diff<-matrix(NA, length(Sam$Vals), 3)
				for (j in 1:(length(Sam$Vals))) {
					Diff[j,1]<-Sam$Vals[j]-xPrime
					Diff[j,2]<-(Diff[j,1])^3
					Diff[j,3]<-(Diff[j,1])^2
				}
				Sk<-((1/NVal)*sum(Diff[,2], na.rm=TRUE))/(((1/NVal)*sum(Diff[,3], na.rm=TRUE))^(3/2))
				SkDev<-2*(sqrt(6/NVal))
				{if(Sk<SkDev | is.nan(Sk)) {Skewness<-FALSE}
					else {Skewness<-TRUE}
				}

				#Calculating confidence intervals and writing results into Results object
				{if (Adapt==TRUE) {
					{if (Skewness==TRUE) {CI<-boot.ci(Boot_Res, conf=ConfLevel, type="basic"); Results[i,1]<-as.numeric(CI$t0); Results[i,2]<-as.numeric(CI$basic[4]); Results[i,3]<-as.numeric(CI$basic[5]); Results[i,4]<-"Basic"; Results[i,5]<-ConfLevel}
					else {CI<-boot.ci(Boot_Res, conf=ConfLevel, type="bca"); Results[i,1]<-as.numeric(CI$t0); Results[i,2]<-as.numeric(CI$bca[4]); Results[i,3]<-as.numeric(CI$bca[5]); Results[i,4]<-"BCA"; Results[i,5]<-ConfLevel}}
					}
					else {CI<-boot.ci(Boot_Res, conf=ConfLevel, type="basic"); Results[i,1]<-as.numeric(CI$t0); Results[i,2]<-as.numeric(CI$basic[4]); Results[i,3]<-as.numeric(CI$basic[5]); Results[i,4]<-"Basic"; Results[i,5]<-ConfLevel}
				}
				
				#Calculating standard deviation including 95% confidence interval of data
				Results[i,6]<-sd(Sam$Vals, na.rm=TRUE)
				Results[i,7]<-Results[i,6]*sqrt((NVal-1)/qchisq((1-((1-ConfLevel)/2)), (NVal-1)))
				Results[i,8]<-Results[i,6]*sqrt((NVal-1)/qchisq((0+((1-ConfLevel)/2)), (NVal-1)))
			}
			else {warning("There were some samples with less than three values. No values calculated!")}}
		}

		{if (Export==TRUE) {write.table(Results, FileNames[x-1], sep="\t")}
		else {assign(VariableNames[x-1], Results, envir=.GlobalEnv)}}
	}
	if (Export==FALSE) {print("The results have been stored in variables"); print(VariableNames)}
}

#--------------------------------------------

#Example

#ConfInt(Test, "Res", Replicates=1000, Export=FALSE)

#ConfInt("Test.txt", "Res", Replicates=1000, Adapt=FALSE, Export=FALSE)

#ConfInt("Test.txt", "Res", ConfLevel=0.6, Replicates=1000, Adapt=FALSE)
#file.show("Res.txt")

#--------------------------------------------
#--------------------------------------------
#Version History
#1.0	Finished Program
#1.1	Adds option to specify coding for missing values
#2.0	Complete rework for more intuitive and logical input file format;
#	possibility added to have more than one response variable;
#	option added not to export data
#2.1	Mean value exported alongside confidence intervals;
#	included security check for suitability of data for accelerated bootstrap
#2.2	Possibility added not to import data but read directly from R variable
#2.2.1	Default NA value switched to NA
#2.2.2	Changed naming of output variables (if Export=FALSE)
#3.0	Functionality added that also calculates the standard deviations and their
#	confidence interval of the data
#--------------------------------------------
#--------------------------------------------
