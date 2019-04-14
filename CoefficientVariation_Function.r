#Calculate the coefficient of variation
#Input data set:
#	-Measurement data for which coefficient of variation should be calculated
#Further reading:
#	McKay, A. T. (1932) "Distribution of the coefficient of variation and the extended...
#		"t" distribution". Journal of the Royal Statistical Society 95 (4): 695-698.
#	Vangel, M. G. (1996) "Confidence intervals for a normal coefficient of variation"...
#		The American Statistician 50 (1): 21-26.

#Author: Manuel Weinkauf (Manuel.Weinkauf@unige.ch)
#Version: 1.1
#Date: 16 November 2016

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.#
#To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/.                   #
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#Creation of test dataset
#Dat<-matrix(c(rep(1, 10), rep(2, 10), sample(4:7, 10, replace=TRUE), sample(2:12, 10, replace=TRUE)), 20, 2)
#colnames(Dat)<-c("Group", "Value")
#Dat2<-cbind(Dat, sample(17:32, 20, replace=TRUE))
#colnames(Dat2)<-c("Group", "Width", "Length")

#**************************************************************************************
#Setting working directory
#setwd("C:/R_TestData")

#########################################################################
# Function to calculate coefficient of variation incl. 95% confidence...#
#   interval                                                            #
# Necessary input variables:                                            #
#   Data: Data matrix, samples in rows, variables in columns. First...  #
#         column is expected to contain groups encoding.                #
#         *matrix*                                                      #
#   conf: Desired level of confidence.                                  #
#         *real*                                                        #
#         default=0.95                                                  #
# Output data: Coefficient of variation and associated confidence...    #
#              interval per group in all parameters                     #
# Input dataset: Data matrix, samples in rows, variables in columns     #
#########################################################################

#Loading packages

CoeffVar<-function (Data, conf=0.95) {
	#Prepare data
	Groups<-as.factor(Data[,1])
	{if (!is.null(colnames(Data))) {Vars<-colnames(Data)[-1]}
	else {Vars<-paste("Var.", 1:(ncol(Data)-1), sep="")}}
	Data<-as.matrix(Data[,-1])
	conf<-1-conf
	
	#Prepare results table
	Res<-list()
	Res$Vangel<-Res$McKay<-array(NA, dim=c(length(levels(Groups)), 3, length(Vars)), dimnames=list(paste("Group.", levels(Groups), sep=""), c("Coeff.Var", "Coeff.Var.LCI", "Coeff.Var.UCI"), Vars))

	#Calculate coefficient of variation and confidence intervals
	for (i in 1:length(Vars)) {
		#Coefficient of variation
		SD<-by(as.vector(Data[,i]), Groups, sd, na.rm=TRUE)
		Mean<-by(as.vector(Data[,i]), Groups, mean, na.rm=TRUE)
		Res$Vangel[,"Coeff.Var",i]<-Res$McKay[,"Coeff.Var",i]<-SD/Mean
		
		#Confidence auxiliary variables
		v<-by(as.vector(Data[,i]), Groups, length)-1
		u1<-qchisq(p=(1-conf/2), df=v)
		u2<-qchisq(p=(conf/2), df=v)
		
		#Confidence intervals
		##McKay
		Res$McKay[,"Coeff.Var.LCI",i]<-Res$McKay[,"Coeff.Var",i]*(((((u1/(v+1))-1)*(Res$McKay[,"Coeff.Var",i]**2))+(u1/v))**(-0.5))
		Res$McKay[,"Coeff.Var.UCI",i]<-Res$McKay[,"Coeff.Var",i]*(((((u2/(v+1))-1)*(Res$McKay[,"Coeff.Var",i]**2))+(u2/v))**(-0.5))
		##Vangel
		Res$Vangel[,"Coeff.Var.LCI",i]<-Res$Vangel[,"Coeff.Var",i]*((((((u1+2)/(v+1))-1)*(Res$Vangel[,"Coeff.Var",i]**2))+(u1/v))**(-0.5))
		Res$Vangel[,"Coeff.Var.UCI",i]<-Res$Vangel[,"Coeff.Var",i]*((((((u2+2)/(v+1))-1)*(Res$Vangel[,"Coeff.Var",i]**2))+(u2/v))**(-0.5))
	}
	
	#Return results
	return(Res)
}

#--------------------------------------------

#Example
#CT1<-CoeffVar(Dat)
#CT2<-CoeffVar(Dat2)

#--------------------------------------------
#--------------------------------------------
#Version History
#1.0	Finished Program
#1.1	Added method by McKay
#--------------------------------------------
#--------------------------------------------
