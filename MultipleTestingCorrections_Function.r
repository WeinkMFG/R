#Correcting the significance levels for the False Discovery Rate and calculating the Bonferroni Correction
#Input data set:
#	-A matrix or vector with one column, containing calculated p-values for each comparison
#Further Reading:
#	Benjamini, Y. and Hochberg, Y. (1995) "Controlling the False Discovery...
#			Rate: A practical and powerful approach to multiple testing", Journal...
#			of the Royal Statistical Society B, 57 (1): 289-300.

#Author: Manuel Weinkauf (Manuel.Weinkauf@unige.ch)
#Version: 1.0
#Date: 5 January 2012

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.#
#To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/.                   #
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#Creation of test dataset
#pVal<-matrix(c(0.7590,0.0004,0.0001,0.0019,0.0095,0.0201,0.0278,0.0298,0.0344,0.3240,0.0459,0.4262,0.5719,0.6528,1.000),15,1)
#write.table(pVal,"pVal.txt",sep="\t")

#pVal2<-c(0.7590,0.0004,0.0001,0.0019,0.0095,0.0201,0.0278,0.0298,0.0344,0.3240,0.0459,0.4262,0.5719,0.6528,1.000)

#**************************************************************************************
#Setting working directory
#setwd("C:/R_TestData")

#########################################################################
# Correction for multiple testing.                                      #
# Necessary input variables:                                            #
#    Input: Name of input matrix (header assumed).                      #
#           *string* or *matrix* or *vector*                            #
#    qo: Desired original significance level.                           #
#        *numeric (real)*, 1<-q<-0                                      #
#        default=0.05                                                   #
#    Import: Do you want to import a .txt file?                         #
#            *logical*                                                  #
#            TRUE=Import data file (Input=*string*).                    #
#            FALSE=Do not import data file (Input=*matrix* or *vector*).#
#            default=TRUE                                               #
# Output data: New q-value of significance.                             #
# Input dataset: Matrix object, one column containing calculated...     #
#                p-values in rows, or vector of calculated p-values.    #
#########################################################################

Benjamini<-function (Input, qo=0.05, Import=TRUE) {
	if (any(qo>1, qo<0)) {stop("qo must have values ranging between 0 and 1!")}
	if (Import==TRUE) {Input<-as.matrix(read.table(Input, header=TRUE, sep="\t"))}
	pList<-sort(Input)
	
	#Setting up results table
	Res<-matrix(NA, length(pList), 3)
	Res[,1]<-pList
	
	#Calculating new p-values
	for (i in (length(pList)):1) {
		qn<-(i/(length(pList)))*qo
		Res[i,2]<-qn
	}
	
	for (i in (length(pList)):1) {
		if (Res[i,1]<=Res[i,2]) {Res[i,3]<-1}
		else {Res[i,3]<-0}
	}
	
	c<-length(pList)
	while (Res[c,3]==0 && c>0) {
		c<-c-1
	}
	
	#Presenting output
	{if (c>0) {CP<-round(Res[c,2],digits=4); writeLines(paste("Level of significance corrected for the False Discovery Rate \nBenjamini and Hochberg (1995) \nq* =", CP, sep=" "))}
	else {writeLines("All null hypotheses must be accepted according to the False Discovery Rate!")}}
	Bon<-round(qo/length(pList), digits=4)
	writeLines(paste("Bonferroni corrected level of significance \nQ* =", Bon, sep=" "))
}

#--------------------------------------------

#Example

#Benjamini(pVal, Import=FALSE)
#Benjamini(pVal2, Import=FALSE)
#Benjamini("pVal.txt")
#Benjamini(pVal, qo=0.1, Import=FALSE)
#Benjamini(pVal, qo=1.1, Import=FALSE)

#--------------------------------------------
#--------------------------------------------
#Version History
#1.0	Finished Program
#--------------------------------------------
#--------------------------------------------
