#Test goodness-of-fit for binomial GLM with partly continuous predictors
#Input data set:
#	- The dataset the GLM should be based on
#Further reading:
#	Stukel, Th. A. (1988) "Generalized logistic models". Journal of the American Statistical
#		Association 83 (402): 426-431.

#Author: Manuel Weinkauf (Manuel.Weinkauf@unige.ch)
#Version: 1.0
#Date: 18 March 2015

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.#
#To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/.                   #
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#Creation of test dataset
#set.seed(1000)
#Suc<-matrix(NA, 20, 4)
#Suc[,1]<-c(rep(0, 10), rep(1, 10))
#Suc[,2]<-sample(0:1, 20, replace=TRUE)
#Suc[,3]<-sample(0:4, 20, replace=TRUE)
#Suc[,4]<-runif(20, min=1.2, max=5.3)
#colnames(Suc)<-c("Success", "Param1", "Param2", "Param3")
#Suc<-as.data.frame(Suc)

#**************************************************************************************
#Setting working directory
#setwd("C:/R_TestData")

#########################################################################
# Function to calculate the goodness-of-fit of the chosen binomial GLM  #
# Necessary input variables:                                            #
#   Data: Data frame to be used for test.                               #
#         *data.frame*                                                  #
#   GLM: Symbolic description of linear model.                          #
#        *string*                                                       #
#   Link: Link function for binomial model.                             #
#         *string* of either "logit", "probit", or "cloglog".           #
# Output data: Results of a likelihood ratio test for goodness of fit...#
#              of the tested model.                                     #
# Input dataset: Data matrix for which binomial linear model was...     #
#                calculated.                                            #
#########################################################################

#Loading packages
require("lmtest")

Stukel<-function (Data, GLM, Link="logit") {
	#Test for valid link function
	if (!any(Link=="logit", Link=="probit", Link=="cloglog")) {stop("'Link' must be either 'logit', 'probit', or 'cloglog'!")}

	#Calculate model which fit should be tested
	Test.Model<-glm(eval(parse(text=GLM)), data=Data, family=binomial(link=Link))
	
	#Calculate null model
	g<-predict(Test.Model)
	za<-zb<-vector(length=length(g), mode="numeric")
	for (i in 1:(length(g))) {
		if (g[i]>=0) {za[i]<-g[i]^2}
		if (g[i]<0) {zb[i]<-g[i]^2}
	}
	Data.Null<-cbind(Data, za, zb)
	GLM.Null<-paste(GLM, "za", "zb", sep="+")
	Null.Model<-glm(eval(parse(text=GLM.Null)), data=Data, family=binomial(link=Link))
	
	#Test for significance of test model (likelihood ratio test)
	Res<-lrtest(Test.Model, Null.Model)
	return(Res)
}

#--------------------------------------------

#Example
#Stukel(Suc, "Success~Param1+Param3")
#Stukel(Suc, "Success~Param1+Param3", Link="probit")
#Stukel(Suc, "Success~Param1+Param2+Param3+Param1*Param2")

#--------------------------------------------
#--------------------------------------------
#Version History
#1.0	Finished Program
#--------------------------------------------
#--------------------------------------------
