#A collection of tools for linear regression
#Input data set:
#	-Dependent and independent variables for linear regression.
#	-Parameters of a fitted linear regression
#Further reading:
#	Helsel, D. R. and Hirsch, R. M. (2002) "Statistical Methods in Water...
#		Resources". IN "Hydrologic Analysis and Interpretation", chap. A3, vol. iv...
#		of "Techniques of Water-Resources Investigations Reports" (United States...
#		Geological Survey), pp. 1-510.
#	Kendall, M. G. (1938) "A New Measurement of Rank Correlation". Biometrika 30 (1-2)...
#		81-93.
#	Richter, S. J. and Stavn, R. H. (2014) "Determining functional relations in multivariate...
#		oceanographic systems: Model II multiple linear regression". Journal of Atmospheric...
#		and Oceanic Technology 31: 1663-72.
#	Sen, P. K. (1968) "Estimates of the Regression Coefficient Based on Kendall's Tau"...
#		Journal of the American Statistical Association 63 (324): 1379-89.
#	Theil, H. (1950) "A Rank-Invariant Method of Linear and Polynomial Regression...
#		Analysis, III". Proceedings of the Koninklijke Nederlandse Akademie van...
#		Wetenschappen 53 (9): 1397-412.
#	https://stackoverflow.com/questions/14069629/plotting-confidence-intervals
#	https://stats.stackexchange.com/questions/101318/understanding-shape-and-calculation-of-confidence-bands-in-linear-regression
#	

#Author: Manuel Weinkauf (Manuel.Weinkauf@unige.ch)
#Version: 2.0.2
#Date: 14 January 2020

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.#
#To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/.                   #
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#Creation of test dataset
#Var1<-c(11.36055,6.98143,6.89490,14.60099,12.52414, 5.40229,6.45010,1.56994,2.15346,22.30509,4.67338, 5.00835,5.32663,5.45830,9.98043,5.67296,6.88813, 3.20948,4.56163,6.55727,4.62544,5.48337,3.87877, 2.80227)
#Var2<-c(2.35529,1.48723,1.15241,3.00136,1.59180, 2.53472,2.20666,0.17953,0.45311,3.28352,2.05125, 0.52108,1.94375,2.16183,2.86684,1.80909,3.10551, 3.32273,1.69203,1.06135,1.05097,0.69146,0.71080, 0.48471)
#Var3<-c(8.40735,5.85045,6.90855,11.30800,10.02730, 3.39939,5.43069,1.08864,1.84174,17.94280,4.85276, 1.19745,4.99449,5.05150,9.11405,8.39094,7.57381, 4.38397,3.77923,4.91364,3.09047,3.66245,2.40817, 1.76638)
#ModII<-cbind(Var1, Var2, Var3)
#colnames(ModIII)<-c("Temperature", "Size")
#Sizes<-matrix(c(1, 1, 1, 4, 4, 4, 7, 7, 7, 10, 10, 10, 2.1, 2.3, 2.2, 3.5, 3.1, 3.2, 4.2, 5.0, 4.8, 6.1, 6.6, 6.2, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3), 12, 3)
#colnames(Sizes)<-c("Salinity", "Length", "Block")

#**************************************************************************************
#Setting working directory
#setwd("C:/R_TestData")

#########################################################################
# Function to perform multiple model II linear regression.              #
# Necessary input variables:                                            #
#   Y: Dependent variable.                                              #
#      *vector*                                                         #
#   X: Independent variables.                                           #
#      *matrix*                                                         #
#   boot.n: Number of bootstrap replications.                           #
#           *numeric (integer)*                                         #
#           default=9999                                                #
#   conf.int: Desired confidence interval.                              #
#             *numeric (real)*                                          #
#             default=0.95                                              #
#   tails: Decide for two-tailed or positive upper or lower tail test.  #
#          *character*, one of "two-tailed", "upper", or "lower"        #
#          default="two.tailed"                                         #
# Output data: Model II linear model.                                   #
# Input dataset: Data for multiple model II linear regression.          #
#########################################################################

Mult.lm2<-function (Y, X, boot.n=9999, conf.int=0.95, tails="two.tailed") {
	#Test data consistency
	if (length(Y)!=nrow(X)) {stop("Dependent and independent variables must be of same length!")}
	if (nrow(X)<2) {stop("With only one independent variable, use package lmodel2 instead!")}
	if (tails!="two.tailed" & tails!="upper" & tails!="lower") {stop("tails must be one of either 'two.tailed', 'upper', or 'lower'!")}
	
	#Link data
	{if (is.null(colnames(X))) {X.names<-paste("X", 1:ncol(X), sep=".")}
	else {X.names<-colnames(X)}}
	Dat<-cbind(X, Y)
	colnames(Dat)<-c(X.names, "Y")
	
	#Data cleaning
	if (any(complete.cases(Dat)==FALSE)) {
		warning("NA's detected and removed!")
		Dat<-Dat[complete.cases(Dat),]
	}
	
	#Calculate base parameters
	Params<-list()
	Params$n<-nrow(Dat)
	Params$means<-apply(Dat, 2, mean)
	Params$sd<-apply(Dat, 2, sd)
	Params$cov.matrix<-cov(Dat)
	Params$cor.matrix<-cor(Dat)
	Params$cov.eigen<-eigen(Params$cov.matrix)
	Params$cor.eigen<-eigen(Params$cor.matrix)
	
	#Calculate regression coefficients
	L.Model<-list()
	L.Model$standardized<-L.Model$unstandardized<-matrix(NA, ncol(X)+1, 4)
	rownames(L.Model$standardized)<-rownames(L.Model$unstandardized)<-c("Intercept", X.names)
	{if (tails=="two.tailed") {P<-"p-value (two-tailed)"}
	else if(tails=="upper") {P<-"p-value (>0)"}
	else {P<-"p-value (<0)"}}
	colnames(L.Model$standardized)<-colnames(L.Model$unstandardized)<-c("Value", "CI.Low", "CI.High", P)
	##Covariance matrix
	Comps<-vector(mode="numeric", length=nrow(Params$cov.eigen$vectors))
	Slopes<-vector(mode="numeric", length=nrow(Params$cov.eigen$vectors)-1)
	for (i in 1:(nrow(Params$cov.eigen$vectors)-1)) {
		Comps[i]<-Params$cov.eigen$vectors[i,ncol(Params$cov.eigen$vectors)]/Params$cov.eigen$vectors[nrow(Params$cov.eigen$vectors),ncol(Params$cov.eigen$vectors)]*Params$means[i]
		Slopes[i]<--Params$cov.eigen$vectors[i,ncol(Params$cov.eigen$vectors)]/Params$cov.eigen$vectors[nrow(Params$cov.eigen$vectors),ncol(Params$cov.eigen$vectors)]
	}
	Comps[length(Comps)]<-Params$means[length(Params$means)]
	L.Model$unstandardized[,"Value"]<-c(sum(Comps), Slopes)
	
	##Correlation matrix
	Comps<-vector(mode="numeric", length=nrow(Params$cor.eigen$vectors))
	Slopes<-vector(mode="numeric", length=nrow(Params$cor.eigen$vectors)-1)
	for (i in 1:(nrow(Params$cor.eigen$vectors)-1)) {
		Comps[i]<-(Params$cor.eigen$vectors[i,ncol(Params$cov.eigen$vectors)]/Params$cor.eigen$vectors[nrow(Params$cor.eigen$vectors),ncol(Params$cor.eigen$vectors)])*(Params$sd[length(Params$sd)]/Params$sd[i])*Params$means[i]
		Slopes[i]<--(Params$cor.eigen$vectors[i,ncol(Params$cov.eigen$vectors)]/Params$cor.eigen$vectors[nrow(Params$cor.eigen$vectors),ncol(Params$cor.eigen$vectors)])*(Params$sd[length(Params$sd)]/Params$sd[i])
	}
	Comps[length(Comps)]<-Params$means[length(Params$means)]
	L.Model$standardized[,"Value"]<-c(sum(Comps), Slopes)
	
	#Bootstrapping evaluation
	##Bootstrap percentiles
	Boot.values<-list()
	Boot.values$standardized<-Boot.values$unstandardized<-matrix(NA, boot.n, ncol(X))
	colnames(Boot.values$unstandardized)<-colnames(Boot.values$standardized)<-X.names

	for(i in 1:boot.n) {
		Boot.index<-sample(1:Params$n, Params$n, replace=TRUE)
		Boot.sample<-Dat[Boot.index,]
		Boot.mat<-list()
		Boot.mat$unstandardized<-eigen(cov(Boot.sample))$vectors
		Boot.mat$standardized<-eigen(cor(Boot.sample))$vectors
		for (j in 1:(nrow(Boot.mat$unstandardized)-1)) {
			Boot.values$unstandardized[i,j]<--Boot.mat$unstandardized[j,ncol(Boot.mat$unstandardized)]/Boot.mat$unstandardized[nrow(Boot.mat$unstandardized),ncol(Boot.mat$unstandardized)]
			Boot.values$standardized[i,j]<--(Boot.mat$standardized[j,ncol(Boot.mat$standardized)]/Boot.mat$standardized[nrow(Boot.mat$standardized),ncol(Boot.mat$standardized)])*(Params$sd[length(Params$sd)]/Params$sd[j])
		}
	}
	
	##Calculate confidence limits of slopes
	for (i in 1:ncol(Boot.values$unstandardized)) {
		L.Model$unstandardized[i+1,"CI.Low"]<-quantile(Boot.values$unstandardized[,i], probs=(1-conf.int)/2, type=8)
		L.Model$unstandardized[i+1,"CI.High"]<-quantile(Boot.values$unstandardized[,i], probs=1-((1-conf.int)/2), type=8)
	}
	for (i in 1:ncol(Boot.values$standardized)) {
		L.Model$standardized[i+1,"CI.Low"]<-quantile(Boot.values$standardized[,i], probs=(1-conf.int)/2, type=8)
		L.Model$standardized[i+1,"CI.High"]<-quantile(Boot.values$standardized[,i], probs=1-((1-conf.int)/2), type=8)
	}
	
	#Calculate p-values
	for (i in 1:ncol(Boot.values$unstandardized)) {
		N<-nrow(Boot.values$unstandardized)
		Boot.H<-length(which(Boot.values$unstandardized[,i]>0))/N
		Boot.L<-length(which(Boot.values$unstandardized[,i]<0))/N
		{if (tails=="two.tailed") {L.Model$unstandardized[i+1,4]<-round(min(c(Boot.H, Boot.L))*2, digits=3)}
		else if (tails=="upper") {L.Model$unstandardized[i+1,4]<-round(1-Boot.H, digits=3)}
		else {L.Model$unstandardized[i+1,4]<-round(1-Boot.L, digits=3)}}
	}
	for (i in 1:ncol(Boot.values$standardized)) {
		N<-nrow(Boot.values$standardized)
		Boot.H<-length(which(Boot.values$standardized[,i]>0))/N
		Boot.L<-length(which(Boot.values$standardized[,i]<0))/N
		{if (tails=="two.tailed") {L.Model$standardized[i+1,4]<-round(min(c(Boot.H, Boot.L))*2, digits=3)}
		else if (tails=="upper") {L.Model$standardized[i+1,4]<-round(1-Boot.H, digits=3)}
		else {L.Model$standardized[i+1,4]<-round(1-Boot.L, digits=3)}}
	}
	
	#Return results
	return(L.Model)
}

#########################################################################
# Bivariate model III robust regression.                                #
# Necessary input variables:                                            #
#   Y: Dependent variable.                                              #
#      *vector*                                                         #
#   X: Independent variables.                                           #
#      *vector*                                                         #
#   Plot: Shall the results be plotted?                                 #
#         *logical*                                                     #
#         TRUE: Plot results                                            #
#         FALSE: Do not plot results                                    #
#         default=TRUE                                                  #
# Output data: Statistics of a Kendall-Theil robust linear regression...#
#              including dip (with 95% confidence interval) and...      #
#              intercept of regression line, Kendall's tau...           #
#              (correlation coefficient), significance p of...          #
#              regression line, and coefficient of determination R2.    #
# Input dataset: Data for model III linear regression.                  #
#########################################################################

Kendall<-function (Y, X, Plot=TRUE) {
	#Test data consistency
	if (length(Y)!=length(X)) {stop("Dependent and independent variables must be of same length!")}
	if (length(X)<3) {stop("At least three data points needed!")}

	#Link data
	Input<-cbind(X, Y)
	colnames(Input)<-c("X", "Y")
	
	#Data cleaning
	if (any(complete.cases(Input)==FALSE)) {
		warning("NA's detected and removed!")
		Input<-Input[complete.cases(Input),]
	}
	
	#Setting up matrix for different slopes
	Slopes<-vector(mode="numeric", length=sum((nrow(Input)-1):1))

	#Setting up temporary-pairs matrix
	n<-nrow(Input)-1
	c<-1
	Temp<-matrix(NA, 2, 2)

	#Setting up results matrix
	Res<-vector(mode="numeric", length=7)
	names(Res)<-c("Dip", "Upper CI Dip", "Lower CI Dip", "Intercept", "Tau", "p", "R2")

	#Calculating regression
	M<-0
	P<-0
	for (j in 1:n){
		for (i in j:n){
			#Calculating slopes
			Temp[1,1]<-Input[j,1]
			Temp[1,2]<-Input[j,2]
			Temp[2,1]<-Input[(i+1),1]
			Temp[2,2]<-Input[(i+1),2]
			{if (Temp[1,1]==Temp[2,1]) {SL<-NA}
			else if (Temp[1,1]<Temp[2,1]) {SL<-(Temp[2,2]-Temp[1,2])/(Temp[2,1]-Temp[1,1])}
			else {SL<-(Temp[1,2]-Temp[2,2])/(Temp[1,1]-Temp[2,1])}
			}
			Slopes[c]<-SL
			c<-c+1
			
			#Calculating significance
			{if (Temp[1,1]>Temp[2,1] & Temp[1,2]>Temp[2,2]) {P<-P+1}
			else if (Temp[1,1]>Temp[2,1] & Temp[1,2]<Temp[2,2]) {M<-M+1}
			else if (Temp[2,1]>Temp[1,1] & Temp[2,2]>Temp[1,2]) {P<-P+1}
			else if (Temp[1,1]==Temp[2,1] | Temp[1,2]==Temp[2,2]) {}
			else {M<-M+1}
			}
		}
	}
	
	#Order slopes
	Slopes<-sort(Slopes)

	#Calculating dip of regression line
	Res[1]<-median(Slopes, na.rm=TRUE)
	Np<-sum((nrow(Input)-1):1)
	N<-nrow(Input)
	Ru<-round(((Np+(1.96*(sqrt((N*(N-1)*(2*N+5))/18))))/2)+1)
	Rl<-round(((Np-(1.96*(sqrt((N*(N-1)*(2*N+5))/18))))/2))
	{if (Ru<=length(Slopes)) {Res[2]<-Slopes[Ru]} else {Res[2]<-"Too few values, unable to calculate"}}
	{if (Rl>=1) {Res[3]<-Slopes[Rl]} else {Res[3]<-"Too few values, unable to calculate"}}
	Res[4]<-median(Input[,2])-(as.numeric(Res[1])*median(Input[,1]))

	#Calculating Tau
	Res[5]<-(P-M)/Np

	#Counting Ties
	XVals<-matrix(c(sort(unique(Input[,1])), rep(0, length(unique(Input[,1])))), length(unique(Input[,1])), 2)
	YVals<-matrix(c(sort(unique(Input[,2])), rep(0, length(unique(Input[,2])))), length(unique(Input[,2])), 2)
	for (j in 1:nrow(Input)) {
		TV<-Input[j,1]
		XVals[match(TV, XVals[,1]),2]<-XVals[match(TV, XVals[,1]),2]+1
		TV<-Input[j,2]
		YVals[match(TV, YVals[,1]),2]<-YVals[match(TV, YVals[,1]),2]+1
	}
	TieList<-c(XVals[,2], YVals[,2])
	TieMin<-min(TieList)
	if (TieMin==0 | TieMin==1) {TieMin<-2}
	TieMax<-max(TieList)
	{if (TieMax>1) {Ties<-TRUE}
		else {Ties<-FALSE}
	}
	if (Ties==TRUE) {
		TieCount<-matrix(0, ((TieMax-TieMin)+1), 3)
		TieCount[,1]<-TieMin:TieMax
		for (j in 1:length(TieList)) {
			TV<-TieList[j]
			TieCount[match(TV, TieCount[,1]),2]<-TieCount[match(TV, TieCount[,1]),2]+1
		}
		for (i in 1:nrow(TieCount)) {
			TieCount[i,3]<-TieCount[i,2]*TieCount[i,1]*(TieCount[i,1]-1)*(2*TieCount[i,1]+5)	
		}
		Corr<-sum(TieCount[,3])
	}

	#Calculating Sigma and Z
	{if (Ties==FALSE) {SigS<-sqrt((N*(N-1)*(2*N+5))/18)}
	else {SigS<-sqrt(((N*(N-1)*(2*N+5))-Corr)/18)}
	}

	{if ((P-M)>0){ZS<-((P-M)-1)/SigS}
	else if ((P-M)<0){ZS<-((P-M)+1)/SigS}
	else {ZS<-0}
	}

	Q<-pnorm(abs(ZS))
	{if (N>10) {Res[6]<-round(2*(1-Q), digits=3)}
	else {Res[6]<-paste("Too few values, consult Helsel and Hirsch (2002, App. B8) for S = ", (P-M), ", n = ", N, sep="")}
	}

	#Plotting data
	if (Plot==TRUE) {
		plot(Input[,1], Input[,2], type="p", main="Kendall-Theil robust line fitting", sub=paste("Significance of correlation p = ", round(Res[6], digits=3), sep=""), xlab="Independent variable", ylab="Dependent variable", col="cornflowerblue", pch=16)
		curve(as.numeric(Res[1])*x+as.numeric(Res[4]), add=TRUE, col="darkgreen", lwd=2)
	}
	
	#Calculating R squared
	ybar<-mean(Input[,2])

	TotVar<-vector()
	for (i in 1:nrow(Input)) {
		TotVar<-append(TotVar, ((Input[i,2]-ybar)^2))
	}
	SS.total<-sum(TotVar)

	TotRes<-vector()
	for (i in 1:nrow(Input)) {
		TotRes<-append(TotRes, ((Input[i,2]-(as.numeric(Res[1])*Input[i,1]+as.numeric(Res[4])))^2))
	}
	SS.residual<-sum(TotRes)

	Res[7]<-1-(SS.residual/SS.total)

	#Return results
	return(Res)
}

#########################################################################
# Function to calculate the confidence band around a linear regression. #
# Necessary input variables:                                            #
#   a: Slope of the regression line.                                    #
#      *numeric (real)*                                                 #
#   b: Intercept of the regression line.                                #
#      *numeric (real)*                                                 #
#   x.original: Original values of the independent variable.            #
#               *vector*                                                #
#   x.fit: Independent variable values for which to fit the...          #
#          confidence band.                                             #
#          *vector*                                                     #
#   y.original: Original values of the dependent variable.              #
#               *vector*                                                #
#   y.fitted: Fitted values for the dependent variable using the model. #
#             *vector*                                                  #
#   conf: Desired level of confidence.                                  #
#         *numeric (real)*                                              #
#         default=0.95                                                  #
# Output data: Confidence band around regression line.                  #
# Input dataset: Coefficients of linear regression.                     #
#########################################################################

Linear.Confidence<-function (a, b, x.original, x.fit, y.original, y.fitted, conf=0.95) {
	#Data consistency tests
	if (length(y.original)!=length(y.fitted)) {stop("y.original and y.fitted must be of same length!")}
	
	#Calculate t-value
	t.val<-qt(p=1-((1-conf)/2), df=length(y.original)-2)
	
	#Calculate standard error of regression
	mse<-sqrt(sum((y.original-y.fitted)^2, na.rm=TRUE)/(length(x.original)-2))
	se.ci<-mse*sqrt((1/length(x.original))+(x.fit-mean(x.original, na.rm=TRUE))^2/sum((x.original-mean(x.original, na.rm=TRUE))^2, na.rm=TRUE))
	se.pi<-mse*sqrt(1+(1/length(x.original))+(x.fit-mean(x.original, na.rm=TRUE))^2/sum((x.original-mean(x.original, na.rm=TRUE))^2, na.rm=TRUE))
	
	#Calculate fitted y-values over entire range
	y.fit<-a*x.fit+b

	#Calculate confidence band
	slope.upper<-y.fit+t.val*se.ci
	slope.lower<-y.fit-t.val*se.ci
	prediction.upper<-y.fit+t.val*se.pi
	prediction.lower<-y.fit-t.val*se.pi
	
	#Coerce data for export and export results
	Res<-list()
	Res$Input<-cbind(x.original, y.original, y.fitted)
	Res$Output<-cbind(x.fit, y.fit, slope.upper, slope.lower, prediction.upper, prediction.lower)
	colnames(Res$Input)<-c("X.data", "Y.data", "Y.model")
	colnames(Res$Output)<-c("X.newdata", "Y.newdata", "CI.upper", "CI.lower", "PI.upper", "PI.lower")
	return(Res)
}

#--------------------------------------------

#Examples
#lm1<-Mult.lm2(Y=ModII[,3], X=ModII[,1:2])
#lm1
#lm2<-Mult.lm2(Y=ModII[,3], X=ModII[,1:2], boot.n=999, tails="upper")
#lm2

#lm3<-Kendall(Y=ModII[,3], X=ModII[,1])
#lm3

#mod<-lm(Sizes[,"Length"] ~ Sizes[,"Salinity"])
#Conf<-Linear.Confidence(a=coef(mod)[2], b=coef(mod)[1], x.original<-Sizes[,"Salinity"], x.fit=seq(from=0, to=30, by=0.5), y.original=Sizes[,"Length"], y.fitted=coef(mod)[2]*Sizes[,"Salinity"]+coef(mod)[1])
#plot(Conf$Input[,"X.data"], Conf$Input[,"Y.data"], xlim=c(0, 30), ylim=c(1, 15))
#lines(Conf$Output[,"X.newdata"], Conf$Output[,"Y.newdata"], col="black", lwd=3)
#lines(Conf$Output[,"X.newdata"], Conf$Output[,"CI.upper"], col="blue", lwd=3)
#lines(Conf$Output[,"X.newdata"], Conf$Output[,"CI.lower"], col="blue", lwd=3)

#--------------------------------------------
#--------------------------------------------
#Version History
#1.0	Finished program
#2.0	Added functions Mult.lm2 and Kendall
#2.0.1	Fixed the Linear.Confidence function for extrapolation and fixed equation errors therin
#2.0.2	Added prediction interval estimation to Linear.Confidence
#
#Note: Function Kendall was taken from Kendall_Theil_Regression_Function.r and the original...
#	version history is given below. As of RegressionTools_Functions.r v. 2.0,...
#	Kendall_Theil_Regression_Function.r is deprecated.
##1.1	Adds option whether or not results should be exported as .txt file
##1.2	Adds option whether or not to import dataset, now calculates R2 value for regression line
##1.3	Added some security checks for input data structure and fixed problems for small datasets
##1.4	A bug was fixed that could under circumstances lead to a miscalculation of the number of ties
##1.4.1	Results are printed onto display by default
##1.5	Added possibility to disable default plotting
##1.6	Removed import/export functionality, streamlined code
#--------------------------------------------
#--------------------------------------------
