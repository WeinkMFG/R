#Function to find parameters of a given curve function

#Further reading: Chrysafi, Loucas and Gordon, Sheldon (2006) "On the curvature function: Where does...
#                 a curve bend the fastest?" Mathematics and Computer Education.

#Author: Manuel Weinkauf (Manuel.Weinkauf@unige.ch)
#Version: 1.0
#Date: 14 August 2014

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.#
#To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/.                   #
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#**************************************************************************************
#Setting working directory
#setwd("C:/R_TestData")

#########################################################################
# Function to find the maximum curvature point of a function            #
# Necessary input variables:                                            #
#    Function: Input function for which point of maximum curvature...   #
#              should be calculated.                                    #
#              *symbolic function*                                      #
#              NOTE: the variable to which to deviate is assumed to...  #
#                    be called "x".                                     #
#    Limes: Range along the x-axis for which results should be plotted. #
#           *numeric (real)* of length 2                                #
#           default=-10:10                                              #
#    N: Number of x-values to interpolate for plotting.                 #
#       *numeric (real)*                                                #
#       default: 100 times x-range                                      #
#    Iter: Precision of iteratively finding the maxima. The smaller...  #
#          this value, the more precise the result.                     #
#          *real*                                                       #
#          default=x-range/10000                                        #
# Output data: The curvature function and point of maximum flecture.    #
# Input data: A symbolic representation of the equation for which the...#
#             point of maximum curvature should be found.               #
#########################################################################

Curvature<-function (Function, Limes=c(-10, 10), N=(diff(Limes)*100), Iter=diff(Limes)/10000) {
	#Calculate first and second derivative
	Deriv1<-D(parse(text=Function), "x")
	Deriv2<-D(Deriv1, "x")
	##Reformat to continuous string
	Temp<-paste(deparse(Deriv1))
	Deriv1<-paste(Temp[1:length(Temp)], collapse="")
	Temp<-paste(deparse(Deriv2))
	Deriv2<-paste(Temp[1:length(Temp)], collapse="")
	
	#Prepare plotting data
	##Original curve
	fx<-function (x) {eval(parse(text=Function))}
	##Curvature function
	Curve.Function<-paste("abs(", Deriv2, ")", "/((1 + (", Deriv1, ")^2)^(3/2)", ")", sep="")
	fc<-function (x) {eval(parse(text=Curve.Function))}
	
	#Find maxima of curve function
	##Function to locate approximate position of maxima
	which.peaks<-function(x) {which(diff(diff(x)>=0)<0)+1}
	
	##Create list of y-values within Limes
	X.Val<-seq(from=Limes[1], to=Limes[2], by=Iter*100)
	Y.Val<-vector(mode="numeric", length=length(X.Val))
	for (i in 1:(length(Y.Val))) {
		x<-X.Val[i]
		Y.Val[i]<-eval(parse(text=Curve.Function))
	}
	
	##Locate approximate position of maxima
	Max.Approx<-X.Val[which.peaks(Y.Val)]
	
	##Iteratively find exact maximum
	Max.Pos<-vector(mode="numeric", length=length(Max.Approx))
	for (i in 1:(length(Max.Approx))) {
		###Calculate start value
		x<-Max.Approx[i]
		Y.Old<-eval(parse(text=Curve.Function))
		
		###Find direction to go
		x<-Max.Approx[i]-Iter
		Y.Left<-eval(parse(text=Curve.Function))
		x<-Max.Approx[i]+Iter
		Y.Right<-eval(parse(text=Curve.Function))
		{if (Y.Left>Y.Old) {Left<-TRUE}
		else if (Y.Right>Y.Old) {Left<-FALSE}
		else {Left<-NA; Max.Pos[i]<-Max.Approx[i]}}
		###Begin iteration
		if (!is.na(Left)) {if (Left==TRUE) {X.Val<-seq(from=Max.Approx[i], to=Max.Approx[i]-Iter*100, by=-Iter)}
		else {X.Val<-seq(from=Max.Approx[i], to=Max.Approx[i]+Iter*100, by=Iter)}}
		Y.Val<-vector(mode="numeric", length=length(X.Val))
		for (j in 1:(length(Y.Val))) {
			x<-X.Val[j]
			Y.Val[j]<-eval(parse(text=Curve.Function))
		}
		Max.Pos[i]<-X.Val[which.peaks(Y.Val)]
	}
	
	#Plot results
	par(mar=c(4,4,1,4))
	curve(fx, from=Limes[1], to=Limes[2], n=N, lwd=2, col="black", xlab="x", ylab="y")
	par(new=TRUE)
	curve(fc, from=Limes[1], to=Limes[2], n=N, lwd=2, col="blue", axes=FALSE, xlab="", ylab="")
	axis(side=4, col.axis="blue", col="blue", line=0)
	mtext("Curvature function", side=4, line=3, col="blue")
	for (i in 1:(length(Max.Pos))) {
		abline(v=Max.Pos[i], lwd=1, lty=2, col="grey50")
	}
	legend("topright", lwd=c(2, 2, 1), lty=c(1, 1, 2), col=c("black", "blue", "grey50"), bg="white", cex=0.7, legend=c("Input function", "Curvature function", "Maximum curvature"))
	
	#Return results equations
	Res<-list()
	Res$Input.Function<-Function
	Res$Curvature.Function<-Curve.Function
	Res$Maximum.Curvature<-Max.Pos
	return(Res)
}

#--------------------------------------------

#Example
##Simple examples
#Curvature("x^2", Limes=c(-2, 2))#Quadratic function
#Curvature("x^3", Limes=c(-2, 2))#Cubic function
#Curvature("x^4", Limes=c(-2, 2))#Higher Polynomial function
#Curvature("3.67^x", Limes=c(-2, 2))#Power function
#Curvature("exp(1)^x", Limes=c(-2, 2))#Exponential function
#Curvature("log(x)", Limes=c(0, 3))#Logarithmic function
#Curvature("sin(x)", Limes=c(0, 7))#Trigonometric function
#Curvature("1/(1+2*exp(1)^(-3*x))", Limes=c(-3, 3))#Logistic function

#Biological examples
#Curvature("1.2/(0.8 + x)", Limes=c(0, 10))#Hyperbolic function
#Curvature("(1.2 * x)/(0.8 + x)", Limes=c(0, 10))#Michaelis-Menten function
#Curvature("(1.2 * x^2)/(0.8^2 + x^2)", Limes=c(0, 10))#Holling type III function
#Curvature("(1.2 * x^2)/(0.8 + -1.2*x + x^2)", Limes=c(-1, 6))#Holling type IV function
#Curvature("1.2*exp(1)^(-0.8*x)", Limes=c(0, 10))#Negative exponential function
#Curvature("1.2*(1-exp(1)^(-0.8*x))", Limes=c(0, 10))#Monomolecular function
#Curvature("1.2*x*exp(1)^(-0.8*x)", Limes=c(0, 10))#Ricker function
#Curvature("1.2/(1+((1.2/1.1)-1)*exp(1)^(-0.6*x))", Limes=c(-6, 10))#Logistic function (population ecology)
#Curvature("1.2*(1-exp(1)^(-0.8*(1-(2/3))*(x-0.1)))^(1/(1-(2/3)))", Limes=c(0, 20))#von Bertalanffy function
#Curvature("(1.2*x)/(0.8+x^1.2)", Limes=c(0, 10))#Shepherd function
#Curvature("(1.2*x)/((0.8+x)^1.2)", Limes=c(0, 10))#Hassel function

#--------------------------------------------
#--------------------------------------------
#Version History
#1.0	Finished Program
#--------------------------------------------
#--------------------------------------------
