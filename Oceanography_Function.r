#Collection of functions for oceanographic analyses.
#Input data set:
#	-CondToSal: Table containing measurements of temperature, pressure (water depth) and conductivity
#
#Further reading:
#	Unesco/ICES/SCOR/IAPSO Joint Panel on Oceanographic Tables and Standards (1981) Background...
#		Papers and Supporting Data on the Practical Salinity Scale 1978. Unesco technical papers...
#		in marine science 37, (Unesco: Paris), 144 pp.
#	Schulz, H. D. and Zabel, M. (eds) (2006) Marine Geochemistry. 2nd ed., 574 pp. (Springer: Berlin,...
#		 Heidelberg). doi:10.1007/3-540-32144-6
#	Kölling, M. and Feseker, T. Alkalinity Titration. http://www.sedgeochem.uni-bremen.de/alkalinity.html


#Author: Manuel Weinkauf (Manuel.Weinkauf@unige.ch)
#Version: 2.0
#Date: 15 February 2017

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.#
#To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/.                   #
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#Creation of test dataset
#set.seed(12)
#T<-runif(10, min=-4, max=42)
#P<-runif(10, min=0, max=100)
#C<-rnorm(10, mean=50, sd=8)
#Conduct<-as.data.frame(cbind(T, P, C))
#colnames(Conduct)<-c("Temperature", "Pressure", "Conductivity")

#Alk<-matrix(c(8.02, 1.32, 0.1, 4.3, 50, 4.95, 22, 0.001, 3.912, 5), 2, 5, byrow=TRUE)
#colnames(Alk)<-c("ph.Start", "V.HCL", "C.HCL", "pH.End", "V.H2O")

#**************************************************************************************
#Setting working directory
#setwd("C:/R_TestData")

#########################################################################
# Function to convert conductivity into salinity.                       #
# Necessary input variables:						#
#   Temp: Temperature data.                                             #
#         *numeric (real)*                                              #
#   Press: Pressure data.                                               #
#         *numeric (real)*                                              #
#   Cond: Conductivity data.                                            #
#         *numeric (real)*                                              #
# Output data: A matrix containing the original data and associated...  #
#              salinity values.                                         #
# Input dataset: Temperature (in °C), pressure (in dbar), and...        #
#                conductivity (in mS/cm).                               #
#########################################################################

CondToSal<-function (Temp, Press, Cond) {
	#Test data integrity
	if (length(Temp)!=length(Press) | length(Temp)!=length(Cond)) {stop("All elemnts need to be of same length!")}
	
	#Calculate auxilliary variables
	Aux<-list()
	Aux$R<-Cond/42.914
	Aux$rT<-0.6766097+0.0200564*Temp+0.0001104259*Temp^2+-0.00000069698*Temp^3+0.0000000010031*Temp^4
	Aux$RP<-1+(0.0000207*Press+-0.000000000637*Press^2+0.000000000000003989*Press^3)/(1+0.03426*Temp+0.0004464*Temp^2+0.4215*Aux$R+-0.003107*Temp*Aux$R)
	Aux$RT<-Aux$R/(Aux$rT*Aux$RP)
	
	#Set up constants
	Const<-list()
	Const$a<-c(0.008, -0.1692, 25.3851, 14.0941, -7.0261, 2.7081)
	Const$b<-c(0.0005, -0.0056, -0.0066, -0.0375, 0.0636, -0.0144)
	
	#Calculate salinites
	a<-b<-list()
	for (i in 1:length(Const$a)) {
		a[[i]]<-Const$a[i]*Aux$RT^((i-1)/2)
	}
	for (i in 1:length(Const$b)) {
		b[[i]]<-Const$b[i]*Aux$RT^((i-1)/2)
	}
	tt<-(Temp-15)/(1+0.0162*(Temp-15))
	Sal<-mapply(sum, a[[1]], a[[2]], a[[3]], a[[4]], a[[5]], a[[6]])+tt*mapply(sum, b[[1]], b[[2]], b[[3]], b[[4]], b[[5]], b[[6]])
	if (any(Temp<2 | Temp>35)) {
		warning("Some temperatures were out of range (-2°C < T < 35°C)! NA's produced.")
		Sal[which(Temp<(-2))]<-NA
		Sal[which(Temp>35)]<-NA
	}
	
	#Return results
	Res<-matrix(NA, length(Temp), 4)
	Res[,1]<-Temp
	Res[,2]<-Press
	Res[,3]<-Cond
	Res[,4]<-Sal
	colnames(Res)<-c("Temperature.degC", "Pressure.dbar", "Conductivity.mS.cm", "Salinity")
	return(Res)
}

#########################################################################
# Function to calculate alkalinity from titration data.                 #
# Necessary input variables:						#
#   ph.Start: Initial pH.                                               #
#             *numeric (real)*                                          #
#   ph.End: Final pH.                                                   #
#           *numeric (real)*                                            #
#   V.H2O: Volume of titrated water.                                    #
#          *numeric (real)*                                             #
#   V.HCl: Volume of added acid.                                        #
#          *numeric (real)*                                             #
#   C.HCl: Concentration of added acid.                                 #
#          *numeric (real)*                                             #
#   f: Fugacity of the water.                                           #
#      *numeric (real)*                                                 #
#      default=0.8                                                      #
# Output data: A matrix containing the original data and associated...  #
#              alkalinity values (in mmol/l).                           #
# Input dataset: Start and end pH of the water, volume of water...      #
#                used in the titration (in ml), volume of acid added...	#
#                to the water (in ml), concentration of the acid (in... #
#                mmol/l).                                               #
#########################################################################

Alkalinity<-function(pH.Start, pH.End, V.H2O, V.HCl, C.HCl, f=0.8) {
	#Test data consistency
	V.lengths<-c(length(pH.Start), length(pH.End), length(V.H2O), length(V.HCl), length(C.HCl))
	if (min(V.lengths)!=max(V.lengths)) {stop("All vectors (except f) must be of equal length!")}
	
	#Calculate alkalinity
	ALK<-(((V.HCl/1000) * C.HCl) - 10**-pH.End * ((V.H2O/1000) + (V.HCl/1000))/f + 10**-pH.Start * (V.H2O/1000)/f)/(V.H2O/1000)
	ALK<-ALK*1000
	
	#Return results
	Res<-matrix(NA, length(pH.Start), 6)
	colnames(Res)<-c("pH.Start", "pH.End", "V.H2O.ml", "V.HCl.ml", "C.HCl.mmol.l", "Alk.mmol.l")
	Res[,"pH.Start"]<-pH.Start
	Res[,"pH.End"]<-pH.End
	Res[,"V.H2O.ml"]<-V.H2O
	Res[,"V.HCl.ml"]<-V.HCl
	Res[,"C.HCl.mmol.l"]<-C.HCl
	Res[,"Alk.mmol.l"]<-ALK
	return(Res)
}

#--------------------------------------------

#Example
#CondToSal(Conduct[,"Temperature"], Conduct[,"Pressure"], Conduct[,"Conductivity"])
#Alkalinity(Alk[,"ph.Start"], Alk[,"pH.End"], Alk[,"V.H2O"], Alk[,"V.HCL"], Alk[,"C.HCL"])

#--------------------------------------------
#--------------------------------------------
#Version History
#1.0	Finished Program
#2.0	Added function Alkalinity
#--------------------------------------------
#--------------------------------------------
