#Bootstrapping of confidence intervals for manual permutated data
#Input data set for manual permutation
#	-A number of data values, obtained from manual Monte-Carlo permutation of specimens
#	-Suitable when several specimens had to be measured together and results are...
#		only means of all specimens measured together, so that only one...
#		value for each sample exists and no classical approach of CI calculation...
#		can be applied
#	-You will have to give the minimum and maximum number of specimens in the bucket,...
#		which refers to the minimal and maximal number of specimens that have been...
#		been weighed together globally, so that you get confidence intervals for the...
#		whole range of occuring sample sizes.
#Input data set for randomized correlation
#	-Several measurement values including errors or confidence intervals obtained with 'ManPerm()'
#	-First column contains variable to correlate with
#	-Second column contains variable that has measurement error, with said error in third and...
#		fourth column.
#Cite as: Weinkauf, M. F. G., Moller, T., Koch, M. C., Kucera, M. (2013) "Calcification...
#		intensity in planktonic Foraminifera reflects ambient conditions irrespective...
#		of environmental stress". Biogeosciences 10 (10): 6639-6655

#Author: Manuel Weinkauf (Manuel.Weinkauf@unige.ch)
#Version: 1.4
#Date: 13 March 2014

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.#
#To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/.                   #
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#Creation of test dataset for ManPerm()
PermDat<-sample(1:100, size=30)


#**************************************************************************************
#Setting working directory
#setwd("C:/R_TestData")

#########################################################################
# Function to apply basic bootstrapping on manual permutation data.     #
# Necessary packages: tcltk                                             #
# Necessary input variables:                                            #
#    Input: Name of input data.                                         #
#           *vector*                                                    #
#    Start: Minimum number of specimens in bucket (i.e. smallest...     #
#           occurring sample size).                                     #
#           *numeric* (integer)                                         #
#           default=1                                                   #
#    Stop: Maximum number of specimens in bucket (i.e. largest...       #
#          occurring sample size).                                      #
#          *numeric* (integer)                                          #
#    Replicates: Number of desired replications.                        #
#                *numeric* (integer)                                    #
#                default=10000                                          #
# Output data: 95% confidence intervals and errors around mean of...    #
#              for different sample sizes.                              #
# Input dataset: Vector containing measurements of several manual...    #
#                Monte-Carlo permutations that are in first...          #
#                approximation a set of individual measurement values.  #
#########################################################################

#Loading necessary packages
require(tcltk)

ManPerm<-function (Input, Start=1, Stop, Replicates=10000) {
	#Test data consistency
	if (Stop<Start) {stop("Max. sample size cannot be smaller than min. sample size!")}
	
	#Creation of function for 'statistic'
	my_mean<-function(data, indices){
		d<-data[indices,]
		mean(d, na.rm=TRUE)
	}

	#Creating output matrix
	RowNum<-Stop-(Start-1)
	CI<-matrix(NA, RowNum, 6)
	colnames(CI)<-c("Number.of.Specimens", "CI.Low", "CI.High", "Mean", "Error.Low", "Error.High")
	CI[,1]<-c(Start:Stop)

	#Bootstrapping
	pb<-tkProgressBar(title="Progress", min=Start, max=Stop, width=300)
	for (i in Start:Stop) {
		SampDataset<-vector(mode="numeric", length=i)
		BootDataset<-vector(mode="numeric", length=Replicates)
	
		for (j in 1:Replicates) {
			for (k in 1:i) {
				Pos<-sample(length(Input), i, replace=TRUE)
				SampDataset[k]<-Input[Pos[k]]
			}
			BootDataset[j]<-mean(SampDataset)
		}
		
		Quant<-quantile(BootDataset, probs=c(0.025, 0.975), type=8)
		
		CI[i,2]<-Quant[1]
		CI[i,3]<-Quant[2]
		CI[i,4]<-mean(BootDataset, na.rm=TRUE)
		Sys.sleep(0.1)
   		setTkProgressBar(pb, i, label=paste(round((i-Start)/(Stop-Start)*100, digits=0), "% of randomisation completed"))
	}

	#Calculations
	for (i in 1:nrow(CI)) {
		CI[i,5]<-CI[i,4]-CI[i,2]
		CI[i,6]<-CI[i,3]-CI[i,4]
	}
	
	#Writing output file
	return(CI)
}

#########################################################################
# Function to calculate correlation of data for which errors or CI's... #
#	were calculated by ManPerm()                                    #
# Necessary input packages: tcltk                                       #
# Necessary input variables:                                            #
#    Input 1: Values, and lower and upper confidence interval of...     #
#             dataset 1.                                                #
#             *list*                                                    #
#    Input 2: Values of dataset 2 (to be correlated with).              #
#             *vector*                                                  #
#    Interval: Are the data to be correlated scalar or ordinal values.  #
#              *logical*                                                #
#              TRUE: All data are scalar                                #
#              FALSE: At least one column is ordinal                    #
#              default=TRUE                                             #
#    Replicates: Number of desired replications.                        #
#                *integer*                                              #
#                default=10000                                          #
# Output data: Set of correlation coefficients and p-values for...      #
#              random correlation of data (regarding variability)...    #
#              with given data, using pairwise deletion.                #
# Input datasets: Values and lower and upper confidence interval of...  #
#                 a datasets. Values of another dataset to correlate... #
#                 with.                                                 #
#########################################################################

#Loading necessary packages
require(tcltk)

RandCorr<-function (Input1, Input2, Interval=TRUE, Replicates=10000) {
	#Check data integrity
	if (!is.list(Input1)) {stop("Input1 must be a list!")}
	if (!is.vector(Input2)) {stop("Input2 must be a vector!")}
	if (length(Input1)<3) {stop("Input1 must have three elements: Value, lower CI, and upper CI!")}
	if (length(Input1)>3) {warning("Input1 has more than three elements, using first three only!")}
	if (length(Input2)!=length(Input1[[1]])) {stop("Input datasets are not of the same length!")}
	if (!all(Input1[[1]]>=Input1[[2]])) {stop("In Input1, lower confidence interval must be smaller than value!")}
	if (!all(Input1[[1]]<=Input1[[3]])) {stop("In Input1, upper confidence interval must be larger than value!")}
	
	#Define method to be used
	Shap1<-shapiro.test(Input1[[1]])
	Shap2<-shapiro.test(Input2)
	{if (Shap1$p.value<0.05 | Shap2$p.value<0.05) {Meth<-"spearman"}
		else if (Interval==FALSE) {Meth<-"spearman"}
		else {Meth<-"pearson"}
	}

	#Creating output matrix
	Res<-matrix(NA, Replicates+1, 2)
	{if (Meth=="pearson") {colnames(Res)<-c("Pearson", "p")}
		else {colnames(Res)<-c("Spearman", "p")}
	}
	rownames(Res)<-c("Original data", paste("Randomisation", 1:Replicates, sep="_"))

	#Testing original data without variability
	Corr<-cor.test(Input2, Input1[[1]], alternative="two.sided", method=Meth, use="pairwise.complete.obs")
	Res[1,1]<-Corr$estimate
	Res[1,2]<-Corr$p.value

	#Setting up resampling matrix
	Resample<-matrix(NA, length(Input1[[1]]), 2)
	Resample[,1]<-Input2

	#Performing randomized correlation
	pb<-tkProgressBar(title="Progress", min=0, max=Replicates, width=300)
	for (i in 2:(Replicates+1)) {
		#Resample values within confidence interval
		for (j in 1:length(Input1[[1]])) {
			Resample[j,2]<-runif(1, min=Input1[[2]], max=Input1[[3]])
		}

		#Perform correlation analysis
		Corr<-cor.test(Resample[,1], Resample[,2], alternative="two.sided", method=Meth, use="pairwise.complete.obs")
		Res[i,1]<-Corr$estimate
		Res[i,2]<-Corr$p.value
		Sys.sleep(0.1)
   		setTkProgressBar(pb, i, label=paste(round(i/Replicates*100, 0),"% of randomisation completed"))
	}

	#Calculating means
	MeanVal<-vector(mode="numeric", length=2)
	MeanVal[1]<-mean(Res[-1,1])
	MeanVal[2]<-mean(Res[-1,2])
	Res<-rbind(Res, MeanVal)

	#Plot histogram
	hist(Res[,2], main=paste("Histogram for", Replicates, "Randomisations"), xlab=expression(paste(italic("p"),"-value")), col="grey")
	mtext(paste("Mean significance =", Res[(Replicates+2),2]), side=4)
	
	#Return results
	return(Res)
}


#--------------------------------------------

#Example

#Permute<-ManPerm(PermDat, Start=1, Stop=50, Replicates=1000)

#CorrDat1<-list()
#CorrDat1$Values<-Permute[,"Mean"]
#CorrDat1$Low.CI<-Permute[,"CI.Low"]
#CorrDat1$Up.CI<-Permute[,"CI.High"]
#CorrDat2<-sample(1:100, size=length(CorrDat1$Values))

#Correlate<-RandCorr(CorrDat1, CorrDat2, Replicates=100)

#--------------------------------------------
#--------------------------------------------
#Version History
#1.0	Finished Program
#1.1	Progress bar added
#1.2	RandCorr allows to define NA strings
#1.2.1	ManPerm allows Start to be equal to Stop
#1.3	Added options not to import or export data
#1.4	Removed import/output functions entirely
#	Generally improved the code
#--------------------------------------------
#--------------------------------------------