#Calculating biodiversity and associated confidence limits and jackknifing for species richness.
#Input data set:
#	-A matrix of species in samples (presence/absence or count data), species in columns and samples in rows.
#	-For the species richness and diversity indices, samples are treated as independent of each other, i.e. as...
#		samples from different time slices or localities. The purpose is to allow a comparison between these...
#		localities.
#	-For jackknifing, samples are treated as replicates on a meaningful scale within one region (e.g. plots). The...
#		purpose is to estimate the species richness in the region based on repliacted sampling of different...
#		localities.
#Further reading:
#	Chao, A., Hwang, W.-H., Chen, Y.-C., and Kuo, C.-Y. (2000) Estimating the number of shared species in...
#		two communities. Statistica Sinica 10 (1): 227-256. http://www.jstor.org/stable/24306714
#	Chao, A. and Shen, T.-J. (2003) Nonparametric estimation of Shannon's index of diversity when there...
#		are unseen species in sample. Environmental and Ecological Statistics 10: 429-443.
#	Hammer, Ø. and Harper, D. (2006) Paleontological Data Analysis. 351 pp., Blackwell Publishing: Malden,...
#		Oxford, and Carlton.
#	Heslop, D., De Schepper, S., and Proske, U. (2011) "Diagnosing the uncertainty of taxa relative...
#		abundances derived from count data". Marine Micropaleontology 79: 114--20.
#	Hill, M. O. (1973) Diversity and evenness: A unifying notation and its consequences. Ecology 54:...
#		427-432. doi:10.2307/1934352
#	Horvitz, D. G. and Thompson, D. J. (1952) A generalization of sampling without replacement from a...
#		finite universe. Journal of the American Statistical Association 47: 663-685.
#	Hurlbert, S. H. (1971) The nonconcept of species diversity: A critique and alternative parameters....
#		Ecology 52: 577-586.
#	Laakso, M. and Taagepera, R. (1979) "Effective" number of parties: A measure with application to West...
#		Europe. Comparative Political Studies 12: 3-25.
#	Legendre, P. and Legendre, L. (2012) Numerical Ecology. 3rd. ed., Developments in Environmental Modelling...
#		no. 24, 990 pp., Elsevier: Amsterdam and Oxford.
#	Leti, G. (1983) Statistica descrittiva. 941 pp., Il Mulino.
#	Lloyd, M. and Ghelardi, R. J. (1964) A table for calculating the 'equitability' component of species...
#		diversity. Journal of Animal Ecology 33: 217-225.
#	Magurran, A. E. (1988) Ecological diversity and its measurement. 179 pp., Croom Helm: Cambridge.
#	Margalef, R. (1958) Information theory in ecology. General Systematics 3: 36-71.
#	Menhinick, E. F. (1964) A comparison of some species-individuals diversity indices applied to samples of...
#		field insects. Ecology 45: 859-861. doi:10.2307/1934933
#	Patten, B. C. (1962) Species diversity in net phytoplankton of Raritan Bay. Journal of Marine Research 20:...
#		57-75.
#	Pesenti, N., Quatto, P., and Ripamonti, E. (2017) Bootstrap confidence intervals for biodiversity measures...
#		based on Gini index and entropy. Quality and Quantity 51: 847-858. doi:10.1007/s11135-016-0443-x
#	Pielou, E. C. (1977) Mathematical Ecology. 385 pp., John Wiley & Sons: New York, London, Sydney, and...
#		Toronto.
#	Shannon, C. E. and Weaver, W. (1949) The Mathematical Theory of Communication. 144 pp., University of...
#		Illinois Press: Urbana and Chicago.
#	Simpson, E. H. (1949) Measurement of diversity. Nature 163: 688. doi:10.1038/163688a0
#	Smith, E. P. and van Belle, G. (1984) Nonparametric estimation of species richness. Biometrics 40 (1): 119-129.
#	https://www.statsdirect.com/help/nonparametric_methods/diversity.htm
#	http://www.dataanalytics.org.uk/Publications/S4E2e%20Support/exercises/Comparing%20shannon%20diversity.htm

#Author: Manuel Weinkauf (Manuel.Weinkauf@unige.ch)
#Version: 1.4.1
#Date: 31 August 2018

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.#
#To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/.                   #
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#Creation of test dataset
#Data.Count<-matrix(sample(0:20, size=400, replace=TRUE, prob=c(0.5, rep(0.5/20, 20))), 40, 10)
#colnames(Data.Count)<-paste("Spec", 1:10, sep=".")
#rownames(Data.Count)<-paste("Sample", 1:40, sep=".")
#Data.Bin<-matrix(sample(x=c(0,1), size=200, replace=TRUE), 20, 10)
#colnames(Data.Bin)<-paste("Spec", 1:10, sep=".")
#rownames(Data.Bin)<-paste("Sample", 1:20, sep=".")
#Data.Jack<-matrix(sample(x=c(0, 1), size=50, replace=TRUE, prob=c(0.8, 0.2)), 5, 10)
#colnames(Data.Jack)<-paste("Spec", 1:10, sep=".")
#rownames(Data.Jack)<-paste("Replicate", 1:5, sep=".")
#Data.Chao<-sample(0:15, size=20, replace=TRUE)
#set.seed(1000)
#Prop1<-matrix(c(1402, 78, 623, 19, 1232, 32), 3, 2, byrow=TRUE)
#colnames(Prop1)<-c("N", "n.obs")
#rownames(Prop1)<-paste("Sam", 1:3, sep=".")
#Prop2<-Prop
#Prop2[,2]<-Prop2[,2]/Prop2[,1]
#Abund<-cbind(sample(40, size=10, replace=TRUE, prob=c(0.5, rep(0.5/39, 39))), sample(40, size=10, replace=TRUE, prob=c(0.5, rep(0.5/39, 39))), sample(40, size=10, replace=TRUE, prob=c(0.5, rep(0.5/39, 39))))
#colnames(Abund)<-paste("Species", 1:3, sep=".")
#rownames(Abund)<-paste("Samples", 1:10, sep=".")

#**************************************************************************************
#Setting working directory
#setwd("C:/R_TestData")

#########################################################################
# Function to calculate species richness indices based on presence...   #
#   absence data.                                                       #
# Necessary input variables:                                            #
#   Dat: Dataset containing species presence/absence or count data...   #
#        (in columns) in samples (in rows).                             #
#        *matrix* or *vector*                                           #
#   Method: Species richness index to calculate. This can be...         #
#           "Richness" (for simple species richness), "Menhinick"...    #
#           (for Menhinick's richness index) or "Margalef" (for...      #
#           Margalef's richness index).                                 #
#           *character*, either of "Richness", "Menhinick", or...       #
#              "Margalef"                                               #
#           default="Richness"                                          #
#   N: Sample size (number of specimens per sample). If Dat is...       #
#      inferred to be count data (at least one value>1) it is...        #
#      calculated automatically. Otherwise, it must be given as a...    #
#      vector of the length nrow(Dat). When N!=NULL for count data it...#
#      will be assumed that N is correct and be used as sample size...  #
#      instead of the internal calculation. Only necessary when...      #
#      Method is "Menhinick" or "Margalef".                             #
#      *vector* of length nrow(Dat) or NULL                             #
#      default=NULL                                                     #
#   Plot: Should results be plotted?                                    #
#         *logical*                                                     #
#         default=TRUE                                                  #
# Output data: Vector of species richness indices.                      #
# Input dataset: Count or presence/absence data of species in columns...#
#                vs. samples in rows.                                   #
#########################################################################

Richness<-function (Dat, Method="Richness", N=NULL, Plot=TRUE) {
	#Test data consistency
	{if (class(Dat)=="numeric") {Length<-1}
	else {Length<-nrow(Dat)}}
	if (max(Dat)>1) {warning("Data seem to be count data, not presence/absence. Diversity indices may be more useful for those data.")}
	if (Method!="Richness" & Method!="Menhinick" & Method!="Margalef") {stop("Method must be either of 'Richness'. 'Menhinick', or 'Margalef'!")}
	if (!is.null(N) & length(N)!=Length) {stop("N must be of length nrow(Dat)!")}
	if (Plot==TRUE & Length==1) {warning("Dat only contains one sample, no plot will be produced.")}
	if (is.data.frame(Dat)==TRUE) {Dat<-as.matrix(Dat)}
	
	#Calculate sample size
	if (max(Dat)>1 & is.null(N)) {N<-apply(Dat, 1, sum)}
	
	#Set up results object
	Res<-vector(mode="numeric", length=Length)
	
	#Calculate richness index
	for (i in 1:Length) {
		{if (class(Dat)=="numeric") {Temp<-Dat}
		else {Temp<-Dat[i,]}}
		Sp.N<-length(which(Temp>0))
		
		{if (Method=="Richness") {Res[i]<-Sp.N}
		else if (Method=="Menhinick") {Res[i]<-Sp.N/sqrt(N[i])}
		else {Res[i]<-(Sp.N-1)/log(N[i], base=exp(1))}}
	}
	
	#Plot results
	if (Plot==TRUE & Length>1) {
		win.graph((Length/2.5), 10, 10)
		plot(Res, type="b", pch=16, lwd=2, col="cornflowerblue", xlab="Sample", ylab="Species richness")
		{if (Method=="Richness") {title("Species Richness Index")}
		else if (Method=="Menhinick") {title("Menhinick's Richness Index")}
		else {title("Margalef's Richness Index")}}
	}
	
	#Return results
	if (!is.null(rownames(Dat))) {names(Res)<-rownames(Dat)}
	return(Res)
}

#########################################################################
# First-order jackknifing to estimate species richness from...          #
#   replicated sampling.                                                #
# Necessary input variables:                                            #
#   Input: Name of input matrix, missing values are not supported.      #
#          *matrix*                                                     #
# Output data: Estimated species richness including 95% confidence...   #
#              interval of estimate.                                    #
# Input dataset: Matrix object with count or presence-absence data of...#
#                species in columns vs. samples in rows. Each sample... #
#                must be a replication on a meaningful scale within...  #
#                the area of interest.                                  #
#########################################################################

Jack<-function (Input) {
	#Clean data from species, that are not present in any sample
	Spec.Present<-apply(Input, 2, sum)
	if (any(Spec.Present==0)) {
		Input<-Input[,-which(Spec.Present==0)]
		warning("Some species do not occur in any sample, missings have been removed!")
	}
	
	#Convert data to presence absence
	if (any(Input>1)) {
		Input<-apply(Input, c(1, 2), function(x) {if (x>0) {x<-1} else {x<-0}})
	}
	
	#Calculate Sample parameters
	#Observed species richness
	S.obs<-ncol(Input)
	
	#Partial estimators
	PE<-vector(mode="numeric", length=nrow(Input))
	for (i in 1:(dim(Input)[1])) {
		T<-which(Input[-i,]==1, arr.ind=TRUE)
		PE[i]<-length(rle(T[,2])$values)
	}
	
	#Pseudo values
	E<-vector(mode="numeric", length=nrow(Input))
	for (i in 1:(nrow(Input))) {
		E[i]<-nrow(Input)*S.obs-(nrow(Input)-1)*PE[i]
	}
	
	#Calculate results
	S.est<-vector(mode="numeric", length=3)
	names(S.est)<-c("S.obs", "S.est", "95%.CI")
	S.est["S.obs"]<-S.obs
	S.est["S.est"]<-mean(E)
	S.est["95%.CI"]<-1.96*(sd(E)/sqrt(nrow(Input)))
	
	#Return results
	return(S.est)
}

#########################################################################
# Second-order jackknifing to estimate species richness from...         #
#   replicated sampling.                                                #
# Necessary input variables:                                            #
#   Input: Name of input matrix, missing values are not supported.      #
#          *matrix*                                                     #
#   Replicates: How often should the random second order bias-...       #
#               correction be repeated?                                 #
#               *numeric (integer)*                                     #
#               default=1000                                            #
# Output data: Estimated species richness including 95% confidence...   #
#              interval of estimate.                                    #
# Input dataset: Matrix object with count or presence-absence data of...#
#                species in columns vs. samples in rows. Each sample... #
#                must be a replication on a meaningful scale within...  #
#                the area of interest.                                  #
#########################################################################

#Load packages
require(tcltk)

Jack2<-function (Input, Replicates=1000) {
	#Clean data from species, that are not present in any sample
	Spec.Present<-apply(Input, 2, sum)
	if (any(Spec.Present==0)) {
		Input<-Input[,-which(Spec.Present==0)]
		warning("Some species do not occur in any sample, missings have been removed!")
	}
	
	#Convert data to presence absence
	if (any(Input>1)) {
		Input<-apply(Input, c(1, 2), function(x) {if (x>0) {x<-1} else {x<-0}})
	}
	
	#Calculate Sample parameters
	#Observed species richness
	S.obs<-ncol(Input)
	
	#Partial estimators
	Rerun<-list()
	Rerun$S.est<-vector(mode="numeric", length=Replicates)
	Rerun$CI95<-vector(mode="numeric", length=Replicates)
	#Set up progress bar
	pb<-tkProgressBar(title="Progress", min=1, max=Replicates, width = 300)
	for (k in 1:Replicates) {
		PE<-vector(mode="numeric", length=nrow(Input))
		for (i in 1:(dim(Input)[1])) {
			T<-which(Input[-i,]==1, arr.ind=TRUE)
			PE[i]<-length(rle(T[,2])$values)
		}
		{if (nrow(Input)%%2==0) {Seq<-sample(nrow(Input), size=nrow(Input), replace=FALSE)}
		else {Seq<-sample(nrow(Input), size=nrow(Input)-1, replace=FALSE)}}
		Seq1<-Seq[seq(from=1, to=(length(Seq)-1), by=2)]
		Seq2<-Seq[seq(from=2, to=length(Seq), by=2)]
		PE2<-vector(mode="numeric", length=length(Seq1))
		for (i in 1:(length(Seq1))) {
			T<-which(Input[-c(Seq1[i], Seq2[i]),]==1, arr.ind=TRUE)
			PE2[i]<-length(rle(T[,2])$values)
		}
	
		#Pseudo values
		E<-vector(mode="numeric", length=nrow(Input))
		for (i in 1:nrow(Input)) {
			E[i]<-nrow(Input)*S.obs-(nrow(Input)-1)*PE[i]
		}
		E2<-vector(mode="numeric", length=length(PE2))
		for (i in 1:length(PE2)) {
			E2[i]<-nrow(Input)*S.obs-(nrow(Input)-1)*PE2[i]
		}
	
		##Calculate rerun results
		Rerun$S.est[k]<-mean(c(E, E2))
		Rerun$CI95[k]<-1.96*(sd(c(E, E2))/sqrt(nrow(Input)))
		
		Sys.sleep(0.1)
		setTkProgressBar(pb, k, label=paste(round((k-1)/(Replicates-1)*100, 0),"% of reruns completed"))
	}
	
	#Calculate results
	S.est<-vector(mode="numeric", length=3)
	names(S.est)<-c("S.obs", "S.est", "95%.CI")
	S.est["S.obs"]<-S.obs
	S.est["S.est"]<-mean(Rerun$S.est)
	S.est["95%.CI"]<-mean(Rerun$CI95)
	
	#Return results
	return(S.est)
}

#########################################################################
# Chao estimate of species richness from single samples.                #
# Necessary input variables:                                            #
#   Input: Name of input vector, missing values are not supported.      #
#          *vector*                                                     #
#   k: Threshold to distinguish rare and abundant species. Default...   #
#      uses suggestion by Chao et al. (2000).                           #
#      *numeric (integer)*                                              #
#      default=10                                                       #
# Output data: Estimated species richness.                              #
# Input dataset: Species richness estimate.                             #
#########################################################################

Chao<-function (Input, k=10) {
	#Test data consistency
	if (!is.vector(Input)) {stop("Input must be vector!")}
	if (max(Input)==1) {stop("Count data needed, not presence absence!")}
	k<-round(k, digits=0)
	Input<-Input[which(Input>0)]
	
	#Calculate parameters
	#Number of rare and abundant species
	S.abun<-length(which(Input>k))
	S.rare<-length(which(Input<=k))
	
	#Number of species per abundance f1
	f1<-length(which(Input==1))
	
	#Sample coverage C
	C<-1-(length(which(Input==1))/sum(Input))
	
	#Rare species coverage C.rare
	C.vec<-vector(length=k, mode="numeric")
	for (i in 1:k) {
		C.vec[i]<-length(which(Input==i))
	}
	C.rare<-1-f1/sum(C.vec)
	if (is.nan(C.rare)) {C.rare<-1}
	
	#Gamma
	Numerator<-Denominator1<-Denominator2<-vector(length=k, mode="numeric")
	for (i in 1:k) {
		Numerator[i]<-i*(i-1)*length(which(Input==i))
		Denominator1[i]<-i*length(which(Input==i))
		Denominator2[i]<-i*length(which(Input==i))-1
	}
	Numerator<-sum(Numerator)
	Denominator1<-sum(Denominator1)
	Denominator2<-sum(Denominator2)
	Gamma<-max(c((S.rare/C)*(Numerator/(Denominator1*Denominator2))-1, 0))
	if(is.infinite(Gamma) | is.nan(Gamma)) {Gamma<-0}
	
	#Calculate species richness estimate
	Res<-vector(mode="numeric", length=2)
	names(Res)<-c("Observed.richness", "Estimated.richness")
	Res[1]<-length(which(Input>0))
	Res[2]<-S.abun+(S.rare/C.rare)+(f1/C.rare)*Gamma
	return(Res)
}

#########################################################################
# Function to calculate diversity indices based on count data. If...    #
#   applicable, functional diversity after Laakso and Taagepera...      #
#   (1979) (for Simpson index) or Leti (1983) (for Shannon--Wiener...   #
#   index) will be calculated as well.                                  #
# Necessary input variables:                                            #
#   Dat: Dataset containing species count data (in columns) in...       #
#        samples (in rows).                                             #
#        *matrix* or *vector*                                           #
#   Method: Diversity index to calculate. This can be "Berger" (for...  #
#           the Berger--Parker index), "Simpson" (for the Simpson...    #
#           index), "Shannon" (for the Shannon--Wiener index) or...     #
#           "Hill" (for Hill's ratio).                                  #
#           *character*, either of "Berger", "Simpson", "Shannon", or...#
#              "Hill"                                                   #
#           default="Shannon"                                           #
#   Correction: Bias correction for the Shannon-Wiener index. Either... #
#               of NULL (for no correction), "Chao1" for simple bias... #
#               correction (Chao et al. 2000), or "Chao2" for...        #
#               bias correction (Chao and Shen 2003).                   #
#               *character* or NULL                                     #
#               default=NULL                                            #
#   CI: Should confidence intervals be calculated? Either NULL (no...   #
#       confidence intervals) or one of "Variance" or "Boot". With...   #
#       "Variance", the calculated confidence intervals are based on... #
#       the estimated variance of the sample (this only works for the...#
#       Simpson and Shannon-Wiener indices). With "Boot", confidence... #
#       intervals are robustly calculated by bootstrapping (works for...#
#       all indices).                                                   #
#       NULL or *vector*, either of "Variance" or "Boot"                #
#       default=NULL                                                    #
#   Boot.N: Number of bootstrap replication in case of CI="Boot".       #
#           *numeric (integer)*                                         #
#           default=999                                                 #
#   Plot: Should results be plotted?                                    #
#         *logical*                                                     #
#         default=TRUE                                                  #
# Output data: List of diversity indices, with confidence intervals...  #
#              desired.                                                 #
# Input dataset: Count data of species in columns vs. samples in rows.  #
#########################################################################

Diversity<-function (Dat, Method="Shannon", Correction=NULL, CI=NULL, Boot.N=999, Plot=TRUE) {
	#Test data consistency
	{if (class(Dat)=="numeric") {Length<-1}
	else {Length<-nrow(Dat)}}
	if (max(Dat)==1) {stop("Data seem to be presence/absence data. Diversity indices cannot be calculated, use species richness indices instead!")}
	if (Method!="Berger" & Method!="Simpson" & Method!="Shannon" & Method!="Hill") {stop("Method must be either of 'Berger'. 'Simpson', 'Shannon', or 'Hill'!")}
	if (!is.null(Correction)) {if (Correction!="Chao1" & Correction!="Chao2") {stop("Correction must be either of NULL, 'Chao1', or 'Chao2'!")}}
	if (!is.null(Correction) & Method!="Shannon") {warning("Bias correction can only be applied to Shannon-Wiener index.")}
	if (!is.null(CI)) {
		if (CI!="Variance" & CI!="Boot") {stop("CI must be either NULL or one of 'Variance' or 'Boot'!")}
		if ((CI=="Variance" & Method=="Berger") | (CI=="Variance" & Method=="Hill")) {stop("Confidence interval variance via variation is only possible for Simpson and Shannon-Wiener indices!")}
		if (CI=="Variance" & (!is.null(Correction) && Correction=="Chao2")) {stop("Confidence interval estimation for Chao2 correction only possible via bootstrapping!")}
	}
	Boot.N<-round(Boot.N, digits=0)
	if (Plot==TRUE & Length==1) {warning("Dat only contains one sample, no plot will be produced.")}
	if (is.data.frame(Dat)==TRUE) {Dat<-as.matrix(Dat)}
	
	#Calculate sample size
	N<-apply(Dat, 1, sum)
	
	#Set up results object
	Res<-Sp.N<-vector(mode="numeric", length=Length)
	if (Method=="Shannon") {Res2<-vector(mode="numeric", length=Length)}
	if (Method=="Simpson" | Method=="Shannon") {FuncDiv<-vector(mode="numeric", length=Length)}
	if (!is.null(CI)) {CI.Low<-CI.Up<-vector(mode="numeric", length=Length)}
	if (!is.null(CI) & Method=="Shannon") {CI.Low2<-CI.Up2<-vector(mode="numeric", length=Length)}
	MinMax<-matrix(NA, Length, 2)
	colnames(MinMax)<-c("Min", "Max")
	rownames(MinMax)<-rownames(Dat)
	
	#Calculate diversity index
	for (i in 1:Length) {
		{if (class(Dat)=="numeric") {Temp<-Dat}
		else {Temp<-Dat[i,]}}
		
		#Berger--Parker index
		{if (Method=="Berger") {
			Res[i]<-max(Temp)/N[i]
			Sp.N[i]<-length(which(Temp>0))
			MinMax[i,1]<-1/length(which(Temp>0))
			MinMax[i,2]<-1
			
			#Calculate confidence interval via bootstrapping
			if (!is.null(CI)) {
				Boot.vals<-vector(mode="numeric", length=Boot.N)
				Size<-sum(Temp)
				Community<-vector(mode="numeric", length=0)
				for (j in 1:length(Temp)) {
					Community<-append(Community, rep(j, Temp[j]))
				}
				for (k in 1:Boot.N) {
					RandCom<-sample(Community, size=Size, replace=TRUE)
					RandAssem<-vector(mode="numeric", length=length(Temp))
					for (j in 1:length(RandAssem)) {
						RandAssem[j]<-length(which(RandCom==j))
					}
					Boot.vals[k]<-max(RandAssem)/Size
				}
				CI.Low[i]<-quantile(Boot.vals, probs=0.025, na.rm=TRUE)
				CI.Up[i]<-quantile(Boot.vals, probs=0.975, na.rm=TRUE)
			}
		}
		#Simpson index
		else if (Method=="Simpson") {
			p<-vector(mode="numeric", length=length(Temp))
			for (j in 1:length(Temp)) {
				p[j]<-Temp[j]/N[i]
			}
			Res[i]<-1-(sum(p**2))
			FuncDiv[i]<-1/(sum(p**2))
			Sp.N[i]<-length(which(Temp>0))
			MinMax[i,1]<-1-(sum(c(1, rep(0, length(p)-1))**2))
			MinMax[i,2]<-1-(sum((rep(1/length(p), length(p)))**2))
			
			#Calculate confidence interval
			if (!is.null(CI)) {
				#Via variance
				{if (CI=="Variance") {
					Var<-(4*N[i]*(N[i]-1)*(N[i]-2)*sum(p**3)+2*N[i]*(N[i]-1)*sum(p**2)-2*N[i]*(N[i]-1)*(2*N[i]-3)*(sum(p**2))**2)/((N[i]*(N[i]-1))**2)
					CI.Low[i]<-Res[i]-(1.96*sqrt(Var))
					CI.Up[i]<-Res[i]+(1.96*sqrt(Var))
				}
				#Via bootstrapping
				else {
					Boot.vals<-vector(mode="numeric", length=Boot.N)
					Size<-sum(Temp)
					Community<-vector(mode="numeric", length=0)
					for (j in 1:length(Temp)) {
						Community<-append(Community, rep(j, Temp[j]))
					}
					for (k in 1:Boot.N) {
						RandCom<-sample(Community, size=Size, replace=TRUE)
						RandAssem<-vector(mode="numeric", length=length(Temp))
						for (j in 1:length(RandAssem)) {
							RandAssem[j]<-length(which(RandCom==j))
						}
						p<-vector(mode="numeric", length=length(RandAssem))
						for (j in 1:length(RandAssem)) {
							p[j]<-RandAssem[j]/Size
						}
						Boot.vals[k]<-1-(sum(p**2))
					}
					CI.Low[i]<-quantile(Boot.vals, probs=0.025, na.rm=TRUE)
					CI.Up[i]<-quantile(Boot.vals, probs=0.975, na.rm=TRUE)
				}
				}
			}
		}
		#Shannon-Wiener index
		else if (Method=="Shannon") {
			p<-vector(mode="numeric", length=length(Temp))
			for (j in 1:length(Temp)) {
				p[j]<-Temp[j]/N[i]
			}
			p<-p[which(p>0)]
			{if (is.null(Correction)) {Res[i]<-(sum(p*log(p, base=exp(1))))*-1}
			else if (Correction=="Chao1") {Res[i]<-((sum(p*log(p, base=exp(1))))*-1)+((Chao(Temp)[2]-1)/(N[i]*2))}
			else {
				C<-1-(length(which(Temp==1))/N[i])
				Res[i]<-(sum((C*p*log(C*p, base=exp(1)))/(1-(1-C*p)^N[i])))*-1
			}}
			FuncDiv[i]<-prod(p**(-p))
			Sp.N[i]<-length(which(Temp>0))
			MinMax[i,1]<-0
			MinMax[i,2]<-log(length(which(Temp>0)), base=exp(1))
			
			#Calculate Pielou's equitability
			Res2[i]<-Res[i]/MinMax[i,2]
			
			#Calculate confidence interval
			if (!is.null(CI)) {
				#Via variance
				{if (CI=="Variance" & (!is.null(Correction) && Correction!="Chao2")) {
					Var<-(((sum(p*log(p**2, base=exp(1)))*-1)-((sum(p*log(p, base=exp(1))))**2)/N[i])/(N[i]**2))+((length(which(Temp>0))-1)/(2*N[i]**2))
					CI.Low[i]<-Res[i]-(1.96*sqrt(Var))
					CI.Up[i]<-Res[i]+(1.96*sqrt(Var))
					CI.Low2[i]<-CI.Low[i]/MinMax[i,2]
					CI.Up2[i]<-CI.Up[i]/MinMax[i,2]
				}
				#Via bootstrapping
				else {
					Boot.vals<-vector(mode="numeric", length=Boot.N)
					Size<-sum(Temp)
					Community<-vector(mode="numeric", length=0)
					for (j in 1:length(Temp)) {
						Community<-append(Community, rep(j, Temp[j]))
					}
					for (k in 1:Boot.N) {
						RandCom<-sample(Community, size=Size, replace=TRUE)
						RandAssem<-vector(mode="numeric", length=length(Temp))
						for (j in 1:length(RandAssem)) {
							RandAssem[j]<-length(which(RandCom==j))
						}
						p<-vector(mode="numeric", length=length(RandAssem))
						for (j in 1:length(RandAssem)) {
							p[j]<-RandAssem[j]/Size
						}
						p<-p[which(p>0)]
						{if (is.null(Correction)) {Boot.vals[k]<-(sum(p*log(p, base=exp(1))))*-1}
						else if (Correction=="Chao1") {Boot.vals[k]<-((sum(p*log(p, base=exp(1))))*-1)+((Chao(RandAssem)[2]-1)/(Size*2))}
						else {
							C<-1-(length(which(RandAssem==1))/Size)
							Boot.vals[k]<-(sum((C*p*log(C*p, base=exp(1)))/(1-(1-C*p)^Size)))*-1
						}}
					}
					CI.Low[i]<-quantile(Boot.vals, probs=0.025, na.rm=TRUE)
					CI.Up[i]<-quantile(Boot.vals, probs=0.975, na.rm=TRUE)
					CI.Low2[i]<-CI.Low[i]/MinMax[i,2]
					CI.Up2[i]<-CI.Up[i]/MinMax[i,2]
				}
				}
			}
		}
		#Hill's ratio
		else {
			p<-vector(mode="numeric", length=length(Temp))
			for (j in 1:length(Temp)) {
				p[j]<-Temp[j]/N[i]
			}
			p<-p[which(p>0)]
			D<-sum(p**2)
			H<-(sum(p*log(p, base=exp(1))))*-1
			Res[i]<-((1/D)-1)/(exp(H)-1)
			Sp.N[i]<-length(which(Temp>0))
			
			#Calculate confidence interval via bootstrapping
			if (!is.null(CI)) {
				Boot.vals<-vector(mode="numeric", length=Boot.N)
				Size<-sum(Temp)
				Community<-vector(mode="numeric", length=0)
				for (j in 1:length(Temp)) {
					Community<-append(Community, rep(j, Temp[j]))
				}
				for (k in 1:Boot.N) {
					RandCom<-sample(Community, size=Size, replace=TRUE)
					RandAssem<-vector(mode="numeric", length=length(Temp))
					for (j in 1:length(RandAssem)) {
						RandAssem[j]<-length(which(RandCom==j))
					}
					p<-vector(mode="numeric", length=length(RandAssem))
					for (j in 1:length(RandAssem)) {
						p[j]<-RandAssem[j]/Size
					}
					p<-p[which(p>0)]
					D<-sum(p**2)
					H<-(sum(p*log(p, base=exp(1))))*-1
					Boot.vals[k]<-((1/D)-1)/(exp(H)-1)
				}
				CI.Low[i]<-quantile(Boot.vals, probs=0.025, na.rm=TRUE)
				CI.Up[i]<-quantile(Boot.vals, probs=0.975, na.rm=TRUE)
			}
		}
		}
	}
	
	#Plot results
	if (Plot==TRUE & Length>1) {
		win.graph((Length/2.5), 10, 10)
		par(mar=c(4.1, 4.1, 4.1, 4.1))
		{if (!is.null(CI)) {
			YLIM<-c(min(c(Res, CI.Low), na.rm=TRUE), max(c(Res, CI.Up), na.rm=TRUE))
			plot(Res, type="b", pch=16, lwd=2, col="cornflowerblue", xlab="Sample", ylim=YLIM, ylab="Diversity")
			Breaks<-sort(unique(c(which(is.nan(CI.Low)), which(is.na(CI.Up)))))
			Poly<-Polx<-list()
			{if (length(Breaks)>0) {
				for (i in 1:(length(Breaks)+1)) {
					{if (i==1) {
						Poly[[i]]<-c(CI.Low[1:(Breaks[1]-1)], rev(CI.Up[1:(Breaks[1]-1)]))
						Polx[[i]]<-c(1:(Breaks[1]-1), rev(1:(Breaks[1]-1)))
					}
					else if (i==(length(Breaks)+1)) {
						Poly[[i]]<-c(CI.Low[(Breaks[i-1]+1):length(CI.Low)], rev(CI.Up[(Breaks[i-1]+1):length(CI.Up)]))
						Polx[[i]]<-c((Breaks[i-1]+1):length(CI.Low), rev((Breaks[i-1]+1):length(CI.Up)))
					}
					else {
						Poly[[i]]<-c(CI.Low[(Breaks[i-1]+1):(Breaks[i]-1)], rev(CI.Up[(Breaks[i-1]+1):(Breaks[i]-1)]))
						Polx[[i]]<-c((Breaks[i-1]+1):(Breaks[i]-1), rev((Breaks[i-1]+1):(Breaks[i]-1)))
					}
					}
				}
			}
			else {
				Poly[[1]]<-c(CI.Low, rev(CI.Up))
				Polx[[1]]<-c(1:Length, Length:1)
			}
			}
			for (i in 1:length(Poly)) {
				polygon(Polx[[i]], Poly[[i]], border=NA, col=rgb(100, 149, 237, alpha=80, maxColorValue=255))
			}
		}
		else {plot(Res, type="b", pch=16, lwd=2, col="cornflowerblue", xlab="Sample", ylab="Diversity")}}
		{if (Method=="Berger") {title("Berger-Parker Index")}
		else if (Method=="Simpson") {title("Simpson Index")}
		else if (Method=="Shannon") {title("Shannon-Wiener index")}
		else {title("Hill's Ratio")}}
		if (Method=="Shannon") {
			par(new=TRUE)
			{if (!is.null(CI)) {
				YLIM<-c(min(c(Res2, CI.Low2), na.rm=TRUE), max(c(Res2, CI.Up2), na.rm=TRUE))
				plot(Res2, type="b", pch=17, lwd=2, lty=2, col="chartreuse4", xlab="", ylim=YLIM, ylab="", axes=FALSE)
				Breaks<-sort(unique(c(which(is.nan(CI.Low2)), which(is.na(CI.Up2)))))
				Poly<-Polx<-list()
				{if (length(Breaks)>0) {
					for (i in 1:(length(Breaks)+1)) {
						{if (i==1) {
							Poly[[i]]<-c(CI.Low2[1:(Breaks[1]-1)], rev(CI.Up2[1:(Breaks[1]-1)]))
							Polx[[i]]<-c(1:(Breaks[1]-1), rev(1:(Breaks[1]-1)))
						}
						else if (i==(length(Breaks)+1)) {
							Poly[[i]]<-c(CI.Low2[(Breaks[i-1]+1):length(CI.Low2)], rev(CI.Up2[(Breaks[i-1]+1):length(CI.Up2)]))
							Polx[[i]]<-c((Breaks[i-1]+1):length(CI.Low2), rev((Breaks[i-1]+1):length(CI.Up2)))
						}
						else {
							Poly[[i]]<-c(CI.Low2[(Breaks[i-1]+1):(Breaks[i]-1)], rev(CI.Up2[(Breaks[i-1]+1):(Breaks[i]-1)]))
							Polx[[i]]<-c((Breaks[i-1]+1):(Breaks[i]-1), rev((Breaks[i-1]+1):(Breaks[i]-1)))
						}
						}
					}
				}
				else {
					Poly[[1]]<-c(CI.Low2, rev(CI.Up2))
					Polx[[1]]<-c(1:Length, Length:1)
				}
				}
				for (i in 1:length(Poly)) {
					polygon(Polx[[i]], Poly[[i]], border=NA, col=rgb(69, 139, 0, alpha=80, maxColorValue=255))
				}
			}
			else{plot(Res2, type="b", pch=17, lwd=2, lty=2, col="chartreuse4", xlab="", ylab="", axes=FALSE)}}
			axis(side=4, line=0)
			mtext("Equitability", side=4, line=3)
			legend("topright", legend=c("Diversity", "Equitability"), lty=c(1, 2), col=c("cornflowerblue", "chartreuse4"), lwd=2, pch=c(16, 17), bty="n", cex=0.8)
		}
	}
	
	#Return results
	names(Res)<-rownames(Dat)
	{if (Method=="Berger") {Name<-"Berger-Parker index"}
	else if (Method=="Simpson") {Name<-"Simpson index"}
	else if (Method=="Shannon") {
		{if (is.null(Correction)) {Name<-"Shannon-Wiener index"}
		else if (Correction=="Chao1") {Name<-"Shannon-Wiener index (corrected after Chao et al. 2000)"}
		else {Name<-"Shannon-Wiener index (corrected after Chao and Shen 2003)"}}
	}
	else {Name<-"Hill's ratio"}}
	
	{if (is.null(CI)) {
		{if (Method=="Hill") {return(list(Method=Name, Species.Richness=Sp.N, Index=Res))}
		else if (Method=="Shannon") {return(list(Method=Name, Species.Richness=Sp.N, Index=Res, Limits=MinMax, Equitability=Res2, Leti.index=FuncDiv))}
		else if (Method=="Simpson") {return(list(Method=Name, Species.Richness=Sp.N, Index=Res, Limits=MinMax, Laakso.Taagepera.index=FuncDiv))}
		else {return(list(Method=Name, Species.Richness=Sp.N, Index=Res, Limits=MinMax))}}
	}
	else if (CI=="Variance") {
		{if (Method=="Hill") {return(list(Method=Name, Species.Richness=Sp.N, Index=Res))}
		else if (Method=="Simpson") {return(list(Method=Name, Species.Richness=Sp.N, Index=Res, CI.Lower=CI.Low, CI.Upper=CI.Up, Limits=MinMax, Laakso.Taagepera.index=FuncDiv))}
		else if (Method=="Shannon") {return(list(Method=Name, Species.Richness=Sp.N, Index=Res, CI.Lower=CI.Low, CI.Upper=CI.Up, Limits=MinMax, Equitability=Res2, Equi.CI.Lower=CI.Low2, Equi.CI.Upper=CI.Up, Leti.index=FuncDiv))}
		else {return(list(Method=Name, Species.Richness=Sp.N, Index=Res, Limits=MinMax))}}
	}
	else {
		{if (Method=="Hill") {return(list(Method=Name, Species.Richness=Sp.N, Index=Res, CI.Lower=CI.Low, CI.Upper=CI.Up))}
		else if (Method=="Shannon") {return(list(Method=Name, Species.Richness=Sp.N, Index=Res, CI.Lower=CI.Low, CI.Upper=CI.Up, Limits=MinMax, Equitability=Res2, Equi.CI.Lower=CI.Low2, Equi.CI.Upper=CI.Up, Leti.index=FuncDiv))}
		else if (Method=="Ssimpson") {return(list(Method=Name, Species.Richness=Sp.N, Index=Res, CI.Lower=CI.Low, CI.Upper=CI.Up, Limits=MinMax, Laakso.Taagepera.index=FuncDiv))}
		else {return(list(Method=Name, Species.Richness=Sp.N, Index=Res, CI.Lower=CI.Low, CI.Upper=CI.Up, Limits=MinMax))}}
	}
	}
}

#########################################################################
# Function to calculate evenness coefficients according to Hurlbert...  #
#   (1971), Patten (1962), or Lloyd and Ghelardi (1964). Compare...     #
#   Legendre and Legendre (2012), p. 256.                               #
# Necessary packages: vegan                                             #
# Necessary input variables:                                            #
#   Div: A diversity object as produced by function Diversity, using... #
#        any supported diversity index except Hill's ratio.             #
#        *list*                                                         #
# Output data: List of evenness in samples.                             #
# Input dataset: List of diversity data as produced by function...      #
#                Diversity.                                             #
#########################################################################

#Load packages
require(vegan)

Evenness<-function (Div) {
	#Test data consistency
	if (Div$Method=="Hill's ratio") {stop("Evenness cannot be calculated for Hill's ratio!")}
	
	#Set up results objects
	Hurlbert<-Patten<-Lloyd<-vector(mode="numeric", length=length(Div$Index))
	
	#Calculate evenness
	for (i in 1:length(Div$Index)) {
		Hurlbert[i]<-(Div$Index[i]-Div$Limits[i,1])/(Div$Limits[i,2]-Div$Limits[i,1])
		Patten[i]<-(Div$Limits[i,2]-Div$Index[i])/(Div$Limits[i,2]-Div$Limits[i,1])
		if (Div$Method=="Shannon-Wiener index") {
			Broken<-bstick(Div$Species.Richness[i])
			M<-(sum(Broken*log(Broken, base=exp(1))))*-1
			Lloyd[i]<-Div$Index[i]/M
		}
	}
	
	#Return results
	{if (Div$Method=="Shannon-Wiener index") {return(list(Hurlbert=Hurlbert, Patten=Patten, Lloyd.Ghelardi=Lloyd))}
	else {return(list(Hurlbert=Hurlbert, Patten=Patten))}}
}

#########################################################################
# Function to apply two-proportions z-test to compare abundances in...  #
#    samples.                                                           #
# Necessary input variables:                                            #
#    Input: Input matrix.                                               #
#           *matrix*                                                    #
#    Obs.col: Column number containing the observed incidences or...    #
#             proportions.                                              #
#             *numeric (integer)*                                       #
#    N.col: Column number containing the sample sizes.                  #
#           *numeric (integer)*                                         #
#    Incidence: Does Obs.col contain incidence counts or proportions... #
#               (0-1)?                                                  #
#               TRUE: Counts are given.                                 #
#               FALSE: Proportions are given.                           #
#               default=TRUE                                            #
#    adj: Adjustment algorithm for p-value (compare p.adjust()).        #
#         *character*                                                   #
#         default="BY"                                                  #
# Output data: P-values whether or not groups are different...          #
#              concerning the observed incidences.                      #
# Input dataset: Matrix object, samples in rows, variables in columns.  #
#########################################################################

Prop.z.test<-function(Input, Obs.col, N.col, Incidence=TRUE, adj="BY") {
	#Read data
	Dat<-Input
	
	#Check data for consistency
	check.integer<-function(x){
		!length(grep("[^[:digit:]]", as.character(x)))
	}
	if (check.integer(Obs.col)==FALSE | check.integer(N.col)==FALSE) {stop("Only integers are allowed as column numbers!")}
	if (Incidence==TRUE) {if (check.integer(Dat[,Obs.col])==FALSE) {stop("Column with observations must contain integers only!")}}
	if (check.integer(Dat[,N.col])==FALSE) {stop("Column with sample size must contain integers only!")}
	if (nrow(Dat)<2) {stop("At least two samples needed for comparison!")}
	if (Incidence==FALSE & max(Dat[,Obs.col])>1) {stop("Proportions must be scaled as fraction (0-1), not per cent!")}
	if (Incidence==TRUE & (any(Dat[,Obs.col]>Dat[,N.col]))) {stop("Incidences cannot be larger than sample sizes!")}
	if (!adj %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")) {stop("adj must be supported by p.adjust()!")}
	
	#Setup results object
	{if (nrow(Dat)==2) {Res<-vector(mode="numeric", length=1); names(Res)<-"p-value"}
	else {
		Res<-matrix(NA, sum(1:(nrow(Dat)-1)), 4)
		colnames(Res)<-c("Sample 1", "Sample 2", "p-value", "adj. p-value")
	}}
	
	#Calculate results
	##Calculate proportions if necessary
	if (Incidence==TRUE) {Dat[,Obs.col]<-Dat[,Obs.col]/Dat[,N.col]}
	
	##Calculate p-values
	C<-1
	for (i in 1:(nrow(Dat)-1)) {
		for (j in (i+1):(nrow(Dat))) {
			PooledProp<-(Dat[i,Obs.col]*Dat[i,N.col]+Dat[j,Obs.col]*Dat[j,N.col])/(Dat[i,N.col]+Dat[j,N.col])
			SE<-sqrt(PooledProp*(1-PooledProp)*((1/Dat[i,N.col])+(1/Dat[j,N.col])))
			z<-(Dat[i,Obs.col]-Dat[j,Obs.col])/SE
			P<-pnorm(abs(z), mean=0, sd=1, lower.tail=FALSE)+pnorm(abs(z)*-1, mean=0, sd=1, lower.tail=TRUE)
			{if (nrow(Dat)>2) {Res[C,1]<-i; Res[C,2]<-j; Res[C,3]<-P; C<-C+1}
			else {Res[1]<-P}}
		}
	}
	if (nrow(Dat)>2) {Res[,4]<-p.adjust(Res[,3], method=adj)}
	print(Res)
}

#########################################################################
# Function to calculate confidence intervals for relative abundances... #
#	of assemblage data based on multinomial models.                 #
# Necessary input variables:                                            #
#    Data: Dataset containing absolute abundances of types...           #
#          (if Type=="Abs") or relative abundances of types (as...      #
#          fraction) and total number of specimens in sample (if...     #
#          Type=="Rel").                                                #
#          *matrix*                                                     #
#    Type: Does dataset contain absolute or relative abundances of...   #
#          types?                                                       #
#          "Abs": Absolute abundances of species provided.              #
#          "Rel": Relative abundances of species and total number of... #
#                 specimens per sample provided.                        #
#          *character*                                                  #
#    Alpha: Desired level of confidence.                                #
#           *real*                                                      #
#           default=0.05                                                #
#    SpecNum: Total number of species in dataset. Can be used to...     #
#             provide information about species richness in sample,...  #
#             if the dataset contains only a subset of the total...     #
#             count data. If NULL, it is assumed that the dataset...    #
#             contains all species and the value is calculated...       #
#             automatically.                                            #
#             *numeric (integer)*                                       #
#             default=NULL                                              #
# Output data: Array containing relative abundances and upper and...    #
#              lower bound of confidence interval of relative...        #
#              abundances.                                              #
# Input dataset: Matrix, samples in rows, types in n columns. If...     #
#	         Type=="Abs" each column must contain absolute...       #
#                abundance counts of one type. If Type=="Rel" the...    #
#                first 1--(n-1) columns contain the relative...         #
#                abundances of the types and the last (n'th) column...  #
#                contains the total absolute count (i.e. all...         #
#                specimens together) of the sample.                     #
#########################################################################

RelConf<-function (Data, Type, Alpha=0.05, SpecNum=NULL) {
	#Check data consistency
	if (Type!="Abs" && Type!="Rel") {stop("Type must be either 'Abs' or 'Rel'!")}
	if (Type=="Rel" && any(Data[,1:(ncol(Data)-1)]>1, na.rm=TRUE)) {stop("Relative abundances must not be >1!")}
	
	#Calculate relative abundances
	{if (Type=="Abs") {
		Sam.Size<-apply(Data, 1, sum)
		Data.Rel<-as.matrix(Data/Sam.Size)
	}
	else {
		Sam.Size<-Data[,ncol(Data)]
		Data.Rel<-as.matrix(Data[,1:(ncol(Data)-1)])
	}}
	
	#Calculate confidence intervals
	##Calculate W per species
	{if (is.null(SpecNum)) {k<-ncol(Data.Rel)}
	else {k<-SpecNum}}
	W<-vector(mode="numeric", length=nrow(Data.Rel))
	for (i in 1:(nrow(Data.Rel))) {
		N<-Sam.Size[i]
		{if (N==0) {W[i]<-0}
		else {W[i]<-((1/N)*qchisq(1-Alpha, k-1))+1}}
	}
	##Calculate confidence intervals
	CI<-array(NA, dim=c(nrow(Data.Rel), ncol(Data.Rel), 3), dimnames<-list(rownames(Data.Rel), colnames(Data.Rel), c("Rel.Abund", "Lower.Bound", "Upper.Bound")))
	for (i in 1:(nrow(Data.Rel))) {
		for (j in 1:(ncol(Data.Rel))) {
			CI[i,j,2]<-(W[i]+(2*Data.Rel[i,j])-1-sqrt(((1-W[i]-(2*Data.Rel[i,j]))^2)-(4*W[i]*Data.Rel[i,j]^2)))/(2*W[i])
			CI[i,j,3]<-(W[i]+(2*Data.Rel[i,j])-1+sqrt(((1-W[i]-(2*Data.Rel[i,j]))^2)-(4*W[i]*Data.Rel[i,j]^2)))/(2*W[i])
			
		}
	}
	CI[,,1]<-Data.Rel
	
	#Eliminate infinite values
	CI[,,2][which(!is.finite(CI[,,2]))]<-NA
	CI[,,3][which(!is.finite(CI[,,3]))]<-NA
	
	#Return results
	return(CI)
}

#--------------------------------------------

#Examples
#Species richness indices
#Richness(Data.Bin)
#Richness(Data.Bin, Method="Margalef", N=sample(1:130, size=20))
#Richness(Data.Bin[,1])
#Richness(Data.Count, Method="Menhinick")

#Species richness estimation
#Jack(Data.Jack)
#Jack2(Data.Jack, Replicates=100)
#Chao(Data.Chao)

#Diversity indices
#Diversity(Data.Count, Method="Berger")
#Diversity(Data.Count, Method="Simpson")
#Diversity(Data.Count, Method="Shannon")
#Diversity(Data.Count, Method="Hill")

#Diversity(Data.Count, Method="Simpson", CI="Variance")
#Div1<-Diversity(Data.Count, Method="Shannon", CI="Variance")
#Div2<-Diversity(Data.Count, Correction="Chao1", Method="Shannon", CI="Variance")
#Div3<-Diversity(Data.Count, Correction="Chao2", Method="Shannon", CI="Variance")#Error
#plot(Div1$Index, Div2$Index)

#Diversity(Data.Count, Method="Berger", CI="Boot", Boot.N=100)
#Diversity(Data.Count, Method="Simpson", CI="Boot", Boot.N=100)
#Div1<-Diversity(Data.Count, Method="Shannon", CI="Boot", Boot.N=100)
#Div2<-Diversity(Data.Count, Correction="Chao1", Method="Shannon", CI="Boot", Boot.N=100)
#Div3<-Diversity(Data.Count, Correction="Chao2", Method="Shannon", CI="Boot", Boot.N=100)
#plot(Div1$Index, Div2$Index)
#plot(Div1$Index, Div3$Index)
#YLIM<-c(min(c(Div1$Index, Div2$Index, Div3$Index)), max(c(Div1$Index, Div2$Index, Div3$Index)))
#plot(1:length(Div1$Index), Div1$Index, type="l", lwd=2, col="blue", main="Shannon-Wiener Index", xlab="Sample", ylim=YLIM, ylab="Diversity")
#lines(1:length(Div2$Index), Div2$Index, lwd=2, col="darkgreen")
#lines(1:length(Div3$Index), Div3$Index, lwd=2, col="firebrick")
#legend("bottomleft", lwd=2, col=c("blue", "darkgreen", "firebrick"), legend=c("Uncorrected", "Chao et al. (2000)", "Chao and Shen (2008)"))
#Diversity(Data.Count, Correction="Chao1", Method="Shannon", CI="Boot", Boot.N=100)
#Diversity(Data.Count, Method="Hill", CI="Boot", Boot.N=100)

#Evenness indices
#Div<-Diversity(Data.Count, Method="Simpson", Plot=FALSE)
#Evenness(Div)
#Div<-Diversity(Data.Count, Method="Shannon", Plot=FALSE)
#Evenness(Div)

#z-test
#Prop.z.test(Prop, Obs.col=2, N.col=1)
#Prop.z.test(Prop2, 2, 1, Incidence=FALSE)
#Prop.z.test(Prop[1:2,], Obs.col=2, N.col=1)

#Multinomial confidence intervals
#CI<-RelConf(Abund, Type="Abs")
#Sam.Size<-apply(Abund, 1, sum)
#RelAbund<-Abund/Sam.Size
#RelAbund<-cbind(RelAbund, Sam.Size)
#CI2<-RelConf(RelAbund, Type="Rel")

#--------------------------------------------
#--------------------------------------------
#Version History
#1.0	Finished Program
#1.1	Added jackknifing functions
#1.2	Added Chao estimator for estimating species richness from single samples
#	Added correction from Chao et al. (2000) to Shannon-Wiener index in Diversity
#	Added correction from Chao and Shen (2003) to Shannon-Wiener index in Diversity
#1.3	Included function RelConf from former MultinomialConfidence_Function.r (which is now deprecated)
#1.4	Included function Prop.z.test from former ProportionComparison_Function.r (which is now deprecated)
#1.4.1	Added choice of p-value adjusted to Prop.z.test
#--------------------------------------------
#--------------------------------------------
