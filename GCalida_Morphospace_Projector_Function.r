#These programs allow to project data-points obtained from specimens of planktonic
#Foraminifera of the Globigerinoides siphonifera/G. calida/G. radians plexus into the principal
#component analysis (PCA) morphospace calculated on the basis of scanning electron microscope
#(SEM) and light microscope (LM) images. Please refer to the final publication (see below)
#for an in depth description of the parameters used in the functions.

#Contents:	Step-by-step guide: Lines 16-48
#		Function SEM.Project: Lines 52-236
#		Function LM.Project: Lines 240-400

#Authors: Agnes K.M. Weiner, Manuel F.G. Weinkauf, Atsushi Kurasawa, Kate F. Darling, Michal Kucera
#Programmers: Manuel F.G. Weinkauf (Manuel.Weinkauf@unige.ch), Michal Kucera (mkucera@marum.de)
#Version: 1.0
#Date: 31 March 2014

######################
# Step-by-step guide #
######################
#This is a detailed guide on how to use this program for users not familiar with R.
#You have to have R installed on your Computer. The term <PATH> in the following refers
#to the full file path, where the respective files are saved on your computer. Note that you
#can either use a slash '/' or a double backslash '\\' as separator between parts of <PATH>,
#but not a single backslash. Furthermore, the whole string, i.e. <PATH>+file name must be
#set in quotation marks. For example "C:/R_Data" and "C:\\R_Data" both work, but "C:\R_Data"
#will lead to an error. In the following, commands to enter in the R console are set in single
#quotation marks (e.g. 'getwd()')---those single quotation marks only serve as a separator
#between the description and the actual command and MUST NOT be typed in the console.
#	1.	Start the R program
#	2.	Read this program file by typing 'source("<PATH>/File S2.r")' into
#		the R console
#	3.	Perform the analysis by using one of the commands defined below, providing all
#		necessary information.
#		3.1 If your data are stored in a .txt file of the name 'File.txt',
#			invoke e.g. 'SEM.Project("<PATH>/File.txt", Import=TRUE)' via the command window
#		3.2 If your files are already read into R into a variable called 'Data',
#			invoke e.g. 'SEM.Project(Data, Colour=FALSE)' via the command window
#		3.3 If the last column of your data contains any coding for the colouring of
#			the points plotted, invoke e.g. 'SEM.Project(Data)'
#		3.4 To provide custom names to be plotted alongside the data points use the 'Names'
#			vector (coded as 'c("Element 1", "Element 2", ... , "Element n")'), e.g. 
#			for a dataset 'Data' with four entries call
#			'SEM.project(Data, Names=c("Carib.", "Pacific", "Atlantic", "Mozambique"))'
#		3.5 The full R functionality is naturally supported. Having a dataset 'Data' with
#			station names in the first column, can be plotted as
#			'SEM.project(Data[,-1], Names=Data[,1])
#	4.	As output a PCA graph is produced, containing the convex hulls of our analysis and
#		projected points of your data. This graph can for instance be saved in various formats
#		via 'File->Save as...'

#####################################################################################################################

###############
# SEM.Project #
###############

#################################################################################
#Function to project data obtained from SEM into PCA morphospace obtained	#
#	by our analysis on the basis of SEM data				#
#										#
#Required input dataset								#
#	A matrix or tab-delimited .txt file containing specimens in rows and	#
#	measurements in 6 columns. Measurements has to occur in the following	#
#	order:									#
#		1-Angle alpha, i.e. angle between plane of spiralisation and	#
#			last whorl						#
#		2-Ratio E(l), i.e. height/width ratio of last chamber in	#
#			lateral view						#
#		3-Fraction PS, i.e. to what degree does the chamber directly in	#
#			front of the aperture cover the aperture. Note that	#
#			this value must be given as fraction (0<x<1) not in	#
#			percent (0<x<100).					#
#		4-Mean chamber elongation E, i.e. mean of height/width ratio of	#
#			all chambers in the last whorl in spiral/umbilical view	#
#		5-Elongation of last chamber E(L), i.e. height/width ratio of	#
#			last chamber in spiral/umbilical view			#
#		6-Mean angle gamma, i.e. mean angle between all chambers in last#
#			whorl in spiral/umbilical view				#
#										#
#Output: A PCA biplot, containing the convex hulls of the three morphotypes we	#
#		defined in the publication, and the specimens from the input	#
#		dataset + the position of type specimens projected into the	#
#		same morphospace.						#
#										#
#Parameters									#
#	Input: The name of the input matrix or file				#
#		*string* if Import=TRUE						#
#		*variable* if Import=FALSE					#
#	Names: Names vector, giving labels for points of Input dataset to be...	#
#		plotted. If not provided, consecutive numbers according to...	#
#		the row number of the respective entry will be plotted.		#
#		*vector*							#
#		default=NULL							#
#	Legend: Should a legend for the colour of test data points be plotted?	#
#		Only used if Colour=TRUE					#
#		*logical*							#
#		default=TRUE							#
#	Import: Should the input dataset be imported from a tab-delimited	#
#		.txt file, or read from an R-object				#
#		If .txt file this is supposed to contain column headings in the	#
#		first row.							#
#		*logical*							#
#		TRUE: Input expected as string, giving the name (or path) of	#
#			the file to be read (for inexperienced users)		#
#		FALSE: Input expected as variable name of matrix in R (for...	#
#			experienced users)					#
#		default=FALSE							#
#	Colour: Should the points of the input dataset be plotted in different	#
#		colours, in that case encoded in the seventh column in some	#
#		meaningful way							#
#		*logical*							#
#		default=TRUE							#
#################################################################################

SEM.Project<-function (Input, Names=NULL, Legend=TRUE, Import=FALSE, Colour=TRUE) {
	#Provide hard-coded original datasets
	##Hard code data from SEM images
	SEMDat<-matrix(c(4.682106452, 28.54989065, 4.050794594, 6.105889097, 5.419854033, 8.755343682, 14.36429531, 1.635293338, 7.258153505, 22.47317379, 1.943448379, 5.226490662, 25.13346166, 3.178328917, 9.094547821, 5.787059671, 13.02456309, 11.78138526, 10.90961968, 17.99352879, 9.358340834, 12.45341827, 2.512877765, 3.314249683, 12.2249818, 10.59464407, 16.25723813, 12.35770285, 7.166992395, 2.355410743, 1.231477474, 4.477558249, 14.83868693, 7.764730322, 8.549152436, 8.435541533, 2.997239721, 3.098937862, 23.70902272, 22.14874516, 0.707141695, 13.65020411, 23.95743137, 1.668521023, 2.278942702, 0.494981839, 16.24184475, 15.70995387, 12.56827235, 24.38065676, 9.506266635, 4.723532545, 32.34503484, 15.02525225, 8.206833049, 13.61545976, 13.01700486, 9.109121408, 5.821927677, 6.958913297, 8.306349571, 7.772012509, 16.6622175, 0.724571015, 0.919202242, 0.788767195, 0.77343943, 0.667827051, 0.712266216, 0.718693303, 0.751729469, 0.698498838, 0.844284764, 0.712759481, 0.715578582, 0.891396807, 0.740001188, 0.922224666, 0.706407341, 0.790383681, 0.781074317, 0.736229279, 0.875907937, 0.694060258, 0.930093799, 0.621044911, 0.702715496, 0.769921937, 0.709415251, 0.931197236, 0.999419766, 0.959098179, 0.795535171, 0.726488973, 0.633954487, 0.591908772, 0.737909377, 0.724087175, 0.756471993, 0.788297058, 0.799793381, 0.915056868, 0.901854197, 0.825518662, 0.907243696, 1.106783336, 0.741154507, 0.856453296, 0.737671111, 0.837356073, 0.718673422, 0.689807318, 0.775886326, 0.929197966, 0.732961339, 0.89521326, 0.738724196, 0.761038624, 0.704894394, 0.776949615, 0.79176613, 0.797899269, 0.807811751, 0.715294432, 0.748358162, 0.762354013, 1, 0.46101774, 1, 1, 1, 1, 0.984292048, 1, 1, 0.30171984, 1, 1, 0.247489824, 0.925561281, 0.968904188, 1, 0.74559322, 0.988659167, 1, 0.401767794, 0.919801642, 0.316417323, 1, 1, 0.946142954, 1, 0.769868198, 0.614503178, 1, 1, 1, 1, 1, 1, 1, 0.962585034, 1, 1, 0.407291842, 0.57181677, 1, 1, 0.695048309, 1, 1, 1, 1, 0.975862069, 0.91663276, 0.750754972, 0.780189928, 1, 0.383495146, 0.539035871, 0.922161172, 0.759656085, 0.617730721, 1, 0.942511841, 0.996376812, 1, 0.987068966, 0.626936532, 0.837194908, 0.885055747, 0.995709771, 0.898106416, 0.854016585, 0.943511388, 0.845100494, 0.839415591, 0.895510512, 0.791166547, 0.901531969, 0.933212039, 0.957913995, 0.867294019, 0.929727432, 0.916124214, 0.957905266, 0.816391415, 0.847375027, 1.008297938, 0.901334932, 0.912204912, 0.913259163, 0.869660726, 0.873259991, 0.895224615, 0.970206828, 0.975298698, 0.882244884, 0.956018595, 0.966618554, 0.845994941, 0.924542166, 0.816703349, 0.790336475, 0.922494229, 0.896207376, 0.842620072, 0.949692894, 1.048516716, 1.087283646, 0.957754858, 1.066808578, 1.110747495, 1.010128838, 0.887580861, 0.963318503, 0.869755891, 0.92251254, 1.000883007, 0.98052164, 0.878181341, 0.947086744, 0.897373597, 0.892070695, 0.929598197, 0.909052081, 0.944778446, 0.941966761, 0.946104761, 0.880812685, 0.94132115, 0.916158736, 0.882587998, 1.025080245, 0.905920883, 0.890488029, 0.801017787, 0.916662497, 0.715973602, 0.925386154, 0.776561764, 0.90273797, 0.894869918, 0.906656193, 0.986026744, 0.832343398, 1.08446816, 0.831911049, 0.9652623, 0.938505664, 0.827895216, 1.15974073, 0.759573527, 0.959277033, 0.911204407, 1.160567587, 0.880781933, 0.958915416, 1.062808843, 1.042960085, 1.012350414, 0.839806915, 0.896148111, 0.693195358, 0.819098755, 0.897594817, 0.891160088, 0.974523227, 0.981120821, 0.819278888, 1.080231987, 1.014650405, 0.970109369, 0.962144778, 1.175126363, 1.332750682, 1.149614367, 0.922633789, 0.891188444, 0.885978006, 0.831432849, 0.988685702, 1.04768025, 0.903157825, 0.931002667, 0.876708427, 0.886854721, 0.970317683, 0.838931997, 0.853070305, 0.924850144, 0.854618288, 0.757481042, 0.863936883, 1.012530117, 82.93040422, 78.13300689, 94.41664705, 74.27988646, 82.98349279, 78.43007123, 77.9193248, 83.50421176, 74.93759597, 91.97674128, 66.43327124, 74.54386286, 78.93182614, 75.88861825, 78.18593552, 66.86586241, 80.41721069, 78.00594149, 72.28674892, 74.31436172, 67.22469806, 68.55962969, 64.19700984, 73.30485016, 70.36434466, 69.30153832, 87.36537417, 76.03586432, 79.13864642, 73.6678195, 71.52233231, 73.13964416, 64.53559649, 80.50214085, 81.77957565, 85.73531082, 78.64280303, 85.58775154, 84.76142201, 84.77278814, 65.06360728, 67.12785481, 86.78658833, 65.58726775, 68.97517785, 70.76642556, 79.55772451, 77.50934969, 76.1505695, 77.14155721, 90.67022408, 81.65951826, 82.39089734, 81.93936668, 77.48352508, 76.54670983, 83.04737651, 74.53435608, 72.59958427, 82.99363041, 79.30630307, 77.81622145, 70.17891925),63,6)
	##Hard code data of drawings of type specimen
	TypeDat<-matrix(c(0.284921721, 29.95689354, 19.70988212, 6.435777062, 7.093455509, 1.335541369, 0.782767685, 0.889163663, 3.472637825, 0.792005173, 1, 0.742857143, 0, 0.961038961, 1, 1.139902895, 1.081459616, 0.932938249, 1.838292105, 0.908821281, 1.083362638, 0.749940858, 0.991323777, 3.160639703, 0.965629734, 71.44619451, 87.49269144, 85.39149987, 61.59074377, 85.00470704), 5, 6)
	##Read test data set and check for consistency
	{if (Import==TRUE) {TestDat<-read.table(Input, header=TRUE, sep="\t")}
	else {TestDat<-Input}}
	{if (Colour==FALSE) {if (dim(TestDat)[2]!=6) {stop("Input data must contain exactly six columns!")}}
	else {if (dim(TestDat)[2]!=7) {stop("Input data must contain exactly seven columns!")}}}
	
	#Prepare PCA
	##Calculate PCA of SEM data
	PCA<-princomp(SEMDat, cor=TRUE)
	##Predict position of holotypes
	Holotypes<-list()
	Holotypes$PCA<-predict(PCA, TypeDat)[,1:2]
	Holotypes$Names<-c("G. radians", "G. siphonifera", "G. calida", "G. adamsi", "G. siphonifera (l)")
	##Predict position of test dataset
	TestSpec<-list()
	{if (dim(TestDat)[1]==1) {TestSpec$PCA<-matrix(c(predict(PCA, TestDat[,1:3])[,1], predict(PCA, TestDat[,1:3])[,2]), 1, 2)}
	else {TestSpec$PCA<-predict(PCA, TestDat[,1:6])[,1:2]}}
	##Create point names
	{if (!is.null(Names)) {
		if (length(Names)!=dim(TestDat)[1]) {stop("Input data and Names vector are not of same length")}
		TestSpec$Names<-Names
	}
	else {TestSpec$Names<-1:(dim(TestDat)[1])}
	}
	##Create point colours
	if (Colour==TRUE) {
		{if (is.factor(TestDat[,7])==TRUE) {Lv<-levels(TestDat[,7])}
		else {Lv<-unique(TestDat[,7])}}
		Col<-rainbow(length(Lv))
		for (i in 1:(dim(TestSpec$PCA)[1])) {
			TestSpec$Colour<-append(TestSpec$Colour,  Col[which(Lv==TestDat[i,7])])
		}
	}
	
	#Prepare data for plotting
	##Main PCA data
	PCA.plot<-list()
	PCA.plot$Comp1<-PCA$scores[,1]
	PCA.plot$Comp2<-PCA$scores[,2]
	PCA.plot$Geno<-as.factor(c("Ib", "IIIb", "Ib", "IIa3", "Ib", "Ib", "Ib", "Ib", "Ib", "IIIb", "IIa5", "Ib", "IIIb", "Ib", "Ia", "IIa5", "IIa4", "Ib", "Ib", "IIIb", "IIa5", "IIIb", "IIa5", "IIa4", "IIa4", "IIa4", "Ia", "Ia", "Ia", "Ib", "IIa5", "IIa5", "IIa5", "Ib", "IIa1", "Ib", "Ib", "Ib", "IIIb", "IIIb", "Ia", "Ia", "IIIb", "Ia", "Ia", "IIa4", "IIa1", "IIa4", "Ib", "IIIb", "IIIb", "Ib", "IIIc", "IIa1", "IIa1", "IIa4", "IIa1", "IIa3", "IIa3", "IIa3", "Ib", "IIa3", "IIa3"))
	PCA.plot$Loadings$Names<-c("Alpha", "El", "PS", "E", "EL", "Gamma")
	PCA.plot$Loadings$Comp1<-PCA$loadings[1:length(PCA.plot$Loadings$Names)]
	PCA.plot$Loadings$Comp2<-PCA$loadings[(length(PCA.plot$Loadings$Names)+1):(length(PCA.plot$Loadings$Names)*2)]
	##Convex hulls
	G1<-list()
	G1$Comp1<-PCA.plot$Comp1[which(PCA.plot$Geno=="Ia")]
	G1$Comp2<-PCA.plot$Comp2[which(PCA.plot$Geno=="Ia")]
	HullI<-chull(G1$Comp1,G1$Comp2)
	HullI<-append(HullI,HullI[1])
	G2<-list()
	G2$Comp1<-PCA.plot$Comp1[which(PCA.plot$Geno!="Ia" & PCA.plot$Geno!="IIIb" & PCA.plot$Geno!="IIIc")]
	G2$Comp2<-PCA.plot$Comp2[which(PCA.plot$Geno!="Ia" & PCA.plot$Geno!="IIIb" & PCA.plot$Geno!="IIIc")]
	HullII<-chull(G2$Comp1,G2$Comp2)
	HullII<-append(HullII,HullII[1])
	G3<-list()
	G3$Comp1<-PCA.plot$Comp1[which(PCA.plot$Geno=="IIIb" | PCA.plot$Geno=="IIIc")]
	G3$Comp2<-PCA.plot$Comp2[which(PCA.plot$Geno=="IIIb" | PCA.plot$Geno=="IIIc")]
	HullIII<-chull(G3$Comp1,G3$Comp2)
	HullIII<-append(HullIII,HullIII[1])
	
	#Plot graph
	Scale=3
	Expand=Scale*1.1
	##Find axis ranges
	X1<-ifelse((min(TestSpec$PCA[,1])<(-5)) ,min(TestSpec$PCA[,1]), -5)
	X2<-ifelse((max(TestSpec$PCA[,1])>4.5) ,max(TestSpec$PCA[,1])+2.5, 4.5)
	Y1<-ifelse((min(TestSpec$PCA[,2])<(-5)) ,min(TestSpec$PCA[,1]), -5)
	Y2<-ifelse((max(TestSpec$PCA[,2])>4.5) ,max(TestSpec$PCA[,1])+1.5, 4.5)
	SX<-c(X1, X2)
	SY<-c(Y1, Y2)
	##Plot test data points
	{if (Colour==TRUE) {plot(TestSpec$PCA, xlim=SX, ylim=SY, xlab="PC 1", ylab="PC 2", main="Projection in Morphospace", pch=16, col=TestSpec$Colour)}
	else {plot(TestSpec$PCA, xlim=SX, ylim=SY, xlab="PC 1", ylab="PC 2", main="Projection in Morphospace", pch=16, col="black")}}
	text(TestSpec$PCA+0.05, adj=c(0, 0), labels=TestSpec$Names)
	#legend("bottomright", title="Type", bg=344, pch=c(16, 16, 15, 15, 15, 15, 17, 17), col=ColRange$Num, legend=ColRange$Levels)
	#Plot loading arrows
	for (i in 1:(length(PCA.plot$Loadings$Comp1))) {
		arrows(x0=0, y0=0, x1=PCA.plot$Loadings$Comp1[i]*Scale, y1=PCA.plot$Loadings$Comp2[i]*Scale, lwd=1, col="red", length=0.3, angle=20)
		text(PCA.plot$Loadings$Comp1[i]*Expand, PCA.plot$Loadings$Comp2[i]*Expand, labels=PCA.plot$Loadings$Names[i], col="red")
	}
	##Plot convex hulls
	lines(G1$Comp1[HullI], G1$Comp2[HullI], lwd=2, col="blue")
	lines(G2$Comp1[HullII], G2$Comp2[HullII], lwd=2, col="green")
	lines(G3$Comp1[HullIII], G3$Comp2[HullIII], lwd=2, col="brown")
	legend("topright", title="Groups", lwd=2, col=c("blue", "brown", "green"), legend=c("G. radians", "G. calida", "G. siphonifera"))
	if (Colour==TRUE && Legend==TRUE) {legend("bottomright", title="Data", pch=16, col=Col, legend=Lv)}
	##Plot holotypes
	points(Holotypes$PCA, cex=2, pch=8, lwd=2)
	text(Holotypes$PCA+0.1, labels=Holotypes$Names, adj=c(0, 0))
}

####
#Examples
##Create test data set
#Example<-matrix(NA, 10, 6)
#set.seed(10);Example[, 1]<-runif(10, 0.4, 35)
#set.seed(0.8);Example[, 2]<-runif(10, 0.5, 1.1)
#set.seed(0.8);Example[, 3]<-runif(10, 0.2, 1)
#set.seed(0.9);Example[, 4]<-runif(10, 0.8, 1.2)
#set.seed(0.9);Example[, 5]<-runif(10, 0.7, 1.3)
#set.seed(76);Example[, 6]<-runif(10, 60, 95)
##Run function
#SEM.Project(Example, Import=FALSE)

##Create test data set
#Example2<-matrix(NA, 10, 7)
#set.seed(10);Example2[, 1]<-runif(10, 0.4, 35)
#set.seed(0.8);Example2[, 2]<-runif(10, 0.5, 1.1)
#set.seed(0.8);Example2[, 3]<-runif(10, 0.2, 1)
#set.seed(0.9);Example2[, 4]<-runif(10, 0.8, 1.2)
#set.seed(0.9);Example2[, 5]<-runif(10, 0.7, 1.3)
#set.seed(76);Example2[, 6]<-runif(10, 60, 95)
#Example2[, 7]<-c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4)
#write.table(Example2, "Test.txt", sep="\t")
##Run function
#SEM.Project(Example2, Import=FALSE, Colour=TRUE)
#SEM.Project("TestSEM.txt", Import=TRUE, Colour=TRUE)

#####################################################################################################################

###############
# LM.Project  #
###############

#################################################################################
#Function to project data obtained from LM into PCA morphospace obtained	#
#	by our analysis on the basis of LM data					#
#										#
#Required input dataset								#
#	A matrix or tab-delimited .txt file containing specimens in rows and	#
#	measurements in 3 columns. Measurements has to occur in the following	#
#	order:									#
#		1-Elongation of last chamber E(L),  i.e. height/width ratio of	#
#			last chamber in spiral/umbilical view			#
#		2-Mean chamber elongation E,  i.e. mean of height/width ratio of#
#			the last three chambers in the last whorl in		#
#			spiral/umbilical view					#
#		3-Mean angle gamma,  i.e. mean angle between last three chambers#
#			in spiral/umbilical view				#
#										#
#Output: A PCA biplot,  containing the convex hulls of the three morphotypes we	#
#		defined in the publication,  and the specimens from the input	#
#		dataset projected into the same morphospace.			#
#										#
#Parameters									#
#	Input: The name of the input matrix or file				#
#		*string* if Import=TRUE						#
#		*variable* if Import=FALSE					#
#	Names: Names vector, giving labels for points of Input dataset to be...	#
#		plotted. If not provided, consecutive numbers according to...	#
#		the row number of the respective entry will be plotted.		#
#		*vector*							#
#		default=NULL							#
#	Legend: Should a legend for the colour of test data points be plotted?	#
#		Only used if Colour=TRUE					#
#		*logical*							#
#		default=TRUE							#
#	Import: Should the input dataset be imported from a tab-delimited	#
#		.txt file,  or read from an R-object				#
#		If .txt file this is supposed to contain column headings in the	#
#		first row.							#
#		*logical*							#
#		TRUE: Input expected as string, giving the name (or path) of	#
#			the file to be read (for inexperienced users)		#
#		FALSE: Input expected as variable name of matrix in R (for...	#
#			experienced users)					#
#		default=FALSE							#
#	Colour: Should the points of the input dataset be plotted in different	#
#		colours,  then encoded in the fourth column in some meaningful	#
#		way								#
#		*logical*							#
#		default=TRUE							#
#################################################################################

LM.Project<-function (Input, Names=NULL, Legend=TRUE, Import=FALSE, Colour=TRUE) {
	#Provide hard-coded original datasets
	##Hard code data from LM images
	LMDat<-matrix(c(1.012392504,1.024032531,1.048823228,0.942134492,0.946311235,0.923710984,0.801811754,0.938420849,0.931917862,0.760135121,1.004219346,0.936759377,0.966001546,1.028089202,1.023043052,0.895167832,0.873221973,1.033911032,0.952685866,0.954959694,0.823066742,0.886429508,0.940122184,0.905431503,1.037403801,0.901441419,0.949339355,0.838878316,1.028921486,0.938858876,0.975048829,0.970928784,0.84333518,0.904645669,0.898135464,0.822275108,0.911791785,0.949823403,0.890464247,0.968573112,0.858303086,1.087479249,0.95149985,1.07906606,0.786797631,0.752374263,0.920816759,0.945373259,0.995255913,0.977443984,0.981272874,0.946913334,0.79687863,0.926424194,0.703815564,0.892662731,0.855017242,0.84453635,0.939330272,0.781626836,0.908251259,0.918297886,0.970940416,0.974688199,1.001017007,0.969516186,0.967331106,0.932207991,0.914887347,0.745578492,0.979476044,0.894673445,0.901586173,0.905706865,0.797146651,0.880124476,1.000165056,0.957285673,0.888030728,0.934564106,0.935661729,0.877306908,0.939479501,1.116576926,1.038799063,0.789569945,0.943015616,1.02263805,1.031264282,1.005735999,0.938014373,0.94766375,0.852468914,0.843211072,0.938102632,0.84373979,1.082623613,1.012127485,0.852486857,0.888780543,0.928009829,0.7557799,0.983756418,1.078336241,0.978992833,1.080956854,0.98436193,1.02356666,0.914064863,1.006239593,0.937276937,0.992594139,0.998937201,0.951744447,1.092629571,0.984441672,1.037255365,1.077121343,1.073985906,1.012956633,1.030252223,0.953059083,0.981129765,0.969194967,0.935440438,1.096200023,1.075475628,0.919703714,0.944145297,0.928752164,1.01732633,0.972645237,0.948753405,0.967279411,0.840694973,0.940589954,0.86673322,0.87515478,0.871960515,0.948002259,0.833152026,1.015906709,0.933947614,0.991922956,0.787067151,0.969962663,0.852294713,0.878365377,0.930253933,0.837902545,0.854609566,0.914595274,0.871076185,0.893419058,0.87579619,0.868524711,0.904961154,0.864724994,0.95515627,0.875969157,0.826798578,0.889869054,0.800155287,0.813006686,0.909288752,0.953051927,0.765543667,0.926318496,0.894902235,0.894782812,0.923784537,0.92180149,0.93622032,0.943869664,0.886580252,0.999122277,0.838154642,0.805455678,1.050289823,1.003383522,0.796932666,0.952708054,0.914860108,0.89719073,0.873321638,0.937371669,0.961287268,0.921579516,0.84892481,0.892936861,0.938374147,0.884766687,0.98093259,0.843275221,0.952375703,0.809102232,0.883279986,0.75366134,0.936659128,0.915639297,0.913217232,1.020722974,0.928784538,0.932705732,1.03130219,0.876072695,0.904424336,0.942482381,0.819010945,0.950313014,0.899833352,0.951640577,0.85894521,0.83419646,0.86693218,0.878407757,0.950526873,0.869372523,0.924387855,0.879088783,0.909332974,0.91671508,0.854799953,0.85063277,0.953834524,0.894039455,0.897533179,0.943828678,0.849279496,0.841018063,0.948614106,0.927895307,0.896610919,0.985228891,0.999228247,1.048604848,0.994561386,0.874043019,0.961749511,1.02694212,0.946228063,0.940157374,1.058458696,0.965443694,0.927656138,0.892237845,0.904272601,0.901613627,0.951941267,0.940352805,0.787914409,0.85582439,0.817501165,0.965084806,0.949172068,0.900497018,90.48708028,80.59423432,86.51361536,91.99830081,81.54349573,85.72566902,81.37434775,89.85391966,100.227489,93.21736958,89.49049708,94.64902723,99.42835541,85.10455172,89.98045314,78.11207078,93.76105473,55.98692556,76.46618358,84.28029288,81.60832678,101.2411281,89.48933301,84.17813715,87.14888533,86.86398926,95.05838001,91.59961799,88.14114166,84.96845607,77.10646381,97.26818765,82.43499137,78.39943487,90.32646034,81.75573334,85.83946151,86.6242728,92.79003844,89.32108133,89.63122636,93.97293609,92.49456333,86.30067307,86.08773122,86.51393856,84.8717222,84.71428522,93.97285678,87.87224338,81.74030389,85.07366405,87.06952037,78.79652987,82.70023062,80.59963066,89.48429498,82.98025738,81.05727533,76.92276451,89.75260782,75.45727603,90.43821946,99.57062514,89.41975581,94.69318805,86.85198559,95.00839776,95.89983162,90.70359656,82.707874,95.41050787,91.89898049,77.87950812,80.90189327,92.72099518,82.05268658,89.32581515,78.14590053,86.07375375,84.18164338,81.32969478,86.65861283,86.46184683,93.33665149,89.52277412,80.60568901,70.57660511,77.97308711,81.3537753,82.19046465,99.8878702,90.2971525,72.44775463,67.7790357,78.72975157,71.61403971,86.9428264,94.70494823,97.60105442,71.76845495,105.5635386,83.50879001,102.5833792,87.36469494,81.3882819,84.79614434,88.0011472,85.61968481,90.30122854,84.27524611,90.49236874,87.37551224,92.2683507,87.97384152,84.61677195,87.1864483,90.23193757,100.0836429,97.27857344,79.15444644,89.95867608,79.60139485,93.54084776,95.79323978,70.41727226,95.68010009,86.18462572), 128, 3)
	##Read test data set and check for consistency
	{if (Import==TRUE) {TestDat<-read.table(Input, header=TRUE, sep="\t")}
	else {TestDat<-Input}}
	{if (Colour==FALSE) {if (dim(TestDat)[2]!=3) {stop("Input data must contain exactly three columns!")}}
	else {if (dim(TestDat)[2]!=4) {stop("Input data must contain exactly four columns!")}}}
	
	#Prepare PCA
	##Calculate PCA of SEM data
	PCA<-princomp(LMDat, cor=TRUE)
	##Predict position of test dataset
	TestSpec<-list()
	{if (dim(TestDat)[1]==1) {TestSpec$PCA<-matrix(c(predict(PCA, TestDat[,1:3])[,1], predict(PCA, TestDat[,1:3])[,2]), 1, 2)}
	else {TestSpec$PCA<-predict(PCA, TestDat[,1:3])[,1:2]}}
	##Create point names
	{if (!is.null(Names)) {
		if (length(Names)!=dim(TestDat)[1]) {stop("Input data and Names vector are not of same length")}
		TestSpec$Names<-Names
	}
	else {TestSpec$Names<-1:(dim(TestDat)[1])}
	}
	##Create point colours
	if (Colour==TRUE) {
		{if (is.factor(TestDat[,4])==TRUE) {Lv<-levels(TestDat[,4])}
		else {Lv<-unique(TestDat[,4])}}
		Col<-rainbow(length(Lv))
		for (i in 1:(dim(TestSpec$PCA)[1])) {
			TestSpec$Colour<-append(TestSpec$Colour,  Col[which(Lv==TestDat[i,4])])
		}
	}
	
	#Prepare data for plotting
	##Main PCA data
	PCA.plot<-list()
	PCA.plot$Comp1<-PCA$scores[, 1]
	PCA.plot$Comp2<-PCA$scores[, 2]
	PCA.plot$Geno<-as.factor(c("Ia", "Ia", "Ia", "Ia", "Ia", "Ib", "Ib", "Ib", "Ib", "Ib", "Ib", "Ib", "Ib", "Ib", "Ib", "Ib", "Ib", "Ib", "Ib", "Ib", "Ib", "Ib", "Ib", "Ib", "Ib", "Ib", "IIa1", "IIa1", "IIa1", "IIa1", "IIa1", "IIa1", "IIa1", "IIa2", "IIa2", "IIa2", "IIa2", "IIa2", "IIa2", "IIa2", "IIa2", "IIa3", "IIa3", "IIa3", "IIa3", "IIa3", "IIa3", "IIa3", "IIa3", "IIa3", "IIa3", "IIa3", "IIa3", "IIa3", "IIa3", "IIa3", "IIa3", "IIa3", "IIa3", "IIa3", "IIa3", "IIa3", "IIa3", "IIa3", "IIa3", "IIa3", "IIa3", "IIa3", "IIa3", "IIa3", "IIa3", "IIa3", "IIa3", "IIa3", "IIa3", "IIa3", "IIa3", "IIa4", "IIa4", "IIa4", "IIa4", "IIa4", "IIa4", "IIa4", "IIa4", "IIa4", "IIa4", "IIa4", "IIa4", "IIa4", "IIa5", "IIa5", "IIb", "IIb", "IIb", "IIb", "IIb", "IIb", "IIb", "IIb", "IIb", "IIb", "IIIb", "IIIb", "IIIb", "IIIb", "IIIb", "IIIb", "IIIb", "IIIb", "IIIb", "IIIb", "IIIb", "IIIb", "IIIb", "IIIb", "IIIb", "IIIb", "IIIc", "IIIc", "IIIc", "IIIc", "IIIc", "IIIc", "IIIc", "IIIc", "IIIc", "IIIc"))
	PCA.plot$Loadings$Names<-c("EL", "E", "Gamma")
	PCA.plot$Loadings$Comp1<-PCA$loadings[1:length(PCA.plot$Loadings$Names)]
	PCA.plot$Loadings$Comp2<-PCA$loadings[(length(PCA.plot$Loadings$Names)+1):(length(PCA.plot$Loadings$Names)*2)]
	##Convex hulls
	G1<-list()
	G1$Comp1<-PCA.plot$Comp1[which(PCA.plot$Geno=="Ia")]
	G1$Comp2<-PCA.plot$Comp2[which(PCA.plot$Geno=="Ia")]
	HullI<-chull(G1$Comp1,G1$Comp2)
	HullI<-append(HullI,HullI[1])
	G2<-list()
	G2$Comp1<-PCA.plot$Comp1[which(PCA.plot$Geno!="Ia" & PCA.plot$Geno!="IIIb" & PCA.plot$Geno!="IIIc")]
	G2$Comp2<-PCA.plot$Comp2[which(PCA.plot$Geno!="Ia" & PCA.plot$Geno!="IIIb" & PCA.plot$Geno!="IIIc")]
	HullII<-chull(G2$Comp1,G2$Comp2)
	HullII<-append(HullII,HullII[1])
	G3<-list()
	G3$Comp1<-PCA.plot$Comp1[which(PCA.plot$Geno=="IIIb" | PCA.plot$Geno=="IIIc")]
	G3$Comp2<-PCA.plot$Comp2[which(PCA.plot$Geno=="IIIb" | PCA.plot$Geno=="IIIc")]
	HullIII<-chull(G3$Comp1,G3$Comp2)
	HullIII<-append(HullIII,HullIII[1])
	
	#Plot graph
	Scale=3
	Expand=Scale*1.1
	##Find axis ranges
	X1<-ifelse((min(TestSpec$PCA[,1])<(-3.5)) ,min(TestSpec$PCA[,1]), -3.5)
	X2<-ifelse((max(TestSpec$PCA[,1])>5) ,max(TestSpec$PCA[,1])+2.5, 5)
	Y1<-ifelse((min(TestSpec$PCA[,2])<(-3)) ,min(TestSpec$PCA[,1]), -3)
	Y2<-ifelse((max(TestSpec$PCA[,2])>3) ,max(TestSpec$PCA[,1])+1.5, 3)
	SX<-c(X1, X2)
	SY<-c(Y1, Y2)
	##Plot test data points
	{if (Colour==TRUE) {plot(TestSpec$PCA, xlim=SX, ylim=SY, xlab="PC 1", ylab="PC 2", main="Projection in Morphospace", pch=16, col=TestSpec$Colour)}
	else {plot(TestSpec$PCA, xlim=SX, ylim=SY, xlab="PC 1", ylab="PC 2", main="Projection in Morphospace", pch=16, col="black")}}
	text(TestSpec$PCA+0.05, adj=c(0, 0), labels=TestSpec$Names)
	#Plot loading arrows
	for (i in 1:(length(PCA.plot$Loadings$Comp1))) {
		arrows(x0=0, y0=0, x1=PCA.plot$Loadings$Comp1[i]*Scale, y1=PCA.plot$Loadings$Comp2[i]*Scale, lwd=1, col="red", length=0.3, angle=20)
		text(PCA.plot$Loadings$Comp1[i]*Expand, PCA.plot$Loadings$Comp2[i]*Expand, labels=PCA.plot$Loadings$Names[i], col="red")
	}
	##Plot convex hulls
	lines(G1$Comp1[HullI], G1$Comp2[HullI], lwd=2, col="blue")
	lines(G2$Comp1[HullII], G2$Comp2[HullII], lwd=2, col="green")
	lines(G3$Comp1[HullIII], G3$Comp2[HullIII], lwd=2, col="brown")
	legend("topright", title="Groups", lwd=2, col=c("blue", "brown", "green"), legend=c("G. radians", "G. calida", "G. siphonifera"))
	if (Colour==TRUE && Legend==TRUE) {legend("bottomright", title="Data", pch=16, col=Col, legend=Lv)}
}

####
#Examples
##Create test data set
#Example<-matrix(NA, 10, 3)
#set.seed(0.9);Example[, 1]<-runif(10, 0.7, 1.3)
#set.seed(0.9);Example[, 2]<-runif(10, 0.8, 1.2)
#set.seed(76);Example[, 3]<-runif(10, 60, 95)
##Run function
#LM.Project(Example, Import=FALSE)

##Create test data set
#Example2<-matrix(NA, 10, 4)
#set.seed(0.9);Example2[, 1]<-runif(10, 0.7, 1.3)
#set.seed(0.9);Example2[, 2]<-runif(10, 0.8, 1.2)
#set.seed(76);Example2[, 3]<-runif(10, 60, 95)
#Example2[, 4]<-c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4)
#write.table(Example2, "TestLM.txt", sep="\t")
##Run function
#LM.Project(Example2, Import=FALSE, Colour=TRUE)
#LM.Project("TestLM.txt", Import=TRUE, Colour=TRUE)
