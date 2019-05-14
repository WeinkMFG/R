# R utility functions

This repository contains ASCII codes for use in R. You can principally open them in any text-editor, and either copy the code into the R console or use source("FilePath") to make the function known to R. Furthermore you can work directly from the console.

However, it is often more convenient to write code in an editor and send command lines directly to R. For that, two possible options are either [Notepad++](https://notepad-plus-plus.org/) with the [NppToR](https://sourceforge.net/projects/npptor/) extension installed on your computer, or the use of [Tinn-R](https://sourceforge.net/projects/tinn-r/) as an editor.

All codes in this repository are distributed under the [Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License](http://creativecommons.org/licenses/by-nc-sa/3.0/) and are also available from my [website](https://sites.google.com/site/weinkaufmanuel/).

# Description of functions

## Coordinate_TestData

A folder containing test data for the file CoordinateTranslation_Function.r.

## DiagramExtraction_Examples

A folder containing test diagrams for the file DiagramExtraction_Function.r.

## Morphometric_TestData

A folder containing test data for the diverse morphometric function files.

## BioDiversity_Function.r v. 1.4.1

A collection of functions for enhanced biodiversity studies, as described in Hammer and Harper (2006) and Legendre and Legendre (2012). It includes some functionalities that go beyond the abilities of packages like vegan, biodiversityR, and picante, but in turn does not contain all functionalities that are provided by those packages. Insofar, it is complementing, not replacing, any of those packages. Functions in this file allow to calculate species richness (including Margalef- and Menhinick corrections), calculate a variety of biodiversity indices (including confidence intervals) and equitability indices, and estimate species richness from replicated sampling via jackknifing. A function to calculate multinomial confidence intervals for for abundance data has also been included.

* Hammer, Ø. and Harper, D. (2006) *Paleontological Data Analysis*. 351 pp. (Blackwell: Malden, Oxford, and Carlton).

* Legendre, P. and Legendre, L. (2012) *Numerical Ecology*. Developments in Environmental Modelling no. xxiv, 3rd ed., 990 pp. (Elsevier: Amsterdam and Oxford).

## Bootstrapping_CI_Function.r v. 3.0

Bootstrapping is just another way to calculate confidence intervals of data means, when standard approaches are not useful due to some reason (e.g. data are not normally distributed). This code enables the calculation of a confidence interval of the mean of several data values to any chosen degree of significance. Although mathematically bootstrapping already works from three measurement values per sample, larger sample sizes are strictly recommended. The code includes a self adaptive approach, that uses accelerated bootstrapping when possible and falls back to basic bootstrapping if the former one is not suitable (Dixon 2002). This self-adaptation can be turned off, to use basic bootstrapping all the time. As of v. 3.0 the standard deviation of the sample values including its 95% confidence interval (Sheskin 2011) is also calculated. If you prefer to use MS Excel, rather than R, you should have a look at [resample.xls](http://woodm.myweb.port.ac.uk/nms/), which basically does nearly the same as this code (though it only allows basic bootstrapping), though it involves much more manual work for large datasets.

* Dixon, Ph. M. (2002) Bootstrap resampling. in El-Shaarawi, A. H. and Piegorsch, W. W. (eds) *Encyclopedia of Environmetrics*, pp. 212–20 (John Wiley & Sons, Ltd.: Chichester).

* Sheskin, D. J. (2011) *Handbook of Parametric and Nonparametric Statistical Procedures*. 5th ed., 1886 pp. (Chapman & Hall/CRC Press: Boca Raton, London, New York). 

## CoefficientVariation_Function.r v. 1.1

The coefficient of variation can be used to assess how variable a population is in any parameter. It is defined as the ratio between the standard deviation and the mean of the population. Many methods to calculate confidence intervals for the coefficient of variation exist, and two such methods (McKay 1932, Vangel 1996) are implemented in the program. The code reads a set of parameters which are somehow grouped and calculates the by-group coefficient of variation including its confidence interval for each parameter.

* McKay, A. T. (1932) Distribution of the coefficient of variation and the extended "t" distribution. *Journal of the Royal Statistical Society* 95 (4): 695-8.

* Vangel, M. G. (1996) Confidence intervals for a normal coefficient of variation. *The American Statistician* 50 (1): 21–6. 

## CoordinateTranslation_Function.r v. 2.1.1

Everybody working in the field of geosciences knows the problem with different systems of geographical coordinates. While many people still prefer the degrees–minute–second scheme, many programs require decimal degrees. Some sources even mix both systems and provide coordinates as degrees–minute.minute. The problem becomes even more immanent if someone compiles data from different sources and wants to convert them into one scheme. This code provides the means to convert a set of coordinates provided in any of the aforementioned systems in any other of those systems. I did not find a way to invoke something like an automated format recognition (and I doubt if that is possible), so that the data encoding must be manually entered in the first or second column of the dataset. 

##  CurveDiscussion_Function.r v. 1.0

Sometimes, one can fit a complex yet explanatory valid function to an observation of data. Oftentimes it is then necessary to work with parameters of that fitted function. While it is comparatively easy to find extremal points of such a function or calculate the goodness of fit, some tasks are harder to perform. The curvature function is a well known equation that allows to find the point of maximum curvature in a function, but it is rather difficult to calculate. Partly so, because it requires complex mathematics with the first and second derivative of the original function, and principally also the resulting curvature function. The R code provided calculates and returns the curvature function of any input function. It also iteratively calculates the position of maximum curvature in the input function. While this iterative approach is principally inferior to the mathematical approach using derivatives, it circumvents some problems and can still principally deliver results with a virtually infinite accuracy.

## DiagramExtraction_Function.r v. 1.3

Most scientists know the problems with the availability of raw data (or better, the lack thereof), especially for older works. Some authors may have presented results that are quite useful for ones own work, but one has only access to the diagrams. Often, in such cases, one wants to reconstruct the data in digitized form, to be able to work with them more fluently. Sure that is possible in many advanced image processing software tools, but those can cause problems, for instance, when *x*- and *y*-axis are not to the same scale. And being able to create a table of *x*–*y*-coordinates in any given program does not necessarily mean, that one can simply export that table (Personal experience of the author: Sometimes the only way is, to manually copy&paste each and every value one by one into a spreadsheet.). One can also perform that task by hand and manually punch in the data into an Excel sheet or something, but that requires a high degree of precision, some manual calculations, and is overall rather time consuming. This R code provides a function to plot a 2D-diagram or geographic map as image, set the scale for both axes separately, digitize the points or lines, and get an automatically exported table with *x*- and *y*-coordinates of all digitized points as .txt file. As image format, .ppm, .jpg, .tif, and .png are all supported. A conversion wrapper for ImageMagic which runs directly in R is included.

## GCalida_Morphospace_Projector_Function.r v. 1.0

This program can be used to project morphometric data extracted from specimens of the planktonic foraminifer genus *Globigerinella* into a predefined morphospace separating the three species *G. siphonifera*, *G. calida*, and *G. radians*. It can thus be used to objectively separate between those species. For further details compare Weiner et al. (2015).

**Cite as:** Weiner, A. K. M., Weinkauf, M. F. G., Kurasawa, A., Darling, K. F., and Kučera, M. (2015) Genetic and morphometric evidence for parallel evolution of the Globigerinella calida morphotype. Marine Micropaleontology 114: 19–35. doi:10.1016/j.marmicro.2014.10.003

## GeographicCalculator_Function.r v. 1.0

This file cointains functions that have to do with geographical calculations. It will probably be expanded in the future, but for now it only contains a function that allows to calculate the corner coordinates of a square of a given size around a given center coordinate.

## GeometricMorphometrics_Functions.r v. 1.11

Geometric morphometrics are one of the two major branches of morphometrics (traditional morphometrics being the other). While traditional morphometrics use only a few selected measurements to describe the morphology of an object, geometric morphometrics try to capture the whole picture of morphology. Geometric morphometrics relies on the definition of relatively few, well chosen landmarks within the object to perform that task and is the newest branch of morphometrics analyses. It tries to capture the whole shape (i.e. size-independent form) of objects using a coherent set of landmarks, instead of doing so using a combination of more or less arbitrarily chosen linear measurements. It can thus describe general shape change independent of size (and is therefore not affected by scaling problems). On the other hand it concentrates on only a few well chosen features of the object, instead of extracting the whole outline regardless of their local explanatory value. It therefore occupies a middle ground between outline analyses and traditional morphometrics. In contrast to both other methods, however, landmarks are not independent of each other after fitting, so that traditional statistic approaches cannot be aplied without modifications to landmark data. The script provided here offers a variety of geometric morphometric analysis tools. Other function files allow extraction of landmarks, reading and writing of .nts and .tps files as intial steps. This function file is complementing, not replacing, the R-packages shapes, Momocs, and geomorph. Apart from R, [MorphoJ](http://www.flywings.org.uk/morphoj_page.htm) is very versatile alternative for landmarks analyses. Most of the functions are based on Claude (2008) and Zelditch et al. (2012), but have been reworked/reassembled to obtain an even higher degree of automation.

* Claude, J. (2008) *Morphometrics with R*. Gentleman, R., Hornik, K., and Parmigiani, G. (eds) Use R! no. ii, 316 pp. (Springer).

* Zelditch, M. L., Swiderski, D. L., and Sheets, H. D. (2012) *Geometric Morphometrics for Biologists: A Primer*. 2nd ed., 478 pp. (London, Waltham, San Diego: Academic Press).

## ManualPermutation_CI_Function.r v. 1.4

Sometimes, for instance when calculating the Measurement Based Weight of organisms (Beer et al. 2010), one only retrieves one data value per sample, which in fact already is the mean of several organisms. Such procedures are often used, when an analysis of each specimen solitarily would be too time consuming. If such a per-specimen analysis is principally possible, however, this code provides the necessary tools to estimate the confidence interval of your data. If you perform a Monte-Carlo permutation, measuring only a random subset of at least one sample several times, those data reflect the variability of your organisms, and can thus be used for confidence interval calculations. The first part of the provided code does exactly that, using the Q8(p) definition for the quantiles (Hyndman and Fan 1996). The second part of the code can be used, to include those calculated confidence intervals in correlation analyses - which would be hard to do using normal approaches. For that, a randomization approach is applied, that in each rerun chooses the value used for the correlation randomly from within the calculated confidence interval.

* Beer, Ch. J., Schiebel, R., and Wilson, P. A. (2010) Technical note: On methodologies for determining the size-normalised weight of planktic Foraminifera. *Biogeosciences* 7: 2193-8.

* Hyndman, R. J. and Fan, Y. (1996) Sample quantiles in statistical packages. *The American Statistician* 50 (4): 361-5.

**Cite as:** Weinkauf, M. F. G., Moller, T., Koch, M. C., and Kučera, M. (2013) Calcification intensity in planktonic Foraminifera reflects ambient conditions irrespective of environmental stress. *Biogeosciences* 10: 6639-55. doi:10.5194/bg-10-6639-2013 

## MorphoFiles_Function.r v. 1.3

This file provides functions that can be used to read and write from morphometric standard file formats (.nts, .tps) into a shapes object in R and vice versa.

## MorphometricExtraction_Functions.r v. 1.1.2

This file provides functions to extract morphometric data from images.

## MultipleTestingCorrections_Function.r v. 1.0

A simple function to recalculate levels of significance for multiple testing. Can be interesting sometimes, but mostly I would recommend the R-function p.adjust().

## Oceanography_Function.r v. 2.0

This is a collection of functions for oceanographic data analysis. It contains (1) a conversion function to translate conductivity values of seawater (under known temperature and pressure) into salinity values corresponding to the practical salinity unit scale. (2) A function to calculate seawater alkalinity from titration with HCl.

## OutlineAnalysis_Functions.r v. 2.2

Outline analyses is a subcategory of geometric morphometrics, one of the two major branches of morphometrics (traditional morphometrics being the other). While traditional morphometrics use only a few selected measurements to describe the morphology of an object, outline analyses tries to capture the whole picture of morphology. Outline analyses use the digitized outline of the object to describe its shape. While outline analyses have been criticised because they neglect features within the object by exclusively focusing on the outline and mostly involve more or less drastic mathematic recalculations, they have some advantages as well. They can capture shape on a more objective level than any other method, because all of the other methods rely on the more or less subjective definition of measurement points. They can describe form (in terms of outer shape) better than any other method. They allow a complete, smooth reconstruction of any possible shape and can model that shape as a whole, not being limited to the reconstruction of the position of a few points. Lastly, in contrast to geometric morphometrics, results can be analysed using traditional multivariate approaches without further parameter tweaking. Outline analyses is particularly problematic in practice, because much of the software is rather old and does hardly run on modern machines. Furthermore, the software is heavily fragmented, often necessitating to switch between four or five programs repeatedly to perform a complete analysis. The program [SHAPE](http://lbm.ab.a.u-tokyo.ac.jp/~iwata/shape/index.html) is one of the few noteworthy programs that run perfectly well on modern computers and are able to perform most tasks (from data acquisition to analysis) in a single framework. The script provided here aims for the same goal in the more adaptable R-environment, containing all functions necessary to extract outlines and perform a Zahn–Roskies Fourier analysis and an elliptic Fourier analysis. Most of the functions are based on Claude (2008), but have been reworked/reassembled to obtain an even higher degree of automation.

* Claude, J. (2008) *Morphometrics with R*. Gentleman, R., Hornik, K., and Parmigiani, G. (eds) Use R! no. ii, 316 pp. (Springer). 

## RegressionTools_Functions.r v. 2.0

Linear regression in general is one of the most abused methods of statistical analyses, given that many people do not understand the fundamental (and partly philosophical) difference between correlation and regression, and simply go for regression (even worse, often far from linear) because a line 'looks neat'. Matters grow worse when one is counting how often model I linear regression is used on datasets that are simply not suitable for that analysis. Beside the obvious prerequisite that a straight line describes the monotonic relationship between dependent and independent variable best, there are several other assumptions made: (a) *x*-values are assumed to be measured basically without error (as far as that is possible in reality), (b) *x*-values were chosen by the experimenter, and (c) *y*-values are normally distributed with the same variance for all values of *x*. If any of those assumptions is violated, a model I linear regression is not suitable, and should be replaced by a model II or model III linear regression method. Given prerequisites for model I linear regression it is mainly of use for laboratory experiments, but especially points (a) and (b) will seldom hold true for data collected in the field or from the fossil record. Model II linear regression is available in R via the lmodel2 package, but only for the bivariate case. Here, I implemented multivariate model II linear regression as described by Richter and Stavn (2014). The Kendall-Theil robust line fitting method (Kendall 1938, Theil 1950, Sen 1968) invoked here is one of the most robust of such model III linear regressions that should be used instead in such cases—and nevertheless still not implemented in any program I know except [KTRLine](http://pubs.usgs.gov/tm/2006/tm4a7/). Further, a function to calculate the confidence band around a linear regression line has been made available.

* Kendall, M. G. (1938) A New Measurement of Rank Correlation. *Biometrika* 30 (1-2): 81-93.

* Richter, S. J. and Stavn, R. H. (2014) Determining functional relations in multivariate oceanographic systems: Model II multiple linear regression. *Journal of Atmospheric and Oceanic Technology* 31: 1663-72.

* Sen, P. K. (1968) Estimates of the Regression Coefficient Based on Kendall's Tau. *Journal of the American Statistical Association* 63 (324): 1379-89.

* Theil, H. (1950) A Rank-Invariant Method of Linear and Polynomial Regression Analysis, III. *Proceedings of the Koninklijke Nederlandse Akademie van Wetenschappen* 53 (9): 1397-412.

## SpiralGrowthAnalysis_Functions v. 1.3

This function provides tools to extract growth data on spirally growing organisms (Foraminifera, Ammonoidea, etc.) and analyse those data using a robust linear regression. Its purpose is to test the deviation of growth from spirality, to estimate to what degree the growth follows an ideal logarithmic spiral.

## StukelsTest_Function v. 1.0

It can sometimes be difficult to calculate the goodness of fit (i.e. the ability of the model to describe the data well) for binomial models. While methods like Pearson's chi-square test work when the data can be grouped into strata, they fail as soon as either any of the strata contains less than *c.*5 cases or at least one of the independent variables is continuous, prohibiting the erection of strata altogether. For such cases, the test invented by Stukel (1988) is a versatile alternative to test, how well the model describes the data.

* Stukel, Th. A. (1988) Generalized logistic models. *Journal of the American Statistical Association* 83 (402): 426–31.

## ToGMT_Function v. 1.1

This file provides a couple of functions to transform data from other sources to be used in the Generic Mapping Tools. Currently works with: World Ocean Atlas .csv data, Ocean Data View spreadsheets (.txt files).

## TraditionalMorphometrics_Functions v. 1.3

Traditional morphometrics are one of the two major branches of morphometrics (geometric morphometrics being the other). Traditional morphometrics use only a few selected measurements to describe the morphology of an object, which is why it is oftentimes neglected nowadays. What it lacks in descriptive power, however, it makes up for in understandability and time-effectiveness of data-gathering. The script provided here includes some of the main approaches used in traditional morphometric analyses. Most of the functions are based on Claude (2008), but have been reworked/reassembled to obtain an even higher degree of automation.

* Claude, J. (2008) Morphometrics with R. Gentleman, R., Hornik, K., and Parmigiani, G. (eds) Use R!, vol. 2, 316 pp. (Springer). 

# Author

Manuel F. G. Weinkauf, Universite de Geneve, Departement des sciences de la Terre, 1205 Geneve, Switzerland; Manuel.Weinkauf@unige.ch















