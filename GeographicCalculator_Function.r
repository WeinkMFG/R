#Set of functions to calculate distances and coordinates on the surface of the Earth
#Further reading: http://www.movable-type.co.uk/scripts/latlong.html	

#Author: Manuel Weinkauf (Manuel.Weinkauf@unige.ch)
#Version: 1.0
#Date: 13 August 2014

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.#
#To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/.                   #
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#**************************************************************************************

#########################################################################
# Function to calculate corner coordinates of rectangle around given... #
#    center point                                                       #
# Necessary packages: maptools, raster                                  #
# Necessary input variables:                                            #
#    Start: Latitude and Longitude of starting point in decimal...      #
#           degrees (-90:90°, -180:180°).                               #
#          *numeric (real)*, length=2                                   #
#    Dist: Distance from starting point into directions N, E, S, W...   #
#          (in this order) in km.                                       #
#          *numeric (real)*, length=4                                   #
# Output data: Map and coordinates of rectangle borders towards N, E,...#
#              S, and W around starting point.                          #
# Input data: Information about position and size of area.              #
#########################################################################

#Loading packages
require(maptools)
data(wrld_simpl)
require(raster)

GeoRect<-function (Start, Dist) {
	#Check data for consistency
	if (length(Start)!=2) {stop("Start coordinates (Lat, Lon) must be provided as vector of length 2")}
	if (length(Dist)!=4) {stop("Distances to N, E, S, and W must be provided as vector of length 4")}
	if (Start[1]<(-90) | Start[1]>90) {stop("Start latitude must lie between -90 and 90")}
	if (Start[2]<(-180) | Start[2]>180) {stop("Start longitude must lie between -180 and 180")}
	
	#Convert Start from degrees to radians
	Start<-Start*(pi/180)
	
	#Calculate parameters of rectangle
	Params<-matrix(NA, 2, 4)
	##Bearings
	Params[1,]<-c(0*(pi/180), 90*(pi/180), 180*(pi/180), 270*(pi/180))
	##Distances
	Params[2,]<-Dist

	#Calculate latitude and longitude of rectangle borders
	Borders<-vector(mode="numeric", length=4)
	names(Borders)<-c("Lat.N", "Lon.E", "Lat.S", "Lon.W")
	##To North
	Borders[1]<-asin((sin(Start[1])*cos(Params[2,1]/6371))+(cos(Start[1])*sin(Params[2,1]/6371)*cos(Params[1,1])))*(180/pi)

	##To East
	Aux<-asin((sin(Start[1])*cos(Params[2,2]/6371))+(cos(Start[1])*sin(Params[2,2]/6371)*cos(Params[1,2])))*(180/pi)
	Borders[2]<-(Start[2]+atan2(sin(Params[1,2])*sin(Params[2,2]/6371)*cos(Start[1]), cos(Params[2,2]/6371)-(sin(Start[1])*sin(Aux))))*(180/pi)

	##To South
	Borders[3]<-asin((sin(Start[1])*cos(Params[2,3]/6371))+(cos(Start[1])*sin(Params[2,3]/6371)*cos(Params[1,3])))*(180/pi)

	##To West
	Aux<-asin((sin(Start[1])*cos(Params[2,4]/6371))+(cos(Start[1])*sin(Params[2,4]/6371)*cos(Params[1,4])))*(180/pi)
	Borders[4]<-(Start[2]+atan2(sin(Params[1,4])*sin(Params[2,4]/6371)*cos(Start[1]), cos(Params[2,4]/6371)-(sin(Start[1])*sin(Aux))))*(180/pi)

	##Plot map
	plot(wrld_simpl, xlim=c(Borders[4]-20, Borders[2]+20), ylim=c(Borders[3]-20, Borders[1]+20), axes=TRUE, col="light yellow")
	box()
	Area<-extent(Borders[4], Borders[2], Borders[3], Borders[1])
	plot(Area, add=TRUE, col="red", lwd=2)
	points(Start[2]*(180/pi), Start[1]*(180/pi), pch=23, col="darkgreen", bg="yellow", cex=2)
	
	##Return coordinates
	return(Borders)
}

#--------------------------------------------

#Example
#GeoRect(c(33, -22), c(300, 150, 150, 300))

#--------------------------------------------
#--------------------------------------------
#Version History
#1.0	Finished Program
#--------------------------------------------
#--------------------------------------------
