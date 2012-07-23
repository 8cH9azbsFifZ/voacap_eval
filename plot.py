#!/opt/local/bin/python2.5
#
#    Copyright: Gerolf Ziegenhain <g@ziegenhain.com> 2012
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

#PYTHONPATH=/opt/local//lib/python2.5/site-packages/ python2.5 ./plot.py 
# FIXME: document with doxygen functions variables etc
import sys
sys.path.append("/opt/local//lib/python2.5/site-packages/")
from numpy import *
import numpy as np # FIXME namespace
from pylab import *
import re
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import os
import datetime 

class voaarea:
	"""This is the voaarea class
	It conducts the actual simulations
	"""
	# Grid parameters
	__min_lon, __max_lon = -180, 180
	__min_lat, __max_lat = -90, 90
	__gridsize = 50 

	# Extremal values
	vsnr_min, vsnr_max = 0, 70
	vrel_min, vrel_max = 0, 1
	
	# Simulation Configuration Parameters
	ssn = 72.8 # FIXME read from file
			# ftp://ftp.ngdc.noaa.gov/STP/SOLAR_DATA/SUNSPOT_NUMBERS/sunspot.predict
			# ssnFetch.py 

	# FIXME: remove my stuff
	def __init__(self,frequency,power,month,hour,
			tx_label="MYCALL", tx_lat=1.1, tx_lon=2.2, 
			filepath="MYHOME/itshfbc/areadata/",
			RSN=27, XNOISE = 145,
			tx_ant_data_file = "samples/sample.00",
			rx_ant_data_file = "samples/sample.00",
			): 
		""" This is the constructor, which
		Input Parameters:
		- RSN: Required SNR for circuit reliability calculation
				(good broadcast station 67, cw 24 or 27, ssb 45, cw as on voacaponline 16)
		- XNOISE: man made noise 
				145 in residental, quiet environment: 155
		- tx_ant_data_file and rx_ant_data_file:
				Standard: isotrope antenna
		What does this contructor?
		- sets standard values
		- simulates if outputfiles (cf. pattern) do not exist
		- plots the data if outputfiles (cf. pattern) do not exist
		"""
		# Save parameters
		self.freq, self.power_watt = float(frequency), int(power)
		self.month, self.utc = int(month), hour
		self.tx_label, self.tx_lat, self.tx_lon = tx_label, tx_lat, tx_lon
		self.RSN, self.XNOISE = RSN, XNOISE
		self.filepath = filepath
		self.tx_ant_data_file, self.rx_ant_data_file = tx_ant_data_file, rx_ant_data_file 
		#tx_ant_data_file = "gz/gzpole" # norcal doublet antenna
		#tx_ant_data_file = "samples/sample.23" # dipole

		# Configure dependent variables
		self.filename = "_%.1f_%d_%d_%d" %(self.freq, self.power_watt, self.utc, self.month)
		self.title = self.tx_label+" %.1fMHz %dW %dh %dM" %(self.freq, self.power_watt, self.utc, self.month)
		self.voafilename = self.filepath+"/"+self.filename+".voa"
		self.vgfilename = self.filepath+"/"+self.filename+".vg1"
		self.snrimgname = "snr"+self.filename+".png"

		# Run the actual code
		if not os.path.isfile (self.voafilename):
			self.make_voacap_file()
		if not os.path.isfile (self.snrimgname):
			self.read_voacap_result()
			self.__plot_snr()
		#self.__plot_rel()

	def __lat_as_string(self,lat):
		""" 
		Function for conversion of format in voacap
		"""
		if lat > 90.0 : lat = 90.0
		if lat < -90.0 : lat = -90.0
		lat_sign = 'N'
		if lat < 0.0:
			lat_sign = 'S'
		return "%5.2f%s" % (abs(lat), lat_sign)

	def __lon_as_string(self,lon):
		""" 
		Function for conversion of format in voacap
		"""
		if lon > 180.0 : lon = 180.0
		if lon < -180.0 : lon = -180.0
		lon_sign = 'E'
		if lon < 0.0:
			lon_sign = 'W'
		return "%6.2f%s" % (abs(lon), lon_sign)

	def make_voacap_file(self):
		"""
		Create a voacap area simulation file
		"""
		#FIXME comments below
		# Some more system configuration - currently I do not use it
		AMIND = 3 #amind = minimum takeoff angle - 3 bad antenna -> better 1
		XLUFP = 90 #xlufp requiered circuit reliability (keep it at >90%)
		PMP = 3 #pmp multipath power tolerance - accept deatul 3db
		DMPX = 0.1 #dmpx - maximal tolerable time delay - accept default 0.10 ms
		PSC1, PSC2, PSC3, PSC4 = 1, 1, 1, 1 # FIXME: what is it?
		tx_ant_design_freq, tx_ant_bearing = 0, 0
		rx_ant_gain, rx_ant_bearing = 0, 0
		tx_ant_power = self.power_watt*0.001 # power (kW)

		# Create the voacap simulation file
		voafile = open (self.voafilename, "w")
		print >>voafile, "Model    :VOACAP"
		print >>voafile, "Colors   :Black    :Blue     :Ignore   :Ignore   :Red      :Black with shading"
		print >>voafile, "Cities   :Receive.cty"
		print >>voafile, "Nparms   :    4"
		print >>voafile, "Parameter:MUF      0"
		print >>voafile, "Parameter:DBU      0"
		print >>voafile, "Parameter:SNRxx    0"
		print >>voafile, "Parameter:REL      0"
		print >>voafile, "Transmit :%10s%10s%18sShort" % (self.__lat_as_string(self.tx_lat), self.__lon_as_string(self.tx_lon), self.tx_label)
		print >>voafile, "Area     :%10.1f%10.1f%10.1f%10.1f" % (self.__min_lon, self.__max_lon, self.__min_lat, self.__max_lat)
		print >>voafile, "Gridsize :  %3d    1" % (self.__gridsize)
		print >>voafile, "Method   :   30" # Short/Long Path smoothing for over 7000km
		#print >>voafile, "Method   :   20" # complete system performance
		print >>voafile, "Coeffs   :CCIR"
		print >>voafile, "Months   :%7.2f" % self.month
		print >>voafile, "Ssns     :%7d" % self.ssn
		print >>voafile, "Hours    :%7d" % self.utc
		print >>voafile, "Freqs    :%7.3f" % self.freq
		print >>voafile, "System   :%5d%10.3f%5d%5d%10.3f%10.3f" % (self.XNOISE, AMIND, XLUFP, self.RSN, PMP, DMPX)
		print >>voafile, "Fprob    :%5.2f%5.2f%5.2f%5.2f" % (PSC1, PSC2, PSC3, PSC4)
		print >>voafile, "Rec Ants :[%21s]  gain=%6.1f%6.1f" % (self.rx_ant_data_file, rx_ant_gain, rx_ant_bearing)
		print >>voafile, "Tx Ants  :[%21s]%7.3f%6.1f%10.4f" % (self.tx_ant_data_file, tx_ant_design_freq, tx_ant_bearing, tx_ant_power)
		voafile.close()
		os.system('voacapl ~/itshfbc area calc '+self.filename+'.voa')


	def read_voacap_result(self):
		""" 
		Read the results generated by voacap.
		"""
		# Abreviations
		gridsize = self.__gridsize
		min_lat, max_lat = self.__min_lat, self.__max_lat
		min_lon, max_lon = self.__min_lon, self.__max_lon

		## Read in the vg output file
		rel_title = "Circuit Reliability"
		# general information positions
		x_first_char, x_last_char = 3, 6
		y_first_char, y_last_char = 0, 3
		snr_first_char, snr_last_char = 86, 92
		rel_first_char, rel_last_char = 98, 104
		# generate grid field
		self.__lats, self.__lons = np.zeros (gridsize**2,float), np.zeros (gridsize**2,float)
		dlat, dlon = (max_lat-min_lat)/float(gridsize-1), (max_lon-min_lon)/float(gridsize-1)
		self.__lats = np.arange (min_lat, max_lat+0.001, dlat)
		self.__lons = np.arange (min_lon, max_lon+0.001, dlon)
		# generate the scalar data fields
		infile = open(self.vgfilename)
		pattern = re.compile(r"[a-z]+")
		self.snr = np.zeros ([gridsize,gridsize],float)
		self.rel = np.zeros ([gridsize,gridsize],float)
		for line in infile: 
			match = pattern.search (line)
			if not match:
				vsnr = float(line[int(snr_first_char):int(snr_last_char)])
				vrel = float(line[int(rel_first_char):int(rel_last_char)])
				vsnr = max (self.vsnr_min, vsnr)
				vsnr = min (self.vsnr_max, vsnr)
				vrel = max (self.vrel_min, vrel)
				vrel = min (self.vrel_max, vrel)
				gx, gy = int(line[x_first_char:x_last_char])-1, int(line[y_first_char:y_last_char])-1
				self.snr[gx][gy] = vsnr
				self.rel[gx][gy] = vrel

	def __plot_basemap(self):
		""" 
		The basis world map for all plots
		"""
		self.__map = Basemap(projection='cyl',lat_0=self.tx_lat,lon_0=self.tx_lon,resolution='i',
			llcrnrlon=self.__min_lon,llcrnrlat=self.__min_lat,urcrnrlon=self.__max_lon,urcrnrlat=self.__max_lat)
		self.__map.drawcoastlines(linewidth=0.25,color="black")
		self.__map.drawcountries(linewidth=0.25,color="black")
		self.__map.drawmapboundary(color="black",linewidth=1.0)

	# SNR
	def __snr_format(self, x, pos):
		"""
		Function for formatted labels in the SNR plot
		"""
		return '%3ddB' % x

	def __plot_snr(self):
		"""
		Plot the SNR map
		"""
		self.__plot_basemap()

		# Interpolate scalar field
		wsnr = np.zeros((self.__gridsize,self.__gridsize),float)
		wsnr, wsnr_lon, wsnr_lat = self.__map.transform_scalar(self.snr,self.__lons,self.__lats,
				self.__gridsize,self.__gridsize,returnxy=True,checkbounds=False,masked=True)
		wsnr = wsnr.filled(self.vsnr_min-1.0)

		c="hot"
		isnr = self.__map.imshow (wsnr, cmap=cm.get_cmap(c),
				extent = (-180,180,-90,90),#FIXME
				origin = "lower",
				norm = Normalize(clip=False,vmin=self.vsnr_min,vmax=self.vsnr_max)
					)
		ylabels = np.arange(self.RSN,self.vsnr_max,10)#(20, 30, 40, 50, 60, 70)
		cb = colorbar(isnr,orientation="horizontal",format = FuncFormatter(self.__snr_format))
		yticks = np.append(ylabels, [self.vsnr_min,self.vsnr_max])
		cb.set_ticks(yticks)
		for t in cb.ax.get_xticklabels():
			t.set_fontsize(8)
			t.set_rotation(45)
		cs = self.__map.contour(wsnr_lon,wsnr_lat,wsnr,ylabels,linewidth=0.25,alpha=0.5,colors="blue")
		#clabel(cs, fmt="%2d", colors="blue",fontsize=6,inline=1,inline_spacing=1)
		
		plt.savefig(self.snrimgname,dpi=100)
		plt.clf()

	# Reliability
	def __rel_format(self, x, pos):
		""" 
		Function for formatted labels in the circuit reliability plot
		"""
		return '%(percent)3d%%' % {'percent':x*100}

	def __plot_rel(self):
		"""
		Plot the circuit reliability
		"""
		self.__plot_basemap()

		# Interpolate scalar field
		wrel = np.zeros((self.__gridsize,self.__gridsize),float)
		wrel, wrel_lon, wrel_lat = self.__map.transform_scalar(self.rel,self.__lons,self.__lats,
				self.__gridsize,self.__gridsize,returnxy=True,checkbounds=False,masked=True)
		wrel = wrel.filled(self.vrel_min-1.0)

		irel = self.__map.imshow (wrel, cmap=cm.get_cmap('hot'),
				norm = Normalize(clip=False,vmin=self.vrel_min,vmax=self.vrel_max)
					)
		ylabels = (0,.5, .7, .9,1)
		cb = colorbar(irel,orientation="horizontal",format = FuncFormatter(self.__rel_format))
		yticks = ylabels
		cb.set_ticks(yticks)
		for t in cb.ax.get_xticklabels():
			t.set_fontsize(8)
			t.set_rotation(45)
		cs = self.__map.contour(wrel_lon,wrel_lat,wrel,ylabels,linewidth=0.25,colors="blue",alpha=0.5)
		
		imgname = "rel"+self.filename+".png"
		plt.savefig(imgname,dpi=100)
		plt.clf()

	def get_filename_snr(self):
		"""
		Return the filename of the created SNR image
		"""
		return self.snrimgname

# Run the code

def calc_qrp(curtime):
	"""
	Calculate the propagation conditions for all QRP frequencies
	"""
	outfile = "qrp."+str(curtime.hour)+"h.png"
	power = 10 #watt
	qrpfrequencies = (3.56, 7.03, 10.116, 14.06, 18.086, 21.06, 24.906, 28.06)
	bands = ("80m", "40m", "30m", "20m", "17m", "15m", "12m", "10m")
	filenames = ""
	# Run the simulations
	for freq in qrpfrequencies:
		a = voaarea(freq,power,curtime.month,curtime.hour)
		filenames = filenames+" "+a.get_filename_snr()
	n = len(qrpfrequencies)
	ff = filenames.split()
	# Image position parameters
	ex_lx, ex_ly = 624, 315
	ex_px, ex_py = 99, 84
	l_x, l_y = 105, 105
	tile = 3
	for i in arange(0,n): # Strip the images
		fex = ff[i].replace(".png","_ex.png")
		label = "-stroke white -fill white -pointsize 18 -draw \"text %d,%d \\\"%s\\\"\" " % (l_x, l_y, bands[i])
		cmd = "convert %s[%dx%d+%d+%d] %s %s" % (ff[i], ex_lx, ex_ly, ex_px, ex_py, label, fex)
		os.system(cmd)
	# Combine the images
	cmd = "montage -tile %dx%d -geometry %dx%d %s %s" % (tile, tile, ex_lx, ex_ly, filenames.replace(".png","_ex.png"), outfile)
	os.system(cmd)
	# Label
	ll_x = (tile-1)*ex_lx+(l_x-ex_px)
	ll_y = (tile-1)*ex_ly+(l_y-ex_py)
	param = "Hour: %s UTC\nMonth: %s\nPower: %dW" % (curtime.hour, curtime.month, power)
	label = "-stroke black -fill black -pointsize 18 -draw \"text %d,%d \\\"%s\\\"\" " % (ll_x, ll_y, param)
	cmd = "mogrify %s %s" % (label, outfile)
	os.system(cmd)


# Run actual code
curtime = datetime.datetime.utcnow()
for h in (range(0,24)):
	calc_qrp(curtime + datetime.timedelta(hours=h) )

	# sort by date for google images
	# for i in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24; do j=`printf %02d $i`;  touch -t 0722${j}01 qrp.${i}h.png;done

