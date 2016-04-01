#!/usr/bin/env python

help = '''
Photometric De-Trending (pdtrend, also Pavlos-Daewon de-Trending) pipeline for astronomical time series. See details in Kim et al. 2008.


Packages to run this pipeline :
	matplotlib : http://matplotlib.sourceforge.net/
	numpy : http://numpy.scipy.org/
	Pycluster : http://bonsai.ims.u-tokyo.ac.jp/~mdehoon/software/cluster/
	hcluster : http://code.google.com/p/scipy-cluster/
	R : http://www.r-project.org/
		You need to install 'nortest' abd 'quadprog' package at R. Use 'install.packages' command in R shell.
	Rpy : http://rpy.sourceforge.net/
	

Usage : ./pdtrending.py -f \'file_list\' [OPTION...]
	
	-f : list of files which you want to de-trend. You can use regular expression. Don't forget for surrounding file_list with \'. i.e. -f \'/home/astro/hht_[0-9]*.lc\'
			
	-c : Correlation matrix file if you have one. i.e. -c /home/astro/corr.dat
			
	-o : Output directory. Default is './pdtrend' under the directory where light-curves are. i.e. -o /home/astro/pdtrend/
			
	-l 0 : Column number of flux in light-curves. Default is 0 (first column).
			
	-i 1 : Not implemented in this version. Number of iteration of de-trending. If no clusters are found, the iteration will stop regardless the number of iteration. Note that png files of determined-trends trends-position will be generated only at the first iteration. Default is 1, which means no iteration but single execution of de-trending. You cannot use this option with -c together.
			
	-d : X and Y coordinates of light-curves. First(second) column should be X(Y) coordinates. The order of stars should be in same order with file_list. i.e. -c /home/astro/xy.dat
			
	-b [True|False] : True if light-curves are in flux based. If they are in magnitude based, then False. Default is True.
			
	-z 20 : Zero magnitude. Default is 20.
			
	-n 0 : Binning light-curves just when calculating correlation values, which can help to increase signal-to-noise ratio. The binning might be useful when low frequency trends are contaminated by high frequency noise source. Default is 0, which means no binning. Note that light-curves will be de-trended in original sampling rate (without binning).
			
	-t 2 : Constrain initial seeds. Default is 2. From 1 to 10 is suitable. If no clusters are found, decrease this number. It might help to find clusters. See Kim et al. 2008 for more details.
	
	-g 0.1 : level of significance, p-value. Default is 0.1. Smaller p-value generally gives more clusters.
			
	-m 5 : Mimimum number of stars in a cluster. Default is 5. According to the number of stars in dataset, you can adjust this value.
			
	-s [True|False] : If True, you can see a dendrogram of hierarchical tree of similarity matrix derived using light-curves. Default is False.
			
	-h : Show this help.
	

mailto : dakim@cfa.harvard.edu
'''

change_log = '''
@ Change log :
		v1.00 at 01/05/2009 : Just change the version number to 1.00 to publish this package.
		
		v4.00 at 11/14/2008 : We gave up multiple linear regression which occassionally destroys intrinsic signals.
		Instead of that, we implemented quadratic program which can constraint free parameters to be bigger than or equal to 0.
		This greatly reduces the signal destruction effects. Another big jump of our algorithm.
		This is our last improvement at this version of De-trending algorithm.
		Next version might start from handling the phase shift of trends.
		
		v3.51 at 09/22/2008 : Save the standard deviation values of light curves after- and before de-trending.
		
		v3.50 at 09/15/2008 : Fisher's tranformation is applied to the distances lists of subsets.
		The transformation convert skewed distribution of correlations to normal distribution;
		which strengthen our normality assumption. Nevertheless, the resulting clusters are almost
		identical because we gather only highly correlated elements and cut outliers by p-value.
		This is the thrid big jump of our algorithm.
		
		v3.03 at 08/25/2008 : Remove Shapiro-Wilk test. Use only Anderson-Darling test.
		
		v3.01 at 08/19/2008 : Save the determined trends light curves as well; as 'trend_[0-9][0-9]'
		
		v3.00 at 08/13/2008 : Clustering algorithm changed to the normality test algorithm.
		Shapiro-Wilk test (Shapiro & Wilk 1965) and Anderson-Darling test (Anderson & Darling 1952)
		is used to test the normality of cluster. Many simulation had performed again.
		New algorithm idenfity the clusters much better than before.
		The second big jump of our de-trending algorithm.
		
		v2.21 at 06/10/2008 : Software is packaged at the first time. Remove not-used function.
		Make simulation light curves for test run. Packaging the light curves separately.
		
		v2.20 at 06/01/2008 : Improvement of combining template set routine to build trend,
		we use weighting according to the standard deviation of normalized light curves.
		
		v2.10 at 05/28/2008 : Improvement of removal of trends routine. Based on the constructed trends,
		Applied multiple regression method.
		
		v2.01 at 05/22/2008 : Save constructed trends light curves as png image. If xy_coords is not null,
		then also save position of each trend on CCD plane as png image.
		
		v2.00 at 05/21/2008 : Change block-clustering algorithm to hierachichal tree clustering algorithm.
		It's one of the biggest changes in our algorithm.
		
		v1.00 was finished at 05/08/2008 by Dae-Won Kim.
'''

import getopt
#import matplotlib
#matplotlib.use('agg')
from pylab import *
from numpy import *
import os
import sys
import glob

from Pycluster import *
from hcluster import *
from clusterlib import *

sys.path.append('./')

def pdetrending(file_list, correlation_matrix='', xy_coords='', output_dir='', \
img=True, show_tree=False, min_template_stars=0, std_save=True, \
flux_based=True, zero_mag=20., column_flux=0, bin_window=0, initial_seed=2, \
n_iteration=1, l_significance=0.1):
	'''
	flie_list :
		list of files. This list should be sorted by brightness of stars.
		Each file has to contain flux list of a star at the first column.
	
	correlation_matrix = '' : 
		Correlation matrix of stars in file_list. If this is '', it will be generated by our routine.
		Newly generated matrix will be saved under the output_dir as 'corr.dat'
	
	xy_coords = '' :
		x and y coordinates of star on CCD plane. This list shoud be sorted by brightness of stars either.
		Therefore total number of lines should be the same with total number of lines of file_list.
		The first column should be x coordinate. The second column should be y coordinate.
		If this coordinates file is provided and img=True, then we save the png file which show where each
		constructed trends are placed on CCD plane. It's named 'trends_position.png' in output_dir.
		
	output_dir = '' : 
		Output directory where de-trended light curves will be saved.
		If this is '', de-trended light curves will be saved under the directory of file_list[0], 'pdtrend'
	
	img = False : 
		Save the image of de-trended light curves as png files. N number of png files will be generated
		under output_dir. N is the total number of stars.
		Top panel is raw, middle is constructed trend and bottom panel is de-trended resules.
		Also constructed trends image file will be saved as 'constructed_trends.png'
	
	show_tree = False :
		If it's true, this routine shows the dendrogram constructed based on the correlation matrix.
	
	min_template_stars = 0:
		Minimum number of template stars per cluster. Empirically, about (total number of star) / 50 is acceptable.
		However, if this value is smaller than 5, it means that constructed trends might suffer from noise
		because the number of template star is too small. If this value is set to 0, then (total number of star) / 50
		is used.
	
	std_save = True:
		Save stadnard deviation values of light curves (before and after de-trending).
	
	flux_based = True:
		If data is in magnitude, change this values to False and set zero_mag to proper values.
		
	zero_mag = 20:
		Only activated when flux_based = False.
	
	column_flux = 0:
		Column of flux in light-curve files.
	
	bin_window = 0:
		Bin light-curves when calculating correlation values. This can help to find clusters when trends are contaminated by noise a lot. 0 means no binning.
		
	initial_seed = 2:
		Constraint initial seeds. For more details, see Kim et al. 2008
		
	n_iteration = 3:
		Number of iteration of de-trending processes. Not implemented in this version.
	'''
	
	#Append current directory for python libraries
	sys.path.append('./')
	
	#Check initial values.
	if file_list == '':
		print '#ERROR : Check the file_list.'
		sys.exit()
	if output_dir == '':
		output_dir = os.path.dirname(file_list[0]) + '/pdtrend/'
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)
	
	#read file into memory.
	if xy_coords != '':
		print '#Start to read xy coordinates..'
		whole_xy = []
		fp = open(xy_coords)
		whole_xy = array([map(float, line.split()) for line in fp])
		max_x = max(whole_xy[:,0])
		max_y = max(whole_xy[:,1])
		fp.close()
		if len(whole_xy) != len(file_list):
			print '#ERROR : The number of xy coords file is not matched with the number of stars.'
			sys.exit()
		
	print '#Start to read file..'
	whole_lc = []
	for i in range(len(file_list)):
		print '	#Reading ' + file_list[i] + '..'
		lc = []
		fp = open(file_list[i])
		for line in fp:
			if line[0] == '#':
				continue
			flux_column = float(line.split()[column_flux])
			if not flux_based:
				flux = math.pow(10, (flux_column - zero_mag) / -2.5)
			else:
				flux = flux_column
			lc.append(flux)
		fp.close()
		whole_lc.append(lc)
	whole_lc = array(whole_lc)
	
	#Check correlation matrix.
	if correlation_matrix == '':
		print '#Calculate correlation matrix..'
		corr_list = cal_correlation_matrix(whole_lc, bin_window)
		
		print '#Save correlation matrix to local file to %s..' % (output_dir + '/corr.dat')
		fp = open(output_dir + '/corr.dat', 'wb')
		for ele in corr_list:
			fp.write(' '.join(map(str, ele)) + '\n')
		fp.close()
	else:
		print '#Read correlation matrix..'
		fp = open(correlation_matrix)
		corr_list = array([map(float, line.split()) for line in fp])
		if corr_list.shape != (len(whole_lc), len(whole_lc)):
			print '#ERROR : Dimension of correlation matrix is not match!. \
Execute this routine without correlation matrix value and generate correlation matrix again.'
			sys.exit()
	
	#select template clusters based on the hierarchical tree structure.
	#for more details, see Kim et al. 2008.
	print '#Now finding clusters for template set..'
	dist_matrix = 1. - corr_list
	print '	#Making Hierarchical tree..'
	tree = treecluster(method='m', distancematrix=dist_matrix.copy())
	if show_tree:
		R = dendrogram(Pycluster_to_hcluster(tree), leaf_font_size=10)
		ylabel('Distance'); xlabel('Index of star'); show()
		
	print '	#Finding clusters in tree..'
	groups = find_group_DW(tree, dist_matrix, len(whole_lc[0]), initial_seed, l_significance)
	print '	#Remove subset of clusters..'
	groups = remove_subset_cluster(groups)
	groups = remove_small(groups, min_template_stars - 1)
	print '	#Total %d clusters found..' % (len(groups))
	#print groups
	for group in groups:
		print '#----------------------------#'
		for ele in group:
			print file_list[ele]
	
	if len(groups) == 0:
		print '#ERROR : No cluster is found; which means this dataset has no \
		reliable template set. You have to check below things. \
		1) Are your data contaminated by other noise sources too much? if yes, retry this algorithm after data binning. \
		2) Do you have enough number of stars? if yes, there is no solution. Our algorithm is \
		developed for wide field survey; which means you should have enough number of stars, \
		~ a few hundreds. \
		3) Do your data have trends? In other words, do your stars just show independent trends with others? \
		Check the light curves by eyes. If the trends appear just in tiny number of stars, \
		then it\'s very hard to find reliable clusters. However, you could retry this algorithm after data binning, \
		which might help to find cluters.'
		sys.exit()
	
	#make trend list of each group.
	print '#Construct trends with identified clusters..'
	trend_set = []
	for group in groups:
		trend_set.append(create_trend(group, whole_lc))
	trend_set = array(trend_set)
		
	if img:
		clf()
		color=['r', 'g', 'b', 'y', 'k']
		#color = ['0.5', '0.6', '0.7', '0.8', '0.9']
		shape=['s', '^', 'o', 'D', '+', 'x', '1', 'h', 'p']
		print '	#Saving constructed trends as png files..'
		for i in range(len(trend_set)):
			ax = subplot(len(trend_set), 1, i + 1)
			#plot(trend_set[i] + 1., color[i%len(color)]+shape[i%len(shape)])
			plot(trend_set[i] + 1., 'b,')
			ylabel('%d' % len(groups[i]))
		xlabel('Time index')
		savefig(output_dir + '/constructed_trends.png')
		
		if xy_coords != '':
			clf()
			print '	#Saving position of trends on CCD as png files..'
			for i in range(len(groups)):
				for ele in groups[i]:
					plot([whole_xy[ele][0]], [whole_xy[ele][1]], color=color[i%len(color)], marker=shape[i%len(shape)])
			axis([min(whole_xy[:,0]), max(whole_xy[:,0]), min(whole_xy[:,1]), max(whole_xy[:,1])])
			xlabel('x-coordinate on CCD plane')
			ylabel('y-coordinate on CCD plane')
			savefig(output_dir + '/trends_position.png')
	
	#Start de-trending.
	print '#Start to de-trend every light curves..'
	#standard deviation list of before and after de-trending
	std_list = []
	for i in range(len(whole_lc)):
		print '	#Now de-trending %dth star.. %s' %((i + 1), file_list[i])
		detrended = get_quadprog(whole_lc[i] / median(whole_lc[i]) - 1. , trend_set)
		detrended += 1.
		detrended *= median(whole_lc[i])
		#detrended = get_quadprog(whole_lc[i] , trend_set)
		print '		#Save the de-trended light curve..'
		fp = open(output_dir + '/' + os.path.basename(file_list[i]), 'w')
		if len(trend_set) == 1:
			detrended = detrended[0]
		for flux in detrended:
			fp.write('%.2f\n' % (round(flux, 2)))
		fp.close()
		if std_save:
			std_list.append([std(whole_lc[i]), std(detrended)])
		
		if img:
			print '		#Save as image file..'
			clf()
			ax = subplot(211)
			title(os.path.basename(file_list[i]))
			plot(whole_lc[i] / median(whole_lc[i]), 'b.')
			text(0.98, 0.9, r'Raw, $\sigma$ : %.1f' % std(whole_lc[i]), horizontalalignment='right', \
			verticalalignment='center', transform=ax.transAxes, fontsize=10, \
			bbox=dict(facecolor='blue', alpha=0.5))
			ylabel('Flux')
			ax = subplot(212)
			plot(detrended / median(whole_lc[i]), 'g.')
			text(0.98, 0.9, r'De-trended, $\sigma$ : %.1f' % std(detrended), horizontalalignment='right', \
			verticalalignment='center', transform=ax.transAxes, fontsize=10, \
			bbox=dict(facecolor='green', alpha=0.5))
			if flux_based:
				ylabel('Normalized Flux')
			else:
				ylabel('Normalized Magnitude')
			xlabel('Time indices')
			savefig(output_dir + '/' + os.path.basename(file_list[i]) + '.png')
	
	if std_save:
		print '#Save standard deviation values..'
		fp = open(output_dir + '/std_list', 'w')
		for i in range(len(std_list)):
			fp.write('%.2f %.2f %s\n' % (round(std_list[i][0], 2), round(std_list[i][1], 2), file_list[i]))
		fp.close()
		
if __name__=='__main__':
	
	file_l = '' # file list of light-curves.
	correlation_matrix = '' # correlation matrix file.
	xy_coords = '' # x and y coordinates of all stars.
	output_d = '' # output directory.
	column_f = 0 # colun number of flux in light-curve files.
	flux_b = True # flux based or magnitude based.
	zero_m = 20 # zero magnitude when light-curves are in magnitude based.
	show_t = False # show dendrogram of hierarchical tree.
	ini_s = 2 # initial seeds constraint.
	bin_w = 0 # binning window size.
	min_t = 5 # minimum number of stars in clusters. 5~10 is good choise.
	n_iter = 1
	l_sig = 0.1
	
	if len(sys.argv) == 1:
		print help
		sys.exit()
	
	#read command line option.
	try:
		optlist, args = getopt.getopt(sys.argv[1:], 'h:f:c:d:b:z:s:o:l:b:t:i:m:g:')
	except getopt.GetoptError, err:
		print help
		sys.exit()
	
	for o, a in optlist:
		if o in ('-h'):
			print help
			sys.exit()
		elif o in ('-f'):
			file_l = a
		elif o in ('-c'):
			correlation_matrix = a
		elif o in ('-d'):
			xy_coords = a
		elif o in ('-b'):
			if a == 'False':
				flux_b = False
			else:
				flux_b = True
		elif o in ('-z'):
			zero_m = float(a)
		elif o in ('-s'):
			show_t = a
		elif o in ('-o'):
			output_d = a
		elif o in ('-l'):
			column_f = int(a)
		elif o in ('-n'):
			bin_w = int(a)
		elif o in ('-t'):
			ini_s = float(a)
		elif o in ('-i'):
			n_iter = int(a)
		elif o in ('-m'):
			min_t = int(a)
		elif o in ('-g'):
			l_sig = float(a)
		else:
			continue
		
	file_list = glob.glob(file_l)
	if len(file_list) == 0:
		print 'There is no file in : ' + file_l
		sys.exit()
	if len(file_list) < 50:
		print 'Too few light curves : %d.\nIf you still want to run the program, modify line#406 in pdtrending.py.' % (len(file_list))
		sys.exit()
	file_list.sort()
	
	pdetrending(file_list, correlation_matrix, xy_coords, show_tree=show_t, \
		min_template_stars=min_t, flux_based=flux_b, zero_mag=zero_m, output_dir=output_d, \
		column_flux=column_f, bin_window=bin_w, initial_seed=ini_s, n_iteration=n_iter, l_significance=l_sig)
	