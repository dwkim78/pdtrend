#!/usr/bin/env python

'''
	Small function set for various usages.
'''

from numpy import *

def sigma_clipping(x, sigma=3., iteration=1):
	'''
	Sigma clipping of x array. Replaced the values over the
	sigma limit with the median values of x array.
	
	x :
		Array of elements.
	
	sigma :
		Threshold.
		
	iteration :
		Number of iteration of sigma clipping. 1~3 is desirable.
	
	return :
		Sigma clipped array of x
	'''
	
	x = array(x)
	xx = x.copy()
	for i in range(iteration):
		median_val = median(x)
		std_val = std(x)
		xx[where(fabs(x - median_val) > std_val * sigma)] = median_val
	return xx

def cmp_length(l, m):
	'''
		Compare the length of two input elements.
		For sort function in length.
	'''
	return len(l) - len(m)

def bin_lc(x, window_size=10, axis=None):
	'''Binning array and return new shortened array
	
	x : 
		original series.
	
	window_size = 10 : 
		default window size.
	
	axis = None : 
		index of x
		
	return :
		list of [bin_lc, std_lc, new_axis if axis != None]
	'''
	x = array(x)
	new_lc = []
	std_lc = []
	for i in range(int(len(x) / window_size) + 1):
		start_index = i * window_size
		end_index = (i + 1) * window_size
		if end_index >= len(x):
			break
		new_lc.append(mean(x[start_index:end_index]))
		std_lc.append(std(x[start_index:end_index]))
	
	if axis != None:
		new_axis = []
		for i in range(int(len(axis) / window_size) + 1):
			start_index = i * window_size
			end_index = (i + 1) * window_size
			if end_index >= len(axis):
				break
			new_axis.append(mean(axis[start_index:end_index]))
	
	if axis != None:
		return [array(new_lc), array(std_lc), array(new_axis)]
	else:
		return [array(new_lc), array(std_lc)]
		
if __name__=='__main__':
	print 'Function set for various usages'
