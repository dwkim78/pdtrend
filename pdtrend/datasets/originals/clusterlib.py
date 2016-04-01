#!/usr/bin/env python

'''
Libraries for PDtrending algorithm. The major part of the libraries
consist of clustering routine.
'''

import sys
import math
from numpy import *

import rpy
from rpy import r

from misclib import *

sys.path.append('./')

def Pycluster_to_hcluster(tree):
	'''
	Transform tree structure of Pycluster to hcluster to use 
	dendrogram function of hcluster which plots dendrogram.
	
	tree :
		Tree from Pycluster.
	
	return :
		hcluster tree.
	'''
	
	hcluster_tree = []
	len_tree = len(tree) + 1
	i = 0
	for node in tree:
		hcluster_node = []
		#make new node for hcluster
		node_left = node.left
		node_right = node.right
		#modify index of Pycluster for hcluster.
		if node_left < 0:
			node_left = (node_left + 1) * -1 + len_tree
		if node_right < 0:
			node_right = (node_right + 1) * -1 + len_tree
		
		#sort
		if node_left > node_right:
			dummy = node_right
			node_right = node_left
			node_left = dummy
		
		#count number of whole elements below each node.
		s_group = []
		simplify_group(find_group_with_node_index(tree, i), s_group)
			
		hcluster_node.append([node_left, node_right, node.distance, len(s_group)])
		hcluster_tree.append(hcluster_node[0])
		i += 1
	hcluster_tree = array(hcluster_tree)
	#print hcluster_tree
	
	return hcluster_tree

def remove_small(groups, size=9):
	'''
	Remove groups smaller than size.
	
	groups :
		List of group.
	
	size :
		Minimum size of group.
	'''
	
	new_groups = []
	
	for group in groups:
		if len(group) > size:
			new_groups.append(group)
			
	return new_groups

def find_group_with_distance(tree, distance):
	'''
	Find group from tree based on the index.
	
	tree :
		Tree structure returned from Pycluster module
		
	distance :
		Maximum distance within group.
		
	return :
		Multiple groups within the distance.
	'''
	
	groups = []
	
	#get node within the distance
	node_list = []
	for i in range(len(tree)):
		if tree[i].distance > distance:
			break
		node_list.append(i)
	
	#make group.
	temp_groups = []
	for node_index in node_list:
		group = find_group_with_node_index(tree, node_index)
		left_group = []
		right_group = []
		simplify_group(group[0], left_group)
		simplify_group(group[1], right_group)
		group = left_group + right_group
		
		temp_groups.append(group)
		
	#check overlapped group and remove.
	groups = remove_subset_cluster(temp_groups)
	
	return groups

def find_group_with_node_index(tree, index):
	'''
	Find group from tree based on the index of node.
	
	tree :
		Tree structure returned from Pycluster module
		
	index :
		Real index of tree.
		
	return :
		Two groups divided by index.
	'''
	
	group = []
	node = tree[index]
	eles = [node.left, node.right]
	for ele in eles:
		if ele >= 0:
			group.append([ele])
		else:
			group.append(find_group_with_link_index(tree, ele))
			
	return group

def find_group_with_link_index(tree, index):
	'''
	Find group from tree based on the index of link.
	Value of link is smaller than -1. -1 means the element
	is linked with 0th node. -2 means it is linked with 1th node and so on.
	
	tree :
		Tree structure returned from Pycluster module
		
	index :
		Linkage index less than -1.
		
	return :
		Two groups divided by index.
	'''
	
	group = []
	node = tree[index * -1 - 1]
	eles = [node.left, node.right]
	for ele in eles:
		if ele >= 0:
			group.append([ele])
		else:
			group.append(find_group_with_link_index(tree, ele))
			
	return group
	
	
def simplify_group(group, s_group):
	'''
	Unify multiple lists inside group into one list.
	
	group :
		Looks like [[6, 4], 3]
	
	s_group :
		simplified group looks like [6, 4, 3]
		
	return :
		None
	'''
	
	if isinstance(group, int):
		s_group.append(group)
		return
	else:
		for ele in group:
			simplify_group(ele, s_group)

def find_group_DW(tree, dist_matrix, mes=0, seed_max=0, l_significance=0.1):
	'''
	For more details, see the Kim et al. 2008.
		
	tree :
		Tree structure returned from Pycluster module.
	
	dist_matrix :
		Distance matrix (= 1. - correlation matrix)
	
	mes = 0 :
		Total number of measurement of each light curve.
	
	seed_max = 0 :
		To get more tighter seed. 1 ~ 10 are good values. '10' gets more tighter clusters than '1'.
		
	return :
		List of clusters.
	'''
	
	r.library('nortest')
	
	clusters = []
	#print tree, len(tree)
	density_list = []
	for i in range(len(dist_matrix) - 1):
		for j in range(i + 1, len(dist_matrix)):
			density_list.append(dist_matrix[i][j])
	density_list_clip = sigma_clipping(density_list, sigma=3.)
	overall_density = (max(density_list_clip) - min(density_list_clip)) / len(dist_matrix)
	#print overall_density, mean(density_list_clip), std(density_list_clip)
	
	#get highly correlated pair of elements.
	initial_seed = []
	for i in range(len(tree)):
		#both left and right element has to be star. not a link to other cluster.
		if tree[i].left >= 0 and tree[i].right >= 0:
			#to get more tight elements.
			if dist_matrix[tree[i].left][tree[i].right] <= median(density_list_clip) / seed_max:
				if mes == 0:
					initial_seed.append(i)
				elif dist_matrix[tree[i].left][tree[i].right] <= (1. - 3. / math.sqrt(mes)):
					initial_seed.append(i)
	#print initial_seed
	
	#start from highly correlated initial pair.
	for i in initial_seed:
		#print tree[i]
		current_node = i
		while current_node < len(tree) - 1:
			cluster_1 = []
			cluster_2 = []
			#find base cluster --> cluster_1
			simplify_group(find_group_with_node_index(tree, current_node), cluster_1)
			#find cluster which will be merged --> cluster_2
			dummy = find_one_side_group(tree, (current_node + 1) * -1)
			current_node = dummy[0]
			simplify_group(dummy[1], cluster_2)
			
			#check the density changes with overall density
			#initial density
			d_1 = []
			for ele_i in range(len(cluster_1) - 1):
				for ele_j in range(ele_i + 1, len(cluster_1)):
					if ele_i != ele_j:
						d_1.append(dist_matrix[cluster_1[ele_i]][cluster_1[ele_j]])
			#density after merged
			d_merge = []
			cluster_3 = hstack([cluster_1, cluster_2])
			for ele_i in range(len(cluster_3) - 1):
				for ele_j in range(ele_i + 1, len(cluster_3)):
					if ele_i != ele_j:
						d_merge.append(dist_matrix[cluster_3[ele_i]][cluster_3[ele_j]])
			
			d_1 = array(d_1)
			d_merge = array(d_merge)
			if len(d_merge) < 8:
				continue
			else:
				#the resulting clusters are almost identical. not use anymore.
				#d_merge = array(d_merge)
				#d_merge = .5 * log((1. + d_merge) / (1. - d_merge))
				
				ad = r.ad_test(d_merge)
				ad_p = ad['p.value']
				p_value = ad_p
				
				#check the level of significance
				#if it's out of normality, the previous cluster is the final cluster.
				if p_value < l_significance:
					#becausd AD test needs at least 8 elements.
					if len(cluster_1) >= 5:
						#print cluster_1
						clusters.append(cluster_1)
					break
				#it's still gaussian, but if there comes outliers into clusters, stop it.
				#the resulting clusters are almost identical. not use anymore.
				#elif len(d_1[where(d_1 > mean(density_list_clip))]) > 0:
				#	if len(cluster_1) >= 5:
				#		clusters.append(cluster_1)
				#	break
	
	return clusters

def find_one_side_group(tree, index):
	'''
		Return group of one side.
		
		tree :
			Tree structure returned from Pycluster module
			
		index :
			Index of tree structure. This routine first find node of index 
			which looks like [index, other_index or node_value]
			and, return the group of [other_index] or [value]
	'''
	
	for i in range(len(tree)):
		if tree[i].left == index:
			group_index = tree[i].right
			break
		if  tree[i].right == index:
			group_index = tree[i].left
			break
	
	if group_index >=0 :
		return [i, [group_index]]
	else:
		return [i, find_group_with_link_index(tree, group_index)]

def cal_correlation_matrix(whole_lc, bin_size=0):
	'''
	Calculate the Pearson correlation matrix.
	
	whole_lc : 
		List of light curves.
	
	bin_size :
		Binning window size. to reduce noise.
		Binning the data only when calculated correlation values.
		Will not change original light curves.
	
	return :
		The Pearson correlation magtrix.
	'''
	
	corr_list = ones((len(whole_lc), len(whole_lc)))
	for i in range(len(whole_lc) - 1):
		#print '	#Correlation values of %dth star is calculating..' % (i + 1)
		for j in range(i, len(whole_lc)):
			if bin_size == 0:
				pear_corr = corrcoef(whole_lc[i], whole_lc[j])[0][1]
			else:
				pear_corr = corrcoef(bin_lc(whole_lc[i], bin_size)[0], bin_lc(whole_lc[j], bin_size)[0])[0][1]
			corr_list[i, j] = pear_corr
			corr_list[j, i] = pear_corr
	
	return corr_list
	
def get_detrened_lc_set(lc, trend_set):
	'''
	Return de-trended lc by using trends set.
	We use Multiple Linear Regression Method to de-trend.
	
	lc :
		Original light curve of flux.
	
	trend_set :
		Set of trend light curves constructed by create_trend routine.
	
	return :
		De-trended light curve.
	'''
	
	#Multiple linear regression method, least square method
	X = ones([len(trend_set[0]), len(trend_set) + 1])
	X[:, 0] = [1]
	X[:, 1:] = transpose(trend_set)
	beta_1 = linalg.inv(dot(transpose(X), X))
	beta_2 = dot(transpose(X), lc)
	beta = dot(beta_1, beta_2)
	
	return lc - dot(X, beta)

def get_quadprog(lc, trend_set):
	'''
	Return de-trended lc by quadratic programming.
	It constraints the free parameters to be bigger than 0.
	See Kim et al. 2008 for more details.
	
	lc :
		Original light curve of flux.
	
	trend_set :
		Set of trend light curves constructed by create_trend routine.
	
	return :
		De-trended light curve.
	'''
	
	r.library('quadprog')
	
	X = transpose(trend_set)
	dmat = r.crossprod(X, X)
	dvec = r.crossprod(lc, X)
	
	results = r.solve_QP(dmat, dvec, r.diag(len(trend_set)))
	#print results['solution'], results['value']
	
	return lc - dot(results['solution'], trend_set)

def get_linprog(lc, trend_set):
	'''
	Return de-trended lc by linear programming.
	It constraints the free parameters to be bigger than 0.
	
	lc :
		Original light curve of flux.
	
	trend_set :
		Set of trend light curves constructed by create_trend routine.
	
	return :
		De-trended light curve.
	'''
	
	r.library('linprog')
	
	X = transpose(trend_set)
	#dmat = r.crossprod(X, X)
	dvec = r.crossprod(lc, X)
	
	results = r.solveLP(dvec, zeros([len(trend_set)]), r.diag(len(trend_set)))
	print results['opt'], results['solution']
	sys.exit()
	
	#return lc - dot(results['solution'], trend_set)
	
def kovacs(lc, trend_set):
	'''
	Apply TFA with determined trends by our algorithm.
	See Kovacs et al. 2005 for more detail.
	
	lc :
		Original light curve of flux
	
	trend_set :
		Set of trend light curves constructed by create_trend routine.
	
	return :
		De-trended light curve.
	'''
	
	dup_lc = lc[::]
	dup_trend = trend_set[::]
	
	dup_lc -= mean(dup_lc)
	for i in range(len(dup_trend)):
		dup_trend[i] -= mean(dup_trend[i])
	
	g = zeros((len(dup_trend), len(dup_trend)))
	h = zeros((len(dup_trend)))
	c = zeros((len(dup_trend)))
	for i in range(len(dup_trend)):
		for j in range(len(dup_trend)):
			g[i, j] = sum(dot(dup_trend[i], dup_trend[j]))
	G = linalg.inv(g)
	
	for i in range(len(dup_trend)):
		h[i] = sum(dot(dup_lc, dup_trend[i]))
	for i in range(len(dup_trend)):
		c[i] = sum(dot(G[i:i + 1], h))

	trend = dot(c, dup_trend)
	detrended = dup_lc - trend
	
	return detrended
	
def create_trend(template_set, whole_lc):
	'''
	Create trend with selected template light curves.
	We use weighted sum of normlized light curves by median values of the lc.
	See Kim et al. 2008 for more details.
	
	template_set :
		Indices of template light curves.
	
	whole_lc : 
		Light curves of every stars.
		
	return :
		Constructed trend.
	'''
	
	trend = []
	trend_lc = []
	
	#normalized by mean value
	normalized = []
	for i in template_set:
		normalized.append((whole_lc[i] - mean(whole_lc[i])) / std(whole_lc[i]))
	
	#weighting by sigma^2
	std_inv_list = []
	for i in range(len(normalized)):
		std_inv_list.append(1. / std(normalized[i])**2.)
	weight_list = std_inv_list / sum(std_inv_list)
	
	for i in range(len(normalized)):
		trend_lc.append(normalized[i] * weight_list[i])
	trend_lc = array(trend_lc)
	
	#make trend..
	for i in range(len(trend_lc[0])):
		#trend.append(median(trend_lc[:, i]))
		trend.append(sum(trend_lc[:, i]))
	
	return array(trend)

def remove_subset_cluster(cluster_list):
	'''
	Remove subset clusters which are included in other clusters.
	
	cluster_list : 
		Initial list of clusters.
		
	return :
		List of clusters after removal of subsets.
	'''
	
	cluster_list.sort(cmp=cmp_length)
	dup = cluster_list[::]
	
	new_cluster_list = []
	#print dup
	
	#check sub-set.
	for i in range(len(dup)):
		if_duplicated = 0
		cluster_d = dup[i]
		len_cluster_d = len(cluster_d)
		#print 'd:', cluster_d
		for j in range(len(cluster_list)):
			if i == j:
				break
			cluster_o = dup[j]
			len_cluster_o = len(cluster_o)
			#print 'o:', cluster_o
			if len_cluster_d <= len_cluster_o:
				duplicated = [ele for ele in cluster_d if ele in cluster_o]
				if len(duplicated) == len_cluster_d:
					if_duplicated = 1
					break
		
		if if_duplicated == 0:
			new_cluster_list.append(cluster_d)
	
	return new_cluster_list

def remove_duplicated_cluster(corr, cluster_list):
	'''
	Make clusters simple by removing duplicated members through all clusters.
	
	corr :
		The Pearson correlation matrix.
	
	cluster_list :
		Initial cluster list.
	
	return :
		list of cluster after removal of dup. clusters.
	'''
	
	dup=corr.copy()
	
	#make a member of cluster list.
	member_list=[]
	for cluster in cluster_list:
		for member in cluster:
			if not member in member_list:
				member_list.append(member)
	
	#check if each member are overlapped through multiple cluster.
	for member in member_list:
		mem_corr_list={}
		cluster_index=0
		for cluster in cluster_list:
			if len(cluster)==2:
				cluster_index+=1
				continue
			#calculate average correlation if overlapped.
			sum_corr=0.
			avg_corr=0.
			if member in cluster:
				for cluster_member in cluster:
					if member!=cluster_member:
						sum_corr+=dup[member, cluster_member]
				avg_corr=sum_corr/(len(cluster)-1)
				mem_corr_list[cluster_index]=avg_corr
			cluster_index+=1
		#print member, mem_corr_list, max(mem_corr_list)
		
		#remain only one highest correlated value, remove else others.
		if len(mem_corr_list)>1:
			max_key=max(mem_corr_list)
			for key, val in mem_corr_list.iteritems():
				#print key, cluster_list[key]
				if key!=max_key:
					cluster_list[key].pop([j for j in range(len(cluster_list[key])) if cluster_list[key][j]==member][0])
	
	return cluster_list
	
if __name__=='__main__':
	print 'Function set for PDtrending algorithm'
