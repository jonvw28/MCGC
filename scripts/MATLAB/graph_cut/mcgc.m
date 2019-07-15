function [finalSeg] = mcgc(data_folder, out_folder, raw_file,...
	chm_file, pr_grid_size, pr_ht_thr, pr_allom, pr_tau, gc_ht_thr,...
	nystrom_fac, max_grp_ratio, sigxy, sigz, allom, r_rat, h_grad_opt,...
	hgrad_w, z_grad_opt, zgrad_w, merge_tol, test_allom, top_frac,...
	alo_share_max, z_olap, min_alo, db_eps, db_pts, min_pts, alo_opt)
%MCGC wrapper for complete MCGC with allometric checking algorithm
%  
% Attempt to apply MCGC graph cut segmentation to LAS data. If this works then
% apply allometric checking as well. All parameters are available to control.
%
%	PLEASE NOTE: For this to work, it is necessary to have Rscript
%				added to the system search path. You will also need to have
%				installed the itcSegment package in R
%
% Syntax
%
%	[finalSeg] = mcgc(data_folder, out_folder, raw_file,...
%		chm_file, pr_grid_size, pr_ht_thr, pr_allom, pr_tau, gc_ht_thr,...
%		nystrom_fac, max_grp_ratio, sigxy, sigz, allom, r_rat, h_grad_opt,...
%		hgrad_w, z_grad_opt, zgrad_w, merge_tol, test_allom, top_frac,...
%		alo_share_max, z_olap, min_alo, db_eps, db_pts, min_pts, alo_opt)
%
%
% Method
%
%	This applies the MCGC normalised cut (using Nystrom extension) [1,2] 
%	followed by allometric testing. For full details of all steps of this
%	process see the documentation of dependancies, or the literature associated 
%	with MCGC.
%
%	Should either step fail, a file will automatically be generated with details
%	of the error. This will end with '_seg_err.txt' if the normalised cut step
%	fails or with '_alo_check_err.txt' if the allometric checking step fails.
%
%	This function gives access to all available arguments across all steps of 
%	the algorithm. For simpler forms of the pipeline see the wrappers which
%	preset many arguments.
%
%
% Inputs
%
%		data_folder		path from the present directory to the folder containing
%						the las data file
%
%		output_folder:	path from the present directory to the folder where the
%						prior shoudl be saved
%
%		raw_file:		A string containing the name of the las file to be
%						segmented - with raw z values (without the .las
%						extension)
%
%		chm_file:		A string containing the name of the las file to be
%						segmented - with chm z values (without the .las 
%						extension)
%
%		pr_grid_size:	Size of grid used to generate CHM and prior (in metres)
%
%		pr_ht_thr:		Minimum height of a prior tree-top for it to be used
%
%		pr_allom:		Allometric look	up table of the co-efficients for 
%						tau = 0.01 to 1	for percentile regression of 
%						R = exp(a) * H^b. Here set the argument to the name of 
%						this file with reference to	the working directory. 
%						(eg '.\data\lut.csv'). The window will then scale with
%						allometry
%
%		pr_tau			Value of tau used for computing	the edge buffer for
%						generating the prior supply as a string, eg '90' for 0.9
%
%		gc_ht_thr		Cut-off height below which points are ignored and not 
%						used in graph cut - to avoid ground returns and very 
%						small vegetation
%
%		nystrom_fac:	This sets the number of points to include in the
%						subsample given by nystrom_fac*(max number of clusters)
%
%		max_grp_ratio:	Factor applied to number of priors to generate maximum
%						number of clusters sought
%
%		sigxy: 			Parameter for significance of planimetric distance in 
%						linkages
%
%		sigz: 			Parameter for significance of vertical distance in 
%						linkages
%
%		allom:			Allometric lookup table for centroid computation. Must 
%						have first 2 columns as height (rounded to nearest
%						metre) and allometric radius respectively.
%
%		r_rat:			Fraction of lookup radius from allom to use for centroid
%						computation (recommended to use 0.5 or 1)
%
%		h_grad_opt:		Sets which option to use for Delta_H comparison:
%							1: 'Uniform adjustment'
%							2: 'Inverse separation weighting'
%							3: 'Weight by Delta difference'
%							4: 'Composite weighting'
%
%		hgrad_w:		Parameter to set significance of Delta_H term
%
%		z_grad_opt:		Sets which option to use for Delta_Z comparison:
%							1: 'Uniform adjustment'
%							2: 'Inverse separation weighting'
%							3: 'Weight by Delta difference'
%							4: 'Composite weighting'
%
%		zgrad_w:		Parameter to set significance of Delta_Z term
%
%		merge_tol:		Should the chm and raw height las files not match up
%						when read in by lasread, this sets the absolute
%						tolerance allowed when matching the x,y coordinates
%						between the two data sources.
%
%		test_allom:		An N-by-2 numeric lookup table for allometry. Used for
%						allometric feasibility checking in processing. The first 
%						column should be tree heights in increments of metres,
%						starting from 0 (ie, 0,1,2,..,N-1). The second column
%						should then be your relevant radius for a tree of this 
%						height (eg 95-th percentile of all trees in your
%						dataset)
%
%		top_frac:		Sets which points are used to compute the 'centre' of a
%						cluster. All points above this proportion of the top 
%						point's height are included (eg 0.98 to include only
%						points at least 98% as high as the highest point)
%
%		alo_share_max:	Maximum proportion of a smaller clusters points that may
%						lie within the allometric radius of a taller tree. If 
%						proportion exceeds this then the smaller tree will be a 
%						candidate for merging. (eg 0.6 if this cutoff should be
%						60% of points)
%
% 		z_olap:			When merging clusters this sets the height overlap 
%						criterion. Here (100*z_olap)-th percentile of heights of
%						the taller tree is compared to (100-100*z_olap)-th 
%						percentile of heights for the smaller tree. The merge
%						occurs if the latter is higher. Eg, if z_olap is 0.1
%						then the merge will be enacted if the the top 10% of 
%						points of the smaller tree are above the bottom 10% of
%						points for the taller tree
%
% 		min_alo:		The minimum proportion of points that must lie within
%						the allometry radius of a cluster for this cluster
%						to be considered valid (eg 0.8 for 80% of points)
%
%		db_eps:			The parameter epsilon to be used in DBScan clustering.
%						If set to zero this switches the DBScan step off
%
%		db_pts:			The number of points that must be in the epsilon 
%						neighbourhood in the DBScan method
%
% 		min_pts:		Minimum number of points needed in each cluster for it
%						not to be rejected
%
% 
%		alo_opt:		Sets the behaviour for remaining points in allometric 
%						feasibility hierarchical clustering
%							0:	Keep these points as their own cluster to also
%									be tested
%							1:	Mark these points are rejected with a 5
%
%
% Outputs
%
%	finalSeg:		The data as supplied (without points below the
%					minimum height threshold removed) with an extra column at 
%					the end with a cluster label (numbered) for each point
%
%					The first two columns will be the x and y co-ordinates,
%					followed by the raw and above ground heights then intensity
%					and original classification in columns 3-6.  Column 7
%					will then contain a numeric label for the crown of each 
%					point which is assigned a label, or 0 if it is not in above
%					crown. Column 8 will contain the label to show what the
%					allometric processing result was with the following
%					options:
%						0: There is an error in the code (please contact me)
%						1: There is no issue with the cluster
%						2: There are too few points in the cluster
%						5: The point for this row was rejected after allometric
%							testing (only in case of allometry rejection)
%						6: The point in this row was thrown away in a pass
%							of DBScan (if enabled)
%
%			The prior that is generated will also be saved in the output folder
%			with the name 'chm_file'_prior.csv
%
%
% References:
%
%       [1]	A Tutorial on Spectral Clustering, U von Luxburg, 
%       	Statistics and Computing, 17 (4), 2007
%
%		[2] Spectral Grouping Using the Nystrom Method,C Fowlkes et al.,
%			IEEE TRANSACTIONS ON PATTERN ANALYSIS AND MACHINE INTELLIGENCE,
%			VOL. 26, NO. 2, FEBRUARY 2004
%
%		[3]	Graph Clustering, Variational Image Segmentation Methods and Hough
%			Transform Scale Detection for Object Measurement in Images, L_ord
%			Calatroni et al., Journal of Mathematical Imaging and Vision, 2017
%
%
% Dependancy Tree
%
%		This function requires mcgc_pipeline which requires mcgc_cut, which also
%		requires nystrom_ext which in turn requires
%		compute_uncon_weights_nystrom. mcgc also requires select_good_tree
%
%		mcgc
%		 ->
%			mcgc_pipeline
%			 ->
%				mcgc_cut
%				 ->
%					nystrom_ext
%					 ->
%						compute_uncon_weights_nystrom
%		 ->
%			select_good_tree
%
%
%
%		Jonathan Williams
%		jonvw28@gmail.com			         
%		10/01/2019		         
%
%		Dependancies heavily influenced by the example code in [2] and adapted
%		from the code used originally in [3].
%
%		This work makes use of the lasread MATLAB function oriinally written by
%		Thomas J. Pingel. See the address below for its source
%		http://www.tpingel.org/code/lasread/lasread.html
%
% 		This work uses the DBSCAN MATLAB software as available from yarpiz.com.
% 		See the directory with DBSCAN to see the license for use of this
%		software



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Function Body %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get mfile location
scrpt_loc = fileparts(which(mfilename));
call_loc = pwd;

% Load dependencies
addpath(genpath([scrpt_loc '\segment']));
addpath(genpath([scrpt_loc '\post_process']));
addpath(genpath([scrpt_loc '\gc_utils']));

% Run the graph cut segmentation
try
	Seg = mcgc_pipeline(data_folder, out_folder, raw_file,...
	chm_file, pr_grid_size, pr_ht_thr, pr_allom, pr_tau, gc_ht_thr,...
	nystrom_fac, max_grp_ratio, sigxy, sigz, allom, r_rat, h_grad_opt,...
	hgrad_w, z_grad_opt, zgrad_w, merge_tol);
catch err
	err_handler(err,[call_loc '\' out_folder '\' raw_file '_seg_err.txt']);
	return
end

% Set up known variables
raw_z_col = 3;
chm_z_col = 4;
class_col = 7;

% Apply allometric checking
try
	finalSeg = select_good_tree(Seg, raw_z_col, chm_z_col,...
						class_col, test_allom, top_frac, alo_share_max,...
						z_olap, min_alo, db_eps, db_pts, min_pts, alo_opt);
catch err
	err_handler(err,[call_loc '\' out_folder '\' raw_file ...
	'_alo_check_err.txt']);
end
end