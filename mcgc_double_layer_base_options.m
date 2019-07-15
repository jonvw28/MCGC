function [LiDAR] = mcgc_double_layer_base_options(save_inter,...
	data_folder, out_folder, raw_file, chm_file, pr_allom, sigxy_1, sigz_1,...
	allom, hgrad_w_1, zgrad_w_1, db_eps_1, db_pts_1, min_pts_1, subSamp,...
	layers_same, sigxy_2, sigz_2, hgrad_w_2, zgrad_w_2, db_eps_2, db_pts_2,...
	min_pts_2)
%MCGC_DOUBLE_LAYER_BASE_OPTIONS wrapper for double layer MCGC to give access to
%only key arguments
%  
% Apply double layer MCGC or subsampled MCGC without needing to supply many
% of the arguments that are presented in the base code
%
%	PLEASE NOTE: For this to work, it is necessary to have Rscript
%				added to the system search path. You will also need to have
%				installed the itcSegment package in R
% Syntax
%
%	[LiDAR] = mcgc_base_options(data_folder, out_folder,...
%		raw_file, chm_file, pr_allom, sigxy, sigz, allom, hgrad_w, zgrad_w,...
%		db_eps, db_pts,	min_pts)
%
%	[LiDAR,sample_idx] = mcgc_base_options(data_folder, out_folder,...
%		raw_file, chm_file, pr_allom, sigxy, sigz, allom, hgrad_w, zgrad_w,...
%		db_eps, db_pts,	min_pts, subSamp)
%
%
% Method
%
%	This applies a double pass of MCGC, using either the full data or subsampled
%	version. This wrapper sets a lot of the arguments to these functions,
%	leaving only a few key arguments open to the user. Those arguments which are
%	automatically set match the settings used with MCGC in literature.
%
%	Full details of the methods for MCGC double layer can be found in the
%	documentation for mcgc_double_layer. For each pass set the values of
%	parameters for *_1 or *_2 respectively. If you want to use the same
%	setting for both passes then set layers_same to 1, otherwise set it to 0.
%
%	save_inter lets you decide if the results of each pass of MCGC should be
%	saved or not.
%
%
% Inputs
%
%		save_inter:		Set to 1 to save results of each pass of MCGC as well as
%						final result, set to 0 to only save final result
%
%		data_folder:	path from the present directory to the folder containing
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
%		pr_allom:		Allometric look	up table of the co-efficients for 
%						tau = 0.01 to 1	for percentile regression of 
%						R = exp(a) * H^b. Here set the argument to the name of 
%						this file with reference to	the working directory. 
%						(eg '.\data\lut.csv'). The window will then scale with
%						allometry
%
%		sigxy_1/2:		Parameter for significance of planimetric distance in 
%						linkages
%
%		sigz_1/2:			Parameter for significance of vertical distance in 
%						linkages
%
%		allom:			Allometric lookup table for centroid computation. Must 
%						have first 2 columns as height (rounded to nearest
%						metre) and allometric radius respectively.
%
%		hgrad_w_1/2:	Parameter to set significance of Delta_H term
%
%		zgrad_w_1/2:	Parameter to set significance of Delta_Z term
%
%		db_eps_1/2:		The parameter epsilon to be used in DBScan clustering.
%						If set to zero this switches the DBScan step off
%
%		db_pts_1/2:		The number of points that must be in the epsilon 
%						neighbourhood in the DBScan method
%
% 		min_pts_1/2:	Minimum number of points needed in each cluster for it
%						not to be rejected
%
%		subSamp:		Sets factor by which to subsample data. A value of 5 
%						means 1/5 of the data are used. Set it to 0 to avoid
%						subsampling
%
%		layers_same:	Set to 1 to use the same parameters for both passes (in
%						which case you don't need to specify the *_2 arguments)
%						Set to 0 to specify using the values of arguments you
%						set for the second pass
%
%
% Outputs
%
%	LiDAR:			The data as supplied with an extra column at 
%					the end with a cluster label (numbered) for each point
%
%			The prior that is generated will also be saved in the output folder
%			with the name 'chm_file'_prior.csv as will the subsampled las files
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
%		This function requires mcgc which uses mcgc_pipeline requiring mcgc_cut,
%		which also requires nystrom_ext which in turn requires
%		compute_uncon_weights_nystrom. mcgc also requires select_good_tree.
%		If using subsampling then mcgc_subsample is also needed which depends on
%		mcgc and is subsequent tree.
%
%		mcgc_double_layer_base_options
%		 ->
%			mcgc_double_layer
%			 ->
%				(mcgc_subsample)
%				 ->
%					mcgc
%					 ->
%						mcgc_pipeline
%						 ->
%							mcgc_cut
%							 ->
%								nystrom_ext
%						 		 ->
%									compute_uncon_weights_nystrom
%					 ->
%						select_good_tree
%
%
%
%		Jonathan Williams
%		jonvw28@gmail.com			         
%		11/01/2019		         
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input Checking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Type checking
if(isnumeric(layers_same)~=1 || size(layers_same,1) ~= 1 ||...
				size(layers_same,2) ~= 1)
	error('layers_same must be a single number')
elseif(isnumeric(subSamp)~=1 || size(subSamp,1) ~= 1 || size(subSamp,2) ~= 1)
	error('subSamp must be a single number')
end

% Sensible values?
if(layers_same~=0 && layers_same~=1)
	error('layers_same must be 1 (copy arguments) or 0 (specify both passes)')
elseif(subSamp~=0 && subSamp < 1)
	error('subSamp must be >= 1. Set to 0 to avoid subsampling')
end


% warn if arguments for second pass are not used
if(layers_same==1 && nargin>16) 
	warning(['Cloning first pass arguments for pass 2, some arguments set '...
			'in function call will be ignored'])
end

% check all arguments are specified for second pass
if(layers_same==0 && nargin<23) 
	error('Not enough arguments specified for second pass')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Function Body %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set parameters which aren't open
pr_grid_size = 0.5;
pr_ht_thr = 5;
pr_tau = '50';
gc_ht_thr = 2;
nystrom_fac = 5;
max_grp_ratio = 2;
r_rat = 0.5;
h_grad_opt = 4;
z_grad_opt = 4;
merge_tol = 9e-3;
test_allom = allom;
top_frac = 0.98;
alo_share_max = 0.6;
z_olap = 0.25;
min_alo = 0.95;
alo_opt = 1;
imp_lim = 2;
sub_db = 0;
out_prec = '%.2f';

% Get mfile location
scrpt_loc = fileparts(which(mfilename));
call_loc = pwd;

% Load dependencies
addpath([scrpt_loc '\scripts\MATLAB\graph_cut']);

% Sort out second pass arguments
if(layers_same==1)
	sigxy_2 = sigxy_1;
	sigz_2 = sigz_1;
	hgrad_w_2 = hgrad_w_1;
	zgrad_w_2 = zgrad_w_1;
	db_eps_2 = db_eps_1;
	db_pts_2 = db_pts_1;
	min_pts_2 = min_pts_1;
end

% Run MCGC
if(subSamp == 0) % full data MCGC

LiDAR = mcgc_double_layer(0, save_inter, out_prec,...
	data_folder, out_folder, raw_file, chm_file, pr_grid_size, pr_ht_thr,...
	pr_ht_thr, pr_allom, pr_tau, gc_ht_thr, gc_ht_thr, nystrom_fac,...
	max_grp_ratio, sigxy_1, sigxy_2, sigz_1, sigz_2, allom, r_rat,...
	h_grad_opt, hgrad_w_1, hgrad_w_2, z_grad_opt, zgrad_w_1, zgrad_w_2,...
	merge_tol, test_allom, top_frac, alo_share_max, z_olap, min_alo,...
	db_eps_1, db_eps_2, db_pts_1, db_pts_2, min_pts_1, min_pts_2, alo_opt,...
	subSamp, subSamp, imp_lim, imp_lim, sub_db);

else % use subsampling
LiDAR= mcgc_double_layer(1, save_inter, out_prec,...
	data_folder, out_folder, raw_file, chm_file, pr_grid_size, pr_ht_thr,...
	pr_ht_thr, pr_allom, pr_tau, gc_ht_thr, gc_ht_thr, nystrom_fac,...
	max_grp_ratio, sigxy_1, sigxy_2, sigz_1, sigz_2, allom, r_rat,...
	h_grad_opt, hgrad_w_1, hgrad_w_2, z_grad_opt, zgrad_w_1, zgrad_w_2,...
	merge_tol, test_allom, top_frac, alo_share_max, z_olap, min_alo,...
	db_eps_1, db_eps_2, db_pts_1, db_pts_2, min_pts_1, min_pts_2, alo_opt,...
	subSamp, subSamp, imp_lim, imp_lim, sub_db);

end

end