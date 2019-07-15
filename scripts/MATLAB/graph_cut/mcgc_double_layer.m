function [LiDAR] = mcgc_double_layer(samp_tog, save_inter, out_prec,...
	data_folder, out_folder, raw_file, chm_file, pr_grid_size, pr_ht_thr_1,...
	pr_ht_thr_2, pr_allom, pr_tau, gc_ht_thr_1, gc_ht_thr_2, nystrom_fac,...
	max_grp_ratio, sigxy_1, sigxy_2, sigz_1, sigz_2, allom, r_rat,...
	h_grad_opt, hgrad_w_1, hgrad_w_2, z_grad_opt, zgrad_w_1, zgrad_w_2,...
	merge_tol, test_allom, top_frac, alo_share_max, z_olap, min_alo,...
	db_eps_1, db_eps_2, db_pts_1, db_pts_2, min_pts_1, min_pts_2, alo_opt,...
	subSamp_1, subSamp_2, imp_lim_1, imp_lim_2, sub_db)
%MCGC_DOUBLE_LAYER wrapper for MCGC to allow double layer application
%  
% Apply MCGC twice with options to alter relavant settings between the layers.
% Available with or without subsampling
%
%	PLEASE NOTE: For this to work, it is necessary to have Rscript
%				added to the system search path. You will also need to have
%				installed the itcSegment package in R
%
% Syntax
%
%	[LiDAR] = mcgc_double_layer(samp_tog, save_inter, out_prec,...
%	data_folder, out_folder, raw_file, chm_file, pr_grid_size, pr_ht_thr_1,...
%	pr_ht_thr_2, pr_allom, pr_tau, gc_ht_thr_1, gc_ht_thr_2, nystrom_fac,...
%	max_grp_ratio, sigxy_1, sigxy_2, sigz_1, sigz_2, allom, r_rat,...
%	h_grad_opt, hgrad_w_1, hgrad_w_2, z_grad_opt, zgrad_w_1, zgrad_w_2,...
%	merge_tol, test_allom, top_frac, alo_share_max, z_olap, min_alo,...
%	db_eps_1, db_eps_2, db_pts_1, db_pts_2, min_pts_1, min_pts_2, alo_opt,...
%	subSamp_1, subSamp_2, imp_lim_1, imp_lim_2, sub_db)
%
%	[LiDAR] = mcgc_double_layer(samp_tog, save_inter, out_prec,...
%	data_folder, out_folder, raw_file, chm_file, pr_grid_size, pr_ht_thr_1,...
%	pr_ht_thr_2, pr_allom, pr_tau, gc_ht_thr_1, gc_ht_thr_2, nystrom_fac,...
%	max_grp_ratio, sigxy_1, sigxy_2, sigz_1, sigz_2, allom, r_rat,...
%	h_grad_opt, hgrad_w_1, hgrad_w_2, z_grad_opt, zgrad_w_1, zgrad_w_2,...
%	merge_tol, test_allom, top_frac, alo_share_max, z_olap, min_alo,...
%	db_eps_1, db_eps_2, db_pts_1, db_pts_2, min_pts_1, min_pts_2, alo_opt)
%
% Method
%
%	This applies the MCGC complete pipeline as provided by the mcgc function.
%	A first application is completed, then all remaining points are used to 
%	complete a sceond pass of MCGC. It is possible to toggle many of the options
%	to alter between the passes. These arguments are marked with an * in the 
%	inputs list. The full explanantion of the steps of this algorithm can be
%	found in the documentation for mcgc and its dependencies or in the
%	literature for MCGC.
%
%	This function can be used with or without the subsampling wrapper. To use
%	the subsampling wrapper, set samp_tog to 1, or set it to 0 to avoid 
%	subsampling. The explanation of this process and its arguments can be found
% 	in the documentation for mcgc_subsample.
%
%	The final results will automatically be saved to the output folder. The
% 	argument out_prec sets the precision to save the data as using a REGEX.
%	It is highly recommended that unless you are fmamiliar with REGEX that you
%	set this to '%.2f' which will store data as float values to cm precision.
%
%	Setting save_inter to 1 will mean the results of each pass of MCGC are
%	independently saved to the output folder in addition to saving the final
%	result. If you only want the final result then set this to 0.
%
%	All other parameters are explained in detail in the relevant dependancies of
%	this function.
%
%	The additonal arguments provided by this function are:
%	samp_tog - set to 1 to use subsampling or 0 to use full dataset
%	out_prec - REGEX setting precision to save output files
%	save_inter - set to 1 to save results of each pass of MCGC or 0 to not
%
%
% Inputs
%
%		samp_tog:		Set to 1 to use subsampling wrapper, or 0 to use full
%						dataset
%
%		save_inter:		Set to 1 to save results of each pass of MCGC as well as
%						final result, set to 0 to only save final result
%
%		out_prec:		Sets precision that las files are saved as. RECOMMENDED
%						TO SET THIS TO '%.2f'		
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
%		pr_grid_size:	Size of grid used to generate CHM and prior (in metres)
%
%	*	pr_ht_thr_1/2:	Minimum height of a prior tree-top for it to be used
%
%		pr_allom:		Allometric look	up table of the co-efficients for 
%						tau = 0.01 to 1	for percentile regression of 
%						R = exp(a) * H^b. Here set the argument to the name of 
%						this file with reference to	the working directory. 
%						(eg '.\data\lut.csv'). The window will then scale with
%						allometry
%
%		pr_tau:			Value of tau used for computing	the edge buffer for
%						generating the prior supply as a string, eg '90' for 0.9
%
%	*	gc_ht_thr_1/2:	Cut-off height below which points are ignored and not 
%						used in graph cut - to avoid ground returns and very 
%						small vegetation
%
%		nystrom_fac:	This sets the number of points to include in the
%						subsample given by nystrom_fac*(max number of clusters)
%
%		max_grp_ratio:	Factor applied to number of priors to generate maximum
%						number of clusters sought
%
%	*	sigxy_1/2:		Parameter for significance of planimetric distance in 
%						linkages
%
%	*	sigz_1/2: 		Parameter for significance of vertical distance in 
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
%	*	hgrad_w_1/2:	Parameter to set significance of Delta_H term
%
%		z_grad_opt:		Sets which option to use for Delta_Z comparison:
%							1: 'Uniform adjustment'
%							2: 'Inverse separation weighting'
%							3: 'Weight by Delta difference'
%							4: 'Composite weighting'
%
%	*	zgrad_w_1/2:	Parameter to set significance of Delta_Z term
%
%		merge_tol:		Sets the absolute tolerances used when comparing
%						coordinates when merging passes of MCGC
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
%	*	db_eps_1/2:		The parameter epsilon to be used in DBScan clustering.
%						If set to zero this switches the DBScan step off
%
%	*	db_pts_1/2:		The number of points that must be in the epsilon 
%						neighbourhood in the DBScan method
%
% 	*	min_pts_1/2:	Minimum number of points needed in each cluster for it
%						not to be rejected
%
%		alo_opt:		Sets the behaviour for remaining points in allometric 
%						feasibility hierarchical clustering
%							0:	Keep these points as their own cluster to also
%									be tested
%							1:	Mark these points are rejected with a 5
%
%	*	subSamp_1/2:	Sets factor by which to subsample data. A value of 5 
%						means 1/5 of the data are used
%
%	*	imp_lim_1/2:		Sets minimum height above ground for points to be
%						included in imputation. It is recommended to set this to
%						match the cut-off height used in MCGC
%
%		sub_db:			toggle to enable use of DBSCAN in subsampled MCGC
%						pipeline or not - 1 to use it, 0 for not. RECOMMENDED
%						TO SET THIS TO 0
%
% Outputs
%
%	LiDAR:			The data as supplied with an extra column at 
%					the end with a cluster label (numbered) for each point
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
%		compute_uncon_weights_nystrom. mcgc also requires select_good_tree
%		If subsampling is used then mcgc_subsample is needed, which has the
%		same dependancy tree as this function 
%
%		mcgc_double_layer
%		 ->
%			(mcgc_subsample)
%			 ->
%				mcgc
%				 ->
%					mcgc_pipeline
%					 ->
%						mcgc_cut
%						 ->
%							nystrom_ext
%						 	 ->
%								compute_uncon_weights_nystrom
%				 ->
%					select_good_tree
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input Checking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	

% Type checking
if(isnumeric(samp_tog)~=1 || size(samp_tog,1) ~= 1 || size(samp_tog,2) ~= 1)
	error('samp_tog must be a single number')
elseif(isnumeric(save_inter)~=1 || size(save_inter,1) ~= 1 ||...
				size(save_inter,2) ~= 1)
	error('save_inter must be a single number')
elseif (isa(out_prec,'char') ~= 1)
	error('out_prec must be a character string')
end

% Sensible values?
if(save_inter~=0 && save_inter~=1)
	error('save_inter must be 1 (save each pass) or 0 (only save final result)')
end

% All arguments present if subsampling used
if(samp_tog==1 && nargin < 46)
	error('Not enough arguments supplied for subsampling, or set samp_tog to 0')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Function Body %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get mfile location
scrpt_loc = fileparts(which(mfilename));
call_loc = pwd;

% Load dependencies
addpath(scrpt_loc);
addpath(genpath([scrpt_loc '\segment']));
addpath(genpath([scrpt_loc '\post_process']));
addpath(genpath([scrpt_loc '\gc_utils']));

if(samp_tog==0)
	% Try first pass
	try
		pass_1 = mcgc(data_folder, out_folder, raw_file, chm_file,...
						pr_grid_size, pr_ht_thr_1, pr_allom, pr_tau,...
						gc_ht_thr_1, nystrom_fac, max_grp_ratio, sigxy_1,...
						sigz_1, allom, r_rat, h_grad_opt, hgrad_w_1,...
						z_grad_opt, zgrad_w_1, merge_tol, test_allom,...
						top_frac, alo_share_max, z_olap, min_alo, db_eps_1,...
						db_pts_1, min_pts_1, alo_opt);
	catch err
		err_handler(err,[call_loc '\' out_folder '\' raw_file...
						'_pass_1_mcgc_err.txt']);
		return
	end
	
	% Save results if requested
	if(save_inter==1)
		dlmwrite([raw_file '_pass_1_result.txt'],pass_1,'precision',out_prec);
	end
	
	% Merge result back into original data
	s = lasread([call_loc '\' data_folder '\' raw_file '.las'],'xyzic','double');
	t = lasread([call_loc '\' data_folder '\' chm_file '.las'],'xyzic','double');	 
	if (sum(s.X==t.X) == size(s.X,1))&&(sum(s.Y==t.Y) == size(s.Y,1))
		LiDAR = [s.X,s.Y,s.Z,t.Z,s.intensity,s.classification];
	else
		[~,las_merge_id]=ismembertol([s.X s.Y s.intensity],[t.X t.Y t.intensity],...
		merge_tol,'ByRows',true,'DataScale',1);
		LiDAR = [s.X,s.Y,s.Z,t.Z(las_merge_id),s.intensity,s.classification];
		clear las_merge_idx
	end
	clear s t
	
	% Add column for tree number
	vec_temp = zeros(size(LiDAR,1),1);
	LiDAR = [LiDAR vec_temp];
	clear vec_temp
	
	[~,merge_idx]=ismembertol(pass_1(pass_1(:,8)==1,1:3),LiDAR(:,1:3),...
								merge_tol,'ByRow',true,'DataScale',1);
	LiDAR(merge_idx,7) = pass_1((pass_1(:,8)==1),7);
	clear pass_1 merge_idx
	
	
	% Create Las files for remaining points
	rem_pts = LiDAR(LiDAR(:,7)==0,:);
	raw_file_2 = [raw_file '_pass_2'];
	chm_file_2 = [chm_file '_pass_2'];
	
	dlmwrite([raw_file_2 '.txt'],rem_pts(:,[1:3,5,6]),'precision',out_prec);
	system(['txt2las -i ' raw_file_2 '.txt -o ' raw_file_2 '.las -parse ' ... 
	'xyzic -odir ' call_loc '\' out_folder]);
	system(['del ' raw_file_2 '.txt']);    

	dlmwrite([chm_file_2 '.txt'],rem_pts(:,[1,2,4,5,6]),'precision',out_prec);
	system(['txt2las -i ' chm_file_2 '.txt -o ' ...
		chm_file_2 '.las -parse xyzic -odir ' call_loc '\' out_folder]);
	system(['del ' chm_file_2 '.txt']);
	
	clear rem_pts
	
	%Try pass two
	try
		pass_2 = mcgc(out_folder, out_folder, raw_file_2, chm_file_2,...
						pr_grid_size, pr_ht_thr_2, pr_allom, pr_tau,...
						gc_ht_thr_2, nystrom_fac, max_grp_ratio, sigxy_2,...
						sigz_2, allom, r_rat, h_grad_opt, hgrad_w_2,...
						z_grad_opt, zgrad_w_2, merge_tol, test_allom,...
						top_frac, alo_share_max, z_olap, min_alo, db_eps_2,...
						db_pts_2, min_pts_2, alo_opt);
	catch err
		err_handler(err,[call_loc '\' out_folder '\' raw_file...
						'_pass_2_mcgc_err.txt']);
		return
	end
	
	% Save results if requested
	if(save_inter==1)
		dlmwrite([raw_file '_pass_2_result.txt'],pass_2,'precision',out_prec);
	end
	
	% ensure no double counting
	shift = max(LiDAR(:,7));
	temp_out = pass_2(:,7);
	pass_2(:,7) = pass_2(:,7) + shift;
	pass_2(temp_out==0,7)=0;
	clear temp_out shift
	
	% Combine results
	[~,merge_idx]=ismembertol(pass_2(pass_2(:,8)==1,1:3),LiDAR(:,1:3),...
								merge_tol,'ByRow',true,'DataScale',1);
	LiDAR(merge_idx,7) = pass_2((pass_2(:,8)==1),7);
	clear pass_2 merge_idx
	
	
	% Relabel trees to go 1:n
	new_label = zeros(size(LiDAR,1),1);
	combos = unique(LiDAR(LiDAR(:,7)~=0,7));
	for i = 1:length(combos)
		ths_com = combos(i);
		new_label(LiDAR(:,7)==ths_com,1)=i;
		clear ths_com
	end
	clear i combos
	LiDAR(:,7)=new_label;
	clear new_label

	dlmwrite([call_loc '\' out_folder '\' raw_file ...
	'_final_double_layer.txt'],LiDAR,'precision',out_prec);
	
	
% Instead use subsampling
elseif(samp_tog==1)
	% Try first pass
	try
		pass_1 = mcgc_subsample(subSamp_1, imp_lim_1, sub_db, out_prec,...
						data_folder, out_folder, raw_file, chm_file,...
						pr_grid_size, pr_ht_thr_1, pr_allom, pr_tau,...
						gc_ht_thr_1, nystrom_fac, max_grp_ratio, sigxy_1,...
						sigz_1, allom, r_rat, h_grad_opt, hgrad_w_1,...
						z_grad_opt, zgrad_w_1, merge_tol, test_allom,...
						top_frac, alo_share_max, z_olap, min_alo, db_eps_1,...
						db_pts_1, min_pts_1, alo_opt);
	catch err
		err_handler(err,[call_loc '\' out_folder '\' raw_file...
						'_pass_1_mcgc_err.txt']);
		return
	end
	
	% Save results if requested
	if(save_inter==1)
		dlmwrite([raw_file '_pass_1_result.txt'],pass_1,'precision',out_prec);
	end
	
	% Merge result back into original data not needed as already merged in
	% subsamp
	
	LiDAR = pass_1;
	clear pass_1
	
	% Create Las files for remaining points
	rem_pts = LiDAR(LiDAR(:,7)==0,:);
	raw_file_2 = [raw_file '_pass_2'];
	chm_file_2 = [chm_file '_pass_2'];
	
	dlmwrite([raw_file_2 '.txt'],rem_pts(:,[1:3,5,6]),'precision',out_prec);
	system(['txt2las -i ' raw_file_2 '.txt -o ' raw_file_2 '.las -parse ' ... 
	'xyzic -odir ' call_loc '\' out_folder]);
	system(['del ' raw_file_2 '.txt']);    

	dlmwrite([chm_file_2 '.txt'],rem_pts(:,[1,2,4,5,6]),'precision',out_prec);
	system(['txt2las -i ' chm_file_2 '.txt -o ' ...
		chm_file_2 '.las -parse xyzic -odir ' call_loc '\' out_folder]);
	system(['del ' chm_file_2 '.txt']);
	
	clear rem_pts
	
	%Try pass two
	try
		pass_2 = mcgc_subsample(subSamp_2, imp_lim_2, sub_db, out_prec,...
						out_folder, out_folder, raw_file_2, chm_file_2,...
						pr_grid_size, pr_ht_thr_2, pr_allom, pr_tau,...
						gc_ht_thr_2, nystrom_fac, max_grp_ratio, sigxy_2,...
						sigz_2, allom, r_rat, h_grad_opt, hgrad_w_2,...
						z_grad_opt, zgrad_w_2, merge_tol, test_allom,...
						top_frac, alo_share_max, z_olap, min_alo, db_eps_2,...
						db_pts_2, min_pts_2, alo_opt);
	catch err
		err_handler(err,[call_loc '\' out_folder '\' raw_file...
						'_pass_2_mcgc_err.txt']);
		return
	end
	
	% Save results if requested
	if(save_inter==1)
		dlmwrite([raw_file '_pass_2_result.txt'],pass_2,'precision',out_prec);
	end
	
	% ensure no double counting
	shift = max(LiDAR(:,7));
	temp_out = pass_2(:,7);
	pass_2(:,7) = pass_2(:,7) + shift;
	pass_2(temp_out==0,7)=0;
	clear temp_out shift
	
	% Combine results
	LiDAR(LiDAR(:,7)==0,7)=pass_2(:,7);
	clear pass_2
	
	% Relabel trees to go 1:n
	new_label = zeros(size(LiDAR,1),1);
	combos = unique(LiDAR(LiDAR(:,7)~=0,7));
	for i = 1:length(combos)
		ths_com = combos(i);
		new_label(LiDAR(:,7)==ths_com,1)=i;
		clear ths_com
	end
	clear i combos
	LiDAR(:,7)=new_label;
	clear new_label

	dlmwrite([call_loc '\' out_folder '\' raw_file ...
	'_final_double_layer.txt'],LiDAR,'precision',out_prec);

else
	error('subSamp must be 0 (no use) or 1 (subsample)');
end
clear pass_num
end