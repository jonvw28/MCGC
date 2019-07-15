function [LiDAR,sample_idx] = mcgc_subsample(subSamp, imp_lim, sub_db,...
	out_prec, data_folder, out_folder, raw_file, chm_file, pr_grid_size,...
	pr_ht_thr, pr_allom, pr_tau, gc_ht_thr, nystrom_fac, max_grp_ratio,...
	sigxy, sigz, allom, r_rat, h_grad_opt, hgrad_w, z_grad_opt, zgrad_w,...
	merge_tol, test_allom, top_frac, alo_share_max, z_olap, min_alo, db_eps,...
	db_pts,	min_pts, alo_opt)
%MCGC_SUBSAMPLE wrapper for MCGC to allow subsampling of data 
%  
% Subsample data before applying MCGC. Then upsample results using k-nearest
% neighbours before applying DBSCAN at the end to ensure connectivity off
% crowns
%
%	PLEASE NOTE: For this to work, it is necessary to have Rscript
%				added to the system search path. You will also need to have
%				installed the itcSegment package in R
%
% Syntax
%
%	[LiDAR,sample_idx] = mcgc_subsample(subSamp, imp_lim, sub_db,...
%		data_folder, out_folder, raw_file, chm_file, pr_grid_size, pr_ht_thr,...
%		pr_allom, pr_tau, gc_ht_thr, nystrom_fac, max_grp_ratio, sigxy, sigz,...
%		allom, r_rat, h_grad_opt, hgrad_w, z_grad_opt, zgrad_w, merge_tol,...
%		test_allom, top_frac, alo_share_max, z_olap, min_alo, db_eps, db_pts,...
%		min_pts, alo_opt,out_prec)
%
%
% Method
%
%	This applies the MCGC complete pipeline as provided by the mcgc function
%	after subsampling. First the data is randomly subsampled by a factor set by
%	subSamp. If subSamp = 5 then 1/5 (20%) of the data will be used. The full
%	pipeline in mcgc will then be applied to this subsampled data. The results
%	of this are then used to impute crowns for the full dataset based on
%	k-nearest neighbours. The value of k is set to subSamp (so if 20% of the 
%	are used, k will be set to 5). This upsampling is constrained by the
%	allometric table used in the allometric checking to ensure imputation of
%	crowns doesn't extend beyond this radius. Additionally imputation only
%	occurs for points which are above a threshold height above gorund, set by
%	imp_lim. This should be set to the same as the threshold for MCGC but
%	is optionally able to be changed. Density dependent arguments are 
%	automatically adjusted for the application of MCGC so should be set as if
%	being used on the complete data.
%
%	After imputation DBSCAN is applied to help avoid artefacts developing in
%	the imputation phase. A final check of number of points in each crown is
%	then also applied before returning the final results in LiDAR, along
%	with the row indices in LiDAR of the subsampled points that are used.
%
%	The toggle sub_db determines if DBSCAN will be used within the subsampled 
%	MCGC step. This is given as a choice as DBSCAN is applied at the end of
%	the imputation step. This is a quality control step so applying this in the
%	usual MCGC pipeline is superfluous. It is recommended to avoid using
%	this intermediate checking (sub_db=0) as the intended use of DBSCAN is as
%	a final quality control step but this is left optional.
%
%	In addition to the outputs this function will additionally save 3 files.
%	Two .las files will be automatically created in the output folder. These
%	will contain only the subsampled points for the raw and chm point clouds.
%	The precision of the output can be set by out_prec. This uses a REGEX
%	formula to specify the precision of the coordinates to be saved. It is
%	heavily recommended that you set this to '%.2f' for float values accurate
%	to cm resolution. Additionally a .csv of the row indices of the subsampled
%	points will be saved in the output folder.
%
%
%	Should the MCGC step fail, a file will automatically be generated with
%	details of the error. This will end with '_sub_mcgc_err.txt'. Should the
%	imputation step fail an error report will be saved with the ending
%	'_knn_err.txt'. If the final DBSCAN and minimum points check fails then
%	a message will be saved to a file ending '_DBSCAN_err.txt'.
%
%	NB: Any error file ending '_seg_err.txt' or '_alo_check_err.txt' are thrown
%	by the function mcgc and should result in a throwing of an error in the
%	MCGC section of this function.
%
%	This function gives access to all available arguments across all steps of 
%	the algorithm. For simpler forms of the pipeline see the wrappers which
%	preset many arguments.
%
%	The additonal arguments presented by the subsampling regime are:
%	subSamp - controlling subsampling factor - setting 5 uses 1/5 of the data)
%	imp_lim - setting minimum height for imputation
%	sub_db - toggle to use DBSCAN in subsampled MCGC (0 - no, 1 - yes)
%	out_prec - REGEX setting precision to save coordinates of subsampled las 
%				files
%
%
% Inputs
%
%		subSamp:		Sets factor by which to subsample data. A value of 5 
%						means 1/5 of the data are used
%
%		imp_lim:		Sets minimum height above ground for points to be
%						included in imputation. It is recommended to set this to
%						match the cut-off height used in MCGC
%
%		sub_db:			toggle to enable use of DBSCAN in subsampled MCGC
%						pipeline or not - 1 to use it, 0 for not. RECOMMENDED
%						TO SET THIS TO 0
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
%		merge_tol:		Sets the absolute tolerances used when comparing
%						coordinates when merging subsampled data back into 
%						complete dataset between the two data sources.
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
%	LiDAR:			The data as supplied with an extra column at 
%					the end with a cluster label (numbered) for each point
%
%	sample_idx:		Row indices in LiDAR for the subsampled datapoints
%
%
%			The prior that is generated will also be saved in the output folder
%			with the name 'chm_file'_prior.csv as will the subsampled las files
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
%		This function requires mcgc which uses mcgc_pipeline requiring mcgc_cut,
%		which also requires nystrom_ext which in turn requires
%		compute_uncon_weights_nystrom. mcgc also requires select_good_tree
%
%		mcgc_subsample
%		 ->
%			mcgc
%			 ->
%				mcgc_pipeline
%				 ->
%					mcgc_cut
%					 ->
%						nystrom_ext
%						 ->
%							compute_uncon_weights_nystrom
%			 ->
%				select_good_tree
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
if(isnumeric(subSamp)~=1 || size(subSamp,1) ~= 1 || size(subSamp,2) ~= 1)
	error('subSamp must be a single number')
elseif(isnumeric(imp_lim)~=1 || size(imp_lim,1) ~= 1 || size(imp_lim,2) ~= 1)
	error('imp_lim must be a single number')		
elseif (isa(out_prec,'char') ~= 1)
	error('out_prec must be a character string')
end

% Sensible values?
if (subSamp <= 1)
	error('subSamp must be >= 1')
elseif (imp_lim <= 0)
	error('imp_lim must be a postive number')
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

%%% Subsample data

% first read in data
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

% Subsample it
subSamp = round(subSamp);
raw_pts = size(LiDAR,1);
[sub_LiDAR, sample_idx] = datasample(LiDAR,round(raw_pts/subSamp),1,'replace',false);

% Save out points and generate subsampled LAS files 
dlmwrite([call_loc '\' out_folder '\' raw_file '_subsamp_raw.txt'],...
	sub_LiDAR(:,[1:3,5,6]),'precision',out_prec);
system(['txt2las -i ' call_loc '\' out_folder '\' raw_file...
	'_subsamp_raw.txt -o ' raw_file '_subsamp_raw.las -parse xyzic -odir '...
	call_loc '\' out_folder]);
system(['del ' call_loc '\' out_folder '\' raw_file '_subsamp_raw.txt']);    

dlmwrite([call_loc '\' out_folder '\' chm_file '_subsamp_chm.txt'],...
	sub_LiDAR(:,[1,2,4,5,6]),'precision',out_prec);
system(['txt2las -i ' call_loc '\' out_folder '\' chm_file...
	'_subsamp_chm.txt -o ' chm_file '_subsamp_chm.las -parse xyzic -odir '...
	call_loc '\' out_folder]);
system(['del ' call_loc '\' out_folder '\' chm_file '_subsamp_chm.txt']);    


%%% Run MCGC
% Update minimum points number
min_sub_pts = min_pts/subSamp;
% Either turn DBscan off or els eupdate the parameters
if(sub_db == 0)
	db_sub_eps = 0;
	db_sub_pts = db_pts;
elseif(sub_db==1)
	db_sub_pts = db_pts/subSamp;
	% Deal with case db_pts becomes too small by increasing search area
	if ((db_sub_pts~=0)&&(db_sub_pts<1))
		db_sub_eps = sqrt(db_eps/db_sub_pts);
		db_sub_pts = 1;
	end
	db_sub_pts = round(db_sub_pts);
else
	error('sub_db must be 0 (off) or 1 (on)')
end


try
	subSeg = mcgc(out_folder,out_folder,[raw_file '_subsamp_raw'],...
		[chm_file '_subsamp_chm'],pr_grid_size,pr_ht_thr,pr_allom,pr_tau,...
		gc_ht_thr,nystrom_fac,max_grp_ratio,sigxy,sigz,allom,r_rat,...
		h_grad_opt,hgrad_w,z_grad_opt,zgrad_w,merge_tol,test_allom, top_frac,...
		alo_share_max,z_olap, min_alo, db_sub_eps, db_sub_pts, min_sub_pts,...
		alo_opt);
catch err
	err_handler(err,[call_loc '\' out_folder '\' raw_file '_sub_mcgc_err.txt']);
	return
end


%%% k-NN upsample
try
	% Add column for tree number
	vec_temp = zeros(size(LiDAR,1),1);
	LiDAR = [LiDAR vec_temp];
	clear vec_temp

	% Add in subsampled data (las files can scramble the order)
	% Keeping only points which pass allometric checking
	[~,merge_idx]=ismembertol(subSeg((subSeg(:,8)==1),1:3),LiDAR(:,1:3),...
								merge_tol,'ByRow',true,'DataScale',1);
	LiDAR(merge_idx,end) = subSeg((subSeg(:,8)==1),7);
	clear subSeg

	% set up k-NN
	knn_X = LiDAR(LiDAR(:,7)~=0,1:3);
	knn_Y = LiDAR(LiDAR(:,7)~=0,7);
	sub_hts = LiDAR(LiDAR(:,7)~=0,4);
	knn_Mdl = fitcknn(knn_X,knn_Y,'NumNeighbors',subSamp);

	% Find indices to be predicted
	ptsleft = find(LiDAR(:,end)==0);
	tall_left = LiDAR(ptsleft,4)>imp_lim;
	LiDAR2predict = ptsleft(tall_left,:);
	clear ptsleft tall_left
	pts_predict = size(LiDAR2predict,1);

	%% Collect a list of current tree tops - to constrain Knn

	% For future reference - what's the tallest height in allom lookup
	max_allom = max(test_allom(:,1));

	% Number of Trees
	trees = sort(unique(knn_Y));
	n_tree = size(trees,1);
	% get a list of all tree 'tops' and their heights
	tops = zeros(n_tree,4);
	for a=1:n_tree
		% extract tree with both z values
		class_a = trees(a);
		cur_tree = [knn_X(knn_Y==class_a,:) sub_hts(knn_Y==class_a)];
		
		% Get centroid of points within top_frac of top height - assumed centre
		[top_ht,~] = max(cur_tree(:,4));
		top_pts = cur_tree(cur_tree(:,4)>=top_frac*top_ht,1:2);
		cur_cen = mean(top_pts,1);
		
		% get allom lookup radius
		top_ht = round(top_ht);
		max_z = min([top_ht max_allom]);
		alo_d = allom((test_allom(:,1)==max_z),2);
		alo_d = alo_d^2;
		
		% Save the tree info
		tops(a,:) = [cur_cen top_ht alo_d];
	end
	clear top_ht cur_tree a class_a top_pts knn_X knn_Y alo_d cur_cen sub_hts
	clear max_z

	% Constrain the size of each block of k-NN fitting to control memory issues
	knn_passes = floor(pts_predict/1000);

	% Apply k-NN by block (no tree updates)
	for knn_i = 1:knn_passes
		% predict trees
		pass = LiDAR(LiDAR2predict(((knn_i-1)*1000+1):((knn_i-1)*1000+1000)),...
		1:3);
		knn_output = predict(knn_Mdl,pass);
		tree_idx = zeros(1000,1);
		
		% Match each point to its crown index
		for b = 1:n_tree
			tree_idx((knn_output==trees(b))) = b;
		end
		clear b
		% get indices for each point in D matrix
		d_inds = sub2ind([1000,n_tree],1:1000,tree_idx');
		
		% Generate 2d Euclidean matrix
		D_mat = repmat(pass(:,1:2),n_tree,1) - repelem(tops(:,1:2),1000,1);
		D_mat = sum(D_mat.^2,2);
		D_mat = reshape(D_mat,1000,n_tree);
		
		% get squared Ds
		D_pts = reshape(D_mat(d_inds),1000,1);
		clear d_inds D_mat
		D_trees = reshape(tops(tree_idx,4),1000,1);
		
		mask = D_trees>D_pts;
		knn_output(~mask) = 0;
		
		LiDAR(LiDAR2predict(((knn_i-1)*1000+1):((knn_i-1)*1000+1000)),7)=...
																	knn_output;
		clear tree_idx pass knn_output mask D_trees D_pts
	end
	clear knn_i

	% Deal with remaining points
	% predict trees
	pass = LiDAR(LiDAR2predict((knn_passes*1000+1):pts_predict),1:3);
	knn_output = predict(knn_Mdl,pass);
	pass_pts = size(pass,1);
	tree_idx = zeros(pass_pts,1);

	% Match each point to its crown
	for b = 1:n_tree
		tree_idx((knn_output==trees(b))) = b;
	end
	clear b
	% get indices for each point in D matrix
	d_inds = sub2ind([pass_pts,n_tree],1:pass_pts,tree_idx');

	% Generate 2d Euclidean matrix
	D_mat = repmat(pass(:,1:2),n_tree,1) - repelem(tops(:,1:2),pass_pts,1);
	D_mat = sum(D_mat.^2,2);
	D_mat = reshape(D_mat,pass_pts,n_tree);

	% get squared Ds
	D_pts = reshape(D_mat(d_inds),pass_pts,1);
	clear d_inds D_mat
	D_trees = reshape(tops(tree_idx,4),pass_pts,1);

	mask = D_trees>D_pts;
	knn_output(~mask) = 0;

	LiDAR(LiDAR2predict((knn_passes*1000+1):pts_predict),7)=knn_output;	
	clear tree_idx pass knn_output mask knn_Mdl LiDAR2predict knn_passes
	clear D_pts D_trees pass_pts
catch err
	err_handler(err,[call_loc '\' out_folder '\' raw_file '_knn_err.txt']);
	return
end



% Apply DBSCAN if requested (ie db_eps ~= 0)
% Get all tree ids
try
	tr_ids = unique(LiDAR(:,7));
		for(db_i=1:length(tr_ids))
			% Rejected points
			if(tr_ids(db_i)==0)
				continue
			end
			this_tree = LiDAR(LiDAR(:,7)==tr_ids(db_i),:);
			[~,tp_idx] = max(this_tree(:,3));
			db_out = DBSCAN(this_tree(:,1:3),db_eps,db_pts); 
				
			% keep cluster with tallest point and add to input
			keep_clust = db_out(tp_idx); 
			new_clust = zeros(size(db_out,1),1);
			new_clust(db_out==keep_clust) = tr_ids(db_i);
			LiDAR(LiDAR(:,end)==tr_ids(db_i),7) = new_clust;
			clear keep_clust new_clust db_out this_tree tp_idx
		end
	clear tr_ids

	% Now drop trees which are too small
	tr_ids = unique(LiDAR(:,end));
	for pt_i = 1:size(tr_ids,1)
		ths_id = tr_ids(pt_i);
		ths_tr = LiDAR(LiDAR(:,7)==ths_id,:);
		if (size(ths_tr,1)< min_pts)
			LiDAR(LiDAR(:,7)==ths_id,7)=0;
		end
	end
	clear pt_i tr_ids ths_id ths_tr
catch err
	err_handler(err,[call_loc '\' out_folder '\' raw_file '_DBSCAN_err.txt']);
	return
end

% Save sample indices in final table of LiDAR data
csvwrite([call_loc '\' out_folder '\' raw_file ...
	'_MCGC_subsamp_ids.txt'],sample_idx);
end