function [ outSeg ] = mcgc_pipeline(data_folder, output_folder, raw_file,...
	chm_file, pr_grid_size, pr_ht_thr, pr_allom, pr_tau, gc_ht_thr,...
	nystrom_fac, max_grp_ratio, sigxy, sigz, allom, r_rat, h_grad_opt,...
	hgrad_w, z_grad_opt, zgrad_w, merge_tol)
%MCGC_PIPELINE Function to call full mcgc pipeline including prior generation
%   
%	PLEASE NOTE: For this to work, it is necessary to have Rscript
%				added to the system search path. You will also need to have
%				installed the itcSegment package in R
%
% Syntax
%
%		outSeg = mcgc_pipeline(data_folder, output_folder, raw_file,...
%				chm_file, pr_grid_size, pr_ht_thr, pr_allom, pr_tau,...
%				gc_ht_thr, nystrom_fac, max_grp_ratio, sigxy, sigz, allom,...
%				r_rat, h_grad_opt, hgrad_w, z_grad_opt, zgrad_w, merge_tol)
%		
% 		This returns the data stored in the LAS file given by 'file' (without
%		.las extension) and its graph cut segmentation in a matrix form, with
%		each point on a row. This takes care of prior generation
%		and uses normalised graph cut using the nystrom extention.
%
%
% Method
%
%		Here the data in LiDAR is used to generate a set of prior tree tops
%		using a sliding window maxima finder from the itcSegment R package.
%		This is then passed to the MCGC implementation in MATLAB. Here the 
%		points are clustered based on a multi class graph cut approach. This
%		uses a normalised cut as outlined in [1]. Details of the clustering
%		method can be found in the documentation of mcgc_cut. This process is
%		acclearated by the nystrom extension as outlined in [2] and adapated
%		from the code base in [3]
%
%		Finally, the resulting segmentation is returned as outSeg, being the 
%		original data with an extra column at the end for cluster identification
%		With each cluster having a unique number here.
%
%		The additional arguments introduced here which are not used by
%		dependancies are:
%			data_folder
%			output_folder
%			raw_file
%			chm_file
%			pr_grid_size
%			pr_ht_thr
%			pr_allom
%			pr_tau
%			merge_tol
%
%      
% Inputs
%
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
%
% Outputs: 
%
%		outSeg:		The data as supplied (without points below the
%					minimum height threshold removed) with an extra column at 
%					the end with a cluster label (numbered) for each point
%
%					The first two columns will be the x and y co-ordinates,
%					followed by the raw and above ground heights then intensity
%					and original classification in columns 3-6. The final column
%					will then contain a numeric label for the cluster of each 
%					point
%
%			The prior that is generated will also be saved in the output folder
%			with the name 'chm_file'_prior.csv
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
%		This function requires mcgc_cut, which requires nystrom_ext which in 
%		turn requires compute_uncon_weights_nystrom.
%
%		mcgc_pipeline
%		 ->
%			mcgc_cut
%			 ->
%				nystrom_ext
%				 ->
%					compute_uncon_weights_nystrom
%
%
%		Jonathan Williams
%		jonvw28@gmail.com			         
%		09/01/2019		         
%
%		Dependancies heavily influenced by the example code in [2] and adapted
%		from the code used originally in [3].
%
%		This work makes use of the lasread MATLAB function oriinally written by
%		Thomas J. Pingel. See the address below for its source
%		http://www.tpingel.org/code/lasread/lasread.html


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input Checking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		     

% Type checking
	
if (isa(data_folder,'char') ~= 1)
	error('data_folder must be a character string')
elseif (isa(output_folder,'char') ~= 1)
	error('output_folder must be a character string')
elseif (isa(raw_file,'char') ~= 1)
	error('raw_file must be a character string')
elseif (isa(chm_file,'char') ~= 1)
	error('chm_file must be a character string')
elseif(isnumeric(pr_grid_size)~=1 || size(pr_grid_size,1) ~= 1 || ...
		size(pr_grid_size,2) ~= 1)
	error('pr_grid_size must be a single number')
elseif(isnumeric(pr_ht_thr)~=1 || size(pr_ht_thr,1) ~= 1 || ...
		size(pr_ht_thr,2) ~= 1)
	error('pr_ht_thr must be a single number')		
elseif(isnumeric(merge_tol)~=1 || size(merge_tol,1) ~= 1 ...
		|| size(merge_tol,2) ~= 1)
		error('merge_tol must be a single number')
end
% Feasibility Checking

% no .las please
if (strcmp(raw_file((end-3):end),'.las')==1)
	error('raw_file must be specified without the .las extension')
end
if (strcmp(chm_file((end-3):end),'.las')==1)
	error('chm_file must be specified without the .las extension')
end
% Sensible values?
if (pr_grid_size <= 0)
	error('pr_grid_size must be a postive number')
elseif (pr_ht_thr <= 0)
	error('pr_ht_thr must be a postive number')
elseif (nargin <20)
	merge_tol = 9e-3;
end
			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Function Body %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get mfile location
scrpt_loc = fileparts(which(mfilename));
call_loc = pwd;

% Load dependencies
addpath(genpath([scrpt_loc '\..\gc_utils']));
addpath(genpath([scrpt_loc '\..\..\utils']));
addpath(genpath([scrpt_loc '\nystrom']));

% Generate and read in prior
system([scrpt_loc '\..\..\..\batch\graph_cut\mcgc_itcseg_prior.bat ' ...
	call_loc '\' data_folder ' ' call_loc '\' output_folder ' ' ...
	chm_file '.las' ' ' pr_allom ' ' num2str(pr_grid_size) ' ' ...
	num2str(pr_ht_thr) ' ' pr_tau]);
priorFile = [call_loc '\' output_folder '\' chm_file '_prior.csv'];
prior = csvread(priorFile,1,0);
    
% read in the LAS data and add cluster column
s = lasread([call_loc '\' data_folder '\' raw_file '.las'],...
	'xyzic','double');
t = lasread([call_loc '\' data_folder '\' chm_file '.las'],...
	'xyzic','double');
if (sum(s.X==t.X) == size(s.X,1))&&(sum(s.Y==t.Y) == size(s.Y,1))
	raw_data = [s.X,s.Y,s.Z,t.Z,s.intensity,s.classification];
else
	[~,las_merge_id]=ismembertol([s.X s.Y s.intensity],[t.x t.Y t.intensity],...
	merge_tol,'ByRows',true,'DataScale',1);
	raw_data = [s.X,s.Y,s.Z,t.Z(las_merge_id),s.intensity,s.classification];
	clear las_merge_id
end
    
% apply cut    
try
	outSeg = mcgc_cut(raw_data,3,4, prior, gc_ht_thr, nystrom_fac,...
		max_grp_ratio, sigxy, sigz, allom, r_rat, h_grad_opt, hgrad_w,...
		z_grad_opt, zgrad_w);
catch err
	err_handler(err,...
		[call_loc '\' output_folder '\' raw_file '_gcseg_err.txt'])
end
end
