function [ outSeg ] = mcgc_cut( LiDAR, raw_z_col, chm_z_col, prior,...
			grndThr, nystrom_fac, max_grp_ratio, sigxy, sigz, allom, r_rat,...
			h_grad_opt, hgrad_w, z_grad_opt, zgrad_w)
%MCGC_cut Applies unconstrained multiclass normalised graph cut for tree 
% segmentation
%
% Syntax
%
%		outSeg = mcgc_cut( LiDAR, raw_z_col, chm_z_col, prior, grndThr,
%							nystrom_fac, max_grp_ratio, sigxy, sigz, allom,
%							r_rat, h_grad_opt, hgrad_w, z_grad_opt, zgrad_w)
%		
% 		This returns the data in LiDAR appended with an extra column labelling
%		each row (ie point) with a cluster generarted via a multiclass
%		normlaised graph cut
%
% Method
%
%		Here the data in LiDAR is clustered based on a multi class graph cut
%		approach. This makes use of a normalised cut as outlined in [1] using
%		the nystrom extension as set out in [2]. First points below a minimum
%		height are removed (to avoid inclusion of ground returns and very small
%		vegetation). Then the relvant adjacency matrices are constructed,
%		details	of this are included in the documentation for 
%		compute_uncon_weights_nystrom.
%
%		Normalised cut is enacted as per the method using L_{sym} in [2]. The 
%		number of clusters sought is flexible. A minimum is set as the number
%		number of prior tree tops (so be selective with these), and a maximum
%		is set as a multiple of this, max_grp_ratio. The relevant eigenvectors
%		for all of these are computed, and the final number is set by the
%		maximum difference between subsequent eigenvalues in this range. k-means
%		clustering is then used on these eigenvalues (per [1]) to produce 
%		clusters.
%
%		The number of points to subsample in the Nystrom extension is set via
%		the parameter nystrom_fac. Here	the maximum number of clusters sought is
%		multiplied by nystrom_fac to get the size of the subsample. For more
%		details see the docs for nystrom_ext or [2]. The implementation of
%		this methods is adapted from the code base for the work in [3].
%
%		Finally, the resulting segmentation is returned as outSeg, being the 
%		original data with an extra column at the end for cluster identification
%		With each cluster having a unique number here
%      
%		The extra parameters introduced here which aren;t used by dependancies
%		are:
%			prior
%			grndThr
%			nystrom_fac
%			max_grp_ratio
%
% Inputs
%
%		LiDAR:			A matrix representation of a point cloud where each 
%						point is in its own row, and the first 2 columns are the
%						x and y co-ordinates respectively
%
%		raw_z_col:		The column containing the raw z value for each point
%
%		chm_z_col:		The column which contains the value of height above the 
%						ground for each point
%
%		prior:			A matrix containing information on suspected tree top 
%						locations. Must be three columns only, being the x, y 
%						and z coordinates of each tree top respectively
%
%		grndThr:		Cut-off height below which points are ignored and not 
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
%       allom:			Allometric lookup table for centroid computation. Must 
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
%
% Outputs: 
%
%		outSeg:		The data as supplied in LiDAR (without points below the
%					minimum height threshold removed) with an extra column at 
%					the end with a cluster label (numbered) for each point
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
%		This function requires nystrom_ext which in turn requires 
%		compute_uncon_weights_nystrom. This function is required by 
%		mcgc_pipeline
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
%		from the code used originally in [3]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input Checking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Type checking
if(ismatrix(LiDAR)~= 1 || isnumeric(LiDAR(:,[1,2,raw_z_col,chm_z_col]))~=1)
	error(['LiDAR must be a matrix, where the first two columns are ' ...
			'numeric values of x,y co-ordinates for the points'])
elseif(ismatrix(prior)~= 1 || isnumeric(prior(:,1:3))~=1)
	error(['prior must be a matrix, where the first three columns are ' ...
			'numeric values of x,y,z co-ordinates for the points'])
elseif(isnumeric(grndThr)~=1 || size(grndThr,1) ~= 1 || size(grndThr,2) ~= 1)
	error('grndThr must be a single number')			
elseif(isnumeric(nystrom_fac)~=1 || size(nystrom_fac,1) ~= 1 ...
		|| size(nystrom_fac,2) ~= 1)
	error('nystrom_fac must be a single number')
elseif(isnumeric(max_grp_ratio)~=1 || size(max_grp_ratio,1) ~= 1 ...
		|| size(max_grp_ratio,2) ~= 1)
	error('max_grp_ratio must be a single number')
end

% Check there are points
if(size(LiDAR,1)==0)
	error('data can not be empty')
end

% Check prior is of form expected
if(size(prior,2)~=3)
	error('prior must only consist of three columns, being x,y & z coordinates')
elseif(size(prior,1)==0)
	error('prior can not be empty')
end

% Check validity of relevant arguments
if (grndThr <= 0)
	error('grndThr must be a postive number')
elseif (nystrom_fac <= 0)
	error('nystrom_fac must be a postive number')
elseif (max_grp_ratio < 1)
	error('max_grp_ratio must be a number bigger than 1')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Function Body %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get mfile location
scrpt_loc = fileparts(which(mfilename));

% Load dependencies
addpath(genpath([scrpt_loc '\nystrom']));
addpath(genpath([scrpt_loc '\..\gc_utils']));

% Filter points by a minimum height threshold
gd = LiDAR(:,chm_z_col) < grndThr;
LiDAR(gd,:)=[];

% Number of trees from prior
numPr = size(prior,1);

% Max Number of Classes
numCl = floor(numPr*max_grp_ratio);

% Run algorithm
disp('Computing Eigenvectors...');
[ncVec,ncEig] = nystrom_ext(LiDAR,raw_z_col,chm_z_col,...
								nystrom_fac*numCl,sigxy,sigz,allom,r_rat,...
								h_grad_opt,hgrad_w,z_grad_opt,zgrad_w,numCl);
disp('Eigenvectors Computed!');

% Use Spectral Gap to select number of clusters
evals = diag(ncEig);

if(numPr == 1) % Special Case
	if(evals(1)>(evals(2)-evals(1)))
		clusts = 1;
	else
		clusts = 2;
	end
else
	[~,clusts]=max(evals(numPr:end) - evals((numPr-1):(end-1)));
	clear evals
	clusts = clusts -1 + numPr;
end

% classify each point
disp('Segmenting Point Cloud...');
label = v2idx(ncVec(:,1:clusts),clusts,0);
disp('Segmentation Complete!');

% Return result with new classification
outSeg = [LiDAR,label];
end