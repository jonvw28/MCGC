function [ outFilt ] = select_good_tree( LiDAR, raw_z_col, chm_z_col,...
						class_col, allom_lkup, top_frac, alo_share_max,...
						z_olap, min_alo, db_eps, db_pts, min_pts, alo_opt)
%SELECT_GOOD_TREE Function to assess quality of MCGC crowns
%
%
% Syntax
%
%		outFilt = select_good_tree_square(LiDAR, raw_z_col, chm_z_col,
%						class_col, allom_lkup, top_frac, alo_share_max, z_olap,
%						min_alo, db_eps, db_pts, min_pts, alo_opt)
%		
% 		This returns the data in LiDAR appended with an extra column labelling
%		each cluster (given in class_col) with the quality of the segmentation
%
% Method
%
%		Here the data in LiDAR is subset by the numeric cluster labels given in
%		the class column (class_col). Each of these is then assessed on various
%		quality controls to decide if the segmentation is believed or not.
%
%		First the top of each cluster is computed as being at the height of 
%		the tallest point above the ground, this data being given in chm_z_col,
%		and at the x-y centroid of all points which are at
%		least top_frac proportion of this height tall. These are then used
%		in combination with the lookup table for tree height (rounded to the 
%		nearest metre) to crown radius supplied in allom_lkup. Starting from the
%		tallest tree, all other	smaller	tree-tops that are within the look up 
%		radius of the tree being tested are marked. Additionally, any other 
%		smaller tree for which at least alo_share_max proportion of their points
%		are within this radius of the tall tree's radius are also marked.
%
%		For each marked tree it is then decided whether to merge this cluster 
%		with the taller tree currently being used for comparison. If the small
%		tree has fewer than min_pts, then it is always merged into the taller
%		tree. Otherwise the distribution of heights of each tree is compared.
%		The (100*z_olap)-th percentile of point heights in the tall tree is 
%		compared to the (100-100*z_olap)-th percentile of points heights in the 
%		short tree. If the latter is taller than the former, then the trees are 
%		merged, otherwise the short tree is merged into the tall tree. 
%		For example, if z_olap is 0.1, a marked tree will be merged if the
%		lowest 10% of points of the taller tree are all below to tallest 10% of
%		points of the smaller tree.
%
%		(Note that in the merge phase, the taller tree when compared, is the
%			version of that tree before any prior merges occured)
%
%		After the merging phase, the function then tests allometric feasibility.
%		Working through the candidate clusters from the previous phase, in
%		descending height order, the location of the tree top is re-computed in 
%		same manner as before. The points in each cluster are then compared
%		to the radius given in the allometric look up table for a tree of this
%		height. A cluster is considered allometrically feasible if at least 
%		min_alo proportion of its points are within this radius, otherwise it is
%		further split to produce new candidate clusters.
%
%		Where necessary, an allometrically infeasible cluster is divided using
%		hierarchical clustering on the raw point cloud into exactly 2 clusters.
%		The cluster containing the highest point is then kept as the updated 
%		version of the 'tree'. The points in the the other cluster are then
%		treated in one of two ways based on the value of alo_opt:
%
%			alo_opt = 0 (default) sees these points appended to the end of 
%			the list of candidate trees, which will in turn be tested for 
%			allometric feasibility.
%			
%			alo_opt = 1 instead marks these points as being rejected in 
%			allometric filtering, resulting in a value of 5 in the added column
%			returned at the end. This enables these points to be pooled for
%			future passes with the graph cut algorithm
%
%		After each pass of this binary hierarchical clustering, the points
%		remaining in the original cluster label are then tested again for 
%		allometric feasibility, and further binary clustering is applied as
%		necessary. Should the number of points in the cluster being tested fall
%		below the minimum (min_pts) then this process is halted as the cluster
%		will later be rejected anyway.
%
%		In order to ensure that trees are made up of joined clusters of points,
%		a DBScan pass is run on the data (unless db_eps is 0) and only those
%		points in the same grouping as the top point of each tree are kept.
%		The remaining points are rejected and re-used in future passes of the
%		algorithm
%
%		Once this is completed, each cluster is then in turn tested for whether
%		it contains enough points, set by a minimum number of points (min_pts),
%		with the exception of any points rejected already if alo_opt = 1. Any 
%		cluster with fewer points is marked as such with a value of 2 in the
%		output column and is ignored from here on. 
%
%		After all of the above, any cluster which has yet to be marked as having 
%		an issue is then marked with a 1 to indicate that it passed all tests.
%
%		Finally, the resulting assesment is returned as outFilt being the 
%		original data with an extra column at the end for assesment status.
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
%		class_col:		The index of the column in LiDAR containing the numeric 
%						cluster label
%
%		allom_lkup:		An N-by-2 numeric lookup table for allometry. The first 
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
% Outputs: 
%
%		outFilt:	The data as supplied in LiDAR with an extra column 
%					indicating the assesment of the cluster for each point
%					(where the assesment is shared across all points in a given
%					cluster). The available options are:
%						0: There is an error in the code (please contact me)
%						1: There is no issue with the cluster
%						2: There are too few points in the cluster
%						5: The point for this row was rejected after allometric
%							testing (only in case of allometry rejection)
%						6: The point in this row was thrown away in a pass
%							of DBScan (if enabled)
%
%
%		Jonathan Williams
%		jonvw28@gmail.com		         
%		09/01/2018		    
%
%     
% This function uses the DBSCAN MATLAB software as available from yarpiz.com.
% See the directory with DBSCAN to see the license for use of this software

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input Checking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% NOT COMPLETE

if(nargin < 12)
	alo_opt = 0;
end

% Careful not to get trapped in infinite loop
if (min_alo > 1)
	error('Minimum ratio of points for allometry must be <= 1')
end

if (z_olap > 1)
	error('z_olap for height comparison must be less than or equal to 1')
end

if (top_frac > 1)
	error('top_frac for tree top estimation must be less than or equal to 1')
end

if(size(allom_lkup,2)~=2)
	error('allometry table must be Nx2')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Function Body %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get mfile location
scrpt_loc = fileparts(which(mfilename));

% Load dependencies
addpath(genpath([scrpt_loc '\..\..\utils']));

% Number of Trees
trees = sort(unique(LiDAR(:,class_col)));
n_tree = size(trees,1);

% Create Column for acceptance or not
temp_col = zeros(size(LiDAR,1),1);
outFilt = [LiDAR temp_col];
clear temp_col

% For future reference - what's the tallest height in allom lookup
max_allom = max(allom_lkup(:,1));


%%%%%%%%%% First merge possible trees based on the allometry area %%%%%%%%%%%%%%

% get a list of all cluster 'tops' and their heights
early_tops = zeros(n_tree,3);
for a=1:n_tree
	class_a = trees(a);
	cur_tree = outFilt(outFilt(:,class_col)==class_a,:);
	
	% Get centroid of points within top_frac of top height - assumed centre
	[top_ht,~] = max(cur_tree(:,chm_z_col));
	top_pts = cur_tree(cur_tree(:,chm_z_col)>=top_frac*top_ht,1:2);
	cur_cen = mean(top_pts,1);
	early_tops(a,:) = [cur_cen top_ht];
end
clear top_ht mx_idx cur_tree a class_a top_pts

% Get list of trees from tallest down
[~,early_tops_ord] = sort(early_tops(:,3),'descend');

% Re-order group ids for ease later
early_tops = early_tops(early_tops_ord,:);
trees = trees(early_tops_ord);
clear early_tops_ord

% Now take allometry radius for each cluster and test all other clusters for
% how much they overlap this and merge if appropriate

% Testing for various overlap issues in x-y space in a matrix - here a vector 
% reset on each loop
% i,j = 0 means no issue
% i,j = 1 means cluster j has top within radius of top of cluster i
% i,j = 2 means cluster j is more than set fraction within radius of cluster i
% i,j = 3 means both of the above

top_tests = zeros(n_tree,1);

% Log which trees have already been merged so can be ignored
skip_idx = zeros(n_tree,1);

% Take each tree - starting from tallest and merge into it any appropriate
% smaller trees, then repeat for all but the smallest tree (as not necessary)
% Thus big trees are hungry, with priority going to the tallest
for b = 1:(n_tree-1)

	% Pick out the tallest tree left
	class_b = trees(b);
	b_top = early_tops(b,:);
	
	% Skip trees already merged
	if (skip_idx(b,1)~=0)
		class_b = skip_idx(b,1);
	end
	
	% Pull out the points
	tree_b = outFilt(outFilt(:,class_col)==class_b,:);
	
	% Lookup allometry radius
	mod_z = round(b_top(3));
	mod_z = min([mod_z max_allom]);
	alo_b = allom_lkup(find(allom_lkup(:,1)==mod_z),2);
	
	% Find tops in here
	alo_tops = rangesearch(early_tops(:,1:2),b_top(:,1:2),alo_b);
	alo_tops = alo_tops{1};

	% mark the tests vector with a 1 for these
	for t = 1:size(alo_tops,2)
		if(alo_tops(1,t)>b)
			top_tests(alo_tops(1,t),1) = 1;
		end
	end
	top_tests(b,1)=0; % As this will have been made to 1
	clear t alo_tops mod_z
	
	% For all the other trees
	for c = (b+1):n_tree
		
		% Skip trees which have been merged
		if (skip_idx(c,1)~=0)
			continue
		end
	
		% Get the smaller tree
		class_c = trees(c);
		tree_c = outFilt(outFilt(:,class_col)==class_c,:);
		
		% Find points within the allometry range
		circ_pts = rangesearch(tree_c(:,1:2),b_top(:,1:2),alo_b);
		
		% If more than allowed mark this
		if (size(circ_pts{1},2) > alo_share_max*size(tree_c,1))
			top_tests(c,1) = top_tests(c,1)+2;
		end
		
		% Deal with clusters that are flagged
		
		% If not enough points to be left in any case, merge it in to tree_b
		if (top_tests(c,1)~=0 && size(tree_c,1) < min_pts)
			outFilt(outFilt(:,class_col)==class_c,class_col) = class_b;
			skip_idx(c,1) = class_b;
			continue
		end
		
		% Now merge the tree if relevant, unless it is low enough to leave out
		if (top_tests(c,1)~=0) 
			% Test the z_olap-th percentile of heights of the tall tree compared
			% to the (1-z_olap)-th percentil of heights of the short tree
			order_b_ht = sort(tree_b(:,chm_z_col));
			b_cutoff = order_b_ht(max([1 round(z_olap*size(order_b_ht,1))]));
			order_c_ht = sort(tree_c(:,chm_z_col));
			c_cutoff = order_c_ht(min([size(order_c_ht,1) ...
									  round((1-z_olap)*size(order_c_ht,1))]));
			
			% If the small tree cutoff is higher than the tall tree's one, we 
			% merge them
			if(c_cutoff>b_cutoff)
				outFilt(outFilt(:,class_col)==class_c,class_col) = class_b;
				skip_idx(c,1) = class_b;
				clear order_b_ht b_cutoff order_c_ht c_cutoff
				continue
			end
			clear order_b_ht b_cutoff order_c_ht c_cutoff
		end
	end
	clear c class_c tree_c circ_pts
	top_tests = zeros(n_tree,1);
end
clear b class_b alo_b b_top top_tests

% Remove merged clusters from list of trees
n_reject = sum(skip_idx~=0);
n_tree = n_tree - n_reject;
trees(skip_idx~=0) = [];

clear skip_idx n_reject

%%%%%%%%%%% Method to deal with allometrically rejected trees set now %%%%%%%%%%

if(alo_opt == 0) 
% Apply recursive binary hierarchical clustering for clusters which have too
% big a proportion of points outside of the allometry radius
	
	% keep track of where we are
	i=1;

	while i<=n_tree
		tree_i = trees(i);
		alo_pts = 0;
		
		% subest current tree
		cur_tree = outFilt(outFilt(:,class_col)==tree_i,:);
		tr_pts = size(cur_tree,1);
		
		% if we'll reject it on number of points we can skip now anyway
		if(tr_pts < min_pts)
			i = i+1;
			continue
		end
		
		% Get centroid of points within top_frac of top height - 
		% assumed centre could have changed
		[max_z,cur_mx_idx] = max(cur_tree(:,chm_z_col));
		top_pts = cur_tree(cur_tree(:,chm_z_col)>=top_frac*max_z,1:2);
		cur_cen = mean(top_pts,1);
		
		% Lookup allometry diameter
		max_z = round(max_z);
		max_z = min([max_z max_allom]);
		alo_d = allom_lkup(find(allom_lkup(:,1)==max_z),2);
		
		% Test if cluster satisfies this up to a percentage
		alo_idx = rangesearch(cur_tree(:,1:2),cur_cen,alo_d);
		alo_pts = size(alo_idx{1},2);
		
		% used later to escape if needed
		break_test = 0;
		
		% If not apply iterative hierarchical clustering
		while (alo_pts<=tr_pts*min_alo && break_test == 0)		
			% cluster
			hc_l = linkage(cur_tree(:,[1,2,raw_z_col]),'ward','euclidean',...
							'savememory','on');
			hc_out = cluster(hc_l,'maxclust',2);
			
			% Keep cluster featuring max point - create new tree for other
			keep_c_id = hc_out(cur_mx_idx); 
			new_clust = zeros(size(hc_out,1),1);
			new_clust(hc_out~=keep_c_id) = max(trees)+1;
			new_clust(hc_out==keep_c_id) = tree_i;
			
			% re-label clusters
			outFilt(outFilt(:,class_col)==tree_i,class_col)=new_clust;
			trees = [trees; max(trees)+1];
			n_tree = n_tree + 1;
			
			% repeat the pre-loop work
			% Subset tree of interest
			cur_tree = outFilt(outFilt(:,class_col)==tree_i,:);
			tr_pts = size(cur_tree,1);
			
			% if we'll reject it on number of points we can skip now anyway
			if(tr_pts < min_pts)
				break_test = 1;
				continue
			end
			
			% Pull out top of the latest version of this tree
			[max_z,cur_mx_idx] = max(cur_tree(:,chm_z_col));
			top_pts = cur_tree(cur_tree(:,chm_z_col)>=top_frac*max_z,1:2);
			cur_cen = mean(top_pts,1);
			
			% Look up radius		
			max_z = round(max_z);
			max_z = min([max_z max_allom]);
			alo_d = allom_lkup(find(allom_lkup(:,1)==max_z),2);
			
			% Check which points sit in here
			alo_idx = rangesearch(cur_tree(:,1:2),cur_cen,alo_d);
			alo_pts = size(alo_idx{1},2);
		end
		i = i+1;
	end
	
elseif (alo_opt == 1)
% Here we still apply recursive binary hierarchical clustering when a cluster
% is allometrically rejected. Only we outright reject the remaining points, 
% rather than treating them as tree of their own

	for i=1:n_tree
		tree_i = trees(i);
		alo_pts = 0;
		
		% subest current tree
		cur_tree = outFilt(outFilt(:,class_col)==tree_i,:);
		tr_pts = size(cur_tree,1);
		
		% if we'll reject it on number of points we can skip now anyway
		if(tr_pts < min_pts)
			continue
		end
		
		% Get centroid of points within top_frac of top height -
		% could have changed from before
		[max_z,cur_mx_idx] = max(cur_tree(:,chm_z_col));
		top_pts = cur_tree(cur_tree(:,chm_z_col)>=top_frac*max_z,1:2);
		cur_cen = mean(top_pts,1);
		
		% Lookup allometry diameter
		max_z = round(max_z);
		max_z = min([max_z max_allom]);
		alo_d = allom_lkup(find(allom_lkup(:,1)==max_z),2);
		
		% Test if cluster satisfies this up to a percentage
		alo_idx = rangesearch(cur_tree(:,1:2),cur_cen,alo_d);
		alo_pts = size(alo_idx{1},2);
		
		% used later to escape if needed
		break_test = 0;
		
		% If not apply iterative hierarchical clustering
		while (alo_pts<=tr_pts*min_alo && break_test == 0)		
			% cluster
			hc_d = pdist(cur_tree(:,[1,2,raw_z_col]));
			hc_l = linkage(hc_d);
			hc_out = cluster(hc_l,'maxclust',2);
			
			% Keep cluster featuring max point - create new tree for other
			keep_c_id = hc_out(cur_mx_idx); 
			new_clust = zeros(size(hc_out,1),1);
			new_clust(hc_out~=keep_c_id) = 0; % marked to be rejected
			new_clust(hc_out==keep_c_id) = tree_i;
			
			% update our quality label column
			clust_flag = zeros(size(hc_out,1),1);
			clust_flag(hc_out~=keep_c_id) = 5; % as allom rejected
			clust_flag(hc_out==keep_c_id) = 0;
			
			% re-label clusters
			outFilt(outFilt(:,class_col)==tree_i,size(outFilt,2))=clust_flag;
			outFilt(outFilt(:,class_col)==tree_i,class_col)=new_clust;
			
			% repeat the pre-loop work
			% Subset tree of interest
			cur_tree = outFilt(outFilt(:,class_col)==tree_i,:);
			tr_pts = size(cur_tree,1);
			
			% if we'll reject it on number of points we can skip now anyway
			if(tr_pts < min_pts)
				break_test = 1;
				continue
			end
			
			% Pull out top of the latest version of this tree
			[max_z,cur_mx_idx] = max(cur_tree(:,chm_z_col));
			top_pts = cur_tree(cur_tree(:,chm_z_col)>=top_frac*max_z,1:2);
			cur_cen = mean(top_pts,1);
			
			% Look up radius		
			max_z = round(max_z);
			max_z = min([max_z max_allom]);
			alo_d = allom_lkup(find(allom_lkup(:,1)==max_z),2);
			
			% Check which point sit in here
			alo_idx = rangesearch(cur_tree(:,1:2),cur_cen,alo_d);
			alo_pts = size(alo_idx{1},2);
		end
	end
	% Tag all allometrically rejected points as being so
	outFilt(outFilt(:,class_col)==0,end)=5;	
end

%%%%%%%%%%%%%%%%%%%%%%%% If enabled, apply DBScan %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if(db_eps~=0)
	for(db_i=1:n_tree)
		% Collect points in tree
		tree_i = trees(db_i);
		this_tree = outFilt(outFilt(:,class_col)==tree_i,:);
		[~,tp_idx] = max(this_tree(:,chm_z_col));
		
		% Apply DBScan in 3D projection
		db_out = DBSCAN(this_tree(:,1:3),db_eps,db_pts); 
			
		% keep cluster with tallest point
		keep_clust = db_out(tp_idx); 
		new_clust = zeros(size(db_out,1),1);
		
		% Relabel clusters
		new_clust(db_out==keep_clust) = tree_i;
		new_clust(db_out~=keep_clust) = 0;
		
		% update our quality label column
		clust_flag = zeros(size(db_out,1),1);
		clust_flag(db_out~=keep_clust) = 6; % as DBScan rejected
		clust_flag(db_out==keep_clust) = 0;
		
		% Add to our data
		outFilt(outFilt(:,class_col)==tree_i,end) = clust_flag;
		outFilt(outFilt(:,class_col)==tree_i,class_col) = new_clust;
		clear keep_clust new_clust db_out clust_flag
		
	end
end

%%%%%%%% Test for number of points  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get a list of all cluster outlines
for d=1:n_tree
	class_d = trees(d);
	cur_tree = outFilt(outFilt(:,class_col)==class_d,:);
	
	% Check if we've already rejected it and skip
	if (sum(cur_tree(:,end))~=0)
		continue
	end
		
	% Check there are enough points - don't bother otherwise and mark this
	if (size(unique(cur_tree(:,1:2),'rows'),1) < min_pts)
		outFilt(outFilt(:,class_col)==class_d,end)=2;
		continue
	end
	% Now tree is good so mark it
	outFilt(outFilt(:,class_col)==class_d,end) = 1;
end
clear cur_tree d class_d

end