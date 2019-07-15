function [idx,X] = v2idx(v,ngroups,method,iterCount)
%V2IDX Function to convert eigenvectors to classes in multiclass graph cut
%
% Syntax
%
%		idx = v2idx(v,ngroups)
%		idx = v2idx(v,ngroups,method)
%		idx = v2idx(v,ngroups,method,iterCount)
%		[idx,X] = v2idx(v,ngroups,method,iterCount)
%		
% 		idx = v2idx(v,ngroups) converts the N x ngroups matrix of eigenvectors
%		in v into a length N vector of classes (1:ngroup) for each point
%		using k-means clustering on the rows of v and return this as idx
%
%		idx = v2idx(v,ngroups,method) achieves the same as above with multiple
%		methods permitted:
%				0: (default) use k-means on rows of v2idx
%				1: First convert entries of v to their sign, then apply k-means
%					to the rows of the resulting matrix
%				2: Assign each point to class based on the column of the 
%					maximum value in its row in v
%
%		idx = v2idx(v,ngroups,method,iterCount) same as above, but can also
%		specify the number of iterations for k-means where relevant
%
%		[idx,X] = v2idx(v,ngroups,method,iterCount) as above, but also returns
%		the matrix used for clustering (either v, or sgn(v), should this be
%		required)
%
% Method
%
%		Convert an N x ngroups matrix, v, of the eigenvectors for the smallest
%		ngroups eigenvectors of a graph adjacency matrix into an indicator
%		vector, idx, where each of the N vertices is assigned to one of ngroups
%		classes. Can also return the form of v eventually used for this if 
%		desired. There are three methods for this conversion:
%
%			0: (default) k-means clustering on the rows of v
%			1: converting entries of v to take only their sign and then applying
%				k-means to the rows
%			2: assigning each vertex based on column containing the maximum 
%				value in its row of v - for hard constraint
%
%		Where k-means is used (0 or 1) it is possible to also specify the number
%		of iterations through iterCount
%      
% Inputs:
%
%		v:			N x ngroups matrix of eigenvectors of graph adjacency matrix
%					with smallest eigenvalues
%
%		ngroups:	number of groups to cluster (needs to match number of 
%					columns of v)
%
%		method:		0 for k-means on rows, 1 for k-means on signs, 2 for maximum 
%					element (see emthod for details)
%
%		iterCount:	number of replicates for k-means where relevant
%
% Outputs: 
%
%		idx:	length N vector with class for each of the N vertices
%
%		X:		eventual form of v used for clustering
%
%		Jonathan Williams
%		jonvw28@gmail.com			         
%		11/07/2017		

% Fill in any unstated arguments
if nargin < 3
    method =0;
end

if nargin <4
    iterCount =20;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input Checking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% type checking
if(ismatrix(v)~= 1 || isnumeric(v)~=1)
	error('v must be a numeric matrix')
elseif(isnumeric(ngroups)~=1 || size(ngroups,1) ~= 1 || size(ngroups,2) ~= 1)
	error('ngroups must be a single number')
elseif(isnumeric(method)~=1 || size(method,1) ~= 1 || size(method,2) ~= 1)
	error('method must be a single number')
elseif(isnumeric(iterCount)~=1 || size(iterCount,1) ~= 1 || ...
		size(iterCount,2) ~= 1)
	error('iterCount must be a single number')
end

% Check feasibility

% Deal with possible non-integer ngroups
ngroups = round(ngroups);
% test feasible input
if(ngroups <= 0)
	error('ngroups must be a positive integer')
end

% Check method selected is valid
if(ismember(method,[0 1 2])==0)
	error('method must be one of 0, 1 or 2')
end

% Deal with possible non-integer maxit
iterCount = round(iterCount);
% test feasible input
if(iterCount <= 0)
	error('iterCount must be a positive integer')
end

% Is ngroup consistent
if(size(v,2)~=ngroups)
	error('ngroup must match the number of columns of v')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Function Body %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opts = statset('Display','off');

if(method == 0)    
    X = v;
    [idx] = kmeans(v(:,1:end),ngroups,'EmptyAction','drop','Replicates',...
        iterCount,'options',opts);
elseif(method==1)
    X = sign(v);
    [idx] = kmeans(X,ngroups,'EmptyAction','drop','Replicates',...
        iterCount,'options',opts);
elseif(method==2)
    X = v;
    [~,idx] = max(v,[],2);
   
end
