function [E,L] = nystrom_ext(data, raw_z_col, chm_z_col, no_samp, sigxy,...
						sigz, allom, r_rat, h_grad_opt, hgrad_w, z_grad_opt,...
						zgrad_w, nVecs)

%NYSTROM_EXT
% Apply Nystrom extension to data, computing Wij and finding the first nVecs
% eigenvectors
%
% Syntax
%
% 	[E,L] = nystrom_ext(data,raw_z_col,chm_z_col,no_samp,sigxy,...
%						sigz,allom,r_rat,h_grad_opt,hgrad_w,z_grad_opt,...
%						zgrad_w,nVecs)
%
%	This takes the points in data, computes matrices of Wij for the
%	nystrom extension (according to the definitions for MCGC) and then returns
%	the first nVecs normalised eigenvectors of L_sym as used in a normalised
%	graph cut, returned in E, and also produces the mathcing set of eigenvalues
%	for the symmetric Laplacian, returned in ascending order in L.
%
%
% Method
%
%	The points in data are used to compute approximations to the eigenvectors
%	of L_sym as defined in the framework of normalised graph cuts. Here the
%	relevant vectors that are kept are 2 to nVecs+1 when ordered by increasing
%	eigenvalue.
%
%	Details of the theory of normalised graph cuts are summarised in [1] and the
%	normalised graph cut approach used here is originally introduced in [2].
%	Details of the nystrom extension applied to normalised graph cuts are
%	available in [3] and [4]. The methods to compute the weights matrix are
%	detailed in the documentation for the function compute_uncon_weights_nystrom
%
%	The only addtionaly parameter to this function is nVecs which determines how
%	many eigenvectors should be returned. Further details on the roles of the
%	remaining parameters can also be found in the documentation for
%	compute_uncon_weights_nystrom.
%
%
% Inputs:
%
%		data: 		A matrix of points where each row is a node and the first
%					two columns are x,y coordinates of this node. this must
%					include raw and aboveground heights of each poitn as well
%
%		raw_z_col:	The column containing the raw z value for each point
%
%		chm_z_col:	The column which contains the value of height above the 
%					ground for each point
%		
%		no_samp:	Number of data points to sample for nystrom extension
%
%		sigxy: 		Parameter for significance of planimetric distance in 
%					linkages
%
%		sigz: 		Parameter for significance of vertical distance in linkages
%
%       allom:		Allometric lookup table for centroid computation. Must have
%					first 2 columns as height (rounded to nearest metre) and
%					allometric radius respectively.
%
%		r_rat:		Fraction of lookup radius from allom to use for centroid
%					computation (recommended to use 0.5 or 1)
%
%		h_grad_opt:	Sets which option to use for Delta_H comparison:
%						1: 'Uniform adjustment'
%						2: 'Inverse separation weighting'
%						3: 'Weight by Delta difference'
%						4: 'Composite weighting'
%
%		hgrad_w:	Parameter to set significance of Delta_H term
%
%		z_grad_opt:	Sets which option to use for Delta_Z comparison:
%						1: 'Uniform adjustment'
%						2: 'Inverse separation weighting'
%						3: 'Weight by Delta difference'
%						4: 'Composite weighting'
%
%		zgrad_w:	Parameter to set significance of Delta_Z term
%
%		nVecs:		Number or eigenvectors and eigenvalues to return
%
%
%	Outputs:
%
%		E:	First nVecs normalised eignenvectors of L_sym (ignoring first
%			eigenvector) in columns of a matrix
%
%		L:	Diagonal maxtrix of first nVecs eigenvalues of L_sym (ignoring
%			first eigenvalue)
%
%	References:
%
%       [1]	A Tutorial on Spectral Clustering, U von Luxburg, 
%       	Statistics and Computing, 17 (4), 2007
%
%		[2] On Spectral Clustering: Analysis and an algorithm, A Y ng et al.,
%			Advances in neural information processing systems. 2002.
%
%		[3] Spectral Grouping Using the Nystrom Method, C Fowlkes et al.,
%			IEEE TRANSACTIONS ON PATTERN ANALYSIS AND MACHINE INTELLIGENCE,
%			VOL. 26, NO. 2, FEBRUARY 2004
%
%		[4]	Graph Clustering, Variational Image Segmentation Methods and Hough
%			Transform Scale Detection for Object Measurement in Images, L_ord
%			Calatroni et al., Journal of Mathematical Imaging and Vision, 2017
%
%
% Dependancy Tree
%
%		This function requires compute_uncon_weights_nystrom. This function is
%		required by mcgc_cut.
%
%		mcgc_cut
%		 ->
%			nystrom_ext
%			 ->
%				compute_uncon_weights_nystrom
%
%
%		Jonathan Williams
%		jonvw28@gmail.com			         
%		09/01/2019	  
%
%		Heavily influenced by the example code in [3] and adapted from the code
%		used originally in [4]



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input Checking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Type checking
if(ismatrix(data)~= 1 || isnumeric(data(:,[1,2,raw_z_col,chm_z_col]))~=1)
	error(['data must be a matrix, where the first three columns are ' ...
			'numeric values of x,y,z co-ordinates for the points'])
elseif(isnumeric(no_samp)~=1 || size(no_samp,1) ~= 1 || size(no_samp,2) ~= 1)
	error('sigz must be a single number')
elseif(ismatrix(allom)~= 1 || isnumeric(allom(:,1:2))~=1)
	error(['allom must be a matrix, where the first two columns are ' ...
			'numeric values of H and R for lookup'])
elseif(isnumeric(sigxy)~=1 || size(sigxy,1) ~= 1 || size(sigxy,2) ~= 1)
	error('sigxy must be a single number')
elseif(isnumeric(sigz)~=1 || size(sigz,1) ~= 1 || size(sigz,2) ~= 1)
	error('sigz must be a single number')
elseif(isnumeric(r_rat)~=1 || size(r_rat,1) ~= 1 || size(r_rat,2) ~= 1)
	error('r_rat must be a single number')
elseif(isnumeric(h_grad_opt)~=1 || size(h_grad_opt,1) ~= 1 ...
		|| size(h_grad_opt,2) ~= 1)
	error('h_grad_opt must be a single number')
elseif(isnumeric(z_grad_opt)~=1 || size(z_grad_opt,1) ~= 1 ...
		|| size(z_grad_opt,2) ~= 1)
	error('z_grad_opt must be a single number')
elseif(isnumeric(hgrad_w)~=1 || size(hgrad_w,1) ~= 1 || size(hgrad_w,2) ~= 1)
	error('hgrad_w must be a single number')
elseif(isnumeric(zgrad_w)~=1 || size(zgrad_w,1) ~= 1 || size(zgrad_w,2) ~= 1)
	error('zgrad_w must be a single number')
elseif(isnumeric(nVecs)~=1 || size(nVecs,1) ~= 1 || size(nVecs,2) ~= 1)
	error('nVecs must be a single number')
end

% Check there are points
if(size(data,1)==0)
	error('data can not be empty')
end

% Check allom is of form expected
if(size(allom,2)~=2)
	error('allom must only consist of two columns, being H and R values')
elseif(size(allom,1)==0)
	error('allom can not be empty')
end

% Check validity of relevant arguments
if (no_samp <= 0)
	error('no_samp must be a postive number')
elseif (nVecs <= 0)
	error('nVecs must be a postive whole number')
elseif(sigxy == 0)
	error('sigxy can not be zero')
elseif(sigz == 0)
	error('sigz can not be zero')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Function Body %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute A B for the chosen methods
[A,B,idx] = compute_uncon_weights_nystrom(data,raw_z_col,chm_z_col,...
											no_samp,sigxy,sigz,allom,...
											r_rat,h_grad_opt,hgrad_w,...
											z_grad_opt,zgrad_w);
%
% Deal with case where A is complete W
%
if (size(A,1) == size(data,1))
	% generate D
	D = sum(A,1);
	D(D(:)<1e-12)=1e-12;
	% Generate L_sym
	D2 = 1./sqrt(D);
	D = diag(D);
	D2 = diag(D2);
	Lsym = D2*(D-A)*D2;
	if not(issymmetric(Lsym))
		Lsym = (Lsym + Lsym')/2;
	end
	% solve for e-vecs
	[U,L,~] = svd(Lsym);
	% order in increasing eigenvalues
	[L,L_ord] = sort(diag(L));
	V(:,:) = U(:,L_ord);
	% Keep only relevant e-vecs
	L = L(2:(nVecs+1));
	L = diag(L);
	E=V(:,2:nVecs+1);
else
	% use Nystrom extension
	%
	n = size(A,1);
	m = size(B,2);
	% Invert A
	[Uinit,Linit,Tinit] = svd(A);
	Ai = Uinit*diag(1./diag(Linit))*Tinit';
	% Row normalisation of W approximation
	d1 = sum([A; B'], 1);
	d2 = sum(B, 1) + sum(B', 1)*Ai*B;
	% Added abs to avoid issues from small negative values
	dhat = sqrt(1./abs([d1 d2]))';
	dhat(dhat(:)<1e-12) = 1e-12;
	A = A.*(dhat(1:n)*dhat(1:n)');
	B = B.*(dhat(1:n)*dhat(n + (1:m))');

	% Invert normalised A
	[Unorm,Lnorm,Tnorm] = svd(A);
	Asi = Tnorm*diag(sqrt(1./diag(Lnorm)))*Unorm';

	% Compute approximations to eigenvectors
	Q = A + Asi*(B*B')*Asi;
	[U,L] = svd(Q);
	V = [A;B']*(Asi*(U*diag(diag(1./sqrt(L)))));
	vecs=V(:,2:nVecs+1);
	
	% Strip L down, compute 1-L and return it (ie move from evals of L to evals
	% of L_sym
	L = diag(L);
	L = L(2:(nVecs+1));
	L = 1 - L;
	L = diag(L);


	% Re-order to match the original data
	E = zeros(n+m,nVecs);
	E(idx,:) = vecs(1:n,:);
	full_idx = 1:(n+m);
	rem_idx = setdiff(full_idx,idx);
	E(rem_idx,:) = vecs(n+(1:m),:);
end

% Normalise the rows as per Ng et al. 2002
normE = sqrt(sum(E.^2,2));
normE(normE==0) = 1;
E = bsxfun(@rdivide,E,normE);
end