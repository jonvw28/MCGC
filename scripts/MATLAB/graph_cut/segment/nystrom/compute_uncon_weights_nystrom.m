function [A,B,idx] = compute_uncon_weights_nystrom(data, raw_z_col,...
						chm_z_col, no_samp, sigxy, sigz, allom, r_rat,...
						h_grad_opt, hgrad_w, z_grad_opt, zgrad_w)
						
%COMPUTE_UNCON_WEIGHTS_NYSTROM
%
% Compute adjacency matrices required for nystrom extension for an undirected
% graph including information on local geometry
%
% Syntax
%
%		[A,B,idx] = compute_uncon_weights_nystrom(data, raw_z_col, chm_z_col,...
%												no_samp, sigxy, sigz, allom,...
%												r_rat, h_grad_opt, hgrad_w,...
%												z_grad_opt, zgrad_w)
%		
%		This returns the matrices A and B of pairwise similarities between 
%		sample points, and between all other points and sample points 
%		respectively. idx gives the indices of the sample points in the original
%		data
%
% Method
%
%	A number of points set by no_samp are selected. For these a full adjacency 
%	matrix is computed and returned as A. For the remaining points a link matrix
%	B is computed which returns the adjacency weight between each remaining
%	point and each point in the sample used for A. The indices of the
%	selected data points in the original data are returned in idx.
%
%	Pairwise adjacency weights are calculated as follows:
%
%	First, a basic similarity based on the distance between each pair of points
%	is computed, $w_{ij}^{base}$. This is separated into a horizontal and a
%	vertical component. Each component has its importance controlled by 
%	separate parameters (sigxy and sigz respectively). Here x_i, y_i and z_i-z_j
%	are the coordinates in the raw point cloud of point i.
%
%	w_{ij}^{base} = exp(- \frac{||(x,y)_{i}-(x,y)_{j}||^{2}_{2}}{sigxy ^2}) *
%			 		exp(-\frac{(z_{i}-z_{j})^2}{sigz ^2}) 
%
%	Next a centroid vector is computed for each point. Here the allometric
%	look-up table in allom is used to find an allometric radius for each point
%	based upon the point's height above ground. A centroid of all points in the
%	dataset within a sphere of this radius (multiplied by r_rat) centred at the
%	point in question is computed (using the raw z values). The vector from the
%	point in question to the centroid is then stored in a vector Delta_{i} (for 
%	the i-th point). This is split into the horziontal (Delta_{H_{i}}) and 
%	vertical (Delta_{Z_{i}}) components. Comparisons between the respective
%	centroid vectors then modify the weight between 2 points as set out in the
%	following.
%
%   First the horizontal components (Delta_{H_{i}} and Delta_{H_{j}}) are
%	compared and the weight is reduced if the angle between these is greater
%	than pi/2. The relative importance of this step is set by hgrad_w. There are
%	4 options set by h_grad_opt which decide the way this reduction is realised
%	as set out below. Dh is the horizontal euclidean distance between points i,j
%	and Kh a scale factor. This is computed by finding the highest point above
%	ground in data and looking up the allometric radius in allom for this height
%
%		1: 'Uniform adjustment'
%				mutliply w_{ij} by exp(-hgrad_w)
%		2: 'Inverse separation weighting'
%				multiply w_{ij} by exp(-hgrad_w * Kh/Dh)
%		3: 'Weight by Delta difference'
%				multiply w_{ij} by exp(-hgrad_w * ||Delta_{H_i}-Delta_{H_j}||_2)
%		4: 'Composite weighting'
%				as 2 with exp(-hgrad_w * Kh/Dh * ||Delta_{H_i}-Delta_{H_j}||_2)
%
%	Similarly the weighting is adjusted based on the vertical components of
%	Delta_i and Delta_j. Here an adjustment is only maded if the point with a
%	taller raw height has postive Delta_Z and the lower point has negative
%	Delta_Z (indicating divergence). In this case there are 4 options set by
%	z_grad_opt with weigthing parameter zgrad_w. Analagously the effects are
%	shown below. Here Dz is the raw height difference between points i,j and
%	Kz is a constant computed from the height of tallest above ground point in
%	data then divided by 2.
%
%		1: Uniform adjustment
%				mutliply w_{ij} by exp(-hgrad_z)
%		2: 'Inverse separation weighting'
%				multiply w_{ij} by exp(-hgrad_z * Kz/Dz)
%		3: 'Weight by Delta difference'
%				multiply w_{ij} by exp(-hgrad_z * |Delta_{Z_i}-Delta_{Z_j}|)
%		4: 'Composite weighting'
%				as 2 with exp(-hgrad_z * Kz/Dz * |Delta_{Z_i}-Delta_{Z_j}|)
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
%
% Outputs: 
%
%		A: 	Pair-wise similarity matrix for selected points
%
%		B: 	Pair-wise similarity matrix for remaining points to selected points
%
%		idx: Indices of the selected data points in the original data
%
%
% Dependancy Tree
%
%		This function is standalone. It is required directly by nystrom_ext
%
%			nystrom_ext
%			 ->
%				compute_uncon_weights_nystrom
%
%
%		Jonathan Williams
%		jonvw28@gmail.com			         
%		09/01/2019	         


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
elseif(isnumeric(h_grad_opt)~=1 || size(h_grad_opt,1) ~= 1 || size(h_grad_opt,2) ~= 1)
	error('h_grad_opt must be a single number')
elseif(isnumeric(z_grad_opt)~=1 || size(z_grad_opt,1) ~= 1 || size(z_grad_opt,2) ~= 1)
	error('z_grad_opt must be a single number')
elseif(isnumeric(hgrad_w)~=1 || size(hgrad_w,1) ~= 1 || size(hgrad_w,2) ~= 1)
	error('hgrad_w must be a single number')
elseif(isnumeric(zgrad_w)~=1 || size(zgrad_w,1) ~= 1 || size(zgrad_w,2) ~= 1)
	error('zgrad_w must be a single number')
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
elseif(sigxy == 0)
	error('sigxy can not be zero')
elseif(sigz == 0)
	error('sigz can not be zero')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Function Body %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of points
R = size(data,1);

% For future reference - what's the tallest height in allom lookup
max_allom = max(allom(:,1));

% For future reference - what's the tallest abovegorund height
max_agh = max(data(:,chm_z_col));
max_agh_r = allom(allom(:,1)==min([round(max_agh) max_allom]),2);
max_agh = max_agh/2;

% Compute V - vector to local centroid
V = zeros(R,3);

% generate lookup clamped Heights
z = round(data(:,chm_z_col));
z = min([z repmat(max_allom,R,1)],[],2);

% get vectors to centroids, going through allowed heights
for j = 0:max_allom
	% pull out points in this band
	idx = z==j;
	ind=find(idx);
	r_alo = allom(j+1,2)*r_rat;
	if(sum(idx)==0) % skip if none
		continue
	end
	nbrs = rangesearch(data(:,[1 2 raw_z_col]),data(idx,[1 2 raw_z_col]),r_alo);
	for n = 1:length(nbrs)
		cent = sum(data(nbrs{n},[1 2 raw_z_col]),1)/size(nbrs{n},2);
		V(ind(n),:) = cent - data(ind(n),[1 2 raw_z_col]);
	end
end


% Take the samples
no_samp = round(no_samp);
% Deal with case that sample size is more than number of data points
if (no_samp >= R)
	sub_data = data;
	A = zeros(R);
	no_samp = R;
	idx = 1:R;
else
	[sub_data, idx] = datasample(data,no_samp,1,'Replace',false);
	A = zeros(no_samp);
end

% Complete A, as upper triangular
for a = 1:(no_samp-1)
	% Difference in x,y
	xy = sub_data(a,1:2); 
	deltaxy = repmat(xy,no_samp-a,1) - sub_data((a+1):no_samp,1:2);
	xyterm = sum(deltaxy.^2,2);
	
	% Difference in z
	z = sub_data(a,raw_z_col); 
	deltaz = repmat(z,no_samp-a,1) - sub_data((a+1):no_samp,raw_z_col);
	zterm = deltaz.^2;
	
	% Local density gradient
	V_a = V(idx(a),:);
	V_all = V(idx((a+1):no_samp),:);
	deltaV = repmat(V_a,no_samp-a,1) - V_all;
	
	% Horizontal
	chi = dot(repmat(V_a(:,1:2),no_samp-a,1),V_all(:,1:2),2);
	if(h_grad_opt == 1) % simple >90deg mask
		chi = (chi < 0);
	elseif(h_grad_opt == 2) % weight mask by inverse distance between points
		chi = (chi < 0);
		h_dist = sqrt(xyterm);
		chi = chi*max_agh_r ./ max([h_dist repmat(1e-15,no_samp-a,1)],[],2);
	elseif(h_grad_opt == 3) % weight mask by magnitude of gradient difference
		chi = (chi < 0);
		gh_dist = sqrt(sum(deltaV(:,1:2).^2,2));
		chi = chi .* gh_dist;
	elseif(h_grad_opt == 4) %  weight mask by inverse distance between points and grad magnitude difference
		chi = (chi < 0);
		h_dist = sqrt(xyterm);
		gh_dist = sqrt(sum(deltaV(:,1:2).^2,2));
		chi = chi * max_agh_r ./max([h_dist repmat(1e-15,no_samp-a,1)],[],2);
		chi = chi .* gh_dist;
	elseif(h_grad_opt == 5) % scale dot product to vary from 0 (0deg) to 2 (180deg)
		chi = chi/ max([sqrt(sum(V_a(1:2).^2)) 1e-15]);
		chi = chi ./ max([sqrt(sum(V_all(:,1:2).^2,2)) repmat(1e-15,no_samp-a,1)],[],2);
		chi = 1 - chi;
	else
		error('h_grad_opt must be 1, 2, 3, 4 or 5 - see help for details')
	end
	hgrad_term = chi*hgrad_w;
	
	% Vertical
	chiz = repmat(sign(V_a(3)),no_samp-a,1).*sign(V_all(:,3));
	% Only opposite directions
	chiz = (chiz < 0);
	% which is taller
	a_tall = (deltaz > 0);
	% Get the sign of the taller point
	tall_sign = zeros(no_samp-a,1);
	tall_sign(a_tall) = sign(V_a(3));
	tall_sign(~a_tall) = V_all(~a_tall,3);
	% now get mask of this (ie different directions)
	chiz = chiz.*tall_sign;
	% And only select opposite direction pairs which diverge
	chiz = (chiz > 0);
	if(z_grad_opt == 1) % simple mask
		
	elseif(z_grad_opt == 2) % weight mask by inverse distance between points
		z_dist = abs(zterm);
		chiz = chiz*max_agh ./ max([z_dist repmat(1e-15,no_samp-a,1)],[],2);
	elseif(z_grad_opt == 3) % weight mask by magnitude of gradient difference
		gz_dist = abs(deltaV(:,3));
		chiz = chiz .* gz_dist;
	elseif(z_grad_opt == 4) %  weight mask by inverse distance between points and grad magnitude difference
		z_dist = abs(zterm);
		gz_dist = abs(deltaV(:,3));
		chiz = chiz*max_agh ./ max([z_dist repmat(1e-15,no_samp-a,1)],[],2);
		chiz = chiz .* gz_dist;
	else
		error('z_grad_opt must be 1, 2, 3 or 4 - see help for details')
	end
	zgrad_term = chiz*zgrad_w;
	
	% Compute exponential weight
	w = exp(-xyterm/sigxy^2) .* exp(-zterm/sigz^2) 	.* exp(-hgrad_term) .*exp(-zgrad_term); 
	A(a,(a+1):no_samp) = w;
end

% Complete the matrix
A = A + A';

% compute B if A is not complete dataset

if (no_samp ~= R)
	% get remaining points
	rem_data = data;
	rem_data(idx,:) = [];
	b_idx = 1:R;
	b_idx(idx)=[];
	B = zeros(no_samp,R-no_samp);

	% Avoid repeated indexing
	Axy = sub_data(:,1:2);
	Az = sub_data(:,raw_z_col);
	AV = V(idx,:);
	 
	% Complete B row-wise
	for b = 1:(R-no_samp)
		% Difference in x,y
		xy = rem_data(b,1:2); 
		deltaxy = repmat(xy,no_samp,1) - Axy;
		xyterm = sum(deltaxy.^2,2);
		
		% Difference in z
		z = rem_data(b,raw_z_col); 
		deltaz = repmat(z,no_samp,1) - Az;
		zterm = deltaz.^2;
		
		% Local density gradient
		V_b = V(b_idx(b),:);
		deltaV = repmat(V_b,no_samp,1) - AV;
		
		% Horizontal
		chi = dot(repmat(V_b(:,1:2),no_samp,1),AV(:,1:2),2);
		if(h_grad_opt == 1) % simple >90deg mask
			chi = (chi < 0);
		elseif(h_grad_opt == 2) % weight mask by inverse distance between points
			chi = (chi < 0);
			h_dist = sqrt(xyterm);
			chi = chi*max_agh_r ./ max([h_dist repmat(1e-15,no_samp,1)],[],2);
		elseif(h_grad_opt == 3) % weight mask by magnitude of gradient difference
			chi = (chi < 0);
			gh_dist = sqrt(sum(deltaV(:,1:2).^2,2));
			chi = chi .* gh_dist;
		elseif(h_grad_opt == 4) %  weight mask by inverse distance between points and grad magnitude difference
			chi = (chi < 0);
			h_dist = sqrt(xyterm);
			gh_dist = sqrt(sum(deltaV(:,1:2).^2,2));
			chi = chi * max_agh_r ./max([h_dist repmat(1e-15,no_samp,1)],[],2);
			chi = chi .* gh_dist;
		elseif(h_grad_opt == 5) % scale dot product to vary from 0 (0deg) to 2 (180deg)
			chi = chi/ max([sqrt(sum(V_b(1:2).^2)) 1e-15]);
			chi = chi ./ max([sqrt(sum(AV(:,1:2).^2,2)) repmat(1e-15,no_samp,1)],[],2);
			chi = 1 - chi;
		else
			error('h_grad_opt must be 1, 2, 3, 4 or 5 - see help for details')
		end
		hgrad_term = chi*hgrad_w;
		
		% Vertical
		chiz = repmat(sign(V_b(3)),no_samp,1).*sign(AV(:,3));
		% Only opposite directions
		chiz = (chiz < 0);
		% which is taller
		a_tall = (deltaz > 0);
		% Get the sign of the taller point
		tall_sign = zeros(no_samp,1);
		tall_sign(a_tall) = sign(V_b(3));
		tall_sign(~a_tall) = AV(~a_tall,3);
		% now get mask of this (ie different directions)
		chiz = chiz.*tall_sign;
		% And only select opposite direction pairs which diverge
		chiz = (chiz > 0);
		if(z_grad_opt == 1) % simple mask
			
		elseif(z_grad_opt == 2) % weight mask by inverse distance between points
			z_dist = abs(zterm);
			chiz = chiz*max_agh ./ max([z_dist repmat(1e-15,no_samp,1)],[],2);
		elseif(z_grad_opt == 3) % weight mask by magnitude of gradient difference
			gz_dist = abs(deltaV(:,3));
			chiz = chiz .* gz_dist;
		elseif(z_grad_opt == 4) %  weight mask by inverse distance between points and grad magnitude difference
			z_dist = abs(zterm);
			gz_dist = abs(deltaV(:,3));
			chiz = chiz*max_agh ./ max([z_dist repmat(1e-15,no_samp,1)],[],2);
			chiz = chiz .* gz_dist;
		else
			error('z_grad_opt must be 1, 2, 3 or 4 - see help for details')
		end
		zgrad_term = chiz*zgrad_w;
		
	w = exp(-xyterm/sigxy^2) .* exp(-zterm/sigz^2) 	.* exp(-hgrad_term) .*exp(-zgrad_term); 
		B(:,b)=w;
	end
else
	B = 0;
end
end