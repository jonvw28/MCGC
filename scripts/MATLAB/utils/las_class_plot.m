function [ ] = las_class_plot( LiDAR, class_col, prior )
%LAS_CLASS_PLOT Function to plot 3d las data coloured by segmentation
% classes
%
% Syntax
%
%		las_class_plot(LiDAR,class_col)
%		las_class_plot(LiDAR,class_col,prior)
%		
% 		las_class_plot(LiDAR,class_col) plots the points in LiDAR columns
%		1 to 3 in x,y,z, coloured by the groupings given by classes numbered in
%		the column given by class_col
%
%		las_class_plot(LiDAR,class_col,prior) addtitonally plots the prior
%		points given in the first three columns of prior
%
% Method
%
%		Points from each cluster in class_col are all plotted in x,y,z space
%		with a randomly assigend colour. If a prior is included then this is 
%		plotted as stars.
%      
% Inputs
%
%		LiDAR: 		A point cloud where the first three columns are the x,y and 
%					z co-ordinates of each point
%
%		class_col: 	the column of LiDAR where class information is included
%
%		prior: 		An optional matrix of prior points, again x,y,z in the first
%			   		three columns.
%
% Outputs: 
%
%		This simply generates a plot
%
%		Jonathan Williams
%		jonvw28@gmail.com			         
%		09/03/2017

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input Checking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		     

if(ismatrix(LiDAR)~= 1 || isnumeric(LiDAR(:,1:3))~=1)
	error(['LiDAR must be a matrix, where the first three columns are ' ...
			'numeric values of x,y,z co-ordinates for the points'])
elseif(isnumeric(class_col)~=1 || size(class_col,1)~=1 || size(class_col,2)~=1)
	error('class_col must be a single number')
elseif(nargin > 2 && (ismatrix(prior)~= 1 || isnumeric(prior(:,1:3))~=1))
	error(['prior must be a matrix, where the first three columns are ' ...
			'numeric values of x,y,z co-ordinates for the points'])
end

% Check class_col makes sense
if(class_col <=0 || class_col > size(LiDAR,2))
	error('class_col must refer to a column in LiDAR')
elseif(isnumeric(LiDAR(:,class_col))~=1)
	error('class column in LiDAR must have classes labelled by numbers')
end
			
			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Function Body %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hold on
% Number of Trees
classes = unique(LiDAR(:,class_col));
N = size(classes,1);
cstring=rand(N,3);

% Plot trees
for i =1:N
    id = find(LiDAR(:,class_col)==classes(i));
    plot3(LiDAR(id,1),LiDAR(id,2),LiDAR(id,3),'.','color',cstring(mod(i,N)+1,:),'markersize',4);
end

% Plot Prior
if nargin > 2
    plot3(prior(:,1),prior(:,2),prior(:,3),'r*','markersize',8);
end
hold off
end