addpath(genpath('../../MCGC'));

mkdir single_layer
mkdir double_layer

file = 'fake_trees_2_tight';

data_folder ='.';
out_folder_1 = 'single_layer';
out_folder_2 = 'double_layer';
pr_allom = 'indomalaya_tropfor_lut.csv';
allom = csvread('indomalaya_tropfor_HR_95.csv');
sig_xy = 4;
sig_z = 2;
hgrad_w = 0.2;
zgrad_w = 0.2;
db_eps=1;
db_pts=5;
min_pts=100;
subSamp = 2;


% test mcgc_base_options with subsampling
rng(123);
try
	Seg1 = mcgc_base_options(data_folder,out_folder_1,file,file,pr_allom,sig_xy,sig_z,allom,hgrad_w,zgrad_w,...
	db_eps, db_pts, min_pts, subSamp);
catch err
	err_handler(err,'seg1_err.txt');
	return
end


% test mcgc_double_layer_base_options with subsampling, layers same
rng(123);
try
	Seg2 = mcgc_double_layer_base_options(1,data_folder,out_folder_2,file,file,pr_allom,sig_xy,sig_z,allom,hgrad_w,zgrad_w,...
	db_eps, db_pts, min_pts, subSamp,1);
catch err
	err_handler(err,'seg2_err.txt');
	return
end


% Plot outputs

figure('Name','Single layer MCGC','NumberTitle','off');
las_class_plot(Seg1,7); % Column 7 as columns are: X,Y,Z,H,intensity,las input classification,MCGC label
title('Single layer MCGC');
figure('Name','Double layer MCGC','NumberTitle','off');
las_class_plot(Seg2,7);
title('Double layer MCGC');