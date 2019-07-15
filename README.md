# MCGC

### *Jonathan Williams*
### *last updated 15/07/2019*

MCGC is a collection of scripts to run the Multi-Class Graph Cut (MCGC) approach to individual tree crown (ITC) detection from Airborne LiDAR data (ALS). We implemented the approach in MATLAB, with some steps requiring the use of other work in R. The exact details of the algorithm are explained in detail in our [manuscript introducing this algorithm][1]. Broad steps are outlined below, but the individual scripts contain explanations of the arguments (though where these are inherited from a dependency descriptions may not be copied for the sake of brevity). I hope to include a more detailed description of the steps here in due course, however the user is highly recommended to refer to our manuscript to understand the motivation and steps of our process, including a description of the key parameters.


# Broad Outline

The key steps and ieas of the MCGC approach are outlined below:

0. *Pre-processing* - (not automated in pipeline) tidies the LiDAR data (.las files) and computes the topographcially corrected point cloud. This uses lasTools and is implemented as a series of batch files.
1. *Prior Generation* - uses an allometrically scaling local window maxima finder to estimate the number of tree tops (this is done in R, using the approach outlined in [Coomes et al (2017)][3] which is implemented as the first step in the R package [itcSegment][4]).
2. *Graph Cut* - computes candidate tree crowns in 3D by the use of a multiclass normalised cut of a graph. The graph is constructed from the point cloud based on geometry, and the local density for each point.
3. *Allometric Checking* - allometry of each candidate crown is checked in terms of overlap with other crowns, extent when compared to allometric limits and then connectivity and minimum size. Crowns are either accepted or rejected

*Optional Step*

4. *Second Layer* - apply the pipeline from 1-3 a second time to those points which were rejected in step 3 of the first pass

It is also possible to complete the full processing on a subsmapled point cloud, with the final reuslts imputed by k-nearest neighbour imputation (within the allowed allometric radius for each crown). This is implemented both for the basic algorithm, and as a part of the double layer approach. This approach speeds up the computation a lot and is highly recommended for large datasets, or where point density is high. A stability and speed analysis is presented in our manuscript.

# Which script do I want?

A brief guide of which function to use for your requirements is outlined below. If the docstring for that function doesn't give a detailed description of any parameter then use the dependency tree in that docstring to find the next level of the dependency tree and look at that docstring. Evenrtually you will be pointed to the function where a parameter is used, where there should be a detailed explanation of its role. In future I hope to include details of the key parameters in this vignette.

| Function 	| Explanation |
|:---------------------:|:-----------------------------------------------------:|
|`mcgc_base_options`| Run single Layer MCGC with most options pre-set |
|`mcgc_double_layer_base_options`| Run double Layer MCGC with most options pre-set |
|`las_class_plot`| Plot point cloud coloured by crown label |

**To get greater control of workings** not recommended unless really desired (and having read the manuscript)

| Function 	| Explanation |
|:---------------------:|:-----------------------------------------------------:|
|`mcgc`| Run single Layer MCGC with access to options but with no subsampling |
|`mcgc_subsample`| As above but incluing subsampling |
|`mcgc_double_layer`| Run double Layer MCGC with access to options with the choice to add subsampling or not |

### Pre-processing

These are found in the batch sub-directory of the scripts directory

| Function 	| Explanation |
|:---------------------:|:-----------------------------------------------------:|
|`gc_preproc.bat`| Remove noise points and generate a topogrpahically corrected point cloud |
|`gc_preproc_batch.bat`| As above but run this for all .las files in a directory |

These are found in the R sub-directory of the scripts directory

| Function 	| Explanation |
|:---------------------:|:-----------------------------------------------------:|
|`gen_allom.R`| Generate allometric lookup tables - use [data][5] from [Jucker et al 2017][6] |
|`utils.R`| Use function `rq_lut` to build an allometric lookup table for a range of tau in quantile regression |

### Useful utilities

These are found in the R sub-directory of the scripts directory

| Function 	| Explanation |
|:---------------------:|:-----------------------------------------------------:|
|`Pts2TreeShp.R`| Convert point cloud with crown labels to a shapefile of crown convex hulls with heights |


# Requirements

The following are requirements of your system and/or environment to be able to run the MCGC package.

### System

* MATLAB (this was developed in R2018a, so reliability in other releases is not guaranteed)
* R - must be on system path so RScript can be used (again developed in R 3.4.4)
* Must be running Windows (to run LASTools)
* A valid installation of [LASTools][2] which must be on the system path
  * The above 2 requirements can be removed with some work on the users behalf (to replace system calls to lastools in the MATLAB code with any other software with the same capability)


### R packages (with version developed using listed)
*NB these should be automatically installed by the scripts if not already present when these are run*

* `quantreg 5.38`
* `raster 2.7-15`
* `rgdal 1.3-6`
* `sp 1.3-1` 


# A worked (artificial) example

I will later add details of pre-processing and how to compute your own allometric look up tables, though the scripts for these processes are listed above and should have relatively clear docstrings (I hope).

First you will need to download a local copy of this repository. Then you should open MATLAB and ensure that the 'scripts' directory is added t the MATLAB path which you can do with the below command in MATLAB:

`addpath(genpath('MCGC/scripts'));`

Once you have done this you can then load in the allometric lookup table for allometric testing as shown for our example:

`allom = csvread('indomalaya_tropfor_HR_95.csv');`

You will then need to set the key variables for the base option. The example script 'example_pipeline.m' shows an example of these. The function call will then look like the below:

`output = mcgc_base_options(data_folder, out_folder, raw_file, chm_file, pr_allom, sigxy, sigz, allom, hgrad_w, zgrad_w, db_eps, db_pts, min_pts, subSamp);`

Here we are using the single layer MCGC. `data_folder` and `output_folder` tell the location data is stored, and where to save output. `raw_file` should be the base las file name as string without .las on the end `chm_file` should be the topographically corrected las file name 9again wihout .las, `pr_allom` should be the file name of the prior look up table (a table of tau, alpha and beta for quantile regression of CD as a power law of H - see the manuscript for more detail), `sigxy` and `sigz` are parameters for the weights in the graph based on distance between points and `hgrad_w` and `hgrad_z` are the parameters for the local density effects in the weights (see manuscript/docstrings for more details on weight formula), `allom` is a pre-loaded table of CD as a function of H for use in segmentation and allometric testing, `db_eps` and `db_pts` are the neighbourhood size and number of points to use in connectivity testing using DBSCAN, `min_pts` is the minimum number of point sin each crown and `subSamp` is the factor by which to subsample the point cloud (can be set to 0 to avoid subsampling)

`output = mcgc_double_layer_base_options(1,data_folder, out_folder, raw_file, chm_file, pr_allom, sig_xy, sig_z, allom, hgrad_w, zgrad_w, db_eps, db_pts, min_pts, subSamp,1);`

This does the same as the baove, but uses a double layer approach. Here only two arguments are added. The very fist argument tells MCGC whether to save intermediate steps between the layers. This is set to 1 here to mean save the results, but setting it to 0 gives only the final results, with no saving after the first layer. The final argument determines if the layers should have the same settings. Here 1 means to use the same settings, but if this were set to 0 then additional setting would have to follow this argument as further arguments, the details of which are in the docstring for this function



## References
[**Three-dimensional Segmentation of Trees Through a Flexible Multi-Class Graph Cut Algorithm (MCGC)**][1] Williams, J; Sch&ouml;nlieb, C-B; Swinfield, T; Lee, J; Cai, X; Qie, L and Coomes, D. ArXiv preprint (2019) (in revisions)
[**Area-based vs tree-centric approaches to mapping forest carbon in Southeast Asian forests from airborne laser scanning data**][3] Coomes, D; Dalponte, M; Jucker, T; Asner, G; Banin, L; Burslem, D; Lewis, S; Nilus, R; Philips, O; Phua, M-H and Qie, L. Remote Sensing of Environment (2017)
[**Allometric equations for integrating remote sensing imagery into forest monitoring programmes**][6] Jucker, T. , Caspersen, J. , Chave, J. , Antin, C. , Barbier, N. , Bongers, F. , Dalponte, M. , Ewijk, K. Y., Forrester, D. I., Haeni, M. , Higgins, S. I., Holdaway, R. J., Iida, Y. , Lorimer, C. , Marshall, P. L., Momo, S. , Moncrieff, G. R., Ploton, P. , Poorter, L. , Rahman, K. A., Schlund, M. , Sonké, B. , Sterck, F. J., Trugman, A. T., Usoltsev, V. A., Vanderwel, M. C., Waldner, P. , Wedeux, B. M., Wirth, C. , Wöll, H. , Woods, M. , Xiang, W. , Zimmermann, N. E. and Coomes, D. A. Global Change Biology (2017)


[1]: https://arxiv.org/abs/1903.08481
[2]: https://rapidlasso.com/lastools/
[3]: https://www.sciencedirect.com/science/article/pii/S0034425717301098#!
[4]: https://cran.r-project.org/web/packages/itcSegment/index.html
[5]: https://figshare.com/articles/Global_Allometric_Database/3413539/1
[6]: https://onlinelibrary.wiley.com/doi/full/10.1111/gcb.13388