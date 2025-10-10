# DICLab2D
An in-development 2D Digital Image Correlation algorithm in Julia Language

HOW TO USE DICLab2D

DICLab2D provides two distinct analysis options: a line probe and a area probe.

The line probe is a method for measuring normal strain along user-defined line segments. It operates on the fundamental principle that the normal strain in a segment is given by the ratio of the change in its length by the initial length. Line probes are defined by two endpoints. DICLab2D tracks the sub-pixel displacement of these points across all deformed images using correlation. The normal strain is then computed directly from the change in the distance between these points. Multiple probes can be specified by an input matrix XY0 with dimensions of 4 rows x nLines columns, where nLines is the number of probes defined. In the XY0, each column vector [x1i y1i x2i y2i]^T contains the integer coordinates of the start (index 1) and end points (index 2) for a single probe i.

The area probe performs the standard full-field DIC analysis for a grid of points of interest (POIs) across a region of interest (ROI). The ROI is defined by an binary image mask, matching the dimensions of the images in the analysis, where pixels inside the ROI are set to 1 (white) and pixels outside the ROI to 0 (black). A padding parameter is required and creates an inner zone around the ROI borders where no points within the subset around POIs are placed to prevent subsets from extending beyond the mask. This padding can be disabled by setting the ROIPadding parameter to zero.

It is crucial that the images in the input folder are arranged in the correct order. The program will always consider the first image as the reference image and the remaining images as the deformed images.

The input parameters must be set in the RunProgram function in the DICLab_Main file. The parameters are given below.

--------------------------------------------------

Probe: the type of probe to be used in the analysis. Options are "Area" or "Line". (::String)

InputFolderPath: the path to the folder containing the images. (::String)

OutputFolderPath: the path to the folder where the output will be saved. (::String)

UpdateStrategy: the update strategy to be used in the analysis. Options are "ICGN" or "BSGN". (::String)

Incremental: a boolean indicating whether the analysis is incremental or not. (::Bool)

SubsetShape: the shape of the subset to be used in the analysis. Options are "square" or "circle". (::String)

SubsetSize: the size of the subset to be used in the analysis. Must always be a odd number. (::Int64)

ShapeFunctionOrder: the order of the shape function to be used in the analysis. Options are 0, 1, 2, 3, or 4. 3rd and 4th orders are only avaiable for BSGN update strategy. (::Int64)

StopCritValue: the value of the stop criterion to be used in the correlation. (::Float64)

MaxIterations: the maximum number of iterations to be used in the correlation. (::Int64)

GaussSTD: the standard deviation of the Gaussian filter to be applied to the images. (::Float64)

GaussWindow: the size of the Gaussian filter to be applied to the images. (::Int64)

MinimumZNCC: the minimum ZNCC value to be considered a valid point of interest. (::Float64)

FilterOutLowQualityPOI: a boolean indicating whether to filter out low quality points of interest or not. (::Bool)

-- Parameters exclusive to the area probe (set to nothing if not using area probe):

ROIMaskPath: the path to the ROI mask image. (::String)

ROIPadding: the padding to be applied to the ROI mask. (::Int64)

GridStep: the step size of the grid to be used in the analysis. (::Int64)

ComputeStrains: a boolean indicating whether to compute strains or not. (::Bool)

StrainWindowSize: the size of the window to be used in the strain computation. Set to nothing if not computing strains. (::Int64)

-- Parameters exclusive to the line probe (set to nothing if not using line probe):

LinesCoordinates: the XY0 matrix with coordinates of the line probes to be used in the analysis. (::Matrix{Float64})

--------------------------------------------------

Fill in the input parameters and execute the RunProgram function to run the analysis.
