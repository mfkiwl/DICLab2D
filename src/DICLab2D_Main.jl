module DICLab2D

    using Images, ImageFiltering
    using Statistics, LinearAlgebra, FFTW
    using DataFrames, CSV, Dates

    include("DICLab2D_Structs.jl")
    include("DICLab2D_AreaProbe.jl")
    include("DICLab2D_LineProbe.jl")
    include("DICLab2D_CorrelateSubset.jl")
    include("DICLab2D_PointsOfInterest.jl")
    include("DICLab2D_InitialGuess.jl")
    include("DICLab2D_Strains.jl")
    include("DICLab2D_ShapeFunctions.jl")
    include("DICLab2D_Interpolations.jl")
    include("DICLab2D_Tables.jl")

    function RunProgram()

        #=

            Input analysis parameters and run this function to start the DICLab2D analysis.

        INPUT:

            Probe: the type of probe to be used in the analysis. Options are "Area" or "Line". (::String)
            InputFolderPath: the path to the folder containing the images. (::String)
            OutputFolderPath: the path to the folder where the output will be saved. (::String)
            UpdateStrategy: the update strategy to be used in the analysis. Options are "ICGN" or "BSGN". (::String)
            Incremental: a boolean indicating whether the analysis is incremental or not. (::Bool)
            SubsetShape: the shape of the subset to be used in the analysis. Options are "square" or "circle". (::String)
            SubsetSize: the size of the subset to be used in the analysis. Must always be a odd number. (::Int64)
            ShapeFunctionOrder: the order of the shape function to be used in the analysis. Options are 0, 1, 2, 3, or 4. (::Int64)
            StopCritValue: the value of the stop criterion to be used in the correlation. (::Float64)
            MaxIterations: the maximum number of iterations to be used in the correlation. (::Int64)
            GaussSTD: the standard deviation of the Gaussian filter to be applied to the images. (::Float64)
            GaussWindow: the size of the Gaussian filter to be applied to the images. (::Int64)
            MinimumZNCC: the minimum ZNCC value to be considered a valid point of interest. (::Float64)
            FilterOutLowQualityPOI: a boolean indicating whether to filter out low quality points of interest or not. (::Bool)

            Parameters exclusive to the area probe (set to nothing if not using area probe):

            ROIMaskPath: the path to the ROI mask image. (::String)
            ROIPadding: the padding to be applied to the ROI mask. (::Int64)
            GridStep: the step size of the grid to be used in the analysis. (::Int64)
            ComputeStrains: a boolean indicating whether to compute strains or not. (::Bool)
            StrainWindowSize: the size of the window to be used in the strain computation. Set to nothing if not computing strains. (::Int64)

            Parameters exclusive to the line probe (set to nothing if not using line probe):

            LinesCoordinates: the coordinates of the line probes to be used in the analysis. (::Matrix{Float64})



        =#

        #INPUT ANALYSIS PARAMETERS HERE:

            Probe = ""

            #General parameters

            InputFolderPath = ""
            OutputFolderPath = ""
            UpdateStrategy = ""
            Incremental = false
            SubsetShape = ""
            SubsetSize = 
            ShapeFunctionOrder = 
            StopCritValue = 1e-4
            MaxIterations = 20
            GaussSTD = 0.4
            GaussWindow = 5
            MinimumZNCC = 0.9
            FilterOutLowQualityPOI = false

            #Area probe parameters

            ROIMaskPath = "DIC Challenge//ROI Masks//Sample14FullMask.png"
            ROIPadding = 0
            GridStep = 5
            ComputeStrains = true
            StrainWindowSize = 5

            #Line module parameters

            LinesCoordinates = nothing

        #

        println("\n","------ Initializing DICLab2D Launcher","\n")

        Input = InputParameters(
            InputFolderPath,
            OutputFolderPath,
            ROIMaskPath,
            ROIPadding,
            UpdateStrategy,
            Incremental,
            SubsetShape,
            SubsetSize,
            GridStep,
            ShapeFunctionOrder,
            StopCritValue,
            MaxIterations,
            GaussSTD,
            GaussWindow,
            StrainWindowSize,
            LinesCoordinates,
            ComputeStrains,
            FilterOutLowQualityPOI,
            MinimumZNCC
        )

        if Probe == "Area"
            
            Data,Info = AreaProbe(Input)

            return Data, Info

        elseif Probe == "Line"

            println("------ Initializing Line Probe Probe","\n")

            LineProbe(Input)

        end

    end

end

#=NOTES:

1. The FilterOutLowQualityPOI parameter is used to filter out low quality points of interest based on the MinimumZNCC value, default value is false. POIs with ZNCC values below the MinimumZNCC will be considered corrupted if this parameter is set to true. Displacements and strains will not be computed for these POIs, and they will be excluded from the analysis for the following deformed images rounds. If this parameter is set to false, all POIs will be considered valid, even if their ZNCC values are below the MinimumZNCC value, with one exception: a POI can be labeled as corrupted if, in the correlation process, the shape function displaces any point of the subset outside the image boundaries. In this case, the POI will be considered corrupted and its displacements and strains will not be computed, even if FilterOutLowQualityPOI is set to false.

2. The MinimumZNCC parameter is also used to evaluate a valid seed point of interest, default value is 0.9.

3. Guidelines for Line probes coordinates:

- Multiple line probes can be defined in any direction;
- Line probes are specified based in their endpoints, arbitrarily referred as point 1 and point 2;
- For a line probe i, endpoints coordinates follow the notation:
                
point 1 (xi1,yi1) o---------o (xi2,yi2) point 2

- Line probe coordinates input are specified as an n×4 LinesCoordinates array, where n represents the number of probes;
                
- LinesCoordinates array is structured as follows:

LinesCoordinates = [x11 y11 x12 y12;x21 y21 x22 y22;...;xn1 yn1 xn2 yn2]

=#
