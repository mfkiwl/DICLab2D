
function AreaProbe(Input)
    
    #=

        This function executes a RGDT area probe analysis
        
    INPUT:

        Input: input parameters struct

    OUTPUT:

        Data: output data struct
        Info: analysis information struct

    =#

    println("------ Initializing Area Probe Analysis","\n")

    InputFolderPath = Input.InputFolderPath
    ROIMaskPath = Input.ROIMaskPath
    ROIPadding = Input.ROIPadding
    UpdateStrategy = Input.UpdateStrategy
    ShapeFunctionOrder = Input.ShapeFunctionOrder
    Incremental = Input.Incremental
    SubsetShape = Input.SubsetShape
    SubsetSize = Input.SubsetSize
    GridStep = Input.GridStep
    GaussSTD = Input.GaussSTD
    GaussWindow = Input.GaussWindow
    ComputeStrains = Input.ComputeStrains
    StrainWindowSize = Input.StrainWindowSize
    MinZNCC = Input.MinZNCC

    ImgFiles = readdir(InputFolderPath,join=true)   #read files in the input folder

    nImages = length(ImgFiles)                      #number of images in the folder
    nRows,nCols = size(load(ImgFiles[1]))           #size of the images in pixels

    println("|  ",nImages," images loaded","\n")

    ROIMask = load(ROIMaskPath)             #loads the ROI mask image              
    ROIMask = Gray.(ROIMask)                #converts to grayscale
    ROIMask = Float64.(ROIMask .> 0.9)      #converts to binary (black & white)

    X0,Y0 = DefinePOI(ROIMask,ROIPadding,SubsetSize,SubsetShape,GridStep)   #initial POI coordinates
        
    nPOI = length(X0)                       #number of POI
    nValidPOI = nPOI                        #number of valid POI
    ValidPOI = trues(nValidPOI)             #valid POI vector

    println("|  ",nPOI," POI per image","\n")

    POIMap = zeros(Int,nRows,nCols)         #initialize POI map

    for i in 1:nPOI

        POIMap[Y0[i],X0[i]] = i             #label every POI in the POI map
            
    end

    Data = [OutputData() for i in 1:nImages, j in 1:nPOI]   #initialize output data array
    Info = AnalysisInfo(nImages,nRows,nCols,nPOI,POIMap)    #assemble the analysis information struct

    for i in 1:nPOI

        Data[1,i] = OutputData(x=X0[i],y=Y0[i])             #initialize initial POI coordinates in the Data array
 
    end

    #local relative coordinates for pixels in the subset

    d = div(SubsetSize,2)

    dX = [i for j in 1:SubsetSize, i in -d:d]
    dY = [i for i in -d:d, j in 1:SubsetSize]

    #run dX and dY by the SubsetMaker function to be shaped

    dX = ShapeSubset(dX,SubsetSize,SubsetShape,d+1,d+1,d)
    dY = ShapeSubset(dY,SubsetSize,SubsetShape,d+1,d+1,d)

    if ComputeStrains == true

        #local relative coordinates for points in the strain window

        dw = div(StrainWindowSize,2)

        dXw = [i for j in 1:StrainWindowSize, i in -dw*GridStep:GridStep:dw*GridStep][:]      
        dYw = [i for i in -dw*GridStep:GridStep:dw*GridStep, j in 1:StrainWindowSize][:]

    end

    println("----- Initializing correlation","\n")

    F = load(ImgFiles[1])                                                               #load reference image
    F = Gray.(F)                                                                        #convert reference image to grayscale
    F = convert(Array{Float64},F)                                                       #convert reference image to a float array
    F = imfilter(F,Kernel.gaussian((GaussSTD,GaussSTD),(GaussWindow,GaussWindow)))      #apply gaussian filter on F

    #compute the partial derivatives of F using convolution

    DiffKernel = [1/12 -8/12 0 8/12 -1/12]          #convolution mask for partial derivatives

    dFdx = imfilter(F,DiffKernel,"replicate")       #convolve along rows
    dFdy = imfilter(F,DiffKernel',"replicate")      #convolve along columns

    G = similar(F)      #initialize deformed image array

    #define warp and parameter update functions
    
    WarpUpdateFunctions = DefineWarpUpdateFunctions(ShapeFunctionOrder,UpdateStrategy)

    for Img in 2:nImages    #deformed images loop

        println("Correlating deformed image ",Img-1," of ",nImages-1,"\n")

        ProcessingQueue = Queue[]                       #initialize the processing queue
        UnprocessedValidPOI = copy(ValidPOI)            #initialize the unprocessed valid POI array

        if Incremental == true && Img > 2

            F = G   #update reference image to the last deformed image

            dFdx = imfilter(F,DiffKernel,"replicate")   #update partial derivatives for F
            dFdy = imfilter(F,DiffKernel',"replicate")

        end

        G = load(ImgFiles[Img])                                                             #load deformed image G
        G = Gray.(G)                                                                        #convert deformed image to grayscale
        G = convert(Array{Float64},G)                                                       #convert deformed image to a float array
        G = imfilter(G,Kernel.gaussian((GaussSTD,GaussSTD),(GaussWindow,GaussWindow)))      #apply gaussian filter on G

        InterpolationTable = Bspline2D(G)               #compute the interpolation coefficients table for G

        while iszero(UnprocessedValidPOI) == false

            #seed POI section
    
                POI = findfirst(UnprocessedValidPOI)   #randomly select a seed POI from the unprocessed valid POI vector

                x0 = X0[POI]        #seed POI x coordinate
                y0 = Y0[POI]        #seed POI y coordinate

                if Img == 2

                    InitialGuess = PCMGuess(F,G,x0,y0,SubsetSize,d)     #compute initial guess

                else

                    InitialGuess = DeformationParameters(u=Data[Img-1,POI].u,v=Data[Img-1,POI].v)   #get initial guess from the previous deformed image
                    
                end
            
                ZNCC = SolvePOI!(Img,POI,x0,y0,F,dFdx,dFdy,d,dX,dY,InterpolationTable,InitialGuess,Input,Data,ProcessingQueue,ValidPOI,nValidPOI,UnprocessedValidPOI,WarpUpdateFunctions)

                if nValidPOI == 0

                    break       #break the loop if there are no valid POIs left to process
                    
                end

                if ZNCC < MinZNCC

                    continue    #restart the loop if the ZNCC is below the minimum threshold

                end

            #

            while iszero(UnprocessedValidPOI) == false

                if length(ProcessingQueue) == 0

                    #if the processing queue is empty, break the loop and return to the seed POI section

                    break
                    
                end

                sort!(ProcessingQueue,by=set->set.ZNCC)         #sort ProcessingQueue in a crescent ZNCC order

                MostReliablePOI = pop!(ProcessingQueue)         #collect and remove the most reliable POI (highest ZNCC) from the processing queue
                POI = MostReliablePOI.POI                       #POI label
                InitialGuess = MostReliablePOI.P                #define the deformation parameters from the POI as initial guess

                x0 = X0[POI]        #POI initial x coordinate
                y0 = Y0[POI]        #POI initial y coordinate

                #define the top, bottom, left, and right neighbors

                Neighbors = [x0 y0-GridStep;x0 y0+GridStep;x0-GridStep y0;x0+GridStep y0]

                for i in 1:4    #neighbors loop

                    x0,y0 = Neighbors[i,1],Neighbors[i,2]       #current neighbor coordinates

                    POI = POIMap[y0,x0]                         #current neighbor POI label

                    if POI !== 0 && UnprocessedValidPOI[POI] == true       #check if current neighbor is a valid POI

                        SolvePOI!(Img,POI,x0,y0,F,dFdx,dFdy,d,dX,dY,InterpolationTable,InitialGuess,Input,Data,ProcessingQueue,ValidPOI,nValidPOI,UnprocessedValidPOI,WarpUpdateFunctions)
                                    
                    end

                end

            end

        end

        if ComputeStrains == true   #check if strains are to be computed

            #collect data for the current deformed image

            X = [i.x for i in Data[Img,:]]
            Y = [i.y for i in Data[Img,:]]
            U = [i.u for i in Data[Img,:]]
            V = [i.v for i in Data[Img,:]]

            #compute strains for current deformed image

            Exx,Eyy,Exy,E1,E2 = CalculateStrains(X0,Y0,X,Y,U,V,dw,dXw,dYw,Input,Info,ValidPOI)

            for i in 1:nPOI

                #write strain data
        
                Data[Img,i].Exx = Exx[i]
                Data[Img,i].Eyy = Eyy[i]
                Data[Img,i].Exy = Exy[i]
                Data[Img,i].E1 = E1[i]
                Data[Img,i].E2 = E2[i]
                            
            end

        end

        if nValidPOI == 0

            break   #break the loop if there are no valid POIs left to process
            
        end

    end

    PrintTable(Input,Data,Info)     #print .CSV file with analysis data

    return Data,Info

end

function SolvePOI!(Img,POI,x0,y0,F,dFdx,dFdy,d,dX,dY,InterpolationTable,InitialGuess,Input,Data,ProcessingQueue,ValidPOI,nValidPOI,UnprocessedValidPOI,WarpUpdateFunctions)
    
    SubsetSize = Input.SubsetSize
    SubsetShape = Input.SubsetShape
    Incremental = Input.Incremental

    if Incremental == true && Img > 2

        xcurrent = Data[Img-1,POI].x    #POI current x coordinate 
        ycurrent = Data[Img-1,POI].y    #POI current y coordinate 
            
        x = round(Int,xcurrent)         #POI shifted x coordinate
        y = round(Int,ycurrent)         #POI shifted y coordinate

    else

        x = x0
        y = y0
            
    end
        
    f = ShapeSubset(F,SubsetSize,SubsetShape,x,y,d)           #gray-level values in current subset
    dfdx = ShapeSubset(dFdx,SubsetSize,SubsetShape,x,y,d)     #gradients in current subset in the x direction
    dfdy = ShapeSubset(dFdy,SubsetSize,SubsetShape,x,y,d)     #gradients in current subset in the y direction

    #correlate POI

    P,ZNCC,nIter,BadPOI = CorrelateSubset(f,x,y,dX,dY,dfdx,dfdy,InterpolationTable,InitialGuess,WarpUpdateFunctions,Input)

    if BadPOI == true   #check if POI is flagged as a corrupted POI

        #write data for corrupted POI

        Data[Img,POI].x = NaN
        Data[Img,POI].y = NaN
        Data[Img,POI].u = NaN
        Data[Img,POI].v = NaN
        Data[Img,POI].ZNCC = ZNCC
        Data[Img,POI].nIter = nIter

        nValidPOI -= 1                          #update valid POI counter
        ValidPOI[POI] = false                   #flag POI as corrupted
        UnprocessedValidPOI[POI] = false        #flag POI as processed

    else

        u = P.u     #u displacement
        v = P.v     #v displacement

        if Incremental == true && Img > 2
            
            u = u + Data[Img-1,POI].u       #accumulated u displacement
            v = v + Data[Img-1,POI].v       #accumulated v displacement

            Data[Img,POI].x = xcurrent+u    #write POI updated x coordinate
            Data[Img,POI].y = ycurrent+v    #write POI updated y coordinate

        else

            Data[Img,POI].x = x+u    #write POI updated x coordinate
            Data[Img,POI].y = y+v    #write POI updated y coordinate

        end

        #write POI data to the output data array

        Data[Img,POI].u = u
        Data[Img,POI].v = v
        Data[Img,POI].ZNCC = ZNCC
        Data[Img,POI].nIter = nIter

        push!(ProcessingQueue,Queue(ZNCC,POI,P))    #add POI to the processing queue
        UnprocessedValidPOI[POI] = false            #flag POI as processed

    end

    return ZNCC

end