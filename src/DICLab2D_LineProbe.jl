
function LineModule(Input)

    InputFolderPath=Input.InputFolderPath
    OutputFolderPath=Input.OutputFolderPath
    SubsetShape=Input.SubsetShape
    SubsetSize=Input.SubsetSize
    GaussSTD=Input.GaussSTD
    GaussWindow=Input.GaussWindow
    XY0=Input.LinesCoordinates
    Incremental=Input.Incremental

    ImgFiles = readdir(InputFolderPath,join=true)    #search for the files in the input folder
    nImages = length(ImgFiles)                       #compute the number of images in the folder
    nLines = size(XY0,1)                             #compute the number of line probes

    println("|  ",nImages, " images loaded","\n")
    println("|  ",nLines," Line probes added","\n")

    #= 
    
    In case Incremental is set to true, it is necessary to create a second subset position array XY to carry the updated subset coordinates. All calculations will be made using the XY array, while XY0 will stay unaltered. If Incremental is set to false, XY will not be updated and will remain unchanged and equal to XY0.
    
    =#

    XY = Float64.(XY0)

    Output = Array{Float64,3}(undef,nImages,21,nLines) #initialize Output array

    #write first lines in the output array

    for i in 1:nLines

        for j in 1:21

            FirstLine = [0 XY0[i,1] XY0[i,1] XY0[i,2] XY0[i,2] XY0[i,3] XY0[i,3] XY0[i,4] XY0[i,4] 0 0 0 0 abs(XY0[i,3]-XY0[i,1]) abs(XY0[i,3]-XY0[i,1]) 0 0 0 0 0 0]

            Output[1,j,i] = FirstLine[j]
                
        end

    end

    d = div(SubsetSize,2) #distance between the subset center to its edges

    dX = [i for j in 1:SubsetSize, i in -d:d]   #relative x coordinates for pixels in subset
    dY = [i for i in -d:d, j in 1:SubsetSize]   #relative y coordinates for pixels in subset

    #run dX and dY by the SubsetMaker function to be shaped

    dX = SubsetMaker(dX,SubsetSize,SubsetShape,d+1,d+1,d)
    dY = SubsetMaker(dY,SubsetSize,SubsetShape,d+1,d+1,d)

    println("----- Initializing correlation","\n")

    F = load(ImgFiles[1])           #load the reference image F
    F = Gray.(F)                    #convert the reference image to grayscale
    F = convert(Array{Float64},F)   #convert the reference image to an Float 64 array
    F = imfilter(F,Kernel.gaussian((GaussSTD,GaussSTD),(GaussWindow,GaussWindow))) #apply a gaussian filter on F

    #compute the partial derivatives of F using convolution

    DiffKernel = [1/12 -8/12 0 8/12 -1/12]          #convolution mask for partial derivatives

    dFdx = imfilter(F,DiffKernel,"replicate")   #convolve along rows
    dFdy = imfilter(F,DiffKernel',"replicate")  #convolve along columns

    G = similar(F)  #initialize a G array

    for img in 2:nImages #images loop
            
        println("Correlating deformed image ",img-1," of ",nImages-1,"\n")

        if Incremental == true && img > 2

            F = G   #update reference image to the last deformed image correlated

            dFdx = imfilter(F,DiffKernel,"replicate")   #update partial derivatives for F
            dFdy = imfilter(F,DiffKernel',"replicate")

            for i in 1:nLines

                XY[i,:] = XY0[i,:]+Output[img-1,10:13,i]    #update XY array

            end

        end

        G = load(ImgFiles[img])     #load the deformed image G 
        G = Gray.(G)
        G = convert(Array{Float64},G)
        G = imfilter(G,Kernel.gaussian((GaussSTD,GaussSTD),(GaussWindow,GaussWindow)))

        InterpolationTable = Bspline2D(G)               #compute de interpolation coefficients table for G

        for line in 1:nLines    #line probes loop
                
            #initial points coordinates

            x1i = XY[line,1]
            y1i = XY[line,2]
            x2i = XY[line,3]
            y2i = XY[line,4]

            #Adaptive subset offset scheme by ZHANG, 2014. If Incremental is set to false, coordinates are already integers and will not be changed

            x1i = round(Int,x1i)
            y1i = round(Int,y1i)
            x2i = round(Int,x2i)
            y2i = round(Int,y2i)

            #gray intensity values for subsets 1 and 2

            f1 = SubsetMaker(F,SubsetSize,SubsetShape,x1i,y1i,d)     
            f2 = SubsetMaker(F,SubsetSize,SubsetShape,x2i,y2i,d)

            #gradients in x direction for subsets 1 and 2

            dfdx1 = SubsetMaker(dFdx,SubsetSize,SubsetShape,x1i,y1i,d)
            dfdx2 = SubsetMaker(dFdx,SubsetSize,SubsetShape,x2i,y2i,d)

            #gradients in y direction for subsets 1 and 2

            dfdy1 = SubsetMaker(dFdy,SubsetSize,SubsetShape,x1i,y1i,d)
            dfdy2 = SubsetMaker(dFdy,SubsetSize,SubsetShape,x2i,y2i,d)

            #initial guess for the deformation parameters P vector for subsets 1 and 2

            if Incremental == true || img == 2
                    
                InitialGuess1 = PCMGuess(F,G,x1i,y1i,SubsetSize,d)
                InitialGuess2 = PCMGuess(F,G,x2i,y2i,SubsetSize,d)

            else

                #deformation parameters temporal transfer (from image to image)
                #in current implementation, only the displacement parameters u and v are transferred

                InitialGuess1[1],InitialGuess1[7] = Output[img-1,10,line],Output[img-1,11,line]
                InitialGuess2[1],InitialGuess2[7] = Output[img-1,12,line],Output[img-1,13,line]

            end

            #call SolveSubset function to correlate subsets 1 and 2

            P1,ZNCC1,nIter1,BadPOI = CorrelateSubset(f1,x1i,y1i,dX,dY,dfdx1,dfdy1,InterpolationTable,InitialGuess1,Input)
            P2,ZNCC2,nIter2,BadPOI = CorrelateSubset(f2,x2i,y2i,dX,dY,dfdx2,dfdy2,InterpolationTable,InitialGuess2,Input)
                
            #get u and v displacements from P1 and P2

            u1 = P1.u
            v1 = P1.v
            u2 = P2.u
            v2 = P2.v

            #if Incremental is set to false, the displacements u and v are already cumulative; otherwise, they are incremental and must be added to the previously accumulated displacements

            if Incremental == true
                    
                u1 = u1 + Output[img-1,10,line]    #accumulated u1 displacement
                v1 = v1 + Output[img-1,11,line]    #accumulated v1 displacement
                u2 = u2 + Output[img-1,12,line]    #accumulated u2 displacement
                v2 = v2 + Output[img-1,13,line]    #accumulated v2 displacement

            end

            #final points coordinates

            x1f = XY0[line,1]+u1
            y1f = XY0[line,2]+v1
            x2f = XY0[line,3]+u2
            y2f = XY0[line,4]+v2

            #distance between points 1 and 2

            dxi = abs(XY0[line,3]-XY0[line,1])
            dyi = abs(XY0[line,2]-XY0[line,4])
            dxf = abs(x2f-x1f)
            dyf = abs(y2f-y1f)

            #strain in x and y direction  

            ϵx = ((dxf-dxi)/dxi)
            ϵy = ((dyf-dyi)/dyi)

            #gather all data from current iteration
                
            LineData = [img-1 x1i x1f y1i y1f x2i x2f y2i y2f u1 v1 u2 v2 dxf dyf ϵx ϵy ZNCC1 ZNCC2 nIter1 nIter2]

            #write data to the Output array

            for j in eachindex(LineData)

                Output[img,j,line] = LineData[j]

            end

        end

    end

    println("----- End of correlation","\n")

    for line in 1:nLines

        println("Printing output file for line probe ",line," of ",nLines,"\n")

        OutputFile = string(OutputFolderPath,"\\DICLab2D_LineProbe_",Dates.format(Dates.now(),"dduyy_HHhMM"),"_LN",line,".CSV")
        OutputData = DataFrame(Output[:,:,line],["img","x1i","x1f","y1i","y1f","x2i","x2f","y2i","y2f","u1","v1","u2","v2","dx12","dy12","Strain_x","Strain_y","ZNCC1","ZNCC2","nIter1","nIter2"])
        CSV.write(OutputFile,OutputData, delim=";")

    end

end