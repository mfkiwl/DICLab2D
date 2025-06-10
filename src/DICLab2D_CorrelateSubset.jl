
function CorrelateSubset(
    
    f::Vector{Float64},
    x::Real,
    y::Real,
    dX::Vector{Int64},
    dY::Vector{Int64},
    dfdx::Vector{Float64},
    dfdy::Vector{Float64},
    InterpolationTable::Matrix{Matrix{Float64}},
    P::DeformationParameters,
    WarpUpdateFunctions::WarpUpdate,
    Input::InputParameters
    )

    #=

        This function correlates a subset around a reference POI.
        
    INPUT:

        f: the light intensity of the subset around the reference POI
        x: the x absolute coordinate of the POI (aka subset center pixel)
        y: the y absolute coordinate of the POI (aka subset center pixel)
        dX: the relative subset pixels x coordinates
        dY: the relative subset pixels y coordinates
        dfdx: the light intensity gradients for the subset in x direction
        dfdy: the light intensity gradients for the subset in y direction
        InterpolationTable: the bspline interpolation coefficients
        P: the initial guess fot the deformation parameters
        WarpUpdateFunctions: the warp update functions (Warp,StructToMat,MatToStruct,VecToStruct,StructToVec)
        Input: input parameters input

    OUTPUT:

        P: the deformation parameters vector (::DeformationParameters)
        ZNCC: the ZNCC coefficient for the subset (::Float64)
        nIterations: the number of iterations for convergence (::Int64)
        BadPOI: the corrupted POI flag (::Bool)

    =#

    UpdateStrategy = Input.UpdateStrategy
    ShapeFunctionOrder = Input.ShapeFunctionOrder
    StopCritValue = Input.StopCritValue
    MaxIterations = Input.MaxIterations
    FilterOutLowQualityPOI = Input.FilterOutLowQualityPOI
    MinZNCC = Input.MinZNCC

    fmean = mean(f)                    #mean of gray intensity values in current reference subset
    fnorm = sqrt(sum(@. (f-fmean)^2))    #normalization value for the current reference subset
    fnormed = @. (f-fmean)/fnorm    

    #pre-compute dfdwdp
                        
    if ShapeFunctionOrder == 0

        dfdWdP = [dfdx dfdy]

    elseif ShapeFunctionOrder == 1

        dfdWdP = @. [dfdx dX*dfdx dY*dfdx dfdy dX*dfdy dY*dfdy]
 
    elseif ShapeFunctionOrder == 2

        dfdWdP = @. [dfdx dX*dfdx dY*dfdx dX^2/2*dfdx dX*dY*dfdx dY^2/2*dfdx dfdy dX*dfdy dY*dfdy dX^2/2*dfdy dX*dY*dfdy dY^2/2*dfdy]

    elseif ShapeFunctionOrder == 3

        dfdWdP = @. [dfdx dX*dfdx dY*dfdx dX^2/2*dfdx dX*dY*dfdx dY^2/2*dfdx dX^3/6*dfdx dX^2/2*dY*dfdx dX*dY^2/2*dfdx dY^3/6*dfdx dfdy dX*dfdy dY*dfdy dX^2/2*dfdy dX*dY*dfdy dY^2/2*dfdy dX^3/6*dfdy dX^2/2*dY*dfdy dX*dY^2/2*dfdy dY^3/6*dfdy]

    elseif ShapeFunctionOrder == 4

        dfdWdP = @. [dfdx dX*dfdx dY*dfdx dX^2/2*dfdx dX*dY*dfdx dY^2/2*dfdx dX^4/24*dfdx dX^3/6*dY*dfdx dX^2/2*dY^2/2*dfdx dX*dY^3/6*dfdx dY^4/24*dfdx dfdy dX*dfdy dY*dfdy dX^2/2*dfdy dX*dY*dfdy dY^2/2*dfdy dX^4/24*dfdy dX^3/6*dY*dfdy dX^2/2*dY^2/2*dfdy dX*dY^3/6*dfdy dY^4/24*dfdy]

    end

    #define warp and auxiliary functions

    Warp = WarpUpdateFunctions.Warp
    StructToMat = WarpUpdateFunctions.StructToMat
    MatToStruct = WarpUpdateFunctions.MatToStruct
    VecToStruct = WarpUpdateFunctions.VecToStruct
    StructToVec = WarpUpdateFunctions.StructToVec

    StopCriterion = (dP) -> sqrt(dP.u^2+dP.v^2)  #stop criterion function
  
    Hinv = inv(dfdWdP'*dfdWdP)      #compute inverse of the Hessian matrix

    nIterations = 0                 #initialize iteration counter
    StopValue = 1.0                 #initialize stop value (1.0 is an arbitrary high value)
    ZNCC = 0.0                      #initialize correlation coefficient value

    BadPOI = false                  #initialize corrupted POI flag
    ContinueIteration = true        #initialize iteration flag

    dXY = Matrix{Float64}(undef,length(dX),2)   #initialize target subset relative pixel positions
    X   = Vector{Float64}(undef,length(dX))     #initialize target subset absolute pixel positions
    Y   = Vector{Float64}(undef,length(dY))     #initialize target subset absolute pixel positions

    while ContinueIteration == true
        
        dXY = Warp(dX,dY,P)     #target subset relative pixel positions
    
        @views begin

            @. X = x+dXY[:,1]         #target subset absolute pixel positions
            @. Y = y+dXY[:,2]

        end

        g,BadPOI = Interpolate(InterpolationTable,X,Y)     #interpolate target subset light intensity values

        if BadPOI == true

            #if the Interpolate function flags a corrupted POI, correlation stops and returns an all zeros deformation parameters vector

            return DeformationParameters(),0.0,0,true
            
        end
        
        gmean = mean(g)                         #mean of gray intensity values in current target subset
        gnorm = sqrt(sum(@. (g-gmean)^2))       #normalization value for the current target subset
        gnormed = @. (g-gmean)/gnorm

        #check if iteration should continue

        if any([StopValue < StopCritValue,nIterations == MaxIterations])

            ZNSSD = sum(@. (fnormed-gnormed)^2)

            ZNCC = 1.0-(ZNSSD/2.0)      #convert ZNSSD to ZNCC correlation coefficient

            #if FilterOutLowQualityPOI is set to true, check if the ZNCC for the POI attends the minimal required value

            if FilterOutLowQualityPOI == true && ZNCC < MinZNCC

                BadPOI = true
                
            end

            ContinueIteration = false   #flag iteration to stop

        else

            J = sum((@. dfdWdP*(f-fmean-((fnorm/gnorm)*(g-gmean))));dims=1)'   #compute the Jacobian

            if UpdateStrategy == "ICGN"

                dP = VecToStruct(-Hinv*J)                     #iterative improvement (struct)

                Pup = StructToMat(P)*inv(StructToMat(dP))     #P updated (matrix)
                P = MatToStruct(Pup)                          #P updated (struct)

            elseif UpdateStrategy == "BSGN"

                dP = -Hinv*J                #iterative improvement dP (vector)

                Pup = StructToVec(P)-dP     #P updated (vector)
                P = VecToStruct(Pup)        #P updated (struct)
                dP = VecToStruct(dP)        #dP (struct)

            end

            StopValue = StopCriterion(dP)   #compute the stop criterion value

        end

        nIterations += 1  #update iteration counter

    end
    
    return P,ZNCC,nIterations,BadPOI

end

#=NOTES:
            


=#